#define TINYMESH_API_EXPORT
#include "restore.h"

#include <iostream>
#include <fstream>
#include <queue>
#include <functional>
#include <unordered_map>

#include <Eigen/LU>

#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"
#include "core/bvh.h"
#include "core/eigen.h"
#include "core/openmp.h"

#include "ops/ops.h"

namespace tinymesh {

namespace {

double solveRigidICP(const std::vector<Vertex *> tgtPatch, const std::vector<Vertex *> srcPatch, const EigenMatrix &hks,
                     EigenMatrix3 &R, EigenVector3 &t) {
    const int nTgtSize = (int)tgtPatch.size();
    const int nSrcSize = (int)srcPatch.size();
    EigenMatrix X(nSrcSize, 3);
    EigenMatrix Y(nSrcSize, 3);

    // Find closest point
    X.setZero();
    Y.setZero();
    for (int i = 0; i < nSrcSize; i++) {
        double minDist = 1.0e20;
        int minId = -1;
        for (int j = 0; j < nTgtSize; j++) {
            const int srcId = srcPatch[i]->index();
            const int tgtId = tgtPatch[j]->index();
            const double dist = (hks.row(srcId) - hks.row(tgtId)).squaredNorm();
            if (dist < minDist) {
                minDist = dist;
                minId = j;
            }
        }

        const Vec3 src = srcPatch[i]->pos();
        const Vec3 tgt = tgtPatch[minId]->pos();
        X.row(i) << src.x(), src.y(), src.z();
        Y.row(i) << tgt.x(), tgt.y(), tgt.z();
    }

    // Solve rigid transformation
    const EigenVector xMean = X.colwise().mean();
    const EigenVector yMean = Y.colwise().mean();
    X.rowwise() -= xMean.transpose();
    Y.rowwise() -= yMean.transpose();

    const EigenMatrix M = Y.transpose() * X;
    EigenMatrix U, V;
    EigenVector sigma;
    eigenSVD(M, U, sigma, V);

    const double detVU = (U * V.transpose()).determinant();

    // Revise sign of determinant
    EigenVector diagH(3);
    diagH << 1.0, 1.0, detVU;

    // Output
    R = U * diagH.asDiagonal() * V.transpose();
    t = yMean - R * xMean;

    // Error
    EigenMatrix Xt = X * R.transpose();
    Xt.rowwise() += t.transpose();

    const double error = (Y - Xt).squaredNorm() / (double)nSrcSize;
    return error;
}

double patchDissimilarity(const std::vector<std::vector<Vertex *>> &patches, const EigenMatrix &hks,
                          const std::vector<int> tgtCands, const std::vector<int> &srcCands) {
    const int nCands = tgtCands.size();
    double dissim = 0.0;
    EigenMatrix3 R_unused;
    EigenVector3 t_unused;
    for (int tgtId : tgtCands) {
        double minDist = 1.0e20;
        for (int srcId : srcCands) {
            const std::vector<Vertex *> tgtPatch = patches[tgtId];
            const std::vector<Vertex *> srcPatch = patches[srcId];
            const double dist = solveRigidICP(tgtPatch, srcPatch, hks, R_unused, t_unused);
            if (dist < minDist) {
                minDist = dist;
            }
        }
        dissim += minDist / (double)nCands;
    }
    return dissim;
}

}  // namespace

void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound) {
    mesh.holeFillMinDihedral_(face, dihedralBound);
}

void holeFillAdvancingFront(Mesh &mesh, Face *face) {
    mesh.holeFillAdvancingFront_(face);
}

void holeFillContextCoherent(Mesh &mesh, int patchRadius, int maxiters) {
    // 0. lock non-hole vertices
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        mesh.vertex(i)->lock();
    }

    // 1. initial hole filling (with advancing front)
    std::vector<Face *> holeFaces;
    for (int i = 0; i < mesh.numFaces(); i++) {
        Face *f = mesh.face(i);
        if (f->isHole()) {
            holeFillAdvancingFront(mesh, f);
        }
    }
    mesh.save("1_initial_hole_fill.ply");

    // 2. Compute HKS and feature line fields
    const EigenSparseMatrix L = getMeshLaplacian(mesh, MeshLaplace::Cotangent);
    const EigenMatrix hks = getHeatKernelSignatures(L);

    const auto flff = getFeatureLineFieldWithFlags(mesh, true);
    const EigenMatrix &flf = std::get<0>(flff);
    const EigenVector &rvFlags = std::get<1>(flff);

    // 3. Collect patches
    const int N = (int)mesh.numVertices();
    std::vector<int> srcIds;
    std::vector<int> tgtIds;
    std::vector<std::vector<Face *>> patchFaces(N);
    for (int i = 0; i < N; i++) {
        Vertex *ctr = mesh.vertex(i);
        std::unordered_set<Face *> visited;
        std::queue<std::pair<Face *, int>> que;
        for (auto it = ctr->f_begin(); it != ctr->f_end(); ++it) {
            que.push(std::make_pair(it.ptr(), 1));
        }

        while (!que.empty()) {
            Face *f = que.front().first;
            int dist = que.front().second;
            que.pop();

            if (visited.count(f) != 0) continue;

            visited.insert(f);
            for (auto it = f->f_begin(); it != f->f_end(); ++it) {
                if (visited.count(it.ptr()) == 0 && dist + 1 <= patchRadius) {
                    que.push(std::make_pair(it.ptr(), dist + 1));
                }
            }
        }

        patchFaces[i].assign(visited.begin(), visited.end());
        if (ctr->isLocked()) {
            srcIds.push_back(i);
        } else {
            tgtIds.push_back(i);
        }
    }
    const int nSrc = (int)srcIds.size();
    const int nTgt = (int)tgtIds.size();

    // Vertex to patch table
    std::vector<std::vector<Vertex *>> patchVerts(N);
    std::vector<std::unordered_set<int>> patchVertSets(N);
    for (int i = 0; i < N; i++) {
        for (Face *f : patchFaces[i]) {
            for (auto it = f->v_begin(); it != f->v_end(); ++it) {
                patchVertSets[i].insert(it->index());
            }
        }

        for (auto vid : patchVertSets[i]) {
            patchVerts[i].push_back(mesh.vertex(vid));
        }
    }

    std::unordered_map<int, std::vector<int>> tgtPatchAssign;
    for (int i = 0; i < nTgt; i++) {
        for (int k : patchVertSets[tgtIds[i]]) {
            if (tgtPatchAssign.count(k) == 0) {
                tgtPatchAssign[k] = std::vector<int>();
            }
            tgtPatchAssign[k].push_back(tgtIds[i]);
        }
    }

    // Compute patch descriptors
    std::vector<EigenVector> descriptors(N);
    for (int i = 0; i < N; i++) {
        EigenVector hksMean(hks.cols());
        EigenVector hksSigma(hks.cols());

        hksMean.setZero();
        hksSigma.setZero();
        for (auto v : patchVerts[i]) {
            const EigenVector f = hks.row(v->index());
            hksMean += f;
            hksSigma += f * f;
        }

        hksMean /= (double)N;
        hksSigma /= (double)N;
        hksSigma = (hksSigma - hksMean * hksMean).array().max(0.0);
        hksSigma = hksSigma.cwiseSqrt();

        hksMean /= hksMean.maxCoeff();
        hksSigma /= hksSigma.maxCoeff();

        descriptors[i].resize(hksMean.rows() + hksSigma.rows());
        descriptors[i] << hksMean, hksSigma;
    }

    // 4. Collect candidate patches
    const int nCands = std::max(5, (int)(0.001 * N));
    std::vector<std::vector<int>> candidates(N, std::vector<int>(nCands, 0));
    for (int i = 0; i < N; i++) {
        std::vector<std::pair<double, int>> dists;
        for (int j = 0; j < N; j++) {
            const double d = (descriptors[i] - descriptors[j]).norm();
            dists.emplace_back(d, j);
        }

        std::sort(dists.begin(), dists.end());
        for (int j = 0; j < nCands; j++) {
            candidates[i][j] = dists[j + 1].second;
        }
    }

    for (int kIter = 0; kIter < maxiters; kIter++) {
        printf("iter = %d\n", kIter);

        // 5. Compute most similar patch for each target patch
        std::vector<int> pairIds(nTgt, 0);
        omp_parallel_for(int i = 0; i < nTgt; i++) {
            double minEps = 1.0e20;
            int minId = -1;
            for (int j = 0; j < nSrc; j++) {
                const int tgtId = tgtIds[i];
                const int srcId = srcIds[j];
                const double eps = patchDissimilarity(patchVerts, hks, candidates[tgtId], candidates[srcId]);
                if (eps < minEps) {
                    minEps = eps;
                    minId = srcId;
                }
            }
            pairIds[i] = minId;
        }

        // 6. Compute rigid registration
        std::vector<EigenMatrix3> rotations(nTgt);
        std::vector<EigenVector3> translations(nTgt);
        std::unordered_map<int, double> dissimilarities(nTgt);
        for (int i = 0; i < nTgt; i++) {
            const std::vector<Vertex *> tgtPatch = patchVerts[tgtIds[i]];
            const std::vector<Vertex *> srcPatch = patchVerts[pairIds[i]];
            solveRigidICP(tgtPatch, srcPatch, hks, rotations[i], translations[i]);

            double dissim = 0.0;
            double sumWgt = 0.0;
            for (Vertex *vt : tgtPatch) {
                Vertex *closest = nullptr;
                double minDist = 1.0e20;
                for (Vertex *vs : srcPatch) {
                    const Vec3 pt = vt->pos();
                    const Vec3 ps = (Vec3)(rotations[i] * vs->pos()) + (Vec3)translations[i];
                    const double dist = length(pt - ps);
                    if (dist < minDist) {
                        minDist = dist;
                        closest = vs;
                    }
                }

                const EigenVector ft_ = flf.row(vt->index());
                const EigenVector fs_ = flf.row(closest->index());

                const Vec3 ft(ft_(0), ft_(1), ft_(2));
                const Vec3 fs(fs_(0), fs_(1), fs_(2));

                double w = 1.0;
                if (rvFlags(closest->index()) != 0.0) {
                    w = (double)tgtPatch.size();
                }
                double l = length(ft - (Vec3)(rotations[i] * fs));
                dissim += w * l * l;
                sumWgt += w;
            }
            dissimilarities[tgtIds[i]] = dissim / sumWgt;
        }

        // Patch BVHs
        std::unordered_map<int, BVH> patchBVHs(nTgt);
        for (int i = 0; i < nTgt; i++) {
            std::vector<Vec3> vertices;
            std::vector<uint32_t> indices;
            std::unordered_map<Vec3, uint32_t> uniqueVertices;
            for (auto f : patchFaces[pairIds[i]]) {
                for (auto it = f->v_begin(); it != f->v_end(); ++it) {
                    const Vec3 p = (Vec3)(rotations[i] * it->pos()) + (Vec3)translations[i];
                    if (uniqueVertices.count(p) == 0) {
                        const uint32_t idx = vertices.size();
                        vertices.push_back(p);
                        uniqueVertices[p] = idx;
                    }
                    indices.push_back(uniqueVertices[p]);
                }
            }
            Assertion(indices.size() % 3 == 0, "Non-triangle face detected!");
            patchBVHs[tgtIds[i]] = BVH(vertices, indices);
        }

        // Update vertex positions
        std::unordered_map<int, Vec3> newPoss;
        std::unordered_map<int, Vec3> offsets;
        for (auto it : tgtPatchAssign) {
            const int tgtId = it.first;
            const std::vector<int> &neighbors = it.second;
            const Vec3 orgPos = mesh.vertex(tgtId)->pos();

            Vec3 newPos(0.0);
            double sumWgt = 0.0;
            for (int k : neighbors) {
                if (patchBVHs.count(k) == 0) continue;

                const BVH &bvh = patchBVHs[k];
                const Vec3 closest = bvh.closestPoint(orgPos);
                const double dissim = dissimilarities[k];
                const double weight = 1.0 / (std::pow(dissim, 5.0) + 1.0e-12);
                newPos += weight * closest;
                sumWgt += weight;
            }
            newPos /= sumWgt;

            if (mesh.vertex(tgtId)->isLocked()) {
                // Boundary vertices
                newPoss[tgtId] = Vec3(0.0);
                offsets[tgtId] = orgPos - newPos;
            } else {
                // Patch vertices
                newPoss[tgtId] = newPos;
                offsets[tgtId] = Vec3(0.0);
            }
        }

        // Smooth offsets
        for (int lapIters = 0; lapIters < 1000; lapIters++) {
            for (int i = 0; i < nTgt; i++) {
                const int tgtId = tgtIds[i];
                const Vertex *v = mesh.vertex(tgtId);

                Vec3 offset(0.0);
                double sumWgt = 0.0;
                for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                    const double weight = it->cotWeight();
                    offset += weight * offsets[it->dst()->index()];
                    sumWgt += weight;
                }
                offsets[tgtId] = offset / sumWgt;
            }
        }

        // Update
        const double invQ = 1.0 / 10.0;
        for (int i = 0; i < nTgt; i++) {
            const int tgtId = tgtIds[i];
            Vertex *v = mesh.vertex(tgtId);

            const Vec3 orgPos = v->pos();
            const Vec3 newPos = orgPos + invQ * (newPoss[tgtId] + offsets[tgtId] - orgPos);
            v->setPos(newPos);
            std::cout << "org: " << orgPos << std::endl;
            std::cout << "new: " << newPos << std::endl;
        }
    }
}

}  // namespace tinymesh
