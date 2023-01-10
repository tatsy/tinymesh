#define TINYMESH_API_EXPORT
#include "restore.h"

#include <iostream>
#include <queue>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include <Eigen/LU>
#include <Eigen/SparseCholesky>

#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"
#include "core/bvh.h"
#include "core/eigen.h"
#include "core/openmp.h"
#include "core/progress.h"

#include "ops/ops.h"

namespace tinymesh {

namespace {

double solveRigidICP(const std::vector<Vertex *> &tgtPatch, const std::vector<Vertex *> &srcPatch,
                     const std::unordered_map<Vertex *, EigenVector> &hks, EigenMatrix3 &R, EigenVector3 &t) {
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
            Vertex *src = srcPatch[i];
            Vertex *tgt = tgtPatch[j];
            const auto itSrc = hks.find(src);
            const auto itTgt = hks.find(tgt);
            const double dist = (itSrc->second - itTgt->second).squaredNorm();
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

}  // namespace

void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound) {
    mesh.holeFillMinDihedral_(face, dihedralBound);
}

void holeFillAdvancingFront(Mesh &mesh, Face *face) {
    mesh.holeFillAdvancingFront_(face);
}

void holeFillContextCoherent(Mesh &mesh, int maxiters) {
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
    const double avgLen = mesh.getMeanEdgeLength();
    const double avgArea = mesh.getMeanFaceArea();
    for (size_t i = 0; i < mesh.numHalfedges(); i++) {
        Halfedge *he = mesh.halfedge(i);
        he->setIsBoundary(false);
    }

    // 2. Compute patch radius
    auto curvatures = getPrincipalCurvatures(mesh);
    const EigenVector &kv_max = std::get<0>(curvatures);
    const EigenVector &kv_min = std::get<1>(curvatures);
    std::vector<double> maxCurves(mesh.numVertices());
    double maxK = 0.0;
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        maxCurves[i] = std::max(std::abs(kv_max(i)), std::abs(kv_min(i)));
        maxK = std::max(maxK, maxCurves[i]);
    }

    std::vector<double> pis = { 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0 };
    std::vector<double> mus = { 0.0, 0.5 * maxK, maxK };
    std::vector<double> vars = { 1.0, 1.0, 1.0 };

    auto gauss = [](double x, double v) -> double { return std::exp(-0.5 * x * x / v) / std::sqrt(2.0 * Pi * v); };

    for (int kIter = 0; kIter < 20; kIter++) {
        // E step
        std::vector<std::vector<double>> gamma(mesh.numVertices(), std::vector<double>(3, 0.0));
        for (int i = 0; i < mesh.numVertices(); i++) {
            double sumWgt = 0.0;
            for (int k = 0; k < 3; k++) {
                gamma[i][k] = pis[k] * gauss(mus[k] - maxCurves[i], vars[k]);
                sumWgt += gamma[i][k];
            }

            for (int k = 0; k < 3; k++) {
                gamma[i][k] /= sumWgt;
            }
        }

        // M step
        pis.assign(3, 0.0);
        for (int i = 0; i < mesh.numVertices(); i++) {
            for (int k = 0; k < 3; k++) {
                pis[k] += gamma[i][k] / (double)mesh.numVertices();
            }
        }

        mus.assign(3, 0.0);
        for (int i = 0; i < mesh.numVertices(); i++) {
            for (int k = 0; k < 3; k++) {
                const double Npi = mesh.numVertices() * pis[k];
                mus[k] += gamma[i][k] * maxCurves[i] / Npi;
            }
        }

        vars.assign(3, 0.0);
        for (int i = 0; i < mesh.numVertices(); i++) {
            for (int k = 0; k < 3; k++) {
                const double Npi = mesh.numVertices() * pis[k];
                const double diff = maxCurves[i] - mus[k];
                vars[k] += gamma[i][k] * diff * diff / Npi;
            }
        }
    }
    printf("GMM[0]: pi=%.3f, mu=%.3f, sigma=%.3f\n", pis[0], mus[0], std::sqrt(vars[0]));
    printf("GMM[1]: pi=%.3f, mu=%.3f, sigma=%.3f\n", pis[1], mus[1], std::sqrt(vars[1]));
    printf("GMM[2]: pi=%.3f, mu=%.3f, sigma=%.3f\n", pis[2], mus[2], std::sqrt(vars[2]));

    const double patchRadiusReal = (1.0 / (avgLen * mus[1]) + 1.0) * avgLen;
    printf("R = %f (avgLen = %f)\n", patchRadiusReal, avgLen);
    const int patchRadius = (int)std::ceil(patchRadiusReal / avgLen);

    // 2. Compute heat kernel signatures
    int hksDims = -1;
    std::unordered_map<Vertex *, EigenVector> hksTable;
    {
        auto start = std::chrono::system_clock::now();

        const EigenSparseMatrix L = getMeshLaplacian(mesh, MeshLaplace::Cotangent);
        const EigenMatrix hks_ = getHeatKernelSignatures(L, 300, 100);
        hksDims = (int)hks_.cols();
        for (size_t i = 0; i < mesh.numVertices(); i++) {
            Vertex *v = mesh.vertex(i);
            hksTable[v] = hks_.row(i);
        }

        auto end = std::chrono::system_clock::now();
        const int64_t elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        printf("[ HKS ] %f sec\n", elapsed / 1000.0);
    }

    // Iteration
    ProgressBar pbar(maxiters);
    for (int kIter = 0; kIter < maxiters; kIter++) {
        const int N = (int)mesh.numVertices();

        // 3. Collect patches
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

        // Inverse table to source/target IDs
        std::unordered_map<int, int> src2glbIds;
        for (int i = 0; i < nSrc; i++) {
            src2glbIds[srcIds[i]] = i;
        }
        std::unordered_map<int, int> tgt2glbIds;
        for (int i = 0; i < nTgt; i++) {
            tgt2glbIds[tgtIds[i]] = i;
        }

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

        // 4. Compute patch descriptors
        std::vector<EigenVector> descriptors(N);
        omp_parallel_for(int i = 0; i < N; i++) {
            EigenVector hksMean(hksDims);
            EigenVector hksSigma(hksDims);

            hksMean.setZero();
            hksSigma.setZero();
            for (auto v : patchVerts[i]) {
                const EigenVector f = hksTable[v];
                hksMean += f;
                hksSigma += f.cwiseProduct(f);
            }

            hksMean /= (double)N;
            hksSigma /= (double)N;
            hksSigma = (hksSigma - hksMean.cwiseProduct(hksMean)).array().max(0.0);
            hksSigma = hksSigma.cwiseSqrt();

            hksMean /= hksMean.maxCoeff();
            hksSigma /= hksSigma.maxCoeff();

            descriptors[i].resize(hksMean.rows() + hksSigma.rows());
            descriptors[i] << hksMean, hksSigma;
        }

        // 5. Collect candidate patches
        const int nCands = std::max(10, (int)(0.001 * nSrc));
        std::vector<std::vector<int>> candidates(nTgt, std::vector<int>(nCands, 0));
        for (int i = 0; i < nTgt; i++) {
            std::vector<std::pair<double, int>> dists;
            for (int j = 0; j < nSrc; j++) {
                const int tgtId = tgtIds[i];
                const int srcId = srcIds[j];
                const double d = (descriptors[tgtId] - descriptors[srcId]).squaredNorm();
                dists.emplace_back(d, srcId);
            }

            std::sort(dists.begin(), dists.end());
            for (int j = 0; j < nCands; j++) {
                candidates[i][j] = dists[j].second;
            }
        }

        const auto flff = getFeatureLineFieldWithFlags(mesh, true);
        const EigenMatrix &flf = std::get<0>(flff);
        const EigenVector &rvFlags = std::get<1>(flff);

        // 6. Compute most similar patch for each target patch
        std::vector<EigenMatrix3> rotations(nTgt);
        std::vector<EigenVector3> translations(nTgt);
        std::vector<int> pairIds(nTgt, 0);
        double newError = 0.0;
        for (int i = 0; i < nTgt; i++) {
            const int tgtId = tgtIds[i];
            EigenMatrix3 R;
            EigenVector3 t;
            double minEps = 1.0e20;
            int minId = -1;
            for (int srcId : candidates[i]) {
                const std::vector<Vertex *> &tgtPatch = patchVerts[tgtId];
                const std::vector<Vertex *> &srcPatch = patchVerts[srcId];
                const double eps = solveRigidICP(tgtPatch, srcPatch, hksTable, R, t);
                if (eps < minEps) {
                    minEps = eps;
                    minId = srcId;
                    rotations[i] = R;
                    translations[i] = t;
                }
            }
            pairIds[i] = minId;
            newError += minEps;
        }
        newError /= nTgt;
        pbar.set_description("error=%7.5f", newError);

        // Compute patch dissimilarities
        std::unordered_map<int, double> dissimilarities;
        for (int i = 0; i < nTgt; i++) {
            double dissim = 0.0;
            double sumWgt = 0.0;
            const std::vector<Vertex *> &tgtPatch = patchVerts[tgtIds[i]];
            const std::vector<Vertex *> &srcPatch = patchVerts[pairIds[i]];
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
                const double l = length(ft - (Vec3)(rotations[i] * fs));
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

        // 7. Update vertex positions
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

        // Smooth offsets by Laplacian
        std::vector<EigenTriplet> tripTgt;
        std::vector<EigenTriplet> tripSrc;
        std::unordered_map<uint32_t, uint32_t> uniqueBoundary;
        int countBoundary = 0;

        for (int i = 0; i < nTgt; i++) {
            const int tgtId = tgtIds[i];
            const Vertex *v = mesh.vertex(tgtId);

            double sumWgt = 0.0;
            for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                Vertex *u = it->dst();
                const double weight = it->cotWeight();
                if (u->isLocked()) {
                    // source
                    if (uniqueBoundary.count(u->index()) == 0) {
                        uniqueBoundary[u->index()] = countBoundary;
                        countBoundary++;
                    }
                    const int j = uniqueBoundary[u->index()];
                    tripSrc.emplace_back(i, j, -weight);
                } else {
                    // target
                    const int j = tgt2glbIds[u->index()];
                    tripTgt.emplace_back(i, j, -weight);
                }
                sumWgt += weight;
            }
            tripTgt.emplace_back(i, i, sumWgt);
        }

        EigenSparseMatrix LL(nTgt, nTgt);
        LL.setFromTriplets(tripTgt.begin(), tripTgt.end());
        EigenSparseMatrix SS(nTgt, countBoundary);
        SS.setFromTriplets(tripSrc.begin(), tripSrc.end());

        EigenMatrix xxSrc(countBoundary, 3);
        xxSrc.setZero();
        for (auto it : uniqueBoundary) {
            const Vec3 &v = offsets[it.first];
            xxSrc.row(it.second) << v.x(), v.y(), v.z();
        }

        EigenMatrix bb = -SS * xxSrc;
        Eigen::SimplicialLDLT<EigenSparseMatrix> solver(LL);
        EigenMatrix xx = solver.solve(bb);

        for (int i = 0; i < nTgt; i++) {
            offsets[tgtIds[i]] = Vec3(xx(i, 0), xx(i, 1), xx(i, 2));
        }

        // Update
        const double invQ = 1.0 / 10.0;
        for (int i = 0; i < nTgt; i++) {
            const int tgtId = tgtIds[i];
            Vertex *v = mesh.vertex(tgtId);

            const Vec3 orgPos = v->pos();
            const Vec3 newPos = orgPos + invQ * (newPoss[tgtId] + offsets[tgtId] - orgPos);
            v->setPos(newPos);
        }

        // 8. Repair mesh
        bool isRepaired = false;

        // Check edge split
        for (int tgtId : tgtIds) {
            bool update = true;
            do {
                update = false;
                if (tgtId >= mesh.numVertices()) {
                    break;
                }

                Vertex *v = mesh.vertex(tgtId);
                for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                    if (it->length() > 1.5 * avgLen) {
                        if (mesh.splitHE(it.ptr())) {
                            isRepaired = true;
                            update = true;
                            break;
                        }
                    }
                }
            } while (update);
        }

        // Check edge flip
        for (int tgtId : tgtIds) {
            bool update = true;
            do {
                update = false;
                if (tgtId >= mesh.numVertices()) {
                    break;
                }

                Vertex *v = mesh.vertex(tgtId);
                for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                    const Vec3 v0 = it->src()->pos();
                    const Vec3 v1 = it->dst()->pos();
                    const Vec3 vl = it->next()->dst()->pos();
                    const Vec3 vr = it->rev()->next()->dst()->pos();
                    const double al = angle(v0, vl, v1);
                    const double ar = angle(v0, vr, v1);
                    if (al + ar > Pi) {
                        if (mesh.flipHE(it.ptr())) {
                            isRepaired = true;
                            update = true;
                            break;
                        }
                    }
                }
            } while (update);
        }

        // Check edge collapse
        for (int tgtId : tgtIds) {
            bool update = true;
            do {
                update = false;
                if (tgtId >= mesh.numVertices()) {
                    break;
                }

                Vertex *v = mesh.vertex(tgtId);
                for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                    const Vec3 v0 = it->src()->pos();
                    const Vec3 v1 = it->dst()->pos();
                    const Vec3 vl = it->next()->dst()->pos();
                    const Vec3 vr = it->rev()->next()->dst()->pos();
                    const double al = angle(v0, vl, v1);
                    const double ar = angle(v0, vr, v1);
                    if (al < Pi / 10.0 || ar < Pi / 10.0) {
                        if (mesh.collapseHE(it.ptr())) {
                            isRepaired = true;
                            update = true;
                            break;
                        }
                    }
                }
            } while (update);
        }

        // Check face collapse
        for (int tgtId : tgtIds) {
            bool update = true;
            do {
                update = false;
                if (tgtId >= mesh.numVertices()) {
                    break;
                }

                Vertex *v = mesh.vertex(tgtId);
                for (auto it = v->f_begin(); it != v->f_end(); ++it) {
                    if (it->index() < 0) continue;

                    if (it->area() < 0.01 * avgArea) {
                        if (mesh.collapseFace(it.ptr())) {
                            isRepaired = true;
                            update = true;
                            break;
                        }
                    }
                }
            } while (update);
        }

        if (isRepaired) {
            // Update HKS table
            for (int i = 0; i < mesh.numVertices(); i++) {
                Vertex *v = mesh.vertex(i);
                if (hksTable.count(v) == 0) {
                    EigenVector hksv = EigenVector::Zero(hksDims);
                    int count = 0;
                    for (auto it = v->v_begin(); it != v->v_end(); ++it) {
                        if (hksTable.count(it.ptr()) != 0) {
                            hksv += hksTable[it.ptr()];
                            count += 1;
                        }
                    }
                    hksv /= (double)count;
                    hksTable[v] = hksv;
                }
            }
        }

        // Verification
        mesh.verify();
        pbar.step();
    }
}

}  // namespace tinymesh
