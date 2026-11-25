#include "filters.h"

#include <unordered_set>

#include "core/debug.h"
#include "core/openmp.h"
#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"
#include "core/constants.h"

#include "core/eigen.h"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>

namespace {

inline double gaussKernel(double x, double sigma) {
    using tinymesh::kTwoPi;

    const double sigma2 = sigma * sigma;
    const double coef = 1.0 / (std::sqrt(kTwoPi) * sigma);
    if (x < 2.0 * sigma) {
        return coef * std::exp(-0.5 * x * x / sigma2);
    }

    if (x < 4.0 * sigma) {
        const double a = (4.0 - std::abs(x) / sigma);
        return coef * (1.0 / (16.0 * std::exp(2.0))) * (a * a * a * a);
    }

    return 0.0;
}

void validateFaces(const tinymesh::Mesh &mesh) {
    const int nf = static_cast<int>(mesh.numFaces());
    for (int i = 0; i < nf; ++i) {
        const auto *f = mesh.face(i);
        if (!f->isHole() && f->numCorners() != 3) {
            tinymesh::Error("Non-triangular face detected!");
        }
    }
}

}  // anonymous namespace

namespace tinymesh {

void denoiseNormalGaussian(Mesh &mesh, double sigma, int iterations) {
    // Check triangularity
    validateFaces(mesh);

    // Average edge length
    double avgEdge = mesh.getMeanEdgeLength();

    for (int it = 0; it < iterations; it++) {
        // Smooth vertex positions
        smoothTaubin(mesh);

        // Compute normals, centroids, and areas
        const int nf = (int)mesh.numFaces();
        std::vector<Vec3> normals(nf);
        std::vector<Vec3> centroids(nf);
        std::vector<double> areas(nf);
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            if (f->isHole()) {
                continue;
            }

            std::vector<Vec3> vs;
            for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
                vs.push_back(vit->pos());
            }

            const Vec3 outer = cross(vs[1] - vs[0], vs[2] - vs[0]);
            normals[i] = normalize(outer);
            centroids[i] = (vs[0] + vs[1] + vs[2]) / 3.0;
            areas[i] = 0.5 * length(outer);
        }

        // Filter
        std::vector<Vec3> newNormals(nf);
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            if (f->isHole()) {
                continue;
            }

            const double selfWeight = gaussKernel(0.0, sigma * avgEdge) * areas[i];
            Vec3 norm = selfWeight * normals[i];
            for (auto fit = f->f_begin(); fit != f->f_end(); ++fit) {
                if (fit->isHole()) {
                    continue;
                }

                const int j = fit->index();
                const double d = length(centroids[i] - centroids[j]);
                const double weight = gaussKernel(d, sigma * avgEdge) * areas[j];
                norm += normals[j] * weight;
            }

            const double l = length(norm);
            if (l != 0.0) {
                newNormals[i] = norm / l;
            } else {
                newNormals[i] = normals[i];
            }
        }

        // Update vertex positions
        const int nv = (int)mesh.numVertices();
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nv; i++) {
            Vertex *v = mesh.vertex(i);
            if (v->isBoundary() || v->isLocked()) {
                continue;
            }

            double sumArea = 0.0;
            Vec3 incr(0.0);
            for (auto fit = v->f_begin(); fit != v->f_end(); ++fit) {
                if (fit->isHole()) {
                    continue;
                }

                const int j = fit->index();
                sumArea += areas[j];
                incr += (areas[j] * dot(newNormals[j], centroids[j] - v->pos())) * newNormals[j];
            }
            incr /= sumArea;
            v->setPos(v->pos() + incr);
        }
    }
}

void denoiseNormalBilateral(Mesh &mesh, double sigmaCenter, double sigmaNormal, int iterations) {
    // Check triangularity
    validateFaces(mesh);

    // Average edge length
    double avgEdge = mesh.getMeanEdgeLength();

    for (int it = 0; it < iterations; it++) {
        // Smooth vertex positions
        smoothTaubin(mesh);

        // Compute normals, centroids, and areas
        const int nf = (int)mesh.numFaces();
        std::vector<Vec3> normals(nf);
        std::vector<Vec3> centroids(nf);
        std::vector<double> areas(nf);
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            if (f->isHole()) {
                continue;
            }

            std::vector<Vec3> vs;
            for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
                vs.push_back(vit->pos());
            }

            const Vec3 outer = cross(vs[1] - vs[0], vs[2] - vs[0]);
            normals[i] = normalize(outer);
            centroids[i] = (vs[0] + vs[1] + vs[2]) / 3.0;
            areas[i] = 0.5 * length(outer);
        }

        // Filter
        std::vector<Vec3> newNormals(nf);
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            if (f->isHole()) {
                continue;
            }

            const double selfWeight =
                gaussKernel(0.0, sigmaCenter * avgEdge) * gaussKernel(0.0, sigmaNormal) * areas[i];
            Vec3 norm = selfWeight * normals[i];
            for (auto fit = f->f_begin(); fit != f->f_end(); ++fit) {
                if (fit->isHole()) {
                    continue;
                }

                const int j = fit->index();
                const double d = length(centroids[i] - centroids[j]);
                const double eta = areas[j];
                const double Wc = gaussKernel(d, sigmaCenter * avgEdge);
                const double nd = length(normals[i] - normals[j]);
                const double Ws = gaussKernel(nd, sigmaNormal);
                const double weight = eta * Wc * Ws;
                norm += normals[j] * weight;
            }

            const double l = length(norm);
            if (l != 0.0) {
                newNormals[i] = norm / l;
            } else {
                newNormals[i] = normals[i];
            }
        }

        // Update vertex positions
        const int nv = (int)mesh.numVertices();
        TINYMESH_OMP_PARALLEL_FOR(int i = 0; i < nv; i++) {
            Vertex *v = mesh.vertex(i);
            if (v->isBoundary() || v->isLocked()) {
                continue;
            }

            double sumArea = 0.0;
            Vec3 incr(0.0);
            for (auto fit = v->f_begin(); fit != v->f_end(); ++fit) {
                if (fit->isHole()) {
                    continue;
                }

                const int j = fit->index();
                sumArea += areas[j];
                incr += (areas[j] * dot(newNormals[j], centroids[j] - v->pos())) * newNormals[j];
            }
            incr /= sumArea;
            v->setPos(v->pos() + incr);
        }
    }
}

void denoiseL0Smooth(Mesh &mesh, double alpha, double beta) {
    const int nv = static_cast<int>(mesh.numVertices());
    std::vector<uint8_t> fixed(nv, 0);
    for (int i = 0; i < nv; ++i) {
        const Vertex *v = mesh.vertex(i);
        fixed[i] = v->isBoundary() || v->isLocked();
    }

    const double avgEdge = mesh.getMeanEdgeLength();
    const double avgDihed = mesh.getMeanDihedralAngle();

    // Collect one halfedge per interior edge
    std::vector<Halfedge *> edges;
    for (size_t i = 0; i < mesh.numHalfedges(); ++i) {
        Halfedge *he = mesh.halfedge(i);
        Halfedge *rev = he->rev();

        if (he->face()->isHole() || rev->face()->isHole()) {
            continue;
        }

        if (he->index() > rev->index()) {
            continue;
        }

        edges.push_back(he);
    }
    const int ne = static_cast<int>(edges.size());

    std::vector<EigenTriplet> tripR;
    for (size_t e = 0; e < edges.size(); e++) {
        Halfedge *he = edges[e];
        Halfedge *rev = he->rev();

        Vertex *vh1 = he->src();
        Vertex *vh2 = he->next()->dst();
        Vertex *vh3 = rev->src();
        Vertex *vh4 = rev->next()->dst();

        const int i1 = vh1->index();
        const int i2 = vh2->index();
        const int i3 = vh3->index();
        const int i4 = vh4->index();
        tripR.emplace_back(e, i1, 1.0);
        tripR.emplace_back(e, i2, -1.0);
        tripR.emplace_back(e, i3, 1.0);
        tripR.emplace_back(e, i4, -1.0);
    }

    EigenSparseMatrix R(ne, nv);
    R.setFromTriplets(tripR.begin(), tripR.end());

    std::vector<EigenTriplet> tripI;
    for (int i = 0; i < nv; i++) {
        tripI.emplace_back(i, i, 1.0);
    }
    EigenSparseMatrix I(nv, nv);
    I.setFromTriplets(tripI.begin(), tripI.end());

    // Store vertex positions to Eigen Matrix
    EigenMatrix pVec, pVecInit;
    pVec.resize(nv, 3);
    pVecInit.resize(nv, 3);
    for (int i = 0; i < nv; i++) {
        const Vec3 p = mesh.vertex(i)->pos();
        pVec(i, 0) = p[0];
        pVec(i, 1) = p[1];
        pVec(i, 2) = p[2];
        pVecInit(i, 0) = p[0];
        pVecInit(i, 1) = p[1];
        pVecInit(i, 2) = p[2];
    }

    // Solve
    const double bmax = 1.0e3 * avgEdge;
    const double lambda = 0.02 * avgEdge * avgEdge * avgDihed;
    const double mu = std::sqrt(2.0);
    alpha = alpha * avgEdge;
    beta = beta * avgEdge;

    // Smooth vertex positions
    smoothTaubin(mesh);

    while (beta < bmax) {
        // Construct sparse matrix D
        std::vector<EigenTriplet> tripD;
        for (int e = 0; e < edges.size(); e++) {
            Halfedge *he = edges[e];
            Halfedge *rev = he->rev();

            Vertex *vh1 = he->src();
            Vertex *vh2 = he->next()->dst();
            Vertex *vh3 = rev->src();
            Vertex *vh4 = rev->next()->dst();

            const int i1 = vh1->index();
            const int i2 = vh2->index();
            const int i3 = vh3->index();
            const int i4 = vh4->index();

            const Vec3 p1 = vh1->pos();
            const Vec3 p2 = vh2->pos();
            const Vec3 p3 = vh3->pos();
            const Vec3 p4 = vh4->pos();

            const double S123 = 0.5 * length(cross(p2 - p1, p2 - p3));
            const double S134 = 0.5 * length(cross(p4 - p1, p4 - p3));
            double l13 = length(p1 - p3);

            const double coef1 =
                (S123 * dot(p4 - p3, p3 - p1) + S134 * dot(p1 - p3, p3 - p2)) / (l13 * l13 * (S123 + S134));
            const double coef2 = S134 / (S123 + S134);
            const double coef3 =
                (S123 * dot(p3 - p1, p1 - p4) + S134 * dot(p2 - p1, p1 - p3)) / (l13 * l13 * (S123 + S134));
            const double coef4 = S123 / (S123 + S134);

            tripD.emplace_back(e, i1, coef1);
            tripD.emplace_back(e, i2, coef2);
            tripD.emplace_back(e, i3, coef3);
            tripD.emplace_back(e, i4, coef4);
        }

        EigenSparseMatrix D(ne, nv);
        D.setFromTriplets(tripD.begin(), tripD.end());

        // Solve sub-problem for delta
        EigenMatrix y = D * pVec;
        EigenMatrix yNorm2 = y.array().square().rowwise().sum().replicate(1, 3);
        EigenMatrix yNorm = yNorm2.array().sqrt().matrix();
        EigenMatrix delta = (yNorm2.array() < (lambda / beta)).select(0.0, y);

        // Solve sub-problem for positions
        EigenSparseMatrix A = I + alpha * R.transpose() * R + beta * D.transpose() * D;
        EigenMatrix bb = pVecInit + beta * D.transpose() * delta;

        // Eliminate fixed vertices from the system
        for (int k = 0; k < A.outerSize(); ++k) {
            for (EigenSparseMatrix::InnerIterator it(A, k); it; ++it) {
                const int row = it.row();
                const int col = it.col();
                if (!fixed[row] && fixed[col]) {
                    bb.row(row) -= it.value() * pVecInit.row(col);
                }

                if (fixed[row] || fixed[col]) {
                    it.valueRef() = 0.0;
                }
            }
        }

        for (int i = 0; i < nv; ++i) {
            if (fixed[i]) {
                A.coeffRef(i, i) = 1.0;
                bb.row(i) = pVecInit.row(i);
            }
        }
        A.prune(0.0);

        Eigen::SimplicialLDLT<EigenSparseMatrix> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            Error("L0 smoothing factorization failed!");
        }

        pVec = solver.solve(bb);
        if (solver.info() != Eigen::Success) {
            Error("L0 smoothing solve failed");
        }

        // Update
        beta *= mu;
        alpha *= 0.5;
    }

    // Update vertex in a mesh
    for (int i = 0; i < nv; i++) {
        if (fixed[i]) {
            continue;
        }

        const Vec3 newPos(pVec(i, 0), pVec(i, 1), pVec(i, 2));
        mesh.vertex(i)->setPos(newPos);
    }
}

}  // namespace tinymesh
