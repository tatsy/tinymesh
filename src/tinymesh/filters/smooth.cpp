#define TINYMESH_API_EXPORT
#include "filters.h"

#include <algorithm>
#include <unordered_map>

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

using ScalarType = double;
using IndexType = int64_t;
using Triplet = Eigen::Triplet<ScalarType, IndexType>;
using SparseMatrix = Eigen::SparseMatrix<ScalarType>;

#include "core/debug.h"
#include "core/openmp.h"
#include "core/face.h"
#include "core/vertex.h"
#include "core/halfedge.h"

namespace tinymesh {

void smoothLaplacian(Mesh &mesh, double epsilon, bool cotangentWeight, int iterations) {
    for (int it = 0; it < iterations; it++) {
        // Volonoi tessellation
        const int nv = (int)mesh.numVertices();
        std::vector<Vec3> centroids(nv);

        // Compute centroids and tangent planes
        omp_parallel_for(int i = 0; i < nv; i++) {
            Vertex *v = mesh.vertex(i);
            std::vector<Vec3> pts;
            Vec3 cent(0.0);
            double sumWgt = 0.0;
            for (auto it = v->ohe_begin(); it != v->ohe_end(); ++it) {
                double weight = 1.0;
                if (cotangentWeight) {
                    weight = it->cotWeight();
                }
                cent += weight * it->dst()->pos();
                sumWgt += weight;
            }
            centroids[i] = cent / sumWgt;
        }

        // Update vertex positions
        omp_parallel_for(int i = 0; i < (int)mesh.numVertices(); i++) {
            Vertex *v = mesh.vertex(i);
            if (v->isBoundary() || v->isLocked()) {
                continue;
            }

            const Vec3 pt = v->pos();
            const Vec3 e = centroids[i] - pt;
            const Vec3 newPos = (1.0 - epsilon) * pt + epsilon * (pt + e);
            v->setPos(newPos);
        }
    }
}

void smoothTaubin(Mesh &mesh, double shrink, double inflate, int iterations) {
    for (int it = 0; it < iterations * 2; it++) {
        // Volonoi tessellation
        const int nv = (int)mesh.numVertices();
        std::vector<Vec3> centroids(nv);

        // Compute centroids and tangent planes
        omp_parallel_for(int i = 0; i < nv; i++) {
            Vertex *v = mesh.vertex(i);

            // Collect surrounding vertices
            std::vector<Vec3> pts;
            for (auto vit = v->v_begin(); vit != v->v_end(); ++vit) {
                pts.push_back(vit->pos());
            }

            // Compute centroids
            Vec3 cent(0.0);
            for (int i = 0; i < (int)pts.size(); i++) {
                cent += pts[i];
            }
            centroids[i] = cent / (double)pts.size();
        }

        // Update vertex positions
        const double epsilon = it % 2 == 0 ? shrink : -inflate;
        omp_parallel_for(int i = 0; i < (int)mesh.numVertices(); i++) {
            Vertex *v = mesh.vertex(i);
            if (v->isBoundary() || v->isLocked()) {
                continue;
            }

            const Vec3 pt = v->pos();
            const Vec3 e = centroids[i] - pt;
            const Vec3 newPos = (1.0 - epsilon) * v->pos() + epsilon * (pt + e);
            v->setPos(newPos);
        }
    }
}

void implicitFairing(Mesh &mesh, double epsilon, int iterations) {
    /*
     * This method is based on the following paper. The normalized version in Sec.5.5 is used.
     * Desbrun et al. "Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow", 1999.
     */

    for (int it = 0; it < iterations; it++) {
        // Indexing vertices
        const int n_verts = (int)mesh.numVertices();
        std::unordered_map<Vertex *, int64_t> v2i;
        for (int i = 0; i < n_verts; i++) {
            v2i.insert(std::make_pair(mesh.vertex(i), i));
        }

        // Compute diffusion matrix
        Eigen::MatrixXd X(n_verts, 3);
        std::vector<Triplet> triplets;

        for (int i = 0; i < n_verts; i++) {
            Vertex *v = mesh.vertex(i);

            // Compute weights
            double sumW = 0.0;
            std::vector<Triplet> tripletsInColumn;
            for (auto he_it = v->ohe_begin(); he_it != v->ohe_end(); ++he_it) {
                const double W = he_it->cotWeight();
                tripletsInColumn.emplace_back(i, he_it->dst()->index(), W);
                if (std::isnan(W) || std::isinf(W)) {
                    Warn("NaN of inf matrix entry is detedted!");
                }
                sumW += W;
            }

            for (const auto &t : tripletsInColumn) {
                triplets.emplace_back(t.row(), t.col(), t.value() / sumW);
            }
            triplets.emplace_back(i, i, -1.0);

            X(i, 0) = v->pos()[0];
            X(i, 1) = v->pos()[1];
            X(i, 2) = v->pos()[2];
        }

        // Solve sparse linear system
        SparseMatrix K(n_verts, n_verts);
        K.setFromTriplets(triplets.begin(), triplets.end());

        SparseMatrix I(n_verts, n_verts);
        I.setIdentity();

        SparseMatrix A = I - epsilon * K;

        Eigen::BiCGSTAB<SparseMatrix> cg;
        cg.setTolerance(1.0e-6);
        cg.setMaxIterations(50);
        cg.compute(A);

        Eigen::MatrixXd Xnext(n_verts, 3);
        Xnext = cg.solve(X);

        for (int i = 0; i < n_verts; i++) {
            const Vec3 newpos = Vec3(Xnext(i, 0), Xnext(i, 1), Xnext(i, 2));
            mesh.vertex(i)->setPos(newpos);
        }
    }
}

}  // namespace tinymesh
