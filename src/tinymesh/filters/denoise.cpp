#define TINYMESH_API_EXPORT
#include "filters.h"

#include "core/debug.h"
#include "core/openmp.h"
#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"

namespace {

inline double K(double x, double sigma) {
    const double sigma2 = sigma * sigma;
    const double coef = 1.0 / (std::sqrt(2.0 * Pi) * sigma);
    if (x < 2.0 * sigma) {
        return coef * std::exp(-x * x / sigma2);
    }

    if (x < 4.0 * sigma) {
        const double a = (4.0 - std::abs(x) / sigma);
        return coef * (1.0 / (16.0 * std::exp(2.0))) * (a * a * a * a); 
    }

    return 0.0;
}

}  // anonymous namespace

namespace tinymesh {

void denoiseNormalGaussian(Mesh &mesh, double sigma, int iterations) {
    // Average edge length
    double avgEdge = 0.0;
    for (int i = 0; i < mesh.num_halfedges(); i++) {
        avgEdge += mesh.halfedge(i)->length();
    }
    avgEdge /= mesh.num_halfedges();

    for (int it = 0; it < iterations; it++) {
        // Compute normals, centroids, and areas
        const int nf = mesh.num_faces();
        std::vector<Vec3> normals(nf);
        std::vector<Vec3> centroids(nf);
        std::vector<double> areas(nf);
        omp_parallel_for(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            std::vector<Vec3> vs;
            for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
                vs.push_back(vit->pos());
            }

            if (vs.size() != 3) {
                FatalError("Mesh is not triangular! Call \"mesh.triangulate()\" first!");
            }

            const Vec3 outer = cross(vs[1] - vs[0], vs[2] - vs[0]);
            normals[i] = normalize(outer);
            centroids[i] = (vs[0] + vs[1] + vs[2]) / 3.0;
            areas[i] = 0.5 * length(outer);
        }

        // Filter
        std::vector<Vec3> newNormals(nf);
        omp_parallel_for(int i = 0; i < nf; i++) {
            Face *f = mesh.face(i);
            Vec3 norm(0.0);
            for (auto fit = f->f_begin(); fit != f->f_end(); ++fit) {
                const int j = fit->index();
                const double d = length(centroids[i] - centroids[j]);
                const double weight = K(d, sigma * avgEdge) * areas[j];
                norm += normals[j] * weight;
            }

            const double l = length(norm);
            if (l != 0.0) {
                newNormals[i] = norm / l;
            }
        }

        // Update vertex positions
        const int nv = mesh.num_vertices();
        omp_parallel_for(int i = 0; i < nv; i++) {
            Vertex *v = mesh.vertex(i);

            double sumArea = 0.0;
            Vec3 incr(0.0);
            for (auto fit = v->f_begin(); fit != v->f_end(); ++fit) {
                const int j = fit->index();
                sumArea += areas[j];
                incr += (areas[j] * dot(newNormals[j], centroids[j] - v->pos())) * newNormals[j];
            }
            incr /= sumArea;
            v->setPos(v->pos() + incr);
        }

        // Smooth vertex positions
        smoothTaubin(mesh);
    }
}

}  // namespace tinymesh