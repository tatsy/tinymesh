#define TINYMESH_API_EXPORT
#include "remesh.h"

#include <random>
#include <algorithm>

#include "core/vec.h"
#include "trimesh/mesh.h"
#include "trimesh/vertex.h"
#include "trimesh/halfedge.h"
#include "trimesh/face.h"
#include "geoproc/smooth.h"

namespace tinymesh {

void remesh(Mesh &mesh, double ratioLower, double ratioUpper, int maxiter) {
    int count;
    double Lavg, Lvar, Lstd;

    std::vector<uint32_t> indices;
    std::random_device randev;
    std::mt19937 rnd(randev());

    for (int k = 0; k < maxiter; k++) {
        printf("*** Original #%d ***\n", k + 1);
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Compute average edge length
        Lavg = 0.0;
        Lvar = 0.0;

        count = 0;
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            Halfedge *he = mesh.halfedge(i);
            const double l = he->length();
            Lavg += l;
            Lvar += l * l;
            count += 1;
        }

        Lavg = Lavg / count;
        Lvar = Lvar / count - Lavg * Lavg;
        Lstd = std::sqrt(Lvar);

        // Split long edges
        indices.clear();
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < mesh.num_halfedges()) {
                Halfedge *he = mesh.halfedge(i);
                const double l = he->length();
                const double p = (l - Lavg) / Lstd;
                if (l >= Lavg * ratioUpper) {
                    mesh.splitHE(he);
                }
            }
        }

        printf("*** After split ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Collapse short edges
        indices.clear();
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            indices.push_back(i);
        }

        std::shuffle(indices.begin(), indices.end(), rnd);

        for (int i : indices) {
            if (i >= 0 && i < mesh.num_halfedges()) {
                Halfedge *he = mesh.halfedge(i);
                const double l = he->length();
                const double p = (l - Lavg) / Lstd;
                if (l <= Lavg * ratioLower) {
                    mesh.collapseHE(he);
                }
            }
        }

        printf("*** After collapse ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Flip edges
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            Halfedge *he = mesh.halfedge(i);
            if (he->face()->isBoundary() || he->rev()->face()->isBoundary()) {
                continue;
            }

            Vertex *v0 = he->src();
            Vertex *v1 = he->dst();
            Vertex *v2 = he->next()->dst();
            Vertex *v3 = he->rev()->next()->dst();
            const int d0 = v0->degree();
            const int d1 = v1->degree();
            const int d2 = v2->degree();
            const int d3 = v3->degree();

            const int score = std::abs(d0 - 6) + std::abs(d1 - 6) + std::abs(d2 - 6) + std::abs(d3 - 6);
            const int after = std::abs(d0 - 1 - 6) + std::abs(d1 - 1 - 6) + std::abs(d2 + 1 - 6) + std::abs(d3 + 1 - 6);
            if (score > after) {
                mesh.flipHE(he);
            }
        }

        // Volonoi tessellation
        for (int l = 0; l < 5; l++) {
            int index;
            const int nv = mesh.num_vertices();
            std::vector<Vec> centroids(nv);
            std::vector<Vec> normals(nv);

            // Compute centroids and tangent planes
            index = 0;
            for (int i = 0; i < mesh.num_vertices(); i++, index++) {
                Vertex *v = mesh.vertex(i);

                // Collect surrounding vertices
                Vec org = v->pos();
                std::vector<Vec> pts;
                for (auto vit = v->v_begin(); vit != v->v_end(); ++vit) {
                    pts.push_back(vit->pos());
                }

                // Compute centroids, tangents, and binormals
                Vec cent(0.0);
                Vec norm(0.0);
                for (int i = 0; i < pts.size(); i++) {
                    const int j = (i + 1) % pts.size();
                    Vec e1 = pts[i] - org;
                    Vec e2 = pts[j] - org;
                    Vec g = (org + pts[i] + pts[j]) / 3.0;

                    cent += g;
                    norm += cross(e1, e2);
                }

                cent /= pts.size();
                const double l = length(norm);

                if (l != 0.0) {
                    centroids[index] = cent;
                    normals[index] = norm / l;
                }
            }

            // Update vertex positions
            index = 0;
            for (int i = 0; i < mesh.num_vertices(); i++, index++) {
                Vertex *v = mesh.vertex(i);
                if (v->isBoundary()) {
                    continue;
                }

                if (length(normals[index]) != 0.0) {
                    const Vec pt = v->pos();
                    Vec e = centroids[index] - pt;
                    e -= normals[index] * dot(e, normals[index]);
                    v->setPos(pt + e);
                }
            }
        }

        // Laplacian smoothing
        smooth(mesh, 1.0);
    }
}

}  // namespace tinymesh