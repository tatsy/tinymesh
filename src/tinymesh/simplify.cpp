#define TINYMESH_API_EXPORT
#include "simplify.h"

#include "mesh.h"
#include "vector.h"
#include "vertex.h"
#include "halfedge.h"
#include "smooth.h"

namespace tinymesh {

void simplify(Mesh &mesh, int maxiter) {
    printf("*** Original ***\n");
    printf("#vert: %d\n", (int)mesh.num_vertices());
    printf("#face: %d\n", (int)mesh.num_faces());

    std::vector<Halfedge*> hes;

    for (int k = 0; k < maxiter; k++) {
        // Compute average edge length
        double L = 0.0;
        int count = 0;
        for (auto it = mesh.he_begin(); it != mesh.he_end(); ++it) {
            L += it->length();
            count += 1;
        }
        L /= count;
        printf("Avg edge length: %f\n", L);

        // Split long edges
        hes.clear();
        for (auto it = mesh.he_begin(); it != mesh.he_end(); ++it) {
            hes.push_back(it.ptr());
        }

        for (Halfedge *he : hes) {
            if (!he || he->index() >= mesh.num_halfedges()) {
                continue;
            }

            if (he->length() > 1.333 * L) {
                mesh.splitHE(he);
            }
        }

        printf("*** After split ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Collapse short edges
        hes.clear();
        for (auto it = mesh.he_begin(); it != mesh.he_end(); ++it) {
            hes.push_back(it.ptr());
        }

        for (Halfedge *he : hes) {
            if (!he || he->index() >= mesh.num_halfedges()) {
                continue;
            }

            if (he->length() < 0.65 * L) {
                if (mesh.collapseHE(he)) {
                }
            }
        }

        printf("*** After collapse ***\n");
        printf("#vert: %d\n", (int)mesh.num_vertices());
        printf("#face: %d\n", (int)mesh.num_faces());

        // Flip edges
        for (auto it = mesh.he_begin(); it != mesh.he_end(); ++it) {
            Vertex *v0 = it->src();
            Vertex *v1 = it->dst();
            Vertex *v2 = it->next()->dst();
            Vertex *v3 = it->rev()->next()->dst();
            const int d0 = v0->degree();
            const int d1 = v1->degree();
            const int d2 = v2->degree();
            const int d3 = v3->degree();

            const int score = std::abs(d0 - 6) + std::abs(d1 - 6) + std::abs(d2 - 6) + std::abs(d3 - 6);
            const int after = std::abs(d0 - 1 - 6) + std::abs(d1 - 1 - 6) + std::abs(d2 + 1 - 6) + std::abs(d3 + 1 - 6);
            if (score > after) {
                mesh.flipHE(it.ptr());
            }
        }

        // Volonoi tessellation
        for (int l = 0; l < 10; l++) {
            int index;
            const int nv = mesh.num_vertices();
            std::vector<Vector> centroids(nv);
            std::vector<Vector> normals(nv);

            // Compute centroids and tangent planes
            index = 0;
            for (auto it = mesh.v_begin(); it != mesh.v_end(); ++it) {
                // Collect surrounding vertices
                Vector org = it->pt();
                std::vector<Vector> pts;
                for (auto vit = it->v_begin(); vit != it->v_end(); ++vit) {
                    pts.push_back(vit->pt());
                }

                // Compute centroids, tangents, and binormals
                Vector cent(0.0);
                Vector norm(0.0);
                for (int i = 0; i < pts.size(); i++) {
                    const int j = (i + 1) % pts.size();
                    Vector e1 = pts[i] - org;
                    Vector e2 = pts[j] - org;
                    Vector g = (org + pts[i] + pts[j]) / 3.0;

                    cent += g;
                    norm += e1.cross(e2);
                }
                cent /= pts.size();
                norm.normalize();

                centroids[index] = cent;
                normals[index] = norm;
                index += 1;
            }

            // Update vertex positions
            index = 0;
            for (auto it = mesh.v_begin(); it != mesh.v_end(); ++it) {
                const Vector pt = it->pt();
                Vector e = centroids[index] - pt;
                e -= normals[index] * e.dot(normals[index]);
                it->setPt(pt + e);
                index += 1;
            }

        }
        // Laplacian smoothing
        //smooth(mesh, 0.1);
    }

    mesh.verify();
}

}  // namespace tinymesh