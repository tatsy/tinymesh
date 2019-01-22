#define TINYMESH_API_EXPORT
#include "simplify.h"

#include "mesh.h"
#include "vector.h"
#include "vertex.h"

namespace tinymesh {

void simplify(Mesh &mesh, int maxiter) {
    int index;
    const int nv = mesh.num_vertices();
    std::vector<Vector> centroids(nv);
    std::vector<Vector> normals(nv);

    for (int k = 0; k < maxiter; k++) {
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

                cent += pts[i];
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
}

}  // namespace tinymesh