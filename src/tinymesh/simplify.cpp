#define TINYMESH_API_EXPORT
#include "simplify.h"

#include "mesh.h"

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
                printf("%f %f %f\n", vit->pt().x, vit->pt().y, vit->pt().z);
                pts.push_back(vit->pt());
            }

            // Compute centroids, tangents, and binormals
            Vector cent(0.0);
            Vector norm(0.0);
            double sumA = 0.0;
            for (int i = 0; i < pts.size(); i++) {
                const int j = (i + 1) % pts.size();
                Vector e1 = pts[i] - org;
                Vector e2 = pts[j] - org;
                double A = 0.5 * e1.cross(e2).length();

                cent += A * pts[i];
                norm += e1.cross(e2);
                sumA += A;
            }
            cent /= sumA;
            norm.normalize();

            centroids[index] = cent;
            normals[index] = norm;
            index += 1;

            printf("%d\n", index);
        }

        // Update vertex positions
        index = 0;
        for (auto it = mesh.v_begin(); it != mesh.v_end(); ++it) {
            const Vector pt = it->pt();
            Vector e = centroids[index] - pt;
            e -= normals[index] * e.dot(normals[index]);
            it->pt() = pt + e;

            index += 1;
        }
    }
}

}  // namespace tinymesh