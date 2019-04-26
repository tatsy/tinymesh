#define TINYMESH_API_EXPORT
#include "smooth.h"

#include <algorithm>

#include "trimesh/face.h"
#include "trimesh/vertex.h"
#include "trimesh/halfedge.h"

namespace tinymesh {

void smooth(Mesh &mesh, double strength) {
    // Clip strength
    strength = std::max(0.0, std::min(strength, 1.0));

    // Simple Laplacian smoothing
    int index;
    const int nv = mesh.num_vertices();
    std::vector<Vec> centroids(nv);

    // Compute centroids and tangent planes
    index = 0;
    for (auto it = mesh.v_begin(); it != mesh.v_end(); ++it) {
        // Collect surrounding vertices
        Vec org = it->pt();
        std::vector<Vec> pts;
        for (auto vit = it->v_begin(); vit != it->v_end(); ++vit) {
            pts.push_back(vit->pt());
        }

        // Compute centroids, tangents, and binormals
        Vec cent(0.0);
        for (int i = 0; i < pts.size(); i++) {
            const int j = (i + 1) % pts.size();
            Vec e1 = pts[i] - org;
            Vec e2 = pts[j] - org;
            Vec g = (org + pts[i] + pts[j]) / 3.0;

            cent += g;
        }
        cent /= pts.size();

        centroids[index] = cent;
        index += 1;
    }

    // Update vertex positions
    index = 0;
    for (auto it = mesh.v_begin(); it != mesh.v_end(); ++it) {
        const Vec pt = it->pt();
        it->setPt((1.0 - strength) * pt + strength * centroids[index]);
        //it->setPt(centroids[index]);
        index += 1;
    }
}

}  // namespace tinymesh