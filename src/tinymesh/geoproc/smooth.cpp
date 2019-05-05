#define TINYMESH_API_EXPORT
#include "smooth.h"

#include <algorithm>

#include "core/openmp.h"
#include "trimesh/face.h"
#include "trimesh/vertex.h"
#include "trimesh/halfedge.h"

namespace tinymesh {

void smooth(Mesh &mesh, double strength) {
    // Volonoi tessellation
    const int nv = mesh.num_vertices();
    std::vector<Vec> centroids(nv);
    std::vector<Vec> normals(nv);

    // Compute centroids and tangent planes
    omp_parallel_for (int i = 0; i < mesh.num_vertices(); i++) {
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
            centroids[i] = cent;
            normals[i] = norm / l;
        }
    }

    // Update vertex positions
    omp_parallel_for (int i = 0; i < mesh.num_vertices(); i++) {
        Vertex *v = mesh.vertex(i);
        if (v->isBoundary()) {
            continue;
        }

        if (length(normals[i]) != 0.0) {
            const Vec pt = v->pos();
            Vec e = centroids[i] - pt;
            e -= normals[i] * dot(e, normals[i]);
            Vec newPos = (1.0 - strength) * v->pos() + strength * (pt + e);
            v->setPos(newPos);
        }
    }
}

}  // namespace tinymesh