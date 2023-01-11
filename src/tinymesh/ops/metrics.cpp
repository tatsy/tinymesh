#define TINYMESH_API_EXPORT
#include "ops.h"

#include <algorithm>

#include "core/bvh.h"
#include "core/vertex.h"
#include "core/mesh.h"

namespace tinymesh {

double getHausdorffDistance(const Mesh &m0, const Mesh &m1) {
    EigenVector dist01 = getPerVertexShortestDistances(m0, m1);
    EigenVector dist10 = getPerVertexShortestDistances(m1, m0);

    const double dist01max = dist01.maxCoeff();
    const double dist10max = dist10.maxCoeff();
    return std::max(dist01max, dist10max);
}

EigenVector getPerVertexShortestDistances(const Mesh &src, const Mesh &dst) {
    const int N = (int)src.numVertices();
    EigenVector distances(N);

    BVH bvh(dst);
    for (size_t i = 0; i < N; i++) {
        const Vec3 query = src.vertex(i)->pos();
        distances(i) = bvh.distance(query);
    }

    return distances;
}

}  // namespace tinymesh
