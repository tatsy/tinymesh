#define TINYMESH_API_EXPORT
#include "restore.h"

#include <queue>
#include <functional>
#include <unordered_map>

#include "core/mesh.h"
#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"
#include "core/eigen.h"

using IndexPair = std::pair<uint32_t, uint32_t>;

struct IndexPairHash : public std::function<size_t(IndexPair)> {
    std::size_t operator()(const IndexPair &k) const {
        return std::get<0>(k) ^ std::get<1>(k);
    }
};

namespace tinymesh {

void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound) {
    mesh.holeFillMinDihedral_(face, dihedralBound);
}

void holeFillAdvancingFront(Mesh &mesh, Face *face) {
    mesh.holeFillAdvancingFront_(face);
}

void holeFillingContextCoherent(Mesh &mesh) {
}

}  // namespace tinymesh
