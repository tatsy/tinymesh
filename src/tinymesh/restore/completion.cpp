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

namespace tinymesh {

void holeFillMinDihedral(Mesh &mesh, Face *face, double dihedralBound) {
    mesh.triangulate(face, dihedralBound);
}

void holeFillAdvancingFront(Mesh &mesh, Face *face) {
}

void holeFillingContextCoherent(Mesh &mesh) {
}

}  // namespace tinymesh
