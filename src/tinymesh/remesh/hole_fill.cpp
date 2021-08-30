#define TINYMESH_API_EXPORT
#include "remesh.h"

#include "core/mesh.h"
#include "core/face.h"

namespace tinymesh {

void holeFill(Mesh& mesh, double dihedralBound) {
    for (int i = 0; i < mesh.numFaces(); i++) {
        Face *f = mesh.face(i);
        if (f->isBoundary()) {
            mesh.triangulate(f, dihedralBound);
        }
    }
}

}  // namespace tinymesh
