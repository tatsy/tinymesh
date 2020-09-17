#define TINYMESH_API_EXPORT
#include "hole_fill.h"

#include "polymesh/mesh.h"
#include "polymesh/face.h"

namespace tinymesh {

void hole_fill(Mesh& mesh, double dihedralBound) {
    for (int i = 0; i < mesh.num_faces(); i++) {
        Face *f = mesh.face(i);
        if (f->isBoundary()) {
            mesh.triangulate(f, dihedralBound);
        }
    }
}

}  // namespace tinymesh
