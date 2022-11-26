#define TINYMESH_API_EXPORT
#include "restore.h"

#include "core/mesh.h"
#include "core/vertex.h"
#include "core/face.h"

namespace tinymesh {

void holeFill(Mesh &mesh, double dihedralBound) {
    std::vector<int> holeFaces;
    for (int i = 0; i < (int)mesh.numFaces(); i++) {
        Face *f = mesh.face(i);

        bool isHole = true;
        for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
            if (!vit->isBoundary()) {
                isHole = false;
                break;
            }
        }

        if (isHole) {
            holeFaces.push_back(i);
        }
    }

    for (auto i : holeFaces) {
        Face *f = mesh.face(i);
        mesh.triangulate(f, dihedralBound);
    }
}

void holeFillAdvancingFront(Mesh &mesh) {
}

void holeFillingContextCoherent(Mesh &mesh) {
}

}  // namespace tinymesh
