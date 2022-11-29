#include <iostream>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace tms = tinymesh;

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("usage: %s [input mesh]\n", argv[0]);
        return 1;
    }

    // Basename of input file
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();

    // Hole filling with minimum dihedral angle [Liepa 2003]
    {
        // tms::Mesh mesh(argv[1]);
        // for (int i = 0; i < mesh.numFaces(); i++) {
        //     tms::Face *f = mesh.face(i);
        //     if (f->isHole()) {
        //         holeFillMinDihedral(mesh, f, Pi / 6.0);
        //     }
        // }
        // const std::string outfile = (dirpath / fs::path((basename + "_fill_minDA" + extension).c_str())).string();
        // mesh.save(outfile);
        // printf("Save: %s\n", outfile.c_str());
    }

    // Hole filling with advancing front algorithm [Zhao et al. 2007]
    {
        tms::Mesh mesh(argv[1]);
        for (int i = 0; i < mesh.numFaces(); i++) {
            tms::Face *f = mesh.face(i);
            if (f->isHole()) {
                holeFillAdvancingFront(mesh, f);
                break;
            }
        }
        const std::string outfile = (dirpath / fs::path((basename + "_fill_adv_front" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }
}
