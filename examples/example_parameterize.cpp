#include <iostream>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace mesh = tinymesh;

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("usage: %s [input mesh]\n", argv[0]);
        return 1;
    }

    // Load
    mesh::Mesh mesh(argv[1]);
    for (int i = 0; i < mesh.numFaces(); i++) {
        mesh.face(i)->setIsStatic(true);
    }

    // Fill holes & remesh
    mesh::holeFill(mesh, Pi / 6.0);
    mesh::abfxx(mesh);

    // Save
    //const fs::path filepath = fs::canonical(fs::path(argv[1]));
    //const fs::path dirpath = filepath.parent_path();
    //const std::string extension = filepath.extension().string();
    //const std::string basename = filepath.stem().string();
    //const std::string outfile = (dirpath / fs::path((basename + "_remesh" + extension).c_str())).string();
    //mesh.save(outfile);
    //printf("Save: %s\n", outfile.c_str());
}
