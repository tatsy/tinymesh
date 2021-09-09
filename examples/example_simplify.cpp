#include <iostream>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("usage: %s [input mesh] [ratio]\n", argv[0]);
        return 1;
    }

    printf("Hello\n");

    const double ratio = argc <= 2 ? 0.1 : atof(argv[2]);

    // Load
    tinymesh::Mesh mesh(argv[1]);

    // Simplify
    const int target = (int)(ratio * mesh.numFaces());
    tinymesh::holeFill(mesh);
    tinymesh::simplifyQEM(mesh, target, 10, true);
    tinymesh::remeshTriangular(mesh, 0.0, 1e6);  // remesh by keeping #faces

    // Save
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();
    const std::string outfile = (dirpath / fs::path((basename + "_simplify" + extension).c_str())).string();
    mesh.save(outfile);
    printf("Save: %s\n", outfile.c_str());
}
