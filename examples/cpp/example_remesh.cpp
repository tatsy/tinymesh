#include <iostream>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace tms = tinymesh;

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("usage: %s [input mesh] [keep angle]\n", argv[0]);
        return 1;
    }

    // Load
    tms::Mesh mesh(argv[1]);

    // Fill holes & remesh
    const double keepAngle = argc > 2 ? std::atof(argv[2]) * Pi / 180.0 : 0.0;
    mesh.fillHoles();
    tms::remeshTriangular(mesh, 0.8, 1.333, keepAngle);

    // Save
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();
    const std::string outfile = (dirpath / fs::path((basename + "_remesh" + extension).c_str())).string();
    mesh.save(outfile);
    printf("Save: %s\n", outfile.c_str());
}
