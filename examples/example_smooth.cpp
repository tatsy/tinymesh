#include <iostream>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace mesh = tinymesh;

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

    // Smooth (Laplacian)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::smoothLaplacian(mesh, 0.5, true, 100);
        const std::string outfile = (dirpath / fs::path((basename + "_laplace_smooth" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }

    // Smooth (Taubin)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::smoothTaubin(mesh, 0.5, 0.53, 100);
        const std::string outfile = (dirpath / fs::path((basename + "_taubin_smooth" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }

    // Smooth (implicit fairing)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::implicitFairing(mesh, 1.0e-2, 1);
        const std::string outfile = (dirpath / fs::path((basename + "_implicit_fair" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }
}
