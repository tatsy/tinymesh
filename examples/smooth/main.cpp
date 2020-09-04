#include <iostream>
#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace mesh = tinymesh;

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "usage: read_write [input mesh] [ratio]" << std::endl;
        return 1;
    }

    const double ratio = argc <= 2 ? 0.1 : atof(argv[2]);

    // Basename of input file
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();

    // Smooth (Laplacian)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::laplace_smooth(mesh, 1.0, 20);
        const std::string outfile =
            (dirpath / fs::path((basename + "_laplace_smooth" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }

    // Smooth (Taubin)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::taubin_smooth(mesh, 1.0, 1.0,20);
        const std::string outfile =
            (dirpath / fs::path((basename + "_taubin_smooth" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }

    // Smooth (implicit fairing)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh);
        mesh::implicit_fair(mesh, 1.0e-4, 1);
        const std::string outfile =
            (dirpath / fs::path((basename + "_implicit_fair" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }
}
