#include <iostream>
#include <random>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace mesh = tinymesh;

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("usage: %s [input mesh] [sigma]\n", argv[0]);
        return 1;
    }

    const double sigma = argc <= 2 ? 0.2 : atof(argv[2]);

    // Basename of input file
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();

    // Random number generator
    std::random_device randev;
    std::mt19937 mt(randev());
    std::uniform_real_distribution<double> dist(-0.01, 0.01);

    // Denoise (Normal Gaussian filter)
    {
        mesh::Mesh mesh(argv[1]);
        mesh::holeFill(mesh, Pi / 6.0);

        // Add noise
        for (int i = 0; i < mesh.num_vertices(); i++) {
            const Vec3 pos = mesh.vertex(i)->pos();
            const Vec3 noise = Vec3(dist(mt), dist(mt), dist(mt));
            mesh.vertex(i)->setPos(pos + noise);
        }
        std::string outfile;
        outfile = (dirpath / fs::path((basename + "_noise" + extension).c_str())).string();
        mesh.save(outfile);

        mesh::denoiseNormalGaussian(mesh, sigma, 10);

        outfile = (dirpath / fs::path((basename + "_denoise" + extension).c_str())).string();
        mesh.save(outfile);
        printf("Save: %s\n", outfile.c_str());
    }
}
