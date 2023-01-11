#include <iostream>
#include <random>

#include "tinymesh/tinymesh.h"

namespace fs = std::filesystem;
namespace tms = tinymesh;

int main(int argc, char **argv) {
    try {
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
        std::uniform_real_distribution<double> dist(-0.5, 0.5);

        // Save noise mesh
        std::string noiseMeshFile;
        {
            tms::Mesh mesh(argv[1]);
            mesh.fillHoles();
            const double Lavg = mesh.getMeanEdgeLength();

            // Add noise
            for (int i = 0; i < (int)mesh.numVertices(); i++) {
                const Vec3 pos = mesh.vertex(i)->pos();
                const Vec3 noise = Vec3(dist(mt), dist(mt), dist(mt));
                mesh.vertex(i)->setPos(pos + Lavg * noise);
            }
            noiseMeshFile = (dirpath / fs::path((basename + "_noise" + extension).c_str())).string();
            mesh.save(noiseMeshFile);
        }

        // Denoise (Normal Gaussian filter)
        {
            tms::Mesh mesh(noiseMeshFile);
            mesh.fillHoles();
            tms::denoiseNormalGaussian(mesh, sigma, 10);

            const std::string outfile =
                (dirpath / fs::path((basename + "_denoise_Gaussian" + extension).c_str())).string();
            mesh.save(outfile);
            printf("Save: %s\n", outfile.c_str());
        }

        // Denoise (Normal bilateral filter)
        {
            tms::Mesh mesh(noiseMeshFile);
            mesh.fillHoles();
            tms::denoiseNormalBilateral(mesh, sigma, 0.1, 10);

            const std::string outfile =
                (dirpath / fs::path((basename + "_denoise_bilateral" + extension).c_str())).string();
            mesh.save(outfile);
            printf("Save: %s\n", outfile.c_str());
        }

        // Denoise (L0 smoothing)
        {
            tms::Mesh mesh(noiseMeshFile);
            mesh.fillHoles();
            tms::denoiseL0Smooth(mesh);

            const std::string outfile = (dirpath / fs::path((basename + "_denoise_l0" + extension).c_str())).string();
            mesh.save(outfile);
            printf("Save: %s\n", outfile.c_str());
        }
    } catch (std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }
}
