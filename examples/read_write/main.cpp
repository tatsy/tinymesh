#include <iostream>
#include <experimental/filesystem>

#include "tinymesh/tinymesh.h"

namespace fs = std::experimental::filesystem;

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "usage: read_write [input mesh]" << std::endl;
        return 1;
    }

    // Load
    tinymesh::Mesh mesh(argv[1]);

    // Save
    const fs::path filepath = fs::canonical(fs::path(argv[1]));
    const fs::path dirpath = filepath.parent_path();
    const std::string extension = filepath.extension().string();
    const std::string basename = filepath.stem().string();
    const std::string outfile = (dirpath / fs::path((basename + "_copy" + extension).c_str())).string();
    mesh.save(outfile);
    printf("Save: %s\n", outfile.c_str());
}
