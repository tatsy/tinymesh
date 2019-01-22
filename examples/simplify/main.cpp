#include <iostream>
#include <tinymesh/mesh.h>

#include "tinymesh/tinymesh.h"

int main(int argc, char **argv) {
    if (argc <= 1) {
        std::cout << "usage: hello_tinymesh [input mesh]" << std::endl;
        return 1;
    }

    // Load
    tinymesh::Mesh mesh(argv[1]);

    // Simplify
    tinymesh::simplify(mesh);

    // Save
    mesh.save("output.obj");
}
