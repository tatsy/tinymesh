#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_config.h"

TEST(MeshTest, MeshLoad) {
    Mesh mesh;
    mesh.load(std::string(DATA_DIRECTORY) + "models/bunny.ply");

    ASSERT_DEATH(mesh.load("mesh_not_found.ply"), "");
}
