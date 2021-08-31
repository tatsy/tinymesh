#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_config.h"

class MeshBasicTest : public ::testing::Test {
protected:
    MeshBasicTest() {}
    virtual ~MeshBasicTest() {}
};

class MeshTest : public MeshBasicTest, public ::testing::WithParamInterface<std::string> {
protected:
    MeshTest() {}
    virtual ~MeshTest() {}

    void SetUp() {
        filename = std::string(DATA_DIRECTORY) + "models/" + GetParam();
    }

    std::string filename;
};

TEST_F(MeshBasicTest, MeshInvaidLoad) {
    Mesh mesh;
    ASSERT_DEATH(mesh.load("mesh_not_found.ply"), "");
}

TEST_P(MeshTest, MeshLoad) {
    Mesh mesh;
    mesh.load(filename);
    ASSERT_TRUE(mesh.verify());

    ASSERT_GT(mesh.numVertices(), 0);
    ASSERT_GT(mesh.numFaces(), 0);
    ASSERT_GT(mesh.numEdges(), 0);
    ASSERT_GT(mesh.numHalfedges(), 0);

    const int nVerts = (int)mesh.numVertices();
    for (int i = 0; i < nVerts; i++) {
        ASSERT_GT(mesh.vertex(i)->degree(), 0);
    }
}

std::vector<std::string> filenames = {
    "box.obj",
    "sphere.obj",
    "torus.obj",
    "plane.obj",
};

INSTANTIATE_TEST_SUITE_P(, MeshTest, ::testing::ValuesIn(filenames));