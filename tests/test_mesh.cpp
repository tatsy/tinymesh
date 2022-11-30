#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

class MeshIOTest : public TinyMeshBaseTestWithParam {
public:
    MeshIOTest() = default;
    virtual ~MeshIOTest() = default;
};

TEST(MeshBaseTest, MeshInvaidLoad) {
    Mesh mesh;
    ASSERT_DEATH(mesh.load("mesh_not_found.ply"), "");
}

TEST_P(MeshIOTest, MeshLoad) {
    std::set_terminate([] { ADD_FAILURE(); });
    ASSERT_NO_FATAL_FAILURE(mesh.load(filename));
    ASSERT_TRUE(mesh.verify());

    EXPECT_GT(mesh.numVertices(), 0);
    EXPECT_GT(mesh.numFaces(), 0);
    EXPECT_GT(mesh.numEdges(), 0);
    EXPECT_GT(mesh.numHalfedges(), 0);

    const int nVerts = (int)mesh.numVertices();
    for (int i = 0; i < nVerts; i++) {
        EXPECT_GT(mesh.vertex(i)->degree(), 0);
    }
}

static std::vector<std::string> filenames = {
    "box.obj", "sphere.obj", "torus.obj", "plane.obj", "fandisk.ply", "bunny.ply",
};

INSTANTIATE_TEST_SUITE_P(, MeshIOTest, ::testing::ValuesIn(filenames));
