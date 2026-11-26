#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

class MeshRemeshTest : public TinyMeshBaseTestWithParam {
public:
    MeshRemeshTest() = default;
    virtual ~MeshRemeshTest() = default;

    void SetUp() override {
        mesh.load(filename);
    }
};

TEST_P(MeshRemeshTest, HoleFill) {
    EXPECT_NO_FATAL_FAILURE(mesh.fillHoles());
    EXPECT_TRUE(mesh.verify());
}

TEST_P(MeshRemeshTest, RemeshTriangular) {
    EXPECT_NO_FATAL_FAILURE(mesh.fillHoles());
    EXPECT_NO_FATAL_FAILURE(remeshTriangular(mesh));
}

TEST_P(MeshRemeshTest, SimplifyQEMOneHalf) {
    EXPECT_NO_FATAL_FAILURE(mesh.fillHoles());
    EXPECT_NO_FATAL_FAILURE(simplifyQEM(mesh, (int)mesh.numFaces() / 2));
}

TEST_P(MeshRemeshTest, SimplifyQEMOneTenth) {
    EXPECT_NO_FATAL_FAILURE(mesh.fillHoles());
    EXPECT_NO_FATAL_FAILURE(simplifyQEM(mesh, (int)mesh.numFaces() / 10));
}

static std::vector<std::string> filenames = {
    "box.obj",
    "torus.obj",
    "fandisk.ply",
    "bunny_mini.ply",
};

INSTANTIATE_TEST_SUITE_P(, MeshRemeshTest, ::testing::ValuesIn(filenames));

TEST(MeshHoleFillTest, ClearsBothHalfedgesOfFilledBoundary) {
    Mesh mesh({ Vec3{ 0, 0, 0 }, Vec3{ 1, 0, 0 }, Vec3{ 0, 1, 0 }, Vec3{ 0, 0, 1 } }, { 0, 2, 1, 0, 1, 3, 0, 3, 2 });

    ASSERT_EQ(mesh.getVertexIndices().size(), 9u);
    mesh.fillHoles();
    EXPECT_EQ(mesh.getVertexIndices().size(), 12u);

    for (int i = 0; i < static_cast<int>(mesh.numHalfedges()); ++i) {
        const Halfedge *he = mesh.halfedge(i);
        ASSERT_NE(he->rev(), nullptr);
        EXPECT_FALSE(he->isBoundary()) << i;
        EXPECT_FALSE(he->rev()->isBoundary()) << i;
    }
    EXPECT_TRUE(mesh.verify());
}
