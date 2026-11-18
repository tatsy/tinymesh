#include "gtest/gtest.h"

#include <format>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

/*
 * Test for a simple mesh
 */
class SimpleMeshFilterTest : public TinyMeshBaseTest {
public:
    SimpleMeshFilterTest() = default;
    virtual ~SimpleMeshFilterTest() = default;

    void SetUp() override {
        const std::vector<Vec3> vertices = {
            Vec3(0.0, 0.0, 3.0), Vec3(1.0, 0.0, 2.0), Vec3(0.0, 1.0, 2.0), Vec3(-1.0, 0.0, 2.0), Vec3(0.0, -1.0, 2.0),
        };
        const std::vector<uint32_t> indices = {
            0, 1, 2, 0, 2, 3, 0, 3, 4, 0, 4, 1,
        };
        mesh.construct(vertices, indices);
    }
};

TEST_F(SimpleMeshFilterTest, SmoothLaplacian) {
    EXPECT_NO_FATAL_FAILURE(smoothLaplacian(mesh));
}
TEST_F(SimpleMeshFilterTest, SmoothTaubin) {
    EXPECT_NO_FATAL_FAILURE(smoothTaubin(mesh));
}
TEST_F(SimpleMeshFilterTest, ImplicitFairingKeepsBoundaryAndMovesCenter) {
    const auto before = mesh.getVertices();
    ASSERT_FALSE(mesh.vertex(0)->isBoundary());

    implicitFairing(mesh, 1.0, 1);

    const auto after = mesh.getVertices();
    ASSERT_EQ(after.size(), before.size());
    for (std::size_t i = 1; i < before.size(); ++i) {
        EXPECT_EQ(after[i], before[i]) << std::format("Boundary vertex {}", i);
    }
    EXPECT_NEAR(after[0].z(), 2.5, 1e-7);
}

TEST_F(SimpleMeshFilterTest, DenoiseNormalGaussian) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalGaussian(mesh));
}
TEST_F(SimpleMeshFilterTest, DenoiseNormalBilateral) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalBilateral(mesh));
}
TEST_F(SimpleMeshFilterTest, DenoiseL0Smoothing) {
    EXPECT_NO_FATAL_FAILURE(denoiseL0Smooth(mesh));
}

/**
 * Test for real meshes
 */
class RealMeshFilterTest : public TinyMeshBaseTestWithParam {
public:
    RealMeshFilterTest() = default;
    virtual ~RealMeshFilterTest() = default;

    void SetUp() override {
        mesh.load(filename);
        mesh.fillHoles();
    }
};

TEST_P(RealMeshFilterTest, SmoothLaplacian) {
    EXPECT_NO_FATAL_FAILURE(smoothLaplacian(mesh));
}
TEST_P(RealMeshFilterTest, SmoothTaubin) {
    EXPECT_NO_FATAL_FAILURE(smoothTaubin(mesh));
}
TEST_P(RealMeshFilterTest, ImplicitFairing) {
    EXPECT_NO_FATAL_FAILURE(implicitFairing(mesh));
}

TEST_P(RealMeshFilterTest, DenoiseNormalGaussian) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalGaussian(mesh));
}
TEST_P(RealMeshFilterTest, DenoiseNormalBilateral) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalBilateral(mesh));
}
TEST_P(RealMeshFilterTest, DenoiseL0Smoothing) {
    EXPECT_NO_FATAL_FAILURE(denoiseL0Smooth(mesh));
}

static std::vector<std::string> filenames = {
    "torus.obj",
    "fandisk.ply",
    "bunny_mini.ply",
};

INSTANTIATE_TEST_SUITE_P(, RealMeshFilterTest, ::testing::ValuesIn(filenames));
