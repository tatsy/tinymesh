#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

class MeshFilterTest : public TinyMeshBaseTestWithParam {
public:
    MeshFilterTest() = default;
    virtual ~MeshFilterTest() = default;

    void SetUp() override {
        mesh.load(filename);
        mesh.fillHoles();
    }
};

// Smoothing

TEST_P(MeshFilterTest, SmoothLaplacian) {
    EXPECT_NO_FATAL_FAILURE(smoothLaplacian(mesh));
}

TEST_P(MeshFilterTest, SmoothTaubin) {
    EXPECT_NO_FATAL_FAILURE(smoothTaubin(mesh));
}

TEST_P(MeshFilterTest, ImplicitFairing) {
    EXPECT_NO_FATAL_FAILURE(implicitFairing(mesh));
}

// Denoising

TEST_P(MeshFilterTest, DenoiseNormalGaussian) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalGaussian(mesh));
}

TEST_P(MeshFilterTest, DenoiseNormalBilateral) {
    EXPECT_NO_FATAL_FAILURE(denoiseNormalBilateral(mesh));
}

TEST_P(MeshFilterTest, DenoiseL0Smoothing) {
    EXPECT_NO_FATAL_FAILURE(denoiseL0Smooth(mesh));
}

static std::vector<std::string> filenames = {
    "torus.obj",
    "fandisk.ply",
    "bunny_mini.ply",
};

INSTANTIATE_TEST_SUITE_P(, MeshFilterTest, ::testing::ValuesIn(filenames));
