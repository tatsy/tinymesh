#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

class MeshOpsTest : public TinyMeshBaseTestWithParam {
public:
    MeshOpsTest() = default;
    virtual ~MeshOpsTest() = default;

    void SetUp() override {
        mesh.load(filename);
        mesh.fillHoles();
    }
};

// Operators

TEST_P(MeshOpsTest, GetLaplacian) {
    ASSERT_NO_FATAL_FAILURE(getMeshLaplacian(mesh, MeshLaplace::Adjacent));
    ASSERT_NO_FATAL_FAILURE(getMeshLaplacian(mesh, MeshLaplace::Cotangent));
    ASSERT_NO_FATAL_FAILURE(getMeshLaplacian(mesh, MeshLaplace::Belkin08));
}

// Features

TEST_P(MeshOpsTest, GetHeatKernelSignatures) {
    auto L = getMeshLaplacian(mesh, MeshLaplace::Cotangent);
    const int K = std::min((int)L.rows(), 200);
    ASSERT_NO_FATAL_FAILURE(getHeatKernelSignatures(L, K));
}

TEST_P(MeshOpsTest, GetPrincipalCurvatures) {
    ASSERT_NO_FATAL_FAILURE(getPrincipalCurvatures(mesh));
}

static std::vector<std::string> filenames = {
    "torus.obj",
    "fandisk.ply",
    "bunny_mini.ply",
};

INSTANTIATE_TEST_SUITE_P(, MeshOpsTest, ::testing::ValuesIn(filenames));
