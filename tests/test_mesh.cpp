#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <vector>
#include <tuple>
#include <random>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

#include "test_utils.h"

class MeshBasicTest : public TinyMeshBaseTestWithParam {
public:
    MeshBasicTest() = default;
    virtual ~MeshBasicTest() = default;
};

TEST(MeshBaseTest, MeshInvaidLoad) {
    Mesh mesh;
    ASSERT_DEATH(mesh.load("mesh_not_found.ply"), "");
}

TEST_P(MeshBasicTest, MeshLoad) {
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

TEST_P(MeshBasicTest, MeshClosetPoint) {
    std::random_device randev;
    std::mt19937 mt(randev());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const Vec3 query(dist(mt), dist(mt), dist(mt));

    std::set_terminate([] { ADD_FAILURE(); });
    ASSERT_NO_FATAL_FAILURE(mesh.load(filename));
    mesh.fillHoles();

    double answer = 1.0e20;
    for (int i = 0; i < mesh.numFaces(); i++) {
        const Triangle t = mesh.face(i)->toTriangle();
        const double dist = t.distance(query);
        if (dist < answer) {
            answer = dist;
        }
    }

    BVH bvh(mesh);

    EXPECT_NEAR(bvh.distance(query), answer, 1.0e-5);
}

static std::vector<std::string> filenames = {
    "box.obj", "sphere.obj", "torus.obj", "plane.obj", "fandisk.ply", "bunny.ply",
};

INSTANTIATE_TEST_SUITE_P(, MeshBasicTest, ::testing::ValuesIn(filenames));
