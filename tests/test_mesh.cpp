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

namespace {

const std::vector<Vec3> validationVertices = {
    Vec3{ 0, 0, 0 },  Vec3{ 1, 0, 0 },  Vec3{ 0, 1, 0 },  Vec3{ 0, 0, 1 },
    Vec3{ -1, 0, 0 }, Vec3{ 0, -1, 0 }, Vec3{ 0, 0, -1 },
};

}  // namespace

TEST(MeshValidationTest, AcceptsManifoldMeshes) {
    struct Case {
        const char *name;
        std::vector<uint32_t> indices;
        size_t holes;
    };

    const Case cases[] = {
        { "single triangle", { 0, 1, 2 }, 1 },
        { "two adjacent triangles", { 0, 1, 2, 0, 2, 3 }, 1 },
        { "tetrahedron", { 0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3 }, 0 },
    };

    for (const auto &c : cases) {
        SCOPED_TRACE(c.name);

        Mesh mesh(validationVertices, c.indices);
        EXPECT_TRUE(mesh.verify());
        EXPECT_EQ(mesh.getVertexIndices().size(), c.indices.size());

        size_t holes = 0;
        for (size_t i = 0; i < mesh.numFaces(); ++i) {
            if (mesh.face(static_cast<int>(i))->isHole()) {
                ++holes;
            }
        }

        EXPECT_EQ(holes, c.holes);
        EXPECT_EQ(mesh.numFaces(), c.indices.size() / 3 + c.holes);
    }
}

TEST(MeshValidationDeathTest, RejectsInvalidMeshes) {
    struct Case {
        const char *name;
        std::vector<uint32_t> indices;
        const char *message;
    };

    const Case cases[] = {
        { "three faces on one edge", { 0, 1, 2, 1, 0, 3, 0, 1, 5 }, "Non-manifold edge" },

        { "inconsistent winding", { 0, 1, 2, 0, 1, 3 }, "Inconsistent winding" },

        { "open bow-tie vertex", { 0, 1, 2, 0, 4, 5 }, "multiple outgoing boundary edges" },

        { "closed bow-tie vertex",
          { 0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3, 0, 5, 4, 0, 4, 6, 4, 5, 6, 5, 0, 6 },
          "disconnected face fans" },

        { "incomplete triangle", { 0, 1 }, "multiple of three" },

        { "invalid vertex index", { 0, 1, 7 }, "out of range" },

        { "repeated vertex", { 0, 1, 1 }, "three distinct vertices" },
    };

    for (const auto &c : cases) {
        SCOPED_TRACE(c.name);

        EXPECT_DEATH(
            {
                std::cout.rdbuf(std::cerr.rdbuf());
                Mesh mesh(validationVertices, c.indices);
            },
            c.message);
    }
}
