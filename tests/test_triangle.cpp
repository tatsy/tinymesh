#include <gtest/gtest.h>

#include "tinymesh/tinymesh.h"

using namespace tinymesh;

TEST(TriangleTest, ClosestPoint) {
    Triangle triangle(Vec3{ 0.0, 0.0, 0.0 }, Vec3{ 1.0, 0.0, 0.0 }, Vec3{ 0.0, 1.0, 0.0 });

    const Vec3 p{ 0.25, 0.25, 0.0 };
    const Vec3 cp = triangle.closestPoint(p);
    EXPECT_DOUBLE_EQ(cp.x(), p.x());
    EXPECT_DOUBLE_EQ(cp.y(), p.y());
    EXPECT_DOUBLE_EQ(cp.z(), p.z());
    EXPECT_DOUBLE_EQ(triangle.distance(p), 0.0);
}
