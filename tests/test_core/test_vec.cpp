#include "gtest/gtest.h"

#include <cstdio>
#include <sstream>
#include <tuple>
#include <algorithm>

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

using VecPair = std::tuple<Vec, Vec>; 

class VecTest : public ::testing::Test {
protected:
    VecTest() {}
    virtual ~VecTest() {}
};

class VecUnaryTest : public VecTest, public ::testing::WithParamInterface<Vec> {
protected:
    Vec v1;

protected:
    VecUnaryTest() {}
    virtual ~VecUnaryTest() {}

    void SetUp() {
        v1 = GetParam();
    }
};

class VecPairwiseTest : public VecTest, public ::testing::WithParamInterface<VecPair> {
protected:
    Vec v1, v2;

protected:
    VecPairwiseTest() {}
    virtual ~VecPairwiseTest() {}

    void SetUp() {
        v1 = std::get<0>(GetParam());
        v2 = std::get<1>(GetParam());
    }
};

TEST_F(VecTest, DefaultInstance) {
    Vec v;
    EXPECT_EQ(0.0, v.x);
    EXPECT_EQ(0.0, v.y);
    EXPECT_EQ(0.0, v.z);
}

TEST_P(VecUnaryTest, Instance) {
    Vec u(v1.x, v1.y, v1.z);
    EXPECT_EQ(v1.x, u.x);
    EXPECT_EQ(v1.y, u.y);
    EXPECT_EQ(v1.z, u.z);

    Vec v = { v1.x, v1.y, v1.z };
    EXPECT_EQ(v1.x, v.x);
    EXPECT_EQ(v1.y, v.y);
    EXPECT_EQ(v1.z, v.z);
}

TEST_P(VecUnaryTest, Copy) {
    Vec v = Vec(v1);
    EXPECT_EQ(v1.x, v.x);
    EXPECT_EQ(v1.y, v.y);
    EXPECT_EQ(v1.z, v.z);
}

TEST_P(VecUnaryTest, Assignment) {
    Vec v;
    v = v1;
    EXPECT_EQ(v1.x, v.x);
    EXPECT_EQ(v1.y, v.y);
    EXPECT_EQ(v1.z, v.z);
}

TEST_P(VecPairwiseTest, PlusOperator) {
    Vec w = v1 + v2;
    EXPECT_EQ(v1.x + v2.x, w.x);
    EXPECT_EQ(v1.y + v2.y, w.y);
    EXPECT_EQ(v1.z + v2.z, w.z);
}

TEST_P(VecPairwiseTest, MinusOperator) {
    Vec w = v1 - v2;
    EXPECT_EQ(v1.x - v2.x, w.x);
    EXPECT_EQ(v1.y - v2.y, w.y);
    EXPECT_EQ(v1.z - v2.z, w.z);
}

TEST_P(VecPairwiseTest, Multiplication) {
    const Vec v = v1 * v2;
    EXPECT_EQ(v1.x * v2.x, v.x);
    EXPECT_EQ(v1.y * v2.y, v.y);
    EXPECT_EQ(v1.z * v2.z, v.z);

    const double a = v2.x;
    const Vec u = v1 * a;
    EXPECT_EQ(v1.x * a, u.x);
    EXPECT_EQ(v1.y * a, u.y);
    EXPECT_EQ(v1.z * a, u.z);

    const double b = v2.y;
    const Vec w = b * v1;
    EXPECT_EQ(b * v1.x, w.x);
    EXPECT_EQ(b * v1.y, w.y);
    EXPECT_EQ(b * v1.z, w.z);
}

TEST_P(VecPairwiseTest, Division) {
    if (v2.x != 0.0 && v2.y != 0.0 && v2.z != 0.0) {
        const Vec v = v1 / v2;
        EXPECT_EQ(v1.x / v2.x, v.x);
        EXPECT_EQ(v1.y / v2.y, v.y);
        EXPECT_EQ(v1.z / v2.z, v.z);
    } else {
        ASSERT_DEATH(v1 / v2, "");
    }

    const double a = v2.x;
    if (a != 0.0) {
        const Vec u = v1 / a;
        EXPECT_DOUBLE_EQ(v1.x / a, u.x);
        EXPECT_DOUBLE_EQ(v1.y / a, u.y);
        EXPECT_DOUBLE_EQ(v1.z / a, u.z);
    } else {
        ASSERT_DEATH(v1 / a, "");
    }
}

TEST_P(VecUnaryTest, NormAndNormalize) {
    const double sqnrm = v1.x * v1.x + 
                         v1.y * v1.y +
                         v1.z * v1.z;
    const double nrm = sqrt(sqnrm);
    EXPECT_EQ(sqrt(sqnrm), length(v1));

    if (nrm != 0.0) {
        const Vec w = normalize(v1);
        EXPECT_EQ(v1.x / nrm, w.x);
        EXPECT_EQ(v1.y / nrm, w.y);
        EXPECT_EQ(v1.z / nrm, w.z);
        EXPECT_FLOAT_EQ(length(w), 1.0);
    } else {
        ASSERT_DEATH(normalize(v1), "");
    }
}

TEST_P(VecUnaryTest, Negation) {
    Vec v = -v1;
    EXPECT_EQ(-v1.x, v.x);
    EXPECT_EQ(-v1.y, v.y);
    EXPECT_EQ(-v1.z, v.z);
}

TEST_P(VecPairwiseTest, Equal) {
    EXPECT_EQ(v1.x == v2.x && v1.y == v2.y && v1.z == v2.z, v1 == v2);
}

TEST_P(VecPairwiseTest, DotAndCross) {
    double dt = dot(v1, v2);
    EXPECT_EQ(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z, dt);
    EXPECT_EQ(dt, dot(v2, v1));

    const Vec w = cross(v1, v2);
    EXPECT_EQ(v1.y * v2.z - v1.z * v2.y, w.x);
    EXPECT_EQ(v1.z * v2.x - v1.x * v2.z, w.y);
    EXPECT_EQ(v1.x * v2.y - v1.y * v2.x, w.z);
}

std::vector<Vec> vectors = {
    Vec(0.0, 1.0, 2.0),
    Vec(-2.0, -1.0, 0.0),
    Vec(3.14, 1.59, 2.65),
    Vec(1.0e8, 1.0e8, 1.0e8),
    Vec(1.0e-8, 1.0e-8, 1.0e-8)
};

INSTANTIATE_TEST_CASE_P(, VecUnaryTest,
    ::testing::ValuesIn(vectors));

INSTANTIATE_TEST_CASE_P(, VecPairwiseTest,
    ::testing::Combine(::testing::ValuesIn(vectors),
::testing::ValuesIn(vectors)));