#include "gtest/gtest.h"

#include "tinymesh/tinymesh.h"
using namespace tinymesh;

class MatrixTest : public ::testing::Test {
protected:
    MatrixTest() {}
    virtual ~MatrixTest() {}
};

class MatrixPairwiseTest : public MatrixTest {
protected:
    MatrixPairwiseTest() {}
    virtual ~MatrixPairwiseTest() {}

    void SetUp() override {
        m1 = Matrix(100, 100);
        m2 = Matrix(100, 100);
    }

    void TearDown() override {

    }

    Matrix m1, m2;
};

TEST_F(MatrixTest, Initialization) {
    Matrix m;
    static const int rows = 123;
    static const int cols = 234;
    
    // Init with zeros
    m = Matrix::zeros(rows, cols);
    ASSERT_EQ(m.rows(), rows);
    ASSERT_EQ(m.cols(), cols);

    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            ASSERT_EQ(m(i, j), 0.0);
        }
    }

    ASSERT_DEATH(m(rows - 1, cols), "");
    ASSERT_DEATH(m(-1, cols - 1), "");

    // Init with constants
    static const double pi = 3.141592653589;
    m = Matrix::constant(rows, cols, pi);
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            ASSERT_EQ(m(i, j), pi);
        }
    }

    // Init with identity
    m = Matrix::identity(rows, cols);
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            if (i == j) {
                ASSERT_EQ(m(i, j), 1.0);
            } else {
                ASSERT_EQ(m(i, j), 0.0);
            }
        }
    }
}

TEST_F(MatrixTest, MatrixOps) {
    Matrix m;
    static const int size = 123;

    m = Matrix::random(size, size);
    Matrix minv = m.inverse();

    // Matrix inversion
    Matrix I = m * minv;
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            if (i == j) {
                ASSERT_NEAR(I(i, j), 1.0, 1.0e-12);
            } else {
                ASSERT_NEAR(I(i, j), 0.0, 1.0e-12);
            }
        }
    }

    I = minv * m;
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            if (i == j) {
                ASSERT_NEAR(I(i, j), 1.0, 1.0e-12);
            } else {
                ASSERT_NEAR(I(i, j), 0.0, 1.0e-12);
            }
        }
    }

    m = Matrix::random(size, size + 1);
    ASSERT_DEATH(m.inverse(), "");

    // Matrix transposition
    m = Matrix::random(size, size + 1);
    Matrix trans = m.T();
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            ASSERT_EQ(m(i, j), trans(j, i));
        }
    }
}