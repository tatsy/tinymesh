#define TINYMESH_API_EXPORT
#include "matrix.h"

#include <cstring>
#include <random>
#include <algorithm>

namespace tinymesh {

Matrix::Matrix()
    : rows_{ 0 }
    , cols_{ 0 } {
    values_ = nullptr;
}

Matrix::Matrix(int rows, int cols)
    : rows_{ rows }
    , cols_{ cols } {
    values_ = std::make_unique<double[]>(rows * cols);
}

Matrix::Matrix(int rows, int cols, double *const m)
    : rows_{ rows }
    , cols_{ cols } {
    values_ = std::make_unique<double[]>(rows * cols);
    std::memcpy(values_.get(), m, sizeof(double)* rows* cols);
}

Matrix::Matrix(const Matrix& m)
    : Matrix{} {
    this->operator=(m);
}

Matrix::Matrix(Matrix&& m) noexcept
    : Matrix{} {
    this->operator=(std::move(m));
}

Matrix::~Matrix() {
}

Matrix& Matrix::operator=(const Matrix& m) {
    this->rows_ = m.rows_;
    this->cols_ = m.cols_;
    this->values_ = std::make_unique<double[]>(rows_ * cols_);
    std::memcpy(values_.get(), m.values_.get(), sizeof(double) * m.rows_ * m.cols_);
    return *this;
}

Matrix& Matrix::operator=(Matrix&& m) noexcept {
    this->rows_ = m.rows_;
    this->cols_ = m.cols_;
    this->values_ = std::move(m.values_);
    return *this;
}

Matrix& Matrix::operator+=(const Matrix& m) {
    if (rows_ != m.rows_ || cols_ != m.cols_) {
        FatalError("Matrix size does not match!");
    }

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            values_[i * cols_ + j] += m.values_[i * cols_ + j];
        }
    }

    return *this;
}

Matrix& Matrix::operator*=(double s) {
    for (int i = 0; i < rows_ * cols_; i++) {
        values_[i] *= s;
    }
    return *this;
}

Matrix& Matrix::operator/=(double s) {
    return this->operator*=(1.0 / s);
}

double Matrix::det() const {
    Matrix lu = factorLU();
    double ret = 1.0;
    for (int i = 0; i < rows_; i++) {
        ret *= lu(i, i);
    }
    return ret;
}

Matrix Matrix::solve(const Matrix& b) const {
    if (rows_ != cols_) {
        FatalError("Matrix is not square. Cannot factorize.");
    }

    if (cols_ != b.rows_) {
        FatalError("Matrix size is invalid");
    }

    int m = rows_;
    int n = b.cols_;

    // LU factorization
    int* order = new int[m];
    Matrix lu = factorLU(order);

    // Compute determinant
    double d = 1.0;
    for (int i = 0; i < m; i++) {
        d *= lu(i, i);
    }

    if (d == 0.0) {
        FatalError("Matrix is singular. Cannot solve.");
    }

    // Reorder entries following pivot selections
    Matrix x(rows_, b.cols_);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            x(i, j) = b(order[i], j);
        }
    }

    // solve for L
    for (int j = 0; j < m; j++) {
        for (int i = j + 1; i < m; i++) {
            for (int k = 0; k < n; k++) {
                x(i, k) = x(i, k) - x(j, k) * lu(i, j);
            }
        }
    }

    // solve for U
    for (int j = m - 1; j >= 0; j--) {
        for (int k = 0; k < n; k++) {
            x(j, k) = x(j, k) / lu(j, j);
            for (int i = 0; i < j; i++) {
                x(i, k) = x(i, k) - x(j, k) * lu(i, j);
            }
        }
    }

    delete[] order;

    return std::move(x);
}

Matrix Matrix::T() const {
    Matrix ret(cols_, rows_);
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
            ret(j, i) = (*this)(i, j);
        }
    }
    return std::move(ret);
}

Matrix Matrix::inverse() const {
    Assertion(rows_ == cols_, "Non square matrix cannot be inverted!");
    Matrix I = Matrix::identity(rows_, cols_);
    return std::move(this->solve(I));
}

Matrix Matrix::identity(int rows, int cols) {
    Matrix ret = Matrix::zeros(rows, cols);
    for (int i = 0; i < std::min(rows, cols); i++) {
        ret(i, i) = 1.0;
    }
    return std::move(ret);
}

Matrix Matrix::zeros(int rows, int cols) {
    return std::move(Matrix::constant(rows, cols, 0.0));
}

Matrix Matrix::ones(int rows, int cols) {
    return std::move(Matrix::constant(rows, cols, 1.0));
}

Matrix Matrix::constant(int rows, int cols, double value) {
    Matrix m(rows, cols);
    std::fill(m.values_.get(), m.values_.get() + rows * cols, value);
    return std::move(m);
}

Matrix Matrix::random(int rows, int cols) {
    std::random_device randev;
    std::mt19937 mt(randev());
    std::uniform_real_distribution<double> rng;

    Matrix m(rows, cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            m(i, j) = rng(mt);
        }
    }
    return std::move(m);
}

Matrix Matrix::factorLU(int* order) const {
    if (rows_ != cols_) {
        FatalError("Matrix is not square. Cannot factorize.");
    }

    int n = rows_;
    Matrix lu = (*this);
    for (int i = 0; i < n; i++) {
        if (order) {
            order[i] = i;
        }
    }

    for (int k = 0; k < n; k++) {
        // Pivot selection
        double maxval = 0.0;
        int pivot = k;
        for (int i = k; i < n; i++) {
            if (maxval < std::abs(lu(i, k))) {
                maxval = std::abs(lu(i, k));
                pivot = i;
            }
        }

        // Swap columns
        if (order) {
            std::swap(order[k], order[pivot]);
        }

        if (pivot != k) {
            for (int j = 0; j < n; j++) {
                double tmp = lu(k, j);
                lu(k, j) = lu(pivot, j);
                lu(pivot, j) = tmp;
            }
        }

        // Eliminate entries
        double iukk = 1.0 / lu(k, k);
        for (int i = k + 1; i < n; i++) {
            double v = lu(i, k) * iukk;
            lu(i, k) = v;
            for (int j = k + 1; j < n; j++) {
                double v = lu(i, j) - lu(i, k) * lu(k, j);
                lu(i, j) = v;
            }
        }
    }

    return std::move(lu);
}

}  // namespace tinymesh

tinymesh::Matrix operator+(const tinymesh::Matrix & m1, const tinymesh::Matrix & m2) {
    tinymesh::Matrix ret = m1;
    ret += m2;
    return std::move(ret);
}

tinymesh::Matrix operator*(const tinymesh::Matrix & m, double s) {
    tinymesh::Matrix ret = m;
    ret *= s;
    return std::move(ret);
}

tinymesh::Matrix operator*(double s, const tinymesh::Matrix & m) {
    tinymesh::Matrix ret = m;
    ret *= s;
    return std::move(ret);
}

tinymesh::Matrix operator*(const tinymesh::Matrix & m1, const tinymesh::Matrix & m2) {
    if (m1.cols() != m2.rows()) {
        FatalError("Matrix sizes are invalid: m1.cols (%d) !=  m2.rows (%d)", m1.cols(), m2.rows());
    }

    tinymesh::Matrix ret(m1.rows(), m2.cols());
    for (int i = 0; i < m1.rows(); i++) {
        for (int j = 0; j < m2.cols(); j++) {
            double v = 0.0;
            for (int k = 0; k < m1.cols(); k++) {
                v += m1(i, k) * m2(k, j);
            }
            ret(i, j) = v;
        }
    }

    return std::move(ret);
}

tinymesh::Matrix operator/(const tinymesh::Matrix & m, double s) {
    tinymesh::Matrix ret = m;
    ret /= s;
    return std::move(ret);
}
