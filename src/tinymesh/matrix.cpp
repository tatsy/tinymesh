#define TINYMESH_API_EXPORT
#include "matrix.h"

#include <algorithm>

namespace tinymesh {

Matrix::Matrix()
    : m_rows{ 0 }
    , m_cols{ 0 } {
    m_values = nullptr;
}

Matrix::Matrix(int rows, int cols)
    : m_rows{ rows }
    , m_cols{ cols } {
    m_values = std::make_unique<double[]>(rows * cols);
}

Matrix::Matrix(int rows, int cols, double *m)
    : m_rows{ rows }
    , m_cols{ cols } {
    m_values = std::make_unique<double[]>(rows * cols);
    std::memcpy(m_values.get(), m, sizeof(double) * rows * cols);
}

Matrix::Matrix(const Matrix &m)
    : Matrix{} {
    this->operator=(m);
}

Matrix::Matrix(Matrix &&m) noexcept
    : Matrix{} {
    this->operator=(std::move(m));
}

Matrix::~Matrix() {
}

Matrix &Matrix::operator=(const Matrix &m) {
    this->m_rows = m.m_rows;
    this->m_cols = m.m_cols;
    this->m_values = std::make_unique<double[]>(m_rows * m_cols);
    std::memcpy(m_values.get(), m.m_values.get(), sizeof(double) * m.m_rows * m.m_cols);
    return *this;
}

Matrix &Matrix::operator=(Matrix &&m) noexcept {
    this->m_rows = m.m_rows;
    this->m_cols = m.m_cols;
    this->m_values = std::move(m.m_values);
    return *this;
}

Matrix &Matrix::operator+=(const Matrix &m) {
    if (m_rows != m.m_rows || m_cols != m.m_cols) {
        FatalError("Matrix size does not match!");
    }

    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            m_values[i * m_cols + j] += m.m_values[i * m_cols + j];
        }
    }

    return *this;
}

Matrix &Matrix::operator*=(double s) {
    for (int i = 0; i < m_rows * m_cols; i++) {
        m_values[i] *= s;
    }
    return *this;
}

double Matrix::det() const {
    Matrix lu = factorLU();
    double ret = 1.0;
    for (int i = 0; i < m_rows; i++) {
        ret *= lu.get(i, i);
    }
    return ret;
}

Matrix Matrix::solve(const Matrix &b) const {
    if (m_rows != m_cols) {
        FatalError("Matrix is not square. Cannot factorize.");
    }

    if (m_cols != b.m_rows) {
        FatalError("Matrix size is invalid");
    }

    int m = m_rows;
    int n = b.m_cols;

    // LU factorization
    int* order = new int[m];
    Matrix LU = factorLU(order);

    // Compute determinant
    double d = 1.0;
    for(int i=0; i<m; i++) {
        d *= LU.get(i, i);
    }

    if (d == 0.0) {
        FatalError("Matrix is singular. Cannot solve.");
    }

    // Reorder entries following pivot selections
    Matrix x(m_rows, b.m_cols);
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            x.set(i, j, b.get(order[i], j));
        }
    }

    // solve for L
    for(int j = 0; j < m; j++) {
        for(int i = j + 1; i < m; i++) {
            for(int k = 0; k < n; k++) {
                x.set(i, k, x.get(i, k) - x.get(j, k) * LU.get(i, j));
            }
        }
    }

    // solve for U
    for(int j = m - 1; j >= 0; j--) {
        for(int k = 0; k < n; k++) {
            x.set(j, k, x.get(j, k) / LU.get(j, j));
            for(int i = 0; i < j; i++) {
                x.set(i, k, x.get(i, k) - x.get(j, k) * LU.get(i, j));
            }
        }
    }

    delete[] order;

    return std::move(x);
}

Matrix Matrix::T() const {
    Matrix ret(m_cols, m_rows);
    for (int i = 0; i < m_rows; i++) {
        for (int j = 0; j < m_cols; j++) {
            ret.set(j, i, get(i, j));
        }
    }
    return std::move(ret);
}

Matrix Matrix::identity(int rows, int cols) {
    Matrix ret = Matrix::zeros(rows, cols);
    for (int i = 0; i < std::min(rows, cols); i++) {
        ret.set(i, i, 1.0);
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
    std::fill(m.m_values.get(), m.m_values.get() + rows * cols, value);
    return std::move(m);
}

Matrix Matrix::factorLU(int *order) const {
    if (m_rows != m_cols) {
        FatalError("Matrix is not square. Cannot factorize.");
    }

    int n = m_rows;
    Matrix LU = (*this);
    for(int i=0; i<n; i++) {
        if (order) {
            order[i] = i;
        }
    }

    for(int k=0; k<n; k++) {
        // Pivot selection
        double maxval = 0.0;
        int pivot  = k;
        for(int i=k; i<n; i++) {
            if(maxval < std::abs(LU.get(i, k))) {
                maxval = std::abs(LU.get(i, k));
                pivot  = i;
            }
        }

        // 行の入れ替え
        if (order) {
            std::swap(order[k], order[pivot]);
        }

        if(pivot != k) {
            for(int j=0; j<n; j++) {
                double tmp = LU.get(k, j);
                LU.set(k, j, LU.get(pivot, j));
                LU.set(pivot, j, tmp);
            }
        }

        // 要素の消去
        double iukk = 1.0 / LU.get(k, k);
        for(int i=k+1; i<n; i++) {
            double v = LU.get(i, k) * iukk;
            LU.set(i, k, v);
            for(int j=k+1; j<n; j++) {
                double v = LU.get(i, j) - LU.get(i, k) * LU.get(k, j);
                LU.set(i, j, v);
            }
        }
    }

    return std::move(LU);
}

}  // namespace tinymesh

tinymesh::Matrix operator+(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2) {
    tinymesh::Matrix ret = m1;
    ret += m2;
    return std::move(ret);
}

tinymesh::Matrix operator*(const tinymesh::Matrix &m, double s) {
    tinymesh::Matrix ret = m;
    ret *= s;
    return std::move(ret);
}

tinymesh::Matrix operator*(double s, const tinymesh::Matrix &m) {
    tinymesh::Matrix ret = m;
    ret *= s;
    return std::move(ret);
}

tinymesh::Matrix operator*(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2) {
    if (m1.cols() != m2.rows()) {
        FatalError("Matrix sizes are invalid: m1.cols (%d) !=  m2.rows (%d)", m1.cols(), m2.rows());
    }

    tinymesh::Matrix ret(m1.rows(), m2.cols());
    for (int i = 0; i < m1.rows(); i++) {
        for (int j = 0; j < m2.cols(); j++) {
            double v = 0.0;
            for (int k = 0; k < m1.cols(); k++) {
                v += m1.get(i, k) * m2.get(k, j);
            }
            ret.set(i, j, v);
        }
    }

    return std::move(ret);
}

