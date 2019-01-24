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

