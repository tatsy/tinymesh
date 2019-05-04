#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MATRIX_H
#define TINYMESH_MATRIX_H

#include <memory>

#include "core/common.h"

namespace tinymesh {

class TINYMESH_API Matrix {
public:
    Matrix();
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, double *const m);
    Matrix(const Matrix &m);
    Matrix(Matrix &&m) noexcept;
    virtual ~Matrix();

    Matrix &operator=(const Matrix &m);
    Matrix &operator=(Matrix &&M) noexcept;

    Matrix &operator+=(const Matrix &m);
    Matrix &operator*=(double s);
    Matrix &operator/=(double s);

    double get(int i, int j) const {
        return m_values[i * m_cols + j];
    }

    void set(int i, int j, double v) {
        m_values[i * m_cols + j] = v;
    }

    double det() const;
    Matrix solve(const Matrix &v) const;
    Matrix T() const;

    static Matrix identity(int rows, int cols);
    static Matrix ones(int rows, int cols);
    static Matrix zeros(int rows, int cols);
    static Matrix constant(int rows, int cols, double value);

    int rows() const { return m_rows; }
    int cols() const { return m_cols; }

private:
    Matrix factorLU(int *order = nullptr) const;

    int m_rows, m_cols;
    std::unique_ptr<double[]> m_values;
};

}  // namespace tinymesh

TINYMESH_API tinymesh::Matrix operator+(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2);
TINYMESH_API tinymesh::Matrix operator*(const tinymesh::Matrix &m, double s);
TINYMESH_API tinymesh::Matrix operator*(double s, const tinymesh::Matrix &m);
TINYMESH_API tinymesh::Matrix operator*(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2);
TINYMESH_API tinymesh::Matrix operator/(const tinymesh::Matrix &m, double s);

#endif  // TINYMESH_MATRIX_H
