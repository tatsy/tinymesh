#ifdef _MSC_VER
#pragma once
#endif

#ifndef TINYMESH_MATRIX_H
#define TINYMESH_MATRIX_H

#include <memory>

#include "common.h"

namespace tinymesh {

class TINYMESH_EXPORTS Matrix {
public:
    Matrix();
    Matrix(int rows, int cols);
    Matrix(int rows, int cols, double *m);
    Matrix(const Matrix &m);
    Matrix(Matrix &&m) noexcept;
    virtual ~Matrix();

    Matrix &operator=(const Matrix &m);
    Matrix &operator=(Matrix &&M) noexcept;

    Matrix &operator+=(const Matrix &m);
    Matrix &operator*=(double s);

    static Matrix ones(int rows, int cols);
    static Matrix zeros(int rows, int cols);
    static Matrix constant(int rows, int cols, double value);

private:
    int m_rows, m_cols;
    std::unique_ptr<double[]> m_values;
};

}  // namespace tinymesh

TINYMESH_EXPORTS tinymesh::Matrix operator+(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2);
TINYMESH_EXPORTS tinymesh::Matrix operator*(const tinymesh::Matrix &m, double s);
TINYMESH_EXPORTS tinymesh::Matrix operator*(double s, const tinymesh::Matrix &m);

#endif  // TINYMESH_MATRIX_H
