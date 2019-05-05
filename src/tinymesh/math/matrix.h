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
    // Public methods
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

    double operator()(int i, int j) const {
        Assertion(i >= 0 && i < rows_ && j >= 0 && j < cols_, "Matrix subscription out of range!");
        return values_[i * cols_ + j];
    }
    double &operator()(int i, int j) {
        Assertion(i >= 0 && i < rows_ && j >= 0 && j < cols_, "Matrix subscription out of range!");
        return values_[i * cols_ + j];
    }

    double det() const;
    Matrix solve(const Matrix &v) const;
    Matrix T() const;
    Matrix inverse() const;

    static Matrix identity(int rows, int cols);
    static Matrix ones(int rows, int cols);
    static Matrix zeros(int rows, int cols);
    static Matrix constant(int rows, int cols, double value);
    static Matrix random(int rows, int cols);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

private:
    Matrix factorLU(int *order = nullptr) const;

    int rows_, cols_;
    std::unique_ptr<double[]> values_ = nullptr;
};

}  // namespace tinymesh

TINYMESH_API tinymesh::Matrix operator+(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2);
TINYMESH_API tinymesh::Matrix operator*(const tinymesh::Matrix &m, double s);
TINYMESH_API tinymesh::Matrix operator*(double s, const tinymesh::Matrix &m);
TINYMESH_API tinymesh::Matrix operator*(const tinymesh::Matrix &m1, const tinymesh::Matrix &m2);
TINYMESH_API tinymesh::Matrix operator/(const tinymesh::Matrix &m, double s);

#endif  // TINYMESH_MATRIX_H
