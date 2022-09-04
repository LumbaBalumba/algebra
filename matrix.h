//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_MATRIX_H
#define ALGEBRA_MATRIX_H

#include "vec.h"

class matrix {
private:
    size_t _size;
    vec* arr;
public:
    matrix() = default;

    matrix(size_t rows, size_t cols);

    matrix(matrix& other);

    matrix(const complex& lambda, size_t rows, size_t cols);

    size_t rows() const;

    size_t cols() const;

    void resize(size_t rows, size_t cols);

    bool real();

    matrix& operator=(const matrix& other);

    matrix operator+(const matrix& other);

    matrix operator-(const matrix& other);

    matrix operator-();

    matrix operator*(const matrix& other);

    matrix operator*(const complex& z);

    matrix operator/(const complex& z);

    matrix& operator+=(const matrix& other);

    matrix& operator-=(const matrix& other);

    matrix& operator*=(const matrix& other);

    matrix& operator*=(const complex& z);

    matrix& operator/=(const complex& z);

    vec operator[](size_t index);

    void row_swap(size_t dest, size_t src);

    void row_add(size_t dest, size_t src, const complex& k);

    void row_add(size_t dest, size_t src);

    void row_mul(size_t dest, const complex& z);

    void col_swap(size_t dest, size_t src);

    void col_add(size_t dest, size_t src, const complex& z);

    void col_add(size_t dest, size_t src);

    void col_mul(size_t dest, const complex& z);

    complex tr();

    matrix transposed();

    matrix conjugate();

    matrix& upper_triangle();

    matrix& lower_triangle();

    complex det();

    complex alg_complement(size_t row, size_t col);

    matrix adjusted();

    matrix inverted();

    friend std::ostream& operator<<(std::ostream& out, const matrix& m);

    friend std::istream& operator>>(std::istream& in, matrix& m);

    vec operator()(vec& v);
};


#endif //ALGEBRA_MATRIX_H