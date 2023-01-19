//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_MATRIX_HPP
#define ALGEBRA_MATRIX_HPP

#include "polynomial.hpp"
#include "vec.hpp"

class matrix {
private:
    size_t _size;
    vec *arr;

public:
    matrix() = default;

    matrix(size_t rows, size_t cols);

    matrix(const matrix &other);

    matrix(const complex &lambda, size_t rows, size_t cols);

    size_t rows() const;

    size_t cols() const;

    void resize(size_t rows, size_t cols);

    bool real();

    matrix &operator=(const matrix &other);

    matrix operator+(const matrix &other) const;

    matrix operator-(const matrix &other) const;

    matrix operator-() const;

    matrix operator*(const matrix &other) const;

    matrix operator*(const complex &z) const;

    matrix operator/(const complex &z) const;

    matrix &operator+=(const matrix &other);

    matrix &operator-=(const matrix &other);

    matrix &operator*=(const matrix &other);

    matrix &operator*=(const complex &z);

    matrix &operator/=(const complex &z);

    vec &operator[](size_t index) const;

    void row_swap(size_t dest, size_t src);

    void row_add(size_t dest, size_t src, const complex &k);

    void row_add(size_t dest, size_t src);

    void row_mul(size_t dest, const complex &z);

    void col_swap(size_t dest, size_t src);

    void col_add(size_t dest, size_t src, const complex &z);

    void col_add(size_t dest, size_t src);

    void col_mul(size_t dest, const complex &z);

    complex tr();

    matrix transposed();

    matrix conjugate();

    matrix upper_triangle();

    complex det();

    complex alg_complement(size_t row, size_t col);

    matrix adjusted();

    matrix inverted();

    friend std::ostream &operator<<(std::ostream &out, const matrix &m);

    friend std::istream &operator>>(std::istream &in, matrix &m);

    vec operator()(vec &v) const;

    size_t rank();

    size_t def();

    [[nodiscard]] polynomial char_pol() const;

    [[nodiscard]] std::vector<complex> eigenvalues() const;

    [[nodiscard]] bool nilpotent() const;

    [[nodiscard]] matrix jordan_form() const;
};

#endif // ALGEBRA_MATRIX_HPP
