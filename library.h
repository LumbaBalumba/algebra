#ifndef ALGEBRA_LIBRARY_H
#define ALGEBRA_LIBRARY_H

#include <iostream>

class complex {
public:
    double re, im;

    complex() = default;

    explicit complex(double x);

    complex(double re, double im);

    complex(const complex& other);

    bool real() const;

    complex conjugate() const;

    double abs() const;

    double arg() const;

    complex operator+(const complex& other) const;

    complex operator-(const complex& other) const;

    complex operator-() const;

    complex operator*(const complex& other) const;

    complex operator/(const complex& other);

    complex& operator+=(const complex& other);

    complex& operator-=(const complex& other);

    complex& operator*=(const complex& other);

    complex& operator/=(const complex& other);

    bool operator==(const complex& other) const;

    bool operator!=(const complex& other) const;

    complex& operator=(complex other);

    complex& operator=(double other);

    friend std::ostream& operator<<(std::ostream& out, const complex& z);

    friend std::istream& operator>>(std::istream& in, complex& z);
};


complex i();

double pi();


class vec {
private:
    size_t _size;
    complex* arr;
public:
    vec() = default;

    explicit vec(size_t _size);

    vec(const vec& other);

    ~vec();

    size_t size() const;

    bool real();

    void fill(const complex& z);

    vec operator+(const vec& other);

    vec operator-(const vec& other);

    vec operator-();

    vec operator*(const complex& z);

    vec operator/(const complex& z);

    vec& operator+=(const vec& other);

    vec& operator-=(const vec& other);

    vec& operator*=(const complex& other);

    vec& operator/=(const complex& other);

    vec& operator=(const vec& other);

    bool operator==(const vec& other);

    bool operator!=(const vec& other);

    complex operator[](size_t index) const;

    friend std::ostream& operator<<(std::ostream& out, const vec& v);

    friend std::istream& operator>>(std::istream&, vec& v);
};


class matrix {
private:
    size_t _size;
    vec* arr;
public:
    matrix() = default;

    matrix(size_t rows, size_t cols);

    matrix(const complex& lambda, size_t rows, size_t cols);

    size_t rows() const;

    size_t cols() const;

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

    vec operate(vec& v);
};

#endif //ALGEBRA_LIBRARY_H
