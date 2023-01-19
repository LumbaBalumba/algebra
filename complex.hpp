//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_COMPLEX_HPP
#define ALGEBRA_COMPLEX_HPP

#include <cmath>
#include <iostream>

class complex {
public:
    double re, im;

    complex() = default;

    complex(double x);

    complex(double re, double im);

    complex(const complex &other);

    [[nodiscard]] bool real() const;

    [[nodiscard]] bool imaginary() const;

    [[nodiscard]] complex conjugate() const;

    [[nodiscard]] double abs() const;

    [[nodiscard]] double arg() const;

    complex operator+(const complex &other) const;

    complex operator-(const complex &other) const;

    complex operator-() const;

    complex operator*(double x) const;

    complex operator/(double x) const;

    complex operator*(const complex &other) const;

    complex operator/(const complex &other) const;

    complex &operator+=(const complex &other);

    complex &operator-=(const complex &other);

    complex &operator*=(const complex &other);

    complex &operator/=(const complex &other);

    bool operator==(const complex &other) const;

    bool operator!=(const complex &other) const;

    complex &operator=(complex other);

    complex &operator=(double other);

    friend std::ostream &operator<<(std::ostream &out, const complex &z);

    friend std::istream &operator>>(std::istream &in, complex &z);

    [[nodiscard]] complex pow(double x) const;
};

complex i();

double pi();

complex sin(const complex &z);

complex cos(const complex &z);

complex exp(const complex &z);

#endif // ALGEBRA_COMPLEX_HPP
