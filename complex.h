//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_COMPLEX_H
#define ALGEBRA_COMPLEX_H

#include <iostream>
#include <cmath>



class complex
{
public:
    double re, im;

    complex() = default;

    complex(double x);

    complex(double re, double im);

    complex(const complex &other);

    bool real() const;

    bool imaginary() const;

    complex conjugate() const;

    double abs() const;

    double arg() const;

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

    complex pow(double x) const;
};


complex i();

double pi();

complex sin(const complex &z);

complex cos(const complex &z);

complex exp(const complex &z);


#endif //ALGEBRA_COMPLEX_H
