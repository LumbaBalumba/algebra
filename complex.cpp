//
// Created by admin on 04.09.22.
//

#include "complex.hpp"

static const double eps = 0.0000001;

double pi() {
    double res = 1;
    double tmp = sqrt(2);
    for (int i = 0; i < 50; ++i) {
        res *= tmp / 2;
        tmp = sqrt(2 + tmp);
    }
    return 2 / res;
}

complex i() {
    return {0, 1};
}

complex complex::conjugate() const {
    return {re, -im};
}

double complex::abs() const {
    return sqrt(re * re + im * im);
}

double complex::arg() const {
    if (abs() < eps) {
        throw std::overflow_error("Zero division error (arg)\n");
    } else if (re > 0 && fabs(im) < eps)
        return 0;
    else if (re > 0 && im > 0)
        return atan(im / re);
    else if (fabs(re) < eps && im > 0)
        return pi() / 2;
    else if (re < 0 && im > 0)
        return pi() + atan(im / re);
    else if (re < 0 && fabs(im) < eps)
        return pi();
    else if (re < 0 && im < 0)
        return pi() + atan(im / re);
    else if (fabs(re) < eps && im < 0)
        return pi() * 3.0 / 2.0;
    else if (re > 0 && im < 0)
        return pi() * 2.0 + atan(re / im);
    return 0;
}

complex::complex(double x)
        : re(x), im(0) {
}

complex::complex(double re, double im)
        : re(re), im(im) {
}

complex::complex(const complex &other) = default;

bool complex::real() const {
    return fabs(im) < eps;
}

bool complex::imaginary() const {
    return fabs(re) < eps;
}

complex complex::operator+(const complex &other) const {
    return {re + other.re, im + other.im};
}

complex complex::operator-(const complex &other) const {
    return {re - other.re, im - other.im};
}

complex complex::operator-() const {
    return {-re, -im};
}

complex complex::operator*(double x) const {
    return {re * x, im * x};
}

complex complex::operator/(double x) const {
    return {re / x, im / x};
}

complex complex::operator*(const complex &other) const {
    return {re * other.re - im * other.im, im * other.re + re * other.im};
}

complex complex::operator/(const complex &other) const {
    if (other.abs() < eps) {
        throw std::overflow_error("Zero division error (/)\n");
    }
    complex res = (*this) * other.conjugate();
    res.re /= other.abs() * other.abs();
    res.im /= other.abs() * other.abs();
    return res;
}

complex &complex::operator+=(const complex &other) {
    return *this = *this + other;
}

complex &complex::operator-=(const complex &other) {
    return *this = *this - other;
}

complex &complex::operator*=(const complex &other) {
    return *this = *this * other;
}

complex &complex::operator/=(const complex &other) {
    return *this = *this / other;
}

complex &complex::operator=(complex other) {
    if (this == &other)
        return *this;
    re = other.re;
    im = other.im;
    return *this;
}

complex &complex::operator=(double other) {
    return *this = complex(other);
}

bool complex::operator==(const complex &other) const {
    return (fabs(re - other.re) < eps) && (fabs(im - other.im) < eps);
}

bool complex::operator!=(const complex &other) const {
    return (fabs(re - other.re) >= eps) || (fabs(im - other.im) >= eps);
}

std::ostream &operator<<(std::ostream &out, const complex &z) {
    if (fabs(z.im) < eps) {
        out << z.re;
    } else if (fabs(z.re) < eps) {
        if (fabs(z.im - 1.0) < eps)
            out << "i";
        else
            out << z.im << "i";
    } else {
        if (z.im >= 0) {
            out << z.re << " + " << z.im << "i";
        } else {
            out << z.re << " " << z.im << "i";
        }
    }
    return out;
}

std::istream &operator>>(std::istream &in, complex &z) {
    in >> z.re;
    in >> z.im;
    return in;
}

complex complex::pow(double x) const {
    if (*this == 0) {
        return 0;
    } else {
        double phi = arg() * x, r = std::pow(abs(), x);
        return {r * std::cos(phi), r * std::sin(phi)};
    }
}

complex operator+(const double left, const complex &right) {
    return {left + right.re, right.im};
}

complex operator-(double left, const complex &right) {
    return {left - right.re, -right.im};
}

complex operator*(double left, const complex &right) {
    return {left * right.re, left * right.im};
}

complex operator/(double left, const complex &right) {
    return left * right.conjugate() / (right.re * right.re + right.im * right.im);
}

complex sin(const complex &z) {
    complex res(0);
    int sgn = -1;
    complex tmp = z;
    unsigned fact = 1;
    for (int i = 1; i <= 10; ++i) {
        if ((i & 1) == 1) {
            sgn *= -1;
            res += tmp * sgn / fact;
            tmp *= z.pow(2);
            fact *= (i + 1) * (i + 2);
        }
    }
    return res;
}

complex cos(const complex &z) {
    complex res(1);
    int sgn = 1;
    complex tmp = z;
    unsigned fact = 1;
    for (int i = 2; i <= 10; ++i) {
        if ((i & 1) == 0) {
            sgn *= -1;
            res += tmp * sgn / fact;
            tmp *= z.pow(2);
            fact *= (i + 1) * (i + 2);
        }
    }
    return res;
}

complex exp(const complex &z) {
    return complex(std::exp(z.re)) * (complex(cos(z.im)) + i() * sin(z.im));
}
