//
// Created by admin on 04.09.22.
//

#include "complex.h"

const double eps = 0.00001;

double pi() {
    return acos(-1);
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
    if(abs() < eps) throw std::overflow_error("Zero division error");
    else if(re > 0 && im == 0) return 0;
    else if(re > 0 && im > 0) return atan(im / re);
    else if(re == 0 && im > 0) return pi() / 2;
    else if(re < 0 && im > 0) return pi() + atan(im / re);
    else if(re < 0 && im == 0) return pi();
    else if(re < 0 && im < 0) return pi() + atan(im / re);
    else if(re == 0 && im < 0) return pi() * 3.0 / 2.0;
    else if(re > 0 && im < 0) return pi() * 2.0 + atan(re / im);
    else return pi() / 4;
}

complex::complex(double x) : re(x), im(0) {}

complex::complex(double re, double im) : re(re), im(im) {}

complex::complex(const complex& other) = default;

bool complex::real() const {
    return fabs(im) < eps;
}

complex complex::operator+(const complex& other) const {
    return {re + other.re, im + other.im};
}

complex complex::operator-(const complex& other) const {
    return {re - other.re, im - other.im};
}

complex complex::operator-() const {
    return {-re, -im};
}

complex complex::operator*(const complex& other) const {
    return {re * other.re - im * other.im, im * other.re + re * other.im};
}

complex complex::operator/(const complex& other) {
    if(other.abs() < eps) throw std::overflow_error("Zero division error");
    complex res = (*this) * other.conjugate();
    res.re /= other.abs();
    res.im /= other.abs();
    return res;
}

complex& complex::operator+=(const complex& other) {
    return *this = *this + other;
}

complex& complex::operator-=(const complex& other) {
    return *this = *this - other;
}

complex& complex::operator*=(const complex& other) {
    return *this = *this * other;
}

complex& complex::operator/=(const complex& other) {
    return *this = *this / other;
}

complex& complex::operator=(complex other) {
    if(this == &other) return *this;
    re = other.re;
    im = other.im;
    return *this;
}

complex& complex::operator=(double other) {
    *this = complex(other);
}

bool complex::operator==(const complex& other) const {
    return (fabs(re - other.re) < eps) && (fabs(im - other.im) < eps);
}

bool complex::operator!=(const complex& other) const {
    return (fabs(re - other.re) >= eps) || (fabs(im - other.im) >= eps);
}

std::ostream& operator<<(std::ostream& out, const complex& z) {
    if(fabs(z.re) < eps) {
        out << z.im << "i";
    } else if(fabs(z.im) < eps) {
        out << z.re;
    } else {
        out << z.re << " + " << z.im << "i";
    }
    return out;
}

std::istream& operator>>(std::istream& in, complex& z) {
    in >> z.re;
    in >> z.im;
    return in;
}