//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_POLYNOMIAL_H
#define ALGEBRA_POLYNOMIAL_H

#include "vec.h"
#include <vector>
#include <algorithm>

class polynomial : public vec
{
public:
    using vec::operator/;
    using vec::operator*;
    using vec::operator*=;
    using vec::operator/=;
    using vec::real;

    polynomial() = default;

    polynomial(size_t size);

    polynomial(const polynomial &other);

    polynomial(const vec &other);

    polynomial(size_t size, complex arr[]);

    polynomial operator+(const polynomial &other);

    polynomial operator-(const polynomial &other);

    polynomial operator*(const polynomial &other);

    polynomial operator/(const polynomial &other);

    polynomial operator%(const polynomial &other);

    size_t deg() const;

    polynomial operator*=(const polynomial &other);

    polynomial operator/=(const polynomial &other);

    polynomial operator%=(const polynomial &other);

    std::vector<complex> roots() const;

    polynomial derivative() const;

    friend std::istream &operator>>(std::istream &in, polynomial &p);

    friend std::ostream &operator<<(std::ostream &out, const polynomial &p);

    complex operator()(const complex &z) const;
};

polynomial Lagrange(const std::vector<std::pair<complex, complex>> &v);

#endif //ALGEBRA_POLYNOMIAL_H
