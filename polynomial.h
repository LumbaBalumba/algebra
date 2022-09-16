//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_POLYNOMIAL_H
#define ALGEBRA_POLYNOMIAL_H

#include "vec.h"

class polynomial : public vec {
public:
    using vec::operator/;
    using vec::operator*;

    polynomial() = default;

    polynomial(size_t size);

    polynomial(const polynomial& other);

    polynomial(const vec& other);

    polynomial operator+(const polynomial& other);

    polynomial operator-(const polynomial& other);

    polynomial operator*(const polynomial& other);

    polynomial operator/(const polynomial& other);

    polynomial operator%(const polynomial& other);

    size_t deg() const;

    polynomial operator*=(const polynomial& other);

    polynomial operator/=(const polynomial& other);


    polynomial operator%=(const polynomial& other);

    complex* roots();

    polynomial derivative() const;

    friend std::istream& operator>>(std::istream& in, polynomial& p);

    friend std::ostream& operator<<(std::ostream& out, const polynomial& p);

    complex operator()(const complex& z) const;
};


#endif //ALGEBRA_POLYNOMIAL_H
