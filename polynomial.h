//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_POLYNOMIAL_H
#define ALGEBRA_POLYNOMIAL_H

#include "vec.h"

class polynomial : public vec {
public:
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

    complex* roots() const;

    polynomial derivative() const;

    friend std::istream& operator>>(std::istream& in, polynomial& p);

    friend std::ostream& operator<<(std::ostream& out, polynomial& p);
};


#endif //ALGEBRA_POLYNOMIAL_H
