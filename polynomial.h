//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_POLYNOMIAL_H
#define ALGEBRA_POLYNOMIAL_H


class polynomial {
    vec v;

    polynomial() = default;

    polynomial(size_t size);

    polynomial(const polynomial& other);

    polynomial(const vec& other);

    polynomial operator=(const polynomial& other);

    polynomial operator+(const polynomial& other);

    polynomial operator*(const polynomial& other);

    polynomial operator-(const polynomial& other);

    polynomial operator-();

    bool operator==(const polynomial& other);

    bool operator!=(const polynomial& other);

    polynomial operator/(const polynomial& other);

    polynomial operator%(const polynomial& other);

    size_t deg();

    polynomial operator*(const complex& other);

    polynomial operator/(const complex& other);

    polynomial operator+=(const polynomial& other);

    polynomial operator-=(const polynomial& other);

    polynomial operator*=(const polynomial& other);

    polynomial operator/=(const polynomial& other);

    polynomial operator%=(const polynomial& other);

    complex* roots();

    polynomial derivative();

    friend std::istream& operator>>(std::istream& in, polynomial& p);

    friend std::ostream& operator<<(std::ostream& out, polynomial& p);
};


#endif //ALGEBRA_POLYNOMIAL_H
