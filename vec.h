//
// Created by admin on 04.09.22.
//

#ifndef ALGEBRA_VEC_H
#define ALGEBRA_VEC_H

#include "complex.h"

class vec
{
private:
    size_t _size;
    complex *arr;
public:
    vec() = default;

    explicit vec(size_t _size);

    vec(const vec &other);

    ~vec();

    size_t size() const;

    void resize(size_t size);

    bool real();

    void fill(const complex &z);

    vec operator+(const vec &other);

    vec operator-(const vec &other);

    vec operator-();

    vec operator*(const complex &z);

    vec operator/(const complex &z);

    vec &operator+=(const vec &other);

    vec &operator-=(const vec &other);

    vec &operator*=(const complex &other);

    vec &operator/=(const complex &other);

    vec &operator=(const vec &other);

    bool operator==(const vec &other);

    bool operator!=(const vec &other);

    complex &operator[](size_t index) const;

    friend std::ostream &operator<<(std::ostream &out, const vec &v);

    friend std::istream &operator>>(std::istream &, vec &v);

    double length();
};


#endif //ALGEBRA_VEC_H
