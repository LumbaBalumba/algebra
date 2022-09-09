//
// Created by admin on 04.09.22.
//

#include "polynomial.h"

polynomial::polynomial(size_t size) : vec(size + 1) {}

polynomial::polynomial(const polynomial& other) : vec(other) {}

polynomial::polynomial(const vec& other) : vec(other) {}

size_t polynomial::deg() const {
    size_t res = 0;
    for(int i = 0; i < size; ++i) {
        if(arr[i].abs() >= eps) res = i;
    }
    return res;
}

polynomial polynomial::operator+(const polynomial& other) {
    if(deg() >= other.deg()) {
        polynomial res(deg());
        for(size_t i = 0; i < deg(); ++i) {
            if(i <= other.deg()) res[i] = (*this)[i] + other[i];
            else res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for(size_t i = 0; i < other.deg(); ++i) {
            if(i <= deg()) res[i] = (*this)[i] + other[i];
            else res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator-(const polynomial& other) {
    if(deg() >= other.deg()) {
        polynomial res(deg());
        for(size_t i = 0; i < deg(); ++i) {
            if(i <= other.deg()) res[i] = (*this)[i] - other[i];
            else res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for(size_t i = 0; i < other.deg(); ++i) {
            if(i <= deg()) res[i] = (*this)[i] - other[i];
            else res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator*=(const polynomial& other) {
    return *this = *this * other;
}

polynomial polynomial::operator/=(const polynomial& other) {
    return *this = *this / other;
}

polynomial polynomial::operator%=(const polynomial& other) {
    return *this = *this % other;
}

polynomial polynomial::derivative() const {
    if(deg() == 0) {
        polynomial res(1);
        res[0] = 0;
        return res;
    }
    polynomial res((vec(deg() - 1)));
    for(size_t i = 0; i < deg() - 1; ++i) {
        res[i] = (*this)[i + 1] * complex(i + 1);
    }
    return res;
}

std::istream& operator>>(std::istream& in, polynomial& p) {
    for(size_t i = 0; i < p.size; ++i) {
        in >> p.arr[i];
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, polynomial& p) {
    for(size_t i = p.deg(); i > 0; --i) {
        out << "(" << p[i] << ")x^" << i << " + ";
    }
    out << p[0];
    return out;
}








