//
// Created by admin on 04.09.22.
//

#include "polynomial.h"

polynomial::polynomial(size_t size)
        : vec(size + 1) {
    fill(complex(0));
}

polynomial::polynomial(const polynomial& other)
        : vec(other) {
}

polynomial::polynomial(const vec& other)
        : vec(other) {
}

size_t polynomial::deg() const {
    size_t res = 0;
    for(int i = 0; i < size(); ++i) {
        if((*this)[i].abs() >= eps)
            res = i;
    }
    return res;
}

polynomial polynomial::operator+(const polynomial& other) {
    if(deg() >= other.deg()) {
        polynomial res(deg());
        for(size_t i = 0; i <= deg(); ++i) {
            if(i <= other.deg())
                res[i] = (*this)[i] + other[i];
            else
                res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for(size_t i = 0; i <= other.deg(); ++i) {
            if(i <= deg())
                res[i] = (*this)[i] + other[i];
            else
                res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator-(const polynomial& other) {
    if(deg() >= other.deg()) {
        polynomial res(deg());
        for(size_t i = 0; i < deg(); ++i) {
            if(i <= other.deg())
                res[i] = (*this)[i] - other[i];
            else
                res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for(size_t i = 0; i < other.deg(); ++i) {
            if(i <= deg())
                res[i] = (*this)[i] - other[i];
            else
                res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator*(const polynomial& other) {
    polynomial res(deg() + other.deg());
    res.fill(complex(0));
    for(size_t i = 0; i <= deg(); ++i) {
        for(size_t j = 0; j <= other.deg(); ++j) {
            res[i + j] += (*this)[i] * other[j];
        }
    }
    return res;
}

polynomial polynomial::operator/(const polynomial& other) {
    if(deg() < other.deg()) {
        polynomial res(0);
        res[0] = 0;
        return res;
    } else {
        polynomial res(deg() - other.deg()), r = *this;
        res.fill(complex(0));
        size_t n = deg(), m = other.deg();
        for(size_t i = n; i >= m; --i) {
            res[i - m] = r[i] / other[m];
            for(size_t j = m + 1; j > 0; --j) {
                r[i - m + j - 1] -= other[j - 1] * res[i - m];
            }
        }
        return res;
    }
}

polynomial polynomial::operator%(const polynomial& other) {
    return (*this) - (*this) / other * other;
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
    for(size_t i = 0; i < p.size(); ++i) {
        in >> p[i];
    }
    return in;
}

std::ostream& operator<<(std::ostream& out, const polynomial& p) {
    out << p[0];
    for(size_t i = 1; i <= p.deg(); ++i) {
        if(p[i].abs() > eps) {
            if(p[i] != complex(1))
                out << " + (" << p[i] << ")x";
            else out << " + x";
        }
        if(i != 1) out << "^" << i;
    }
    return out;
}

complex polynomial::operator()(const complex& z) const {
    complex res(0);
    for(size_t i = 0; i <= deg(); ++i) {
        if((*this)[i] != complex(0)) {
            complex tmp(1);
            for(size_t j = 0; j < i; ++j) {
                tmp *= z;
            }
            tmp *= (*this)[i];
            res += tmp;
        }
    }
    return res;
}

complex* polynomial::roots() const{
    auto* res = new complex[deg()];

    return res;
}

