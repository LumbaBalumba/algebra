//
// Created by admin on 04.09.22.
//

#include "polynomial.hpp"

polynomial::polynomial(size_t size) : vec(size + 1) {
    fill(0);
}

polynomial::polynomial(const polynomial &other) : vec(other) {
}

polynomial::polynomial(const vec &other) : vec(other) {
}

polynomial::polynomial(size_t size, complex arr[]) : vec(size) {
    for (int i = 0; i < size; ++i) {
        (*this)[i] = arr[i];
    }
}

size_t polynomial::deg() const {
    size_t res = 0;
    for (int i = 0; i < size(); ++i) {
        if ((*this)[i] != 0)
            res = i;
    }
    return res;
}

polynomial polynomial::operator+(const polynomial &other) {
    if (deg() >= other.deg()) {
        polynomial res(deg());
        for (size_t i = 0; i <= deg(); ++i) {
            if (i <= other.deg())
                res[i] = (*this)[i] + other[i];
            else
                res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for (size_t i = 0; i <= other.deg(); ++i) {
            if (i <= deg())
                res[i] = (*this)[i] + other[i];
            else
                res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator-(const polynomial &other) {
    if (deg() >= other.deg()) {
        polynomial res(deg());
        for (size_t i = 0; i < deg(); ++i) {
            if (i <= other.deg())
                res[i] = (*this)[i] - other[i];
            else
                res[i] = (*this)[i];
        }
        return res;
    } else {
        polynomial res(other.deg());
        for (size_t i = 0; i < other.deg(); ++i) {
            if (i <= deg())
                res[i] = (*this)[i] - other[i];
            else
                res[i] = other[i];
        }
        return res;
    }
}

polynomial polynomial::operator*(const polynomial &other) {
    polynomial res(deg() + other.deg());
    res.fill(0);
    for (size_t i = 0; i <= deg(); ++i) {
        for (size_t j = 0; j <= other.deg(); ++j) {
            res[i + j] += (*this)[i] * other[j];
        }
    }
    return res;
}

polynomial polynomial::operator/(const polynomial &other) {
    if (deg() < other.deg()) {
        polynomial res(0);
        res[0] = 0;
        return res;
    } else {
        polynomial res(deg() - other.deg()), r = *this;
        res.fill(0);
        size_t n = deg(), m = other.deg();
        for (size_t i = n; i >= m; --i) {
            res[i - m] = r[i] / other[m];
            for (size_t j = m + 1; j > 0; --j) {
                r[i - m + j - 1] -= other[j - 1] * res[i - m];
            }
        }
        return res;
    }
}

polynomial polynomial::operator%(const polynomial &other) {
    return (*this) - (*this) / other * other;
}

polynomial polynomial::operator*=(const polynomial &other) {
    return *this = *this * other;
}

polynomial polynomial::operator/=(const polynomial &other) {
    return *this = *this / other;
}

polynomial polynomial::operator%=(const polynomial &other) {
    return *this = *this % other;
}

polynomial polynomial::derivative() const {
    if (deg() == 0) {
        polynomial res(1);
        res[0] = 0;
        return res;
    }
    polynomial res((vec(deg() - 1)));
    for (size_t i = 0; i < deg() - 1; ++i) {
        res[i] = (*this)[i + 1] * ((double) i + 1);
    }
    return res;
}

std::istream &operator>>(std::istream &in, polynomial &p) {
    for (size_t i = 0; i < p.size(); ++i) {
        in >> p[i];
    }
    return in;
}

std::ostream &operator<<(std::ostream &out, const polynomial &p) {
    out << p[0];
    for (size_t i = 1; i <= p.deg(); ++i) {
        if (p[i] != 0) {
            if (p[i] != 1) {
                out << " + (" << p[i] << ")x";
            } else {
                out << " + x";
            }
        }
        if (i != 1)
            out << "^" << i;
    }
    return out;
}

complex polynomial::operator()(const complex &z) const {
    complex res(0);
    for (size_t i = 0; i <= deg(); ++i) {
        complex tmp(1);
        for (size_t j = 0; j < i; ++j) {
            tmp *= z;
        }
        tmp *= (*this)[i];
        res += tmp;
    }
    return res;
}

std::vector<complex> polynomial::roots() const {
    std::vector<complex> res;
    polynomial tmp = (*this);
    bool flag = true;
    while (tmp.deg() > 0 && flag) {
        const size_t sz = 7;
        complex common_roots[sz] = {0, 1, -1, i(), -i(), 2, -2};
        for (size_t i = 0; i < sz;) {
            if (tmp(common_roots[i]) == 0) {
                polynomial d(1);
                d[0] = -common_roots[i];
                d[1] = 1;
                res.push_back(common_roots[i]);
                tmp /= d;
                // std::cout << tmp << std::endl;
            } else {
                ++i;
            }
        }
        switch (tmp.deg()) {
            case 0:
                flag = false;
                break;
            case 1: {
                res.push_back(-tmp[0] / tmp[1]);
                flag = false;
                break;
            }
            case 2: {
                complex a = tmp[2], b = tmp[1], c = tmp[0];
                complex D = b.pow(2) - a * c * 4;
                res.push_back((-b + D.pow(0.5)) / (a * 2));
                res.push_back((-b - D.pow(0.5)) / (a * 2));
                flag = false;
                break;
            }
            case 3: {
                complex a = tmp[3], b = tmp[2], c = tmp[1], d = tmp[0];
                complex p = (a * c * 3 - b * b) / (a * a * 3);
                complex q = (b * b * b * 2 - a * b * c * 9 + a * a * d * 27) / (a * a * a * 27);
                complex Q = (p / 3).pow(3) + (q / 2).pow(2);
                complex alpha = (-q / 2 + Q.pow(0.5)).pow(1.0 / 3.0);
                complex beta = (-q / 2 - Q.pow(0.5)).pow(1.0 / 3.0);
                res.push_back(alpha + beta - b / (a * 3));
                res.push_back(-(alpha + beta) / 2 + i() * (alpha - beta) / 2 * sqrt(3) - b / (a * 3));
                res.push_back(-(alpha + beta) / 2 - i() * (alpha - beta) / 2 * sqrt(3) - b / (a * 3));
                flag = false;
                break;
            }
            case 4: {
                complex A = tmp[4], B = tmp[3], C = tmp[2], D = tmp[1], E = tmp[0];
                complex alpha = -B * B * 3 / 8 / A / A + C / A,
                        beta = B * B * B / 8 / A / A / A - B * C / 2 / A / A + D / A,
                        gamma =
                        -B * B * B * B * 3 / 256 / A / A / A / A + B * B * C / 16 / A / A / A - B * D / 4 / A / A +
                        E / A;
                if (beta == 0) {
                    res.push_back(-B / 4 / A + ((-alpha + (alpha * alpha - gamma * 4).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / 4 / A + ((-alpha - (alpha * alpha - gamma * 4).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / 4 / A - ((-alpha + (alpha * alpha - gamma * 4).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / 4 / A - ((-alpha - (alpha * alpha - gamma * 4).pow(0.5)) / 2).pow(0.5));
                } else {
                    complex P = -alpha * alpha / 12 - gamma;
                    complex Q = -alpha * alpha * alpha / 108 + alpha * gamma / 3 - beta * beta / 8;
                    complex R = -Q / 2 + (Q * Q / 4 + P * P * P / 27).pow(0.5);
                    complex U = R.pow(1.0 / 3.0);
                    complex y = -alpha * 5 / 6 + U;
                    if (U == 0) {
                        y -= Q.pow(1.0 / 3.0);
                    } else {
                        y -= P / 3 / U;
                    }
                    complex W = (alpha + y * 2).pow(0.5);
                    res.push_back(-B / 4 / A + (W + (-alpha * 3 - y * 2 - beta * 2 / W).pow(0.5)) / 2);
                    res.push_back(-B / 4 / A + (-W + (-alpha * 3 - y * 2 + beta * 2 / W).pow(0.5)) / 2);
                    res.push_back(-B / 4 / A + (W - (-alpha * 3 - y * 2 - beta * 2 / W).pow(0.5)) / 2);
                    res.push_back(-B / 4 / A + (-W - (-alpha * 3 - y * 2 + beta * 2 / W).pow(0.5)) / 2);
                }
                flag = false;
                break;
            }
            default: {

                // IDK LOL
            }
        }
    }
    return res;
}

static polynomial l(size_t index, const std::vector<std::pair<complex, complex>> &v) {
    polynomial res(0);
    res[0] = 1;
    for (int i = 0; i < v.size(); ++i) {
        if (i != index) {
            polynomial tmp(1);
            tmp[0] = -v[i].first;
            tmp[1] = 1;
            tmp /= v[index].first - v[i].first;
            res *= tmp;
        }
    }
    return res;
}

polynomial Lagrange(const std::vector<std::pair<complex, complex>> &v) {
    polynomial res(v.size() - 1);
    for (size_t i = 0; i < v.size(); ++i) {
        res += l(i, v) * v[i].second;
    }
    return res;
}
