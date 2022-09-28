//
// Created by admin on 04.09.22.
//

#include "polynomial.h"

polynomial::polynomial(size_t size)
        : vec(size + 1) {
    fill(0);
}

polynomial::polynomial(const polynomial &other)
        : vec(other) {
}

polynomial::polynomial(const vec &other)
        : vec(other) {
}

polynomial::polynomial(size_t size, complex arr[]) : vec(size) {
    for (int i = 0; i < size; ++i) {
        (*this)[i] = arr[i];
    }
}

size_t polynomial::deg() const {
    size_t res = 0;
    for (int i = 0; i < size(); ++i) {
        if ((*this)[i].abs() >= eps)
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
        if (p[i].abs() > eps) {
            if (p[i] != 1)
                out << " + (" << p[i] << ")x";
            else out << " + x";
        }
        if (i != 1) out << "^" << i;
    }
    return out;
}

complex polynomial::operator()(const complex &z) const {
    complex res(0);
    for (size_t i = 0; i <= deg(); ++i) {
        if ((*this)[i] != 0) {
            complex tmp(1);
            for (size_t j = 0; j < i; ++j) {
                tmp *= z;
            }
            tmp *= (*this)[i];
            res += tmp;
        }
    }
    return res;
}

std::vector<complex> polynomial::roots() const {
    std::vector<complex> res(deg());
    polynomial tmp = (*this);
    bool flag = true;
    while (res.size() < deg() && flag) {
        switch (deg()) {
            case 0:
                flag = false;
                break;
            case 1: {
                res.push_back(-tmp[0] / tmp[1]);
                flag = false;
                break;
            }
            case 2: {
                complex D = (tmp[1].pow(2) - tmp[0] * tmp[2] * 4.0).pow(0.5);
                res.push_back((-tmp[1] + D) / 2 / tmp[0]);
                res.push_back((-tmp[1] - D) / 2 / tmp[0]);
                flag = false;
                break;
            }
            case 3: {
                complex p = (tmp[3] * tmp[1] * 3.0 - tmp[2].pow(2)) / (tmp[3].pow(2) * 3);
                complex q = (tmp[2].pow(3) * 2 - tmp[3] * tmp[2] * tmp[1] * 9 +
                             tmp[3].pow(2) * tmp[0]) * 27 / (tmp[3].pow(3) * 27);
                complex Q = (p / 3).pow(3) + (q / 2).pow(2);
                complex alpha = (-q / 2 + Q.pow(0.5)).pow(1.0 / 3.0);
                complex beta = (-q / 2 - Q.pow(0.5)).pow(1.0 / 3.0);
                res.push_back(alpha + beta);
                res.push_back(-(alpha + beta) / 2 + i() * (alpha - beta) / 2 * complex(3).pow(0.5));
                res.push_back(-(alpha + beta) / 2 - i() * (alpha - beta) / 2 * complex(3).pow(0.5));
                flag = false;
                break;
            }
            case 4: {
                complex A = tmp[4], B = tmp[3], C = tmp[2], D = tmp[1], E = tmp[0];
                complex alpha = -(complex(3) * B.pow(2)) / (complex(8) * A.pow(2)) + C / A;
                complex beta = B.pow(3) / (complex(8) * A.pow(3)) - B * C / (complex(2) * A.pow(2)) + D / A;
                complex gamma =
                        -(complex(3) * B.pow(4)) / (complex(256) * A.pow(4)) + B.pow(2) * C / (complex(16) * A.pow(3)) -
                        B * D / (complex(4) * A.pow(2)) + E / A;
                if (beta == complex(0)) {
                    res.push_back(-B / (complex(4) * A) +
                                  ((-alpha + (alpha.pow(2) - complex(4) * gamma).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / (complex(4) * A) +
                                  ((-alpha - (alpha.pow(2) - complex(4) * gamma).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / (complex(4) * A) -
                                  ((-alpha + (alpha.pow(2) - complex(4) * gamma).pow(0.5)) / 2).pow(0.5));
                    res.push_back(-B / (complex(4) * A) -
                                  ((-alpha - (alpha.pow(2) - complex(4) * gamma).pow(0.5)) / 2).pow(0.5));
                } else {
                    complex P = -alpha.pow(2) / complex(12) - gamma;
                    complex Q = -alpha.pow(3) / complex(108) + alpha * gamma / complex(3) - beta.pow(2) / complex(8);
                    complex R = -Q / 2 + (Q.pow(2) / complex(4) + P.pow(3) / complex(27)).pow(2);
                    complex U = R.pow(1.0 / 3.0);
                    complex y{};
                    if (U == complex(0)) {
                        y = -complex(5.0 / 6.0) * alpha + U - Q.pow(1.0 / 3.0);
                    } else {
                        y = -complex(5.0 / 6.0) * alpha + U - P / (complex(3) * U);
                    }
                    complex W = (alpha + y * complex(2));
                    res.push_back(-B / (complex(4) * A) +
                                  (W + (-(complex(3) * alpha + complex(2) * y + complex(2) * beta / W)).pow(0.5)) /
                                  complex(2));
                    res.push_back(-B / (complex(4) * A) +
                                  (-W + (-(complex(3) * alpha + complex(2) * y - complex(2) * beta / W)).pow(0.5)) /
                                  complex(2));
                    res.push_back(-B / (complex(4) * A) +
                                  (W - (-(complex(3) * alpha + complex(2) * y + complex(2) * beta / W)).pow(0.5)) /
                                  complex(2));
                    res.push_back(-B / (complex(4) * A) +
                                  (-W - (-(complex(3) * alpha + complex(2) * y - complex(2) * beta / W)).pow(0.5)) /
                                  complex(2));
                }
                flag = false;
                break;
            }
            default: {
                const size_t sz = 7;
                complex common_roots[sz] = {0, 1, -1, i(), -i(), 2, -2};
                for (auto &common_root: common_roots) {
                    if (tmp(common_root) == 0) {
                        polynomial d(1);
                        d[0] = 1;
                        d[1] = -common_root;
                        res.push_back(common_root);
                        tmp /= d;
                    }
                }
                //IDK LOL
            }
        }
    }
    return res;
}

