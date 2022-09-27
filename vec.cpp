//
// Created by admin on 04.09.22.
//

#include "vec.h"

vec::vec(size_t _size) : _size(_size), arr(new complex[_size]) {}

vec::vec(const vec &other) : _size(other._size), arr(new complex[_size]) {
    for (size_t j = 0; j < _size; ++j) {
        arr[j] = other.arr[j];
    }
}

vec::~vec() {
    delete[] arr;
}

size_t vec::size() const {
    return _size;
}

void vec::resize(size_t size) {
    arr = (complex *) realloc(arr, sizeof(complex) * size);
    _size = size;
}

bool vec::real() {
    for (size_t i = 0; i < _size; ++i) {
        if (!arr[i].real()) return false;
    }
    return true;
};

void vec::fill(const complex &z) {
    for (size_t j = 0; j < _size; ++j) {
        arr[j] = z;
    }
}


vec vec::operator+(const vec &other) {
    if (_size != other._size) {
        throw std::out_of_range("Incorrect vector size");;
    }
    vec res(other);
    for (size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] + other.arr[j];
    }
    return res;
}

vec vec::operator-(const vec &other) {
    if (_size != other._size) {
        throw std::out_of_range("Incorrect vector size");
    }
    vec res(other);
    for (int j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] + other.arr[j];
    }
    return res;
}

vec vec::operator-() {
    vec res(_size);
    for (size_t i = 0; i < _size; ++i) {
        res.arr[i] = -arr[i];
    }
    return res;
}

vec vec::operator*(const complex &z) {
    vec res(*this);
    for (size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] * z;
    }
    return res;
}

vec vec::operator/(const complex &z) {
    vec res(*this);
    for (size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] / z;
    }
    return res;
}

vec &vec::operator=(const vec &other) {
    if (this == &other) return *this;
    _size = other._size;
    arr = new complex[_size];
    for (size_t j = 0; j < _size; ++j) {
        arr[j] = other.arr[j];
    }
    return *this;
}

vec &vec::operator+=(const vec &other) {
    return *this = *this + other;
}

vec &vec::operator-=(const vec &other) {
    return *this = *this - other;
}

vec &vec::operator*=(const complex &other) {
    return *this = *this * other;
}

vec &vec::operator/=(const complex &other) {
    return *this = *this / other;
}

bool vec::operator==(const vec &other) {
    for (size_t j = 0; j < _size; ++j) {
        if (arr[j] != other.arr[j]) return false;
    }
    return true;
}

bool vec::operator!=(const vec &other) {
    for (int j = 0; j < _size; ++j) {
        if (arr[j] == other.arr[j]) return false;
    }
    return true;
}

complex &vec::operator[](size_t index) const {
    return arr[index];
}

std::ostream &operator<<(std::ostream &out, const vec &v) {
    for (size_t i = 0; i < v.size(); ++i) {
        out << v[i] << " ";
    }
    return out;
}

std::istream &operator>>(std::istream &in, vec &v) {
    for (size_t i = 0; i < v.size(); ++i) {
        in >> v.arr[i];
    }
    return in;
}

double vec::length() {
    if (!real()) return 0;
    double res = 0;
    for (int i = 0; i < _size; ++i) {
        res += arr[i].re * arr[i].re;
    }
    return sqrt(res);
}