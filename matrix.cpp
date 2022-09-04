//
// Created by admin on 04.09.22.
//

#include "matrix.h"

matrix::matrix(size_t rows, size_t cols) : _size(rows), arr(new vec[cols]) {}

matrix::matrix(matrix& other) : _size(other.rows()), arr(new vec[other.cols()]) {
    for(int i = 0; i < _size; ++i) {
        arr[i] = other.arr[i];
    }
}

matrix::matrix(const complex& lambda, size_t rows, size_t cols) : _size(rows), arr(new vec[cols]) {
    for(size_t j = 0; j < rows; ++j) {
        arr[j].fill(complex(0));
        if(j < cols) arr[j][j] = lambda;
    }
}

size_t matrix::rows() const {
    return _size;
}

size_t matrix::cols() const {
    return arr[0].size();
}

bool matrix::real() {
    for(size_t i = 0; i < _size; ++i) {
        if(!arr[i].real()) return false;
    }
    return true;
}

matrix matrix::operator+(const matrix& other) {
    if(this->rows() != other.rows() || this->cols() != other.cols()) throw std::out_of_range("Incorrect matrix size");
    matrix res(rows(), cols());
    for(size_t j = 0; j < _size; ++j) {
        res[j] = this->arr[j] + other.arr[j];
    }
    return res;
}

matrix matrix::operator-(const matrix& other) {
    if(this->rows() != other.rows() || this->cols() != other.cols()) throw std::out_of_range("Incorrect matrix size");
    matrix res(rows(), cols());
    for(int j = 0; j < _size; ++j) {
        res[j] = this->arr[j] - other.arr[j];
    }
    return res;
}

matrix matrix::operator-() {
    matrix res(*this);
    for(size_t i = 0; i < rows(); ++i) {
        for(size_t j = 0; j < cols(); ++j) {
            res[i][j] = -res[i][j];
        }
    }
    return res;
}

matrix matrix::operator*(const complex& z) {
    matrix res(rows(), cols());
    for(size_t j = 0; j < _size; ++j) {
        res[j] = this->arr[j] * z;
    }
    return res;
}

matrix matrix::operator/(const complex& z) {
    matrix res(rows(), cols());
    for(size_t j = 0; j < _size; ++j) {
        res[j] = this->arr[j] / z;
    }
    return res;
}

matrix matrix::operator*(const matrix& other) {
    if(cols() != other.rows()) throw std::out_of_range("Incorrect matrix size");
    matrix res(rows(), other.cols());
    for(size_t i = 0; i < rows(); ++i) {
        for(size_t j = 0; j < other.cols(); ++j) {
            complex z(0);
            for(int k = 0; k < cols(); ++k) {
                z += this->arr[i][k] * other.arr[k][j];
            }
            res.arr[i][j] = z;
        }
    }
    return res;
}

matrix& matrix::operator=(const matrix& other) {
    if(this == &other) return *this;
    _size = other._size;
    arr = new vec[_size];
    for(size_t i = 0; i < _size; ++i) {
        arr[i] = other.arr[i];
    }
    return *this;
}

matrix& matrix::operator+=(const matrix& other) {
    return *this = *this + other;
}

matrix& matrix::operator-=(const matrix& other) {
    return *this = *this - other;
}

matrix& matrix::operator*=(const matrix& other) {
    return *this = *this * other;
}

matrix& matrix::operator*=(const complex& z) {
    return *this = *this * z;
}

matrix& matrix::operator/=(const complex& z) {
    return *this = *this / z;
}

vec matrix::operator[](size_t index) {
    return arr[index];
}

void matrix::row_swap(size_t dest, size_t src) {
    vec tmp(arr[dest]);
    arr[dest] = arr[src];
    arr[src] = tmp;
}

void matrix::row_add(size_t dest, size_t src, const complex& k) {
    arr[dest] += arr[src] * k;
}

void matrix::row_add(size_t dest, size_t src) {
    arr[dest] += arr[src];
}

void matrix::row_mul(size_t dest, const complex& z) {
    arr[dest] *= z;
}

void matrix::col_swap(size_t dest, size_t src) {
    for(size_t j = 0; j < rows(); ++j) {
        complex z = arr[j][dest];
        arr[j][dest] = arr[j][src];
        arr[j][src] = z;
    }
}

void matrix::col_add(size_t dest, size_t src, const complex& z) {
    for(size_t j = 0; j < rows(); ++j) {
        arr[j][dest] += arr[j][src] * z;
    }
}

void matrix::col_add(size_t dest, size_t src) {
    for(size_t j = 0; j < rows(); ++j) {
        arr[j][dest] += arr[j][src];
    }
}

void matrix::col_mul(size_t dest, const complex& z) {
    for(size_t j = 0; j < rows(); ++j) {
        arr[j][dest] *= z;
    }
}

complex matrix::tr() {
    if(rows() != cols()) throw std::out_of_range("Incorrect matrix size");
    complex res(0);
    for(size_t j = 0; j < rows(); ++j) {
        res += arr[j][j];
    }
    return res;
}

matrix matrix::transposed() {
    matrix res(cols(), rows());
    for(size_t i = 0; i < rows(); ++i) {
        for(size_t j = 0; j < cols(); ++j) {
            res.arr[j][i] = arr[i][j];
        }
    }
    return res;
}

matrix matrix::conjugate() {
    matrix res(cols(), rows());
    for(size_t i = 0; i < rows(); ++i) {
        for(size_t j = 0; j < cols(); ++j) {
            res.arr[j][i] = arr[i][j].conjugate();
        }
    }
    return res;
}

matrix& matrix::upper_triangle() {
    matrix res(*this);
    for(size_t i = 0; i < rows() - 1; ++i) {
        if(res[i][i] == complex(0)) {
            bool flag = false;
            for(size_t j = i + 1; j < rows() && !flag; ++j) {
                if(res[i][j] != complex(0)) {
                    res.row_swap(i, j);
                    flag = true;
                }
            }
            if(!flag) throw std::exception();
        }
        for(size_t j = i + 1; j < rows(); ++j) {
            res.row_add(j, i, complex(-1) / res[i][i] * res[i][j]);
        }
    }
}

matrix& matrix::lower_triangle() {
    matrix res(*this);
    for(size_t i = rows(); i >= 0; --i) {
        if(res[i][i] == complex(0)) {
            bool flag = false;
            for(size_t j = i - 1; j >= 0 && !flag; --j) {
                if(res[i][j] != complex(0)) {
                    res.row_swap(i, j);
                    flag = true;
                }
            }
            if(!flag) throw std::exception();
        }
        for(size_t j = i - 1; j >= 0; --j) {
            res.row_add(j, i, complex(-1) / res[i][i] * res[i][j]);
        }
    }
}

complex matrix::det() {
    if(rows() != cols()) throw std::out_of_range("Incorrect matrix size");
    complex res(1);
    matrix tmp = this->upper_triangle();
    for(size_t i = 0; i < rows(); ++i) {
        res *= tmp[i][i];
    }
    return res;
}

complex matrix::alg_complement(size_t row, size_t col) {
    if(rows() != cols()) throw std::out_of_range("Incorrect matrix size");
    matrix res(*this);
    for(size_t i = 0; i < rows(); ++i) {
        for(size_t j = 0; j < cols(); ++j) {
            if(i < row && j < col) res[i][j] = (*this)[i][j];
            else if(i < row && j > col) res[i][j - 1] = (*this)[i][j];
            else if(i > row && j < col) res[i - 1][j] = (*this)[i][j];
            else if(i > row && j > col) res[i - 1][j - 1] = (*this)[i][j];
        }
    }
    return (row + col) % 2 == 0 ? res.det() : -res.det();
}

matrix matrix::adjusted() {
    if(rows() != cols()) throw std::out_of_range("Incorrect matrix size");
    matrix res(cols(), rows());
    for(size_t i = 0; i < cols(); ++i) {
        for(size_t j = 0; j < rows(); ++j) {
            res[i][j] = alg_complement(j, i);
        }
    }
    return res;
}

matrix matrix::inverted() {
    if(rows() != cols()) throw std::exception();
    return adjusted() / det();
}

std::ostream& operator<<(std::ostream& out, const matrix& m) {
    for(size_t i = 0; i < m.rows(); ++i) {
        out << m.arr[i] << std::endl;
    }
    return out;
}

std::istream& operator>>(std::istream& in, matrix& m) {
    for(size_t i = 0; i < m.rows(); ++i) {
        in >> m.arr[i];
    }
    return in;
}

vec matrix::operator()(vec& v) {
    if(v.size() != cols()) throw std::out_of_range("Incorrect matrix and vector size");
    vec res(rows());
    for(size_t i = 0; i < rows(); ++i) {
        res[i] = 0;
        for(size_t j = 0; j < cols(); ++j) {
            res[i] += v[j] * (*this)[i][j];
        }
    }
    return res;
}

