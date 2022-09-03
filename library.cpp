#include "library.h"

#include <cmath>

const double eps = 0.00001;

double pi() {
    return acos(-1);
}

complex i() {
    return {0, 1};
}

complex complex::conjugate() const {
    return {re, -im};
}

double complex::abs() const {
    return sqrt(re * re + im * im);
}

double complex::arg() const {
    if(abs() < eps) throw std::overflow_error("Zero division error");
    else if(re > 0 && im == 0) return 0;
    else if(re > 0 && im > 0) return atan(im / re);
    else if(re == 0 && im > 0) return pi() / 2;
    else if(re < 0 && im > 0) return pi() + atan(im / re);
    else if(re < 0 && im == 0) return pi();
    else if(re < 0 && im < 0) return pi() + atan(im / re);
    else if(re == 0 && im < 0) return pi() * 3.0 / 2.0;
    else if(re > 0 && im < 0) return pi() * 2.0 + atan(re / im);
    else return pi() / 4;
}

complex::complex(double x) : re(x), im(0) {}

complex::complex(double re, double im) : re(re), im(im) {}

complex::complex(const complex& other) = default;

bool complex::real() const {
    return fabs(im) < eps;
}

complex complex::operator+(const complex& other) const {
    return {re + other.re, im + other.im};
}

complex complex::operator-(const complex& other) const {
    return {re - other.re, im - other.im};
}

complex complex::operator-() const {
    return {-re, -im};
}

complex complex::operator*(const complex& other) const {
    return {re * other.re - im * other.im, im * other.re + re * other.im};
}

complex complex::operator/(const complex& other) {
    if(other.abs() < eps) throw std::overflow_error("Zero division error");
    complex res = (*this) * other.conjugate();
    res.re /= other.abs();
    res.im /= other.abs();
    return res;
}

complex& complex::operator+=(const complex& other) {
    return *this = *this + other;
}

complex& complex::operator-=(const complex& other) {
    return *this = *this - other;
}

complex& complex::operator*=(const complex& other) {
    return *this = *this * other;
}

complex& complex::operator/=(const complex& other) {
    return *this = *this / other;
}

complex& complex::operator=(complex other) {
    if(this == &other) return *this;
    re = other.re;
    im = other.im;
    return *this;
}

complex& complex::operator=(double other) {
    *this = complex(other);
}

bool complex::operator==(const complex& other) const {
    return (fabs(re - other.re) < eps) && (fabs(im - other.im) < eps);
}

bool complex::operator!=(const complex& other) const {
    return (fabs(re - other.re) >= eps) || (fabs(im - other.im) >= eps);
}

std::ostream& operator<<(std::ostream& out, const complex& z) {
    if(fabs(z.re) < eps) {
        out << z.im << "i";
    } else if(fabs(z.im) < eps) {
        out << z.re;
    } else {
        out << z.re << " + " << z.im << "i";
    }
    return out;
}

std::istream& operator>>(std::istream& in, complex& z) {
    in >> z.re;
    in >> z.im;
    return in;
}

vec::vec(size_t _size) : _size(_size), arr(new complex[_size]) {}

vec::vec(const vec& other) : _size(other._size), arr(new complex[_size]) {
    for(size_t j = 0; j < _size; ++j) {
        arr[j] = other.arr[j];
    }
}

vec::~vec() {
    delete[] arr;
}

size_t vec::size() const {
    return _size;
}

bool vec::real() {
    for(size_t i = 0; i < _size; ++i) {
        if(!arr[i].real()) return false;
    }
    return true;
};

void vec::fill(const complex& z) {
    for(size_t j = 0; j < _size; ++j) {
        arr[j] = z;
    }
}


vec vec::operator+(const vec& other) {
    if(_size != other._size) {
        throw std::out_of_range("Incorrect vector size");;
    }
    vec res(other);
    for(size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] + other.arr[j];
    }
}

vec vec::operator-(const vec& other) {
    if(_size != other._size) {
        throw std::out_of_range("Incorrect vector size");
    }
    vec res(other);
    for(int j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] + other.arr[j];
    }
}

vec vec::operator-() {
    vec res(_size);
    for(size_t i = 0; i < _size; ++i) {
        res.arr[i] = -arr[i];
    }
    return res;
}

vec vec::operator*(const complex& z) {
    vec res(*this);
    for(size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] * z;
    }
}

vec vec::operator/(const complex& z) {
    vec res(*this);
    for(size_t j = 0; j < _size; ++j) {
        res.arr[j] = arr[j] / z;
    }
}

vec& vec::operator=(const vec& other) {
    if(this == &other) return *this;
    _size = other._size;
    arr = new complex[_size];
    for(size_t j = 0; j < _size; ++j) {
        arr[j] = other.arr[j];
    }
    return *this;
}

vec& vec::operator+=(const vec& other) {
    return *this = *this + other;
}

vec& vec::operator-=(const vec& other) {
    return *this = *this - other;
}

vec& vec::operator*=(const complex& other) {
    return *this = *this * other;
}

vec& vec::operator/=(const complex& other) {
    return *this = *this / other;
}

bool vec::operator==(const vec& other) {
    for(size_t j = 0; j < _size; ++j) {
        if(arr[j] != other.arr[j]) return false;
    }
    return true;
}

bool vec::operator!=(const vec& other) {
    for(int j = 0; j < _size; ++j) {
        if(arr[j] == other.arr[j]) return false;
    }
    return true;
}

complex vec::operator[](size_t index) const {
    return arr[index];
}

std::ostream& operator<<(std::ostream& out, const vec& v) {
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i] << " ";
    }
    return out;
}

std::istream& operator>>(std::istream& in, vec& v) {
    for(size_t i = 0; i < v.size(); ++i) {
        in >> v.arr[i];
    }
    return in;
}

double vec::length() {
    if(!real()) return 0;
    double res = 0;
    for(int i = 0; i < _size; ++i) {
        res += arr[i].re * arr[i].re;
    }
    return sqrt(res);
}

matrix::matrix(size_t rows, size_t cols) : _size(rows), arr(new vec[cols]) {}

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

vec matrix::operate(vec& v) {
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

