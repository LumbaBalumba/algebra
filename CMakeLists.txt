cmake_minimum_required(VERSION 3.22)
project(algebra)

set(CMAKE_CXX_STANDARD 20)

add_library(algebra algebra.cpp complex.cpp vec.cpp matrix.cpp polynomial.cpp)
