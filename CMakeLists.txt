cmake_minimum_required(VERSION 3.9)
project(algebra)

set(CMAKE_CXX_STANDARD 23)

add_library(algebra algebra.cpp complex.cpp vec.cpp matrix.cpp polynomial.cpp)
add_library(algebra_shared SHARED algebra.cpp complex.cpp vec.cpp matrix.cpp polynomial.cpp)