# cmake_minimum_required(VERSION <specify CMake version here>)
project(untitled)

set(CMAKE_CXX_STANDARD 11)

add_library(untitled SHARED library.cpp library.h)