cmake_minimum_required(VERSION 3.31)
project(wildcard_FSPBWT)

set(CMAKE_CXX_STANDARD 20)

add_executable(wildcard_FSPBWT main.cpp
        wildcard-FSPBWT.h
        tools.h)
