cmake_minimum_required(VERSION 3.31)
project(wildcard_FSPBWT)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXE_LINKER_FLAGS "-static")
add_executable(wildcard_FSPBWT main.cpp
        wildcard-FSPBWT.h
        tools.h)
