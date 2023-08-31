cmake_minimum_required(VERSION 3.16)
project(tpSCAN)

set(CMAKE_CXX_STANDARD 14)

add_executable(main tpSCAN.hpp tpSCAN.cpp main.cpp)