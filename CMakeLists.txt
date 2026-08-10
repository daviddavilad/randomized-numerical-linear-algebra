cmake_minimum_required(VERSION 3.21)
project(randomized_numerical_linear_algebra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(scratch src/main.cpp)
target_include_directories(scratch PRIVATE include)
target_compile_options(scratch PRIVATE -Wall -Wextra -Wpedantic)