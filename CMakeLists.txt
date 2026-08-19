cmake_minimum_required(VERSION 3.21)
project(randomized_numerical_linear_algebra LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)

add_library(rnla src/matrix.cpp src/linalg.cpp src/random.cpp src/test_matrices.cpp src/rsvd.cpp)
target_include_directories(rnla PUBLIC include)
target_link_libraries(rnla PUBLIC ${ACCELERATE_FRAMEWORK})
target_compile_options(rnla PRIVATE -Wall -Wextra -Wpedantic)

add_executable(scratch src/main.cpp)
target_link_libraries(scratch PRIVATE rnla)
target_compile_options(scratch PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(tests tests/test_linalg.cpp)
target_link_libraries(tests PRIVATE rnla)
target_compile_options(tests PRIVATE -Wall -Wextra -Wpedantic)
add_test(NAME linalg COMMAND tests)

add_executable(bench_spectrum benchmarks/bench_spectrum.cpp)
target_link_libraries(bench_spectrum PRIVATE rnla)
target_compile_options(bench_spectrum PRIVATE -Wall -Wextra -Wpedantic)