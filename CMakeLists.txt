cmake_minimum_required(VERSION 3.16)
project(numerical_methods_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

add_executable(q2 src/q2.cpp)
add_executable(q3 src/q3.cpp)
add_executable(q4 src/q4.cpp)

# rk4.cpp is kept as a reference file; q4.cpp already includes RK4 in many versions.
