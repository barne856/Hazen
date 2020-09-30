cmake_minimum_required(VERSION 3.15)

project(HAZEN)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(./inc)


set(TEST_SRC
./src/tests/vectors_example.cpp)

add_executable(test ${TEST_SRC})