cmake_minimum_required(VERSION 3.15)

project(HAZEN)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(./inc)

set(LIB_SRC
./src/FrictionMethods.cpp
./src/HydraulicShapes.cpp
./src/HydraulicUtil.cpp
./src/Headloss.cpp
)

set(TEST_SRC
./src/tests/hydraulic_util_test.cpp)

add_library(HAZEN SHARED ${LIB_SRC})

add_executable(test ${TEST_SRC})

target_link_libraries(test HAZEN)