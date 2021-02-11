cmake_minimum_required(VERSION 3.15)

project(HAZEN)

include_directories("./ext/Eigen")

add_definitions(-D_USE_MATH_DEFINES)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(./inc)

set(LIB_SRC
./src/Hazen/FrictionMethods.cpp
./src/Hazen/Headloss.cpp
./src/Hazen/HydraulicComponents.cpp
./src/Hazen/HydraulicLinks.cpp
./src/Hazen/HydraulicNetwork.cpp
./src/Hazen/HydraulicShapes.cpp
./src/Hazen/Hydraulics.cpp
)

set(TEST_SRC
./src/tests/channel_series_test.cpp)

add_library(HAZEN STATIC ${LIB_SRC})
