cmake_minimum_required(VERSION 3.15)

project(HAZEN)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

include_directories(./inc)

set(LIB_SRC
./src/FrictionMethods.cpp
./src/Headloss.cpp
./src/HydraulicComponents.cpp
./src/HydraulicLinks.cpp
./src/HydraulicNetwork.cpp
./src/HydraulicShapes.cpp
./src/HydraulicUtil.cpp
)

set(TEST_SRC
./src/tests/hydraulic_util_test.cpp)

add_library(HAZEN SHARED ${LIB_SRC})

add_executable(test ${TEST_SRC})

target_link_libraries(test HAZEN)