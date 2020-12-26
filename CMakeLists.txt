cmake_minimum_required(VERSION 3.15)

project(HAZEN)

include_directories("./ext" "C:/Program Files (x86)/Intel/oneAPI/mkl/2021.1.1/include")
# link_directories("C:/Program Files (x86)/Intel/oneAPI/mkl/2021.1.1/lib/intel64")

add_definitions(-D_USE_MATH_DEFINES)
# add_definitions(-D_USE_MATH_DEFINES -DEIGEN_USE_MKL_ALL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

add_library(HAZEN SHARED ${LIB_SRC})

# target_link_libraries(HAZEN
# mkl_intel_lp64
# mkl_sequential
# mkl_core
# )

# target_compile_options(HAZEN PUBLIC "/arch:AVX2" "/openmp" "/fp:fast")

add_executable(test ${TEST_SRC})

target_link_libraries(test HAZEN)