# HAZEN - Hydraulic AnalyZer ENgine - Version 1.0

## Introduction
HAZEN is a steady state one-dimensional hydraulic computation library written in C++. Open channel flow, pressure driven flow, and hydraulic structures such as weirs and orifices can all be modeled. Branching flows that diverge, combine, or networks with multiple outfalls and inflows can all be modeled seamlessly. There are no dependencies other than the C++ standard library. HAZEN can be used to construct hydraulic networks of a size bounded only by the host computer's available memory and computational power.

## Building the Library
The library can be build with cmake and easily added into any project that currently uses cmake as its build system.

To build the library, run cmake --build into a directory of your choosing. For example from the root directory of the source code:

mkdir build
cd build
cmake ..
cd ..
cmake --build ./build

will build a shared library (.dll on Windows) into the build folder (or a subdirectory therein).

## Tutorial / Examples

## Documentation
The full documentation for the APIs can be built using Doxygen. For a thourough description of the computations performed for each link in a network and the implementation details of the network flow solver, see the Theory PDF in the project root directory. 

## TODO
1. Tutorial / Examples
2. Optimize network solver for large numbers of links
3. Additional Hydraulic Components:
    * Manhole
    * Pump/Special Loss
    * Bar Screen
    * Flow Control Valve
    * in-line minor loss K vlaues (Reducers/Meters/Bends/etc.)
    * Three-way minor loss (Tee/Wye)
    * Four-way minor loss (Cross)
    * Baffle
    * Granular Filter
    * Diffuser
    * Slide Gate
    * Stop Valves (Plug and Gate)
    * Check Valve
    * Flume