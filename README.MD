#  TSP

## What's this?
This is an implementation of the "traveling salesman problem".
It generates a problem based on some parameters, and tries to
solve it in realtime using an evolutionary algorithm.

### Libraries
TSP uses the libraries glm and GLFW for windowing/drawing and linear algebra.
These libraries are not written by me, and can be found at:
Library | Github
--------|----------------------------
glm:    | https://github.com/g-truc/glm
GLFW:   | https://github.com/glfw/glfw

## How do I use it?
### Building
This project is configured via CMake.
CMake can be downloaded from https://cmake.org/download/

To build the project, run 
cmake -G [build system] -A [architecture] -B [output folder] to generate a 
project within your favorite build system.

If you're, for example, using visual studio 2022 to build the project for an x64 architecture,
and you want the project in a folder called "build"
the full command would be cmake -G "Visual studio 17 2022" -A x64 -B build.

Consult the cmake documentation for help:
https://cmake.org/runningcmake/

### Usage
tsp --[paramter1] [value] --[parameter2] [value2] ...

Parameter   | Description
------------|----------------------------------------------------------------------------
seed        | An integer that's entered to the random number generator
points      | An integer that determines the number of "cities" generated in the problem
population  | An integer that determines how many solutions are generated per generation
generations | An integer that determines how many generations to generate
mutation    | A floating point number between 0 and 1 which determines the chance that a solution is mutated for the next generation
kill        | A floating point number between 0 and 1 which determines the fraction of each generation which is killed
mix         | A floating point number between 0 and 1 which determines how much solutions are mixed with eachother when reproducing


