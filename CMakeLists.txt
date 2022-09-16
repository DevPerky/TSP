cmake_minimum_required(VERSION 3.24)

project(TSP)

set(CMAKE_BUILD_TYPE Release)

add_executable(tsp
    src/main.cpp
    src/tsp.cpp
)

add_subdirectory(glfw)
add_subdirectory(glm)
find_package(OpenGL REQUIRED)

target_link_libraries(tsp glfw)
target_link_libraries(tsp OpenGL::GL)
target_include_directories(tsp PUBLIC glfw)
target_include_directories(tsp PUBLIC glm)