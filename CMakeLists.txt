cmake_minimum_required(VERSION 3.27.0)
project(map_scripts VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)

include_directories(../Vector2D ../timeit ./inc)

add_library(tectonicSimulation STATIC src/tectonicSimulation.cpp src/plate.cpp)
#add_executable(test_simulation test.cpp)
#target_link_libraries(test_simulation tectonicSimulation)

add_subdirectory(dev/algorithms/boundary_identification)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
