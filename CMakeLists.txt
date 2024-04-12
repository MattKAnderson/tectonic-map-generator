cmake_minimum_required(VERSION 3.27.0)
project(map_scripts VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

include(FetchContent)
FetchContent_Declare(
  vector2d
  GIT_REPOSITORY https://github.com/MattKAnderson/Vector2D
  GIT_TAG main
)
FetchContent_Declare(
    cpp-timeit
    GIT_REPOSITORY https://github.com/MattKAnderson/cpp-timeit
    GIT_TAG main
)
FetchContent_MakeAvailable(vector2d)
FetchContent_MakeAvailable(cpp-timeit)

include_directories(inc)

add_library(tectonicSimulation STATIC src/tectonicSimulation.cpp src/plate.cpp)
target_link_libraries(tectonicSimulation PUBLIC Vector2D::Vector2D)
# add_executable(test_simulation test/test.cpp)
# target_link_libraries(test_simulation tectonicSimulation timeit::timeit)

add_subdirectory(dev/algorithms/boundary_identification)
add_subdirectory(dev/algorithms/space_partitioning)
add_subdirectory(dev/algorithms/voronoi_diagram)
add_subdirectory(dev/algorithms/poisson_disk_sampling)
add_subdirectory(dev/algorithms/polygon_center)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)
