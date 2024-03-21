add_library(boundary boundary.cpp)
target_link_libraries(boundary Vector2D::Vector2D)

add_executable(boundary_visual_test test_visually_boundary.cpp)
target_link_libraries(boundary_visual_test boundary tectonicSimulation timeit::timeit)
target_include_directories(boundary_visual_test PRIVATE ./)


add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)