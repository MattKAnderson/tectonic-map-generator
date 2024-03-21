add_library(space_partition space_partition.cpp)
target_link_libraries(space_partition Vector2D::Vector2D)

add_library(dev_space_partition tmp.cpp)
target_link_libraries(dev_space_partition Vector2D::Vector2D)

add_executable(stochastic_space_partition_visual_test test.cpp)
target_link_libraries(stochastic_space_partition_visual_test 
                      space_partition dev_space_partition timeit::timeit)
target_include_directories(stochastic_space_partition_visual_test PRIVATE ./)


add_definitions(-DCMAKE_EXPORT_COMPILE_COMMANDS=ON)