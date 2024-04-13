add_subdirectory(point_in_polygon)

add_library(polygon_center polygon_center.cpp)
target_include_directories(polygon_center PUBLIC ./)
target_link_libraries(polygon_center PUBLIC Vector2D)

add_executable(test_polygon_center test.cpp)
target_link_libraries(test_polygon_center polygon_center)