add_library(voronoi_diagram voronoi_diagram.cpp geometry.cpp LineClipper.cpp)
target_include_directories(voronoi_diagram PUBLIC ./)
 
add_executable(test_voronoi_diagram test.cpp)
target_link_libraries(test_voronoi_diagram voronoi_diagram Vector2D)

add_executable(test_fortunes_alg test_fortunes_alg.cpp)
target_link_libraries(test_fortunes_alg voronoi_diagram Vector2D)