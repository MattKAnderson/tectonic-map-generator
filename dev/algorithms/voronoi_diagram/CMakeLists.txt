add_library(voronoi_diagram voronoi_diagram.cpp)
target_include_directories(voronoi_diagram PUBLIC ./)
 
add_executable(test_voronoi_diagram test.cpp)
target_link_libraries(test_voronoi_diagram voronoi_diagram)