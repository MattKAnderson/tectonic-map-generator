add_subdirectory(algorithms)

add_executable(voronoi_from_poisson_disk voronoi_from_poisson_disk.cpp)
target_link_libraries(voronoi_from_poisson_disk voronoi_diagram 
                                                poisson_disk_sampler
                                                Vector2D)