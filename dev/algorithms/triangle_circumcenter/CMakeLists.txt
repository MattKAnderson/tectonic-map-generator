add_executable(test_triangle_circumcenter test_triangle_circumcenter.cpp
                                          triangle_circumcenter.cpp)
target_include_directories(test_triangle_circumcenter PUBLIC ./)
target_link_libraries(test_triangle_circumcenter timeit::timeit)