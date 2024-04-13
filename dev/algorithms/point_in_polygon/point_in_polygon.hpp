/*
 *   Need to solve the problem of determining if a point lies within a polygon
 *   that is defined by a series of vertices
 * 
 *   Can implement the Ray-casting algorithm defined on wikipedia
 */
#pragma once
#include <vector>
#include <Coordinate.hpp>


bool in_polygon(std::vector<Coordinate>& vertices, Coordinate c);
