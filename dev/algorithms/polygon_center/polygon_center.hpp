
/*
 *   Goal is to produce a method that determines a polygon center on a 
 *   surface that is not euclidean, e.g. a taurus surface or sphere
 * 
 *   Initial idea is to use a monte carlo approach
 */

// core algorithm take a series of points and find the approximate point that
// minimizes distance to all these points
#include <iostream>
#include <random>
#include <vector>
#include <Coordinate.hpp>
#include <Vector2D.hpp>


Coordinate find_toroidal_centroid(
    std::vector<Coordinate>& points, int xmax, int ymax, int iters, double k
);
