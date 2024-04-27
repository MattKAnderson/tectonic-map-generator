/*
 *   Geometry functions for implementing fortunes algorithm
 * 
 */
#pragma once
#include <cmath>
#include <vector>
#include <Coordinate.hpp>

double euclidean_distance(RealCoordinate& c1, RealCoordinate& c2);

double toroidal_distance(
    Coordinate& c1, Coordinate& c2, int xsize, int ysize
);

bool is_between(double x, double a, double b);

RealCoordinate triangle_centroid(
    const RealCoordinate& a, const RealCoordinate& b, const RealCoordinate& c
);

RealCoordinate triangle_circumcenter(
    const RealCoordinate& a, const RealCoordinate& b, const RealCoordinate& c
);

RealCoordinate lines_intercept(
    const RealCoordinate& l1a, const RealCoordinate& l1b, 
    const RealCoordinate& l2a, const RealCoordinate& l2b
);

RealCoordinate midpoint(const RealCoordinate& a, const RealCoordinate& b); 

std::vector<double> quadratic_roots(double a, double b, double c);

double parabola_x_from_y(double directrix, const RealCoordinate& focus, double y);

RealCoordinate parabola_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
);

double parabolae_y_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
);