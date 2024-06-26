/*
 *   Geometry functions for implementing fortunes algorithm
 * 
 */
#pragma once
#include <cmath>
#include <vector>
#include <Coordinate.hpp>


double euclidean_distance(const RealCoordinate& c1, const RealCoordinate& c2);

//double toroidal_distance(
//    Coordinate& c1, Coordinate& c2, int xsize, int ysize
//);

//bool is_between(double x, double a, double b);

//RealCoordinate triangle_centroid(
//    const RealCoordinate& a, const RealCoordinate& b, const RealCoordinate& c
//);

RealCoordinate triangle_circumcenter(
    const RealCoordinate& a, const RealCoordinate& b, const RealCoordinate& c
);

//RealCoordinate lines_intercept(
//    const RealCoordinate& l1a, const RealCoordinate& l1b, 
//    const RealCoordinate& l2a, const RealCoordinate& l2b
//);

//bool line_segments_intersect(
//    const RealCoordinate& l1a, const RealCoordinate& l1b,
//    const RealCoordinate& l2a, const RealCoordinate& l2b
//);

//bool in_polygon(std::vector<RealCoordinate>& vertices, const RealCoordinate& c);

//bool num_is_between(double x, double a, double b);

//RealCoordinate midpoint(const RealCoordinate& a, const RealCoordinate& b); 

std::vector<double> quadratic_roots(double a, double b, double c);

double parabola_x_from_y(double directrix, const RealCoordinate& focus, double y);

RealCoordinate parabola_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
);

double parabolae_y_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
);