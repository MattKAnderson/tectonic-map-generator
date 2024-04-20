/*
 *   Geometry functions for implementing fortunes algorithm
 * 
 */
#include <cmath>
#include <vector>
#include <Coordinate.hpp>

double euclidean_distance(RealCoordinate& c1, RealCoordinate& c2);

double toroidal_distance(
    Coordinate& c1, Coordinate& c2, int xsize, int ysize
);

bool is_between(double x, double a, double b);

RealCoordinate triangle_centroid(
    RealCoordinate& a, RealCoordinate& b, RealCoordinate& c
);

RealCoordinate triangle_centroid(
    RealCoordinate& a, RealCoordinate& b, RealCoordinate& c
);

RealCoordinate lines_intercept(
    RealCoordinate& l1a, RealCoordinate& l1b, 
    RealCoordinate& l2a, RealCoordinate& l2b
);

RealCoordinate midpoint(RealCoordinate& a, RealCoordinate& b); 

RealCoordinate triangle_circumcircle_center(
    RealCoordinate& a, RealCoordinate& b, RealCoordinate& c
);

std::vector<double> quadratic_roots(double a, double b, double c);


double parabola_x_from_y(double directrix, RealCoordinate& focus,double y);

RealCoordinate parabola_intercept(
    double directrix, RealCoordinate& focus1, RealCoordinate& focus2
);

double parabolae_y_intercept(
    double directrix, RealCoordinate& focus1, RealCoordinate& focus2
);