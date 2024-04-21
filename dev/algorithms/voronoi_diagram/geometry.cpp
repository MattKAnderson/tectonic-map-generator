#include <geometry.hpp>


double euclidean_distance(RealCoordinate& c1, RealCoordinate& c2) {
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return std::sqrt(dx * dx + dy*dy);
}


double toroidal_distance(
    Coordinate& c1, Coordinate& c2, int xsize, int ysize
) {
    double dx = std::abs(c1.x - c2.x);
    if (dx > xsize / 2) {
        dx = xsize - dx;
    }
    double dy = std::abs(c1.y - c2.y);
    if (dy > ysize / 2) {
        dy = ysize - dy;
    }
    return std::sqrt(dx * dx + dy * dy);
}


RealCoordinate triangle_centroid(
    RealCoordinate& a, RealCoordinate& b, RealCoordinate& c
) {
    return {(a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0};
}


RealCoordinate triangle_circumcenter(
    RealCoordinate& p1, RealCoordinate& p2, RealCoordinate& p3
) {
    double a = p1.x - p2.x;
    double b = p1.x - p2.y;
    double c = p2.x - p3.x;
    double d = p2.y - p3.y;
    double u = 0.5 * (p1.x * p1.x + p1.y * p1.y - p2.x * p2.x - p2.y * p2.y);
    double v = 0.5 * (p2.x * p2.x + p2.y * p2.y - p3.x * p3.x - p3.y * p3.y);

    double inv_denominator = 1.0 / (a * d - b * c);
    double x = (d * u - b * v)  * inv_denominator;
    double y = (a * v - c * u) * inv_denominator;

    return {x, y};
}


RealCoordinate lines_intercept(
    RealCoordinate& l1a, RealCoordinate& l1b, 
    RealCoordinate& l2a, RealCoordinate& l2b
) {
    double m1 = (l1b.y - l1a.y) / (l1b.x - l1a.x);
    double b1 = l1a.y - m1 * l1a.x;
    double m2 = (l2b.y - l2a.y) / (l2b.x -l2a.x);
    double b2 = l2a.y - m2 * l2a.x;
    double x0 = (b2 - b1) / (m1 - m2);
    double y0 = m1 * x0 + b1;
    return {x0, y0};
}


RealCoordinate midpoint(RealCoordinate& a, RealCoordinate& b) {
    return {(a.x + b.x) / 2.0, (a.y + b.y) / 2.0};
}


std::vector<double> quadratic_roots(double a, double b, double c) {
    double inner = b * b - 4.0 * a * c;
    if (inner < 0) {
        return {};
    }
    else if (inner == 0) {
        return {(-b / (2.0 * a))};
    }
    else {
        double inv_denominator = 0.5 / a;
        double root_inner = std::sqrt(inner);
        return {
            (-b + root_inner) * inv_denominator,
            (-b - root_inner) * inv_denominator
        };
    }
}


bool is_between(double x, double a, double b) {
    return std::min(a, b) < x && x < std::max(a, b);
}


double parabola_x_from_y(double directrix, RealCoordinate& focus,double y) {
    return (y - focus.y) * (y - focus.y) / (2.0 * (focus.x - directrix))
           + (focus.x + directrix) * 0.5;
}


RealCoordinate parabola_intercept(
    double directrix, RealCoordinate& focus1, RealCoordinate& focus2
) {
    /*
     Will return the intersection point of the two parabolas that is 
     between the two focuses (on y-axis).
     */
    auto [xf1, yf1] = focus1;
    auto [xf2, yf2] = focus2;
    double xd = directrix;

    double a = xf2 - xf1;
    double b = 2.0 * (yf2 * (xf1 - xd) - yf1 * (xf2 - xd));
    double c = (yf1 * yf1) * (xf2 - xd) - (yf2 * yf2) * (xf1 - xd) 
             - (xf2 - xf1) * (xf1 - xd) * (xf2 - xd);

    // There will always exist 2 intercepts for 2 parabolas defined 
    // in this way. 
    std::vector<double> roots = quadratic_roots(a, b, c);
    double y = is_between(roots[0], yf1, yf2) ? roots[0] : roots[1];
    return {parabola_x_from_y(directrix, focus1, y), y};
}


double parabolae_y_intercept(
    double directrix, RealCoordinate& focus1, RealCoordinate& focus2
) {
    auto [xf1, yf1] = focus1;
    auto [xf2, yf2] = focus2;
    double xd = directrix;

    double a = xf2 - xf1;
    double b = 2.0 * (yf2 * (xf1 - xd) - yf1 * (xf2 - xd));
    double c = (yf1 * yf1) * (xf2 - xd) - (yf2 * yf2) * (xf1 - xd) 
             - (xf2 - xf1) * (xf1 - xd) * (xf2 - xd);

    std::vector<double> roots = quadratic_roots(a, b, c);
    return is_between(roots[0], yf1, yf2) ? roots[0] : roots[1];
}