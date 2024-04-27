#include <geometry.hpp>
#include <iostream>

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
    const RealCoordinate& p1, const RealCoordinate& p2, const RealCoordinate& p3
) {
    double a = p1.x - p2.x;
    double b = p1.y - p2.y;
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
    const RealCoordinate& l1a, const RealCoordinate& l1b, 
    const RealCoordinate& l2a, const RealCoordinate& l2b
) {
    double m1 = (l1b.y - l1a.y) / (l1b.x - l1a.x);
    double b1 = l1a.y - m1 * l1a.x;
    double m2 = (l2b.y - l2a.y) / (l2b.x -l2a.x);
    double b2 = l2a.y - m2 * l2a.x;
    double x0 = (b2 - b1) / (m1 - m2);
    double y0 = m1 * x0 + b1;
    return {x0, y0};
}


bool on_line_segment(
    RealCoordinate& l1, RealCoordinate& l2, RealCoordinate& p
) {
     
}


RealCoordinate midpoint(const RealCoordinate& a, const RealCoordinate& b) {
    return {(a.x + b.x) / 2.0, (a.y + b.y) / 2.0};
}


std::vector<double> quadratic_roots(double a, double b, double c) {
    double inner = b * b - 4.0 * a * c;
    const double eps = 1e-10;
    if (inner < 0) {
        return {};
    }
    else if (0 - eps < inner && inner < 0 + eps) {
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


double parabola_x_from_y(double directrix, const RealCoordinate& focus, double y) {
    return (y - focus.y) * (y - focus.y) / (2.0 * (focus.x - directrix))
           + (focus.x + directrix) * 0.5;
}


RealCoordinate parabola_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
) {
    const double eps = 1e-8;
    auto [xf1, yf1] = focus1;
    auto [xf2, yf2] = focus2;
    double xd = directrix;

    if (xf1 - xd < 0 + eps && xf1 - xd > 0 - eps) {
        return {parabola_x_from_y(directrix, focus2, yf1), yf1};
    }
    if (xf2 - xd < 0 + eps && xf2 - xd > 0 - eps) {
        return {parabola_x_from_y(directrix, focus1, yf2), yf2};
    }

    double y; 
    double a = xf2 - xf1;
    if (a == 0) {
        y = 0.5 * (yf1 + yf2);
        return {parabola_x_from_y(directrix, focus1, y), y};
    }
    double b = 2.0 * (yf2 * (xf1 - xd) - yf1 * (xf2 - xd));
    double c = (yf1 * yf1) * (xf2 - xd) - (yf2 * yf2) * (xf1 - xd) 
             - (xf2 - xf1) * (xf1 - xd) * (xf2 - xd);

    std::vector<double> roots = quadratic_roots(a, b, c);
    if (roots.size() == 0) {
        std::cout << "ran into case with no roots" << std::endl;
        std::cout << "f1: (" << xf1 << ", " << yf1 << ")\n"
                  << "f2: (" << xf2 << ", " << yf2 << ")\n"
                  << "d:  " << directrix << std::endl;
    }
    if (roots.size() == 1) {
        y = roots[0];
    }
    if (xf1 > xf2) {
        y = std::min(roots[0], roots[1]);
    }
    else {
        y = std::max(roots[0], roots[1]);
    }
    
    return {parabola_x_from_y(directrix, focus1, y), y};
}

/* 
 Warn: This method is designed such that the order the focuses are given
 matters and determines which of the 2 possible y values to return.
 */
double parabolae_y_intercept(
    double directrix, const RealCoordinate& focus1, const RealCoordinate& focus2
) {
    auto [xf1, yf1] = focus1;
    auto [xf2, yf2] = focus2;
    double xd = directrix;

    double a = xf2 - xf1;
    if (a == 0) {
        return 0.5 * (yf1 + yf2);
    }
    double b = 2.0 * (yf2 * (xf1 - xd) - yf1 * (xf2 - xd));
    double c = (yf1 * yf1) * (xf2 - xd) - (yf2 * yf2) * (xf1 - xd) 
             - (xf2 - xf1) * (xf1 - xd) * (xf2 - xd);

    std::vector<double> roots = quadratic_roots(a, b, c);
    if (roots.size() == 1) {
        return roots[0];
    }
    if (xf1 > xf2) {
        return std::min(roots[0], roots[1]);
    }
    else {
        return std::max(roots[0], roots[1]);
    }
}