#include <triangle_circumcenter.hpp>

RealCoordinate lin_alg_solution(
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

RealCoordinate perp_bisector_solution(
    const RealCoordinate& p1, const RealCoordinate& p2, const RealCoordinate& p3
) {
    double x1 = 0.5 * (p2.x + p1.x);
    double y1 = 0.5 * (p2.y + p1.y);
    double m1 = - (p2.x - p1.x) / (p2.y - p1.y); 

    double x2 = 0.5 * (p3.x + p2.x);
    double y2 = 0.5 * (p3.y + p2.y);
    double m2 = - (p3.x - p2.x) / (p3.y - p2.y);

    double y_int = (y1 - m1 * (y2 / m2 - x2 +    x1)) / (1 - (m1 / m2));
    double x_int = (y_int - y1) / m1 + x1;

    return {x_int, y_int};
}