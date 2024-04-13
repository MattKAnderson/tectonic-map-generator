#include <point_in_polygon.hpp>

bool between(double x, double a, double b) {
    return std::min(a, b) <= x && x <= std::max(a, b);
}

// assumes flat space
bool in_polygon(std::vector<Coordinate>& vertices, Coordinate c) {
    
    int intersect_count = 0;
    for (int i = 1; i < vertices.size(); i++) {
        Coordinate v1 = vertices[i - 1];
        Coordinate v2 = vertices[i];
        if (v2.y == v1.y) {
            continue;
        }
        if (v2.x == v1.x) {
            intersect_count += v2.x < c.x;
            continue; 
        }
        double m = (v2.y - v1.y) / (v2.x - v1.x);
        double b = v2.y - m * v2.x;
        double x_intersect = (c.y - b) / m;
        if (x_intersect < c.x && between(x_intersect, v1.x, v2.x)) {
            ++intersect_count; 
        }

    }
    return intersect_count % 2;
}