#include <point_in_polygon.hpp>
#include <iostream>

bool between(double x, double a, double b) {
    return std::min(a, b) < x && x < std::max(a, b);
}

// assumes flat space
bool in_polygon(std::vector<Coordinate>& vertices, Coordinate c) {

    if (vertices.size() < 3) {
        return false;
    }

    int intersect_count = 0;
    if (vertices[0].y == c.y
        && vertices[0].x < c.x 
        && between(c.y, vertices[1].y, vertices[vertices.size() - 2].y)
    ) {
        intersect_count = 1;
    }
    
    for (int i = 1; i < vertices.size(); ++i) {
        Coordinate v1 = vertices[i - 1];
        Coordinate v2 = vertices[i];
        if (v2.y == c.y) {
            if (v2.x <= c.x && i != vertices.size() - 1) {
                intersect_count += between(c.y, v1.y, vertices[i + 1].y);
                ++i;
            }
        }
        else if (v2.x == v1.x) {
            intersect_count += v2.x < c.x && between(c.y, v1.y, v2.y);
        }
        else if (v2.y != v1.y) {
            double m = static_cast<double>(v2.y - v1.y) / (v2.x - v1.x);
            double b = v2.y - m * v2.x;
            double x_intersect = (c.y - b) / m;
            if (x_intersect < c.x && between(x_intersect, v1.x, v2.x)) {
                ++intersect_count;
            }
        }
    }
    return intersect_count % 2;
}