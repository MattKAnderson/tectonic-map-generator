#include <iostream>
#include <vector>
#include <Coordinate.hpp>
#include <point_in_polygon.hpp>

int main() {
    Coordinate sample(3, 7);

    std::vector<Coordinate> polygon = {
        {2, 2}, {3, 6}, {8, 2}
    };

    bool is_in = in_polygon(polygon, sample);

    if (in_polygon(polygon, sample)) {
        std::cout << "The sample is in the polygon\n";
    }
    else {
        std::cout << "The sample was not in the polygon\n";
    }
}