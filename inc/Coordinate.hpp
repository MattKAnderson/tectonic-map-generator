#pragma once
#include <vector>
#include <limits>

struct Coordinate {
    int x, y;
    bool operator==(const Coordinate& other) const {
        return x == other.x && y == other.y;
    };
};

struct CoordinateHash {
    std::size_t operator()(const Coordinate& c) const {
        return std::hash<long long>{}(
            c.x 
            * (std::numeric_limits<int>::max() + static_cast<long long>(1))
            + c.y
        );
    }
};