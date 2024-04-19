#pragma once
#include <vector>
#include <limits>
#include <geometry.hpp>


struct RealCoordinate {
    double x, y;
    bool operator==(const RealCoordinate& other) const {
        return x == other.x && y == other.y;
    }
};


struct Coordinate {
    int x, y;
    bool operator==(const Coordinate& other) const {
        return x == other.x && y == other.y;
    };
};


namespace std {
template <> 
struct hash<RealCoordinate> {
    size_t operator()(const RealCoordinate& c) const {
        size_t hx = hash<double>()(c.x);
        size_t hy = hash<double>()(c.y);
        return hx ^ (hy << 1);
    }
};

template <>
struct hash<Coordinate> {
    size_t operator()(const Coordinate& c) const {
        return hash<long long>{}(
            c.x 
            * (numeric_limits<int>::max() + static_cast<long long>(1))
            + c.y
        );
    }
};
} // namespace std