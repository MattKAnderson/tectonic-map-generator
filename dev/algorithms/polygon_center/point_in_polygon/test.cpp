#include <iostream>
#include <vector>
#include <Coordinate.hpp>
#include <point_in_polygon.hpp>

int main() {
    Coordinate sample(3, 7);

    std::vector<Coordinate> polygon = {
        {2, 2}, {3, 6}, {8, 2}
    };

    std::vector<Coordinate> vertices = {
        {147, 189}, {139, 248}, {64, 268}, {0, 245}, {6, 172},
        {38, 108}, {116, 51}, {233, 16}, {424, 0}, {621, 23},
        {755, 55}, {786, 124}, {808, 207}, {764, 289}, {652, 313},
        {521, 304}, {607, 227}, {689, 223}, {723, 204}, {709, 132},
        {626, 127}, {562, 157}, {531, 231}, {508, 255}, {428, 281},
        {416, 254}, {440, 170}, {476, 142}, {474, 101}, {395, 94},
        {326, 130}, {312, 168}, {332, 213}, {354, 247}, {336, 270},
        {270, 272}, {226, 256}, {227, 207}, {227, 140}, {197, 116},
        {145, 119}, {123, 131}, {116, 145}, {129, 161}, {149, 172},
        {147, 189}
    };
    Coordinate sample2(606, 142);
    
    if (in_polygon(vertices, sample2)) {
        std::cout << "The sample is in the polygon\n";
    }
    else {
        std::cout << "The sample was not in the polygon\n";
    }
}