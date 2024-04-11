#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <Coordinate.hpp>

template<typename T>
std::string to_string(const T& t) {
    return std::to_string(t);
}

template <>
std::string to_string<Coordinate>(const Coordinate& c) {
    return "[" + std::to_string(c.x) + ", " + std::to_string(c.y) + "]";
}

template <typename T>
void matrix_output(std::string filename, std::vector<std::vector<T>>& mat) {
    std::ofstream outfile;
    outfile.open(filename);
    outfile << "{\"matrix\": [";
    for (int i = 0; i < mat.size() - 1; i++) {
        outfile << "[";
        for (int k = 0; k < mat[i].size() - 1; k++) {
            outfile << to_string(mat[i][k]) << ", ";
        }
        outfile << to_string(mat[i].back()) << "], ";
    }
    outfile << "[";
    for (int k = 0; k < mat.back().size() - 1; k++) {
        outfile << to_string(mat.back()[k]) << ", ";
    }
    outfile << to_string(mat.back().back()) << "]]}";
};

template <typename T>
void vector_output(std::string filename, std::vector<T>& vec) {
    std::ofstream outfile;
    outfile.open(filename);
    outfile << "{\"vector\": [";
    for (int i = 0; i < vec.size() - 1; i++) {
        outfile << to_string(vec[i]) << ", ";
    }
    outfile << to_string(vec.back()) << "]}";
}