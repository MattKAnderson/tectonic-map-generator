#pragma once
#include <algorithm>
#include <queue>
#include <vector>
#include <random>
#include <Coordinate.hpp>


struct Node {
    RealCoordinate coord;
    std::vector<Node*> edges;
};


class VoronoiDiagram {
public:
    VoronoiDiagram(int seed = 0);
    void generate(int xsize, int ysize, int nseeds);
    void generate_from_seeds(
        std::vector<Coordinate>& seeds, int xsize, int ysize
    );
    void lloyd_iteration();
    void nesting_iteration(int nseeds);
    void nesting_iteration_from_seeds(std::vector<Coordinate>& seeds);
    std::vector<std::vector<int>> get_diagram();

private:
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    std::vector<Node> vertices;
    std::vector<Node> delauney_triangulation;
    int nregions;

    void fortunes_algorithm(
        std::vector<RealCoordinate>& seeds, int xsize, int ysize
    );
};
