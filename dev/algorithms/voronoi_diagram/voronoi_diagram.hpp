#pragma once
#include <vector>
#include <random>
#include <Coordinate.hpp>

class VoronoiDiagram {
public:
    VoronoiDiagram(int seed = 0);
    void generate(int xsize, int ysize, int nseeds);
    void iteration();
    std::vector<std::vector<int>> get_diagram();

private:
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    int nregions;
};


std::vector<std::vector<int>> generate_voronoi_diagram(int xsize, int ysize, int nseeds);

std::vector<std::vector<int>> relax_voronoi_diagram(std::vector<std::vector<int>>& diagram);