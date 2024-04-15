#include <iostream>
#include <FileIO.hpp>
#include <voronoi_diagram.hpp>
#include <PoissonDiskSampler.hpp>

int main() {
    
    int seed = 72335649;
    int xsize = 1024;
    int ysize = xsize;
    int k = 30;
    double radius = 32.0;

    PoissonDiskSampler sampler(seed);
    std::vector<Coordinate> seeds = sampler.fill_grid(xsize, radius, k);
    VoronoiDiagram voronoi_diagram(seed);
    voronoi_diagram.generate_from_seeds(seeds, xsize, ysize);
    for (int i = 0; i < 4; i++) {
        voronoi_diagram.lloyd_iteration();
    }
    std::vector<std::vector<int>> diagram = voronoi_diagram.get_diagram();
    matrix_output("poisson_voronoi_diagram.json", diagram);

    // nesting is still not a great approach. It would be nice to group cells
    // together into larger regions. Maybe I need to try a BFS approach for 
    // "aggregating" voronio cells together
    voronoi_diagram.nesting_iteration(16);
    diagram = voronoi_diagram.get_diagram();
    matrix_output("nested_poisson_voronoi_diagram.json", diagram);
}