#include <iostream>
#include <voronoi_diagram.hpp>
#include <FileIO.hpp>

int main() {

    int seed = 1284197328;
    int xsize = 1024;
    int ysize = xsize;
    int nseeds = 32;

    VoronoiDiagram voronoi_diagram(seed);
    voronoi_diagram.generate(xsize, ysize, nseeds);
    auto initial_diagram = voronoi_diagram.get_diagram();
    std::cout << "Created initial diagram" << std::endl;

    voronoi_diagram.iteration();
    auto iterated_diagram = voronoi_diagram.get_diagram();
    std::cout << "Completed first iteration" << std::endl;

    voronoi_diagram.iteration();
    voronoi_diagram.iteration();
    auto iterated_x3_diagram = voronoi_diagram.get_diagram();

    for (int i = 3; i < 20; i++) {
        voronoi_diagram.iteration();
    }
    auto iterated_x20_diagram = voronoi_diagram.get_diagram();

    matrix_output("initial_diagram.json", initial_diagram);
    matrix_output("iterated_diagram.json", iterated_diagram);
    matrix_output("iterated_x3_diagram.json", iterated_x3_diagram);
    matrix_output("iterated_x20_diagram.json", iterated_x20_diagram);
}