#include <PoissonDiskSampler.hpp>
#include <FileIO.hpp>

int main() {
    int seed = 1823983;
    int size = 1024;
    int k = 30;
    double fill_radius = 0.15;
    double grid_radius = 128.0;

    PoissonDiskSampler sampler(seed);

    std::vector<Vector2D<double>> unit_fill = sampler.fill(fill_radius, k);

    std::vector<Coordinate> grid_fill = sampler.fill_grid(1024, grid_radius, k);

    vector_output("unit_fill.json", unit_fill);
    vector_output("grid_fill.json", grid_fill);
}