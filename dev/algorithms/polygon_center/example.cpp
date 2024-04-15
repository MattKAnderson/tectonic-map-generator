#include <algorithm>
#include <chrono>
#include <limits>
#include <random>
#include <vector>
#include <point_in_polygon.hpp>
#include <polygon_center.hpp>
#include <Coordinate.hpp>
#include <FileIO.hpp>


std::vector<Coordinate> sample_polygon(
    std::vector<Coordinate>& vertices, int nsamples, int xlim, int ylim,
    std::mt19937_64& rng
);

std::vector<Coordinate> translate(
    std::vector<Coordinate> vertices, int dx, int dy, int xlim, int ylim
);

int main() {

    int seed = 1711181;
    int nsamples = 3000;
    int xlim = 1024;
    int ylim = 1024;
    int iters = 20000;
    int k = 1.0;
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
    std::mt19937_64 rng(seed);
    std::vector<Coordinate> polygon_samples;
    polygon_samples = sample_polygon(vertices, nsamples, xlim, ylim, rng);

    CentroidFinder centroid_finder(seed);
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Coordinate> vcentroid = {
        centroid_finder.toroidal_centroid(
            polygon_samples, xlim, ylim, iters, k, true
        )
    };
    auto end = std::chrono::high_resolution_clock::now();
    std::vector<Coordinate> visit_history;
    visit_history = centroid_finder.get_visit_history();
    std::vector<double> cost_history;
    cost_history = centroid_finder.get_cost_history();
    vector_output("example_vertices.json", vertices);
    vector_output("example_polygon_samples.json", polygon_samples);
    vector_output("example_centroid.json", vcentroid);
    vector_output("example_visit_history.json", visit_history);
    vector_output("example_cost_history.json", cost_history);
    auto elapsed_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
        end - start
    ).count();
    std::cout << "Took " << elapsed_time / 1000000.0 
              << " ms to find the centroid\n"; 

    std::vector<Coordinate> samples_2 = translate(
        polygon_samples, 400, 0, xlim, ylim
    );
    vcentroid = {
        centroid_finder.toroidal_centroid(
            samples_2, xlim, ylim, iters, k, true
        )
    };
    visit_history = centroid_finder.get_visit_history();
    cost_history = centroid_finder.get_cost_history();
    vector_output("example_translated_samples.json", samples_2);
    vector_output("example_translated_centroid.json", vcentroid);
    vector_output("example_translated_visit_history.json", visit_history);
    vector_output("example_translated_cost_history.json", cost_history);
}

std::vector<Coordinate> sample_polygon(
    std::vector<Coordinate>& vertices, int nsamples, int xlim, int ylim,
    std::mt19937_64& rng
) {
    std::uniform_int_distribution<int> y_domain(0, ylim - 1);
    std::uniform_int_distribution<int> x_domain(0, xlim - 1);
    std::vector<Coordinate> polygon_samples;
    polygon_samples.reserve(nsamples);
    for (int i = 0; i < nsamples; i++) {
        Coordinate sample = {x_domain(rng), y_domain(rng)};
        while (!in_polygon(vertices, sample)) {
            sample = {x_domain(rng), y_domain(rng)};
        }
        polygon_samples.push_back(sample);
    }
    return polygon_samples; 
}

std::vector<Coordinate> translate(
    std::vector<Coordinate> vertices, int dx, int dy, int xlim, int ylim
) {
    for (auto it = vertices.begin(); it != vertices.end(); it++) {
        it->x = it->x + dx;
        if (it->x >= xlim) {
            it->x -= xlim;
        }
        it->y = it->y + dy;
        if (it->y >= ylim) {
            it->y -= ylim;
        }
    }
    return vertices;
}