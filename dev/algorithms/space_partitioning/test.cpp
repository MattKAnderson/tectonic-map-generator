#include <fstream>
#include <iostream>
#include <timeit.hpp>

#include "space_partition.hpp"
#include "tmp.hpp"

template <typename T>
void matrix_output(std::string filename, std::vector<std::vector<T>>& mat);

template <typename T>
void vector_output(std::string filename, std::vector<T>& vec);

int main() {
    
    int xsize = 1024;
    int ysize = xsize;
    int num_seeds = 16;

    std::mt19937_64 rng(325);
    std::geometric_distribution dist(0.33);

    //space_partition::StochasticSpacePartitioner partitioner(rng, dist);

    dev::SpacePartitioner partitioner(31552);

    //partitioner.random_walk_partition(512, 512, 8, 0.03);
    partitioner.voronoi_partition(xsize, ysize, num_seeds);

    dev::SpacePartitioner partitioner2(8172);

    std::vector<std::vector<int>> region_map, border_map;
    std::vector<std::vector<dev::Coordinate>> borderlines;
    std::vector<dev::Coordinate> traversal_history;
    region_map = partitioner.get_region_map();
    border_map = partitioner.get_border_map();
    borderlines = partitioner.get_borderlines();
    traversal_history = partitioner.get_traversal_history();


    timeIt::print(
        [](decltype(partitioner)& partitioner, int a, int b, int c) {
            partitioner.voronoi_partition(a, b, c);
        },
        partitioner2,
        xsize,
        ysize,
        num_seeds
    );

/*
    for (auto& [field, stat] : stats) {
        std::cout << field << ": " << stat << " ms" << std::endl;
    }
*/
/*
    m = partitioner.iterative_growth_partition(2048, 2048, 16);

    timeIt::print(
        [](decltype(partitioner)& partitioner, int a, int b, int c) { 
            partitioner.iterative_growth_partition(a, b, c); 
        },
        partitioner,
        2048, 
        2048, 
        16
    );
*/
    std::cout << "Writing region partition to: region_partitions.json" << std::endl;
    matrix_output("region_partitions.json", region_map);

    std::cout << "Writing boundary map to: border_map.json" << std::endl;
    matrix_output("border_map.json", border_map);

    std::cout << "Writing borderlines to: borderlines.json" << std::endl;
    matrix_output("borderlines.json", borderlines);

    std::cout << "Writing traversal history to: traversal_history.json" << std::endl;
    vector_output("traversal_history.json", traversal_history);

}

template<typename T>
std::string to_string(const T& t) {
    return std::to_string(t);
}

template <>
std::string to_string<dev::Coordinate>(const dev::Coordinate& c) {
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