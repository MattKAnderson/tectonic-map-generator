#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <timeit.hpp>

#include "boundary.hpp"
#include "tectonicSimulation.hpp"

template <typename T>
void matrix_output(std::string filename, std::vector<std::vector<T>>& mat);
void output_edge_list(std::string filename, std::vector<edge>& edges);
void output_labelled_edges(std::string filename, std::vector<std::pair<int, edge>>& labelled_edges);

int main() {
    tectonicSimulation sim(12, 1024, 1024, true, 9999);

    sim.plate_assignment();

    timeIt::print(
        []() { tectonicSimulation sim(12, 1024, 1024, true, 9999); sim.plate_assignment(); }
    );

    auto plate_assignments = sim.plate_idmap();

    //std::cout << "before return boundaries" << std::endl;
    //auto plate_boundaries = return_boundaries(plate_assignments);

    std::cout << "Writing results\n"; 
    matrix_output<int>("plate_partitioning.json", plate_assignments);
    //matrix_output<int>("plate_boundaries.json", plate_boundaries);

    auto edge_list = boundary_outlines(plate_assignments);
    output_edge_list("boundary_edges.json", edge_list);

    auto labelled_edges = label_edges(edge_list);
    output_labelled_edges("labelled_boundary_edges.json", labelled_edges);
}

template <typename T>
void matrix_output(std::string filename, std::vector<std::vector<T>>& mat) {
    std::ofstream outfile;
    outfile.open(filename);
    outfile << "{\"matrix\": [";
    std::cout << "writing to " << filename << std::endl;
    for (int i = 0; i < mat.size() - 1; i++) {
        outfile << "[";
        for (int k = 0; k < mat[i].size() - 1; k++) {
            outfile << std::to_string(mat[i][k]) << ", ";
        }
        outfile << std::to_string(mat[i].back()) << "], ";
    }
    outfile << "[";
    for (int k = 0; k < mat.back().size() - 1; k++) {
        outfile << std::to_string(mat.back()[k]) << ", ";
    }
    outfile << std::to_string(mat.back().back()) << "]]}";
};

void output_edge_list(std::string filename, std::vector<edge>& edges) {
    std::ofstream outfile;
    outfile.open(filename);
    outfile << " {\"edges\": [";
    for (int i = 0; i < edges.size() - 1; i++) {
        outfile << "[[" << edges[i].p1.x << ", " << edges[i].p1.y << "], [" << edges[i].p2.x << ", " << edges[i].p2.y << "]], ";
    }
    outfile << "[[" << edges.back().p1.x << ", " << edges.back().p1.y << "], [" << edges.back().p2.x << ", " << edges.back().p2.y << "]]]}"; 
}

void output_labelled_edges(std::string filename, std::vector<std::pair<int, edge>>& labelled_edges) {
    std::ofstream outfile;
    outfile.open(filename);
    outfile << "{\"labelled_edges\": [";
    for (int i = 0; i < labelled_edges.size() - 1; i++) {
        auto e = labelled_edges[i].second;
        outfile << "{\"label\": " << labelled_edges[i].first << " , \"vertices\": [[";
        outfile << e.p1.x << ", " << e.p1.y << "], [" << e.p2.x << ", " << e.p2.y << "]]}, ";
    }
    auto e = labelled_edges.back().second;
    outfile << "{\"label\": " << labelled_edges.back().first << " , \"vertices\": [[";
    outfile << e.p1.x << ", " << e.p1.y << "], [" << e.p2.x << ", " << e.p2.y << "]]}]}";
}