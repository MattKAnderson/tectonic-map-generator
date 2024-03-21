#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <utility>
#include <timeit.hpp>

#include "tectonicSimulation.hpp"



template <typename T>
void matrix_output(std::string filename, std::vector<std::vector<T>>& mat);
void time_plate_assignment();
void check_plate_areas(tectonicSimulation& sim);

int main() {
    
    tectonicSimulation sim(16, 1024, 1024);

    sim.plate_assignment();

    auto assignment_result = sim.plate_idmap();
    std::string filename = "plate_partitioning.json";
    matrix_output<int>(filename, assignment_result);
    
}

void check_plate_areas(tectonicSimulation& sim) {
    auto plate_areas = sim.plate_areas();
    std::cout << "plate areas length: " << plate_areas.size() << "\n";
    auto cmp = [](std::pair<int, int>& a, std::pair<int, int>& b) { return a.second < b.second; };
    std::sort(plate_areas.begin(), plate_areas.end(), cmp);
    double total_area = sim.get_sizex() * sim.get_sizey();
    for (auto pa : plate_areas) {
        std::cout << "id:" << std::setw(4) << pa.first << "  area: " << std::setw(8) << pa.second;
        std::cout << "  percent total: " << std::setw(4) << std::setprecision(2) << pa.second / total_area << "\n";  
    }
}

void time_plate_assignment() {
    std::vector<Vector2D<int>> sim_sizes = {
        {200, 200}, {400, 400}, {800, 800}, {1600, 1600}, {3200, 3200}
    };
    
    for (auto& size : sim_sizes) {
        tectonicSimulation sim(10, size.x, size.y);
        timeIt::print([](tectonicSimulation& s) { s.plate_assignment(); }, sim);
    }

    auto timestats = timeIt::stats(10, []() { tectonicSimulation s(10, 1000, 1000); s.plate_assignment(); });

    std::cout << "min:    " << timestats["min"] << " ms\n";
    std::cout << "max:    " << timestats["max"] << " ms\n";
    std::cout << "mean:   " << timestats["mean"] << " ms\n";
    std::cout << "center: " << timestats["center"] << " ms\n";    
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

