#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <Coordinate.hpp>
#include <Vector2D.hpp>


class CentroidFinder {
public:
    CentroidFinder(int seed);
    Coordinate toroidal_centroid(
        std::vector<Coordinate>& points, int xmax, int ymax, 
        int iters, int k, bool collect_metrics = false
    );
    std::vector<Coordinate> get_visit_history();
    std::vector<double> get_total_dist_history();

private:
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> die{0.0, 1.0};
    std::vector<Coordinate> visit_history;
    std::vector<double> total_dist_history;

    Coordinate toroidal_mc(
        std::vector<Coordinate>& points, int xmax, int ymax, int iters, int k,
        Coordinate point
    );
    Coordinate toroidal_mc_record_metrics(
        std::vector<Coordinate>& points, int xmax, int ymax, int iters, int k,
        Coordinate point
    );
};

