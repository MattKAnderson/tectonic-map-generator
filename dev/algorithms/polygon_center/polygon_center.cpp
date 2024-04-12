#include <polygon_center.hpp>


double toroidal_distance(Coordinate c1, Coordinate c2, int xmax, int ymax) {

    double dx = std::abs(c2.x - c1.x);
    if (dx > xmax / 2) {
        dx = xmax - dx;
    }
    double dy = std::abs(c2.y - c1.y);
    if (dy > ymax / 2) {
        dy = ymax - dy;
    }
    return std::sqrt(dx * dx + dy * dy);
}


double distance_sum(
    Coordinate c, std::vector<Coordinate>& points, int xmax, int ymax
) {
    double sum = 0.0;
    for (int i = 0; i < points.size(); i++) {
        sum += toroidal_distance(c, points[i], xmax, ymax);
    }
    return sum;
}


// The chance to go ahead with the new sample based on the distances
double chance(double dist, double last_dist, double k) {
    return std::exp(k * (last_dist - dist));
}


std::vector<Coordinate> adjacent_grid_cells(
    Coordinate& c, int xmax, int ymax
) {
    std::vector<Coordinate> adjacents;
    if (c.x > 0) {
        adjacents.emplace_back(c.x - 1, c.y);
    }
    if (c.x + 1 < xmax) {
        adjacents.emplace_back(c.x + 1, c.y);
    }
    if (c.y > 0) {
        adjacents.emplace_back(c.x, c.y - 1);
    }
    if (c.y + 1 < ymax) {
        adjacents.emplace_back(c.x, c.y + 1);
    }
    return adjacents;
}


Coordinate find_toroidal_centroid(
    std::vector<Coordinate>& points, int xmax, int ymax, int iters, double k
) {
    // get this as a parameter later
    std::mt19937_64 rng(129129);
    std::uniform_real_distribution<double> die(0, 1.0);

    if (points.size() == 0) {
        std::cout << "No points given" << std::endl;
        return {0, 0};
    }
    else if (points.size() == 1) {
        std::cout << "Only 1 point was given" << std::endl;
        return points[0];
    }

    // initialize by taking midpoints between existing points
    int start_x = (points[0].x + points[1].x) / 2;
    int start_y = (points[0].y + points[1].y) / 2;
    Coordinate candidate(start_x, start_y);
    double total_dist = distance_sum(candidate, points, xmax, ymax);

    // dumb implementation to start that can only move one grid
    // cell at a time
    // just iterate a certain number of times. Can determine exit 
    // conditions later, can save agrgegate distances when this
    // gets converted to a class later... 

    for (int i = 0; i < iters; i++) {
        std::vector<Coordinate> candidates;
        candidates = adjacent_grid_cells(candidate, xmax, ymax);
        
        std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
        Coordinate prospect = candidates[dist(rng)];

        double prospect_dist = distance_sum(prospect, points, xmax, ymax);
        if (chance(prospect_dist, total_dist, k) > die(rng)) {
            //std::cout << "prospect dist: " << prospect_dist << "  total dist: " << total_dist
            //          << "  chance: " << chance(prospect_dist, total_dist, k) << std::endl; 
            candidate = prospect;
            total_dist = prospect_dist;
        }
    }

    return candidate;
}