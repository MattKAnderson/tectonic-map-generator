#include <polygon_center.hpp>


double sq_toroidal_distance(Coordinate c1, Coordinate c2, int xmax, int ymax) {

    double dx = std::abs(c2.x - c1.x);
    if (dx > xmax / 2) {
        dx = xmax - dx;
    }
    double dy = std::abs(c2.y - c1.y);
    if (dy > ymax / 2) {
        dy = ymax - dy;
    }
    return dx * dx + dy * dy;
}


double toroidal_distance_sum(
    Coordinate c, std::vector<Coordinate>& points, int xmax, int ymax
) {
    double sum = 0.0;
    for (int i = 0; i < points.size(); i++) {
        sum += sq_toroidal_distance(c, points[i], xmax, ymax);
    }
    return sum;
}


// The chance to go ahead with the new sample based on the distances
double chance(double dist, double last_dist, double k) {
    return std::exp(k * (last_dist - dist));
}


std::vector<Coordinate> toroidal_adjacent_grid_cells(
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

CentroidFinder::CentroidFinder(int seed) {
    rng = std::mt19937_64(seed);
}

Coordinate CentroidFinder::toroidal_centroid(
    std::vector<Coordinate>& points, int xmax, int ymax, int iters, int k, 
    bool collect_metrics
) {
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

    visit_history = {};
    cost_history = {};
    if (collect_metrics) {
        return toroidal_mc_record_metrics(
            points, xmax, ymax, iters, k, candidate
        );
    }
    else {
        return toroidal_mc(points, xmax, ymax, iters, k, candidate);
    }
}

std::vector<Coordinate> CentroidFinder::get_visit_history() {
    return visit_history;
}

std::vector<double> CentroidFinder::get_cost_history() {
    return cost_history;
}

Coordinate CentroidFinder::toroidal_mc(
        std::vector<Coordinate>& points, int xmax, int ymax, int iters, int k,
        Coordinate point
) {
    double total_dist = toroidal_distance_sum(point, points, xmax, ymax);
    
    for (int i = 0; i < iters; i++) {
        std::vector<Coordinate> candidates;
        candidates = toroidal_adjacent_grid_cells(point, xmax, ymax);
        
        std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
        Coordinate candidate = candidates[dist(rng)];

        double candidate_dist = toroidal_distance_sum(
            candidate, points, xmax, ymax
        );
        if (chance(candidate_dist, total_dist, k) > die(rng)) {
            point = candidate;
            total_dist = candidate_dist;
        }
    }

    return point;
}

Coordinate CentroidFinder::toroidal_mc_record_metrics(
        std::vector<Coordinate>& points, int xmax, int ymax, int iters, int k,
        Coordinate point
) {
    double total_dist = toroidal_distance_sum(point, points, xmax, ymax);
    visit_history.push_back(point);
    cost_history.push_back(total_dist);

    for (int i = 0; i < iters; i++) {
        std::vector<Coordinate> candidates;
        candidates = toroidal_adjacent_grid_cells(point, xmax, ymax);
        
        std::uniform_int_distribution<int> dist(0, candidates.size() - 1);
        Coordinate candidate = candidates[dist(rng)];

        double candidate_dist = toroidal_distance_sum(
            candidate, points, xmax, ymax
        );
        if (chance(candidate_dist, total_dist, k) > die(rng)) {
            point = candidate;
            total_dist = candidate_dist;
            visit_history.push_back(point);
            cost_history.push_back(total_dist);
        }
    }

    return point;
}