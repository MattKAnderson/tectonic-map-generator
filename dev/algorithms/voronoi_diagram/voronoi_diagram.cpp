#include <voronoi_diagram.hpp>


std::vector<Coordinate> generate_seeds(
    int nseeds, int xsize, int ysize, std::mt19937_64& rng
) {
    std::uniform_int_distribution random_x(0, xsize - 1);
    std::uniform_int_distribution random_y(0, ysize - 1);
    std::vector<Coordinate> seeds;
    seeds.reserve(nseeds);
    for (int i = 0; i < nseeds; i++) {
        int x = random_x(rng);
        int y = random_y(rng);
        seeds.emplace_back(x, y);
    }
    return seeds;
}


double euclidean_distance(Coordinate c1, Coordinate c2) {
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    return std::sqrt(dx * dx + dy*dy);
}


double toroidal_distance(
    Coordinate c1, Coordinate c2, int xsize, int ysize
) {
    double dx = std::abs(c1.x - c2.x);
    if (dx > xsize / 2) {
        dx = xsize - dx;
    }
    double dy = std::abs(c1.y - c2.y);
    if (dy > ysize / 2) {
        dy = ysize - dy;
    }
    return std::sqrt(dx * dx + dy * dy);
}


int find_closest_id(Coordinate point, std::vector<Coordinate>& points) {
    int closest_id = 0;
    double shortest_dist = euclidean_distance(points[0], point);
    for (int i = 1; i < points.size(); i++) {
        double dist = euclidean_distance(points[i], point);
        if (dist < shortest_dist) {
            closest_id = i;
            shortest_dist = dist;
        }
    }
    return closest_id;
}


std::vector<std::vector<int>> generate_diagram_from_seeds(
    std::vector<Coordinate>& seeds, int xsize, int ysize
) {
    std::vector<std::vector<int>> diagram(xsize, std::vector<int>(ysize, -1));
    for (int i = 0; i < seeds.size(); ++i) {
        diagram[seeds[i].x][seeds[i].y] = i;
    }
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            if (diagram[i][j] == -1) {
                diagram[i][j] = find_closest_id({i, j}, seeds);
            }
        }
    }
    return diagram;
}


std::vector<Coordinate> compute_voronoi_cell_centroids(
    std::vector<std::vector<int>>& diagram, int ncells
) {
    std::vector<int> pixel_counts(ncells, 0);
    std::vector<int> x_sums(ncells, 0);
    std::vector<int> y_sums(ncells, 0);
    for (int i = 0; i < diagram.size(); ++i) {
        for (int j = 0; j < diagram[i].size(); ++j) {
            int cell = diagram[i][j];
            x_sums[cell] += i;
            y_sums[cell] += j;
            ++pixel_counts[cell];
        }
    }
    std::vector<Coordinate> centroids;
    centroids.reserve(ncells);
    for (int i = 0; i < ncells; i++) {
        if (pixel_counts[i] == 0) {
            continue;
        }
        centroids.emplace_back(
            x_sums[i] / pixel_counts[i],
            y_sums[i] / pixel_counts[i]
        );
    }
    return centroids;
}


VoronoiDiagram::VoronoiDiagram(int seed) {
    rng = std::mt19937_64(seed);
}


void VoronoiDiagram::generate(int xsize, int ysize, int nseeds) {
    nregions = nseeds;
    std::vector<Coordinate> seeds = generate_seeds(nseeds, xsize, ysize, rng);
    diagram = generate_diagram_from_seeds(seeds, xsize, ysize);
}


void VoronoiDiagram::generate_from_seeds(
    std::vector<Coordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    diagram = generate_diagram_from_seeds(seeds, xsize, ysize);
}


void VoronoiDiagram::lloyd_iteration() {
    std::vector<Coordinate> new_seeds;
    new_seeds = compute_voronoi_cell_centroids(diagram, nregions);
    diagram = generate_diagram_from_seeds(
        new_seeds, diagram.size(), diagram[0].size()
    );
}


void VoronoiDiagram::nesting_iteration(int nseeds) {
    std::vector<Coordinate> seeds = generate_seeds(
        nseeds, diagram.size(), diagram[0].size(), rng
    );
    nesting_iteration_from_seeds(seeds);
}


void VoronoiDiagram::nesting_iteration_from_seeds(
    std::vector<Coordinate>& seeds
) {
    std::vector<Coordinate> centroids;
    centroids = compute_voronoi_cell_centroids(diagram, nregions);
    std::vector<int> region_id_map(centroids.size(), -1);
    for (int i = 0; i < centroids.size(); i++) {
        region_id_map[i] = find_closest_id(centroids[i], seeds);
    }
    for (int i = 0; i < diagram.size(); i++) {
        for (int j = 0; j < diagram[0].size(); j++) {
            diagram[i][j] = region_id_map[diagram[i][j]];
        }
    }
    nregions = seeds.size();
}


std::vector<std::vector<int>> VoronoiDiagram::get_diagram() {
    return diagram;
}

