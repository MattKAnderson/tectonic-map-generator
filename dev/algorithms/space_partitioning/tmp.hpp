#include <algorithm>
#include <concepts>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Vector2D.hpp>

#include <fstream>

namespace dev {

struct Coordinate {
    int x, y;
    bool operator==(const Coordinate& other) const {
        return x == other.x && y == other.y;
    };
};

struct OrderedCoordinate {
    Coordinate coord;
    int order;
    bool operator<(const OrderedCoordinate& other) const {
        return order < other.order;
    };
};

struct CoordinateHash {
    std::size_t operator()(const Coordinate& c) const {
        return std::hash<long long>{}(
            c.x * (INT_MAX + static_cast<long long>(1)) + c.y
        );
    };
};

struct GridVertex {
    Coordinate coord;
    GridVertex* edges[4] = {}; // +x, -x, +y, -y
    int visit_number;
    int edge_traversal_number[4] = {};
    GridVertex(Coordinate c, int visit_number) : 
        coord(c), visit_number(visit_number) {};
};

/*
 *  Class for partitioning a grid into regions using the desired method. 
 *  Results stored in the object.  
 */
class SpacePartitioner {
public:
    SpacePartitioner(int seed);

    void seed_growth_partition(int xsize, int ysize, int num_seeds);
    void random_walk_partition(
        int xsize, int ysize, int ncuts, double alpha=0.05, 
        double max_drift_change=0.18
    );
    void voronoi_partition(int xsize, int ysize, int num_seeds);

    std::vector<std::vector<Coordinate>> get_borderlines();
    std::vector<std::vector<int>> get_region_map();
    std::vector<std::vector<int>> get_border_map();
    std::vector<Coordinate> get_traversal_history();
    std::unordered_map<Coordinate, GridVertex*, CoordinateHash> get_graph();

private:
    std::uniform_real_distribution<double> random_angle_ {
        -std::numbers::pi, std::numbers::pi
    };
    std::uniform_real_distribution<double> random_chance_threshold_ {
        0.0, 1.0
    };

    int xsize_;
    int ysize_;
    std::unordered_map<Coordinate, GridVertex*, CoordinateHash> graph_;
    std::vector<std::vector<int>> region_map_;
    std::vector<std::vector<int>> border_map_;
    std::mt19937_64 rng_;
    std::vector<std::vector<Coordinate>> borderlines_;    
    std::uniform_int_distribution<int> random_x_;
    std::uniform_int_distribution<int> random_y_;
    std::normal_distribution<double> random_drift_change_;
    double max_drift_change_;

    void initialize_containers();
    void random_walks_on_graph(int nwalks);
    void convert_borderlines_to_border_map();
    void convert_graph_to_borderlines();
    void label_areas_between_boundaries();
    void convert_region_map_to_graph();
    void ensure_nodes_and_edge(
        Coordinate a, Coordinate b, int edge_a, int edge_b
    );
    void travel_until_intersect(
        Coordinate pos, GridVertex* node, Vector2D<double> drift
    );
    void dfs(
        std::unordered_set<GridVertex*>& visited, 
        GridVertex* node, 
        int line_index, 
        int exhausted_edge
    );
    std::vector<Coordinate> follow_line(
        std::unordered_set<GridVertex*>& visited,
        GridVertex* node,
        int start_edge
    );
    void follow_branching_lines(
        std::unordered_set<GridVertex*>& visited,
        GridVertex* node
    );
    void naive_voronoi(int num_seeds);
    Coordinate closest_seed(
        std::vector<Coordinate>& seed_location,
        Coordinate c
    );
    double euclidean_distance(Coordinate& c1, Coordinate& c2);
    double toroidal_distance(Coordinate& c1, Coordinate& c2);
    std::vector<Coordinate> adjacent_coords(Coordinate coord);
    std::vector<OrderedCoordinate> unassigned_adjacent_grid_cells(
        Coordinate c
    );

    inline int above(int y) {
        return y - 1 < 0 ? ysize_ - 1 : y - 1;
    }
    inline int below(int y) {
        return y + 1 < ysize_ ? y + 1 : 0;
    }
    inline int left(int x) {
        return x - 1 < 0 ? xsize_ - 1 : x - 1;
    }
    inline int right(int x) {
        return x + 1 < xsize_ ? x + 1 : 0;
    }
};

} // namespace dev