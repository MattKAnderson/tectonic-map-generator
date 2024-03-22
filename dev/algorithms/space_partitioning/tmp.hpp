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

// TODO - produce a Readme.md for the dev work with some of the results.
//      - explore adding noise to voronoi diagram
//      - explore random voronoi cell aggregation, see if this is an 
//        effective means for varying sizes (the weighted distance
//        approach was not very effective)
//      - explore clustered/non-uniform seed placement, this would 
//        naturally lead to larger areas
//      - explore starting with fewer cells and implementing random
//        partitioning within a cell, i.e take cell after applying noise
//        and start a random walk with a drift vector from a random 
//        point on the boundary of the cell to "carve out" a new area
//      - explore iterative carving, slight change on the random lines 
//        approach, don't start new lines randomly anywhere start at 
//        midpoints between existing junctions with some random adjustment 
//      - snaky portions and pointy edges are the bane of it looking 
//        natural. The cocos and Juan de Fuca plates are thin and wedged
//        between other larger plates on Earth. Maybe there is some post
//        processing that can be done to eliminate portions of regions 
//        that start snaking by breaking the snaking part off into a new 
//        region
//      - some amount of the unbalanced voronoi could be used, especially 
//        if fully contained regions are 'absorbed' into the parent region
//        these are easy to identify by checking for disjoint loops in the
//        boundary graph. An alpha value ~0.05 seems to be good for a little
//        but not too much distortion of the base plot.
namespace dev {

struct Coordinate {
    int x, y;
    bool operator==(const Coordinate& other) const {
        return x == other.x && y == other.y;
    };
};

struct Edge {
    Coordinate start, end;
};

struct OrderedEdge {
    Edge edge;
    int order;
    bool operator<(const OrderedEdge& other) const {
        return order < other.order;
    }
    OrderedEdge(Coordinate start, Coordinate end, int order) :
        edge({start, end}), order(order) {};
};

struct OrderedCoordinate {
    Coordinate coord;
    int order;
    bool operator<(const OrderedCoordinate& other) const {
        return order < other.order;
    }
};

struct CoordinateHash {
    std::size_t operator()(const Coordinate& c) const {
        return std::hash<long long>{}(
            c.x * (INT_MAX + static_cast<long long>(1)) + c.y
        );
    }
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
    void unbalanced_voronoi_partition(int xsize, int ysize, int num_seeds, double alpha);

    std::vector<std::vector<Coordinate>> get_borderlines();
    std::vector<std::vector<int>> get_region_map();
    std::vector<std::vector<int>> get_border_map();
    std::vector<Edge> get_traversal_history();
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
    // A better algo would be Fortune's algorithm, if there's time for that
    void naive_voronoi(int num_seeds);
    void unbalanced_voronoi(int num_seeds, double alpha);
    std::vector<double> normally_distributed_seed_weights(
        int num_seeds, double alpha, double min_weight = 1.0
    );
    std::vector<Coordinate> place_seeds_and_get_coords(int num_seeds);
    Coordinate weighted_closest_seed(
        std::vector<Coordinate>& seeds,
        std::vector<double>& weights,
        Coordinate c
    );
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

    inline int above(int x) {
        return x - 1 < 0 ? xsize_ - 1 : x - 1;
    }
    inline int below(int x) {
        return x + 1 < xsize_ ? x + 1 : 0;
    }
    inline int left(int y) {
        return y - 1 < 0 ? ysize_ - 1 : y - 1;
    }
    inline int right(int y) {
        return y + 1 < ysize_ ? y + 1 : 0;
    }
};

} // namespace dev