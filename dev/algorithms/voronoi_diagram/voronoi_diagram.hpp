#pragma once
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>
#include <queue>
#include <vector>
#include <random>
#include <unordered_map>
#include <Coordinate.hpp>
#include <geometry.hpp>


struct VertexNode {
    RealCoordinate coord;
    std::vector<VertexNode*> connected;
    VertexNode(const RealCoordinate& coord): coord(coord) {}
    bool operator==(const VertexNode& other) const { return coord == other.coord; }
};

struct RegionNode {
    std::vector<RealCoordinate> vertices;
    std::vector<RegionNode*> adjacent;
    RegionNode() {}
    RealCoordinate centroid();
};

namespace std {
template <>
struct hash<VertexNode> {
    size_t operator()(const VertexNode& n) const {
        return hash<RealCoordinate>()(n.coord);
    }
};   
} // namespace std

class VoronoiDiagram {
public:
    VoronoiDiagram(int seed = 0);
    void generate(int xsize, int ysize, int nseeds);
    void generate_from_seeds(
        std::vector<RealCoordinate>& seeds, int xsize, int ysize
    );
    void lloyd_iteration();
    void nesting_iteration(int nseeds);
    void nesting_iteration_from_seeds(std::vector<RealCoordinate>& seeds);
    std::vector<std::vector<int>> get_diagram();
    std::vector<RealCoordinate> get_seeds();
    std::vector<VertexNode*> consume_vertices();
    std::vector<RegionNode*> consume_region_graph();

private:
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    std::vector<RealCoordinate> seeds;
    std::vector<VertexNode*> vertices;
    std::vector<RegionNode*> regions;
    int nregions;

    void fortunes_algorithm(
        std::vector<RealCoordinate>& seeds, int xsize, int ysize
    );
    void grid_algorithm(
        std::vector<RealCoordinate>& seeds, int xsize, int ysize
    );
    std::vector<RealCoordinate> generate_seeds(
        int nseeds, int xsize, int ysize
    );
};

namespace Impl {

class BoundaryRay {
public:
    BoundaryRay() {};
    void clip_to_bbox(double xmin, double xmax, double ymin, double ymax);
    VertexNode* v[2] = {nullptr, nullptr}; // left, right vertices
private:
    void clip_infinite_ray(double xmin, double xmax, double ymin, double ymax);
    RealCoordinate* clip_infinite_ray(
        double x_int, double y_int, double x0, double y0, double ymin, double ymax
    );
    void clip_vertex_to_bbox(
        RealCoordinate* v, RealCoordinate* other, double xmin, double xmax, 
        double ymin, double ymax
    );
};

/*
class Region {
public:
    Region(const RealCoordinate& seed): seed(seed) {}
    void add_ray(BoundaryRay* ray);
    void add_adjacent(Region* region);
    void prune_rays();
    void draw_rays_on_bounds(double min, double max);
    std::vector<RealCoordinate> get_bounds();
    RealCoordinate seed;
    std::vector<BoundaryRay*> rays;
    std::vector<Region*> adjacent;
};
*/
struct HalfEdge;

struct Region {
    RealCoordinate seed;
    HalfEdge* an_edge = nullptr;
};

struct HalfEdge {
    Region* region = nullptr;
    HalfEdge* prev = nullptr;
    HalfEdge* next = nullptr;
    HalfEdge* twin = nullptr;
    VertexNode* origin = nullptr;
};

struct Event;

struct Arc {
    Arc(const RealCoordinate& focus, Region* region): focus(focus), region(region) {};
    RealCoordinate focus;
    Region* region = nullptr;
    bool red = true;
    bool active = true;
    Arc* left = nullptr;
    Arc* right = nullptr;
    Arc* parent = nullptr;
    Arc* lower = nullptr;
    Arc* upper = nullptr;
    HalfEdge* upper_edge = nullptr;
    HalfEdge* lower_edge = nullptr;
    //Impl::BoundaryRay* lower_ray = nullptr;
    //Impl::BoundaryRay* upper_ray = nullptr;
};

class BeachLine {
public:
    BeachLine() {}
    BeachLine(const RealCoordinate& r1, Region* region): head(new Arc(r1, region)) {}; 
    BeachLine(const BeachLine& other) = delete;
    BeachLine& operator=(const BeachLine& other);
    ~BeachLine();
    Arc* find_intersected_arc(const RealCoordinate& c);
    void set_head(const RealCoordinate& focus, Region* region);
    Arc* get_head();
    Arc* get_lowest();
    void insert_arc_above(Arc* arc, Arc* new_arc);
    void insert_arc_below(Arc* arc, Arc* new_arc);
    void remove_arc(Arc* arc);
    void reserve(int n);

private:
    Arc* head;
    std::vector<Arc*> closed_regions;
    void insert_balance(Arc* arc);
    void delete_balance(Arc* arc);
    void rotate_left(Arc* arc);
    void rotate_right(Arc* arc);
    void flip_colors(Arc* arc);
    bool is_red(Arc* arc);
};

class FortunesAlgorithm {
public:
    FortunesAlgorithm() {}
    void compute(std::vector<RealCoordinate>& seeds, double min, double max);
    //std::vector<Node*> vertex_graph();
    //std::vector<Node*> region_graph();
    std::vector<VertexNode*> consume_vertex_graph();
    std::vector<RegionNode*> consume_region_graph();
private:
    int num_seeds;
    double min;
    double max;
    std::vector<Impl::BoundaryRay*> rays;
    std::vector<Impl::Region*> regions;
    std::vector<VertexNode*> vertices;
    std::vector<HalfEdge*> half_edges;
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> event_queue;
    BeachLine beach_line;

    void site_event(const RealCoordinate& focus);
    void intersection_event(const Event& event);
    void bound_DCEL();
    void flush_beachline();
    void add_vertices_for_bounds_corners();
    void connect_vertices_on_bounds();
    bool compare_bounds_vertices(VertexNode* va, VertexNode* vb);
    std::vector<RegionNode*> region_graph_from_regions();

    std::vector<RealCoordinate> order_region_vertices(
        std::vector<RealCoordinate>& vertex_pairs,
        const RealCoordinate& focus
    );
    void connect_border_vertices(
        std::vector<RealCoordinate>& vertex_pairs,
        const RealCoordinate& focus
    );
};

struct Event {
    RealCoordinate coord;
    RealCoordinate* intersect_point = nullptr;
    Arc* associated_arc = nullptr;
    Event(
        const RealCoordinate& coord, const RealCoordinate& intersect, 
        Arc* associated_arc
    );
    Event(const RealCoordinate& coord);
    Event(const Event& other);
    ~Event();
    Event& operator=(const Event& other);
    bool operator< (const Event& other) const;
    bool operator> (const Event& other) const;
    bool operator== (const Event& other) const;
};
/*
inline Arc::Arc(const RealCoordinate& focus, Region* region):
    focus(focus), region(region), red(true), active(true),
    left(nullptr), right(nullptr), parent(nullptr), lower(nullptr),
    upper(nullptr), upper_edge(nullptr), lower_edge(nullptr) 
{
        std::cout << "This constructor was called" << std::endl;
        if (left != nullptr) {
            std::cout << "left was not nullptr" << std::endl;
        }
}*/

inline BeachLine::~BeachLine() {
    std::vector<Arc*> stack = {head};
    while (!stack.empty()) {
        Arc* arc = stack.back(); stack.pop_back();
        if (arc->left) {
            stack.push_back(arc->left);
        }
        if (arc->right) {
            stack.push_back(arc->right);
        }
        delete arc;
    }
    for (Arc* arc : closed_regions) {
        delete arc;
    }
}

inline BeachLine& BeachLine::operator=(const BeachLine& other) {
    head = new Arc(other.head->focus, other.head->region);
    return *this;
}

inline Arc* BeachLine::get_head() {
    return head;
}

inline Event::Event(const RealCoordinate& coord): coord(coord) {}

inline Event::Event(
    const RealCoordinate& coord, const RealCoordinate& intersect,
    Arc* associated_arc
): coord(coord), intersect_point(new RealCoordinate(intersect)), 
   associated_arc(associated_arc) {}

inline Event::Event(const Event& other) {
    coord = other.coord;
    associated_arc = other.associated_arc;
    if (other.intersect_point) {
        intersect_point = new RealCoordinate(*other.intersect_point);
    }
}

inline Event::~Event() {
    delete intersect_point;
}

inline Event& Event::operator=(const Event& other) {
    if (this == &other) {
        return *this;
    }
    coord = other.coord;
    associated_arc = other.associated_arc;
    delete intersect_point;
    if (other.intersect_point) {
        intersect_point = new RealCoordinate(*other.intersect_point);
    } 
    else {
        intersect_point = nullptr;
    }
    return *this;
}

inline bool Event::operator<(const Event& other) const {
    return (
        coord.x < other.coord.x
        || (coord.x == other.coord.x && coord.y < other.coord.y)
    );
}

inline bool Event::operator>(const Event& other) const {
    return (
        coord.x > other.coord.x 
        || (coord.x == other.coord.x && coord.y > other.coord.y)
    );
}

inline bool Event::operator==(const Event& other) const {
    return coord.x == other.coord.x && coord.y == other.coord.y;
}
} // Impl