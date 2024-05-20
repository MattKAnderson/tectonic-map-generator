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
#include <LineClipper.hpp>

#include <chrono>

// forward dec
namespace Impl { class FortunesAlgorithm; }

struct VertexNode {
    RealCoordinate coord;
    std::vector<VertexNode*> connected;
    VertexNode() { connected.reserve(3); }
    VertexNode(const RealCoordinate& coord): coord(coord) { connected.reserve(3); }
    bool operator==(const VertexNode& other) const { return coord == other.coord; }
};

class VertexGraph {
public:
    VertexGraph(VertexNode* graph, std::vector<VertexNode*>& refs);
    ~VertexGraph();
    VertexNode* get_head();
    std::vector<VertexNode*> get_vertices(); 
private:
    VertexNode* data;
    std::vector<VertexNode*> refs;
};

struct RegionNode {
    std::vector<RealCoordinate> vertices;
    std::vector<RegionNode*> adjacent;
    RegionNode() {}
    RealCoordinate centroid();
};

class RegionGraph {
public: 
    RegionGraph(RegionNode* graph, std::vector<RegionNode*>& refs);
    ~RegionGraph();
    RegionNode* get_head();
    std::vector<RegionNode*> get_regions();
private:
    RegionNode* data;
    std::vector<RegionNode*> refs;
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
    ~VoronoiDiagram() { delete generator; }
    void generate(int xsize, int ysize, int nseeds);
    void generate_from_seeds(
        std::vector<RealCoordinate>& seeds, int xsize, int ysize
    );
    std::vector<std::vector<int>> get_diagram();
    std::vector<RealCoordinate> get_seeds();
    VertexGraph get_vertex_graph();
    RegionGraph get_region_graph();

private:
    int nregions;
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    std::vector<RealCoordinate> seeds;
    Impl::FortunesAlgorithm* generator;

    std::vector<RealCoordinate> generate_seeds(
        int nseeds, int xsize, int ysize
    );
};

inline VertexGraph::VertexGraph(
    VertexNode* data, std::vector<VertexNode*>& refs
): data(data), refs(refs) {}

inline VertexGraph::~VertexGraph() { if (data) { delete[] data; } }

inline RegionGraph::RegionGraph(
    RegionNode* data, std::vector<RegionNode*>& refs
): data(data), refs(refs) {}

inline RegionGraph::~RegionGraph() { if (data) { delete[] data; } }

namespace Impl {

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
    //VertexNode* origin = nullptr;
    int origin_id = -1;
};

struct Arc {
    Arc() {}
    Arc(const RealCoordinate& focus, Region* region): focus(focus), region(region) {};
    RealCoordinate focus;
    Region* region = nullptr;
    bool red = true;
    int event_id = -1;
    Arc* left = nullptr;
    Arc* right = nullptr;
    Arc* parent = nullptr;
    Arc* lower = nullptr;
    Arc* upper = nullptr;
    HalfEdge* upper_edge = nullptr;
    HalfEdge* lower_edge = nullptr;
};

struct Event {
    RealCoordinate coord;
    RealCoordinate intersect_point;
    Arc* associated_arc = nullptr;
    Event(
        const RealCoordinate& coord, const RealCoordinate& intersect, 
        Arc* associated_arc
    );
    Event(const RealCoordinate& coord);
};

class EventManager {
public:
    EventManager() {};
    const Event& get(int id);
    int create(const RealCoordinate& coord);
    int create(
        const RealCoordinate& coord, const RealCoordinate& intersect, 
        Arc* associated_arc
    );
    void remove(int id);

private:
    std::vector<Event> events;
    std::vector<int> available_stack;
};

class EventQueue {
public:
    void set_event_manager(EventManager* em) { this->em = em; }
    void insert(int event_id);
    void remove(int event_id);
    void reserve_space(int nevents);
    void print_ordered_x();
    int consume_next();
    bool empty();
private:
    EventManager* em = nullptr;
    std::vector<int> event_id_heap;
    std::vector<int> id_to_location{8}; // for fast delete
    int lchild(int id);
    int rchild(int id);
    int parent(int id);
    void heapify(int id);
    void up_heapify(int id);
    void down_heapify(int id);
    void swap(int ida, int idb);
    bool compare(int ida, int idb);
    bool compare_event_id(int event_ida, int event_idb);
};

class BeachLine {
public:
    BeachLine();
    BeachLine(const BeachLine& other) = delete;
    BeachLine& operator=(const BeachLine& other);
    ~BeachLine();
    Arc* new_arc(const RealCoordinate& r1, Region* region);
    Arc* find_intersected_arc(const RealCoordinate& c);
    void set_head(Arc* arc);
    Arc* get_head();
    Arc* get_lowest();
    void insert_arc_above(Arc* arc, Arc* new_arc);
    void insert_arc_below(Arc* arc, Arc* new_arc);
    void remove_arc(Arc* arc);
    void reserve(int n);

private:
    const static int POOL_ALLOC_SIZE = 512; 
    Arc* head = nullptr;
    std::vector<Arc*> arc_pools;
    int next_index = 0;
    std::vector<Arc*> available_arcs;
    //std::vector<Arc*> closed_regions;
    void insert_balance(Arc* arc);
    void delete_balance(Arc* arc);
    Arc* insert_local_balance(Arc* arc);
    Arc* rotate_left(Arc* arc);
    Arc* rotate_right(Arc* arc);
    void flip_colors(Arc* arc);
    bool is_red(Arc* arc);
};

class FortunesAlgorithm {
public:
    FortunesAlgorithm() {}
    ~FortunesAlgorithm();
    void compute(std::vector<RealCoordinate>& seeds, double min, double max);
    void crop(double min, double max);
    void bound();
    VertexGraph get_vertex_graph();
    RegionGraph get_region_graph();
private:
    int num_seeds;
    //double min;
    //double max;
    EventManager event_manager;
    std::vector<VertexNode*> vertices;
    std::vector<HalfEdge*> half_edges;
    // std::vector<Region*> region_refs; TODO use this
    HalfEdge* internal_half_edges = nullptr;
    //HalfEdge* boundary_half_edges = nullptr;
    int next_half_edge_index = 0;
    VertexNode* internal_vertices = nullptr;
    //VertexNode* exterior_vertices = nullptr;
    std::vector<HalfEdge*> additional_half_edges;
    std::vector<VertexNode*> additional_vertices;
    int next_vertex_index = 0;
    EventQueue event_queue;
    BeachLine beach_line;
    Region* regions = nullptr;
    int next_region_id = 0;

    HalfEdge* new_interior_edge(Region* region);
    VertexNode* new_interior_vertex(const RealCoordinate& c);
    void initialize();
    void site_event(const RealCoordinate& focus);
    void intersection_event(const Event& event);
    void bound_DCEL();
    void crop_DCEL();
    void connect_DCEL_exterior(
        std::list<std::pair<RealCoordinate, HalfEdge*>>& exterior,
        double xmax, double xmin, double ymax, double ymin
    );
    RealCoordinate clip_infinite_edge(
        HalfEdge*, double xmax, double xmin, double ymax, double ymin
    );
    void clip_to_bbox(double min, double max);
    std::vector<RegionNode*> region_graph_from_regions();
};

inline BeachLine::BeachLine() {
    arc_pools.push_back(new Arc[POOL_ALLOC_SIZE]);
}

inline BeachLine::~BeachLine() {
    for (Arc* pool : arc_pools) {
        delete[] pool;
    }
}

inline void BeachLine::set_head(Arc* arc) {
    head = arc;
}

inline Arc* BeachLine::new_arc(const RealCoordinate& focus, Region* region) {
    Arc* arc = nullptr;
    if (!available_arcs.empty()) {
        arc = available_arcs.back(); available_arcs.pop_back();
        *arc = Arc(focus, region);
    }
    else {
        if (next_index < POOL_ALLOC_SIZE) {
            arc = &arc_pools.back()[next_index++];
        }
        else {
            arc_pools.push_back(new Arc[POOL_ALLOC_SIZE]);
            next_index = 1;
            arc = &arc_pools.back()[0]; // can remove this subscript right
        }
        arc->focus = focus;
        arc->region = region;
    }
    return arc;
}

inline BeachLine& BeachLine::operator=(const BeachLine& other) {
    //head = new Arc(other.head->focus, other.head->region);
    return *this; // this is hacky, TODO: change
}

inline Arc* BeachLine::get_head() {
    return head;
}

inline Event::Event(const RealCoordinate& coord): coord(coord) {}

inline Event::Event(
    const RealCoordinate& coord, const RealCoordinate& intersect,
    Arc* associated_arc
): coord(coord), intersect_point(intersect), 
   associated_arc(associated_arc) {}

inline FortunesAlgorithm::~FortunesAlgorithm() {
    if (internal_half_edges) { delete[] internal_half_edges; }
    if (internal_vertices) { delete[] internal_vertices; }
    if (regions) { delete[] regions; }
    for (VertexNode* vertices : additional_vertices) {
        delete[] vertices;
    }
    for (HalfEdge* edges : additional_half_edges) {
        delete[] edges;
    }
}

inline HalfEdge* FortunesAlgorithm::new_interior_edge(Region* region) {
    HalfEdge* edge = &internal_half_edges[next_half_edge_index++];
    edge->region = region;
    return edge;
}

inline VertexNode* FortunesAlgorithm::new_interior_vertex(
    const RealCoordinate& c
) {
    VertexNode* vertex = &internal_vertices[next_vertex_index++];
    vertex->coord = c;
    return vertex;
}

} // Impl