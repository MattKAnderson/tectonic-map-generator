#pragma once
#include <algorithm>
#include <iostream>
#include <memory>
#include <queue>
#include <vector>
#include <random>
#include <unordered_map>
#include <Coordinate.hpp>
#include <geometry.hpp>


struct Node {
    RealCoordinate coord;
    std::vector<Node*> edges;
    Node(const RealCoordinate& coord): coord(coord) {};
    bool operator==(const Node& other) const { return coord == other.coord; }
};


namespace std {
template <>
struct hash<Node> {
    size_t operator()(const Node& n) const {
        return hash<RealCoordinate>()(n.coord);
    }
};   
}


struct LineSegment {
    RealCoordinate a, b;
};


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
    std::vector<Node*> consume_vertices();

private:
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    std::vector<RealCoordinate> seeds;
    std::vector<Node*> vertices;
    std::vector<Node*> regions;
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


class BoundaryRay {
public:
    BoundaryRay() {}
    ~BoundaryRay();
    BoundaryRay(const BoundaryRay& other);
    BoundaryRay& operator=(const BoundaryRay& other);
    void add_vertex(const RealCoordinate v);
    RealCoordinate* a = nullptr;
    RealCoordinate* b = nullptr;
};


class BeachLineItem {
public:
    BeachLineItem* left = nullptr;
    BeachLineItem* right = nullptr;
    BeachLineItem* parent = nullptr;
    BeachLineItem() {}
    BeachLineItem(const RealCoordinate& coord);
    BeachLineItem(const RealCoordinate& coord, BoundaryRay* ray);
    BeachLineItem(const BeachLineItem&) = delete;
    ~BeachLineItem();
    BeachLineItem& operator=(const BeachLineItem&) = delete;
    void set_coord(const RealCoordinate& coord);
    RealCoordinate& coord(); 
    bool has_coord();
    void publish_vertex(const RealCoordinate& coord);
    BoundaryRay* associated_ray = nullptr;
private:
    RealCoordinate* coord_ = nullptr; // focus for regions
};


// TODO: refactor and update the destructor to delete all the used memory
class BeachLine {
public:
    BeachLine(const RealCoordinate& first_p); 
    ~BeachLine();
    BeachLineItem* find_intersected_region(const RealCoordinate& loc);
    BeachLineItem* get_upper_edge(BeachLineItem* node);
    BeachLineItem* get_lower_edge(BeachLineItem* node);
    BeachLineItem* get_upper_region(BeachLineItem* node);
    BeachLineItem* get_lower_region(BeachLineItem* node);
    BeachLineItem* split_region(
        BeachLineItem* region, const RealCoordinate& focus
    );
    BeachLineItem* close_region(
        BeachLineItem* region, const RealCoordinate& coord
    );
    bool is_behind(double directrix, const RealCoordinate& loc);
    void flush();
    std::vector<Node*> vertex_graph();
    //std::vector<Node*> consume_vertex_graph();
    //std::vector<Node*> consume_region_graph();
private:
    BeachLineItem* head = nullptr;
    //std::vector<BeachLineItem*> finished_regions;
    std::vector<BoundaryRay*> rays;
    //std::unordered_map<RealCoordinate, Node*> vertex_map;
    //std::vector<Node*> vertex_graph_;
    //std::unordered_map<RealCoordinate, Node*> region_map;
    //std::vector<Node*> region_graph;
    BeachLineItem* add_subtree(
        BeachLineItem* region, const RealCoordinate& focus
    );
    void record_edge(const RealCoordinate& a, const RealCoordinate& b);
    void link_regions(
        const RealCoordinate& a, const RealCoordinate& b,
        const RealCoordinate& c
    ); 
};


struct FortunesAlgoEvent {
    RealCoordinate coord;
    RealCoordinate* intersect_point = nullptr;
    BeachLineItem* associated_region = nullptr;
    FortunesAlgoEvent(RealCoordinate& coord);
    FortunesAlgoEvent(RealCoordinate& coord, RealCoordinate& intersect);
    FortunesAlgoEvent(const FortunesAlgoEvent& other);
    ~FortunesAlgoEvent();
    FortunesAlgoEvent& operator=(const FortunesAlgoEvent& other);
    bool operator< (const FortunesAlgoEvent& other) const;
    bool operator> (const FortunesAlgoEvent& other) const;
    bool operator== (const FortunesAlgoEvent& other) const;
    friend std::ostream& operator<<(std::ostream& os, const FortunesAlgoEvent& e);
};


inline BoundaryRay::~BoundaryRay() {
    delete a;
    delete b;
}


inline BoundaryRay::BoundaryRay(const BoundaryRay& other) {
    if (other.a != nullptr) {
        a = new RealCoordinate(*other.a);
    }
    if (other.b != nullptr) {
        b = new RealCoordinate(*other.b);
    }
}


inline BoundaryRay& BoundaryRay::operator=(const BoundaryRay& other) {
    if (this == &other) {
        return *this;
    }
    BoundaryRay ray(other);
    std::swap(this->a, ray.a);
    std::swap(this->b, ray.b);
}


inline FortunesAlgoEvent::FortunesAlgoEvent(RealCoordinate& coord): 
    coord(coord) {}


inline FortunesAlgoEvent::FortunesAlgoEvent(RealCoordinate& coord, RealCoordinate& intersect): 
    coord(coord), intersect_point(new RealCoordinate(intersect)) {}
    

inline FortunesAlgoEvent::FortunesAlgoEvent(const FortunesAlgoEvent& other) {
    coord = other.coord;
    associated_region = other.associated_region;
    if (other.intersect_point) {
        intersect_point = new RealCoordinate(*other.intersect_point);
    }
}
    

inline FortunesAlgoEvent::~FortunesAlgoEvent() {
    delete intersect_point;
}
    

inline FortunesAlgoEvent& FortunesAlgoEvent::operator=(
    const FortunesAlgoEvent& other
) {
    if (this != &other) {
        coord = other.coord;
        associated_region = other.associated_region;
        intersect_point = nullptr;
        if (other.intersect_point) {
            intersect_point = new RealCoordinate(*other.intersect_point);
        }
    }
    return *this;
}
    

inline bool FortunesAlgoEvent::operator< (
    const FortunesAlgoEvent& other
) const {
    return (
        coord.x < other.coord.x 
        || (coord.x == other.coord.x && coord.y > other.coord.y)
    );
}
    

inline bool FortunesAlgoEvent::operator> (
    const FortunesAlgoEvent& other
) const {
    return (
        coord.x > other.coord.x
        || coord.x == other.coord.x && coord.y > other.coord.y
    );
}


inline bool FortunesAlgoEvent::operator== (
    const FortunesAlgoEvent& other
) const {
    return coord.x == other.coord.x && coord.y == other.coord.y;
}
    

inline std::ostream& operator<<(std::ostream& os, const FortunesAlgoEvent& e) {
    os << "coord: (" << e.coord.x << ", " << e.coord.y << ")";
    if (e.intersect_point == nullptr) {
        os << ", intersect: nullptr";
    }
    else {
        os << ", intersect: (" << e.intersect_point->x << ", " << e.intersect_point->y << ")";
    }
    if (e.associated_region == nullptr) {
        os << ", associated_region: nullptr";
    }
    else {
        os << ", associated_region: exists";
    }
    return os;
}


inline BeachLineItem::BeachLineItem(const RealCoordinate& coord): 
    coord_(new RealCoordinate(coord)) {};
    

inline BeachLineItem::BeachLineItem(const RealCoordinate& coord, BoundaryRay* ray):
    coord_(new RealCoordinate(coord)), associated_ray(ray) {};


inline BeachLineItem::~BeachLineItem() {
    delete coord_;
}


inline void BeachLineItem::set_coord(const RealCoordinate& coord) {
    delete coord_;
    coord_ = new RealCoordinate(coord);
}


inline RealCoordinate& BeachLineItem::coord() {
    return *coord_;
}


inline bool BeachLineItem::has_coord() {
    return coord_ != nullptr;
}


inline void BeachLineItem::publish_vertex(const RealCoordinate& v) {
    if (associated_ray != nullptr) {
        associated_ray->add_vertex(v);
    }
    else {
        std::cout << "Attempting to publish vertex, but there is no "
                  << "associated ray!" << std::endl;
    }
}


inline BeachLine::BeachLine(const RealCoordinate& first_p): 
    head(new BeachLineItem(first_p)) {}; 


inline BeachLine::~BeachLine() {
    if (head == nullptr) {
        return;
    }
    std::vector<BeachLineItem*> item_stack;
    item_stack.push_back(head);
    while (!item_stack.empty()) {
        BeachLineItem* top = item_stack.back(); 
        item_stack.pop_back();
        if (top->left != nullptr) {
            item_stack.push_back(top->left);
        }
        if (top->right != nullptr) {
            item_stack.push_back(top->right);
        }
        delete top;
    }
}

