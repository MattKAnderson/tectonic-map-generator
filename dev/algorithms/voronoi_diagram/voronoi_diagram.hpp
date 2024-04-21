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

private:
    std::mt19937_64 rng;
    std::vector<std::vector<int>> diagram;
    std::vector<Node> vertices;
    std::vector<Node> regions;
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


struct BeachLineItem {
    RealCoordinate coord; // focus for regions, base for edges
    BeachLineItem* left = nullptr;
    BeachLineItem* right = nullptr;
    BeachLineItem* parent = nullptr;
    BeachLineItem(const RealCoordinate& coord);
    BeachLineItem(const BeachLineItem&) = delete;
    BeachLineItem& operator=(const BeachLineItem&) = delete;
};


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
    void flush();
    std::vector<Node> get_vertex_graph();
    std::vector<Node> get_region_graph();
private:
    BeachLineItem* head = nullptr;
    std::vector<BeachLineItem*> finished_regions;
    std::vector<LineSegment> line_segments;
    std::unordered_map<RealCoordinate, Node*> vertex_map;
    std::vector<Node> vertex_graph;
    std::unordered_map<RealCoordinate, Node*> region_map;
    std::vector<Node> region_graph;
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
};


inline FortunesAlgoEvent::FortunesAlgoEvent(RealCoordinate& coord): 
    coord(coord) {}

inline FortunesAlgoEvent::FortunesAlgoEvent(RealCoordinate& coord, RealCoordinate& intersect): 
    coord(coord), intersect_point(new RealCoordinate(intersect)) {}
    
inline FortunesAlgoEvent::FortunesAlgoEvent(const FortunesAlgoEvent& other) {
    coord = other.coord;
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

inline BeachLineItem::BeachLineItem(const RealCoordinate& coord): 
    coord(coord) {};
    
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

