#pragma once
#include <algorithm>
#include <iostream>
#include <memory>
#include <queue>
#include <vector>
#include <random>
#include <unordered_set>
#include <Coordinate.hpp>
#include <geometry.hpp>


struct Node {
    RealCoordinate coord;
    std::vector<Node*> edges;
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
    std::vector<Node> delauney_triangulation;
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


struct FortunesAlgoEvent {
    RealCoordinate coord;
    RealCoordinate* intersect_point = nullptr;
    FortunesAlgoEvent(RealCoordinate& coord);
    FortunesAlgoEvent(RealCoordinate& coord, RealCoordinate& intersect);
    FortunesAlgoEvent(const FortunesAlgoEvent& other);
    ~FortunesAlgoEvent();
    FortunesAlgoEvent& operator=(const FortunesAlgoEvent& other);
    bool operator< (const FortunesAlgoEvent& other) const;
    bool operator> (const FortunesAlgoEvent& other) const;
    bool operator== (const FortunesAlgoEvent& other) const;
};


struct BeachLineItem {
    RealCoordinate* coord = nullptr;
    BeachLineItem* left = nullptr;
    BeachLineItem* right = nullptr;
    BeachLineItem* parent = nullptr;
    BeachLineItem() {};
    BeachLineItem(const RealCoordinate& coord);
    ~BeachLineItem();
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
    BeachLineItem* close_region(BeachLineItem* region);
    BeachLineItem* get_region_2_upper(BeachLineItem* region);
    BeachLineItem* get_region_2_lower(BeachLineItem* region);
    std::vector<Node> get_vertex_graph();
    std::vector<Node> get_region_graph();
private:
    BeachLineItem* head = nullptr;
    std::vector<Node> vertices;
    std::vector<Node> regions;
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
    coord(new RealCoordinate(coord)) {};
    
inline BeachLineItem::~BeachLineItem() { 
    delete coord; 
};

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

