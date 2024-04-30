#pragma once
#include <algorithm>
#include <iostream>
#include <iomanip>
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
    std::vector<Node*> consume_region_graph();

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

namespace Impl {

class BoundaryRay {
public:
    BoundaryRay(const RealCoordinate& r1, const RealCoordinate& r2);
    void clip_to_bbox(double xmin, double xmax, double ymin, double ymax);
    RealCoordinate r1;
    RealCoordinate r2;
    RealCoordinate* v[2] = {nullptr, nullptr}; // left, right vertices
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

class FortunesAlgorithm {
public:
    FortunesAlgorithm() {}
    void compute(std::vector<RealCoordinate>& seeds, double min, double max);
    std::vector<Node*> vertex_graph();
    std::vector<Node*> region_graph();
private:
    int num_seeds;
    double min;
    double max;
    std::vector<BoundaryRay*> rays;
};

struct BeachLineItem {
    BeachLineItem(const RealCoordinate& coord, int direction): 
        coord(coord), direction(direction) {}
    BeachLineItem(const RealCoordinate& coord): coord(coord) {}
    RealCoordinate coord;
    int direction;
    BeachLineItem* left = nullptr;
    BeachLineItem* right = nullptr;
    BeachLineItem* parent = nullptr;
    BoundaryRay* ray = nullptr;

};

struct Event {
    RealCoordinate coord;
    RealCoordinate* intersect_point = nullptr;
    BeachLineItem* associated_region = nullptr;
    Event(
        const RealCoordinate& coord, const RealCoordinate& intersect, 
        BeachLineItem* associated_region
    );
    Event(const RealCoordinate& coord);
    Event(const RealCoordinate& coord, const RealCoordinate& intersect);
    Event(const Event& other);
    ~Event();
    Event& operator=(const Event& other);
    bool operator< (const Event& other) const;
    bool operator> (const Event& other) const;
    bool operator== (const Event& other) const;
};

class BeachLine {
public:
    BeachLine(const RealCoordinate& f1, const RealCoordinate& f2);
    BeachLine(const BeachLine& other) = delete;
    ~BeachLine();
    BeachLineItem* get_upper_edge(BeachLineItem* region);
    BeachLineItem* get_lower_edge(BeachLineItem* region);
    BeachLineItem* get_upper_region(BeachLineItem* edge);
    BeachLineItem* get_lower_region(BeachLineItem* edge);
    BeachLineItem* find_intersected_region(const RealCoordinate& c);
    BeachLineItem* head = nullptr;
    std::vector<BeachLineItem*> closed_regions;
};

inline BoundaryRay::BoundaryRay(
    const RealCoordinate& r1, const RealCoordinate& r2
): r1(r1), r2(r2) {}

inline Event::Event(const RealCoordinate& coord): coord(coord) {}

inline Event::Event(
    const RealCoordinate& coord, const RealCoordinate& intersect,
    BeachLineItem* associated_region
): coord(coord), intersect_point(new RealCoordinate(intersect)), 
   associated_region(associated_region) {}

inline Event::Event(const RealCoordinate& coord, const RealCoordinate& intersect):
coord(coord), intersect_point(new RealCoordinate(intersect)) {}

inline Event::Event(const Event& other) {
    coord = other.coord;
    associated_region = other.associated_region;
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
    associated_region = other.associated_region;
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

inline BeachLine::BeachLine(const RealCoordinate& f1, const RealCoordinate& f2) {
    double eps = 1e-10;
    double x;
    if (f1.x < f2.x + eps && f1.x > f2.x - eps) {
        x = std::numeric_limits<double>::min();
    } 
    else {
        x = parabola_x_from_y(f2.x, f1, f2.y);
    }
    int direction = f2.y > f1.y ? 0 : 1;
    head = new BeachLineItem{{x, f2.y}, direction};
    head->right = new BeachLineItem{{x, f2.y}, (direction ^ 1)};
    head->right->parent = head;
    head->left = new BeachLineItem{f1};
    head->left->parent = head;
    head->right->right = new BeachLineItem{f1};
    head->right->right->parent = head->right;
    head->right->left = new BeachLineItem{f2};
    head->right->left->parent = head->right;
}

inline BeachLine::~BeachLine() {
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
    for (auto* region : closed_regions) {
        delete region;
    }
}
} // namespace Impl

namespace Impl2 {
/*
class BoundaryRay {
public:
    BoundaryRay(const RealCoordinate& r1, const RealCoordinate& r2);
    void clip_to_bbox(double xmin, double xmax, double ymin, double ymax);
    RealCoordinate r1;
    RealCoordinate r2;
    RealCoordinate* v[2] = {nullptr, nullptr}; // left, right vertices
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
*/

class FortunesAlgorithm {
public:
    FortunesAlgorithm() {}
    void compute(std::vector<RealCoordinate>& seeds, double min, double max);
    std::vector<Node*> vertex_graph();
    std::vector<Node*> region_graph();
private:
    int num_seeds;
    double min;
    double max;
    std::vector<Impl::BoundaryRay*> rays;
};

struct Event;

struct Arc {
    RealCoordinate focus;
    Arc* left = nullptr;
    Arc* right = nullptr;
    Arc* parent = nullptr;
    Arc* lower = nullptr;
    Arc* upper = nullptr;
    bool active = true;
    Impl::BoundaryRay* lower_ray = nullptr;
    Impl::BoundaryRay* upper_ray = nullptr;
};

class BeachLine {
public:
    BeachLine(const RealCoordinate& r1): head(new Arc{r1}) {};
    ~BeachLine();
    Arc* find_intersected_arc(const RealCoordinate& c);
    Arc* get_head();
    void insert_arc_above(Arc* arc, Arc* new_arc);
    void insert_arc_below(Arc* arc, Arc* new_arc);
    void remove_arc(Arc* arc);
    void reserve(int n);

private:
    Arc* head;
    std::vector<Arc*> closed_regions;
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
    std::cout << "About to delete closed regions" << std::endl;
    for (Arc* arc : closed_regions) {
        delete arc;
    }
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
}