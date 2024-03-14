#include <vector>
#include <Vector2D.hpp>

struct edge {
    edge(Vector2D<int> p1, Vector2D<int> p2, int group_a, int group_b): p1(p1), p2(p2), group_a(group_a), group_b(group_b) {};
    Vector2D<int> p1;
    Vector2D<int> p2;
    int group_a;
    int group_b;
};

struct hash_int_Vector2D {
    size_t operator()(const Vector2D<int>& p) const {
        return std::hash<long long>{}(p.x * (INT_MAX + static_cast<long long>(1)) + p.y);
    }
};

std::vector<std::vector<int>> return_boundaries(std::vector<std::vector<int>>& group_map);

std::vector<edge> boundary_outlines(std::vector<std::vector<int>>& group_map);

std::vector<std::pair<int, edge>> label_edges(std::vector<edge>& edges);