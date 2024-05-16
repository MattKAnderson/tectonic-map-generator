#include <chrono>
#include <iostream>
#include <limits>
#include <unordered_set>
#include <voronoi_diagram.hpp>
#include <FileIO.hpp>

template<>
std::string to_json<std::pair<RealCoordinate, RealCoordinate>>(
    const std::pair<RealCoordinate, RealCoordinate>& p
) {
    return "[[" + to_json(p.first.x) + ", " + to_json(p.first.y) + "], ["
           + to_json(p.second.x) + ", " + to_json(p.second.y) + "]]";
}


std::vector<std::pair<RealCoordinate, RealCoordinate>> graph_to_line_segments(
    std::vector<VertexNode*> graph
) {
    std::cout << "Graph size: " << graph.size() << std::endl;
    std::vector<std::pair<RealCoordinate, RealCoordinate>> line_segments;
    std::unordered_set<VertexNode*> visited;
    std::vector<VertexNode*> stack;
    stack.push_back(graph[0]);
    while (stack.size() != 0) {
        VertexNode* node = stack.back(); stack.pop_back();
        visited.insert(node);
        for (VertexNode* edge : node->connected) {
            if (visited.find(edge) != visited.end()) {
                continue;
            }
            if (edge == nullptr) {
                std::cout << "Null edge detected" << std::endl;
            }
            line_segments.emplace_back(node->coord, edge->coord);
            stack.push_back(edge);
        }
    }
    std::cout << "Number of line segments: " << line_segments.size() << std::endl;
    return line_segments;
}

void check_vertex_graph(std::vector<VertexNode*>& graph);

int main() {

    int seed = 182310;
    int xsize = 4096;
    int ysize = xsize;
    int nseeds = 10000000;

    VoronoiDiagram voronoi_diagram(seed);
    auto t1 = std::chrono::high_resolution_clock::now();
    voronoi_diagram.generate(xsize, ysize, nseeds);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "finished generate function. ";
    std::cout << "Took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0 << " ms\n";

    std::vector<RealCoordinate> seeds = voronoi_diagram.get_seeds();
    std::vector<VertexNode*> vgraph = voronoi_diagram.consume_vertices();
    std::vector<RegionNode*> rgraph = voronoi_diagram.consume_region_graph();
    std::cout << "Finished getting vertices and region graphs" << std::endl;
    check_vertex_graph(vgraph);
    typedef std::vector<std::pair<RealCoordinate, RealCoordinate>> ls_vector;
    ls_vector vertex_line_segments = graph_to_line_segments(vgraph);
    //ls_vector region_adjacency_segments = graph_to_line_segments(rgraph);
    
    std::cout << "Writing results" << std::endl;
    vector_output("vertex_line_segments.json", vertex_line_segments);
    //vector_output(
    //    "region_adjacency_line_segments.json", region_adjacency_segments
    //);
    vector_output("diagram_seeds.json", seeds);
    std::cout << "Done" << std::endl;
}

void check_vertex_graph(std::vector<VertexNode*>& graph) {

    std::vector<VertexNode*> vertices_with_improper_edge_counts;

    double xmin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymin = std::numeric_limits<double>::max();
    double ymax = std::numeric_limits<double>::min();
    for (VertexNode* vn : graph) {
        xmin = std::min(xmin, vn->coord.x);
        xmax = std::max(xmax, vn->coord.y);
        ymin = std::min(ymin, vn->coord.y);
        ymax = std::max(ymax, vn->coord.y);
    }
    for (VertexNode* vn : graph) {
        RealCoordinate& c = vn->coord;
        if (c.x == xmin || c.x == xmax || c.y == ymin || c.y == ymax) {
            if (vn->connected.size() != 2 && vn->connected.size() != 3) {
                vertices_with_improper_edge_counts.push_back(vn);
            }
        }
        else if (vn->connected.size() != 3) {
            vertices_with_improper_edge_counts.push_back(vn);
        }
    }
    int violation_count = vertices_with_improper_edge_counts.size();
    if (violation_count > 0) {
        std::cout << "WARNING: vertices with improper edge counts detected!\n"
                  << "         Number of violations " << violation_count << "\n";
    }
}
