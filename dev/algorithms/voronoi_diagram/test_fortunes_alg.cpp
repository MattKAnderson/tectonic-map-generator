#include <chrono>
#include <iostream>
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
    std::vector<Node*> graph
) {
    std::cout << "Graph size: " << graph.size() << std::endl;
    std::vector<std::pair<RealCoordinate, RealCoordinate>> line_segments;
    std::unordered_set<Node*> visited;
    std::vector<Node*> stack;
    stack.push_back(graph[0]);
    while (stack.size() != 0) {
        Node* node = stack.back(); stack.pop_back();
        visited.insert(node);
        for (Node* edge : node->edges) {
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
    return line_segments;
}


int main() {

    int seed = 1213468;
    int xsize = 4096;
    int ysize = xsize;
    int nseeds = 1000;

    VoronoiDiagram voronoi_diagram(seed);
    auto t1 = std::chrono::high_resolution_clock::now();
    voronoi_diagram.generate(xsize, ysize, nseeds);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "finished generate function" << std::endl;
    std::cout << "Took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0 << " ms\n";

    std::vector<RealCoordinate> seeds = voronoi_diagram.get_seeds();
    std::vector<Node*> vgraph = voronoi_diagram.consume_vertices();
    std::vector<Node*> rgraph = voronoi_diagram.consume_region_graph();
    std::cout << "Finished getting vertices and region graphs" << std::endl;
    typedef std::vector<std::pair<RealCoordinate, RealCoordinate>> ls_vector;
    ls_vector vertex_line_segments = graph_to_line_segments(vgraph);
    ls_vector region_adjacency_segments = graph_to_line_segments(rgraph);
    
    std::cout << "Writing results" << std::endl;
    vector_output("vertex_line_segments.json", vertex_line_segments);
    vector_output(
        "region_adjacency_line_segments.json", region_adjacency_segments
    );
    vector_output("diagram_seeds.json", seeds);
    std::cout << "Done" << std::endl;
}