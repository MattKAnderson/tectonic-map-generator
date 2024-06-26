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

std::vector<std::pair<RealCoordinate, RealCoordinate>> region_adjacency_line_segments(
    RegionGraph& graph
) {
    std::vector<std::pair<RealCoordinate, RealCoordinate>> line_segments;
    std::unordered_set<RegionNode*> visited;
    for (RegionNode* region : graph.get_regions()) {
        for (RegionNode* adj_region : region->adjacent) {
            if (visited.find(adj_region) != visited.end()) { continue; }
            line_segments.emplace_back(region->centroid(), adj_region->centroid());
        }
        visited.insert(region);
    }
    return line_segments;
}

std::string to_json(RegionNode* p) {
    std::string centroid = to_json(p->centroid());
    std::string vertices_string = "[";
    for (const RealCoordinate& vertex : p->vertices) {
        vertices_string += to_json(vertex) + ", ";
    }
    vertices_string.pop_back(); 
    vertices_string.pop_back();
    vertices_string += "]";
    return "{\"centroid\": " + centroid + ", \"vertices\": " + vertices_string + "}";
}

int main() {

    int seed = 12490; 
    int xsize = 4096;
    int ysize = xsize;
    int nseeds = 1000;

    VoronoiDiagram voronoi_diagram(seed);
    auto t1 = std::chrono::high_resolution_clock::now();
    voronoi_diagram.generate(xsize, ysize, nseeds);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "finished generate function. ";
    std::cout << "Took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0 << " ms\n";

    std::vector<RealCoordinate> seeds = voronoi_diagram.get_seeds();
    VertexGraph vgraph = voronoi_diagram.get_vertex_graph();
    RegionGraph rgraph = voronoi_diagram.get_region_graph();
    std::cout << "Finished getting vertices and region graphs" << std::endl;
    typedef std::vector<std::pair<RealCoordinate, RealCoordinate>> ls_vector;
    ls_vector vertex_line_segments = graph_to_line_segments(vgraph.get_vertices());
    ls_vector region_adjacency_segments = region_adjacency_line_segments(rgraph);
    std::vector<RegionNode*> regions = rgraph.get_regions();
    std::cout << "Writing results" << std::endl;
    vector_output("vertex_line_segments.json", vertex_line_segments);
    vector_output(
        "region_adjacency_line_segments.json", region_adjacency_segments
    );
    vector_output("regions.json", regions);
    vector_output("diagram_seeds.json", seeds);
    std::cout << "Done" << std::endl;
}