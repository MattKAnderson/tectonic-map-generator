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


int main() {

    int seed = 1213468;
    int xsize = 4096;
    int ysize = xsize;
    int nseeds = 25000;

    VoronoiDiagram voronoi_diagram(seed);
    auto t1 = std::chrono::high_resolution_clock::now();
    voronoi_diagram.generate(xsize, ysize, nseeds);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "finished generate function" << std::endl;
    std::cout << "Took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0 << " ms\n";

    std::vector<Node*> vgraph = voronoi_diagram.consume_vertices();
    std::vector<std::pair<RealCoordinate, RealCoordinate>> line_segments;
    std::unordered_set<Node*> visited;
    std::vector<Node*> stack;
    stack.push_back(vgraph[0]);
    //std::cout << "Starting while loop" << std::endl;
    int count = 0;
    while (stack.size() != 0) {
        Node* node = stack.back(); stack.pop_back();
        visited.insert(node);
        //std::cout << "num edges: " << node->edges.size() << std::endl;
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
    std::vector<RealCoordinate> seeds = voronoi_diagram.get_seeds();
    std::cout << "Outputting results" << std::endl;
    vector_output("line_segments.json", line_segments);
    vector_output("diagram_seeds.json", seeds);

}