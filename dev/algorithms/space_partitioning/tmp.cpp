#include "tmp.hpp"
#include <array>

namespace dev {
namespace Impl {
    int find_parent(std::vector<int>& parent_labels, int label);
    void join_labels(
        std::vector<int>& parent_labels, int label_a, int label_b
    );
    std::vector<int> relabel(std::vector<int>& parent_labels);
}


SpacePartitioner::SpacePartitioner(int seed) {
    rng_ = std::mt19937_64(seed);
}

std::vector<std::vector<Coordinate>> SpacePartitioner::get_borderlines() {
    return borderlines_;
}

std::vector<std::vector<int>> SpacePartitioner::get_region_map() {
    return region_map_;
}

std::vector<std::vector<int>> SpacePartitioner::get_border_map() {
    return border_map_;
}

std::unordered_map<Coordinate, GridVertex*, CoordinateHash> SpacePartitioner::
get_graph() {
    return graph_;
}

std::vector<Coordinate> SpacePartitioner::get_traversal_history() {
    std::vector<OrderedCoordinate> coords;
    for (auto& [coord, vertex] : graph_) {
        coords.emplace_back(coord, vertex->visit_number);
    }
    std::sort(coords.begin(), coords.end());
    std::vector<Coordinate> ordered_coords;
    ordered_coords.reserve(coords.size());
    std::transform(
        coords.begin(), coords.end(), std::back_inserter(ordered_coords),
        [](OrderedCoordinate& oc) { return oc.coord; }
    );
    return ordered_coords;
}

void SpacePartitioner::seed_growth_partition(
    int xsize, int ysize, int num_seeds
) {
    xsize_ = xsize;
    ysize_ = ysize;
    random_x_ = std::uniform_int_distribution(0, xsize - 1);
    random_y_ = std::uniform_int_distribution(0, ysize - 1);
    initialize_containers();
    // TODO: continue 
    std::cout << "WARNING: seed_growth_partition method not implemented" 
              << std::endl;
}

void SpacePartitioner::random_walk_partition(
    int xsize, int ysize, int ncuts, double alpha, double max_drift_change
) {
    xsize_ = xsize;
    ysize_ = ysize;
    max_drift_change_ = max_drift_change;
    random_drift_change_ = std::normal_distribution<double>(0.0, alpha);
    random_x_ = std::uniform_int_distribution(0, xsize);
    random_y_ = std::uniform_int_distribution(0, ysize);
    //std::cout << "initializing containers" << std::endl;
    initialize_containers();
    //std::cout << "doing random walk on graph" << std::endl;
    random_walks_on_graph(ncuts);
    //std::cout << "converting graph to borderlines" << std::endl;
    convert_graph_to_borderlines();
    //std::cout << "convert borderlines to border map" << std::endl;
    convert_borderlines_to_border_map();
    //std::cout << "labelling areas between boundaries" << std::endl;
    label_areas_between_boundaries();
    //std::cout << "done" << std::endl;
}

void SpacePartitioner::initialize_containers() {
    region_map_ = std::vector<std::vector<int>>(
        xsize_, std::vector<int>(ysize_, -1)
    );
    border_map_ = std::vector<std::vector<int>>(
        2 * xsize_, std::vector<int>(ysize_, -1)
    );
    borderlines_ = {};
    graph_ = {};
}

void SpacePartitioner::travel_until_intersect(
    Coordinate pos, GridVertex* node, Vector2D<double> drift
) {
    while (true) {
        graph_[pos] = node;
        double absx = std::abs(drift.x);
        double absy = std::abs(drift.y);
        double chance_x = absx / (absx + absy);

        double random_drift_change = random_drift_change_(rng_);
        if (random_drift_change > max_drift_change_) {
            random_drift_change = max_drift_change_;
        }
        drift.rotate(random_drift_change);

        Coordinate new_pos; 
        int edge_index = 0;
        int comp_edge_index = 1;
        if (chance_x > random_chance_threshold_(rng_)) {
            if (drift.x > 0) {
                new_pos = {(pos.x + 1 < xsize_) * (pos.x + 1), pos.y};
            }
            else {
                std::swap(edge_index, comp_edge_index);
                new_pos = {pos.x - 1 + (pos.x - 1 < 0) * xsize_, pos.y};
            }
        }
        else {
            edge_index = 2;
            comp_edge_index = 3;
            if (drift.y > 0) {
                new_pos = {pos.x, (pos.y + 1 < ysize_) * (pos.y + 1)};
            }
            else {
                std::swap(edge_index, comp_edge_index);
                new_pos = {pos.x, pos.y - 1 + (pos.y - 1 < 0) * ysize_};
            }
        }
        
        if (auto search = graph_.find(new_pos); search != graph_.end()) {
            if (search->second->edges[comp_edge_index] == node) {
                // not allowed to re-cross the same edge in reverse
                // can happen when drift is straddling one of the axes
                // back and forth
                continue;
            }
            search->second->edges[comp_edge_index] = node;
            node->edges[edge_index] = search->second;
            break;
        }
        else {
            GridVertex* new_node = new GridVertex(new_pos, graph_.size());
            new_node->edges[comp_edge_index] = node;
            node->edges[edge_index] = new_node;
            node = new_node;
            pos = new_pos;
        }
    }
}

void SpacePartitioner::random_walks_on_graph(int nwalks) {
    for (int i = 0; i < nwalks; i++) {
        Coordinate pos(random_x_(rng_), random_y_(rng_));
        while (graph_.find(pos) != graph_.end()) {
            pos = {random_x_(rng_), random_y_(rng_)};
        }
        GridVertex* node = new GridVertex(pos, graph_.size());
        Vector2D<double> drift(0.0, 1.0);
        drift.rotate(random_angle_(rng_));
        
        travel_until_intersect(pos, node, drift);
        travel_until_intersect(pos, node, -drift);
    }
}

void SpacePartitioner::convert_borderlines_to_border_map() {
    for (int i = 0; i < borderlines_.size(); i++) {
        for (int j = 0; j < borderlines_[i].size() - 1; j++) {
            Coordinate& a = borderlines_[i][j];
            Coordinate& b = borderlines_[i][j + 1];
            int x, y;
            if (a.x != b.x) {
                y = a.y;
                x = std::abs(a.x - b.x) > 1 ? xsize_ - 1 : std::min(a.x, b.x); 
                x = 2 * x + 1;
            }
            else {
                x = 2 * a.x;
                y = std::abs(a.y - b.y) > 1 ? ysize_ - 1 : std::min(a.y, b.y); 
            }
            border_map_[x][y] = i;
        }
    }
}
/*
void SpacePartitioner::convert_graph_to_borderlines() {
    
    // for every edge we can know the edge traversal number
    // this number can be used to know the order of the node 
    // sequence to define the order of points in the line
    // Can do a DFS on the graph to extract the lines, any
    // 3-way or 4-way junction spawns a new label
    // visited nodes coordinate is remembered in a hashmap
    // the initial node will probably be in the middle of a
    // line so constructing this first line can be handled 
    // as a preliminary special case
    // dead-end errors can be checked for
    std::unordered_set<GridVertex*> visited;

    // if a high alpha and max_drift_change values are used then
    // it is likely that the traversals are disjoint on the graph
    // and not all can be reached in a single DFS 
    for (auto& [coord, node] : graph_) {
        if (visited.find(node) == visited.end()) {
            std::cout << "starting new search at node: (" << coord.x << ", " << coord.y << ")\n";
            int start_index = borderlines_.size();
            borderlines_.push_back({}); 
            dfs(visited, node, start_index, -1);
            
            // if the DFS starts at a junction, remove the
            // first list which only contains the start coord
            // else, do the same but join the next two lines
            int edge_count = 0;
            for (int i = 0; i < 4; ++i) {
                if (node->edges[i] != nullptr) { ++edge_count; }
            }
            if (edge_count > 2) {
                borderlines_[start_index] = std::move(borderlines_.back());
            }
            else {
                std::reverse(
                    borderlines_[start_index + 1].begin(),
                    borderlines_[start_index + 1].end()
                );
                //borderlines_[start_index + 1].pop_back();
                borderlines_[start_index + 1].insert(
                    borderlines_[start_index + 1].end(),
                    borderlines_[start_index + 2].begin(),
                    borderlines_[start_index + 2].end()
                );
                borderlines_[start_index + 2] = std::move(borderlines_.back());
                borderlines_.pop_back();
                borderlines_[start_index] = std::move(borderlines_.back());
                borderlines_.pop_back();
            }
        }
    }
}
*/
void SpacePartitioner::convert_graph_to_borderlines() {

    std::unordered_set<GridVertex*> visited;
    for (auto& [c, node] : graph_) {
        if (visited.find(node) != visited.end()) {
            continue;
        }
        visited.insert(node);
        int edge_count = 0;
        std::array<int, 4> edge_ids = {};
        for (int i = 0; i < 4; i++) {
            if (node->edges[i] != nullptr) {
                edge_ids[edge_count++] = i;
            }
        }
        auto v1 = node->edges[edge_ids[0]]->coord;
        auto v2 = node->edges[edge_ids[1]]->coord;

        if (edge_count == 2) {
            auto line1 = follow_line(visited, node, edge_ids[0]);
            auto line2 = follow_line(visited, node, edge_ids[1]);

            follow_branching_lines(visited, graph_[line1.back()]);;
            follow_branching_lines(visited, graph_[line2.back()]);
            std::reverse(line1.begin(), line1.end());
            line1.pop_back();
            line1.reserve(line1.size() + line2.size());
            line1.insert(line1.end(), line2.begin(), line2.end());
            borderlines_.push_back(std::move(line1));
        }
        else if (edge_count > 2) {
            follow_branching_lines(visited, node);
        }
        else {
            std::cout << "Error: node (" << c.x << ", " << c.y << ") had "
                      << "less than 2 edges\n";
        } 
    }
}

// 0 - > 1, 1 -> 0, 2 -> 3, 3 -> 2
int complementary_edge(int e) {
    return (e == 0) + (e > 1) * (2 + (e == 2));
}

std::vector<Coordinate> SpacePartitioner::follow_line(
    std::unordered_set<GridVertex*>& visited,
    GridVertex* node,
    int edge
) { 
    std::vector<Coordinate> line = {node->coord};
    int edge_count;
    do {
        node = node->edges[edge];
        line.push_back(node->coord);
        visited.insert(node);
        int comp_edge = complementary_edge(edge);
        edge_count = 0;
        for (int i = 0; i < 4; i++) {
            if (node->edges[i] != nullptr && i != comp_edge) { 
                edge = i;
                ++edge_count; 
            }
        }
    } while (edge_count < 2);
    return line;    
}

void SpacePartitioner::follow_branching_lines(
    std::unordered_set<GridVertex*>& visited,
    GridVertex* node
) {
    int edge_count = 0;
    for (int i = 0; i < 4; i++) {
        if (node->edges[i] != nullptr) {
            if (visited.find(node->edges[i]) == visited.end()) {
                borderlines_.push_back(follow_line(visited, node, i));
                GridVertex* final_vertex = graph_[borderlines_.back().back()];
                follow_branching_lines(visited, final_vertex);
            }
            ++edge_count;
        }
    }
    if (edge_count < 2) {
        auto c = node->coord;
        std::cout << "WARNING: node (" << c.x << ", " << c.y << ") was not "
                  << "a junction, but it was treated as one\n";
    }
}

void SpacePartitioner::dfs(
    std::unordered_set<GridVertex*>& visited,
    GridVertex* node,
    int line_index,
    int exhausted_edge 
) {
    visited.insert(node);
    borderlines_[line_index].push_back(node->coord);
    std::array<int, 4> current_edge_ids = {};
    int current_edge_count = 0;
    for (int i = 0; i < 4; ++i) { 
        if (node->edges[i] != nullptr && i != exhausted_edge) {
            current_edge_ids[current_edge_count] = i;
            ++current_edge_count;
        }
    }

    if (current_edge_count == 0) {
        Coordinate& c = node->coord;
        std::cout << "Hit a dead end on (" << c.x << ", " << c.y << ")"
                  << std::endl;
    }
    else if (current_edge_count == 1) {
        auto& coord = node->edges[current_edge_ids[0]]->coord;
        //std::cout << "following coord: (" << coord.x << ", " << coord.y << ")" << std::endl;
        GridVertex* next_node = node->edges[current_edge_ids[0]];
        if (visited.find(next_node) == visited.end()) {
            int e = current_edge_ids[0];
            int complementary_edge = (e == 0) + (e > 1) * (2 + (e == 2));
            dfs(visited, next_node, line_index, complementary_edge);
        }
        else {
            borderlines_[line_index].push_back(next_node->coord);
        }
    }
    else {
        std::cout << "junction at: (" << node->coord.x << ", " << node->coord.y << ")\n";
        for (int i = 0; i < current_edge_count; i++) {

        }
        for (int i = 0; i < current_edge_count; i++) {
            auto& coord = node->edges[current_edge_ids[i]]->coord;
            //std::cout << "splitting coord: (" << coord.x << ", " << coord.y << ")" << std::endl;
            GridVertex* next_node = node->edges[current_edge_ids[i]];
            if (visited.find(next_node) == visited.end()) {
                line_index = borderlines_.size();
                borderlines_.push_back({node->coord});
                int e = current_edge_ids[i];
                int complementary_edge = (e == 0) + (e > 1) * (2 + (e == 2));
                dfs(visited, next_node, line_index, complementary_edge);
            }
        }
    }
}

void SpacePartitioner::label_areas_between_boundaries() {

    /*
     * function will expect the boundary_map is filled
     * edge locations in border_map:
     *    above: (2x, y)
     *    left:  (2x + 1, y)
     * 
     * looking above and the left, the boundary conditions
     * don't need to be considered. First row treated special
     * final column needs to check for union with first column
     * final row needs to check for union with first row
     */

    std::vector<int> label_parents(1, 0);
    region_map_[0][0] = 0;
    for (int j = 1; j < ysize_; j++) {
        if (border_map_[1][j] == -1) {
            region_map_[0][j] = region_map_[0][j - 1];
        }
        else {
            label_parents.push_back(label_parents.size());
            region_map_[0][j] = label_parents.back();
        }
    }
    if (border_map_[1][0] == -1 && region_map_[0][0] != region_map_[0].back()) {
        Impl::join_labels(
            label_parents, region_map_[0][0], region_map_[0].back()
        );
    }
    for (int i = 1; i < xsize_; i++) {
        if (border_map_[2 * i][0] == -1) {
            region_map_[i][0] = region_map_[i - 1][0];
        }
        else {
            label_parents.push_back(label_parents.size());
            region_map_[i][0] = label_parents.back();
        }
        for (int j = 1; j < ysize_; j++) {
            if (border_map_[2 * i][j] == -1) {
                region_map_[i][j] = region_map_[i - 1][j];
                if (border_map_[2 * i + 1][j] == -1 
                    && region_map_[i][j] != region_map_[i][j - 1]
                ) {
                    Impl::join_labels(
                        label_parents, region_map_[i][j], region_map_[i][j - 1]
                    );
                }
            }
            else if (border_map_[2 * i + 1][j] == -1) {
                region_map_[i][j] = region_map_[i][j - 1];
            }
            else {
                label_parents.push_back(label_parents.size());
                region_map_[i][j] = label_parents.back();
            }
        }
        if (border_map_[2 * i + 1][0] == -1 
            && region_map_[i][0] != region_map_[i].back()
        ) {
            Impl::join_labels(
                label_parents, region_map_[i][0], region_map_[i].back()
            );
        }
    }
    for (int j = 0; j < ysize_; j++) {
        if (border_map_[0][j] == -1 
            && region_map_[0][j] != region_map_.back()[j]
        ) {
            Impl::join_labels(
                label_parents, region_map_[0][j], region_map_.back()[j]
            );
        }
    }

    std::vector<int> label_map = Impl::relabel(label_parents);
    for (int i = 0; i < xsize_; i++) {
        for (int j = 0; j < ysize_; j++) {
            region_map_[i][j] = label_map[region_map_[i][j]];
        }
    }
}

std::vector<Coordinate> SpacePartitioner::adjacent_coords(Coordinate c) {
    int above = c.y - 1 < 0 ? ysize_ - 1 : c.y - 1;
    int below = c.y + 1 < ysize_ ? c.y + 1 : 0;
    int left = c.x - 1 < 0 ? xsize_ - 1 : c.x - 1;
    int right = c.x + 1 < xsize_ ? c.x + 1 : 0;

    return {
        {left, above}, {c.x, above}, {right, above}, {left, c.y}, 
        {right, c.y}, {left, below}, {c.x, below}, {left, below}
    };    
}

std::vector<OrderedCoordinate> SpacePartitioner::
unassigned_adjacent_grid_cells(Coordinate coord) {
    
    int affiliated_plate = region_map_[coord.x][coord.y];
    int size_x = region_map_.size();
    int size_y = region_map_[0].size();

    std::vector<OrderedCoordinate> unassigned_adjacent;
    auto adjacent_cells = adjacent_coords(coord); 
    for (Coordinate& c : adjacent_cells) {
        if (region_map_[c.x][c.y] == -1) {
            auto adj_adj_cells = adjacent_coords(c);
            int adj_count = 0;
            for (Coordinate& aac : adj_adj_cells) {
                adj_count += (
                    region_map_[aac.x][aac.y] == affiliated_plate
                );
            }
            unassigned_adjacent.emplace_back(c, adj_count);
        }
    }
    return unassigned_adjacent;
}

void SpacePartitioner::voronoi_partition(
    int xsize, int ysize, int num_seeds
) {
    random_x_ = std::uniform_int_distribution(0, xsize);
    random_y_ = std::uniform_int_distribution(0, ysize);
    xsize_ = xsize;
    ysize_ = ysize;
    initialize_containers();
    naive_voronoi(num_seeds);
    convert_region_map_to_graph();
    convert_graph_to_borderlines();
    convert_borderlines_to_border_map();
}

void SpacePartitioner::naive_voronoi(int num_seeds) {
    // naive is compute distance from each pixel to every seed, complexity: (elems * num seed)
    std::vector<Coordinate> seed_locations;
    seed_locations.reserve(num_seeds);
    for (int i = 0; i < num_seeds; i++) {
        int x = random_x_(rng_);
        int y = random_y_(rng_);
        region_map_[x][y] = i;
        seed_locations.emplace_back(x, y);
    }

    for (int i = 0; i < xsize_; i++) {
        for (int j = 0; j < ysize_; j++) {
            if (region_map_[i][j] == -1) {
                Coordinate c = closest_seed(seed_locations, {i, j});
                region_map_[i][j] = region_map_[c.x][c.y];
            }
        }
    }
}

Coordinate SpacePartitioner::closest_seed(
    std::vector<Coordinate>& seed_locations, Coordinate c
) {
    int closest_id = 0;
    double shortest_dist = toroidal_distance(seed_locations[0], c);
    for (int i = 1; i < seed_locations.size(); i++) {
        double dist = toroidal_distance(seed_locations[i], c);
        if (dist < shortest_dist) {
            closest_id = i;
            shortest_dist = dist;
        }
    }
    return seed_locations[closest_id];
}

double SpacePartitioner::euclidean_distance(
    Coordinate& c1, Coordinate& c2
) {
    double dx = c2.x - c1.x;
    double dy = c2.y - c1.y;
    return std::sqrt(dx * dx + dy * dy);
}

double SpacePartitioner::toroidal_distance(
    Coordinate& c1, Coordinate& c2
) {
    double dx = std::abs(c2.x - c1.x);
    if (dx > xsize_ / 2) {
        dx = xsize_ - dx;
    }
    double dy = std::abs(c2.y - c1.y);
    if (dy > xsize_ / 2) {
        dy = ysize_ - dy;
    }
    return std::sqrt(dx * dx + dy * dy);
}

void SpacePartitioner::convert_region_map_to_graph() {
    // require iterating over the region map and populating the nodes and 
    // edges as boundary lines are encountered. Consider upper and left
    // boundary for each cell (uniquely looking at each edge per grid cell)

    // assumes the graph is empty

    std::cout << "Converting the region map to a graph" << std::endl;

    for (int i = 0; i < xsize_; i++) {
        for (int j = 0; j < ysize_; j++) {
            if (region_map_[i][j] != region_map_[i][above(j)]) {
                ensure_nodes_and_edge({i, j}, {right(i), j}, 0, 1);
            }
            if (region_map_[i][j] != region_map_[left(i)][j]) {
                ensure_nodes_and_edge({i, j}, {i, below(j)}, 3, 2);
            }
        }
    }

    std::cout << "Size of graph: " << graph_.size() << std::endl;
}

// ensure the nodes and the edge between them exists
void SpacePartitioner::ensure_nodes_and_edge(
    Coordinate a, Coordinate b, int edge_a, int edge_b
) {
    GridVertex* node_a;
    GridVertex* node_b;
    if (auto search = graph_.find(a); search != graph_.end()) {
        node_a = search->second;
    }
    else {
        node_a = new GridVertex(a, graph_.size());
        graph_.emplace(a, node_a);
    }
    if (auto search = graph_.find(b); search != graph_.end()) {
        node_b = search->second;
    }
    else {
        node_b = new GridVertex(b, graph_.size());
        graph_.emplace(b, node_b);
    }
    node_a->edges[edge_a] = node_b;
    node_b->edges[edge_b] = node_a;
}

namespace Impl {
int find_parent(std::vector<int>& parent_labels, int label) {
    if (parent_labels[label] == label) {
        return label;
    }
    else {
        int parent = find_parent(parent_labels, parent_labels[label]);
        parent_labels[label] = parent;
        return parent;
    }
}

void join_labels(
    std::vector<int>& parent_labels, int label_a, int label_b
) {
    //std::cout << "Joining labels " << label_a << "  &  " << label_b << std::endl;
    int parent_a = find_parent(parent_labels, label_a);
    int parent_b = find_parent(parent_labels, label_b);
    parent_labels[parent_a] = parent_b;
}

std::vector<int> relabel(std::vector<int>& parent_labels) {
    std::unordered_map<int, int> parent_map;
    std::vector<int> relabelled(parent_labels.size());
    for (int i = 0;  i < parent_labels.size(); i++) {   
        int parent = find_parent(parent_labels, i);
        auto search = parent_map.find(parent);
        if (search == parent_map.end()) {
            int new_label = parent_map.size() + 1;
            relabelled[i] = new_label;
            parent_map[parent] = new_label;
        }
        else {
            relabelled[i] = search->second;
        }
    }
    return relabelled;
}
} // namespace Impl
} // namespace dev