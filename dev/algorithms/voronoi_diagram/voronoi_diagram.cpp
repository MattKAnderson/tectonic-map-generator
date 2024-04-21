#include <voronoi_diagram.hpp>

/*
double compute_edge_y_val(BeachlineNode* node) {
    RealCoordinate& upper_focus = find_left_child(node)->coord;
    RealCoordinate& lower_focus = find_right_child(node)->coord;

}
*/


VoronoiDiagram::VoronoiDiagram(int seed) {
    rng = std::mt19937_64(seed);
}


void VoronoiDiagram::generate(int xsize, int ysize, int nseeds) {
    nregions = nseeds;
    std::vector<RealCoordinate> seeds = generate_seeds(nseeds, xsize, ysize);
    grid_algorithm(seeds, xsize, ysize);
}


void VoronoiDiagram::generate_from_seeds(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    grid_algorithm(seeds, xsize, ysize);
}


void VoronoiDiagram::lloyd_iteration() {
    std::vector<RealCoordinate> new_seeds;
    new_seeds = compute_voronoi_cell_centroids(diagram, nregions);
    grid_algorithm(new_seeds, diagram.size(), diagram[0].size());
}


void VoronoiDiagram::nesting_iteration(int nseeds) {
    std::vector<RealCoordinate> seeds = generate_seeds(
        nseeds, diagram.size(), diagram[0].size()
    );
    nesting_iteration_from_seeds(seeds);
}


std::vector<RealCoordinate> compute_voronoi_cell_centroids(
    std::vector<std::vector<int>>& diagram, int ncells
) {
    std::vector<int> pixel_counts(ncells, 0);
    std::vector<int> x_sums(ncells, 0);
    std::vector<int> y_sums(ncells, 0);
    for (int i = 0; i < diagram.size(); ++i) {
        for (int j = 0; j < diagram[i].size(); ++j) {
            int cell = diagram[i][j];
            x_sums[cell] += i;
            y_sums[cell] += j;
            ++pixel_counts[cell];
        }
    }
    std::vector<RealCoordinate> centroids;
    centroids.reserve(ncells);
    for (int i = 0; i < ncells; i++) {
        if (pixel_counts[i] == 0) {
            continue;
        }
        centroids.emplace_back(
            static_cast<double>(x_sums[i]) / pixel_counts[i],
            static_cast<double>(y_sums[i]) / pixel_counts[i]
        );
    }
    return centroids;
}


void VoronoiDiagram::nesting_iteration_from_seeds(
    std::vector<RealCoordinate>& seeds
) {
    std::vector<RealCoordinate> centroids;
    centroids = compute_voronoi_cell_centroids(diagram, nregions);
    std::vector<int> region_id_map(centroids.size(), -1);
    for (int i = 0; i < centroids.size(); i++) {
        region_id_map[i] = find_closest_id(centroids[i], seeds);
    }
    for (int i = 0; i < diagram.size(); i++) {
        for (int j = 0; j < diagram[0].size(); j++) {
            diagram[i][j] = region_id_map[diagram[i][j]];
        }
    }
    nregions = seeds.size();
}


std::vector<std::vector<int>> VoronoiDiagram::get_diagram() {
    return diagram;
}


FortunesAlgoEvent intersection_event(
    BeachLineItem* r1, BeachLineItem* r2, BeachLineItem* r3
) {
    RealCoordinate intersect = triangle_circumcenter(
        r1->coord, r2->coord, r3->coord
    );
    double dist = euclidean_distance(r1->coord, intersect);
    RealCoordinate known_at = {intersect.x, intersect.y + dist};
    FortunesAlgoEvent event(known_at);
    event.intersect_point = new RealCoordinate(intersect);
    return event;
}


void VoronoiDiagram::fortunes_algorithm(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    std::priority_queue<FortunesAlgoEvent, 
                        std::vector<FortunesAlgoEvent>, 
                        std::greater<FortunesAlgoEvent>> event_queue;

    for (RealCoordinate& seed : seeds) {
        event_queue.push(FortunesAlgoEvent(seed));
    }
    BeachLine beach_line(event_queue.top().coord); 
    event_queue.pop();

    while (!event_queue.empty()) {
        const FortunesAlgoEvent& next_event = event_queue.top();
        if (next_event.intersect_point == nullptr) {
            BeachLineItem* new_region = beach_line.split_region(
                beach_line.find_intersected_region(next_event.coord),
                next_event.coord
            );
            
            BeachLineItem* l_region = beach_line.get_lower_edge(new_region);
            l_region = beach_line.get_lower_region(l_region);
            BeachLineItem* ll_region = beach_line.get_lower_edge(ll_region);
            if (ll_region != nullptr) {
                ll_region = beach_line.get_lower_region(ll_region);
                event_queue.push(intersection_event(new_region, l_region, ll_region));
            }
            
            BeachLineItem* u_region = beach_line.get_upper_edge(new_region);
            u_region = beach_line.get_upper_region(u_region);
            BeachLineItem* uu_region = beach_line.get_upper_edge(u_region);
            if (uu_region != nullptr) {
                uu_region = beach_line.get_upper_region(uu_region);
                event_queue.push(intersection_event(new_region, u_region, uu_region));
            }

        }
        else if (next_event.associated_region->parent != nullptr) {
            beach_line.close_region(
                next_event.associated_region, *next_event.intersect_point
            );
        }
    }
    beach_line.flush();
    vertices = beach_line.get_vertex_graph();
    regions = beach_line.get_region_graph();
}



int find_closest_id(RealCoordinate point, std::vector<RealCoordinate>& points) {
    int closest_id = 0;
    double shortest_dist = euclidean_distance(points[0], point);
    for (int i = 1; i < points.size(); i++) {
        double dist = euclidean_distance(points[i], point);
        if (dist < shortest_dist) {
            closest_id = i;
            shortest_dist = dist;
        }
    }
    return closest_id;
}


void VoronoiDiagram::grid_algorithm(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    for (int i = 0; i < seeds.size(); ++i) {
        diagram[seeds[i].x][seeds[i].y] = i;
    }
    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            if (diagram[i][j] == -1) {
                diagram[i][j] = find_closest_id(
                    {static_cast<double>(i), static_cast<double>(j)}, seeds
                );
            }
        }
    }
    
}


std::vector<RealCoordinate> VoronoiDiagram::generate_seeds(
    int nseeds, int xsize, int ysize
) {
    std::uniform_real_distribution random_x(0, xsize - 1);
    std::uniform_real_distribution random_y(0, ysize - 1);
    std::vector<RealCoordinate> seeds;
    seeds.reserve(nseeds);
    for (int i = 0; i < nseeds; i++) {
        int x = random_x(rng);
        int y = random_y(rng);
        seeds.emplace_back(x, y);
    }
    return seeds;
}


BeachLineItem* BeachLine::find_intersected_region(const RealCoordinate& loc) {
    if (head == nullptr || head->left == nullptr) {
        return head;
    }

    BeachLineItem* node = head;
    while (node->left != nullptr) {
        RealCoordinate& focus_a = get_upper_region(node)->coord;
        RealCoordinate& focus_b = get_lower_region(node)->coord;
        double y = parabolae_y_intercept(loc.x, focus_a, focus_b);
        if (loc.y > y) {
            node = node->left;
        }
        else if (loc.y < y) {
            node = node->right;
        }
        else {
            // TODO: proper handling of the degenerate case
            std::cout << "Hit degenerate case, treating as loc.y < y\n";
            node = node->right;
        }
    }
    return node;
}


BeachLineItem* BeachLine::get_upper_edge(BeachLineItem* node) {
    BeachLineItem* parent = node->parent;
    while (parent != nullptr && parent->right != node) {
        node = parent;
        parent = parent->parent;
    }
    return parent;
}


BeachLineItem* BeachLine::get_lower_edge(BeachLineItem* node) {
    BeachLineItem* parent = node->parent;
    while (parent != nullptr && parent->left != node) {
        node = parent;
        parent = parent->parent;
    }
    return parent;
}


BeachLineItem* BeachLine::get_upper_region(BeachLineItem* node) {
    node = node->left;
    while (node->right != nullptr) { node = node->right; }
    return node;
}


BeachLineItem* BeachLine::get_lower_region(BeachLineItem* node) {
    node = node->right;
    while (node->left != nullptr) { node = node->left; }
    return node;
}


BeachLineItem* add_subtree(
    BeachLineItem* region, const RealCoordinate& focus
) {
    double x = parabola_x_from_y(focus.x, region->coord, focus.y);
    RealCoordinate bound_start(x, focus.y);
    BeachLineItem* left_edge = new BeachLineItem(bound_start);
    BeachLineItem* right_edge = new BeachLineItem(bound_start);
    BeachLineItem* new_region = new BeachLineItem(focus);
    BeachLineItem* right_region = new BeachLineItem(region->coord);
    left_edge->left = region;
    region->parent = left_edge;
    left_edge->right = right_edge;
    right_edge->parent = left_edge;
    right_edge->left = new_region;
    new_region->parent = right_edge;
    right_edge->right = right_region;
    right_region->parent = right_edge;
    return left_edge;
}


BeachLineItem* BeachLine::split_region(
    BeachLineItem* region, const RealCoordinate& focus
) {
    region_graph.emplace_back(focus);
    region_map.emplace(focus, &region_graph.back());
    if (region == head) {
        head = add_subtree(region, focus);
        return head;
    }
    
    BeachLineItem* parent = head->parent;
    if (parent->left == region) {
        parent->left = add_subtree(region, focus);
        return parent->left;
    }   
    else {
        parent->right = add_subtree(region, focus);
        return parent->right;
    } 
}


BeachLineItem* BeachLine::close_region(
    BeachLineItem* region, const RealCoordinate& coord
) {
    BeachLineItem* parent = region->parent;
    BeachLineItem* grand_parent = parent->parent;
    BeachLineItem* other_child;
    other_child = (region == parent->left) ? parent->right : parent->left;
    other_child->parent = grand_parent;
    if (grand_parent->left == parent) {
        grand_parent->left = other_child;
    }
    else {
        grand_parent->right = other_child;
    }
    record_edge(grand_parent->coord, coord);
    record_edge(parent->coord, coord);
    link_regions(
        region->coord,
        get_lower_region(grand_parent)->coord,
        get_upper_region(grand_parent)->coord
    );
    grand_parent->coord = coord;
    delete parent;
    region->parent = nullptr; 
    finished_regions.push_back(region);
    return grand_parent;
}


Node* ensure_node_in_graph(
    const RealCoordinate& c, 
    std::vector<Node>& graph,
    std::unordered_map<RealCoordinate, Node*>& map    
) {
    if (auto it = map.find(c); it == map.end()) {
        graph.emplace_back(c);
        map.emplace(c, &graph.back());
        return &graph.back();
    }
    else {
        return it->second;
    }
    
}


void BeachLine::record_edge(
    const RealCoordinate& a, const RealCoordinate& b
) {
    Node *va = ensure_node_in_graph(a, vertex_graph, vertex_map);
    Node *vb = ensure_node_in_graph(b, vertex_graph, vertex_map);
    va->edges.push_back(vb);
    vb->edges.push_back(va);
}


void ensure_region_link(Node* a, Node* b) {
    if (std::find(a->edges.begin(), a->edges.end(), b) == a->edges.end()) {
        a->edges.push_back(b);
        b->edges.push_back(a);
    }
}


void BeachLine::link_regions(
    const RealCoordinate& a, const RealCoordinate& b,
    const RealCoordinate& c
) {
    Node *ra = ensure_node_in_graph(a, region_graph, region_map);
    Node *rb = ensure_node_in_graph(b, region_graph, region_map);
    Node *rc = ensure_node_in_graph(c, region_graph, region_map); 
    ensure_region_link(ra, rb);
    ensure_region_link(ra, rc);
    ensure_region_link(rb, rc);
}


void BeachLine::flush() {
// TODO: This is the method to remove the remaining 
//       items from the beachline and create the 
//       associated region and vertex items
}


std::vector<Node> BeachLine::get_vertex_graph() {
    return std::move(vertex_graph);
}


std::vector<Node> BeachLine::get_region_graph() {
    return std::move(region_graph);
}