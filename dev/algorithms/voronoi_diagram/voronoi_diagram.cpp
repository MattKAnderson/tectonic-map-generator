#include <voronoi_diagram.hpp>

/*
double compute_edge_y_val(BeachlineNode* node) {
    RealCoordinate& upper_focus = find_left_child(node)->coord;
    RealCoordinate& lower_focus = find_right_child(node)->coord;

}
*/
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


VoronoiDiagram::VoronoiDiagram(int seed) {
    rng = std::mt19937_64(seed);
}


void VoronoiDiagram::generate(int xsize, int ysize, int nseeds) {
    nregions = nseeds;
    seeds = generate_seeds(nseeds, xsize, ysize);
    //grid_algorithm(seeds, xsize, ysize);
    fortunes_algorithm(seeds, xsize, ysize);
}


void VoronoiDiagram::generate_from_seeds(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    //grid_algorithm(seeds, xsize, ysize);
    this->seeds = seeds;
    fortunes_algorithm(seeds, xsize, ysize);
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


std::vector<RealCoordinate> VoronoiDiagram::get_seeds() {
    return seeds;
}


std::vector<Node*> VoronoiDiagram::consume_vertices() {
    return std::move(vertices);
}


void BoundaryRay::add_vertex(const RealCoordinate v) {
    if (a == nullptr) {
        a = new RealCoordinate(v);
    }
    else if (b == nullptr) {
        b = new RealCoordinate(v);
    }
    else {
        std::cout << "Erorr: Tried to add a third point to line!!" << std::endl;
    }
}


FortunesAlgoEvent intersection_event(
    BeachLineItem* r1, BeachLineItem* r2, BeachLineItem* r3
) {
    RealCoordinate intersect = triangle_circumcenter(
        r1->coord(), r2->coord(), r3->coord()
    );
    double dist = euclidean_distance(r1->coord(), intersect);
    RealCoordinate known_at = {intersect.x + dist, intersect.y};

    auto c1 = r1->coord();
    auto c2 = r2->coord();
    auto c3 = r3->coord();
    //std::cout << "    Creating intersection from: (" << c1.x << ", " << c1.y << "), (" << c2.x << ", " << c2.y << "), (" << c3.x << ", " << c3.y << ")\n"; 
    //std::cout << "    Creating intersection event at: (" << known_at.x << ", " << known_at.y << ")" << std::endl;
    //std::cout << "    Intersection is: (" << intersect.x << ", " << intersect.y << ")" << std::endl;
    FortunesAlgoEvent event(known_at);
    event.intersect_point = new RealCoordinate(intersect);
    event.associated_region = r2; // it's the mid region that gets squeezed
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
        //std::cout << "[" << seed.x << ", " << seed.y << "],\n";
    }
    BeachLine beach_line(event_queue.top().coord); 
    event_queue.pop();

    while (!event_queue.empty()) {
        const FortunesAlgoEvent next_event = event_queue.top();
        auto coord = next_event.coord;
        double directrix = next_event.coord.x;
        if (next_event.intersect_point == nullptr) {
            //std::cout << "Handling SITE:\n"
            //          << "    location: (" << coord.x << ", " << coord.y << ")\n";
            BeachLineItem* new_region = beach_line.split_region(
                beach_line.find_intersected_region(next_event.coord),
                next_event.coord
            );
            event_queue.pop();
            
            BeachLineItem* l_region = beach_line.get_lower_edge(new_region);
            l_region = beach_line.get_lower_region(l_region);
            BeachLineItem* ll_region = beach_line.get_lower_edge(l_region);
            if (ll_region != nullptr) {
                ll_region = beach_line.get_lower_region(ll_region);
                FortunesAlgoEvent ie = intersection_event(
                    new_region, l_region, ll_region
                );
                //if (next_event.coord == RealCoordinate(882.0, 848.0)) {
                //    std::cout << "Intersection event is: (" << ie.intersect_point->x 
                //              << ", " << ie.intersect_point->y << ")\n";
                //}
                if (coord.y > ie.intersect_point->y) {
                    //std::cout << "    Published\n";
                    event_queue.push(ie);
                }
            }
            
            BeachLineItem* u_region = beach_line.get_upper_edge(new_region);
            u_region = beach_line.get_upper_region(u_region);
            BeachLineItem* uu_region = beach_line.get_upper_edge(u_region);
            if (uu_region != nullptr) {
                uu_region = beach_line.get_upper_region(uu_region);
                FortunesAlgoEvent ie = intersection_event(
                    new_region, u_region, uu_region
                );
                //if (next_event.coord == RealCoordinate(882.0, 848.0)) {
                //    std::cout << "Intersection event is: (" << ie.intersect_point->x 
                //              << ", " << ie.intersect_point->y << ")\n";
                //}
                if (coord.y < ie.intersect_point->y) {
                    //std::cout << "    Published\n";
                    event_queue.push(ie);
                }
            }

        }
        else if (next_event.associated_region->parent != nullptr) {
            auto inter = *next_event.intersect_point;
            //std::cout << "Handling INTERSECTION:\n"
            //          << "    location: (" << inter.x << ", " << inter.y << ")\n";
            
            BeachLineItem* l_region = beach_line.get_lower_edge(next_event.associated_region);
            l_region = beach_line.get_lower_region(l_region);
            BeachLineItem* u_region = beach_line.get_upper_edge(next_event.associated_region);
            u_region = beach_line.get_upper_region(u_region);
            beach_line.close_region(
                next_event.associated_region, *next_event.intersect_point
            );
            event_queue.pop();
            BeachLineItem* ll_region = beach_line.get_lower_edge(l_region);
            BeachLineItem* uu_region = beach_line.get_upper_edge(u_region);
            
            if (ll_region != nullptr) {
                ll_region = beach_line.get_lower_region(ll_region);
                if (ll_region->coord() != u_region->coord()) {
                    FortunesAlgoEvent ie = intersection_event(
                        u_region, l_region, ll_region
                    );
                    RealCoordinate& fu = u_region->coord();
                    RealCoordinate& fl = l_region->coord();
                    double midpoint_x = 0.5 * (fu.x + fl.x);
                    double arc_dir = (fu.y - fl.y);//* (midpoint_x - inter.x);
                    double int_dir = (ie.intersect_point->x - inter.x);
                    bool towards_intersect = (arc_dir * int_dir) >= 0;
                    RealCoordinate& fll = ll_region->coord();
                    //double midpoint_x2 = 0.5 * (fuu.x + fu.x);
                    double arc_dir_2 = (fl.y - fll.y);
                    RealCoordinate l_edge = parabola_intercept(
                        directrix, fl, fll
                    );
                    double int_dir_2 = (ie.intersect_point->x - l_edge.x);
                    bool arc_2_towards_intersect = (arc_dir_2 * int_dir_2) >= 0;
                    bool pushed = false;
                    //if (next_event.intersect_point->x > 848.29 && next_event.intersect_point->x < 848.3) {
                    //    std::cout << "Event in question: " 
                    //              << "    towards_intersect:   " << towards_intersect << "\n"
                    //              << "    towards_intersect_2: " << arc_2_towards_intersect << "\n"
                    //              << "    arc_dir_2: " << arc_dir_2 << "\n"
                    //              << "    int_dir_2: " << int_dir_2 << "\n";
                    //}
                    if (towards_intersect && arc_2_towards_intersect && ie.coord.x > directrix) {
                        if (*ie.intersect_point != *next_event.intersect_point) {
                            event_queue.push(ie);
                            pushed = true;
                        }
                    }
                
                //double dist = euclidean_distance(u_region->coord(), *ie.intersect_point); 
                //RealCoordinate u_edge = parabola_intercept(
                //    directrix, u_region->coord(), l_region->coord()
                //);
                //RealCoordinate l_edge = parabola_intercept(
                //    directrix, l_region->coord(), ll_region->coord()
                //);
                //BeachLineItem* l_edge_ptr = beach_line.get_lower_edge(l_region);
                //bool pushed = false;
                //if (*next_event.intersect_point != *ie.intersect_point) {
                //if (u_region->coord() != ll_region->coord()) {
                //if (
                //    (ie.intersect_point->x > l_edge.x && l_edge.x - l_edge_ptr->coord().x > 0)
                //) {
                //    if (ie.coord.x > directrix) {
                //        event_queue.push(ie);
                //        pushed = true;
                //   }
                //}
                //}
                //}
                //if (ie.intersect_point->x + dist > directrix) {
                //    event_queue.push(ie);
                //}
                //if (*next_event.intersect_point != *ie.intersect_point) {
                //    if (beach_line.is_behind(directrix, *ie.intersect_point)) {
                //        if (!pushed) {
                //        std::cout << "The beachline was behind of intersection, but it was not "
                //                  << "pushed. Int: (" << ie.intersect_point->x << ", " 
                //                  << ie.intersect_point->y << ")\n";
                //        }
                //    }
                //    else if (pushed) {
                //        std::cout << "The beachline was ahead of the prospective intersection "
                //                  << "at: (" << ie.intersect_point->x << ", " 
                //                  << ie.intersect_point->y << ")\n"; 
                        
                //    }
                }
            }
            if (uu_region != nullptr) {
                uu_region = beach_line.get_upper_region(uu_region);
                if (uu_region->coord() != l_region->coord()) {
                    FortunesAlgoEvent ie = intersection_event(
                        l_region, u_region, uu_region
                    );
                    RealCoordinate& fu = u_region->coord();
                    RealCoordinate& fl = l_region->coord();
                    double midpoint_x = 0.5 * (fu.x + fl.x);
                    double arc_dir = (fu.y - fl.y);// * (midpoint_x - inter.x);
                    double int_dir = (ie.intersect_point->x - inter.x);
                    bool towards_intersect = (arc_dir * int_dir) >= 0;

                    RealCoordinate& fuu = uu_region->coord();
                    //double midpoint_x2 = 0.5 * (fuu.x + fu.x);
                    double arc_dir_2 = (fuu.y - fu.y);
                    RealCoordinate u_edge = parabola_intercept(
                        directrix, fuu, fu
                    );
                    double int_dir_2 = (ie.intersect_point->x - u_edge.x);
                    bool arc_2_towards_intersect = (arc_dir_2 * int_dir_2) >= 0;
                    bool pushed = false;
                    if (towards_intersect && arc_2_towards_intersect && ie.coord.x > directrix) {
                        if (*ie.intersect_point != *next_event.intersect_point) {
                            event_queue.push(ie);
                            pushed = true;
                        }
                    }
                
                //double dist = euclidean_distance(l_region->coord(), *ie.intersect_point);
                //RealCoordinate l_edge = parabola_intercept(
                //    directrix, u_region->coord(), l_region->coord()
                //);
                //double l_edge_x_midpoint = 0.5 * (u_region->coord().x + l_region->coord().y);
                //BeachLineItem* l_edge_ptr = beach_line.get_lower_edge(u_region);
                //RealCoordinate u_edge = parabola_intercept(
                //    directrix, uu_region->coord(), u_region->coord()
                //);
                //BeachLineItem* u_edge_ptr = beach_line.get_upper_edge(u_region);
                //bool pushed = false;
                //if (*next_event.intersect_point != *ie.intersect_point) {
                //if (l_region->coord() != uu_region->coord()) {
                //if (
                //    (ie.intersect_point->x > u_edge.x && u_edge.x - u_edge_ptr->coord().x > 0)
                // ) {
                //    if (ie.coord.x > directrix) {
                //        event_queue.push(ie);
                //        pushed = true;
                //    }
                //}
                //}
                //}
                //if (!beach_line.is_behind(directrix, *ie.intersect_point)) {
                //    if (pushed) {
                //        std::cout << "The beachline was ahead of the prospective intersection "
                //                  << "at: (" << ie.intersect_point->x << ", " 
                //                  << ie.intersect_point->y << ")\n"; 
                //    }
                //}
                //else {
                //    if (!pushed) {
                //        std::cout << "The beachline was behind of intersection, but it was not "
                //                  << "pushed. Int: (" << ie.intersect_point->x << ", " 
                //                  << ie.intersect_point->y << ")\n";
                //    }
                //}
                //if (ie.intersect_point->x + dist > directrix) {
                //   event_queue.push(ie);
                //}
                //if (*next_event.intersect_point != *ie.intersect_point) {
                //if (beach_line.is_behind(directrix, *ie.intersect_point)) {
                //    std::cout << "    Published\n";
                //    event_queue.push(ie);
                //}
                //}
                }
            }
        }
        else {
            auto inter = *next_event.intersect_point;
            //std::cout << "Removing invalidated INTERSECTION\n"
            //          << "    location: (" << inter.x << ", " << inter.y << ")\n";
            event_queue.pop();
        }
    }
    beach_line.flush();
    vertices = beach_line.vertex_graph();
    //regions = beach_line.consume_region_graph();
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
    std::uniform_real_distribution<double> random_x(0, xsize - 1);
    std::uniform_real_distribution<double> random_y(0, ysize - 1);
    std::vector<RealCoordinate> seeds;
    seeds.reserve(nseeds);
    for (int i = 0; i < nseeds; i++) {
        //int x = random_x(rng);
        //int y = random_y(rng);
        double x = random_x(rng);
        double y = random_y(rng);
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
        RealCoordinate& focus_a = get_upper_region(node)->coord();
        RealCoordinate& focus_b = get_lower_region(node)->coord();
        double y = parabolae_y_intercept(loc.x, focus_a, focus_b);
        //if (loc == RealCoordinate(882, 848)) {
        //    std::cout << "    Computing edge location:\n"
        //              << "         focus_a: (" << focus_a.x << ", " << focus_a.y << ")\n"
        //              << "         focus_b: (" << focus_b.x << ", " << focus_b.y << ")\n"
        //              << "         y_int: " << y << std::endl;
        //}
        if (loc.y > y) {
            node = node->left;
        }
        else if (loc.y < y) {
            node = node->right;
        }
        else {
            // TODO: proper handling of the degenerate case
            std::cout << "Hit degenerate case, treating as loc.y < y\n";
            std::cout << "    focus_a: (" << focus_a.x << ", " << focus_a.y << ")\n"
                      << "    focus_b: (" << focus_b.x << ", " << focus_b.y << ")\n"
                      << "    loc    : (" << loc.x << ", " << loc.y << ")\n"
                      << "    y_inter:  " << y << std::endl;
            node = node->right;
        }
    }
    //if (loc == RealCoordinate(882.0, 848.0)) {
    //    std::cout << "    Intersected region with focus: (" << node->coord().x << ", " << node->coord().y << ")\n";
    //}
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


BeachLineItem* BeachLine::add_subtree(
    BeachLineItem* region, const RealCoordinate& focus
) {
    double x = parabola_x_from_y(focus.x, region->coord(), focus.y);
    RealCoordinate bound_start(x, focus.y);
    rays.push_back(new BoundaryRay());
    BeachLineItem* left_edge = new BeachLineItem(bound_start, rays.back());
    BeachLineItem* right_edge = new BeachLineItem(bound_start, rays.back());
    BeachLineItem* new_region = new BeachLineItem(focus);
    BeachLineItem* right_region = new BeachLineItem(region->coord());
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
    //region_graph.push_back(new Node(focus));
    //region_map.emplace(focus, region_graph.back());
    if (region == head) {
        head = add_subtree(region, focus);
        return head->right->left;
    }
    
    BeachLineItem* parent = region->parent;
    if (parent->left == region) {
        parent->left = add_subtree(region, focus);
        parent->left->parent = parent;
        return parent->left->right->left;
    }   
    else {
        parent->right = add_subtree(region, focus);
        parent->right->parent = parent;
        return parent->right->right->left;
    } 
}


BeachLineItem* BeachLine::close_region(
    BeachLineItem* region, const RealCoordinate& coord
) {
    //std::cout << "Starting to close region" << std::endl;
    BeachLineItem* parent = region->parent;
    BeachLineItem* other_child = nullptr;
    BeachLineItem* comp_edge = nullptr;
    if (parent->left == region) {
        comp_edge = get_upper_edge(region);
        other_child = parent->right;
    }
    else {
        comp_edge = get_lower_edge(region);
        other_child = parent->left;
    }
    if (parent->parent->left == parent) {
        parent->parent->left = other_child;
        other_child->parent = parent->parent;
    }
    else {
        parent->parent->right = other_child;
        other_child->parent = parent->parent;
    }
    
    parent->publish_vertex(coord);
    comp_edge->publish_vertex(coord);
    rays.push_back(new BoundaryRay());
    rays.back()->add_vertex(coord);
    comp_edge->associated_ray = rays.back();
    //link_regions(
    //    region->coord(),
    //    get_lower_region(grand_parent)->coord(),
    //    get_upper_region(grand_parent)->coord()
    //);
    delete parent;
    region->parent = nullptr; 
    //finished_regions.push_back(region);
    return comp_edge;
}


bool BeachLine::is_behind(double directrix, const RealCoordinate& loc) {
    RealCoordinate c = {directrix, loc.y};
    BeachLineItem* region = find_intersected_region(c);
    double beach_line_x = parabola_x_from_y(directrix, region->coord(), loc.y);
    double eps = 1e-10;
    //std::cout << "    Computed beachline x: " << beach_line_x << std::endl;
    return loc.x > beach_line_x + eps;
}


Node* ensure_node_in_graph(
    const RealCoordinate& c, 
    std::vector<Node*>& graph,
    std::unordered_map<RealCoordinate, Node*>& map    
) {
    if (auto it = map.find(c); it == map.end()) {
        graph.push_back(new Node(c));
        map.emplace(c, graph.back());
        return graph.back();
    }
    else {
        return it->second;
    }
    
}


/*
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
*/

/*void BeachLine::link_regions(
    const RealCoordinate& a, const RealCoordinate& b,
    const RealCoordinate& c
) {
    Node *ra = ensure_node_in_graph(a, region_graph, region_map);
    Node *rb = ensure_node_in_graph(b, region_graph, region_map);
    Node *rc = ensure_node_in_graph(c, region_graph, region_map); 
    ensure_region_link(ra, rb);
    ensure_region_link(ra, rc);
    ensure_region_link(rb, rc);
}*/


void BeachLine::flush() {
// TODO: This is the method to remove the remaining 
//       items from the beachline and create the 
//       associated region and vertex items
}


std::vector<Node*> BeachLine::vertex_graph() {
    //std::cout << "building vertex graph" << std::endl;
    std::vector<Node*> graph;
    std::unordered_map<RealCoordinate, Node*> vertex_map;
    for (BoundaryRay* ray : rays) {
        if (ray->b == nullptr) {
            continue;
        }
        Node* vertex_a = nullptr;
        Node* vertex_b = nullptr;
        const RealCoordinate& va = *ray->a;
        const RealCoordinate& vb = *ray->b; 
        //std::cout << ray->a->x << std::endl;
        if (auto s = vertex_map.find(va); s == vertex_map.end()) {
            vertex_a = new Node(va);
            graph.push_back(vertex_a);
            vertex_map.emplace(va, graph.back());
            vertex_a = graph.back();
        }
        else {
            vertex_a = s->second;
        }
        if (auto s = vertex_map.find(vb); s == vertex_map.end()) {
            graph.push_back(new Node(vb));
            vertex_map.emplace(vb, graph.back());
            vertex_b = graph.back();
        }
        else {
            vertex_b = s->second;
        }
        vertex_a->edges.push_back(vertex_b);
        vertex_b->edges.push_back(vertex_a);
    }
    return graph;
}


//std::vector<Node*> BeachLine::consume_vertex_graph() {
//    return std::move(vertex_graph_);
//}


//std::vector<Node*> BeachLine::consume_region_graph() {
//    return std::move(region_graph);
//}