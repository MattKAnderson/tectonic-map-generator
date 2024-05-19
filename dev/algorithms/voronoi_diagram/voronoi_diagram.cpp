#include <voronoi_diagram.hpp>

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
    
    //auto t1 = std::chrono::high_resolution_clock::now();
    seeds = generate_seeds(nseeds, xsize, ysize);
    //auto t2 = std::chrono::high_resolution_clock::now();
    //double et = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    //std::cout << "Time spent generating seeds: " << et / 1e6 << " ms\n";
    //grid_algorithm(seeds, xsize, ysize);
    
    Impl::FortunesAlgorithm generator;
    generator.compute(seeds, 0.0, xsize);
    //t1 = std::chrono::high_resolution_clock::now();
    vertices = generator.consume_vertex_graph();
    //t2 = std::chrono::high_resolution_clock::now();
    //et = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    //std::cout << "Time spent creating vertex graph: " << et / 1e6 << " ms\n";
    //t1 = std::chrono::high_resolution_clock::now();
    regions = generator.consume_region_graph();
    //t2 = std::chrono::high_resolution_clock::now();
    //et = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
    //std::cout << "Time spent creating region graph: " << et / 1e6 << " ms\n";
}

void VoronoiDiagram::generate_from_seeds(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    //grid_algorithm(seeds, xsize, ysize); 
    this->seeds = seeds;
    Impl::FortunesAlgorithm generator;
    generator.compute(seeds, 0.0, xsize);
    //std::cout << "about to produce vertex graph" << std::endl;
    vertices = generator.consume_vertex_graph();
    regions = generator.consume_region_graph();
    
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

std::vector<VertexNode*> VoronoiDiagram::consume_vertices() {
    return std::move(vertices);
}

std::vector<RegionNode*> VoronoiDiagram::consume_region_graph() {
    return std::move(regions);
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

namespace Impl {

Arc* BeachLine::find_intersected_arc(const RealCoordinate& c) {
    Arc* node = head;
    double y;
    //std::cout << "c: (" << c.x << ", " << c.y << ")" << std::endl;
    // TODO: handle case when c.y == y
    //       it is super rare but I've seen it at 100k seeds
    while (true) {
        if (node->left != nullptr) {
            if (node->focus.x < node->upper->focus.x && c.y > node->upper->focus.y) {
                node = node->left;
                continue;
            }
            if (node->focus.x <= node->upper->focus.x || c.y > node->focus.y) {
                y = parabolae_y_intercept(c.x, node->upper->focus, node->focus);
                if (c.y > y) {
                    node = node->left;
                    continue;
                } 
            }
        }
        if (node->right != nullptr) {
            if (node->focus.x < node->lower->focus.x && c.y < node->lower->focus.y) {
                node = node->right;
                continue;
            }
            if (node->focus.x <= node->lower->focus.x || c.y < node->focus.y) {
                y = parabolae_y_intercept(c.x, node->focus, node->lower->focus);
                if (c.y < y) {
                    node = node->right;
                    continue;
                }
            }
        }
        return node;
    }
}   

Arc* BeachLine::get_lowest() {
    Arc* lowest = head;
    while (lowest->right != nullptr) { lowest = lowest->right; }
    return lowest;
}

void BeachLine::insert_arc_above(Arc* arc, Arc* new_arc) {
    if (arc->left == nullptr) {
        arc->left = new_arc;
        new_arc->parent = arc;
        if (arc->upper != nullptr) {
            arc->upper->lower = new_arc;
        }
    }
    else {
        arc->upper->right = new_arc;
        new_arc->parent = arc->upper;
        arc->upper->lower = new_arc;
    }
    new_arc->upper = arc->upper;
    new_arc->lower = arc;
    arc->upper = new_arc;
    //insert_balance(new_arc);
}

void BeachLine::insert_balance(Arc* arc) {
    arc->red = true;
    arc = arc->parent;
    Arc* parent = arc->parent;
    while (parent != nullptr) {
        if (parent->left == arc) { parent->left = insert_local_balance(arc); }
        else { parent->right = insert_local_balance(arc); }
        arc = parent;
        parent = parent->parent;
    }
    arc = insert_local_balance(arc);
    arc->red = false;
    head = arc;
}

Arc* BeachLine::insert_local_balance(Arc* arc) {
    if (is_red(arc->right) && !is_red(arc->right)) { arc = rotate_left(arc); }
    if (is_red(arc->left) && is_red(arc->left->left)) { arc = rotate_right(arc); }
    if (is_red(arc->left) && is_red(arc->right)) { flip_colors(arc); }
    return arc;
}

void BeachLine::delete_balance(Arc* arc) {
    // TODO
}

Arc* BeachLine::rotate_left(Arc* arc) {
    Arc* rchild = arc->right;
    rchild->parent = arc->parent;
    arc->parent = rchild;
    arc->right = rchild->left;
    if (arc->right) { arc->right->parent = arc; }
    arc->right->parent = arc;
    rchild->left = arc;
    rchild->red = arc->red;
    arc->red = true;
    return rchild;
}

Arc* BeachLine::rotate_right(Arc* arc) {
    Arc* lchild = arc->left;
    lchild->parent = arc->parent;
    arc->parent = lchild;
    arc->left = lchild->right;
    if (arc->left) { arc->left->parent = arc; }
    lchild->right = arc;
    lchild->red = arc->red;
    arc->red = true;
    return lchild;
}

void BeachLine::flip_colors(Arc* arc) {
    arc->red = true;
    arc->left->red = false;
    arc->right->red = false;
}

bool BeachLine::is_red(Arc* arc) {
    return arc != nullptr && arc->red;
}

/*
void BeachLine::insert_arc_below(Arc* arc, Arc* new_arc) {
    if (arc->right == nullptr) {
        arc->right = new_arc;
        new_arc->parent = arc;
        if (arc->lower != nullptr) {
            arc->lower->upper = new_arc;
        }
    }
    else {
        arc->lower->left = new_arc;
        new_arc->parent = arc->lower;
        arc->lower->upper = new_arc;
    }
    new_arc->lower = arc->lower;
    new_arc->upper = arc;
    arc->lower = new_arc;
    // TODO: balance op
}*/

void BeachLine::remove_arc(Arc* arc) {
    // Handle the doubly linked list portion of this
    // Data structure
    if (arc->upper){
        arc->upper->lower = arc->lower;
    }
    if (arc->lower) {
        arc->lower->upper = arc->upper;
    }

    if (arc->left == nullptr) {
        if (arc->right == nullptr) {
            // leaf node case
            if (arc == head) { head = nullptr; }
            else if (arc->parent->left == arc) { arc->parent->left = nullptr; }
            else { arc->parent->right = nullptr; }

        }
        else {
            // single child on right case
            arc->right->parent = arc->parent;
            if (arc == head) { head = arc->right; }
            else if (arc->parent->left == arc) { arc->parent->left = arc->right; }
            else {arc->parent->right = arc->right; }
        } 
    }
    else if (arc->right == nullptr) {
        // single child on left case
        arc->left->parent = arc->parent;
        if (arc == head) { head = arc->left; }
        else if (arc->parent->left == arc) { arc->parent->left = arc->left; }
        else {arc->parent->right = arc->left; }
    }
    else {
        // two children case
        Arc* upper = arc->upper;

        // connect upper's child to upper's parent
        if (upper->parent->left == upper) { upper->parent->left = upper->left; }
        else { upper->parent->right = upper->left; }
        if (upper->left) { upper->left->parent = upper->parent; }

        // put upper in place of arc
        upper->parent = arc->parent;
        upper->left = arc->left;
        upper->right = arc->right;
        upper->red = arc->red;
        if (arc->left) { arc->left->parent = upper; }
        if (arc->right) { arc->right->parent = upper; }
        if (arc == head) { head = upper; }
        else if (arc->parent->left == arc) { arc->parent->left = upper; }
        else { arc->parent->right = upper; }
    }
    available_arcs.push_back(arc);
    //closed_regions.push_back(arc);
    // TODO: delete_balance(y);
}

//void BeachLine::reserve(int n) {
//    closed_regions.reserve(n);
//}

Event new_intersection_event(RealCoordinate& intersect, Arc* closing_region) {
    double dist = euclidean_distance(closing_region->focus, intersect);
    Event event({intersect.x + dist, intersect.y}, intersect, closing_region);
    return event;
}

void FortunesAlgorithm::site_event(const RealCoordinate& focus) {
    regions[next_region_id] = {focus, nullptr};
    Arc* arc = beach_line.find_intersected_arc(focus);
    Arc* new_arc = beach_line.new_arc(focus, &regions[next_region_id]); //new Arc{focus, &regions[next_region_id]};
    Arc* split_arc = beach_line.new_arc(arc->focus, arc->region); //new Arc{arc->focus, arc->region};
    ++next_region_id;

    beach_line.insert_arc_above(arc, new_arc);
    beach_line.insert_arc_above(new_arc, split_arc);

    new_arc->upper_edge = new HalfEdge{new_arc->region};
    half_edges.push_back(new_arc->upper_edge);
    new_arc->lower_edge = new_arc->upper_edge;
    new_arc->region->an_edge = new_arc->upper_edge;
    split_arc->upper_edge = arc->upper_edge;
    arc->upper_edge = new HalfEdge{arc->region};
    half_edges.push_back(arc->upper_edge);
    split_arc->lower_edge = arc->upper_edge;
    arc->upper_edge->twin = new_arc->upper_edge;
    new_arc->lower_edge->twin = arc->upper_edge;
    arc->region->an_edge = arc->upper_edge;
    new_arc->region->an_edge = new_arc->upper_edge;

    Arc* upper_upper = split_arc->upper;
    Arc* lower_lower = arc->lower;
    if (upper_upper) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, split_arc->focus, upper_upper->focus
        );
        if (focus.y < intersect.y) {
            double dist = euclidean_distance(split_arc->focus, intersect); 
            int event_id = event_manager.create(
                {intersect.x + dist, intersect.y}, intersect, split_arc
            );
            //auto t1 = std::chrono::high_resolution_clock::now();
            if (split_arc->event_id != -1) {
                event_queue.remove(split_arc->event_id);
            }
            split_arc->event_id = event_id;
            event_queue.insert(event_id);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //queue_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
        }
    } 
    if (lower_lower) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, arc->focus, lower_lower->focus
        );
        if (focus.y > intersect.y) {
            double dist = euclidean_distance(arc->focus, intersect);
            int event_id = event_manager.create(
                {intersect.x + dist, intersect.y}, intersect, arc
            );
            //auto t1 = std::chrono::high_resolution_clock::now();
            if (arc->event_id != -1) {
                event_queue.remove(arc->event_id);
            }
            arc->event_id = event_id;
            event_queue.insert(event_id);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //queue_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
        }
    }
}

double ray_direction(
    const RealCoordinate& intersection,
    const RealCoordinate& focus_1,
    const RealCoordinate& focus_2,
    const RealCoordinate& focus_3
) {
    if (focus_1.x == focus_2.x) {
        //similar for f1.y == f2.y ???
        return (focus_3.x - focus_1.x) * (intersection.x - focus_1.x);
    }
    //else if (focus_1.y == focus_2.y) {
    //
    //}
    else if (focus_1.x > focus_2.x) {
        return (intersection.y - focus_1.y) * (focus_2.y - focus_1.y);
    }
    else {
        return (intersection.y - focus_2.y) * (focus_1.y - focus_2.y);
    }
}

void FortunesAlgorithm::intersection_event(const Event& event) {
    const RealCoordinate& intersect = event.intersect_point;
    Arc* arc = event.associated_arc;
    Arc* u_arc = arc->upper;
    Arc* l_arc = arc->lower;

    VertexNode* new_intersect = new VertexNode(intersect);
     
    vertices.push_back(new VertexNode(intersect));

    arc->upper_edge->origin = vertices.back();
    arc->upper_edge->prev = arc->lower_edge;
    arc->lower_edge->next = arc->upper_edge;
    HalfEdge* new_upper_half_edge = new HalfEdge{u_arc->region};
    HalfEdge* new_lower_half_edge = new HalfEdge{l_arc->region};
    half_edges.push_back(new_upper_half_edge);
    half_edges.push_back(new_lower_half_edge);
    new_upper_half_edge->twin = new_lower_half_edge;
    new_upper_half_edge->origin = vertices.back();
    new_upper_half_edge->prev = u_arc->lower_edge;
    u_arc->lower_edge->next = new_upper_half_edge;
    u_arc->lower_edge = new_upper_half_edge;
    new_lower_half_edge->twin = new_upper_half_edge;
    new_lower_half_edge->next = l_arc->upper_edge;
    l_arc->upper_edge->prev = new_lower_half_edge;
    l_arc->upper_edge->origin = vertices.back();
    l_arc->upper_edge = new_lower_half_edge;

    beach_line.remove_arc(arc);
    event_queue.remove(arc->event_id);

    const RealCoordinate& fu = u_arc->focus;
    const RealCoordinate& fl = l_arc->focus;
    Arc* uu_arc = u_arc->upper;
    Arc* ll_arc = l_arc->lower;
    if (uu_arc && fl != uu_arc->focus) {
        const RealCoordinate& fuu = uu_arc->focus;
        RealCoordinate new_intersect = triangle_circumcenter(fl, fu, fuu);
        double known_at_x = new_intersect.x + euclidean_distance(new_intersect, fu);
        if (new_intersect != intersect 
            && event.coord.x < known_at_x 
            && (fu.y - fl.y) * (new_intersect.x - intersect.x) >= 0.0
        ) {
            RealCoordinate u_edge = parabola_intercept(event.coord.x, fuu, fu);
            if ((fuu.y - fu.y) * (new_intersect.x - u_edge.x) >= 0.0) {
                int event_id = event_manager.create(
                    {known_at_x, new_intersect.y}, new_intersect, u_arc
                );
                if (u_arc->event_id != -1) {
                    event_queue.remove(u_arc->event_id);
                }
                u_arc->event_id = event_id;
            //auto t1 = std::chrono::high_resolution_clock::now();
                event_queue.insert(event_id);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //queue_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
            }
        }

    }
    if (ll_arc && u_arc->focus != ll_arc->focus) {
        const RealCoordinate& fll = ll_arc->focus;
        RealCoordinate new_intersect = triangle_circumcenter(fu, fl, fll);
        double known_at_x = new_intersect.x + euclidean_distance(new_intersect, fl);
        if (new_intersect != intersect
            && event.coord.x < known_at_x
            && (fu.y - fl.y) * (new_intersect.x  - intersect.x) >= 0.0
        ) {
            RealCoordinate l_edge = parabola_intercept(event.coord.x, fl, fll);
            if ((fl.y - fll.y) * (new_intersect.x - l_edge.x) >= 0.0) {
                int event_id = event_manager.create(
                    {known_at_x, new_intersect.y}, new_intersect, l_arc
                );
                if (l_arc->event_id != -1) {
                    event_queue.remove(l_arc->event_id);
                }
                l_arc->event_id = event_id;
            //auto t1 = std::chrono::high_resolution_clock::now();
                event_queue.insert(event_id);
            //event_queue.push(event_id);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //queue_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count());
            }
        }
    }
}

bool compare_x_lower(const RealCoordinate& a, const RealCoordinate& b) {
    return a.x < b.x;
}

bool compare_x_greater(const RealCoordinate& a, const RealCoordinate& b) {
    return a.x > b.x;
}

bool compare_y_lower(const RealCoordinate& a, const RealCoordinate& b) {
    return a.y < b.y;
}

bool compare_y_greater(const RealCoordinate& a, const RealCoordinate& b) {
    return a.y > b.y;
}


void FortunesAlgorithm::compute(std::vector<RealCoordinate>& seeds, double min, double max) {
    using namespace Impl; 
    this->min = min;
    this->max = max;
    num_seeds = seeds.size();
    if (regions != nullptr) {
        delete[] regions;
    }
    regions = new Region[num_seeds];
    half_edges = {};
    half_edges.reserve(num_seeds * 3 - 6);
    beach_line = BeachLine();
    event_manager = EventManager();
    event_queue = EventQueue();
    event_queue.set_event_manager(&event_manager);
    int event_id;
    for (const RealCoordinate& seed : seeds) { 
        event_id = event_manager.create(seed);
        event_queue.insert(event_id);
    }
    //event_queue.print_ordered_x();
    event_id = event_queue.consume_next();
    const RealCoordinate& s1 = event_manager.get(event_id).coord;
    regions[next_region_id] = {s1, nullptr};
    beach_line.set_head(beach_line.new_arc(s1, &regions[next_region_id]));
    ++next_region_id;
    event_manager.remove(event_id);
    auto t1 = std::chrono::high_resolution_clock::now();
    //std::cout << "About to start main loop" << std::endl;
    while (!event_queue.empty()) {
        
        //auto t1a = std::chrono::high_resolution_clock::now();
        event_id = event_queue.consume_next();
        //auto t2a = std::chrono::high_resolution_clock::now();
        const Event& event = event_manager.get(event_id);
        //queue_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2a - t1a).count());

        if (event.associated_arc == nullptr) { 
            //auto t1 = std::chrono::high_resolution_clock::now();
            //std::cout << "Site x: (" << event.coord.x << ", " << event.coord.y << ")" << std::endl;
            site_event(event.coord); 
            //auto t2 = std::chrono::high_resolution_clock::now();
            //site_event_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        }
        else {
            //auto t1 = std::chrono::high_resolution_clock::now();
            //std::cout << "Int:  (" << event.coord.x <<  ", " << event.coord.y << ")" << std::endl;
            intersection_event(event); 
            //auto t2 = std::chrono::high_resolution_clock::now();
            //int_event_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count());
        }
        event_manager.remove(event_id);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    double et = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    std::cout << "Time spent on main loop: " << et / 1e6 << " ms" << std::endl;
    //t1 = std::chrono::high_resolution_clock::now();
    bound_DCEL();
    //t2 = std::chrono::high_resolution_clock::now();
    //et = std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count();
    //std::cout << "Time spent bounding DCEL: " << et / 1e6 << " ms\n";


    //double find_time = std::accumulate(find_arc_times.begin(), find_arc_times.end(), 0.0);
    //double ins_time = std::accumulate(insert_arc_times.begin(), insert_arc_times.end(), 0.0);
    //double del_time = std::accumulate(delete_arc_times.begin(), delete_arc_times.end(), 0.0);
    //double new_time = std::accumulate(new_allocation_times.begin(), new_allocation_times.end(), 0.0);
    //double site_time = std::accumulate(site_event_times.begin(), site_event_times.end(), 0.0);
    //double int_time = std::accumulate(int_event_times.begin(), int_event_times.end(), 0.0);
    //double queue_time = std::accumulate(queue_times.begin(), queue_times.end(), 0.0);
    //double site_new_int_time = std::accumulate(site_new_int_times.begin(), site_new_int_times.end(), 0.0);
    //std::cout << "Number of triangle circumcenter ops: " << op_counter << std::endl;
    //std::cout << "Aggregate time spent finding: " << find_time / 1e6 << " ms\n"
    //          << "Per call time spent finding:  " << find_time / find_arc_times.size() << " ns\n\n"
    //          << "Aggregate time spent inserting: " << ins_time / 1e6 << " ms\n"
    //          << "Per call time spent inserting:  " << ins_time / insert_arc_times.size() << " ns\n\n"
    //         << "Aggregate time spent deleting:  " << del_time / 1e6 << " ms\n" 
    //          << "Per call time spent deleting:   " << del_time / delete_arc_times.size() << " ns\n\n"
    //          << "New allocations time: " << new_time / 1e6 << " ms\n"
    //std::cout          << "Aggregate time spent on site events:    " << site_time / 1e6 << " ms\n"
    //          << "Avg Per call time spent on site events: " << site_time / site_event_times.size() << " ns\n"
    //          << "Aggregate time spent on site new int events:    " << site_new_int_time / 1e6 << " ms\n"
    //          << "Avg Per call time spent on site new int events: " << site_new_int_time / site_new_int_times.size() << " ns\n\n"
    //          << "Aggregate time spent on int events:     " << int_time / 1e6 << " ms\n"
    //          << "Avg Per call time spent on int events:  " << int_time / int_event_times.size() << " ns\n"
    //          << "Aggregate time spent on queue:    " << queue_time / 1e6 << " ms\n"
    //          << "Avg Per call time spent on queue: " << queue_time / queue_times.size() << " ns\n";
;
}

RealCoordinate clip_infinite_edge(
    HalfEdge* edge, double xmax, double xmin, double ymax, double ymin
) {
    const auto& [x0, y0] = edge->twin->origin->coord;
    const auto& [rx1, ry1] = edge->region->seed;
    const auto& [rx2, ry2] = edge->twin->region->seed;
    double y_int, x_int;
    if (rx1 == rx2) {
        double rx3 = edge->next->twin->region->seed.x;
        x_int = (rx3 - rx1) * (x0 - rx1) > 0 ? xmax : xmin;
        y_int = 0.5 * (ry1 + ry2);
    }
    else if (ry1 == ry2) {
        // TODO figure this case out properly
        x_int = 0.5 * (rx1 + rx2);
        y_int = y0 > ry1 ? ymin : ymax;
    }
    else {
        RealCoordinate mid = {0.5 * (rx1 + rx2), 0.5 * (ry1 + ry2)};
        if (x0 - mid.x == 0.0) { std::cout << "It was zero" << std::endl; }
        double m = (y0 - mid.y) / (x0 - mid.x);
        double b = y0 - m * x0;
        x_int = ry2 > ry1 ? xmax : xmin;
        y_int = m * x_int + b;
        if (y_int < ymin) {
            y_int = ymin;
            x_int = (y_int - b) / m;
        }
        else if (y_int > ymax) {
            y_int = ymax;
            x_int = (y_int - b) / m;
        }
    }
    return {x_int, y_int};
}

bool is_before_on_bbox_exterior(
    const RealCoordinate& a, const RealCoordinate& b, double xmax, double xmin,
    double ymax, double ymin
) { 
    if (a.x == xmin) {
        if (b.x == xmin) { return a.y > b.y; }
        else { return a.y != ymax; }
    }
    if (a.y == ymin) {
        if (b.y == ymin) { return a.x < b.x; }
        else if (b.x == xmin) { return false; }
        else { return true; }
    }
    if (a.x == xmax) {
        if (b.x == xmax) { return a.y < b.y; }
        else if (b.y == ymax) { return true; }
        else { return false; }
    }
    return b.y == ymax && a.x > b.x;
}

bool is_corner_node(
    VertexNode* node, double xmax, double xmin, double ymax, double ymin
) {
    if (
        node->coord == RealCoordinate(xmin, ymin)
        || node->coord == RealCoordinate(xmin, ymax)
        || node->coord == RealCoordinate(xmax, ymin)
        || node->coord == RealCoordinate(xmax, ymax)
    ) { 
        return true; 
    }
    return false;
}

void FortunesAlgorithm::bound_DCEL() {
    double xmin = std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double xmax = std::numeric_limits<double>::min();
    double ymax = std::numeric_limits<double>::min();
    for (VertexNode* v : vertices) {
        xmin = std::min(xmin, v->coord.x);
        ymin = std::min(ymin, v->coord.y);
        xmax = std::max(xmax, v->coord.x);
        ymax = std::max(ymax, v->coord.y);
    }

    std::list<VertexNode*> exterior;
    std::unordered_map<VertexNode*, HalfEdge*> node_to_edge_in;
    Arc* lower = beach_line.get_lowest();
    Arc* upper = lower->upper;
    while (upper != nullptr) {
        HalfEdge* edge = upper->lower_edge->origin ? lower->upper_edge : upper->lower_edge;
        RealCoordinate v = clip_infinite_edge(edge, xmax, xmin, ymax, ymin);
        vertices.push_back(new VertexNode(v));
        exterior.push_back(vertices.back());
        edge->origin = vertices.back();
        node_to_edge_in.insert({vertices.back(), edge->twin});
        lower = upper;
        upper = upper->upper;
    }

    RealCoordinate corners[4] = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
    int corner_to_place = 0;
    for (auto it = exterior.begin(); it != exterior.end(); ++it) {
        for (int i = corner_to_place; i < 4; i++) {
            if (corners[i] == (*it)->coord) {
                ++corner_to_place;
                break;
            }
            else if (is_before_on_bbox_exterior(corners[i], (*it)->coord, xmax, xmin, ymax, ymin)) {
                vertices.push_back(new VertexNode(corners[i]));
                exterior.insert(it, vertices.back());
                HalfEdge* dummy_edge = new HalfEdge;
                dummy_edge->region = node_to_edge_in[(*it)]->twin->region;
                node_to_edge_in.insert({vertices.back(), dummy_edge});
                ++corner_to_place;
            }
            else {
                break;
            }
        }
    }

    HalfEdge* prev_edge = nullptr;
    HalfEdge* first_corner_edge = nullptr;
    HalfEdge* edge_in = nullptr;
    for (VertexNode* vertex : exterior) {
        edge_in = node_to_edge_in[vertex];
        if (edge_in->origin == nullptr) {
            HalfEdge* new_edge = new HalfEdge{
                edge_in->region, prev_edge, nullptr, nullptr, vertex
            };
            HalfEdge* twin_edge = new HalfEdge{
                nullptr, nullptr, nullptr, new_edge, nullptr
            };
            new_edge->twin = twin_edge;
            if (first_corner_edge == nullptr) { first_corner_edge = new_edge; }
            if (prev_edge) { 
                prev_edge->next = new_edge;
                prev_edge->twin->origin = vertex;
            }
            prev_edge = new_edge;
            half_edges.push_back(new_edge);
            half_edges.push_back(twin_edge);
        }
        else {
            HalfEdge* new_edge = new HalfEdge{
                edge_in->region, edge_in, nullptr, nullptr, vertex
            };
            HalfEdge* twin_edge = new HalfEdge{
                nullptr, nullptr, nullptr, new_edge, nullptr
            };
            new_edge->twin = twin_edge;
            if (prev_edge) { 
                prev_edge->next = edge_in->twin; 
                prev_edge->twin->origin = vertex;
                edge_in->twin->prev = prev_edge;    
            }
            prev_edge = new_edge;
            edge_in->next = prev_edge;
            half_edges.push_back(prev_edge);
            half_edges.push_back(twin_edge);
        }
    }
    edge_in = node_to_edge_in[*exterior.begin()];
    if (edge_in->origin == nullptr) {
        first_corner_edge->prev = prev_edge;
        prev_edge->next = first_corner_edge;
        prev_edge->twin->origin = *exterior.begin();
    }
    else {
        HalfEdge* next = edge_in->twin;
        prev_edge->next = next;
        prev_edge->twin->origin = *exterior.begin();
        next->prev = prev_edge;
    }
}

bool outside_bbox(const RealCoordinate& c, double min, double max) {
    return c.x < min || c.x > max || c.y < min || c.y > max;
}

bool inside_bbox(const RealCoordinate& c, double min, double max) {
    return c.x >= min && c.x <= max && c.y >= min && c.y <= max;
}

void FortunesAlgorithm::clip_to_bbox(double min, double max) {

    for (auto* edge : half_edges) {
        const RealCoordinate& v = edge->origin->coord;
        if (inside_bbox(v, min, max)) { continue; }

    }

}

std::vector<VertexNode*> FortunesAlgorithm::consume_vertex_graph() {
    //std::cout << "Consuming vertex graph" << std::endl;
    //std::cout << "Size of regions: " << regions.size() << std::endl;
    for (HalfEdge* half_edge : half_edges) {
        half_edge->origin->connected.push_back(half_edge->twin->origin);
    }
    return std::move(vertices);
}

std::vector<RegionNode*> FortunesAlgorithm::consume_region_graph() {
    return std::move(region_graph_from_regions());
}

std::vector<RegionNode*> FortunesAlgorithm::region_graph_from_regions() {
    std::vector<RegionNode*> nodes;
    std::unordered_map<Region*, RegionNode*> node_map;
    nodes.reserve(num_seeds);
    for (int i = 0; i < num_seeds; i++) {
        nodes.push_back(new RegionNode);
        node_map.insert({&regions[i], nodes.back()});
    }
    //std::cout << "Finished iterating over regions" << std::endl;
    for (int i = 0; i < num_seeds; i++) {
        
        RegionNode* this_region = node_map[&regions[i]];
        HalfEdge* edge_ptr = regions[i].an_edge;
        if (edge_ptr == nullptr) {
            std::cout << "edge_ptr was nullptr" << std::endl;
        }
        while (edge_ptr->next != regions[i].an_edge) {
            this_region->vertices.push_back(edge_ptr->origin->coord);
            if (edge_ptr->twin->region != nullptr) {
                this_region->adjacent.push_back(node_map[edge_ptr->twin->region]);
            }
            edge_ptr = edge_ptr->next;
        }
        this_region->vertices.push_back(edge_ptr->origin->coord);
    }
    return std::move(nodes);
}

const Event& EventManager::get(int id) {
    return events[id];
}

int EventManager::create(const RealCoordinate& coord) {
    int id;
    if (available_stack.size()) {
        id = available_stack.back(); available_stack.pop_back();
        events[id].coord = coord;
        events[id].associated_arc = nullptr;
        // intersect_point will contain garbage, it should not be accessed
        // given that associated_arc is set to nullptr
    }   
    else {
        id = events.size();
        events.emplace_back(coord);
    }
    return id;
}

int EventManager::create(
    const RealCoordinate& coord, const RealCoordinate& intersect, 
    Arc* associated_arc
) {
    int id;
    if (available_stack.size()) {
        id = available_stack.back(); available_stack.pop_back();
        events[id].coord = coord;
        events[id].associated_arc = associated_arc;
        events[id].intersect_point = intersect;
    }
    else {
        id = events.size();
        events.emplace_back(coord, intersect, associated_arc);
    }
    return id;
}

void EventManager::remove(int id) {
    available_stack.push_back(id);
}

void EventQueue::insert(int event_id) {
    event_id_heap.push_back(event_id);
    if (event_id >= id_to_location.size()) {
        int new_size = id_to_location.size() * 2;
        while (new_size <= event_id) { new_size *= 2; }
        id_to_location.resize(new_size, -1);
    }
    id_to_location[event_id] = event_id_heap.size() - 1;
    up_heapify(event_id_heap.size() - 1);
}

void EventQueue::remove(int event_id) {
    int heap_id = id_to_location[event_id];
    if (heap_id == -1) { return; } 
    id_to_location[event_id] = -1;
    event_id_heap[heap_id] = event_id_heap.back();
    id_to_location[event_id_heap.back()] = heap_id;
    event_id_heap.pop_back();
    if (compare_event_id(event_id, event_id_heap[heap_id])) {
        up_heapify(heap_id);
    }
    else {
        down_heapify(heap_id);
    }
}

int EventQueue::consume_next() {
    int event_id = event_id_heap[0];
    swap(0, event_id_heap.size() - 1);
    event_id_heap.pop_back();
    id_to_location[event_id] = -1;
    down_heapify(0);
    return event_id;
}

bool EventQueue::empty() {
    return event_id_heap.size() == 0;
}

int EventQueue::lchild(int id) {
    return 2 * id + 1;
}

int EventQueue::rchild(int id) {
    return 2 * id + 2;
}

int EventQueue::parent(int id) {
    return (id - 1 - (id & 1 ^ 1)) >> 1;
}

void EventQueue::up_heapify(int id) {
    int parent_id = parent(id);
    while (id > 0 && compare(parent_id, id)) {
        swap(parent_id, id);
        id = parent_id;
        parent_id = parent(id);
    }
}

void EventQueue::down_heapify(int id) {
    while (id < event_id_heap.size()) {
        int id_of_lowest = id;
        int lc_id = lchild(id);
        int rc_id = rchild(id);
        int size = event_id_heap.size();
        if (lc_id < size && compare(id_of_lowest, lc_id)) { id_of_lowest = lc_id; }
        if (rc_id < size && compare(id_of_lowest, rc_id)) { id_of_lowest = rc_id; }
        if (id_of_lowest != id) {
            swap(id_of_lowest, id);
            id = id_of_lowest;
        }
        else {
            break;
        }
    }
}

void EventQueue::swap(int ida, int idb) {
    int event_id_a = event_id_heap[ida];
    int event_id_b = event_id_heap[idb];
    id_to_location[event_id_a] = idb;
    id_to_location[event_id_b] = ida;
    event_id_heap[ida] = event_id_b;
    event_id_heap[idb] = event_id_a;
 }

// this might be really inefficient
bool EventQueue::compare(int ida, int idb) {
    int event_ida = event_id_heap[ida];
    int event_idb = event_id_heap[idb];
    return compare_event_id(event_ida, event_idb);
}

bool EventQueue::compare_event_id(int event_ida, int event_idb) {
    const RealCoordinate& ca = em->get(event_ida).coord;
    const RealCoordinate& cb = em->get(event_idb).coord;
    return ca.x > cb.x || (ca.x == cb.x && ca.y > cb.y);
}

void EventQueue::print_ordered_x() {
    std::queue<int> index_queue;
    index_queue.push(0);
    while (!index_queue.empty()) {
        int id = index_queue.front(); index_queue.pop();
        std::cout << em->get(event_id_heap[id]).coord.x << "\n";
        int lc_id = lchild(id);
        int rc_id = rchild(id);
        if (lc_id < event_id_heap.size()) { index_queue.push(lc_id); }
        if (rc_id < event_id_heap.size()) { index_queue.push(rc_id); }
    }
}

} // namespace Impl