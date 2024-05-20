#include <voronoi_diagram.hpp>

VertexNode* VertexGraph::get_head() { 
    if (refs.size() == 0) { return nullptr; }
    else { return refs[0]; }
}

std::vector<VertexNode*> VertexGraph::get_vertices() {
    return refs;
}

RegionNode* RegionGraph::get_head() {
    if (refs.size() == 0) { return nullptr; }
    else { return refs[0]; }
}

std::vector<RegionNode*> RegionGraph::get_regions() {
    return refs;
}

VoronoiDiagram::VoronoiDiagram(int seed) {
    rng = std::mt19937_64(seed);
    generator = new Impl::FortunesAlgorithm;
}

void VoronoiDiagram::generate(int xsize, int ysize, int nseeds) {
    seeds = generate_seeds(nseeds, xsize, ysize);
    generate_from_seeds(seeds, xsize, ysize);
}

void VoronoiDiagram::generate_from_seeds(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    this->seeds = seeds;
    generator->compute(seeds, 0.0, xsize);
    //generator->bound();
    generator->crop(0.0, xsize);
}

std::vector<std::vector<int>> VoronoiDiagram::get_diagram() {
    return diagram;
}

std::vector<RealCoordinate> VoronoiDiagram::get_seeds() {
    return seeds;
}

VertexGraph VoronoiDiagram::get_vertex_graph() {
    return generator->get_vertex_graph();
}

RegionGraph VoronoiDiagram::get_region_graph() {
    return generator->get_region_graph();
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
    // TODO: delete_balance(y);
}

Event new_intersection_event(RealCoordinate& intersect, Arc* closing_region) {
    double dist = euclidean_distance(closing_region->focus, intersect);
    Event event({intersect.x + dist, intersect.y}, intersect, closing_region);
    return event;
}

void FortunesAlgorithm::site_event(const RealCoordinate& focus) {
    regions[next_region_id] = {focus, nullptr};
    Arc* arc = beach_line.find_intersected_arc(focus);
    Arc* new_arc = beach_line.new_arc(focus, &regions[next_region_id]); 
    Arc* split_arc = beach_line.new_arc(arc->focus, arc->region); 
    ++next_region_id;

    beach_line.insert_arc_above(arc, new_arc);
    beach_line.insert_arc_above(new_arc, split_arc);

    new_arc->upper_edge = new_interior_edge(new_arc->region); 
    half_edges.push_back(new_arc->upper_edge);
    new_arc->lower_edge = new_arc->upper_edge;
    new_arc->region->an_edge = new_arc->upper_edge;
    split_arc->upper_edge = arc->upper_edge;
    arc->upper_edge = new_interior_edge(arc->region); 
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
            if (split_arc->event_id != -1) {
                event_queue.remove(split_arc->event_id);
            }
            split_arc->event_id = event_id;
            event_queue.insert(event_id);
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
            if (arc->event_id != -1) {
                event_queue.remove(arc->event_id);
            }
            arc->event_id = event_id;
            event_queue.insert(event_id);
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
     
    vertices.push_back(new_interior_vertex(intersect));

    arc->upper_edge->origin_id = vertices.size() - 1;
    arc->upper_edge->prev = arc->lower_edge;
    arc->lower_edge->next = arc->upper_edge;
    HalfEdge* new_upper_half_edge = new_interior_edge(u_arc->region); 
    HalfEdge* new_lower_half_edge = new_interior_edge(l_arc->region); 
    half_edges.push_back(new_upper_half_edge);
    half_edges.push_back(new_lower_half_edge);
    new_upper_half_edge->twin = new_lower_half_edge;
    new_upper_half_edge->origin_id = vertices.size() - 1;
    new_upper_half_edge->prev = u_arc->lower_edge;
    u_arc->lower_edge->next = new_upper_half_edge;
    u_arc->lower_edge = new_upper_half_edge;
    new_lower_half_edge->twin = new_upper_half_edge;
    new_lower_half_edge->next = l_arc->upper_edge;
    l_arc->upper_edge->prev = new_lower_half_edge;
    l_arc->upper_edge->origin_id = vertices.size() - 1;
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
                event_queue.insert(event_id);
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
                event_queue.insert(event_id);
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

void FortunesAlgorithm::initialize() {
    if (regions != nullptr) {
        delete[] regions;
    }
    regions = new Region[num_seeds];
    half_edges = {};
    half_edges.reserve(2 * (num_seeds * 3 - 6));
    //if (boundary_half_edges != nullptr) {
    //    delete[] boundary_half_edges;
    //}
    //boundary_half_edges = nullptr;
    if (internal_half_edges != nullptr) {
        delete[] internal_half_edges;
    }
    internal_half_edges = new HalfEdge[6 * (num_seeds) - 12];
    if (internal_vertices != nullptr) {
        delete[] internal_vertices;
    }
    internal_vertices = new VertexNode[2 * num_seeds - 5];
    //if (exterior_vertices != nullptr) {
    //    delete[] exterior_vertices;
    //}
    //exterior_vertices = nullptr;
    vertices = {};
    vertices.reserve(3 * num_seeds - 6);
    beach_line = BeachLine();
    event_manager = EventManager();
    event_queue = EventQueue();
    event_queue.set_event_manager(&event_manager);
    next_half_edge_index = 0;
    next_vertex_index = 0;
    next_region_id = 0;
}

void FortunesAlgorithm::compute(std::vector<RealCoordinate>& seeds, double min, double max) {
    //this->min = min;
    //this->max = max;
    num_seeds = seeds.size();
    initialize();

    int event_id;
    for (const RealCoordinate& seed : seeds) { 
        event_id = event_manager.create(seed);
        event_queue.insert(event_id);
    }
    event_id = event_queue.consume_next();
    const RealCoordinate& s1 = event_manager.get(event_id).coord;
    regions[next_region_id] = {s1, nullptr};
    beach_line.set_head(beach_line.new_arc(s1, &regions[next_region_id]));
    ++next_region_id;
    event_manager.remove(event_id);
    
    while (!event_queue.empty()) {
        event_id = event_queue.consume_next();
        const Event& event = event_manager.get(event_id);
        if (event.associated_arc == nullptr) { 
            site_event(event.coord); 
        }
        else {
            intersection_event(event); 
        }
        event_manager.remove(event_id);
    }
}

RealCoordinate FortunesAlgorithm::clip_infinite_edge(
    HalfEdge* edge, double xmax, double xmin, double ymax, double ymin
) {
    const auto& [x0, y0] = vertices[edge->twin->origin_id]->coord;
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

void FortunesAlgorithm::bound() {
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

    std::list<std::pair<RealCoordinate, HalfEdge*>> exterior;
    Arc* lower = beach_line.get_lowest();
    Arc* upper = lower->upper;
    while (upper != nullptr) {
        HalfEdge* edge = upper->lower_edge->origin_id != -1 ? lower->upper_edge : upper->lower_edge;
        RealCoordinate v = clip_infinite_edge(edge, xmax, xmin, ymax, ymin);
        exterior.emplace_back(v, edge); 
        lower = upper;
        upper = upper->upper;
    } 
    connect_DCEL_exterior(exterior, xmax, xmin, ymax, ymin);
    /*
    exterior_vertices = new VertexNode[exterior.size() + 4];
    boundary_half_edges = new HalfEdge[exterior.size() * 2 + 8];
    int next_vertex_id = 0;
    int next_half_edge_id = 0;

    RealCoordinate corners[4] = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
    int corner_to_place = 0;
    for (auto it = exterior.begin(); it != exterior.end(); ++it) {
        for (int i = corner_to_place; i < 4; i++) {
            if (corners[i] == (*it).first) {
                ++corner_to_place;
                break;
            }
            else if (is_before_on_bbox_exterior(corners[i], (*it).first, xmax, xmin, ymax, ymin)) {
                exterior.insert(it, {corners[i], nullptr});
                ++corner_to_place;
            }
            else {
                break;
            }
        }
    }
    for (int i = corner_to_place; i < 4; i++) {
        exterior.emplace_back(corners[i], nullptr);
    }

    HalfEdge* prev_edge = nullptr;
    HalfEdge* first_corner_edge = nullptr;
    int first_origin_id = -1;
    for (auto [vertex, edge_out] : exterior) {
        VertexNode* vertex_node = &exterior_vertices[next_vertex_id++];
        vertex_node->coord = vertex;
        vertices.push_back(vertex_node);
        HalfEdge* new_edge = &boundary_half_edges[next_half_edge_id++];
        HalfEdge* twin_edge = &boundary_half_edges[next_half_edge_id++];
        if (edge_out == nullptr) {
            //new_edge->region = nullptr;
            new_edge->prev = prev_edge;
            new_edge->origin_id = vertices.size() - 1;
            new_edge->twin = twin_edge;
            twin_edge->twin = new_edge;
            if (first_corner_edge == nullptr) { 
                first_corner_edge = new_edge; 
                first_origin_id = vertices.size() - 1;    
            }
            if (prev_edge) {
                prev_edge->next = new_edge;
                prev_edge->twin->origin_id = vertices.size() - 1;
            }
            prev_edge = new_edge;
            half_edges.push_back(new_edge);
            half_edges.push_back(twin_edge);
        }
        else {
            edge_out->origin_id = vertices.size() - 1;
            new_edge->region = edge_out->twin->region;
            new_edge->prev = edge_out->twin;
            new_edge->origin_id = vertices.size() - 1;
            new_edge->twin = twin_edge;
            twin_edge->twin = new_edge;
            if (prev_edge) { 
                prev_edge->next = edge_out; 
                prev_edge->twin->origin_id = vertices.size() - 1;
                edge_out->prev = prev_edge;    
            }
            prev_edge = new_edge;
            edge_out->twin->next = prev_edge;
            half_edges.push_back(prev_edge);
            half_edges.push_back(twin_edge);
        }
    }
    HalfEdge* edge_out = exterior.begin()->second;
    if (edge_out == nullptr) {
        first_corner_edge->prev = prev_edge;
        prev_edge->next = first_corner_edge;
        prev_edge->twin->origin_id = first_origin_id;
    }
    else {
        HalfEdge* next = edge_out;
        prev_edge->next = next;
        prev_edge->twin->origin_id = first_origin_id;
        next->prev = prev_edge;
    }
    */
}

bool outside_bbox(const RealCoordinate& c, double min, double max) {
    return c.x < min || c.x > max || c.y < min || c.y > max;
}

bool inside_bbox(const RealCoordinate& c, double min, double max) {
    return c.x >= min && c.x <= max && c.y >= min && c.y <= max;
}

class crop_compare {
public:
    crop_compare(double min, double max): min(min), max(max) {}
    bool operator()(
        const std::pair<RealCoordinate, HalfEdge*>& a,
        const std::pair<RealCoordinate, HalfEdge*>& b
    ) {
        return is_before_on_bbox_exterior(a.first, b.first, max, min, max, min);
    }
private:
    double min, max;
};

void FortunesAlgorithm::crop(double min, double max) {
    LineClipper clipper(min, max);
    crop_compare comp(min, max);
    //int inside_count = 0;
    std::list<std::pair<RealCoordinate, HalfEdge*>> exterior;
    for (HalfEdge* edge : half_edges) {
        if (edge->origin_id == -1) {
            const RealCoordinate& twin_c = vertices[edge->twin->origin_id]->coord;
            if (inside_bbox(twin_c, min, max)) { 
                //++inside_count;
                RealCoordinate v = clip_infinite_edge(edge, max, min, max, min);
                exterior.emplace_back(v, edge);
            }
            continue;
        }
        const RealCoordinate& c = vertices[edge->origin_id]->coord;
        if (outside_bbox(c, min, max)) {
            if (edge->twin->origin_id == -1) { continue; }
            const RealCoordinate& twin_c = vertices[edge->twin->origin_id]->coord;
            if (!clipper.CohenSutherlandClip(c, twin_c)) { continue; } 
            vertices[edge->origin_id]->coord = clipper.get_clipped_a();
            exterior.emplace_back(clipper.get_clipped_a(), edge);
        }
        else if (c.x == min || c.x == max || c.y == min || c.y == max) {
            exterior.emplace_back(c, edge);
        }
        //++inside_count;
    }
    exterior.sort(comp);

    connect_DCEL_exterior(exterior, max, min, max, min);

    std::vector<HalfEdge*> cropped_half_edges;
    std::vector<VertexNode*> cropped_vertices; 
    std::vector<int> vertex_id_map(vertices.size(), -1);
    // std::vector<Region*> cropped_regions; TODO
    for (int i = 0; i < vertices.size(); i++) {
        VertexNode* vertex = vertices[i];
        if (inside_bbox(vertex->coord, min, max)) {
            vertex_id_map[i] = cropped_vertices.size();
            cropped_vertices.push_back(vertex);
        }
    }
    for (HalfEdge* edge : half_edges) {
        if (edge->origin_id == -1) { continue; }
        if (inside_bbox(vertices[edge->origin_id]->coord, min, max)) {
            cropped_half_edges.push_back(edge);
            edge->origin_id = vertex_id_map[edge->origin_id];
            if (edge->origin_id >= cropped_vertices.size()) {
                std::cout << "Origin id out of bounds: " << edge->origin_id << std::endl;
            }
            else if (edge->origin_id < 0 ) {
                std::cout << "Origin id out of bounds: " << edge->origin_id << std::endl;
            }
            if (edge->region) {
                edge->region->an_edge = edge;
            }
        }
    }
    vertices = cropped_vertices;
    half_edges = cropped_half_edges;
}

void FortunesAlgorithm::connect_DCEL_exterior(
    std::list<std::pair<RealCoordinate, HalfEdge*>>& exterior,
    double xmax, double xmin, double ymax, double ymin
) {

    VertexNode* exterior_vertices = new VertexNode[exterior.size() + 4];
    HalfEdge* boundary_half_edges = new HalfEdge[exterior.size() * 2 + 8];
    additional_vertices.push_back(exterior_vertices);
    additional_half_edges.push_back(boundary_half_edges);

    int next_vertex_id = 0;
    int next_half_edge_id = 0;

    RealCoordinate corners[4] = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
    int corner_to_place = 0;
    for (auto it = exterior.begin(); it != exterior.end(); ++it) {
        for (int i = corner_to_place; i < 4; i++) {
            if (corners[i] == (*it).first) {
                ++corner_to_place;
                break;
            }
            else if (is_before_on_bbox_exterior(corners[i], (*it).first, xmax, xmin, ymax, ymin)) {
                exterior.insert(it, {corners[i], nullptr});
                ++corner_to_place;
            }
            else {
                break;
            }
        }
    }
    for (int i = corner_to_place; i < 4; i++) {
        exterior.emplace_back(corners[i], nullptr);
    }

    HalfEdge* prev_edge = nullptr;
    HalfEdge* first_corner_edge = nullptr;
    int first_origin_id = -1;
    for (auto [vertex, edge_out] : exterior) {
        VertexNode* vertex_node = &exterior_vertices[next_vertex_id++];
        vertex_node->coord = vertex;
        vertices.push_back(vertex_node);
        HalfEdge* new_edge = &boundary_half_edges[next_half_edge_id++];
        HalfEdge* twin_edge = &boundary_half_edges[next_half_edge_id++];
        if (edge_out == nullptr) {
            //new_edge->region = nullptr;
            new_edge->prev = prev_edge;
            new_edge->origin_id = vertices.size() - 1;
            new_edge->twin = twin_edge;
            twin_edge->twin = new_edge;
            if (first_corner_edge == nullptr) { 
                first_corner_edge = new_edge; 
                first_origin_id = vertices.size() - 1;    
            }
            if (prev_edge) {
                prev_edge->next = new_edge;
                prev_edge->twin->origin_id = vertices.size() - 1;
            }
            prev_edge = new_edge;
            half_edges.push_back(new_edge);
            half_edges.push_back(twin_edge);
        }
        else {
            edge_out->origin_id = vertices.size() - 1;
            new_edge->region = edge_out->twin->region;
            new_edge->prev = edge_out->twin;
            new_edge->origin_id = vertices.size() - 1;
            new_edge->twin = twin_edge;
            twin_edge->twin = new_edge;
            if (prev_edge) { 
                prev_edge->next = edge_out; 
                prev_edge->twin->origin_id = vertices.size() - 1;
                edge_out->prev = prev_edge;    
            }
            prev_edge = new_edge;
            edge_out->twin->next = prev_edge;
            half_edges.push_back(prev_edge);
            half_edges.push_back(twin_edge);
        }
    }
    HalfEdge* edge_out = exterior.begin()->second;
    if (edge_out == nullptr) {
        first_corner_edge->prev = prev_edge;
        prev_edge->next = first_corner_edge;
        prev_edge->twin->origin_id = first_origin_id;
    }
    else {
        HalfEdge* next = edge_out;
        prev_edge->next = next;
        prev_edge->twin->origin_id = first_origin_id;
        next->prev = prev_edge;
    }
}

/*
void FortunesAlgorithm::clip_to_bbox(double min, double max) {

    for (auto* edge : half_edges) {
        const RealCoordinate& v = edge->origin->coord;
        if (inside_bbox(v, min, max)) { continue; }

    }
}
*/

VertexGraph FortunesAlgorithm::get_vertex_graph() {
    std::vector<VertexNode*> vertex_refs; vertex_refs.reserve(vertices.size());
    VertexNode* new_vertices = new VertexNode[vertices.size()];
    for (int i = 0; i < vertices.size(); i++) {
        vertex_refs.push_back(&new_vertices[i]);
        new_vertices[i].coord = vertices[i]->coord;
    }
    for (HalfEdge* half_edge : half_edges) {
        if (half_edge->twin->origin_id == -1) {
            std::cout << "Caught -1 id" << std::endl;
        }
        else if (half_edge->twin->origin_id > vertices.size()) {
            std::cout << "Caught out of bounds id" << std::endl;
        }
        new_vertices[half_edge->origin_id].connected.push_back(
            &new_vertices[half_edge->twin->origin_id]
        );
    }
    return VertexGraph(new_vertices, vertex_refs);
}

RegionGraph FortunesAlgorithm::get_region_graph() {
    std::vector<RegionNode*> nodes; nodes.reserve(num_seeds);
    RegionNode* nodes_array = new RegionNode[num_seeds];
    int next_id = 0;
    std::unordered_map<Region*, RegionNode*> node_map;
    for (int i = 0; i < num_seeds; i++) {
        nodes_array[i].adjacent.reserve(8);
        nodes_array[i].vertices.reserve(8);
        nodes.push_back(&nodes_array[i]);
        node_map.insert({&regions[i], nodes.back()});
    }
    for (int i = 0; i < num_seeds; i++) {
        RegionNode* this_region = node_map[&regions[i]];
        HalfEdge* edge_ptr = regions[i].an_edge;
        while (edge_ptr->next != regions[i].an_edge) {
            this_region->vertices.push_back(vertices[edge_ptr->origin_id]->coord);
            if (edge_ptr->twin->region != nullptr) {
                this_region->adjacent.push_back(node_map[edge_ptr->twin->region]);
            }
            edge_ptr = edge_ptr->next;
        }
        this_region->vertices.push_back(vertices[edge_ptr->origin_id]->coord);
    }
    return RegionGraph(nodes_array, nodes);
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