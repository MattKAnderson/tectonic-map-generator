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
    seeds = generate_seeds(nseeds, xsize, ysize);
    //grid_algorithm(seeds, xsize, ysize);
    
    Impl::FortunesAlgorithm generator;
    generator.compute(seeds, 0.0, xsize);
    vertices = generator.consume_vertex_graph();
    //regions = generator.region_graph();
}

void VoronoiDiagram::generate_from_seeds(
    std::vector<RealCoordinate>& seeds, int xsize, int ysize
) {
    nregions = seeds.size();
    //grid_algorithm(seeds, xsize, ysize); 
    this->seeds = seeds;
    Impl::FortunesAlgorithm generator;
    generator.compute(seeds, 0.0, xsize);
    std::cout << "about to produce vertex graph" << std::endl;
    vertices = generator.consume_vertex_graph();
    //regions = generator.region_graph();
    
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
/*
std::vector<Node*> VoronoiDiagram::consume_region_graph() {
    return std::move(regions);
}
*/
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

/*
void BoundaryRay::clip_to_bbox(
    double xmin, double xmax, double ymin, double ymax
) {
    if (v[0] == nullptr && v[1] == nullptr) {
        std::cout << "Both were nullptrs... shouldnt happen" << std::endl;
    }
    if (v[0] == nullptr) {
        if (v[1]->x < xmin || v[1]->x > xmax || v[1]->y < ymin || v[1]->y > ymax) {
            delete v[1];
            v[1] = nullptr;
        }
        else {
            v[0] = clip_infinite_ray(xmin, ymin, v[1]->x, v[1]->y, ymin, ymax);
        }
    }
    else if (v[1] == nullptr) {
        if (v[0]->x < xmin || v[0]->x > xmax || v[0]->y < ymin || v[0]->y > ymax) {
            delete v[0];
            v[0] = nullptr;
        }
        else {
            v[1] = clip_infinite_ray(xmax, ymax, v[0]->x, v[0]->y, ymin, ymax);        
        }
    }
    else {
        clip_vertex_to_bbox(v[0], v[1], xmin, xmax, ymin, ymax);
        if (v[0] == nullptr) {
            delete v[1];
            v[1] = nullptr;
        }
        else {
            clip_vertex_to_bbox(v[1], v[0], xmin, xmax, ymin, ymax);
        }
    }
}
*/
/*
RealCoordinate* BoundaryRay::clip_infinite_ray(
    double x_int, double y_int, double x0, double y0, double ymin, double ymax
) {
    if (r1.x == r2.x) {
        y_int = r1.y;
    }
    else if (r1.y == r2.y) {
        x_int = r1.x;
    }
    else {
        RealCoordinate mp = {0.5 * (r1.x + r2.x), 0.5 * (r1.y + r2.y)};
        double m = (y0 - mp.y) / (x0 - mp.x);
        double b = y0 - m * x0;
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
    return new RealCoordinate{x_int, y_int};
}
*/
/*
void BoundaryRay::clip_vertex_to_bbox(
    RealCoordinate* v, RealCoordinate* other, double xmin, double xmax, 
    double ymin, double ymax
) {
    double m = (v->y - other->y) / (v->x - other->x);
    double b = v->y - m * v->x;
    if (v->x < xmin && other->x > xmin) {
        double y_int = m * xmin + b;
        if (y_int >= ymin && y_int <= ymax) {
            *v = {xmin, y_int};
            return;
        }
    }
    if (v->x > xmax && other->x < xmax) {
        double y_int = m * xmax + b;
        if (y_int >= ymin && y_int <= ymax) {
            *v = {xmax, y_int};
            return;
        }
    }
    if (v->y < ymin && other->y > ymin) {
        double x_int = (ymin - b) / m;
        if (x_int >= xmin && x_int <= xmax) {
            *v = {x_int, ymin};
            return;
        }
    }
    if (v->y > ymax && other->y < ymax) {
        double x_int = (ymax - b) / m;
        if (x_int >= xmin && x_int <= xmax) {
            *v = {x_int, ymax};
            return;
        }
    }
    if (v->x > xmax || v->x < xmin || v->y > ymax || v->y < ymin) {
        if (v == this->v[0]) {
            delete this->v[0];
            this->v[0] = nullptr;
        }
        else {
            delete this->v[1];
            this->v[1] = nullptr;
        }
    }
}
*/

Arc* BeachLine::find_intersected_arc(const RealCoordinate& c) {
    Arc* node = head;
    double y;
    //std::cout << "c: (" << c.x << ", " << c.y << ")" << std::endl;
    // TODO: handle case when c.y == y
    //       it is super rare but I've seen it at 100k seeds

    while (true) {
        if (node->left != nullptr) {
            y = parabolae_y_intercept(c.x, node->upper->focus, node->focus);
            if (c.y > y) {
                node = node->left;
                continue;
            } 
        }
        if (node->right != nullptr) {
            y = parabolae_y_intercept(c.x, node->focus, node->lower->focus);
            if (c.y < y) {
                node = node->right;
                continue;
            }
        }
        return node;
    }
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
    // TODO: insert_balance(new_arc);
}

void BeachLine::insert_balance(Arc* arc) {
    // TODO 
}

void BeachLine::delete_balance(Arc* arc) {
    // TODO
}

void BeachLine::rotate_left(Arc* arc) {
    Arc* rchild = arc->right;
    rchild->parent = arc->parent;
    if (arc->parent->left == arc) { arc->parent->left = rchild; }
    else { arc->parent->right = rchild; }
    arc->parent = rchild;
    arc->right = rchild->left;
    rchild->left = arc;
    rchild->red = arc->red;
    arc->red = true;
}

void BeachLine::rotate_right(Arc* arc) {
    Arc* lchild = arc->left;
    lchild->parent = arc->parent;
    if (arc->parent->left == arc) { arc->parent->left = lchild; }
    else { arc->parent->right = lchild; }
    arc->parent = lchild;
    arc->left = lchild->right;
    lchild->right = arc;
    lchild->red = arc->red;
    arc->red = true;
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
    arc->active = false;
    closed_regions.push_back(arc);
    // TODO: delete_balance(y);
}

void BeachLine::reserve(int n) {
    closed_regions.reserve(n);
}

Event new_intersection_event(RealCoordinate& intersect, Arc* closing_region) {
    double dist = euclidean_distance(closing_region->focus, intersect);
    Event event({intersect.x + dist, intersect.y});
    event.intersect_point = new RealCoordinate(intersect);
    event.associated_arc = closing_region;
    return event;
}

void FortunesAlgorithm::site_event(const RealCoordinate& focus) {
    regions.push_back(new Region{focus});
    Arc* arc = beach_line.find_intersected_arc(focus);
    Arc* new_arc = new Arc{focus, regions.back()};
    Arc* split_arc = new Arc{arc->focus, arc->region};

    beach_line.insert_arc_above(arc, new_arc);
    beach_line.insert_arc_above(new_arc, split_arc);

    new_arc->upper_edge = new HalfEdge{new_arc->region};
    new_arc->lower_edge = new_arc->upper_edge;
    split_arc->upper_edge = arc->upper_edge;
    arc->upper_edge = new HalfEdge{arc->region};
    split_arc->lower_edge = arc->upper_edge;
    arc->upper_edge->twin = new_arc->upper_edge;
    new_arc->lower_edge->twin = arc->upper_edge;
    
    Arc* upper_upper = split_arc->upper;
    Arc* lower_lower = arc->lower;
    if (upper_upper) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, split_arc->focus, upper_upper->focus
        );
        if (focus.y < intersect.y) {
           event_queue.push(new_intersection_event(intersect, split_arc)); 
        }
    } 
    if (lower_lower) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, arc->focus, lower_lower->focus
        );
        if (focus.y > intersect.y) {
            event_queue.push(new_intersection_event(intersect, arc));
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
    const RealCoordinate& intersect = *event.intersect_point;
    Arc* arc = event.associated_arc;
    Arc* u_arc = arc->upper;
    Arc* l_arc = arc->lower;
    
    vertices.push_back(new VertexNode(intersect));

    /*
     *  Double check this logic tomorrow
     */
    arc->upper_edge->origin = vertices.back();
    arc->upper_edge->prev = arc->lower_edge;
    arc->lower_edge->next = arc->upper_edge;
    HalfEdge* new_upper_half_edge = new HalfEdge{u_arc->region};
    HalfEdge* new_lower_half_edge = new HalfEdge{l_arc->region};
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
                event_queue.push({{known_at_x, new_intersect.y}, new_intersect, u_arc});
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
                event_queue.push({{known_at_x, new_intersect.y}, new_intersect, l_arc});
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
    //rays = {};
    regions = {};
    regions.reserve(num_seeds);
    //rays.reserve(3 * num_seeds);
    event_queue = {};
    for (const RealCoordinate& seed : seeds) { event_queue.emplace(seed); }
    const RealCoordinate s1 = event_queue.top().coord; event_queue.pop();
    regions.push_back(new Region{s1});
    beach_line = BeachLine(s1, regions.back());
    while (!event_queue.empty()) {
        const Event event = event_queue.top(); event_queue.pop();
        if (event.intersect_point == nullptr) { site_event(event.coord); }
        else if (event.associated_arc->active) { intersection_event(event); }
    }

    std::cout << "Before getting size of vertices" << std::endl;
    std::cout << "Size of vertices: " << vertices.size() << std::endl;

    //for (BoundaryRay* ray : rays) {
    //    ray->clip_to_bbox(min, max, min, max);
    //}
    //add_vertices_for_bounds_corners(); 
    //connect_vertices_on_bounds();
    //for (Region* region : regions) {
    //    region->prune_rays();
    //    region->draw_rays_on_bounds(min, max);
    //}
}

std::vector<VertexNode*> FortunesAlgorithm::consume_vertex_graph() {
    std::cout << "Consuming vertex graph" << std::endl;
    std::cout << "Size of regions: " << regions.size() << std::endl;
    for (Region* region : regions) {
        HalfEdge* edge_ptr = region->an_edge;
        if (region->an_edge == nullptr) {
            continue;
        }
        int counter = 0;
        while (edge_ptr->next != region->an_edge && counter != 1000) {
            edge_ptr->origin->connected.push_back(edge_ptr->next->origin);
            ++counter;
        }
        if (counter == 1000) {
            std::cout << "Hit counter of 1000" << std::endl;
        }
        edge_ptr->origin->connected.push_back(edge_ptr->next->origin);
    }
    return std::move(vertices);
}

std::vector<RegionNode*> FortunesAlgorithm::consume_region_graph() {
    return region_graph_from_regions();
}

std::vector<RegionNode*> FortunesAlgorithm::region_graph_from_regions() {
    std::vector<RegionNode*> nodes;
    std::unordered_map<Region*, RegionNode*> node_map;
    nodes.reserve(regions.size());
    for (Region* region : regions) {
        nodes.emplace_back();
        node_map.insert({region, nodes.back()});
    }
    for (Region* region : regions) {
        RegionNode* this_region = node_map[region];
        HalfEdge* edge_ptr = region->an_edge;
        while (edge_ptr->next != region->an_edge) {
            this_region->vertices.push_back(edge_ptr->origin->coord);
            this_region->adjacent.push_back(node_map[edge_ptr->twin->region]);
            edge_ptr = edge_ptr->next;
        }
        this_region->vertices.push_back(edge_ptr->origin->coord);
    }
    return std::move(nodes);
}

/*
void FortunesAlgorithm::add_vertices_for_bounds_corners() {
    // need to iterate through to ensure that a duplicate is not being added here
    // this would be more efficient to do from the connect vertices on bounds method
    vertices.push_back(new VertexNode({min, min}));
    vertices.push_back(new VertexNode({min, max}));
    vertices.push_back(new VertexNode({max, max}));
    vertices.push_back(new VertexNode({max, min}));
}

bool FortunesAlgorithm::compare_bounds_vertices(
    VertexNode* va, VertexNode* vb
) {
    const RealCoordinate& a = va->coord;
    const RealCoordinate& b = vb->coord;
    if (a.x == min) { return (b.x != min) || (a.y < b.y); }
    else if (a.y == max) { return (b.x != min) || (b.y == max && b.x > a.x); } 
    else if (a.x == max) { return (b.y == min && b.x != min) || (b.x == max && a.y > b.y); }
    else { return b.y == min && b.x < a.x && b.x != min; }
}

void FortunesAlgorithm::connect_vertices_on_bounds() {
    std::vector<VertexNode*> bounds_vertices;
    for (VertexNode* vertex : vertices) {
        const RealCoordinate& c = vertex->coord;
        if (c.x == min || c.x == max || c.y == min || c.y == max) {
            bounds_vertices.push_back(vertex);
        }
    }
    std::sort(bounds_vertices.begin(), bounds_vertices.end(), compare_bounds_vertices);
    for (int i = 0; i < bounds_vertices.size() - 1; i++) {
        bounds_vertices[i]->connected.push_back(bounds_vertices[i + 1]);
        bounds_vertices[i + 1]->connected.push_back(bounds_vertices[i]);
    }
    bounds_vertices[0]->connected.push_back(bounds_vertices.back());
    bounds_vertices.back()->connected.push_back(bounds_vertices[0]);
}
*/
/*
void Region::add_ray(BoundaryRay* ray) {
    rays.push_back(ray);
}

void Region::add_adjacent(Region* region) {
    adjacent.push_back(region);
}

void Region::prune_rays() {
    std::vector<BoundaryRay*> pruned_rays;
    for (BoundaryRay* ray : rays) {
        if (ray->v[0] != nullptr && ray->v[1] != nullptr) {
            pruned_rays.push_back(ray);
        }
    }
    rays = std::move(pruned_rays);
}

void delete_if_exists_add_if_not(
    std::vector<VertexNode*>& nodes,
    std::vector<RealCoordinate>& direction,
    VertexNode* node,
    VertexNode* prev
) {
    if (auto s = std::find(nodes.begin(), nodes.end(), node); s == nodes.end()) {
        nodes.push_back(node);
        direction.emplace_back(
            node->coord.x - prev->coord.x, node->coord.y - prev->coord.y
        );
    }
    else {
        int index = std::distance(s, nodes.begin());
        nodes[index] = nodes.back();
        nodes.pop_back();
        direction[index] = direction.back();
        direction.pop_back();
    }
}

void Region::draw_rays_on_bounds(double min, double max) {
    std::vector<VertexNode*> boundary_nodes;
    std::vector<RealCoordinate> direction;
    for (BoundaryRay* ray : rays) {
        const RealCoordinate& v0 = ray->v[0]->coord;
        const RealCoordinate& v1 = ray->v[1]->coord;
        if (v0.x == min || v0.x == max || v0.y == min || v0.y == max) {
            delete_if_exists_add_if_not(
                boundary_nodes, direction, ray->v[0], ray->v[1]
            );
        }   
        if (v1.x == min || v1.x == max || v1.y == min || v1.y == max) {
            delete_if_exists_add_if_not(
                boundary_nodes, direction, ray->v[1], ray->v[0]
            );
        }
    }
    if (boundary_nodes.size() == 0) { return; }

}

VertexNode* clone_if_not_in_map(
    std::unordered_map<VertexNode*, VertexNode*>& map,
    VertexNode* node
) {
    auto s = map.find(node);
    if (auto s = map.find(node); s == map.end()) {
        VertexNode* clone = new VertexNode(node->coord);
        map.insert({node, clone});
        return clone;
    }
    else {
        return s->second;
    }
}

std::vector<RealCoordinate> Region::get_bounds() {
    std::vector<RealCoordinate> bounds_vertices;
    std::unordered_map<VertexNode*, VertexNode*> clone_map;
    for (BoundaryRay* ray : rays) {
        if (ray->v[0] == nullptr || ray->v[1] == nullptr) {
            continue;
        }
        VertexNode* v0 = clone_if_not_in_map(clone_map, ray->v[0]);
        VertexNode* v1 = clone_if_not_in_map(clone_map, ray->v[1]);
        v0->connected.push_back(v1);
        v1->connected.push_back(v0);
    }
    VertexNode* node = clone_map.begin()->second;
    VertexNode* end = node;
    VertexNode* next = nullptr;
    const RealCoordinate& c1 = node->connected[0]->coord;
    const RealCoordinate& c2 = node->connected[1]->coord;
    if (c1.x > c2.x || (c1.x == c2.x && c1.y > c2.y)) {
        next = node->connected[0];
    }
    else {
        next = node->connected[1];
    }
    while (node != end) {
        bounds_vertices.push_back(node->coord);
        if (next->connected[0] == node) {
            node = next;
            next = next->connected[1];
        }
        else {
            node = next;
            next = next->connected[0];
        }
    }
    return bounds_vertices;
}

RegionNode* clone_if_not_in_map(
    std::unordered_map<Region*, RegionNode*>& map,
    Region* region
) {
    auto s = map.find(region);
    if (s == map.end()) {
        RegionNode* node = new RegionNode;
        node->vertices = region->get_bounds();
        map.insert({region, node});
        return node;
    }
    else {
        return s->second;
    }
}

void add_adjacency_if_dne(RegionNode* node_a, RegionNode* node_b) {
    auto adj_list = node_a->adjacent;
    if (std::find(adj_list.begin(), adj_list.end(), node_b) != adj_list.end()) {
        node_a->adjacent.push_back(node_b);
        node_b->adjacent.push_back(node_a);
    }
}

std::vector<RegionNode*> FortunesAlgorithm::region_graph_from_regions() {
    std::vector<RegionNode*> nodes;
    std::unordered_map<Region*, RegionNode*> region_to_node;
    for (Region* region : regions) {
        RegionNode* rn = clone_if_not_in_map(region_to_node, region);
        for (Region* adj_region : region->adjacent) {
            RegionNode* arn = clone_if_not_in_map(region_to_node, adj_region);
            add_adjacency_if_dne(rn, arn);
        }
    }
    return nodes;
}
*/
/*
std::vector<VertexNode*> FortunesAlgorithm::vertex_graph() {
    std::vector<VertexNode*> graph;
    std::vector<RealCoordinate> x_min_points, x_max_points, y_min_points, y_max_points;
    std::unordered_map<RealCoordinate, VertexNode*> vertex_map;
    for (BoundaryRay* ray : rays) {
        if (ray->v[0] == nullptr || ray->v[1] == nullptr) {
            continue;
        }
        for (int i = 0; i < 2; i++) {
            if (ray->v[i]->x == min) { x_min_points.push_back(*ray->v[i]); }
            else if (ray->v[i]->x == max) { x_max_points.push_back(*ray->v[i]); }
            if (ray->v[i]->y == min) { y_min_points.push_back(*ray->v[i]); }
            else if (ray->v[i]->y == max) { y_max_points.push_back(*ray->v[i]); }
        }
        VertexNode* vna = nullptr;
        VertexNode* vnb = nullptr;
        const RealCoordinate& va = *ray->v[0];
        const RealCoordinate& vb = *ray->v[1]; 
        if (auto s = vertex_map.find(va); s == vertex_map.end()) {
            vna = new VertexNode(va);
            graph.push_back(vna);
            vertex_map.emplace(va, graph.back());
        }
        else {
            vna = s->second;
        }
        if (auto s = vertex_map.find(vb); s == vertex_map.end()) {
            graph.push_back(new VertexNode(vb));
            vertex_map.emplace(vb, graph.back());
            vnb = graph.back();
        }
        else {
            vnb = s->second;
        }
        vna->connected.push_back(vnb);
        vnb->connected.push_back(vna);
    }
    // need to handle cases where a vertex does exist at one of these corners 
    // this is rare but possible
    graph.push_back(new VertexNode({min, min}));
    vertex_map.emplace(RealCoordinate{min, min}, graph.back());
    graph.push_back(new VertexNode({min, max}));
    vertex_map.emplace(RealCoordinate{min, max}, graph.back());
    graph.push_back(new VertexNode({max, min}));
    vertex_map.emplace(RealCoordinate{max, min}, graph.back());
    graph.push_back(new VertexNode({max, max}));
    vertex_map.emplace(RealCoordinate{max, max}, graph.back());
    std::sort(x_min_points.begin(), x_min_points.end(), compare_y_lower);
    std::sort(x_max_points.begin(), x_max_points.end(), compare_y_greater);
    std::sort(y_min_points.begin(), y_min_points.end(), compare_x_greater);
    std::sort(y_max_points.begin(), y_max_points.end(), compare_x_lower);
    std::vector<RealCoordinate> boundary_points;
    boundary_points.reserve(
        x_min_points.size() + y_min_points.size() + x_max_points.size()
        + y_max_points.size() + 4
    );
    boundary_points.insert(boundary_points.end(), x_min_points.begin(), x_min_points.end());
    boundary_points.push_back({min, max});
    boundary_points.insert(boundary_points.end(), y_max_points.begin(), y_max_points.end());
    boundary_points.push_back({max, max});
    boundary_points.insert(boundary_points.end(), x_max_points.begin(), x_max_points.end());
    boundary_points.push_back({max, min});
    boundary_points.insert(boundary_points.end(), y_min_points.begin(), y_min_points.end());
    boundary_points.push_back({min, min});
    VertexNode* start = vertex_map[{min, min}];
    for (RealCoordinate& point : boundary_points) {
        VertexNode* next = vertex_map[point];
        start->connected.push_back(next);
        next->connected.push_back(start);
        start = next;
    }
    return graph; 
}
*/
/*
std::vector<Node*> FortunesAlgorithm::region_graph() {
    std::vector<Node*> graph;
    graph.reserve(num_seeds);
    std::unordered_map<RealCoordinate, Node*> region_map;
    for (BoundaryRay* ray : rays) {
        if (ray->v[0] == nullptr || ray->v[1] == nullptr) {
            continue;
        }
        Node* region_a = nullptr;
        Node* region_b = nullptr;
        if (auto s = region_map.find(ray->r1); s == region_map.end()) {
            region_a = new Node(ray->r1);
            graph.push_back(region_a);
            region_map.emplace(ray->r1, region_a);
        }
        else {
            region_a = s->second;
        }
        if (auto s = region_map.find(ray->r2); s == region_map.end()) {
            region_b = new Node(ray->r2);
            graph.push_back(region_b);
            region_map.emplace(ray->r2, region_b);
        }
        else {
            region_b = s->second;
        }
        region_a->edges.push_back(region_b);
        region_b->edges.push_back(region_a);
    }
    return graph;
}
*/
} // namespace Impl