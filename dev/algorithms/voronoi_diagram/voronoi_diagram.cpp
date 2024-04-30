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
    vertices = generator.vertex_graph();
    regions = generator.region_graph();
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
    vertices = generator.vertex_graph();
    regions = generator.region_graph();
    
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

std::vector<Node*> VoronoiDiagram::consume_region_graph() {
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

Arc* BeachLine::find_intersected_arc(const RealCoordinate& c) {
    Arc* node = head;
    double y;

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

void site_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& event_queue,
    std::vector<Impl::BoundaryRay*>& rays,
    BeachLine& beach_line,
    const RealCoordinate& focus
) {
    Arc* region = beach_line.find_intersected_arc(focus);
    Arc* new_region = new Arc{focus};
    Arc* split_region = new Arc{region->focus};

    beach_line.insert_arc_above(region, new_region);
    beach_line.insert_arc_above(new_region, split_region);

    rays.push_back(new Impl::BoundaryRay(focus, region->focus));
    split_region->upper_ray = region->upper_ray;
    split_region->lower_ray = rays.back();
    region->upper_ray = rays.back();
    new_region->upper_ray = rays.back();
    new_region->lower_ray = rays.back();

    Arc* upper_upper = split_region->upper;
    Arc* lower_lower = region->lower;
    if (upper_upper) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, split_region->focus, upper_upper->focus
        );
        if (focus.y < intersect.y) {
           event_queue.push(new_intersection_event(intersect, split_region)); 
        }
    } 
    if (lower_lower) {
        RealCoordinate intersect = triangle_circumcenter(
            focus, region->focus, lower_lower->focus
        );
        if (focus.y > intersect.y) {
            event_queue.push(new_intersection_event(intersect, region));
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
    else if (focus_1.y == focus_2.y) {

    }
    else if (focus_1.x > focus_2.x) {
        return (intersection.y - focus_1.y) * (focus_2.y - focus_1.y);
    }
    else {
        return (intersection.y - focus_2.y) * (focus_1.y - focus_2.y);
    }
}

void intersection_event(
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>>& event_queue,
    std::vector<Impl::BoundaryRay*>& rays,
    BeachLine& beach_line,
    const Event& event
) {
    const RealCoordinate& intersect = *event.intersect_point;
    Arc* arc = event.associated_arc;
    Arc* u_arc = arc->upper;
    Arc* l_arc = arc->lower;
    
    int uv = ray_direction(intersect, arc->focus, u_arc->focus, l_arc->focus) > 0;
    int lv = ray_direction(intersect, arc->focus, l_arc->focus, u_arc->focus) > 0;

    arc->upper_ray->v[uv] = new RealCoordinate(intersect);
    arc->lower_ray->v[lv] = new RealCoordinate(intersect);

    beach_line.remove_arc(arc);

    const RealCoordinate& fu = u_arc->focus;
    const RealCoordinate& fl = l_arc->focus;
    rays.push_back(new BoundaryRay(fu, fl));
    rays.back()->v[(fl.y > fu.y)] = new RealCoordinate(intersect);
    u_arc->lower_ray = rays.back();
    l_arc->upper_ray = rays.back();
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

void FortunesAlgorithm::compute(std::vector<RealCoordinate>& seeds, double min, double max) {
    using namespace Impl; 
    num_seeds = seeds.size();
    rays = {};
    rays.reserve(3 * num_seeds);
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> event_queue; 
    for (const RealCoordinate& seed : seeds) { event_queue.emplace(seed); }
    
    const RealCoordinate s1 = event_queue.top().coord; event_queue.pop();
    BeachLine beach_line(s1);

    while (!event_queue.empty()) {
        const Event event = event_queue.top(); event_queue.pop();
        if (event.intersect_point == nullptr) {
            site_event(event_queue, rays, beach_line, event.coord);
        }
        else if (event.associated_arc->active) {
            intersection_event(event_queue, rays, beach_line, event);
        }
    }
    for (BoundaryRay* ray : rays) {
        ray->clip_to_bbox(min, max, min, max);
    }
}


std::vector<Node*> FortunesAlgorithm::vertex_graph() {
    std::vector<Node*> graph;
    graph.reserve(2 * num_seeds);
    std::unordered_map<RealCoordinate, Node*> vertex_map;
    //int counter = 0;
    for (BoundaryRay* ray : rays) {
        if (ray->v[0] == nullptr || ray->v[1] == nullptr) {
            //++counter;
            continue;
        }
        Node* vertex_a = nullptr;
        Node* vertex_b = nullptr;
        const RealCoordinate& va = *ray->v[0];
        const RealCoordinate& vb = *ray->v[1]; 
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


std::vector<Node*> FortunesAlgorithm::region_graph() {
    std::vector<Node*> graph;
    graph.reserve(num_seeds);
    std::unordered_map<RealCoordinate, Node*> region_map;
    for (BoundaryRay* ray : rays) {
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

} // namespace Impl
