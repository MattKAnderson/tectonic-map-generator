#include <algorithm>
#include <concepts>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <vector>
#include <unordered_map>
#include <utility>
#include <Vector2D.hpp>

#include <fstream>

namespace space_partition {
namespace Impl {
struct GridCell {
    int x, y;
    bool operator==(const GridCell& other) const {
        return x == other.x && y == other.y;
    }
};

struct GridCellHash {
    std::size_t operator()(const GridCell& gc) const {
        return std::hash<long long>{}(
            gc.x * (INT_MAX + static_cast<long long>(1)) + gc.y
        );
    }
};

struct GridCellAdjacencyCounter {
    GridCell cell;
    int adjacency_count;
    bool operator<(const GridCellAdjacencyCounter& other) const {
        return adjacency_count < other.adjacency_count; 
    };
};

struct Node {
    GridCell cell;
    Node* edges[4] = {}; // +x, -x, +y, -y
    Node(GridCell g) : cell(g) {};
};

std::vector<GridCell> adjacent_ids(GridCell& cell, size_t size_x, size_t size_y);

std::vector<GridCellAdjacencyCounter> unassigned_adjacent_grid_cells(
    GridCell& cell, std::vector<std::vector<int>>& region_affiliation
);

template<typename URBG, typename Distribution>
concept URBGCompatibleDistribution = requires (URBG urbg, Distribution d) {
    requires std::uniform_random_bit_generator<URBG>;
    requires std::integral<typename Distribution::result_type>;
    {d(urbg)} -> std::same_as<typename Distribution::result_type>;
};

int find_parent(std::vector<int>& parent_labels, int label);

void join_labels(std::vector<int>& parent_labels, int label_a, int label_b);

std::vector<int> relabel(std::vector<int>& parent_labels);
} // namespace Impl

template<typename URBG, typename Distribution>
requires Impl::URBGCompatibleDistribution<URBG, Distribution>
class StochasticSpacePartitioner {
public:
    StochasticSpacePartitioner(URBG rng, 
                               Distribution distribution);
    std::vector<std::vector<int>> iterative_growth_partition(int xsize, 
                                                             int ysize, 
                                                             int num_regions);
    std::vector<std::vector<int>> slice_areas_partition(int xsize, 
                                                        int ysize,
                                                        int ncuts);
    
private:
    URBG rng_;
    Distribution distribution_; 
    
    void travel_until_intersect(    
        std::unordered_map<Impl::GridCell, Impl::Node*, Impl::GridCellHash>& graph,
        std::normal_distribution<double>& random_drift_change,
        Impl::GridCell pos,
        Vector2D<double> drift,
        Impl::Node* node,
        size_t xsize, 
        size_t ysize
    );
};


template<typename URBG, typename Distribution>
requires Impl::URBGCompatibleDistribution<URBG, Distribution>
StochasticSpacePartitioner<URBG, Distribution>::StochasticSpacePartitioner(
   URBG rng, Distribution distribution
) {
    rng_ = rng;
    distribution_ = distribution;
}


template<typename URBG, typename Distribution>
requires Impl::URBGCompatibleDistribution<URBG, Distribution>
std::vector<std::vector<int>> StochasticSpacePartitioner<URBG, Distribution>::
iterative_growth_partition(
    int xsize, int ysize, int num_regions
) {
    std::vector<std::vector<int>> region_affiliation(
        xsize, std::vector<int>(ysize, -1)
    );

    /*
     * Data structure for holding the collection of grid cells forming the 
     * perimeter of a region. This is used for caching the grid cells on
     * the exterior of a region that have non-assigned adjacent grid cells.
     * Enables efficient random selection of a grid cell to add the region
     */
    std::vector<std::vector<Impl::GridCell>> region_perimeters(num_regions);

    std::uniform_int_distribution<int> random_x (0, xsize - 1);
    std::uniform_int_distribution<int> random_y (0, ysize - 1);

    // random placement of the initial point for each region
    for (int i = 0; i < num_regions; i++) {
        int x = random_x(rng_);
        int y = random_y(rng_);

        while (region_affiliation[x][y] != -1) {
            x = random_x(rng_);
            y = random_y(rng_);
        }

        region_affiliation[x][y] = i;
        region_perimeters[i] = {{x, y}};
    } 

    /*
     *  A discrete distribution will pick which region to grow on each 
     *  iteration of the algorithm. These weights are randomly assigned
     *  using the weight_distribution. The weights and indices are saved
     *  so that the discrete distribution can be rebuilt every time a 
     *  region can no longer grow (rebuilt without said region). 
     */
    std::vector<int> growable_region_indices(num_regions);
    std::vector<int> growable_region_weights(num_regions);
    for (int i = 0; i < growable_region_weights.size(); i++) {
        growable_region_indices[i] = i;
        growable_region_weights[i] = distribution_(rng_);
        if (growable_region_weights[i] < 1) {
            growable_region_weights[i] = 1;
        }
    }
    std::discrete_distribution<int> region_selector(
        growable_region_weights.begin(), growable_region_weights.end()
    );

    while (!growable_region_indices.empty()) {
        int growable_index = region_selector(rng_);
        int region_index = growable_region_indices[growable_index];

        /*
        if (region_perimeters[region_index].size() == 0) {
            // no longer growable, remove from list of growable regions. This 
            // has to be done reactively because another region growing can 
            // make the final point on the 
            growable_region_indices[growable_index] = growable_region_indices.back();
            growable_region_indices.pop_back();
            growable_region_weights[growable_index] = growable_region_weights.back();
            growable_region_weights.pop_back();
            continue;
        }
        */

        std::uniform_int_distribution<int> cell_selector {
            0, region_perimeters[region_index].size() - 1
        };
        int perimeter_index = cell_selector(rng_);

        std::vector<Impl::GridCellAdjacencyCounter> candidates(
            Impl::unassigned_adjacent_grid_cells(
                region_perimeters[region_index][perimeter_index],
                region_affiliation
            )
        );

        if (candidates.size() == 0) {
            // this cell no longer has unassigned adjacent cells. Remove it from the
            // perimeter. Must be done re-actively, due to other regions growing
            auto& perimeter = region_perimeters[region_index];
            perimeter[perimeter_index] = perimeter.back();
            perimeter.pop_back();
            if (region_perimeters[region_index].size() == 0) {
                // the region is no longer growable, so remove it from the list
                // and re-build the weight distribution without this region
                std::swap(
                    growable_region_indices[growable_index],
                    growable_region_indices.back()
                );
                growable_region_indices.pop_back();
                std::swap(
                    growable_region_weights[growable_index],
                    growable_region_weights.back()
                );
                growable_region_weights.pop_back();              
                std::discrete_distribution<int> region_selector(
                    growable_region_weights.begin(), 
                    growable_region_weights.end()
                );  
            }
            continue;
        }
        else if (candidates.size() == 1) {
            Impl::GridCell& cell = candidates[0].cell;
            region_perimeters[region_index][perimeter_index] = {
                cell.x,
                cell.y
            };
            region_affiliation[cell.x][cell.y] = region_index;
        }
        else {
            int highest_adjacency = std::max_element(
                candidates.begin(), 
                candidates.end()          
            )->adjacency_count;         
            std::vector<Impl::GridCellAdjacencyCounter> filtered_candidates;
            std::copy_if(
                candidates.begin(), 
                candidates.end(), 
                std::back_inserter(filtered_candidates),
                [highest_adjacency](Impl::GridCellAdjacencyCounter& a) { 
                    return a.adjacency_count == highest_adjacency;
                }
            );  
            std::uniform_int_distribution<int> candidate_selector {
                0, filtered_candidates.size() - 1
            };
            int candidate_index = candidate_selector(rng_);
            int x = filtered_candidates[candidate_index].cell.x;
            int y = filtered_candidates[candidate_index].cell.y;

            region_affiliation[x][y] = region_index;
            region_perimeters[region_index].emplace_back(x, y); 
        }
    } // end while

    return region_affiliation;
}

template<typename URBG, typename Distribution>
requires Impl::URBGCompatibleDistribution<URBG, Distribution>
std::vector<std::vector<int>> StochasticSpacePartitioner<URBG, Distribution>::
slice_areas_partition(int xsize, int ysize, int ncuts) {
    std::vector<std::vector<int>> region_affiliation(
        xsize, std::vector<int>(ysize, -1)
    );

    // the matrix of grid cell corners is bigger than the grid by 1 in both
    // dimensions. This is why xsize/ysize are included in the distributions
    std::uniform_int_distribution<int> random_x (0, xsize);
    std::uniform_int_distribution<int> random_y (0, ysize);
    
    std::uniform_real_distribution<double> random_angle (
        -std::numbers::pi, std::numbers::pi
    );

    // TODO: tune the parameter here
    std::normal_distribution<double> random_drift_change (0.0, 0.05);


    /*
     *   Probably best to use an unordered map for a partial representation 
     *   of the graph. Only nodes that get visited get added to the unordered
     *   map. Only edges that are travelled are added. Use the GridCell as the 
     *   key, add a hash method for it. 
     * 
     *   In addition to initializing the regions, this method could probably
     *   also dump out labelled boundaries if this is updated. The work of 
     *   identifying the vertices will already be completed naturally
     */
    // std::cout << "At definition of unordered_map for the graph" << std::endl;
    std::unordered_map<Impl::GridCell, Impl::Node*, Impl::GridCellHash> graph;
    for (int i = 0; i < ncuts; i++) {

        Impl::GridCell pos(random_x(rng_), random_y(rng_));
        while (graph.find(pos) != graph.end())  {
            pos = {random_x(rng_), random_y(rng_)};
        }

        Vector2D drift(0.0, 1.0);
        drift.rotate(random_angle(rng_));

        Impl::Node* node = new Impl::Node(pos);

        // std::cout << "starting travel_until_intersect for cut " << i << std::endl;

        travel_until_intersect(
            graph, random_drift_change, pos, drift, node, xsize, ysize
        );
        travel_until_intersect(
            graph, random_drift_change, pos, -drift, node, xsize, ysize
        );

    }
    // std::cout << "Completed cut tracing" << std::endl;

    /*
     * so presumably now graph is filled with nodes connected to other nodes 
     * such that there are no dangling ends. This representation of the 
     * boundaries needs to be converted into labelling regions. 
     * 
     * Basically need to iterate left to right, top to bottom on the grid
     * monotonically incrementing labels when the adjacent is unlabelled,
     * or this cell is separated by an edge
     * 
     * Separation can be identified for all edges of a cell by getting the
     * upper left and lower right nodes and inspecting the appropriate edge.
     * 
     * Easier to just iterate through all nodes and output edges into a 
     * 2xsize x ysize matrix representing the boundaries of cells. Then
     * assigning labels to areas can proceed without continual lookups to
     * hashmap. 
     */

    //std::cout << "Number of vertices touched in cuts: " << graph.size() << std::endl;

    // std::vector<std::pair<Impl::GridCell, Impl::GridCell>> edges;
    std::vector<std::vector<int>> boundary_map(
        2 * xsize, std::vector<int>(ysize, 0)
    );
    for (auto& [cell, node] : graph) {
        int x = (cell.x - 1) + (cell.x - 1 < 0) * (xsize);
        int y = (cell.y - 1) + (cell.y - 1 < 0) * (ysize);
        boundary_map[2 * x][y] = (node->edges[0] != nullptr);
        boundary_map[2 * x + 1][y] = (node->edges[2] != nullptr);
        /*
        if (node->edges[0] != nullptr) {
            edges.push_back({{cell.x + 1, cell.y}, {cell.x, cell.y}});
        }
        if (node->edges[2] != nullptr) {
            edges.push_back({{cell.x, cell.y + 1}, {cell.x, cell.y}});
        }
        if (node->edges[1] != nullptr) {
            edges.push_back({{cell.x-1, cell.y}, {cell.x, cell.y}});
        }
        if (node->edges[3] != nullptr) {
            edges.push_back({{cell.x, cell.y - 1}, {cell.x, cell.y}});
        }*/
    }
    /*
    std::cout << "writing to diagnostic file" << std::endl;
    std::ofstream outfile;
    outfile.open("space_partition_edges.json");
    outfile << "{\"edges\": [";
    for (int i = 0; i < edges.size() - 1; i++) {
        auto& ep = edges[i];
        outfile << "[[" << ep.first.x << ", " << ep.first.y << "], ["
                << ep.second.x << ", " << ep.second.y << "]], ";
    }
    auto& ep = edges.back();
    outfile << "[[" << ep.first.x << ", " << ep.first.y << "], ["
            << ep.second.x << ", " << ep.second.y << "]]]}";
    outfile.close();
    */
    // std::cout << "Completed writing traces to boundary map" << std::endl;

    std::vector<int> parent_labels;
    for (int x = 0; x < xsize; x++) {
        for (int y = 0; y < ysize; y++) {
            int up = (x - 1) + (x - 1 < 0) * (xsize);
            int left = (y - 1) + (y - 1 < 0) * (ysize); 
            if (boundary_map[2 * x][left] == 0
                && region_affiliation[x][left] != -1) 
            {
                region_affiliation[x][y] = region_affiliation[x][left];
            }
            if (boundary_map[2 * up + 1][y] == 0
                && region_affiliation[up][y] != -1) 
            {
                if (region_affiliation[x][y] != -1) {
                    Impl::join_labels(
                        parent_labels, 
                        region_affiliation[x][y], 
                        region_affiliation[up][y]
                    );
                }
                else {
                    region_affiliation[x][y] = region_affiliation[up][y];
                }
            }
            if (region_affiliation[x][y] == -1) {
                parent_labels.push_back(parent_labels.size());
                region_affiliation[x][y] = parent_labels.back();
            }

        }
    }

    // std::cout << "Completed labelling the regions" << std::endl;

    std::vector<int> relabel_map = Impl::relabel(parent_labels);
    for (int x = 0; x < xsize; x++) {
        for (int y = 0; y < ysize; y++) {
            region_affiliation[x][y] = relabel_map[region_affiliation[x][y]];
        }
    }

    // std::cout << "Completed monotonic re-labelling of regions" << std::endl;

    return region_affiliation;    
};

template<typename URBG, typename Distribution>
requires Impl::URBGCompatibleDistribution<URBG, Distribution>
void StochasticSpacePartitioner<URBG, Distribution>::travel_until_intersect(
    std::unordered_map<Impl::GridCell, Impl::Node*, Impl::GridCellHash>& graph,
    std::normal_distribution<double>& random_drift_change,
    Impl::GridCell pos,
    Vector2D<double> drift,
    Impl::Node* node,
    size_t xsize, 
    size_t ysize
) {
    // std::cout << "Travel start: (" << pos.x << ", " << pos.y 
    //          << "), Travel drift: (" << drift.x << ", " << drift.y << ")" 
    //          << std::endl; 
    std::uniform_real_distribution<double> random_chance_threshold (0.0, 1.0);
    while (true) {
        graph[pos] = node;
        double absx = std::abs(drift.x);
        double absy = std::abs(drift.y);
        double chance_x = absx / (absx + absy);

        drift.rotate(random_drift_change(rng_));

        
        int edge_index, comp_edge_index;
        if (chance_x > random_chance_threshold(rng_)) {
            int step = (drift.x > 0) - (drift.x < 0);
            if (step == 0) {
                std::cout << "Zero X step detected!" << std::endl;
            }
            if (step > 0) {
                edge_index = 0;
                comp_edge_index = 1;
                pos = {(pos.x + 1 < xsize) * (pos.x + 1), pos.y};
            }
            else {
                edge_index = 1;
                comp_edge_index = 0;
                pos = {pos.x - 1 + (pos.x - 1 < 0) * xsize, pos.y};
            }
        }
        else {
            int step = (drift.y > 0) - (drift.y < 0);
            if (step == 0) {
                std::cout << "Zero Y step detected!" << std::endl;
            }
            if (step > 0) {
                edge_index = 2;
                comp_edge_index = 3;
                pos = {pos.x, (pos.y + 1 < ysize) * (pos.y + 1)};
            }
            else {
                edge_index = 3;
                comp_edge_index = 2;
                pos = {pos.x, pos.y - 1 + (pos.y - 1 < 0) * ysize};
            }
        }

        // branchless implementation for picking a dimension and direction 
        // to step in while handling world wrap conditions
        // comparable to what the compiler creates with the if-statements
        // but the compiler actually came up with a better solution when
        // the if statements were used.
        /*
        int go_x = (chance_x > random_chance_threshold(rng_));
        int go_y = go_x == 0;
        int step = go_x * ((drift.x > 0) - (drift.x < 0)) 
                 + go_y * ((drift.y > 0) - (drift.y < 0));
        int edge_index = (step < 0) + 2 * go_y;
        int comp_edge_index = (step > 0) + 2 * go_y;
        int x = pos.x;
        int y = pos.y;
        int new_x = go_y * x + go_x * ((step > 0) * (x + step < xsize) 
                  * (x + step) + (step < 0) * (x + step + (x + step < 0) * xsize)); 
        int new_y = go_x * y + go_y * ((step > 0) * (y + step < ysize) 
                  * (y + step) + (step < 0) * (y + step + (y + step < 0) * ysize));
        pos = {new_x, new_y};        
        */

        if (auto search = graph.find(pos); search != graph.end()) {
            // node has been visited before
            //if (search->second->edges[comp_edge_index] == node) {
            //    std::cout << "Falling back alond existing edge!" << std::endl;
            //}
            search->second->edges[comp_edge_index] = node;
            node->edges[edge_index] = search->second;
            node = search->second;
            break;
        }
        else {
            Impl::Node* new_node = new Impl::Node(pos);
            new_node->edges[comp_edge_index] = node;
            node->edges[edge_index] = node;
            node = new_node;
        }
        //std::cout << "Next pos: (" << pos.x << ", " << pos.y << std::endl;
    }
    //std::cout << "Ending on the position: (" << node->cell.x << ", "
    //          << node->cell.y << ")" << std::endl;
}

} // namespace space_partition 