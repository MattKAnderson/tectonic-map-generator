#include <algorithm>
#include <chrono>
#include <unordered_set>
#include "tectonicSimulation.hpp"

tectonicSimulation::tectonicSimulation(int n_plates, int size_x, int size_y, bool custom_seed, int seed) :
    n_plates(n_plates), size_x(size_x), size_y(size_y), custom_seed(custom_seed), seed(seed)
{
    seed_rng();
    plates = std::vector<plate>(n_plates);
    fractional_plate_movements = std::vector<Vector2D<double>>(n_plates, {0.0, 0.0});
    plate_affiliation = std::vector<std::vector<int>>(size_x, std::vector<int>(size_y, -1));
    heights = std::vector<std::vector<double>>(size_x, std::vector<double>(size_y, 0.0));
    //densities = std::vector<std::vector<double>>(size_x, std::vector<double>(size_y, base_oceanic_crust_density));
    
    updated_plate_affiliation = std::vector<std::vector<int>>(size_x, std::vector<int>(size_y));
    updated_heights = std::vector<std::vector<double>>(size_x, std::vector<double>(size_y));
    updated_densities = std::vector<std::vector<double>>(size_x, std::vector<double>(size_y));
}

int tectonicSimulation::get_sizex() {
    return size_x;
}

int tectonicSimulation::get_sizey() {
    return size_y;
}

std::vector<std::vector<int>> tectonicSimulation::plate_idmap() {
    return plate_affiliation;
}

std::vector<std::pair<int, int>> tectonicSimulation::plate_areas() {
    std::vector<std::pair<int, int>> areas(plates.size());
    for (int i = 0; i < areas.size(); i++) {
        areas[i] = {i, 0};
    }
    for (int i = 0; i < plate_affiliation.size(); i++) {
        for (int k = 0; k < plate_affiliation[i].size(); k++) {
            areas[plate_affiliation[i][k]].second++;
        }
    }
    return areas;
}

void tectonicSimulation::plate_assignment() {

    // Data structure for holding the collection of grid cells forming the perimeter of a 
    // plate. Plate index corresponds to 'plates' datastructure 
    std::vector<std::vector<Vector2D<int>>> plate_perimeters(plates.size()); 

    std::uniform_int_distribution<int> random_x {0, plate_affiliation.size() - 1};
    std::uniform_int_distribution<int> random_y (0, plate_affiliation[0].size() - 1);

    // randomly place the initial points of the plates
    for (int i = 0; i < plates.size(); i++) {
        int x = random_x(rng);
        int y = random_y(rng);

        while (plate_affiliation[x][y] != -1) {
            x = random_x(rng);
            y = random_y(rng);
        }

        plate_affiliation[x][y] = i;
        plate_perimeters[i] = {{x, y}};
    }

    /*
     *  Core algorithm for assigning grid cells to plates. 
     *
     *  The method is to recall the perimeter of each plate in an array. Each iteration, 
     *  one of the non-empty plate perimeter arrays is randomly selected (i.e. which plate 
     *  to grow is randomly selected), a random perimeter grid cell is then selected from 
     *  this array. The adjacent grid cells to this grid cell are checked to see which are
     *  currently unassigned. All unassigned adjacent grid cells are candidates. From the 
     *  candidates, any with multiple adjacencies to the same plate are preferentially 
     *  chosen. One candidate from the group of candidates with maximum adjacency is chosen
     *  at random. In this way the plate should grow in a roughly blob manner without 'fingers'
     *
     */

    /*
     *  
     * Data structures to support biased random selection of a growable plate. Together forms 
     * a pair of the plate index and a weight to use when randomly selecting which plate to 
     * grow. When a plate is no longer growable, the pair is swapped to the back and popped. 
     * 
     * Using the geometric distribution, p=0.33, for the growth weights, clamped to a min of
     * 1, results in a distribution of plate areas similar to Earth.
     */
    
    std::geometric_distribution<int> growth_weight_selector(0.33);
    std::vector<int> growable_plate_indices(plate_perimeters.size());
    std::vector<int> growable_plate_weights(plate_perimeters.size());
    for (int i = 0; i < growable_plate_indices.size(); i++) {
        growable_plate_indices[i] = i;
        growable_plate_weights[i] = growth_weight_selector(rng);
        if (growable_plate_weights[i] < 1) {
            growable_plate_weights[i] = 1;
        }
        std::cout << "plate " << i << " weight: " << growable_plate_weights[i] << "\n";
    } 

    while (!growable_plate_indices.empty()) {
        std::discrete_distribution<int> plate_selector(growable_plate_weights.begin(), growable_plate_weights.end());
        int growable_plate_index = plate_selector(rng);
        int plate_index = growable_plate_indices[growable_plate_index];

        if (plate_perimeters[plate_index].size() == 0) { std::cout << "caught it" << std::endl; }
        std::uniform_int_distribution<int> cell_selector {0, plate_perimeters[plate_index].size() - 1};
        int perim_index = cell_selector(rng);

        auto candidates = unassigned_adjacent_grid_cells(plate_perimeters[plate_index][perim_index]);

        // Check that there exist candidates. If there are none, remove the grid cell from the 
        // perimeter array and continue to the next iteration. 
        if (candidates.size() == 0) {
            plate_perimeters[plate_index][perim_index] = plate_perimeters[plate_index].back();
            plate_perimeters[plate_index].pop_back();
            if (plate_perimeters[plate_index].size() == 0) {
                growable_plate_indices[growable_plate_index] = growable_plate_indices.back();
                growable_plate_indices.pop_back();
                growable_plate_weights[growable_plate_index] = growable_plate_weights.back();
                growable_plate_weights.pop_back();
            }
            continue;
        }

        int highest_adjacency = std::max_element(
            candidates.begin(), candidates.end(), 
            [](std::pair<Vector2D<int>, int>& a, std::pair<Vector2D<int>, int>& b) {
                return a.second < b.second;
            }
        )->second;

        std::vector<std::pair<Vector2D<int>, int>> filtered_candidates;
        std::copy_if(
            candidates.begin(), candidates.end(), std::back_inserter(filtered_candidates),
            [highest_adjacency](std::pair<Vector2D<int>, int>& a) { return a.second == highest_adjacency; }
        );

        std::uniform_int_distribution<int> candidate_selector {0, filtered_candidates.size() - 1};
        int candidate_index = candidate_selector(rng);
        int x = filtered_candidates[candidate_index].first.x;
        int y = filtered_candidates[candidate_index].first.y;

        plate_affiliation[x][y] = plate_index;
        plate_perimeters[plate_index].emplace_back(x, y); 
    }
}

/*
 *  advance() advances the tectonic plate simulation one timestep. 
 *
 *  At each timestep, each plate is moved. The fractional step for each plate in both the x and y
 *  dimensions is computed and accumulated. When the fractional step exceeds 1.0 in a dimension, the
 *  plate is moved by 1 pixel in that direction. Convergent boundary interactions are computed using 
 *  the velocities of the pixel under question and the adjacent pixels in the direction of motion 
 *  belonging to different plates. Divergent boundary interactions are handled retroactively by 
 *  scanning the grid for vacancies.
 * 
 *  Plate interaction situations
 *    1.  Subduction  (oceanic-oceanic and oceanic-continental)
 *            In this case one plate begins to subduct under the other. The subducting plate is 
 *            destroyed as it moves under the other and it's material recycled into the mantle. 
 *            The plate on top experiences magmatism, 100-200 km from the subduction boundary that
 *            adds material to the plate (mix continental and oceanic crust). Additionally 
 *            relamination, a process of partial metlting of the lighter components of the subducting
 *            plate and their subsequent re-application to the bottom of the overtop plate, adds more
 *            light material to the overtop plate (reducing it's average density). This process 
 *            is responsible for continental crust creation. In terms of simulation, this boundary
 *            will act as sink of the subducting plate and lower the density and add thickness to the
 *            overtop plate. Adding material to the overtop plate will happen with 2 components, a 
 *            probabilistic volcanism component and a consistent re-lamination component. 
 * 
 *            Subduction applies to oceanic-oceanic and oceanic-continental interactions
 * 
 *    2.  Oceanic rift  (any empty cell)
 *            In this case, two plates move apart from each other. This results in creation of new
 *            oceanic crust between them. There is also a chance of volcanism along this boundary
 * 
 *    3.  Transform fault
 *            In this case, two plates slide past each other. No new crust creation or destruction. 
 *            This interaction won't be recognized in the simulation
 *   
 *    4.  Rift valley
 *            In this case a continental plate is in the process of splitting. This won't be handled
 *            yet
 * 
 *    5.  Obduction zones  (oceanic-continental)
 *            In this case, plate consisting of continental crust is pushed under oceanic crust. 
 *            Presumably because the obducting plate had consisted of some oceanic plate that was 
 *            subducting, but all the oceanic crust was consumed. This situation is gravitationally
 *            unstable and should result in buckling that causes the oceanic crust to subduct under
 *            the continental crust. This can be modelled by probabilistically flipping the 
 *            obduction interaction into a subduction interaction. This necessitates remembering the
 *            interaction type at grid cells. 
 * 
 *    6.  Orogenic zone  (continental-continental)
 *            Collision of two continental plates. Increases height at grid cell. Density unaffected.
 *            Increases irregularity of terrain. 
 * 
 */
void tectonicSimulation::advance() {

    // iterate over every grid cell, check the adjacent cells in direction of motion for interactions
    // evaluate interactions and if the associated plate's fractional movement is over 1 in a 
    // dimension, move the pixel, unless the interaction is a destruction interaction. 

    for (int i = 0; i < plates.size(); i++) {
        fractional_plate_movements[i] += plates[i].drift;      
    }

    for (int i = 0; i < size_x; i++) {
        for (int k = 0; k < size_y; k++) {

            Vector2D<double> drift = plates[plate_affiliation[i][k]].drift;
            int y_dir = static_cast<int>(drift.y > 0) - static_cast<int>(drift.y < 0);
            int x_dir = static_cast<int>(drift.x > 0) - static_cast<int>(drift.x < 0);

            int x_adj = i + x_dir;
            int y_adj = k + y_dir;
            x_adj = x_adj * (x_adj < size_x) * (x_adj > -1) + (size_x - 1) * (x_adj < 0);
            y_adj = y_adj * (y_adj < size_y) * (y_adj > -1) + (size_y - 1) * (y_adj < 0);

            resolve_movement(i, k, x_adj, k);
            resolve_movement(i, k, i, y_adj);
            resolve_movement(i, k, x_adj, y_adj);

            // TODO: workout movement when fractional movement is over 1

        }
    }

}

void tectonicSimulation::resolve_movement(int x1, int y1, int x2, int y2) {

    // if they are on the same plate, check if they are to be moved
    // if they are on different plates, compare motion (direction this grid cell is in), 
    //   if the other plate has greater speed in the same direction then there is a rift 
    //   forming, do nothing on this for now
    //   else, (this plate has greater speed, or the directions are opposite)

}

void tectonicSimulation::seed_rng() {
    if (!custom_seed) {
        seed = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now().time_since_epoch()
        ).count();
    }
    rng.seed(seed);
}

// TODO parametrize the speed distribution
void tectonicSimulation::set_plates_in_motion() {
    std::uniform_real_distribution<double> plate_speed_distribution(0.1, 1.0); 
    for (plate& p : plates) {
        p.update_drift({plate_speed_distribution(rng), plate_speed_distribution(rng)});
    }
}

std::vector<std::pair<Vector2D<int>, int>> tectonicSimulation::unassigned_adjacent_grid_cells(Vector2D<int>& cell) {
    int affiliated_plate = plate_affiliation[cell.x][cell.y];
    std::vector<std::pair<Vector2D<int>, int>> unassigned_adjacent;

    auto adjacent_cells = adjacent_ids(cell); 

    for (Vector2D c : adjacent_cells) {
        if (plate_affiliation[c.x][c.y] == -1) {
            auto adj_adj_cells = adjacent_ids(c);
            int adj_count = 0;
            for (Vector2D aac : adj_adj_cells) {
                adj_count += (plate_affiliation[aac.x][aac.y] == affiliated_plate);
            }
            unassigned_adjacent.emplace_back(c, adj_count);
        }
    }

    return unassigned_adjacent;
}

std::vector<Vector2D<int>> tectonicSimulation::adjacent_ids(Vector2D<int>& cell) {   
    int above = cell.y - 1 < 0 ? size_y - 1 : cell.y - 1;
    int below = cell.y + 1 < size_y ? cell.y + 1 : 0;
    int left = cell.x - 1 < 0 ? size_x - 1 : cell.x - 1;
    int right = cell.x + 1 < size_x ? cell.x + 1 : 0;

    return {
        {left, above}, {cell.x, above}, {right, above}, {left, cell.y}, {right, cell.y},
        {left, below}, {cell.x, below}, {left, below}
    };    
}
