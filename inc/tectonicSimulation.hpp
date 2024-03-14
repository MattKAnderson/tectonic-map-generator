#include <random>

#include "Vector2D.hpp"
#include "plate.hpp"

class tectonicSimulation {
public:
    tectonicSimulation(int n_plates, int size_x, int size_y, bool custom_seed=false, int seed=0);

    int get_sizex();
    int get_sizey();
    std::vector<std::vector<int>> plate_idmap();
    std::vector<std::pair<int, int>> plate_areas();
    void plate_assignment();
    void advance(); 

    //static const double base_oceanic_crust_density = 3.0;
    //static const double base_continental_crust_density = 2.7;
    //static const double base_upper_mantle_density = 3.35;

private:
    int n_plates;
    int size_x;
    int size_y;
    int seed;
    bool custom_seed;
    std::vector<plate> plates;
    std::vector<Vector2D<double>> fractional_plate_movements;
    std::vector<std::vector<int>> plate_affiliation; 
    std::vector<std::vector<int>> updated_plate_affiliation;
    std::vector<std::vector<int>> interaction;
    std::vector<std::vector<int>> updated_interaction;
    std::vector<std::vector<double>> heights;
    std::vector<std::vector<double>> updated_heights;
    std::vector<std::vector<double>> densities;
    std::vector<std::vector<double>> updated_densities;

    std::mt19937_64 rng;

    void seed_rng();
    void set_plates_in_motion();
    void resolve_movement(int x1, int y1, int x2, int y2);
    //void plate_assignment();
    std::vector<std::pair<Vector2D<int>, int>> unassigned_adjacent_grid_cells(Vector2D<int>& cell);
    std::vector<Vector2D<int>> adjacent_ids(Vector2D<int>& cell);
};