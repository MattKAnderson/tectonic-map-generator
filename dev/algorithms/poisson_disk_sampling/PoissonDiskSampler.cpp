#include <PoissonDiskSampler.hpp>


PoissonDiskSampler::PoissonDiskSampler(int seed) {
    rng = std::mt19937_64(seed);
}


std::vector<Vector2D<double>> PoissonDiskSampler::sample_annulus(
    Vector2D<double> origin, double r1, double r2, double nsamples
) {
    double xmin = origin.x - r2 < 0.0 ? 0.0 : origin.x - r2;
    double xmax = origin.x + r2 > 1.0 ? 1.0 : origin.x + r2;
    double ymin = origin.y - r2 < 0.0 ? 0.0 : origin.y - r2;
    double ymax = origin.y + r2 > 1.0 ? 1.0 : origin.y + r2;
    std::uniform_real_distribution<double> x_dist(xmin, xmax);
    std::uniform_real_distribution<double> y_dist(ymin, ymax);
    
    std::vector<Vector2D<double>> samples;
    samples.reserve(nsamples);

    while (samples.size() < nsamples) {
        Vector2D<double> c = {x_dist(rng), y_dist(rng)};
        double r_sq = (c - origin).VlenSq();
        if (r_sq > r1 * r1 && r_sq < r2 * r2) {
            samples.push_back(c);
        }
    }

    return samples;
}

// the 4 corners of the 5x5 grid don't need to be searched, but leaving
// unrolling the loop to avoid these corners as an optimization to do later
bool does_not_overlap(
    std::vector<std::vector<int>>& grid,
    std::vector<Vector2D<double>>& existing,
    Vector2D<double>& candidate, 
    int x,
    int y,
    double radius 
) {
    int x_start = x - 2 < 0 ? 0 : x - 2;
    int x_end = x + 2 < grid.size() ? x + 3 : grid.size();
    int y_start = y - 2 < 0 ? 0 : y - 2;
    int y_end = y + 2 < grid[0].size() ? y + 3 : grid[0].size();
    for (int i = x_start; i < x_end; ++i) {
        for (int j = y_start; j < y_end; ++j) {
            if (grid[i][j] == -1) {
                continue;
            }
            double sq_dist = (candidate - existing[grid[i][j]]).VlenSq();
            if (sq_dist < radius * radius) {
                return false;
            }
        }
    }
    return true;
}


std::vector<Vector2D<double>> PoissonDiskSampler::fill(double radius, int k) {
    
    std::uniform_real_distribution<double> domain_dist(0, 1.0);
    std::vector<Vector2D<double>> samples;
    std::vector<Vector2D<double>> active;
    samples.emplace_back(domain_dist(rng), domain_dist(rng));
    active.push_back(samples.back());

    double cell_length = radius * (std::sqrt(2) / 2);
    double inv_cell_length = std::ceil(1.0 / cell_length);
    int grid_size = std::ceil(inv_cell_length);
    std::vector<std::vector<int>> grid(
        grid_size, std::vector<int>(grid_size, -1)
    );
    int idx = active.back().x * inv_cell_length;
    int idy = active.back().y * inv_cell_length;
    grid[idx][idy] = 0;

    //int counter = 0;

    while (!active.empty()) {

        //if (counter++ == 10) { break; }
        
        std::uniform_int_distribution<int> progenitor_selector(
            0, active.size() - 1
        );
        int active_id = progenitor_selector(rng);
        std::vector<Vector2D<double>> candidates = sample_annulus(
            active[active_id], radius, 2 * radius, k
        );

        Vector2D<double>* new_sample = nullptr;
        for (Vector2D<double>& candidate : candidates) {
            idx = candidate.x * inv_cell_length;
            idy = candidate.y * inv_cell_length;
            if (does_not_overlap(grid, samples, candidate, idx, idy, radius)) {
                new_sample = &candidate;
                break;
            }
        }

        if (new_sample != nullptr) {
            grid[idx][idy] = samples.size();
            active.push_back(*new_sample);
            samples.push_back(*new_sample);
        }
        else {
            active[active_id] = active.back();
            active.pop_back();
        }
    }
    std::cout << "Number of samples: " << samples.size() << std::endl;
    return samples;
}


std::vector<Coordinate> PoissonDiskSampler::fill_grid(
    int size, double radius, int k
) {
    double unit_radius = radius / size;
    std::vector<Vector2D<double>> unit_fill = fill(unit_radius, k);

    std::vector<Coordinate> samples;
    samples.reserve(unit_fill.size());
    for (Vector2D<double>& sample : unit_fill) {
        samples.emplace_back(sample.x * size, sample.y * size);
    }

    return samples;
}

