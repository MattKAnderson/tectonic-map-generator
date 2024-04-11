#pragma once
#include <cmath>
#include <iostream>
#include <vector>
#include <random>
#include <Coordinate.hpp>
#include <Vector2D.hpp>

// implement the fast algorithm that fills the space with points
// this can be used to initialize a voronoi diagram by first
// identifying the average desired area of cells, then setting 
// the poisson disk radius such that the disk's area is equal to
// or a little less than the desired area


// could combine this with the nested voronoi diagram approach tried 
// earlier today (Apr 10)

// ALSO -- need to retry the nested voronoi diagram approach. I was not
//         computing the midpoints correctly on the torus -- the approach
//         used would have been correct for euclidean space

// ALSO pt.2 -- This means I was not computing the cell centroids correctly
//              when implementing Llyod's algorithm either. Which is probably
//              why it was never converging...


class PoissonDiskSampler {
public:
    PoissonDiskSampler(int seed);
    std::vector<Vector2D<double>> fill(double radius, int k);
    std::vector<Coordinate> fill_grid(
        int size, double radius, int k
    );

private:
    std::mt19937_64 rng;
    std::vector<Vector2D<double>> sample_annulus(
        Vector2D<double> origin, double r1, double r2, double nsamples
    );
};
