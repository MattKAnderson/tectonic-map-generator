#include <iostream>
#include <vector>
#include <Coordinate.hpp>
#include <polygon_center.hpp>
#include <FileIO.hpp>

/*
 *   This initial test reveals that the initial implementation is working correctly
 *   with an appropriate choice of k parameter. This implementation could be further
 *   improved by:
 *       -  allowing more than 1 grid cell to be traversed in a single step
 *       -  initializing k parameter automatically by testing how much the distance
 *          total can change from the initial point
 *       -  as the result is converged towards, a smaller step size, and a larger k
 *          parameter is used to home in on convergence
 *       -  in the final iters, gather samples and return the most frequent sample as 
 *          the centroid (at the end the simulation may jostle around the centroid, but
 *          it will spend most it's time at the centroid as this is the minimum total
 *          distance point)
 */


int main() {

    int seed = 452653;
    int xmax = 1024;
    int ymax = 1024;
    int iters = 1000;
    double k = 2.5;

    std::vector<Coordinate> points = {{
        {546, 345}, {591, 288}, {521, 333}, {531, 341}, {555, 321},
        {524, 322}, {511, 291}, {544, 301}, {546, 319}, {575, 297}
    }};

    CentroidFinder centroid_finder(seed);

    Coordinate centroid = centroid_finder.toroidal_centroid(
        points, xmax, ymax, iters, k, true
    );

    std::vector<Coordinate> wrapped_centroid = {centroid};
    auto visit_history = centroid_finder.get_visit_history();
    auto cost_history = centroid_finder.get_cost_history();

    std::cout << "The centroid found is: (" << centroid.x << ", " 
              << centroid.y << ")\n";

    vector_output("visit_history.json", visit_history);
    vector_output("cost_history.json", cost_history);
    vector_output("sample_points.json", points);
    vector_output("centroid.json", wrapped_centroid);
}