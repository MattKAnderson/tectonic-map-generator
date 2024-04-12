#include <iostream>
#include <vector>
#include <Coordinate.hpp>
#include <polygon_center.hpp>


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
 * 
 *    Also TODO: wrap this into a class, and add data structure to track the visit 
 *               coordinates and total distance, so that visualizations/movies can
 *               be made from how the simulation is working.
 */


int main() {

    int xmax = 1024;
    int ymax = 1024;
    int iters = 1000;
    double k = 0.75;

    std::vector<Coordinate> points = {{
        {546, 345}, {591, 288}, {521, 333}, {531, 341}, {555, 321},
        {524, 322}, {511, 291}, {544, 301}, {546, 319}, {575, 297}
    }};

    Coordinate centroid = find_toroidal_centroid(
        points, xmax, ymax, iters, k
    );

    std::cout << "The centroid found is: (" << centroid.x << ", " 
              << centroid.y << ")\n";
}