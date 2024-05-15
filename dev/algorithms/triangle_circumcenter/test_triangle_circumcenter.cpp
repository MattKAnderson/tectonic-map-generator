#include <chrono>
#include <iostream>
#include <random>
#include <vector>
#include <Coordinate.hpp>
#include <triangle_circumcenter.hpp>

std::ostream& operator<<(std::ostream& s, const RealCoordinate& c) {
    s << "(" << c.x << ", " << c.y << ")";
    return s;
}

template <typename Func>
std::vector<double> time_solution(
    Func f, std::mt19937_64& rng, std::uniform_real_distribution<double>& dist, 
    int nsamples
);

int main() {

    unsigned int seed = 9123;    
    unsigned int test_samples = 1000000;

    std::mt19937_64 rng(seed);
    std::uniform_real_distribution dist(-1.0, 1.0);

    RealCoordinate vertex_a(dist(rng), dist(rng));
    RealCoordinate vertex_b(dist(rng), dist(rng));
    RealCoordinate vertex_c(dist(rng), dist(rng));

    RealCoordinate las_result = lin_alg_solution(vertex_a, vertex_b, vertex_c);
    RealCoordinate pbs_result = perp_bisector_solution(vertex_a, vertex_b, vertex_c);

    std::cout << "las result: " << las_result << "\n";
    std::cout << "pbs result: " << pbs_result << "\n";

    std::vector<double> las_times;
    std::vector<double> pbs_times;
    las_times = time_solution(lin_alg_solution, rng, dist, test_samples);
    pbs_times = time_solution(perp_bisector_solution, rng, dist, test_samples);

    double las_average = std::accumulate(las_times.begin(), las_times.end(), 0.0) / las_times.size();
    double pbs_average = std::accumulate(pbs_times.begin(), pbs_times.end(), 0.0) / pbs_times.size();

    std::cout << "las average time: " << las_average << " ns\n";
    std::cout << "pbs average time: " << pbs_average << " ns\n";
    std::cout << "las / pbs ratio:  " << las_average / pbs_average << "\n";
    std::cout << "pbs / las ratio:  " << pbs_average / las_average << "\n";
}

template <typename Func>
std::vector<double> time_solution(
    Func f, std::mt19937_64& rng, std::uniform_real_distribution<double>& dist, 
    int nsamples
) {
    std::vector<double> timings;
    for (int i = 0; i < nsamples; i++) {
        RealCoordinate a(dist(rng), dist(rng));
        RealCoordinate b(dist(rng), dist(rng));
        RealCoordinate c(dist(rng), dist(rng));
        auto start = std::chrono::high_resolution_clock::now();
        RealCoordinate result = f(a, b, c);
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        timings.push_back(elapsed);
    }
    return timings;
}