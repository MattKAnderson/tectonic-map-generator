#include <iostream>
#include <random>

int main() {

    // can play around with this parameter later, but it does seem like sampling the normal
    // distribution for drift change will produce the behaviour I would like to see
    std::mt19937_64 rng(0);
    std::normal_distribution<double> dist(0.0, 0.1);

    for (int i = 0; i < 10; i++) {
        std::cout << dist(rng) << "\n";
    }

}