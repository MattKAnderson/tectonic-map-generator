#include <chrono>
#include <iostream>
#include <random>


template<typename ValueType>
class GenericDistributionHolder {
public:
    virtual ~GenericDistributionHolder() {};
    virtual ValueType get_sample(std::mt19937_64& g) = 0;
};

template <typename DistributionType, typename ValueType>
class DistributionHolder : public GenericDistributionHolder<ValueType> {
public:
    DistributionHolder(DistributionType distribution) : distribution_(distribution) {};
    ValueType get_sample(std::mt19937_64& g) { return distribution_(g); };

private:
    DistributionType distribution_;
};


template <typename ValueType>
class GenericDistribution {
public:

    template <typename DistributionType>
    GenericDistribution(DistributionType distribution) : 
        distribution(new DistributionHolder<DistributionType, ValueType>(distribution)) {};

    ~GenericDistribution() { delete distribution; };

    ValueType operator()(std::mt19937_64& g) { return distribution->get_sample(g); };

private:
    GenericDistributionHolder<ValueType>* distribution;
};


// test
int main() {

    std::mt19937_64 rng;
    rng.seed(
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now().time_since_epoch()
        ).count()
    );

    GenericDistribution<int> my_distribution(
        std::uniform_int_distribution(1, 6)
    );

    std::cout << "Here are 10 random numbers!\n";
    for (int i = 0; i < 10; i++) {
        std::cout << my_distribution(rng) << "\n";
    }
};