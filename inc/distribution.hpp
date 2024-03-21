/*
 *   A wrapper for the standard library random number generator distributions. 
 *   Defines an abstract base class for the distributions, so that a general
 *   reference to a distribution can be made. 
 * 
 *   Written By: Matt Anderson 
 *   Last updated: 15/03/2024
 */

#include <concepts>
#include <random>

template<typename ResultType>
requires std::unsigned_integral<ResultType>
class GenericUniformRandomBitGeneratorHolder {
public:
    virtual ~GenericDistributionHolder() {};
    virtual ResultType operator()() = 0;
};

template <typename URBG, typename ResultType>
requires std::uniform_random_bit_generator<URBG> 
         && std::unsigned_integral<ResultType>
class UniformRandomBitGeneratorHolder : 
      GenericUniformRandomBitGeneratorHolder<ResultType> {
public:
    UniformRandomBitGeneratorHolder(URBG generator) : generator_(generator) {};
    ResultType operator()() { return generator_(); }
private:
    URBG generator_;
};

template <typename ResultType>
requires std::unsigned_integral<ResultType>
class UniformRandomBitGenerator {
public:
    template <typename URBG>
    requires std::uniform_random_bit_generator<URBG>
    UniformRandomBitGenerator(URBG urbg) : urbg_(new UniformRandomBitGeneratorHolder<URBG, ResultType>(urbg)) {};
    GenericUniformRandomBitGeneratorHolder<ResultType>* operator()() { return urbg_->(); };
    
private:
    GenericUniformRandomBitGeneratorHolder<ResultType>* urbg_;
};


template<typename ValueType>
class GenericDistributionHolder {
public:
    virtual ~GenericDistributionHolder() {};
    virtual ValueType get_sample(UniformRandomBitGenerator& g) = 0;
};

template <typename DistributionType, typename ValueType>
class DistributionHolder : public GenericDistributionHolder<ValueType> {
public:
    DistributionHolder(DistributionType distribution) : distribution_(distribution) {};
    ValueType get_sample(UniformRandomBitGenerator& g) { return distribution_(g); };

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


/*

template <typename T>
class GenericDistribution {
public:
    virtual ~GenericDistribution() {};

    virtual T operator()() = 0;
};

template <typename IntType = int>
class UniformIntDistribution : public GenericDistribution {
public:
    UniformIntDistribution(IntType a);
    UniformIntDistribution(IntType a, IntType b) {
        dist = std::uniform_int_distribution<IntType>(a, b);
    }

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::uniform_int_distribution<IntType> dist;
};

template<typename RealType = double>
class UniformRealDistribution : public GenericDistribution {
public:
    UniformRealDistribution(RealType a);
    UniformRealDistribution(RealType a, RealType b);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::uniform_real_distribution<RealType> dist;
};

class BernoulliDistribution : public GenericDistribution<bool> {
public:
    BernoulliDistribution(double p);
    
    template <class Generator>
    bool operator()(Generator& g);

private:
    std::bernoulli_distribution dist;
};

template <typename IntType = int>
class BinomialDistribution : public GenericDistribution {
public:
    BinomialDistribution(IntType t, double p);

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::binomial_distribution<IntType> dist;
};

template <typename IntType = int>
class NegativeBinomialDistribution : public GenericDistribution {
public:
    NegativeBinomialDistribution(IntType k, double p);

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::negative_binomial_distribution<IntType> dist;
};

template <typename IntType = int>
class GeometricDistribution : public GenericDistribution {
public:
    GeometricDistribution(double p);

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::geomertric_distribution<IntType> dist;
};

template <typename IntType = int>
class PoissonDistribution : public GenericDistribution {
public:
    PoissonDistribution(double mean);

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::poisson_distribution<IntType> dist;
};

template<typename RealType = double>
class ExponentialDistribution : public GenericDistribution {
public:
    ExponentialDistribution(RealType lambda);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::uniform_real_distribution<RealType> dist;
};

template<typename RealType = double>
class GammaDistribution : public GenericDistribution {
public:
    GammaDistribution(RealType alpha);
    GammaDistribution(RealType alpha, RealType beta);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::gamma_distribution<RealType> dist;
};

template<typename RealType = double>
class WeibullDistribution : public GenericDistribution {
public:
    WeibullDistribution(RealType a);
    WeibullDistribution(RealType a, RealType b);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::weibull_distribution<RealType> dist;
};

template<typename RealType = double>
class ExtremeValueDistribution : public GenericDistribution {
public:
    ExtremeValueDistribution(RealType a);
    ExtremeValueDistribution(RealType a, RealType b);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::extreme_value_distribution<RealType> dist;
};

template<typename RealType = double>
class NormalDistribution : public GenericDistribution {
public:
    NormalDistribution(RealType mean);
    NormalDistribution(RealType mean, RealType stddev);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::normal_distribution<RealType> dist;
};

template<typename RealType = double>
class LogNormalDistribution : public GenericDistribution {
public:
    LogNormalDistribution(RealType m);
    LogNormalDistribution(RealType m, RealType s);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::lognormal_distribution<RealType> dist;
};

template<typename RealType = double>
class ChiSquaredDistribution : public GenericDistribution {
public:
    ChiSquaredDistribution(RealType n);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::chi_squared_distribution<RealType> dist;
};

template<typename RealType = double>
class CauchyDistribution : public GenericDistribution {
public:
    CauchyDistribution(RealType a);
    CauchyDistribution(RealType a, RealType b);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::cauchy_distribution<RealType> dist;
};

template<typename RealType = double>
class FisherFDistribution : public GenericDistribution {
public:
    FisherFDistribution(RealType m);
    FisherFDistribution(RealType m, RealType n);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::fisher_f_distribution<RealType> dist;
};

template<typename RealType = double>
class StudentTDistribution : public GenericDistribution {
public:
    StudentTDistribution(RealType n);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::student_t_distribution<RealType> dist;
};

template<typename IntType = int>
class DiscreteDistribution : public GenericDistribution {
public:
    template <class InputIt>
    DiscreteDistribution(InputIt first, InputIt last);
    DiscreteDistribution(std::initializer_list<double> weights);
    template <class UnaryOperation>
    DiscreteDistribution(std::size_t count, double min, double max, 
                         UnaryOperation unary_op);

    template <class Generator>
    IntType operator()(Generator& g);

private:
    std::discrete_distribution<IntType> dist;
};

template<typename RealType = double>
class PiecewiseConstantDistribution : public GenericDistribution {
public:
    template<class InputIt1, class InputIt2>
    PiecewiseConstantDistribution(InputIt1 first_1, InputIt1 last_1, 
                                  InputIt2 first_w);
    template<class UnaryOperation>
    PiecewiseConstantDistribution(std::initializer_list<RealType> ilist_i,
                                  UnaryOperation fw);
    template<class UnaryOperation>
    PiecewiseConstantDistribution(std::size_t nw, RealType xmin, RealType xmax,
                                  UnaryOperation fw);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::piecewise_constant_distribution<RealType> dist;
};

template<typename RealType = double>
class PiecewiseLinearDistribution : public GenericDistribution {
public:
    template<class InputIt1, class InputIt2>
    PiecewiseLinearDistribution(InputIt1 first_1, InputIt1 last_1, 
                                InputIt2 first_w);
    template<class UnaryOperation>
    PiecewiseLinearDistribution(std::initializer_list<RealType> ilist_i,
                                UnaryOperation fw);
    template<class UnaryOperation>
    PiecewiseLinearDistribution(std::size_t nw, RealType xmin, RealType xmax,
                                UnaryOperation fw);

    template <class Generator>
    RealType operator()(Generator& g);

private:
    std::piecewise_linear_distribution<RealType> dist;
};

template <typename IntType>
UniformIntDistribution<IntType>::UniformIntDistribution(IntType a) {
    dist = std::uniform_int_distribution<IntType>(a);
}

template <typename IntType>
UniformIntDistribution<IntType>::UniformIntDistribution(IntType a, IntType b) {
    dist = std::uniform_int_distribution<T>(a, b);
}

template <typename IntType>
template <class Generator>
IntType UniformIntDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
UniformRealDistribution<RealType>::UniformRealDistribution(RealType a) {
    dist = std::uniform_real_distribution<RealType>(a);
}

template <typename RealType>
UniformRealDistribution<RealType>::UniformRealDistribution(RealType a, 
                                                           RealType b ) {
    dist = std::uniform_real_distribution<RealType>(a, b);
}

template <typename RealType>
template <class Generator>
RealType UniformRealDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

BernoulliDistribution::BernoulliDistribution(double p) {
    dist = std::bernoulli_distribution(p);
}

template <class Generator>
bool BernoulliDistribution::operator()(Generator& g) {
    return dist(g);
}

template <typename IntType>
BinomialDistribution<IntType>::BinomialDistribution(IntType t, double p) {
    dist = std::binomial_distribution<IntType>(t, p);
}

template <typename IntType>
template <class Generator>
IntType BinomialDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename IntType>
NegativeBinomialDistribution<IntType>::NegativeBinomialDistribution(IntType k,
                                                                    double p) {
    dist = std::negative_binomial_distribution<IntType>(k, p);
}

template <typename IntType>
template <class Generator>
IntType NegativeBinomialDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename IntType>
GeometricDistribution<IntType>::GeometricDistribution(double p) {
    dist = std::geometric_distribution<IntType>(p);
}

template <typename IntType>
template <class Generator>
IntType GeometricDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename IntType>
PoissonDistribution<IntType>::PoissonDistribution(double mean) {
    dist = std::poisson_distribution<IntType>(mean);
}

template <typename IntType>
template <class Generator>
IntType PoissonDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
ExponentialDistribution<RealType>::ExponentialDistribution(RealType lambda) {
    dist = std::exponential_distribution<RealType>(lambda);
}

template <typename RealType>
template <class Generator>
RealType ExponentialDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
GammaDistribution<RealType>::GammaDistribution(RealType alpha) {
    dist = std::gamma_distribution<RealType>(alpha);
}

template <typename RealType>
GammaDistribution<RealType>::GammaDistribution(RealType alpha, RealType beta) {
    dist = std::gamma_distribution<RealType>(alpha, beta);
}

template <typename RealType>
template <class Generator>
RealType GammaDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
WeibullDistribution<RealType>::WeibullDistribution(RealType a) {
    dist = std::weibull_distribution<RealType>(a);
}

template <typename RealType>
WeibullDistribution<RealType>::WeibullDistribution(RealType a, RealType b) {
    dist = std::weibull_distribution<RealType>(a, b);
}

template <typename RealType>
template <class Generator>
RealType WeibullDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
ExtremeValueDistribution<RealType>::ExtremeValueDistribution(RealType a) {
    dist = std::extreme_value_distribution<RealType>(a);
}

template <typename RealType>
ExtremeValueDistribution<RealType>::ExtremeValueDistribution(RealType a, 
                                                             RealType b) {
    dist = std::extreme_value_distribution<RealType>(a, b);
}

template <typename RealType>
template <class Generator>
RealType ExtremeValueDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
NormalDistribution<RealType>::NormalDistribution(RealType mean) {
    dist = std::normal_distribution<RealType>(mean);
}

template <typename RealType>
NormalDistribution<RealType>::NormalDistribution(RealType mean, 
                                                 RealType stddev) {
    dist = std::normal_distribution<RealType>(mean, stddev);
}

template <typename RealType>
template <class Generator>
RealType NormalDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
LogNormalDistribution<RealType>::LogNormalDistribution(RealType m) {
    dist = std::log_normal_distribution<RealType>(m);
}

template <typename RealType>
LogNormalDistribution<RealType>::LogNormalDistribution(RealType m, 
                                                       RealType s) {
    dist = std::log_normal_distribution<RealType>(m, s);
}

template <typename RealType>
template <class Generator>
RealType LogNormalDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
ChiSquaredDistribution<RealType>::ChiSquaredDistribution(RealType n) {
    dist = std::chi_squared_distribution<RealType>(n);
}

template <typename RealType>
template <class Generator>
RealType ChiSquaredDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
CauchyDistribution<RealType>::CauchyDistribution(RealType a) {
    dist = std::cauchy_distribution<RealType>(a);
}

template <typename RealType>
CauchyDistribution<RealType>::CauchyDistribution(RealType a, 
                                                             RealType b) {
    dist = std::cauchy_distribution<RealType>(a, b);
}

template <typename RealType>
template <class Generator>
RealType CauchyDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
FisherFDistribution<RealType>::FisherFDistribution(RealType m) {
    dist = std::fisher_f_distribution<RealType>(m);
}

template <typename RealType>
FisherFDistribution<RealType>::FisherFDistribution(RealType m, RealType n) {
    dist = std::fisher_f_distribution<RealType>(m, n);
}

template <typename RealType>
template <class Generator>
RealType FisherFDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
StudentTDistribution<RealType>::StudentTDistribution(RealType n) {
    dist = std::student_t_distribution<RealType>(n);
}

template <typename RealType>
template <class Generator>
RealType StudentTDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename IntType>
template <class InputIt>
DiscreteDistribution<IntType>::DiscreteDistribution(InputIt first, InputIt last) {
    dist = std::discrete_distribution<IntType>(first, last);
}

template <typename IntType>
DiscreteDistribution<IntType>::DiscreteDistribution(
    std::initializer_list<double> weights 
) {
    dist = std::discrete_distribution<IntType>(weights);
}

template <typename IntType>
template <class UnaryOperation>
DiscreteDistribution<IntType>::DiscreteDistribution(std::size_t count, 
                                                     double xmin,
                                                     double xmax,
                                                     UnaryOperation unary_op) {
    dist = std::discrete_distribution<IntType>(count, xmin, xmax, unary_op);
}

template <typename IntType>
template <class Generator>
IntType DiscreteDistribution<IntType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
template <class InputIt1, class InputIt2>
PiecewiseConstantDistribution<RealType>::PiecewiseConstantDistribution(
    InputIt1 first_i, InputIt1 last_i, InputIt2 first_w
) {
    dist = std::piecewise_constant_distribution<RealType>(first_i, last_i, first_w);
}

template <typename RealType>
template <class UnaryOperation>
PiecewiseConstantDistribution<RealType>::PiecewiseConstantDistribution(
    std::initializer_list<RealType> ilist_i, UnaryOperation fw
) {
    dist = std::piecewise_constant_distribution<RealType>(ilist_i, fw);
}

template <typename RealType>
template <class UnaryOperation>
PiecewiseConstantDistribution<RealType>::PiecewiseConstantDistribution(
    std::size_t nw, RealType xmin, RealType xmax, UnaryOperation fw
) {
    dist = std::piecewise_constant_distribution<RealType>(nw, xmin, xmax, fw);
}

template <typename RealType>
template <class Generator>
RealType PiecewiseConstantDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

template <typename RealType>
template <class InputIt1, class InputIt2>
PiecewiseLinearDistribution<RealType>::PiecewiseLinearDistribution(
    InputIt1 first_i, InputIt1 last_i, InputIt2 first_w
) {
    dist = std::piecewise_linear_distribution<RealType>(first_i, last_i, first_w);
}

template <typename RealType>
template <class UnaryOperation>
PiecewiseLinearDistribution<RealType>::PiecewiseLinearDistribution(
    std::initializer_list<RealType> ilist_i, UnaryOperation fw
) {
    dist = std::piecewise_linear_distribution<RealType>(ilist_i, fw);
}

template <typename RealType>
template <class UnaryOperation>
PiecewiseLinearDistribution<RealType>::PiecewiseLinearDistribution(
    std::size_t nw, RealType xmin, RealType xmax, UnaryOperation fw
) {
    dist = std::piecewise_linear_distribution<RealType>(nw, xmin, xmax, fw);
}

template <typename RealType>
template <class Generator>
RealType PiecewiseLinearDistribution<RealType>::operator()(Generator& g) {
    return dist(g);
}

*/