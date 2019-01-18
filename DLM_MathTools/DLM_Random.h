
#ifndef DLM_RANDOM_H
#define DLM_RANDOM_H

#include <random>

class DLM_Random{
public:
    DLM_Random(const unsigned& seed);
    ~DLM_Random();
    double Uniform(const double& from=0, const double& to=1);
    double Gauss(const double& mean, const double& sigma);
    double Cauchy(const double& mean, const double& sigma);
    double Stable(const double& stability, const double& skewness, const double& scale, const double& location);
    double StableR(const double& stability, const double& skewness, const double& scale, const double& location, const unsigned short& dim);
    double StableR(const double* stability, const double* skewness, const double* scale, const double* location, const unsigned short& dim);
private:
    const unsigned SEED;
    std::mt19937_64* MT_RanGen;
    std::uniform_real_distribution<>* RealDist;
    std::exponential_distribution<>* ExpDist;
    std::normal_distribution<>* NormDist;
    std::cauchy_distribution<>* CauchyDist;
};

#endif

