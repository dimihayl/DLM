
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
    double Exponential(const double& lambda);
    double Stable(const double& stability=2, const double& location=0, const double& scale=1, const double& skewness=0);
    double GaussR(const unsigned short& dim, const double& mean, const double& sigma);
    double CauchyR(const unsigned short& dim, const double& mean, const double& sigma);
    double StableR(const unsigned short& dim, const double& stability=2, const double& location=0, const double& scale=1, const double& skewness=0);
    double GaussR(const unsigned short& dim, const double* mean, const double* sigma);
    double CauchyR(const unsigned short& dim, const double* mean, const double* sigma);
    double StableR(const unsigned short& dim, const double* stability, const double* location, const double* scale, const double* skewness);
    double StableR(const unsigned short& dim, const double& stability, const double* location, const double* scale, const double* skewness);
    double GaussDiffR(const unsigned short& dim, const double& mean=0, const double& sigma=1);
    double CauchyDiffR(const unsigned short& dim, const double& mean=0, const double& sigma=1);
    double StableDiffR(const unsigned short& dim,const double& stability=2, const double& location=0, const double& scale=1, const double& skewness=0);
    double StableNolan(const unsigned short& dim,const double& stability=2, const double& location=0, const double& scale=1, const double& skewness=0);
    double GaussDiffR(  const unsigned short& dim,
                        const double& mean1, const double& sigma1,
                        const double& mean2, const double& sigma2);
    double CauchyDiffR( const unsigned short& dim,
                        const double& mean1, const double& sigma1,
                        const double& mean2, const double& sigma2);
    double StableDiffR( const unsigned short& dim,
                        const double& stability1, const double& location1, const double& scale1, const double& skewness1,
                        const double& stability2, const double& location2, const double& scale2, const double& skewness2);
    double GaussDiffR(const unsigned short& dim, const double* mean, const double* sigma);
    double CauchyDiffR(const unsigned short& dim, const double* mean, const double* sigma);
    double StableDiffR(const unsigned short& dim,const double* stability, const double* location, const double* scale, const double* skewness);
    double StableDiffR(const unsigned short& dim,const double& stability, const double* location, const double* scale, const double* skewness);

    double StableDiffR( const unsigned short& dim,
                        const double* stability1, const double* location1, const double* scale1, const double* skewness1,
                        const double* stability2, const double* location2, const double* scale2, const double* skewness2);
    double StableDiffR( const unsigned short& dim,
                        const double& stability1, const double* location1, const double* scale1, const double* skewness1,
                        const double& stability2, const double* location2, const double* scale2, const double* skewness2);
private:
    const unsigned SEED;
    std::mt19937_64* MT_RanGen;
    std::uniform_real_distribution<>* RealDist;
    std::exponential_distribution<>* ExpDist;
    std::normal_distribution<>* NormDist;
    std::cauchy_distribution<>* CauchyDist;
    double* STAB_TEMP;
    double* SKEW_TEMP;
    double* SCAL_TEMP;
    double* LOCA_TEMP;
};

#endif

