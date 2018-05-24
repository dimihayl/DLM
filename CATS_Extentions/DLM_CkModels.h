
#ifndef DLM_CKMODELS_H
#define DLM_CKMODELS_H

double CoulombEta(const double& Momentum, const double& RedMass, const double& Q1Q2);
double CoulombPenetrationFactor(const double& eta);

double LednickyAsInStar(const double& Momentum, const double& GaussR, const double& ScattLenSin, const double& EffRangeSin,
                        const double& Norm, const double& lambda, const double& ares, const double& RadRes);

double GeneralLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const double& ScattLenTri, const double& EffRangeTri,
                       const bool& SinOnly, const bool& QS);

double GeneralCoulombLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const bool& QS, const double& RedMass, const double& Q1Q2);
double GeneralCoulombLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const double& ScattLenTri, const double& EffRangeTri,
                       const bool& QS, const double& RedMass, const double& Q1Q2);

double Flat_Residual(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_Identical_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double LednickyCoulomb_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar);
double LednickyCoulomb_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar);
double LednickyCoulomb_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double LednickyCoulomb_Identical_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_SingletTriplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_gauss_Sigma0(const double &Momentum, const double* SourcePar, const double* PotPar);
double pXi_pheno(const double &Momentum, const double* SourcePar, const double* PotPar);

#endif
