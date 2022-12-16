#ifndef DLM_CKMODELS_H
#define DLM_CKMODELS_H

#include <complex>

class CATSparameters;

double LednickyAsInStar(const double &Momentum, const double &GaussR, const double &ScattLenSin, const double &EffRangeSin,
                        const double &Norm, const double &lambda, const double &ares, const double &RadRes);

double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const double &ScattLenSin, const double &EffRangeSin,
                       const double &ScattLenTri, const double &EffRangeTri,
                       const bool &SinOnly, const bool &QS, const bool &InverseScatLen = false);
double GeneralLednicky2channel(const double &Momentum, const double &GaussR,
                               const double &ScattLenSin, const double &EffRangeSin,
                               const double &ScattLenTri, const double &EffRangeTri,
                               const bool &QS, const bool &InverseScatLen,
                               const double &Weight1, const double &Weight2);
double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const std::complex<double> &ScattLenSin, const double &EffRangeSin,
                       const std::complex<double> &ScattLenTri, const double &EffRangeTri,
                       const bool &SinOnly, const bool &QS, const bool &InverseScatLen = false);
double GeneralCoulombLednicky(const double &Momentum, const double &GaussR,
                              const double &ScattLenSin, const double &EffRangeSin,
                              const bool &QS, const double &RedMass, const double &Q1Q2);
double GeneralCoulombLednicky(const double &Momentum, const double &GaussR,
                              const double &ScattLenSin, const double &EffRangeSin,
                              const double &ScattLenTri, const double &EffRangeTri,
                              const bool &QS, const double &RedMass, const double &Q1Q2);
double GeneralLednickySill_twochannels(const double &Momentum, const double &GaussR, const double &MassR, const double &Gamma1, const double &Gamma2, const double &m11, const double &m12, const double &m21, const double &m22);

double GeneralLednickySillConstraints_twochannels(const double &Momentum, const double &GaussR, const double &MassR, const double &GammaTilde1, const double &GammaTot, const double &m11, const double &m12, const double &m21, const double &m22);

// Lednicky + Coulomb for complex scattering length and effective Range

double GeneralCoulombLednickyAvg(const double &Momentum, const double &GaussR,
                                 const std::complex<double> &ScattLenSin, const double &EffRangeSin,
                                 const bool &QS, const double &RedMass, const double &Q1Q2);

double Flat_Residual(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Identical_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickyCoulomb_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickyCoulomb_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickyCoulomb_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickyCoulomb_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_SingletTriplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Identical_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar);
double ComplexLednicky_Singlet_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickySilltwochannels_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickySilltwochannelsConstraints_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar);

double LednickySilltwochannelsERE_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar);

double Lednicky_2channel(const double &Momentum, const double *SourcePar, const double *PotPar);

double Lednicky_gauss_Sigma0(const double &Momentum, const double *SourcePar, const double *PotPar);
double pXi_pheno(const double &Momentum, const double *SourcePar, const double *PotPar);
double LednickySingletScatAmplitude(const double &kStar,
                                    const double *SourcePar,
                                    const double *PotPar);

double Lednicky_gauss_pAp_test(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_pAp_CATS(const double &Momentum, const double *SourcePar, const double *PotPar);

double Lednicky_gauss_pAL(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_pAL_varup(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_pAL_varlow(const double &Momentum, const double *SourcePar, const double *PotPar);

// double Lednicky_gauss_pAL(const double &Momentum,  const double* SourcePar,const std::complex<double>& ScattLenSin, const double& EffRangeSin);
double Lednicky_gauss_LAL(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_LAL_varup(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_LAL_varlow(const double &Momentum, const double *SourcePar, const double *PotPar);

// double Lednicky_gauss_LAL(const double &Momentum,  const double* SourcePar,const std::complex<double>& ScattLenSin, const double& EffRangeSin);

double Lednicky_gauss_pAL_v2(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_gauss_LAL_v2(const double &Momentum, const double *SourcePar, const double *PotPar);

// 2-channel Lednicky for pAL
double Lednicky_CC_pAL(const double &Momentum, const double *SourcePar, const double *PotPar);
double Lednicky_CC_LAL(const double &Momentum, const double *SourcePar, const double *PotPar);

// AD-hoc function to read

void SetLedniIntegral_SourceFunction(double (*AS)(double *), CATSparameters &Pars);
void SetLedniIntegral_SourceClass(void *context, const unsigned &numparameters = 0);
void RemoveLedniIntegral_SourceFunction();
void RemoveLedniIntegral_SourceClass();

#endif
