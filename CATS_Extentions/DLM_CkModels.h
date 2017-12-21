
#ifndef DLM_CKMODELS_H
#define DLM_CKMODELS_H

double Flat_Residual(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_SingletTriplet(const double& Momentum, const double* SourcePar, const double* PotPar);
double Lednicky_gauss_Sigma0(const double &Momentum, const double* SourcePar, const double* PotPar);
double pXi_pheno(const double &Momentum, const double* SourcePar, const double* PotPar);

#endif
