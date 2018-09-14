
#ifndef DLM_SOURCE_H
#define DLM_SOURCE_H

double GaussSource(double* Pars);
double GaussSourceTF1(double* x, double* Pars);
double GaussSourceTheta(double* Pars);
double CauchySource(double* Pars);
double CauchySourceTheta(double* Pars);

double DoubleGaussSource(double* Pars);
double GaussCauchySource(double* Pars);

//a monte-carlo out-side-long Gaussian source. Works very slowly!
double GaussOSL_MC(double* Pars);

double Gauss_Exp_Approx(double* Pars);
double Gauss_Exp(double* Pars);
double GaussExpSimple(double* Pars);
double GaussExpTotSimple(double* Pars);
double GaussExpTotIdenticalSimple(double* Pars);
#endif
