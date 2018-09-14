
#ifndef DLM_INTEGRATION_H
#define DLM_INTEGRATION_H

double DLM_INT_Trapez(const double& a, const double& b, const unsigned& N);
double DLM_INT_Simpson(const double& a, const double& b, const unsigned& N);
double DLM_INT_TrapezWiki(const double& a, const double& b, const unsigned& N);
double DLM_INT_SimpsonWiki(const double& a, const double& b, const unsigned& N);

void DLM_INT_SetFunction(double (*f)(const double&));
//sometimes one could want to fit function of many parameters, but only integrate over one of them.
//to achieve this one specifies which parameter n is to be integrated over
void DLM_INT_SetFunction(double (*f)(double*), double* par, const unsigned& n);

#endif

