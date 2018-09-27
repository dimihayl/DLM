
#ifndef DLM_INTEGRATION_H
#define DLM_INTEGRATION_H

//note that all functions with "Wiki" in their name have being "stolen" from Wikipedia

double DLM_INT_Trapez(const double& a, const double& b, const unsigned& N);
double DLM_INT_Simpson(const double& a, const double& b, const unsigned& N);
double DLM_INT_TrapezWiki(const double& a, const double& b, const unsigned& N);
double DLM_INT_SimpsonWiki(const double& a, const double& b, const unsigned& N);
double DLM_INT_aSimpsonWiki(const double& a, const double& b, const double& epsilon=1e-6, const int& maxRecursionDepth=128);

void DLM_INT_SetFunction(double (*f)(const double&));
//sometimes one could want to fit function of many parameters, but only integrate over one of them.
//to achieve this one specifies which parameter n is to be integrated over
void DLM_INT_SetFunction(double (*f)(double*), double* par, const unsigned& n);

#endif

