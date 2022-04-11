#ifndef DLM_ROOTFINDER_H
#define DLM_ROOTFINDER_H

#include <iostream>

double NewtonRapson(double (*FUN)(const double&, const double*),const double&  xMin, const double& xMax,
                      const double* pars=NULL, const double& EpsilonX=0, const unsigned& maxIter=1024);



#endif
