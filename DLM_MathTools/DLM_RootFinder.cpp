
#include "DLM_RootFinder.h"

#include <math.h>
#include <stdio.h>

double NewtonRapson(double (*FUN)(const double&, const double*),const double&  xMin, const double& xMax,
                    const double* pars, const double& EpsilonX, const unsigned& maxIter){
  double DeltaX;
  double xVal=(xMax+xMin)*0.5;
  double fVal;
  double DeltaF;
  double EpsX = fabs(EpsilonX);
  double xmax,xmin;
  if(xMin<xMax){
    xmin = xMin;
    xmax = xMax;
  }
  else{
    xmin = xMax;
    xmax = xMin;
  }
  unsigned maxIt = maxIter;
  if(maxIt<8) maxIt = 8;
  //we have quadratic convergance, so we would expect the precision to go a maxIter^2
  //to have some buffer for fluctuations, we only take ^1.75 in this estimator
  if(EpsX==0){
    EpsX = (xmax-xmin)*pow(double(maxIt),-1.75);
    if(EpsX<=0) EpsX = 1e-307;
  }

  for(unsigned iIter=0; iIter<maxIt; iIter++){
//printf("Call (0) : %f\n",xVal);
      fVal = (FUN)(xVal,pars);
      if(xVal==xVal+EpsX){EpsX=xVal*(1e-15);}
//printf("Call (1) : %f\n",xVal+EpsX);
      DeltaF = (FUN)(xVal+EpsX,pars) - fVal;
      if(!DeltaF){
          DeltaX=1;
          {printf("\033[1;33mWARNING (DLM_RootFinder):\033[0m DeltaX is fishy with the NewtonRapson solver! Please contact the developers!\n");}
      }
      else DeltaX = -fVal*EpsX/DeltaF;
      xVal += DeltaX;
      int counter = 0;
      while(xVal<xmin || xVal>xmax){
        counter++;
        xVal -= DeltaX*pow(2., -counter);
      }
      //counter=0;
      double fValNew = (FUN)(xVal,pars);
      while(fabs(fValNew)>fabs(fVal) && counter<16){
          counter++;
          xVal -= DeltaX*pow(2., -counter);
          if(xVal<xmin || xVal>xmax){
            {printf("\033[1;33mWARNING (DLM_RootFinder):\033[0m xVal is fishy with the NewtonRapson solver! Please contact the developers!\n");}
          }        
          fValNew = (FUN)(xVal,pars);
      }
      if(counter==16 && fabs(fValNew)>fabs(fVal)){
        {printf("\033[1;33mWARNING (DLM_RootFinder):\033[0m The backtracking of the NewtonRapson root-finder failed!\n");}
      }
      if(fabs(fValNew)<fabs(DeltaF)){
          return xVal;
      }
  }
  {printf("\033[1;33mWARNING (DLM_RootFinder):\033[0m The NewtonRapson root-finder failed!\n");}
  return xVal;
}
