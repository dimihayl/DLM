#include "DLM_Integration.h"
#include <iostream>
#include <stdio.h>
#include <math.h>

double (*DLM_INT_TEMP_FUN1)(const double& x);
double (*DLM_INT_TEMP_FUN2)(double* par);
double* DLM_INT_TEMP_PAR;
unsigned DLM_INT_TEMP_N;

void DLM_INT_SetFunction(double (*f)(const double&)){
    DLM_INT_TEMP_FUN1 = f;
    DLM_INT_TEMP_FUN2 = NULL;
}
void DLM_INT_SetFunction(double (*f)(double*), double* par, const unsigned& n){
    DLM_INT_TEMP_FUN1 = NULL;
    DLM_INT_TEMP_FUN2 = f;
    DLM_INT_TEMP_N = n;
    DLM_INT_TEMP_PAR = par;
}
double FUNCTION(const double& x){
    if(DLM_INT_TEMP_FUN1) return DLM_INT_TEMP_FUN1(x);
    else if(DLM_INT_TEMP_FUN2 && DLM_INT_TEMP_PAR){
//printf("x=%f; n=%u\n",x,DLM_INT_TEMP_N);
        DLM_INT_TEMP_PAR[DLM_INT_TEMP_N] = x;
        return DLM_INT_TEMP_FUN2(DLM_INT_TEMP_PAR);
    }
    else{
        printf("\033[1;31mERROR!\033[0m DLM_Integration: NULL pointer to the function!\n");
    }
    return 0;
}


double DLM_INT_Trapez(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	for(unsigned i=1; i<N; i++){ //interm. steps
		result += h*FUNCTION(a+i*h);
	}
	result += 0.5*h*FUNCTION(a);  //first step
	result += 0.5*h*FUNCTION(b);//last step
	return result;
}

double DLM_INT_Simpson(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	for(unsigned i=1; i<N; i++){ //interm. steps
		result += (1+i%2)*2./3.*h*FUNCTION(a+i*h);
	}
	result += h/3.*FUNCTION(a);  //first step
	result += h/3.*FUNCTION(b);//last step
	return result;
}


double DLM_INT_TrapezWiki(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	//double temp;
	for(unsigned i=0; i<N; i++){ //interm. steps
	    result += h*(FUNCTION(a+i*h) + FUNCTION(a+(i+1)*h))/2.;
	}
	return result;
}

double DLM_INT_SimpsonWiki(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	for(unsigned i=0; i<N; i++){ //interm. steps
		result += h*(FUNCTION(a+i*h)+4*FUNCTION(a+i*h+0.5*h)+FUNCTION(a+(i+1)*h))/6.;
	}
	//result += h/3.*Function1(a, 5);  //first step
	//result += h/3.*Function1(b-h, 5);//last step
	return result;
}

//
// Recursive auxiliary function for adaptiveSimpsons() function below
//
double DLM_INT_adaptiveSimpsonsAuxWiki(const double& a, const double& b, const double& epsilon,
                         const double& S, const double& fa, const double&fb, const double& fc, const double& bottom) {
  double c = (a + b)/2, h = b - a;
  double d = (a + c)/2, e = (c + b)/2;
//printf(" FUNCTION(%f)=%f\n",d,FUNCTION(d));
//printf("  FUNCTION(%f)=%f\n",e,FUNCTION(e));
  double fd = FUNCTION(d), fe = FUNCTION(e);
  double Sleft = (h/12)*(fa + 4*fd + fc);
  double Sright = (h/12)*(fc + 4*fe + fb);
  double S2 = Sleft + Sright;

  if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)   // magic 15 comes from error analysis
    return S2 + (S2 - S)/15;
  return DLM_INT_adaptiveSimpsonsAuxWiki(a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1) +
         DLM_INT_adaptiveSimpsonsAuxWiki(c, b, epsilon/2, Sright, fc, fb, fe, bottom-1);
}

//
// Adaptive Simpson's Rule
//
double DLM_INT_aSimpsonWiki(const double& a, const double& b,  // interval [a,b]
                           const double& epsilon,  // error tolerance
                           const int& maxRecursionDepth) {   // recursion cap
  double c = (a + b)/2, h = b - a;
  double fa = FUNCTION(a), fb = FUNCTION(b), fc = FUNCTION(c);
  double S = (h/6)*(fa + 4*fc + fb);
//printf("a=%f; b=%f\n",a,b);
  return DLM_INT_adaptiveSimpsonsAuxWiki(a, b, epsilon, S, fa, fb, fc, maxRecursionDepth);
}


