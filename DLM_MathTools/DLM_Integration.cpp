#include "DLM_Integration.h"
#include <iostream>
#include <stdio.h>

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
        DLM_INT_TEMP_PAR[DLM_INT_TEMP_N] = x;
        return DLM_INT_TEMP_FUN2(DLM_INT_TEMP_PAR);
    }
    else{
        printf("\033[1;31mERROR!\033[0m DLM_Integration: NULL pointer to the function!\n");
    }
}


double DLM_INT_Trapez(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	for(int i=1; i<N; i++){ //interm. steps
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
	for(int i=1; i<N; i++){ //interm. steps
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
	for(int i=0; i<N; i++){ //interm. steps
	    result += h*(FUNCTION(a+i*h) + FUNCTION(a+(i+1)*h))/2.;
	}
	return result;
}

double DLM_INT_SimpsonWiki(const double& a, const double& b, const unsigned& N){
    if(!N) return 0;
	double result = 0;
	double h = (b-a)/double(N);
	for(int i=0; i<N; i++){ //interm. steps
		result += h*(FUNCTION(a+i*h)+4*FUNCTION(a+i*h+0.5*h)+FUNCTION(a+(i+1)*h))/6.;
	}
	//result += h/3.*Function1(a, 5);  //first step
	//result += h/3.*Function1(b-h, 5);//last step
	return result;
}
