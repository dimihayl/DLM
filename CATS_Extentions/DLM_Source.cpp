#include "DLM_Source.h"

#include "math.h"

const double PI = 3.141592653589793;

double GaussSourceTF1(double* x, double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = *x;
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return Pars[0]*4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

double GaussSource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
//same as GaussSource, but we assume that the user calls the function wanting to sample from theta as well.
//since Integral(dTheta) = 2 for a flat theta distribution and the whole Source function needs to be normalized to 1,
//in order to preserve this we should divide the whole function by 2
double GaussSourceTheta(double* Pars){
    return 0.5*GaussSource(Pars);
}

//double CauchySource(const double& Radius, const double& Size){
double CauchySource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 2.97*2.*Size*Radius*Radius/PI*pow(Radius*Radius+0.25*2.97*2.97*Size*Size,-2.);
}

double CauchySourceTheta(double* Pars){
    return 0.5*CauchySource(Pars);
}
