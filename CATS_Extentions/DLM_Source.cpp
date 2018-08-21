#include "DLM_Source.h"
#include "DLM_CRAB_PM.h"

#include "math.h"

const double PI = 3.141592653589793;

//!TEST
//#include <fstream>
//#include <stdio.h>

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
//printf(" G-> r=%.2f, s=%.2f => %.2f\n",Radius,Size,4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size)));
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

double DoubleGaussSource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];

    double& Size1 = Pars[3];
    double& Size2 = Pars[4];
    double& Weight1 = Pars[5];

//printf(" G1=%.1f x %e\n",Weight1, 4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1)));
//printf(" G2=%.1f x %e\n",1.-Weight1, 4.*PI*Radius*Radius*pow(4.*PI*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2)));

    return      Weight1 *4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*4.*PI*Radius*Radius*pow(4.*PI*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2));

}

double GaussCauchySource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];

    double& Size1 = Pars[3];
    double& Size2 = Pars[4];
    double& Weight1 = Pars[5];

//printf(" G1=%.1f x %e\n",Weight1, 4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1)));
//printf(" G2=%.1f x %e\n",1.-Weight1, 4.*PI*Radius*Radius*pow(4.*PI*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2)));

    return      Weight1 *4.*PI*Radius*Radius*pow(4.*PI*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*2.*Size2*Radius*Radius/PI*pow(Radius*Radius+0.25*Size2*Size2,-2.);

}

//a monte-carlo out-side-long Gaussian source. Works very slowly!
//pars[3] = R_OUT
//pars[4] = R_SIDE
//pars[5] = R_LONG
//pars[6] = TAU
//pars[7] = Temperature
double GaussOSL_MC(double* Pars){
    static DLM_CRAB_PM* PM1 = NULL;
    static double Old_rOut;
    static double Old_rSide;
    static double Old_rLong;
    static double Old_Tau;
    static double Old_Temp;
    //bool SourceChanged = false;
    double& Radius = Pars[1];
    double& rOut = Pars[3];
    double& rSide = Pars[4];
    double& rLong = Pars[5];
    double& Tau = Pars[6];
    double& Temp = Pars[7];
    static TH1F* hRO=NULL;
    static TH1F* hRS=NULL;
    static TH1F* hRL=NULL;
    static TH1F* hTau=NULL;
    if(!PM1 || Old_rOut!=rOut || Old_rSide!=rSide || Old_rLong!=rLong || Old_Tau!=Tau || Old_Temp!=Temp){
        if(PM1) delete PM1;
        PM1 = new DLM_CRAB_PM();
        if(hRO) {delete hRO;}
        if(hRS) {delete hRS;}
        if(hRL) {delete hRL;}
        if(hTau) {delete hTau;}
        //SourceChanged = true;
        const unsigned NumBins = 1024;
        const double MeanPR = 0; //SigmaPR = ProtonRS;
        double CurrentX;
        hRO = new TH1F("hRO","hRO",NumBins,MeanPR-5*rOut,MeanPR+5*rOut);
        hRS = new TH1F("hRS","hRS",NumBins,MeanPR-5*rSide,MeanPR+5*rSide);
        hRL = new TH1F("hRL","hRL",NumBins,MeanPR-5*rLong,MeanPR+5*rLong);
        for(unsigned uBin=1; uBin<=NumBins; uBin++){
            CurrentX = hRO->GetBinCenter(uBin);
            hRO->SetBinContent(uBin, 0.5*(1+erf((CurrentX-MeanPR)/(rOut*sqrt(2)))));
            hRS->SetBinContent(uBin, 0.5*(1+erf((CurrentX-MeanPR)/(rSide*sqrt(2)))));
            hRL->SetBinContent(uBin, 0.5*(1+erf((CurrentX-MeanPR)/(rLong*sqrt(2)))));
        }

        const double MeanPT = 0; //SigmaPT=SigmaPR*ProtonTSscale;
        hTau = new TH1F("hTau","hTau",NumBins,MeanPT-5*Tau,MeanPT+5*Tau);
        for(unsigned uBin=1; uBin<=NumBins; uBin++){
            CurrentX = hTau->GetBinCenter(uBin);
            hTau->SetBinContent(uBin, 0.5*(1+erf((CurrentX-MeanPT)/(Tau*sqrt(2)))));
        }

        PM1->SetNumEvents(4000);
        PM1->SetNumPerEvent(50);
        PM1->SetNumSpecies(1);
        PM1->SetParticleMass(0,938);
        PM1->SetParticleSpecies(0,2212);
        PM1->SetTemperature(Temp);

        PM1->SetShapeRx(0,hRO);
        PM1->SetShapeRy(0,hRS);
        PM1->SetShapeRz(0,hRL);
        PM1->SetShapeTau(0,hTau);

        PM1->RunPhasemaker(11);
    }

    Old_rOut = rOut;
    Old_rSide = rSide;
    Old_rLong = rLong;
    Old_Tau = Tau;
    Old_Temp = Temp;
if(rLong>0.5){
//printf("For r=%.2f, NOW=%.2f; OLD=%.2f\n", Radius, PM1->GetIntegratedSource(0)->GetBinContent(PM1->GetIntegratedSource(0)->FindBin(Radius), GaussSource(Pars)));

//printf("For r=%.2f, NOW=%.2f; OLD=%.2f\n", Radius, PM1->GetIntegratedSource(0)->GetBinContent(PM1->GetIntegratedSource(0)->FindBin(Radius)), GaussSource(Pars));
//4.*PI*Radius*Radius*pow(4.*PI*rOut*rOut,-1.5)*exp(-(Radius*Radius)/(4.*rOut*rOut))
//printf(" rOut=%.2f; rSide=%.2f; rLong=%.2f\n\n",rOut,rSide,rLong);
}

    return PM1->GetIntegratedSource(0)->GetBinContent(PM1->GetIntegratedSource(0)->FindBin(Radius));

}
