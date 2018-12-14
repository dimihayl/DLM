#include "DLM_Source.h"
#include "DLM_CRAB_PM.h"
#include "DLM_Integration.h"

#include "math.h"

const double PI = 3.141592653589793;

//!TEST
//#include <fstream>
//#include <stdio.h>

double GaussSourceTF1(double* x, double* Pars){
    Pars[1] = *x;
    return GaussSource(Pars);
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

//the Ansatz is: single particle are emitted according to a Gaussian source with R1 and R2
//after they are emitted, they decay following an exponential law. The momentum is assumed to be unchanged
//the convolution can be analytically derived with, e.g. Mathematica, but the solution is literally two pages long
//and unsurprisingly numerically unstable. I have found a workaround in which I just start with an ordinary Gaussian,
//I shift it by some amount (t*k/m), this shift is dynamically modified to grow with r (starting from zero) by introducing some atan dependence,
//and finally the power law of r is changed depending on the radius and tau, we start off with r^2 and transition to ln^2(r+1)
//N.B. This function will most likely fail completely if t*k/m is comparable to the size of the system!
double Gauss_Exp_Approx(double* Pars){
    double& MOM = Pars[0];
    double& RAD = Pars[1];
    double& SIG1 = Pars[3];//size
    double& SIG2 = Pars[4];//size
    double SIG = sqrt(0.5*(SIG1*SIG1+SIG2*SIG2));
    double& TAU1 = Pars[5];
    double& MASS1 = Pars[6];
    double& TAU2 = Pars[7];
    double& MASS2 = Pars[8];
    double MASS = 0.5*(MASS1+MASS2);
    double TKM1 = MASS1?TAU1*MOM/MASS1:0;
    double TKM2 = MASS2?TAU2*MOM/MASS2:0;
    double TKM = TKM1+TKM2;
    //double rMOM1 = sqrt(pow(MMass,4.)-2.*pow(MMass*Mass[0],2.)+pow(Mass[0],4.)-2.*pow(MMass*Mass[1],2.)-2.*pow(Mass[0]*Mass[1],2.)+pow(Mass[1],4.))/(2.*Mass[0]);
    double OldRad = RAD;
    static double oldMOM = 0;
    static double oldSIG = 0;
    static double oldTKM = 0;
    static double oldMASS = 0;
    static double NORM = 1;
    //this function is not normalized by default. This is done here. However we save information about the last time a normalization
    //was performed, if nothing has changed, this means NORM already has the correct value and we skip this rather expensive step
    if( (oldMOM!=MOM || oldSIG!=SIG || oldTKM!=TKM || oldMASS!=MASS) && NORM!=1.23456789 ){
        //this value is just a flag, which basically tells Gauss_Exp_Approx that we are currently computing NORM, hence
        //we are interested in the true value of the integral, without normalization. Thus later on the Result is not normalized
        NORM = 1.23456789;
        DLM_INT_SetFunction(Gauss_Exp_Approx,Pars,1);
        NORM = DLM_INT_aSimpsonWiki(0,SIG*8.+TKM*8.);
        RAD = OldRad;
        oldMOM = MOM;
        oldSIG = SIG;
        oldTKM = TKM;
        oldMASS = MASS;
    }
    RAD = TKM?log(RAD*TKM/80.+1.)*80./TKM-TKM*atan(1.5*RAD/TKM)*2./PI:RAD;
    if(RAD<0) RAD=0;
    if(MASS<=0) {printf("\033[1;33mWARNING:\033[0m Gauss_Exp_Approx got MASS=0\n"); return 0;}
    double Result = NORM!=1.23456789?GaussSource(Pars)/NORM:GaussSource(Pars);
    RAD = OldRad;
    return Result;
}

//the exact solution of the same thing obtained with Mathematica. Stable only if Sigma and Tau are very similar (up to factor 2)
double Gauss_Exp_Exact(double* Pars){
    double& MOM = Pars[0];
    double& RAD = Pars[1];
    double& SIG1 = Pars[3];//size
    double& SIG2 = Pars[4];//size
    double SIG = sqrt(0.5*(SIG1*SIG1+SIG2*SIG2));
    double& TAU1 = Pars[5];
    double& MASS1 = Pars[6];
    double& TAU2 = Pars[7];
    double& MASS2 = Pars[8];
    double TKM1 = MASS1?TAU1*MOM/MASS1:0;
    double TKM2 = MASS2?TAU2*MOM/MASS2:0;
    const double E = exp(1);
    if(TKM1==0 && TKM2==0){
        return GaussSource(Pars);
    }
    else if(TKM1==0){
        return (pow(E,-pow(RAD,2)/(4.*pow(SIG,2)) - RAD/TKM2)*(4*pow(E,RAD/TKM2)*SIG*TKM2*(pow(SIG,2) + pow(TKM2,2)) -
               2*pow(E,pow(RAD,2)/(4.*pow(SIG,2)))*SIG*TKM2*(2*pow(SIG,2) + TKM2*(-RAD + 2*TKM2)) +
               pow(E,pow(RAD,2)/(4.*pow(SIG,2)) + pow(SIG,2)/pow(TKM2,2))*sqrt(PI)*(-4*pow(SIG,4) + 2*pow(SIG,2)*(RAD - 3*TKM2)*TKM2 + RAD*pow(TKM2,3))*
                (erf(RAD/(2.*SIG) - SIG/TKM2) + erf(SIG/TKM2))))/(sqrt(PI)*pow(TKM2,5));
    }
    else if(TKM2==0){
        return (pow(E,-pow(RAD,2)/(4.*pow(SIG,2)) - RAD/TKM1)*(4*pow(E,RAD/TKM1)*SIG*TKM1*(pow(SIG,2) + pow(TKM1,2)) -
               2*pow(E,pow(RAD,2)/(4.*pow(SIG,2)))*SIG*TKM1*(2*pow(SIG,2) + TKM1*(-RAD + 2*TKM1)) +
               pow(E,pow(RAD,2)/(4.*pow(SIG,2)) + pow(SIG,2)/pow(TKM1,2))*sqrt(PI)*(-4*pow(SIG,4) + 2*pow(SIG,2)*(RAD - 3*TKM1)*TKM1 + RAD*pow(TKM1,3))*
                (erf(RAD/(2.*SIG) - SIG/TKM1) + erf(SIG/TKM1))))/(sqrt(PI)*pow(TKM1,5));
    }
//!shouldn't I take TKM=TKM1+TKM2 for the computation, verify!
    else if(TKM1==TKM2){
static bool WARNING = true;
if(WARNING) printf("WARNING: shouldn't I take TKM=TKM1+TKM2 for the computation, verify!\n");
WARNING=false;
        return (pow(E,(pow(SIG,2) - RAD*TKM1)/pow(TKM1,2))*((-2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM1,2))*TKM1) + (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM1,2)))*pow(SIG,3))/TKM1 -
               (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2))))*pow(SIG,3))/TKM1 + sqrt(PI)*pow(SIG,2)*erf(SIG/TKM1) +
               (2*sqrt(PI)*pow(SIG,4)*erf(SIG/TKM1))/pow(TKM1,2) - (sqrt(PI)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM1)*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1)))/
                (pow(TKM1,3)*abs(-RAD/(2.*SIG) + SIG/TKM1)) - (pow(2*pow(SIG,2) - RAD*TKM1,3)*
                  (-abs(-2*pow(SIG,2) + RAD*TKM1) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*sqrt(PI)*SIG*TKM1*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1))))/
                (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*pow(SIG,2)*pow(TKM1,4)*pow(abs(-RAD/(2.*SIG) + SIG/TKM1),3))))/(sqrt(PI)*pow(SIG,2)*TKM1);
    }
    else{
        return (pow(E,(pow(SIG,2) - RAD*TKM1)/pow(TKM1,2))*((-2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM1,2))*TKM1) + (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM1,2)))*pow(SIG,3))/TKM1 -
        (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2))))*pow(SIG,3))/TKM1 + sqrt(PI)*pow(SIG,2)*erf(SIG/TKM1) +
        (2*sqrt(PI)*pow(SIG,4)*erf(SIG/TKM1))/pow(TKM1,2) - (sqrt(PI)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM1)*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1)))/
         (pow(TKM1,3)*abs(RAD/(2.*SIG) - SIG/TKM1)) - (pow(2*pow(SIG,2) - RAD*TKM1,3)*
           (-abs(-2*pow(SIG,2) + RAD*TKM1) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*sqrt(PI)*SIG*TKM1*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1))))/
         (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*pow(SIG,2)*pow(TKM1,4)*pow(abs(RAD/(2.*SIG) - SIG/TKM1),3))) +
     pow(E,(pow(SIG,2) - RAD*TKM2)/pow(TKM2,2))*((2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM2,2))*TKM2) - (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM2,2)))*pow(SIG,3))/TKM2 +
        (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2))))*pow(SIG,3))/TKM2 - sqrt(PI)*pow(SIG,2)*erf(SIG/TKM2) -
        (2*sqrt(PI)*pow(SIG,4)*erf(SIG/TKM2))/pow(TKM2,2) + (sqrt(PI)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM2)*erf(abs(-2*pow(SIG,2) + RAD*TKM2)/(2.*SIG*TKM2)))/
         (pow(TKM2,3)*abs(RAD/(2.*SIG) - SIG/TKM2)) + (pow(2*pow(SIG,2) - RAD*TKM2,3)*
           (-abs(-2*pow(SIG,2) + RAD*TKM2) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2)))*sqrt(PI)*SIG*TKM2*erf(abs(-2*pow(SIG,2) + RAD*TKM2)/(2.*SIG*TKM2))))/
         (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2)))*pow(SIG,2)*pow(TKM2,4)*pow(abs(RAD/(2.*SIG) - SIG/TKM2),3))))/(sqrt(PI)*pow(SIG,2)*(TKM1 - TKM2));
    }
}

double Gauss_Exp(double* Pars){
    double& MOM = Pars[0];
    //double& RAD = Pars[1];
    double& SIG1 = Pars[3];//size
    double& SIG2 = Pars[4];//size
    double SIG = sqrt(0.5*(SIG1*SIG1+SIG2*SIG2));
    double& TAU1 = Pars[5];
    double& MASS1 = Pars[6];
    double& TAU2 = Pars[7];
    double& MASS2 = Pars[8];
    double TKM1 = MASS1?TAU1*MOM/MASS1:0;
    double TKM2 = MASS2?TAU2*MOM/MASS2:0;
    if(SIG/TKM1>3 || SIG/TKM2>3){
        return Gauss_Exp_Approx(Pars);
    }
    else{
        return Gauss_Exp_Exact(Pars);
    }
}


double GaussExpSimple_Approx(double* Pars){
    double& MOM = Pars[0];
    double& RAD = Pars[1];
    double& SIG = Pars[3];//size
    double& TKM = Pars[4];
    //double rMOM1 = sqrt(pow(MMass,4.)-2.*pow(MMass*Mass[0],2.)+pow(Mass[0],4.)-2.*pow(MMass*Mass[1],2.)-2.*pow(Mass[0]*Mass[1],2.)+pow(Mass[1],4.))/(2.*Mass[0]);
    double OldRad = RAD;
    static double oldMOM = 0;
    static double oldSIG = 0;
    static double oldTKM = 0;
    static double NORM = 1;
    //this function is not normalized by default. This is done here. However we save information about the last time a normalization
    //was performed, if nothing has changed, this means NORM already has the correct value and we skip this rather expensive step
    if( (oldMOM!=MOM || oldSIG!=SIG || oldTKM!=TKM ) && NORM!=1.23456789 ){
        //this value is just a flag, which basically tells Gauss_Exp_Approx that we are currently computing NORM, hence
        //we are interested in the true value of the integral, without normalization. Thus later on the Result is not normalized
        NORM = 1.23456789;
        DLM_INT_SetFunction(GaussExpSimple_Approx,Pars,1);
        NORM = DLM_INT_aSimpsonWiki(0,SIG*8.+TKM*8.);
        RAD = OldRad;
        oldMOM = MOM;
        oldSIG = SIG;
        oldTKM = TKM;
//printf("NORMALIZING! NORM=%.2e\n",NORM);
    }
    RAD = TKM?log(RAD*TKM/80.+1.)*80./TKM-TKM*atan(1.5*RAD/TKM)*2./PI:RAD;
    if(RAD<0) RAD=0;
    double Result = NORM!=1.23456789?GaussSource(Pars)/NORM:GaussSource(Pars);
    RAD = OldRad;
//printf(" RES=%.2e\n",Result);
    return Result;
}
double GaussExpSimple_Exact(double* Pars){
    //double& MOM = Pars[0];
    double& RAD = Pars[1];
    double& SIG = Pars[3];//size
    double& TKM = Pars[4];
    const double E = exp(1);
    if(TKM==0){
        return GaussSource(Pars);
    }
    else{
        return (pow(E,(pow(SIG,2) - RAD*TKM)/pow(TKM,2))*((-2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM,2))*TKM) + (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM,2)))*pow(SIG,3))/TKM -
               (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2))))*pow(SIG,3))/TKM + sqrt(PI)*pow(SIG,2)*erf(SIG/TKM) +
               (2*sqrt(PI)*pow(SIG,4)*erf(SIG/TKM))/pow(TKM,2) - (sqrt(PI)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM)*erf(abs(-2*pow(SIG,2) + RAD*TKM)/(2.*SIG*TKM)))/
                (pow(TKM,3)*abs(-RAD/(2.*SIG) + SIG/TKM)) - (pow(2*pow(SIG,2) - RAD*TKM,3)*
                  (-abs(-2*pow(SIG,2) + RAD*TKM) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2)))*sqrt(PI)*SIG*TKM*erf(abs(-2*pow(SIG,2) + RAD*TKM)/(2.*SIG*TKM))))/
                (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2)))*pow(SIG,2)*pow(TKM,4)*pow(abs(-RAD/(2.*SIG) + SIG/TKM),3))))/(sqrt(PI)*pow(SIG,2)*TKM);
    }
}
double GaussExpSimple(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    double& SIG = Pars[3];//size
    double& TKM = Pars[4];
    if(TKM==0){
        return GaussSource(Pars);
    }
    else if(SIG/TKM>3){
        return GaussExpSimple_Approx(Pars);
    }
    else{
        return GaussExpSimple_Exact(Pars);
    }
}
//the particles can have different TKM and different weights
//total of 8 parameters
double GaussExpTotSimple(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    //double& SIG = Pars[3];//size
    double& TKMA = Pars[4];//t*p/m of particle A
    double& TKMB = Pars[5];//t*p/m of particle B
    //double TKM=TKMA+TKMB;
    const double oldTKMA = TKMA;
    const double oldTKMB = TKMB;
    //fraction of primaries for A, should be between 0 and 1
    if(Pars[6]<0) Pars[6]=0; if(Pars[6]>1) Pars[6]=1;
    const double& primA = Pars[6];
    //fraction of primaries for B, should be between 0 and 1
    if(Pars[6]<0) Pars[7]=0; if(Pars[7]>1) Pars[7]=1;
    const double& primB = Pars[7];
    double Result=0;
    TKMA=0;                 Result += primA*primB*GaussExpSimple(Pars);//both primary
    TKMA=oldTKMB;           Result += primA*(1.-primB)*GaussExpSimple(Pars);//B comes from resonance
    TKMA=oldTKMA;           Result += (1.-primA)*primB*GaussExpSimple(Pars);//A comes from resonance
    TKMA=oldTKMA+oldTKMB;   Result += (1.-primA)*(1.-primB)*GaussExpSimple(Pars);//both come from resonance

    TKMA = oldTKMA;
    TKMB = oldTKMB;

    return Result;
}
//same as GaussExpTotSimple but for identical particles
//total of 6 parameters
double GaussExpTotIdenticalSimple(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    //double& SIG = Pars[3];//size
    double& TKM = Pars[4];
    //double TKM=TKMA+TKMB;
    const double oldTKM = TKM;
    //fraction of primaries
    const double& prim = Pars[5];
    double Result=0;
    TKM=0;          Result += prim*prim*GaussExpSimple(Pars);//both primary
//printf("Result = %.2e; GES=%.2e*%.2e\n",Result,prim*prim,GaussExpSimple(Pars));
    TKM=oldTKM;     Result += 2*prim*(1.-prim)*GaussExpSimple(Pars);//1 comes from resonance
//printf(" Result = %.2e; GES=%.2e*%.2e\n",Result,2*prim*(1.-prim),GaussExpSimple(Pars));
    TKM=2*oldTKM; Result += (1.-prim)*(1.-prim)*GaussExpSimple(Pars);//both come from resonance
//printf("  Result = %.2e; GES=%.2e*%.2e\n",Result,(1.-prim)*(1.-prim),GaussExpSimple(Pars));
//printf("  Result = %.2e; GES=%.2e*%.2e\n",Result,(1.-prim)*(1.-prim),GaussExpSimple(Pars));
//printf("  Result = %.2e; GES=%.2e*%.2e\n",Result,(1.-prim)*(1.-prim),GaussExpSimple(Pars));
//printf("  Should: %.2e\n",(1.-prim)*(1.-prim)*GaussExpSimple(Pars));
    TKM = oldTKM;

    return Result;
}

//introducing the momentum dependence. Pars[4] is now just t/m
//the momentum is taken assuming two body decay of a resonance to primary+pion
double GaussExpTotIdenticalSimple_2body(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    //double& SIG = Pars[3];//size
    double& rTAU = Pars[4];
    //5 = weight of the primaries
    double& rMASS = Pars[6];
    double& p1MASS = Pars[7];
    double& p2MASS = Pars[8];
    double oldrTAU = rTAU;
    rTAU = rTAU*sqrt(pow(rMASS,4.)-2.*pow(rMASS*p1MASS,2.)+pow(p1MASS,4.)-2.*pow(rMASS*p2MASS,2.)-2.*pow(p1MASS*p2MASS,2.)+pow(p2MASS,4.))/(2.*p1MASS)/rMASS;
    double Result = GaussExpTotIdenticalSimple(Pars);
    rTAU = oldrTAU;
    return Result;
}

//introducing the momentum dependence. Pars[4] is now just t/m
//the momentum is taken assuming two body decay of a resonance to primary+pion
double GaussExpTotSimple_2body(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    //double& SIG = Pars[3];//size
    double& rTAU = Pars[4];
    double& prim = Pars[5];//fraction of primaries
    double& rMASS = Pars[6];
    double& p1MASS = Pars[7];
    double& p2MASS = Pars[8];
    double& rTAU_2 = Pars[9];
    double& prim_2 = Pars[10];
    double& rMASS_2 = Pars[11];
    double& p1MASS_2 = Pars[12];
    double& p2MASS_2 = Pars[13];
    double PARS[8];
    double& TKMA = PARS[4];
    double& TKMB = PARS[5];
    TKMA = rTAU*sqrt(pow(rMASS,4.)-2.*pow(rMASS*p1MASS,2.)+pow(p1MASS,4.)-2.*pow(rMASS*p2MASS,2.)-2.*pow(p1MASS*p2MASS,2.)+pow(p2MASS,4.))/(2.*p1MASS)/rMASS;
    TKMB = rTAU_2*sqrt(pow(rMASS_2,4.)-2.*pow(rMASS_2*p1MASS_2,2.)+pow(p1MASS_2,4.)-2.*pow(rMASS_2*p2MASS_2,2.)-2.*pow(p1MASS_2*p2MASS_2,2.)+pow(p2MASS_2,4.))/(2.*p1MASS_2)/rMASS_2;
    PARS[0] = Pars[0]; PARS[1] = Pars[1]; PARS[2] = Pars[2]; PARS[3] = Pars[3];
    PARS[4] = TKMA; PARS[5] = TKMB; PARS[6] = prim; PARS[7] = prim_2;
    return GaussExpTotSimple(PARS);
}

//Gauss including kT dependence, convoluted further with an exponential. The latter can be used to model resonances. In case we have
//two resonances, a nice approximation is just to add their tau/mass
//the parameters are the following:
//[3] = Number of kT bins NkT
//[4] = Number of resonances related to particle 1
//[5] = Number of resonances related to particle 2

//[6 : 5+NkT] = the radii of the different kT bins
//[6+NkT : 3+2*NkT] = the weights of each kT bin
//
//double GaussExpTotSimple_2body_kT(double* Pars){

//}

double MemberSourceForwarder(void* context, double* Pars){return static_cast<MemberSource*>(context)->Eval(Pars);}

double MemberSource::Eval(double* Pars){
    return 0;
}
double MemberSource::Eval(const double& Momentum, const double Radius, const double& Angle){
        PARS[0] = Momentum;
        PARS[1] = Radius;
        PARS[2] = Angle;
        return Eval(PARS);
}

double MS_Gauss::Eval(double* Pars){
    double& Radius = Pars[1];
    return 4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

MS_GaussExp_mT_Simple::MS_GaussExp_mT_Simple(){
    Num_mT = 0;
    Mean_mT = NULL;
    Weight_mT = NULL;
    Linear_mT = 0;
    Slope_mT = 0;
    Mass = new double [2];
    MassR = new double [2];
    MassD = new double [2];
    Tau = new double [2];
    Weight_R = new double [2];
    Parameters = new double [14];
}
MS_GaussExp_mT_Simple::~MS_GaussExp_mT_Simple(){
    delete[]Mass; Mass=NULL;
    delete[]MassR; MassR=NULL;
    delete[]MassD; MassD=NULL;
    delete[]Tau; Tau=NULL;
    delete[]Weight_R; Weight_R=NULL;
    delete[]Parameters; Parameters=NULL;
}

void MS_GaussExp_mT_Simple::SetNum_mT(const unsigned& nmt){
    if(!nmt){
        printf("\033[1;33mWARNING:\033[0m You must have at least one mT bin\n");
        return;
    }
    if(Num_mT==nmt) return;
    Num_mT=nmt;

    if(Mean_mT)delete[]Mean_mT;Mean_mT=new double [nmt];
    if(Weight_mT)delete[]Weight_mT;Weight_mT=new double [nmt];
}
void MS_GaussExp_mT_Simple::SetMean_mT(const unsigned& umt, const double& mmt){
    if(umt>=Num_mT){
        printf("\033[1;33mWARNING:\033[0m Current number of mT bins is %u and you attempt to set bin nr. %u\n",Num_mT,umt);
        return;
    }
    Mean_mT[umt] = mmt;
}
void MS_GaussExp_mT_Simple::SetWeight_mT(const unsigned& umt, const double& wmt){
    if(umt>=Num_mT){
        printf("\033[1;33mWARNING:\033[0m Current number of mT bins is %u and you attempt to set bin nr. %u\n",Num_mT,umt);
        return;
    }
    Weight_mT[umt] = wmt;
}
void MS_GaussExp_mT_Simple::SetLinear_mT(const double& lin){
    Linear_mT = lin;
}
void MS_GaussExp_mT_Simple::SetSlope_mT(const double& slope){
    Slope_mT = slope;
}
void MS_GaussExp_mT_Simple::SetMass(const unsigned short& particle, const double& mass){
    if(particle>1){
        printf("\033[1;33mWARNING:\033[0m MS_GaussExp_mT_Simple can be used only for particle 0 or 1\n");
        return;
    }
    if(mass<0){
        printf("\033[1;33mWARNING:\033[0m You are setting a negative mass!\n");
        return;
    }
    Mass[particle] = mass;
}
void MS_GaussExp_mT_Simple::SetMassR(const unsigned short& particle, const double& mass){
    if(particle>1){
        printf("\033[1;33mWARNING:\033[0m MS_GaussExp_mT_Simple can be used only for particle 0 or 1\n");
        return;
    }
    if(mass<0){
        printf("\033[1;33mWARNING:\033[0m You are setting a negative mass!\n");
        return;
    }
    MassR[particle] = mass;
}
void MS_GaussExp_mT_Simple::SetMassD(const unsigned short& particle, const double& mass){
    if(particle>1){
        printf("\033[1;33mWARNING:\033[0m MS_GaussExp_mT_Simple can be used only for particle 0 or 1\n");
        return;
    }
    if(mass<0){
        printf("\033[1;33mWARNING:\033[0m You are setting a negative mass!\n");
        return;
    }
    MassD[particle] = mass;
}
void MS_GaussExp_mT_Simple::SetTau(const unsigned short& particle, const double& tau){
    if(particle>1){
        printf("\033[1;33mWARNING:\033[0m MS_GaussExp_mT_Simple can be used only for particle 0 or 1\n");
        return;
    }
    if(tau<0){
        printf("\033[1;33mWARNING:\033[0m You are setting a negative lifetime!\n");
        return;
    }
    Tau[particle] = tau;
}
void MS_GaussExp_mT_Simple::SetResonanceWeight(const unsigned short& particle, const double& weight){
    if(particle>1){
        printf("\033[1;33mWARNING:\033[0m MS_GaussExp_mT_Simple can be used only for particle 0 or 1\n");
        return;
    }
    if(weight<0){
        printf("\033[1;33mWARNING:\033[0m You are setting a negative weight!\n");
        return;
    }
    if(weight>1){
        printf("\033[1;33mWARNING:\033[0m You are setting a bigger than 1!\n");
        return;
    }
    Weight_R[particle] = weight;
}

/*
double GaussExpTotSimple_2body(double* Pars){
    //double& MOM = Pars[0];
    //double& RAD = Pars[1];
    //double& SIG = Pars[3];//size
    double& rTAU = Pars[4];
    double& prim = Pars[5];//fraction of primaries
    double& rMASS = Pars[6];
    double& p1MASS = Pars[7];
    double& p2MASS = Pars[8];
    double& rTAU_2 = Pars[9];
    double& prim_2 = Pars[10];
    double& rMASS_2 = Pars[11];
    double& p1MASS_2 = Pars[12];
    double& p2MASS_2 = Pars[13];
*/
double MS_GaussExp_mT_Simple::Eval(double* Pars){
    double Result=0;
    Parameters[1]=Pars[1];
    Parameters[4] = Tau[0];
    Parameters[5] = 1.-Weight_R[0];
    Parameters[6] = MassR[0];
    Parameters[7] = Mass[0];
    Parameters[8] = MassD[0];
    Parameters[9] = Tau[1];
    Parameters[10] = 1.-Weight_R[1];
    Parameters[11] = MassR[1];
    Parameters[12] = Mass[1];
    Parameters[13] = MassD[1];
    double Rad_mT;
    for(unsigned umt=0; umt<Num_mT; umt++){
        if(!Weight_mT[umt]) continue;
        Rad_mT = Linear_mT+Slope_mT*Mean_mT[umt];
        Parameters[3] = Rad_mT;
        Result += Weight_mT[umt]*GaussExpTotSimple_2body(Parameters);
    }
if(Result!=Result || Result>1 || Result<0){
    /*
    printf("Result=%e\n",Result);
    for(unsigned uPar=0; uPar<14; uPar++){
        printf("   Parameters[%u]=%f\n",uPar,Parameters[uPar]);
    }
    */
    Result = 0;
}


    return Result;
}
