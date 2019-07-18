#include "DLM_Source.h"
#include "DLM_CRAB_PM.h"
#include "DLM_Integration.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_MathFunctions.h"
#include "DLM_Histo.h"
#include "CATSconstants.h"

#include "math.h"

//!TEST
//#include <fstream>
//#include <stdio.h>

double GaussSourceTF1(double* x, double* Pars){
    double& Radius = *x;
    double& Size = Pars[0];
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

//2-particle
double GaussSource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
//printf(" G-> r=%.2f, s=%.2f => %.2f\n",Radius,Size,4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size)));
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
//same as GaussSource, but we assume that the user calls the function wanting to sample from theta as well.
//since Integral(dTheta) = 2 for a flat theta distribution and the whole Source function needs to be normalized to 1,
//in order to preserve this we should divide the whole function by 2
double GaussSourceTheta(double* Pars){
    return 0.5*GaussSource(Pars);
}

//double CauchySource(const double& Radius, const double& Size){
//single particle
double CauchySource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 2.97*Size*sqrt(2)*Radius*Radius/Pi*pow(Radius*Radius+0.125*2.97*2.97*Size*Size,-2.);
}

double CauchySource_v2(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 12.*Radius*Size*sqrt(2)/(Pi*Pi*sqrt(Radius*Radius+4.*Size*Size)*(Radius*Radius+6.*Size*Size))*
            atan(sqrt(2.+0.5*pow(Radius/Size,2.))*Radius/Size/sqrt(2));
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

//printf(" G1=%.1f x %e\n",Weight1, 4.*Pi*Radius*Radius*pow(4.*Pi*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1)));
//printf(" G2=%.1f x %e\n",1.-Weight1, 4.*Pi*Radius*Radius*pow(4.*Pi*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2)));

    return      Weight1 *4.*Pi*Radius*Radius*pow(4.*Pi*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*4.*Pi*Radius*Radius*pow(4.*Pi*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2));

}

double GaussCauchySource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];

    double& Size1 = Pars[3];
    double& Size2 = Pars[4];
    double& Weight1 = Pars[5];

//printf(" G1=%.1f x %e\n",Weight1, 4.*Pi*Radius*Radius*pow(4.*Pi*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1)));
//printf(" G2=%.1f x %e\n",1.-Weight1, 4.*Pi*Radius*Radius*pow(4.*Pi*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2)));

    return      Weight1 *4.*Pi*Radius*Radius*pow(4.*Pi*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*2.*Size2*Radius*Radius/Pi*pow(Radius*Radius+0.25*Size2*Size2,-2.);

}


double LevyIntegral3D_2particle(double* Pars){
    //just a dummy, we do not expect angular dependence, but need the memory for the integration variable
    double& IntVar = Pars[2];
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double Dim = 3;
    //double RadSzRatio = Radius/Scale;
    if(Radius==0) return 0;
    if(IntVar==0) return 0;
//printf("IntVar=%f; Radius=%f; Scale=%f; Stability=%f\n",IntVar,Radius,Scale,Stability);
//printf(" DLM_Bessel1(Dim*0.5-1.,Radius*IntVar)=%f\n",DLM_Bessel1(Dim*0.5-1.,Radius*IntVar));
//double RetVal = pow(Radius*IntVar,Dim*0.5)*DLM_Bessel1(Dim*0.5-1.,Radius*IntVar)*exp(-pow(Scale*2./Stability*IntVar,Stability));
//printf(" RetVal=%f\n",RetVal);
    //return IntVar*IntVar/sqrt(IntVar*RadSzRatio)*DLM_Bessel1(0.5,IntVar*RadSzRatio)*exp(-pow(IntVar,Stability));
    return pow(Radius*IntVar,Dim*0.5)*DLM_Bessel1(Dim*0.5-1.,Radius*IntVar)*exp(-pow(Scale*2./Stability*IntVar,Stability));
//return RetVal;
}

//so I take the definition from something called
//"Multivariate stable densities and distribution functions: general and elliptical case", equation 10
//the normalization is properly done, however for whatever reason if I simulate R=sqrt(X2+Y2+Z2) with X,Y,Z being Levy,
//and they are simulated as according to Wikipedia, I do get some difference, which is corrected for by
//dividing the Scale /= sqrt(Stability) (of equation 10).
//Another funny thing I noticed, if you simulate R=sqrt(X2+Y2+Z2) and correct with Scale *= 2./sqrt(Stability),
//you end up with the distribution of R=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2).
//Effectively to modify this function we need to only make Scale *= 2./Stability
//THIS IS ALL IMPLEMENTED IN HERE, SO THAT WE GET THE DISTRIBUTION CORRESPONDING TO R=sqrt((X1-X2)^2+(Y1-Y2)^2+(Z1-Z2)^2)
double LevySource3D_2particle(double* Pars){
    double& Radius = Pars[1];
    //double& Scale = Pars[3];
    const double Dim = 3;
    double& Scale = Pars[3];
    double& Stability = Pars[4];
//static unsigned NumFunctionCalls=0;
//NumFunctionCalls++;
//if(NumFunctionCalls%1000==0)
//printf("Function call Nr. %u\n",NumFunctionCalls);
//if(Stability==1.01){
//printf("Radius=%f; Scale=%f; Stability=%f\n",Pars[1],Pars[3],Pars[4]);
//}

    if(Stability==1){
        return 2.97*2.*Scale*sqrt(2)*Radius*Radius/Pi*pow(Radius*Radius+0.5*2.97*2.97*Scale*Scale,-2.);
    }
    else if(Stability==2){
        return 4.*Pi*Radius*Radius*pow(4.*Pi*Scale*Scale,-1.5)*exp(-(Radius*Radius)/(4.*Scale*Scale));
    }

    DLM_INT_SetFunction(LevyIntegral3D_2particle,Pars,2);
    if(Radius==0) return 0;
    //return DLM_INT_aSimpsonWiki(0.,16.+Radius,1e-8,128)*pow(2.*Pi*pow(Scale*Scale/3.,3.),-1.5)*4.*Pi*Radius*Radius;
    unsigned NSteps;
    if(Radius>108) NSteps = 1024;
    else if(Radius>36) NSteps = 512;
    else if(Radius>12) NSteps = 256;
    else if(Radius>4) NSteps = 128;
    else NSteps = 64;

    //x2+y2+z2=r2 => r2=3x2 x2=r2/3 x2y2z2=r6/27 = r2
    //return DLM_INT_SimpsonWiki(0.,16.,NSteps)*pow(2.*Pi*pow(Scale*Scale/3.,3.),-1.5)*4.*Pi*Radius*Radius;
    double ReturnVal = DLM_INT_SimpsonWiki(0.,16.,NSteps)*2./(pow(2.,Dim*0.5)*exp(gammln(Dim*0.5)));
//if(ReturnVal!=ReturnVal)
//printf("ReturnVal=%f\n",ReturnVal);
    return ReturnVal;
}

double LevyIntegral3D_single(double* Pars){
    //just a dummy, we do not expect angular dependence, but need the memory for the integration variable
    double& IntVar = Pars[2];
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double Dim = 3;
    if(Radius==0) return 0;
    if(IntVar==0) return 0;
    return pow(Radius*IntVar,Dim*0.5)*DLM_Bessel1(Dim*0.5-1.,Radius*IntVar)*exp(-pow(Scale/sqrt(Stability)*IntVar,Stability));
}
//using this function, we get the 3D Levy distribution (single particle)
double LevySource3D_single(double* Pars){
    double& Radius = Pars[1];
    const double Dim = 3;
    double& Scale = Pars[3];
    double& Stability = Pars[4];

    if(Stability==1){
        return 2.97*Scale*sqrt(2)*Radius*Radius/Pi*pow(Radius*Radius+0.125*2.97*2.97*Scale*Scale,-2.);
    }
    else if(Stability==2){
        return 4.*Pi*Radius*Radius*pow(2.*Pi*Scale*Scale,-1.5)*exp(-(Radius*Radius)/(2.*Scale*Scale));
    }

    DLM_INT_SetFunction(LevyIntegral3D_single,Pars,2);
    if(Radius==0) return 0;
    unsigned NSteps;
    if(Radius>108) NSteps = 1024;
    else if(Radius>36) NSteps = 512;
    else if(Radius>12) NSteps = 256;
    else if(Radius>4) NSteps = 128;
    else NSteps = 64;

    double ReturnVal = DLM_INT_SimpsonWiki(0.,16.,NSteps)*2./(pow(2.,Dim*0.5)*exp(gammln(Dim*0.5)));
    return ReturnVal;
}

double LevyIntegral3D(double* Pars){
    //just a dummy, we do not expect angular dependence, but need the memory for the integration variable
    double& IntVar = Pars[2];
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double Dim = 3;
    if(Radius==0) return 0;
    if(IntVar==0) return 0;
    return pow(Radius*IntVar,Dim*0.5)*DLM_Bessel1(Dim*0.5-1.,Radius*IntVar)*exp(-pow(Scale*IntVar,Stability));
}
//using this function, we get the 3D Levy distribution (single particle)
//what is a bit worrying is that this guy give slightly different results compared to my DLM_Random::Stable, but whatever
//This implementation over here is written to be numerically unstable, so maybe that is the reason...
double LevySource3D(double* Pars){
    double& Radius = Pars[1];
    const double Dim = 3;
    double& Scale = Pars[3];
    double& Stability = Pars[4];

    if(Stability==1){
        return 2.97*Scale*sqrt(2)*Radius*Radius/Pi*pow(Radius*Radius+0.125*2.97*2.97*Scale*Scale,-2.);
    }
    else if(Stability==2){
        return 4.*Pi*Radius*Radius*pow(4.*Pi*Scale*Scale,-1.5)*exp(-(Radius*Radius)/(4.*Scale*Scale));
    }

    DLM_INT_SetFunction(LevyIntegral3D,Pars,2);
    if(Radius==0) return 0;
    unsigned NSteps;
    if(Radius>108) NSteps = 1024;
    else if(Radius>36) NSteps = 512;
    else if(Radius>12) NSteps = 256;
    else if(Radius>4) NSteps = 128;
    else NSteps = 64;

    double ReturnVal = DLM_INT_SimpsonWiki(0.,16.,NSteps)*2./(pow(2.,Dim*0.5)*exp(gammln(Dim*0.5)));
    return ReturnVal;
}

/*
double LevySource_A(double* Pars){
    const unsigned MaxIter = 64;
    const unsigned MinIter = 4;
    const double Epsilon = 1e-10;
    double& Radius = Pars[1];
    double& Size = Pars[3];
    double& Stability = Pars[4];
    double Dim = round(Pars[5]);

    double Result=0;
    double Increase;
    double kIter;
    //double GAMMA;

    if(Stability<=0 || Stability>2) return 0;
    else if(Radius/Size>3){
        double k1 = pow(2.,Stability)*sin(Pi*Stability*0.5)*exp(gammln((Stability+2.)*0.5))*exp(gammln((Stability+Dim)*0.5))/
                    (exp(gammln(Dim*0.5))*Pi*Stability*0.5);
        Result = Stability*k1*pow(Size/Radius,Stability)/Radius;
    }
    else if(Stability<1){

    }
    //Cauchy
    else if(Stability==1){

    }
    else if(Stability<2){
        for(unsigned uIter=0; uIter<MaxIter; uIter++){
            kIter = uIter;
            Increase=(uIter%2==0?1:-1)*exp(gammln((2.*kIter+Dim)/Stability)-gammln((2.*kIter+Dim)*0.5))*pow(0.5*Radius/Size,2.*kIter+Dim-1.)/
                        (factrl(uIter));
            Result+=Increase;
            if(fabs(Increase/Result)<Epsilon && uIter>=MinIter) {
                printf("uIter=%u (%f)\n",uIter,Radius);
                printf(" Increase=%f (%f)\n",Increase,Result);

                break;
            }
        }
        printf(" Result(%f)=%f\n",Radius,Result);
        Result *= 2./(Stability*Size*exp(gammln(Dim*0.5)));
        printf(" Result(%f)=%f\n",Radius,Result);
    }
    //Gauss
    else{

    }
    return Result;
}
*/

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
//4.*Pi*Radius*Radius*pow(4.*Pi*rOut*rOut,-1.5)*exp(-(Radius*Radius)/(4.*rOut*rOut))
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
    RAD = TKM?log(RAD*TKM/80.+1.)*80./TKM-TKM*atan(1.5*RAD/TKM)*2./Pi:RAD;
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
               pow(E,pow(RAD,2)/(4.*pow(SIG,2)) + pow(SIG,2)/pow(TKM2,2))*sqrt(Pi)*(-4*pow(SIG,4) + 2*pow(SIG,2)*(RAD - 3*TKM2)*TKM2 + RAD*pow(TKM2,3))*
                (erf(RAD/(2.*SIG) - SIG/TKM2) + erf(SIG/TKM2))))/(sqrt(Pi)*pow(TKM2,5));
    }
    else if(TKM2==0){
        return (pow(E,-pow(RAD,2)/(4.*pow(SIG,2)) - RAD/TKM1)*(4*pow(E,RAD/TKM1)*SIG*TKM1*(pow(SIG,2) + pow(TKM1,2)) -
               2*pow(E,pow(RAD,2)/(4.*pow(SIG,2)))*SIG*TKM1*(2*pow(SIG,2) + TKM1*(-RAD + 2*TKM1)) +
               pow(E,pow(RAD,2)/(4.*pow(SIG,2)) + pow(SIG,2)/pow(TKM1,2))*sqrt(Pi)*(-4*pow(SIG,4) + 2*pow(SIG,2)*(RAD - 3*TKM1)*TKM1 + RAD*pow(TKM1,3))*
                (erf(RAD/(2.*SIG) - SIG/TKM1) + erf(SIG/TKM1))))/(sqrt(Pi)*pow(TKM1,5));
    }
//!shouldn't I take TKM=TKM1+TKM2 for the computation, verify!
    else if(TKM1==TKM2){
static bool WARNING = true;
if(WARNING) printf("WARNING: shouldn't I take TKM=TKM1+TKM2 for the computation, verify!\n");
WARNING=false;
        return (pow(E,(pow(SIG,2) - RAD*TKM1)/pow(TKM1,2))*((-2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM1,2))*TKM1) + (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM1,2)))*pow(SIG,3))/TKM1 -
               (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2))))*pow(SIG,3))/TKM1 + sqrt(Pi)*pow(SIG,2)*erf(SIG/TKM1) +
               (2*sqrt(Pi)*pow(SIG,4)*erf(SIG/TKM1))/pow(TKM1,2) - (sqrt(Pi)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM1)*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1)))/
                (pow(TKM1,3)*abs(-RAD/(2.*SIG) + SIG/TKM1)) - (pow(2*pow(SIG,2) - RAD*TKM1,3)*
                  (-abs(-2*pow(SIG,2) + RAD*TKM1) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*sqrt(Pi)*SIG*TKM1*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1))))/
                (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*pow(SIG,2)*pow(TKM1,4)*pow(abs(-RAD/(2.*SIG) + SIG/TKM1),3))))/(sqrt(Pi)*pow(SIG,2)*TKM1);
    }
    else{
        return (pow(E,(pow(SIG,2) - RAD*TKM1)/pow(TKM1,2))*((-2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM1,2))*TKM1) + (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM1,2)))*pow(SIG,3))/TKM1 -
        (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2))))*pow(SIG,3))/TKM1 + sqrt(Pi)*pow(SIG,2)*erf(SIG/TKM1) +
        (2*sqrt(Pi)*pow(SIG,4)*erf(SIG/TKM1))/pow(TKM1,2) - (sqrt(Pi)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM1)*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1)))/
         (pow(TKM1,3)*abs(RAD/(2.*SIG) - SIG/TKM1)) - (pow(2*pow(SIG,2) - RAD*TKM1,3)*
           (-abs(-2*pow(SIG,2) + RAD*TKM1) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*sqrt(Pi)*SIG*TKM1*erf(abs(-2*pow(SIG,2) + RAD*TKM1)/(2.*SIG*TKM1))))/
         (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM1,2)/(4.*pow(SIG,2)*pow(TKM1,2)))*pow(SIG,2)*pow(TKM1,4)*pow(abs(RAD/(2.*SIG) - SIG/TKM1),3))) +
     pow(E,(pow(SIG,2) - RAD*TKM2)/pow(TKM2,2))*((2*pow(SIG,3))/(pow(E,pow(SIG,2)/pow(TKM2,2))*TKM2) - (4*(-1 + pow(E,-pow(SIG,2)/pow(TKM2,2)))*pow(SIG,3))/TKM2 +
        (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2))))*pow(SIG,3))/TKM2 - sqrt(Pi)*pow(SIG,2)*erf(SIG/TKM2) -
        (2*sqrt(Pi)*pow(SIG,4)*erf(SIG/TKM2))/pow(TKM2,2) + (sqrt(Pi)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM2)*erf(abs(-2*pow(SIG,2) + RAD*TKM2)/(2.*SIG*TKM2)))/
         (pow(TKM2,3)*abs(RAD/(2.*SIG) - SIG/TKM2)) + (pow(2*pow(SIG,2) - RAD*TKM2,3)*
           (-abs(-2*pow(SIG,2) + RAD*TKM2) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2)))*sqrt(Pi)*SIG*TKM2*erf(abs(-2*pow(SIG,2) + RAD*TKM2)/(2.*SIG*TKM2))))/
         (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM2,2)/(4.*pow(SIG,2)*pow(TKM2,2)))*pow(SIG,2)*pow(TKM2,4)*pow(abs(RAD/(2.*SIG) - SIG/TKM2),3))))/(sqrt(Pi)*pow(SIG,2)*(TKM1 - TKM2));
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
    RAD = TKM?log(RAD*TKM/80.+1.)*80./TKM-TKM*atan(1.5*RAD/TKM)*2./Pi:RAD;
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
               (4*(-1 + pow(E,-pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2))))*pow(SIG,3))/TKM + sqrt(Pi)*pow(SIG,2)*erf(SIG/TKM) +
               (2*sqrt(Pi)*pow(SIG,4)*erf(SIG/TKM))/pow(TKM,2) - (sqrt(Pi)*pow(SIG,3)*(2*pow(SIG,2) - RAD*TKM)*erf(abs(-2*pow(SIG,2) + RAD*TKM)/(2.*SIG*TKM)))/
                (pow(TKM,3)*abs(-RAD/(2.*SIG) + SIG/TKM)) - (pow(2*pow(SIG,2) - RAD*TKM,3)*
                  (-abs(-2*pow(SIG,2) + RAD*TKM) + pow(E,pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2)))*sqrt(Pi)*SIG*TKM*erf(abs(-2*pow(SIG,2) + RAD*TKM)/(2.*SIG*TKM))))/
                (8.*pow(E,pow(-2*pow(SIG,2) + RAD*TKM,2)/(4.*pow(SIG,2)*pow(TKM,2)))*pow(SIG,2)*pow(TKM,4)*pow(abs(-RAD/(2.*SIG) + SIG/TKM),3))))/(sqrt(Pi)*pow(SIG,2)*TKM);
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
printf("r=%f -> %f\n",Pars[1],GaussExpTotSimple(PARS));
for(int i=0; i<14; i++){
    printf(" Pars[%i]=%f\n",i,Pars[i]);
}
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

/*
void MS_Gauss::SetParameter(const unsigned& WhichPar, const double& Value){
    printf("MS_Gauss::SetParameter is a DUMMY at the moment!\n");
}
double MS_Gauss::Eval(double* Pars){
    double& Radius = Pars[1];
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}
*/

MS_GaussExp_mT_Simple::MS_GaussExp_mT_Simple(){
    Num_mT = 0;
    Mean_mT = NULL;
    Weight_mT = NULL;
    Linear_mT = 0;
    Slope_mT = 0;
    FunctionValue = NULL;
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
    if(FunctionValue)delete[]FunctionValue;FunctionValue=NULL;
    if(Mean_mT)delete[]Mean_mT;Mean_mT=NULL;
    if(Weight_mT)delete[]Weight_mT;Weight_mT=NULL;
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
    if(FunctionValue)delete[]FunctionValue;FunctionValue=NULL;
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
void MS_GaussExp_mT_Simple::SetCustomFunction(const unsigned& umt, const double& value){
    if(umt>=Num_mT){
        printf("\033[1;33mWARNING:\033[0m Current number of mT bins is %u and you attempt to set bin nr. %u\n",Num_mT,umt);
        return;
    }
    if(!FunctionValue){
        FunctionValue = new double [Num_mT];
    }
    FunctionValue[umt] = value;
}
void MS_GaussExp_mT_Simple::RemoveCustomFunction(){
    if(FunctionValue){
        delete[]FunctionValue;
        FunctionValue=NULL;
    }
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
void MS_GaussExp_mT_Simple::SetParameter(const unsigned& WhichPar, const double& Value){
    printf("MS_GaussExp_mT_Simple::SetParameter is a DUMMY at the moment!\n");
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

//note that the radius cannot be fitted, as it is taken from the mT slope
double MS_GaussExp_mT_Simple::Eval(double* Pars){
    double Result=0;
    Parameters[1] = Pars[1];
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
    double& Rad_mT = Parameters[3];
    for(unsigned umt=0; umt<Num_mT; umt++){
        if(!Weight_mT[umt]) continue;
        if(FunctionValue) Rad_mT = FunctionValue[umt];
        else Rad_mT = Linear_mT+Slope_mT*Mean_mT[umt];
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
double MS_GaussExp_mT_Simple::EvalROOT(double* x, double* Pars){
    double Result=0;
    Parameters[1] = x[0];
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
    double& Rad_mT = Parameters[3];
    for(unsigned umt=0; umt<Num_mT; umt++){
        if(!Weight_mT[umt]) continue;
        if(FunctionValue) Rad_mT = FunctionValue[umt];
        else Rad_mT = Linear_mT+Slope_mT*Mean_mT[umt];
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
unsigned MS_GaussExp_mT_Simple::GetNumPars(){
    return 14;
}

DLM_StableDistribution::DLM_StableDistribution(const unsigned& numgridpts):NumGridPts(numgridpts){
    Histo = NULL;
    RanGen = new DLM_Random(NumGridPts);
    //NumIter = 65536;
    NumIter = 131072*4;
    Generated = false;
}
DLM_StableDistribution::~DLM_StableDistribution(){
    if(Histo) {delete Histo; Histo=NULL;}
    if(RanGen) {delete RanGen; RanGen=NULL;}
}
void DLM_StableDistribution::SetStability(const double& val){
    if(Stability==val) return;
    if(val<=0||val>2){
        Stability=2;
    }
    else{
        Stability=val;
    }
    Generated = false;
}
void DLM_StableDistribution::SetLocation(const double& val){
    if(Location==val) return;
    if(val<0||val>2){
        Location=2;
    }
    else{
        Location=val;
    }
    Generated = false;
}
void DLM_StableDistribution::SetScale(const double& val){
    if(Scale==val) return;
    if(val<0||val>2){
        Scale=2;
    }
    else{
        Scale=val;
    }
    Generated = false;
}
void DLM_StableDistribution::SetSkewness(const double& val){
    if(Skewness==val) return;
    if(val<0||val>2){
        Skewness=2;
    }
    else{
        Skewness=val;
    }
    Generated = false;
}
void DLM_StableDistribution::SetNumIter(const unsigned& val){
    NumIter = val;
    if(NumIter<100) NumIter=100;
}
void DLM_StableDistribution::Generate(const double& stability, const double& location, const double& scale, const double& skewness){
//printf("Hello there, stability=%e\n",stability);
    if(!Histo){
        Histo = new DLM_Histo1D<double>(NumGridPts,0,64);
    }
    double RanVal;
    if(RanGen) delete RanGen;
    RanGen = new DLM_Random(NumGridPts);
//printf("Hi again\n");
    for(unsigned uBin=0; uBin<NumIter; uBin++){
//printf(" uBin=%u\n",uBin);
        RanVal = RanGen->StableDiffR(3,stability,location,scale,skewness);
//printf(" RanVal=%f\n",RanVal);
        Histo->AddAt(RanVal,1.);
    }
    Histo->Scale(1./double(NumIter));
    Histo->ScaleToBinWidth();
    Generated = true;
//printf("Nothing to see here!\n");
}
void DLM_StableDistribution::SetParameter(const unsigned& WhichPar, const double& Value){
    switch(WhichPar){
    case 0 : SetStability(Value); break;
    case 1 : SetLocation(Value); break;
    case 2 : SetScale(Value); break;
    case 3 : SetSkewness(Value); break;
    default : break;
    }
}
double DLM_StableDistribution::Eval(double* Pars){
    double& rVal=Pars[1];
    SetStability(Pars[3]);
    SetLocation(Pars[4]);
    SetScale(Pars[5]);
    SetSkewness(Pars[6]);

//if(fabs(Stability<1e-128))
//for(int i=0; i<7; i++){
//    printf("%i = %e\n",i,Pars[i]);
//}

    if(rVal<=0) return 0;
    if(!Generated){
        Generate(Stability,Location,Scale,Skewness);
    }
//printf("I have been called!\n");
//printf( "rVal=%f Histo->Eval(rVal)=%f\n",rVal,Histo->Eval(rVal));
    return Histo->Eval(rVal);
}
unsigned DLM_StableDistribution::GetNumPars(){
    return 4;
}

DLM_CleverLevy::DLM_CleverLevy(){
    NumPtsStability = 64;
    MinStability=1;
    MaxStability=2;
    NumPtsScale = 128;
    MinScale=0.25;
    MaxScale=4;
    NumPtsRad = 512;
    MinRad = 0;
    MaxRad = 64;
    Histo = NULL;
    Type = 1;
}
DLM_CleverLevy::~DLM_CleverLevy(){
    if(Histo) {delete Histo;Histo=NULL;}
}
void DLM_CleverLevy::InitStability(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsStability==numPts&&MinStability==minVal&&MaxStability==maxVal) return;
    Reset();
    NumPtsStability=numPts;
    MinStability=minVal;
    MaxStability=maxVal;
}
void DLM_CleverLevy::InitScale(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsScale==numPts&&MinScale==minVal&&MaxScale==maxVal) return;
    Reset();
    NumPtsScale=numPts;
    MinScale=minVal;
    MaxScale=maxVal;
}
void DLM_CleverLevy::InitRad(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsRad==numPts&&MinRad==minVal&&MaxRad==maxVal) return;
    Reset();
    NumPtsRad=numPts;
    MinRad=minVal;
    MaxRad=maxVal;
}
void DLM_CleverLevy::InitType(const int& type){
    if(type<0 || type>1) Type = 1;
    Type = type;
}
double DLM_CleverLevy::RootEval(double* x, double* Pars){
    double PARS[5];
    PARS[1] = *x;
    PARS[3] = Pars[0];
    PARS[4] = Pars[1];
    return Eval(PARS);
}
double DLM_CleverLevy::Eval(double* Pars){
    if(!Histo) {Init();}
    if(!Histo) return -1;
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double RSS[3] = {Radius,Scale,Stability};
    int RadBin = Histo->GetBin(0,Radius);
    int ScaleBin = Histo->GetBin(1,Scale);
    int StabilityBin = Histo->GetBin(2,Stability);
    unsigned WhichBin[3];
    double BIN_PARS[5];
//printf("  Called with: r=%f; σ=%f; α=%f\n",Radius,Scale,Stability);
//printf("  RadBin=%u; ScaleBin=%u; StabilityBin=%u\n",RadBin,ScaleBin,StabilityBin);
    for(int iBin0=RadBin-1; iBin0<=RadBin+1; iBin0++){
        if(iBin0<0||iBin0>=int(Histo->GetNbins(0))) continue;
        WhichBin[0] = iBin0;
        BIN_PARS[1] = Histo->GetBinCenter(0,iBin0);
        for(int iBin1=ScaleBin-1; iBin1<=ScaleBin+1; iBin1++){
            if(iBin1<0||iBin1>=int(Histo->GetNbins(1))) continue;
            WhichBin[1] = iBin1;
            BIN_PARS[3] = Histo->GetBinCenter(1,iBin1);
            for(int iBin2=StabilityBin-1; iBin2<=StabilityBin+1; iBin2++){
                if(iBin2<0||iBin2>=int(Histo->GetNbins(2))) continue;
                WhichBin[2] = iBin2;
                BIN_PARS[4] = Histo->GetBinCenter(2,iBin2);
                if(Histo->GetBinContent(WhichBin)>=0.99e6){
                    #pragma omp critical
                    {
//printf("   Setting up bin %i %i %i\n",iBin0,iBin1,iBin2);
                    if(Type==0) {Histo->SetBinContent(WhichBin,LevySource3D_single(BIN_PARS));}
                    else if(Type==1) {Histo->SetBinContent(WhichBin,LevySource3D_2particle(BIN_PARS));}
                    else {Histo->SetBinContent(WhichBin,LevySource3D(BIN_PARS));}
                    }
                }
            }
        }
    }
/*
printf("  I will try to evaluate at r=%f; σ=%f; α=%f\n",Radius,Scale,Stability);
printf("  The REAL value is: %f\n",LevySource3D_2particle(Pars));
printf("  The relevant bins are: %i(%f) %i(%f) %i(%f)\n",
       RadBin,Histo->GetBinCenter(0,RadBin),
       ScaleBin,Histo->GetBinCenter(1,ScaleBin),
       StabilityBin,Histo->GetBinCenter(2,StabilityBin));
WhichBin[0] = RadBin;
WhichBin[1] = ScaleBin;
WhichBin[2] = StabilityBin;
printf("  The value in this bin is: %f\n",Histo->GetBinContent(WhichBin));
double PARSTMP[5];
PARSTMP[1] = Histo->GetBinCenter(0,RadBin);
PARSTMP[3] = Histo->GetBinCenter(1,ScaleBin);
PARSTMP[4] = Histo->GetBinCenter(2,StabilityBin);
printf("  The value in bin SHOULD be: %f\n",LevySource3D_2particle(PARSTMP));
printf("  The Eval value is: %f\n",Histo->Eval(RSS));
printf("---------------------------------------------------\n");
*/

//static unsigned NumFunctionCalls=0;
//NumFunctionCalls++;
//if(NumFunctionCalls%10000==0)
//printf("Function call Nr. %u\n",NumFunctionCalls);

    return Histo->Eval(RSS);
}
unsigned DLM_CleverLevy::GetNumPars(){
    return 2;
}
void DLM_CleverLevy::Reset(){
    if(Histo) {delete Histo;Histo=NULL;}
}
void DLM_CleverLevy::Init(){
    Reset();
    Histo  = new DLM_Histo<double>();
    Histo->SetUp(3);
    if(NumPtsRad==1) {Histo->SetUp(0,NumPtsRad,MinRad,MaxRad);}
    else{
        double BinWidth = (MaxRad-MinRad)/double(NumPtsRad-1);
        Histo->SetUp(0,NumPtsRad,MinRad-BinWidth*0.5,MaxRad+BinWidth*0.5);
    }
    if(NumPtsScale==1) {Histo->SetUp(1,NumPtsScale,MinScale,MaxScale);}
    else{
        double BinWidth = (MaxScale-MinScale)/double(NumPtsScale-1);
        Histo->SetUp(1,NumPtsScale,MinScale-BinWidth*0.5,MaxScale+BinWidth*0.5);
    }
    if(NumPtsStability==1) {Histo->SetUp(2,NumPtsStability,MinStability,MaxStability);}
    else{
        double BinWidth = (MaxStability-MinStability)/double(NumPtsStability-1);
        Histo->SetUp(2,NumPtsStability,MinStability-BinWidth*0.5,MaxStability+BinWidth*0.5);
    }
    Histo->Initialize();
    Histo->AddToAll(1e6);
}




DLM_CleverMcLevyReso::DLM_CleverMcLevyReso(){
    NumPtsStability = 64;
    MinStability=1;
    MaxStability=2;
    NumPtsScale = 128;
    MinScale=0.25;
    MaxScale=4;
    NumPtsRad = 512;
    MinRad = 0;
    MaxRad = 64;
    Histo = NULL;
    Type = 1;
    NumResonances = new unsigned [2];
    ResoWeight = new double* [2];
    ResoMass = new double* [2];
    ResoTau = new double* [2];
    ChildMass0 = new double* [2];
    ChildMass1 = new double* [2];
    SmearResoMomentum = new double* [2];
    SmearResoMass = new bool* [2];
    for(unsigned uParticle=0; uParticle<2; uParticle++){
        NumResonances[uParticle]=0;
        ResoWeight[uParticle]=NULL;
        ResoMass[uParticle]=NULL;
        ResoTau[uParticle]=NULL;
        ChildMass0[uParticle]=NULL;
        ChildMass1[uParticle]=NULL;
        SmearResoMomentum[uParticle]=NULL;
        SmearResoMass[uParticle]=NULL;
    }
    NumMcIter = 1000000;
}
DLM_CleverMcLevyReso::~DLM_CleverMcLevyReso(){
    if(Histo) {delete Histo;Histo=NULL;}
}
void DLM_CleverMcLevyReso::InitStability(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsStability==numPts&&MinStability==minVal&&MaxStability==maxVal) return;
    Reset();
    NumPtsStability=numPts;
    MinStability=minVal;
    MaxStability=maxVal;
}
void DLM_CleverMcLevyReso::InitScale(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsScale==numPts&&MinScale==minVal&&MaxScale==maxVal) return;
    Reset();
    NumPtsScale=numPts;
    MinScale=minVal;
    MaxScale=maxVal;
}
void DLM_CleverMcLevyReso::InitRad(const unsigned& numPts, const double& minVal, const double& maxVal){
//printf(" NPR=%u(%u); MinR=%.3f(%.3f); MaxR=%.3f(%.3f)\n",numPts,NumPtsRad,minVal,MinRad,maxVal,MaxRad);
    if(NumPtsRad==numPts&&MinRad==minVal&&MaxRad==maxVal) return;
    Reset();
    NumPtsRad=numPts;
    MinRad=minVal;
    MaxRad=maxVal;
}
void DLM_CleverMcLevyReso::InitType(const int& type){
    if(type<0 || type>1) Type = 1;
    Type = type;
}
void DLM_CleverMcLevyReso::InitReso(const unsigned& whichparticle, const unsigned& numreso){
    if(whichparticle>=2) {printf("\033[1;33mWARNING:\033[0m You can call InitReso only for particle 0 or 1\n"); return;}
    if(NumResonances[whichparticle]==numreso) return;
    NumResonances[whichparticle] = numreso;
    if(ResoWeight[whichparticle]){delete[]ResoWeight[whichparticle];ResoWeight[whichparticle]=new double[NumResonances[whichparticle]];}
    if(ResoMass[whichparticle]){delete[]ResoMass[whichparticle];ResoMass[whichparticle]=new double[NumResonances[whichparticle]];}
    if(ResoTau[whichparticle]){delete[]ResoTau[whichparticle];ResoTau[whichparticle]=new double[NumResonances[whichparticle]];}
    if(ChildMass0[whichparticle]){delete[]ChildMass0[whichparticle];ChildMass0[whichparticle]=new double[NumResonances[whichparticle]];}
    if(ChildMass1[whichparticle]){delete[]ChildMass1[whichparticle];ChildMass1[whichparticle]=new double[NumResonances[whichparticle]];}
    if(SmearResoMomentum[whichparticle]){delete[]SmearResoMomentum[whichparticle];SmearResoMomentum[whichparticle]=new double[NumResonances[whichparticle]];}
    if(SmearResoMass[whichparticle]){delete[]SmearResoMass[whichparticle];SmearResoMass[whichparticle]=new bool[NumResonances[whichparticle]];}
    if(Histo) Histo->SetBinContentAll(1e6);
}
void DLM_CleverMcLevyReso::SetUpReso(const unsigned& whichparticle, const unsigned& whichreso, const double& weight, const double& mass, const double& tau, const double& mass0, const double& mass1, const double& momSmear, const bool& massSmear){
    if(whichparticle>=2) {printf("\033[1;33mWARNING:\033[0m You can call SetUpReso only for particle 0 or 1\n"); return;}
    if(whichreso>=NumResonances[whichparticle]) {printf("\033[1;33mWARNING:\033[0m Only %u number of resonances are currently allowed. Change using InitReso\n",NumResonances[whichparticle]); return;}
    if(!ResoWeight[whichparticle]){ResoWeight[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ResoMass[whichparticle]){ResoMass[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ResoTau[whichparticle]){ResoTau[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ChildMass0[whichparticle]){ChildMass0[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ChildMass1[whichparticle]){ChildMass1[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!SmearResoMomentum[whichparticle]){SmearResoMomentum[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!SmearResoMass[whichparticle]){SmearResoMass[whichparticle]=new bool[NumResonances[whichparticle]];}
    if(ResoWeight[whichparticle][whichreso]==weight&&ResoMass[whichparticle][whichreso]==mass&&ResoTau[whichparticle][whichreso]==tau&&
       ChildMass0[whichparticle][whichreso]==mass0&&ChildMass1[whichparticle][whichreso]==mass1&&
       SmearResoMomentum[whichparticle][whichreso]==momSmear&&SmearResoMass[whichparticle][whichreso]==massSmear){return;}
    if(Histo) Histo->SetBinContentAll(1e6);
    ResoWeight[whichparticle][whichreso] = weight;
    ResoMass[whichparticle][whichreso] = mass;
    ResoTau[whichparticle][whichreso] = tau*FmToNu;
    ChildMass0[whichparticle][whichreso] = mass0;
    ChildMass1[whichparticle][whichreso] = mass1;
    SmearResoMomentum[whichparticle][whichreso] = momSmear;
    SmearResoMass[whichparticle][whichreso] = massSmear;
//printf("ResoWeight=%f\n",ResoWeight[whichparticle][whichreso]);
}
void DLM_CleverMcLevyReso::InitNumMcIter(const unsigned& numiter){
    if(NumMcIter==numiter) return;
    if(Histo) Histo->SetBinContentAll(1e6);
    NumMcIter=numiter;
}

double DLM_CleverMcLevyReso::RootEval(double* x, double* Pars){
    double PARS[5];
    PARS[1] = *x;
    PARS[3] = Pars[0];
    PARS[4] = Pars[1];
    return Eval(PARS);
}
double DLM_CleverMcLevyReso::Eval(double* Pars){

//static unsigned NumFunctionCalls=0;
//NumFunctionCalls++;
//if(NumFunctionCalls%10000==0)
//printf("DLM_CleverMcLevyReso::Eval: Function call Nr. %u\n",NumFunctionCalls);
//printf("r=%.3f, R=%.3f, alpha=%.3f\n",Pars[1],Pars[3],Pars[4]);

    if(!Histo) {Init();}
    if(!Histo) return -1;
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double RSS[3] = {Radius,Scale,Stability};
    //int RadBin = Histo->GetBin(0,Radius);
    int ScaleBin = Histo->GetBin(1,Scale);
    int StabilityBin = Histo->GetBin(2,Stability);
    unsigned WhichBin[3];
    double par_stability;
    double par_scale;

    WhichBin[0] = 0;
    for(int iBin1=ScaleBin-1; iBin1<=ScaleBin+1; iBin1++){
        if(iBin1<0||iBin1>=int(Histo->GetNbins(1))) continue;
        WhichBin[1] = iBin1;
        par_scale = Histo->GetBinCenter(1,iBin1);
        for(int iBin2=StabilityBin-1; iBin2<=StabilityBin+1; iBin2++){
            if(iBin2<0||iBin2>=int(Histo->GetNbins(2))) continue;
            WhichBin[2] = iBin2;
            par_stability = Histo->GetBinCenter(2,iBin2);
            if(Histo->GetBinContent(WhichBin)>=0.99e6){
                for(unsigned iBin0=0; iBin0<Histo->GetNbins(0); iBin0++){
                    WhichBin[0] = iBin0;
                    #pragma omp critical
                    {
                    Histo->SetBinContent(WhichBin,0);
                    }
                }
                double RAD;
                //this spares us a renormalization of the whole histogram
                double DiffVal = 1./double(NumMcIter);
                //the size of the r-bin
                double BinSize;
                unsigned TotBin;
                //we use the same seed every time to make our fake function as smooth as possible
                //between the different stability-scale bins
                DLM_Random RanGen(11);
                double RanVal;
                double ResoMomentum;
                double TKM;
                //#pragma omp for
                for(unsigned uIter=0; uIter<NumMcIter; uIter++){
//if(uIter%10000==0)
//printf("   %u/%u\n",uIter,NumMcIter);
                    //we simulate random 'core' distance RAD
                    if(Type==0) {RAD = RanGen.StableR(3,par_stability,0,par_scale,0);}
                    else if(Type==1) {RAD = RanGen.StableDiffR(3,par_stability,0,par_scale,0);}
                    else {RAD = RanGen.StableNolan(3,par_stability,0,par_scale,0);}
//printf("par_scale=%f\n",par_scale);
//printf("par_stability=%f\n",par_stability);
//printf("RAD = %f\n",RAD);
//static unsigned RESO=0;
//static unsigned ATTEMPTS=0;
                    //after that we throw a random dice to decide IF each particle is primary or not
                    for(unsigned uParticle=0; uParticle<2; uParticle++){
                        if(!ResoWeight||!ResoWeight[uParticle]) continue;
                        RanVal = RanGen.Uniform(0,1);
                        double CummulativeWeight=0;
                        int WhichReso = -1;
                        for(unsigned uReso=0; uReso<NumResonances[uParticle]; uReso++){
//printf(" RanVal=%f; CW=%f; CW+RW=%f\n",RanVal,CummulativeWeight,CummulativeWeight+ResoWeight[uParticle][uReso]);
                            if(RanVal>=CummulativeWeight&&RanVal<CummulativeWeight+ResoWeight[uParticle][uReso]){
                                WhichReso = int(uReso);
//printf("WhichReso=%i from %.2f<%.2f<%.2f\n",WhichReso,CummulativeWeight,RanVal,CummulativeWeight+ResoWeight[uParticle][uReso]);
//RESO++;
//ATTEMPTS++;
                                break;
                            }
//if(uReso==NumResonances[uParticle]-1){
//ATTEMPTS++;
//printf("WhichReso=%i from %.2f>%.2f\n",WhichReso,RanVal,CummulativeWeight+ResoWeight[uParticle][uReso]);
//}
                            CummulativeWeight+=ResoWeight[uParticle][uReso];
                        }
//printf("RESO/ATTAMPTS = %.3f\n",double(RESO)/double(ATTEMPTS));
                        //if we have a primary particle, we do nothing
                        if(WhichReso==-1) continue;
                        //we make the radius bigger, according to an exponential decay law, in case we have a resonance


                        double SmearedResoMass = ResoMass[uParticle][WhichReso];
                        if(SmearResoMass[uParticle][WhichReso]&&ResoTau[uParticle][WhichReso]){
                            double Gamma = 1./ResoTau[uParticle][WhichReso];
                            do{
                                SmearedResoMass = RanGen.Cauchy(ResoMass[uParticle][WhichReso],Gamma/sqrt(2));
                            }
                            while(SmearedResoMass<ChildMass0[uParticle][WhichReso]+ChildMass1[uParticle][WhichReso]);
                        }

                        //we compute the effective momentum of the resonance
                        ResoMomentum = sqrt(pow(SmearedResoMass,4.)-2.*pow(SmearedResoMass*ChildMass0[uParticle][WhichReso],2.)+
                                            pow(ChildMass0[uParticle][WhichReso],4.)-2.*pow(SmearedResoMass*ChildMass1[uParticle][WhichReso],2.)-
                                            2.*pow(ChildMass0[uParticle][WhichReso]*ChildMass1[uParticle][WhichReso],2.)+
                                            pow(ChildMass1[uParticle][WhichReso],4.))/(2.*ChildMass0[uParticle][WhichReso]);
//printf("ResoMomentum=%f\n",ResoMomentum);
                        double SmearedResoMomentum = ResoMomentum;
                        if(SmearResoMomentum[uParticle][WhichReso]>0){
                            do{
                                SmearedResoMomentum = RanGen.Gauss(ResoMomentum,ResoMomentum*SmearResoMomentum[uParticle][WhichReso]);
//printf("Smear %.2f%% from %.3f to %.3f\n",SmearResoMomentum[uParticle][WhichReso]*100,ResoMomentum,SmearedResoMomentum);
                            }
                            while(SmearedResoMomentum<0);
                        }

                        //! is this the correct def. or should it be 1/TKM???
                        //! I think it should be M/(p*tau), which is the opposite as done before??????? Check for consistency!!!
                        //! also make sure that the ResoTau gets here with the correct units (natural, i.e. 1/MeV)
//printf("Do not trust the results before you resolve this issue!\n");
//printf("ResoTau[uParticle][WhichReso] = %f\n",ResoTau[uParticle][WhichReso]);
                        //TKM = ResoTau[uParticle][WhichReso]*SmearedResoMomentum/ResoMass[uParticle][WhichReso];
                        TKM = SmearedResoMass/(ResoTau[uParticle][WhichReso]*SmearedResoMomentum+1e-64);
                        //we increase RAD following an exponential
//printf(" TKM = %f\n",TKM);
                        RAD += RanGen.Exponential(TKM)*NuToFm;
//printf(" eRAD = %f\n",RAD);
                    }

                    WhichBin[0] = Histo->GetBin(0,RAD);
                    TotBin = Histo->GetTotBin(WhichBin);
                    BinSize = Histo->GetBinSize(0,WhichBin[0]);
                    //we also normalize to the bin size of the radius, so that we get a
                    //pdf (for these particular stability and scale) that is properly normalized
                    if(Histo->GetBinContent(TotBin)>0.99e6){
                        Histo->SetBinContent(TotBin,0);
                    }
                    //BinSize==0, happens when we are on the edges, i.e. outside of our histo
                    //#pragma omp critical
                    {
                    if(BinSize!=0) Histo->Add(TotBin,DiffVal/BinSize);
                    }

//if(Histo->GetBinContent(TotBin)>1e6){
//if(BinSize==0){
//printf("We have problem with BinSize --> WhichBin[0]=%u (DiffVal=%f, GetBinContent=%f)\n",WhichBin[0],DiffVal,Histo->GetBinContent(TotBin));
//}
//if(DiffVal!=DiffVal){
//printf("We have problem with DiffVal\n");
//}
//printf(" fdsojhvkjfdshbgvkjsfdbgvjkds\n");
                }
            }
        }
    }
    double RETVAL = Histo->Eval(RSS);
    //if(RETVAL>0.99e6) return 0;//this happens in case we lack statistics, so most likely in a bin that should be zero

    //if(RETVAL>0.99e6){
    //    printf("RETVAL>0.99e6 = %e\n",RETVAL);
    //    printf(" At: %f; %f; %f\n",RSS[0],RSS[1],RSS[2]);
        //for(int iBin1=ScaleBin-1; iBin1<=ScaleBin+1; iBin1++){
        //    if(iBin1<0||iBin1>=int(Histo->GetNbins(1))) continue;
        //    WhichBin[1] = iBin1;
        //    par_scale = Histo->GetBinCenter(1,iBin1);
        //    for(int iBin2=StabilityBin-1; iBin2<=StabilityBin+1; iBin2++){
        //        WhichBin[0] = Histo->GetBin(0,RSS[0]);
        //        unsigned TotBin = Histo->GetTotBin(WhichBin);
        //        double BinSize = Histo->GetBinSize(0,WhichBin[0]);
        //        printf("TotBin=%u\n",TotBin);
        //        printf("BinSize=%f\n",BinSize);
        //        printf("BC=%f\n",Histo->GetBinContent(TotBin));
        //    }
        //}
    //}
    //if(RETVAL!=RETVAL){
    //    printf("What happened!?\n");
    //}

//printf("RETVAL=%f\n",RETVAL);

    return RETVAL;
}
unsigned DLM_CleverMcLevyReso::GetNumPars(){
    return 2;
}
void DLM_CleverMcLevyReso::Reset(){
//printf("RESET\n");
    if(Histo) {delete Histo;Histo=NULL;}
}
void DLM_CleverMcLevyReso::Init(){
//printf(" Time to init new histo!\n");
    Reset();
    Histo  = new DLM_Histo<double>();
    Histo->SetUp(3);
    if(NumPtsRad==1) {Histo->SetUp(0,NumPtsRad,MinRad,MaxRad);}
    else{
        double BinWidth = (MaxRad-MinRad)/double(NumPtsRad-1);
        Histo->SetUp(0,NumPtsRad,MinRad-BinWidth*0.5,MaxRad+BinWidth*0.5);
    }
    if(NumPtsScale==1) {Histo->SetUp(1,NumPtsScale,MinScale,MaxScale);}
    else{
        double BinWidth = (MaxScale-MinScale)/double(NumPtsScale-1);
        Histo->SetUp(1,NumPtsScale,MinScale-BinWidth*0.5,MaxScale+BinWidth*0.5);
    }
    if(NumPtsStability==1) {Histo->SetUp(2,NumPtsStability,MinStability,MaxStability);}
    else{
        double BinWidth = (MaxStability-MinStability)/double(NumPtsStability-1);
        Histo->SetUp(2,NumPtsStability,MinStability-BinWidth*0.5,MaxStability+BinWidth*0.5);
    }
    Histo->Initialize();
    Histo->AddToAll(1e6);
}
