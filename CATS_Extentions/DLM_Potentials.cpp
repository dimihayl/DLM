
#include "DLM_Potentials.h"
#include "DLM_StefanoPotentials.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include <omp.h>

//int DlmPot=0;
//int DlmPotFlag=0;
DLM_StefanoPotentials*** fV18pot=NULL;
unsigned NumThreads_DLMPOT;

void CleanUpV18Pot(){
    if(fV18pot){
        for(unsigned uThread=0; uThread<NumThreads_DLMPOT; uThread++){
            if(fV18pot[uThread]){
                for(unsigned iPot=0; iPot<30; iPot++){
                    if(fV18pot[uThread][iPot]){delete fV18pot[uThread][iPot]; fV18pot[uThread][iPot]=NULL;}
                }
                delete [] fV18pot[uThread];
                fV18pot[uThread] = NULL;
            }
        }
        delete [] fV18pot;
        fV18pot = NULL;
    }

}

double ZeroPotential(double* Radius){
    return 0;
}

double DoubleGaussSum(double* Pars){
    return Pars[2]*exp(-pow(Pars[0]/Pars[3],2))+Pars[4]*exp(-pow(Pars[0]/Pars[5],2));
}

//V0*exp(-r^2/β0^2)+V1*exp(-r^2/β1^2)+V2*exp(-r^2/β2^2)
//[0] - r; [1] = k; [2]=V0; [3]=μ0; [4]=V1; [5]=μ1; [6]=V2; [7]=μ2
double TripleGaussSum(double* Pars){
    return Pars[2]*exp(-pow(Pars[0]/Pars[3],2))+Pars[4]*exp(-pow(Pars[0]/Pars[5],2))+Pars[6]*exp(-pow(Pars[0]/Pars[7],2));
}

/*
double CustomUsmaniStefano1(const double& Radius,const int& ipart){

    const double mju = 0.98;
    const double vZero = -49.2;
    double PotVal;

    if(ipart==0){
        PotVal = mju*mju*vZero/Power(CosH(mju*Radius),2);
    }
    else{
        PotVal = 0;
        //Printf("KUREC!");
        //PotVal = mju*mju*vZero/Power(CosH(mju*Radius),2);
    }

    //if(ipart==0) hCustomUsmaniStefano1->SetBinContent(hCustomUsmaniStefano1->FindBin(Radius),PotVal);

    return PotVal;
}

double CustomUsmaniStefano2(const double& Radius,const int& ipart){

    const double mju=2.08;
    const double vZero = -57.8;
    double PotVal;
    if(ipart==0){
        PotVal = mju*mju*vZero*(Exp(-mju*Radius)-2*Exp(-2*mju*Radius));
    }
    else{
        PotVal = 0;
        //Printf("KUREC!");
        //PotVal = mju*mju*vZero*(Exp(-mju*Radius)-2*Exp(-2*mju*Radius));
    }

    //if(ipart==0) hCustomUsmaniStefano2->SetBinContent(hCustomUsmaniStefano2->FindBin(Radius),PotVal);

    return PotVal;
}
*/
//f(x) = (x > 0 ? a*exp(-b*x*x)+c*(1-exp(-d*x*x))**1*(exp(-e*x)/x)**1 : a) #  for 1S0 effective
//g(x) = (x > 0 ? a*exp(-b*x*x)+f*exp(-g*x*x)+c*(1-exp(-d*x*x))**1*(exp(-e*x)/x)**2 : a+f+c*d)  # for 3S1
//Other pars --> used for the dummy I=1 case, in which we use a shifted 3P1 of the AV18, OtherPars[0] contains info about the shift
double LatticePots_pXi(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){

    double parA;
    double parB;
    double parC;
    double parD;
    double parE;
    double parF;
    double parG;
    double Result;
    const double& rad=Radius[0];
    const double rad2=rad*rad;

    //I=0; 1S0
    if(IsoSpin==0&&Spin==0&&AngMom==0){
        switch(DlmPotFlag){
        case 9 :
            parA = 74.71316;
            parB = 5.63683;
            parC = -95.357;
            parD = 4.2328;
            parE = 1.2523;
            parF = 0;
            parG = 0;
            break;
        case 10 :
            parA = 59.9704;
            parB = 6.22966;
            parC = -124.526;
            parD = 2.86352;
            parE = 1.43602;
            parF = 0;
            parG = 0;
            break;
        case 11 :
            parA = 102.822;
            parB = 5.35978;
            parC = -111.948;
            parD = 4.91408;
            parE = 1.37556;
            parF = 0;
            parG = 0;
            break;
        case 12 :
            parA = 155.796;
            parB = 4.80141;
            parC = -111.387;
            parD = 8.28641;
            parE = 1.39018;
            parF = 0;
            parG = 0;
            break;
        default :
            return 0;
        }
 //f(x) = (x > 0 ? a*exp(-b*x*x)+c*(1-exp(-d*x*x))**1*(exp(-e*x)/x)**1 : a) #  for 1S0 effective
        Result = parA*exp(-parB*rad2)+parC*(1.-exp(-parD*rad2))*exp(-parE*rad)/rad;
    }
    //I=0; 3S1
    else if(IsoSpin==0&&Spin==1&&AngMom==0){
        switch(DlmPotFlag){
        case 9 :
            parA = 263.001;
            parB = 46.6953;
            parC = 50.7773;
            parD = 12.7769;
            parE = 1.36322;
            parF = -53.3663;
            parG = 1.00358;
            break;
        case 10 :
            parA = 233.389;
            parB = 48.4379;
            parC = 60.8502;
            parD = 11.8011;
            parE = 1.60436;
            parF = -49.9276;
            parG = 0.932905;
            break;
        case 11 :
            parA = 282.933;
            parB = 45.1205;
            parC = 54.009;
            parD = 12.7096;
            parE = 1.54701;
            parF = -47.2561;
            parG = 0.871518;
            break;
        case 12 :
            parA = 143.473;
            parB = 72.8712;
            parC = 25.0747;
            parD = 32.9093;
            parE = 0.625839;
            parF = -47.8783;
            parG = 0.776256;
            break;
        default :
            return 0;
        }
        Result = parA*exp(-parB*rad2)+parF*exp(-parG*rad2)+parC*(1.-exp(-parD*rad2))*pow(exp(-parE*rad)/rad,2.);
    }
    //I=1; 1S0 (1-6 is 9-14t)
    else if(IsoSpin==1&&Spin==0&&AngMom==0){
        double ShiftedRad;
        switch(DlmPotFlag){
        case 1 :
            parA = 887.15;
            parB = 86.6343;
            parC = 47.8335;
            parD = 44.4054;
            parE = 0.0521;
            parF = -107.262;
            parG = 0.8034;
            break;
        case 2 :
            parA = 1246.97;
            parB = 69.2017;
            parC = 50.681;
            parD = 35.2924;
            parE = 0.1242;
            parF = -106.289;
            parG = 0.809274;
            break;
        case 3 :
            parA = 1117.45;
            parB = 73.714;
            parC = 50.7074;
            parD = 38.0843;
            parE = 0.2032;
            parF = -89.8129;
            parG = 0.78483;
            break;
        case 4 :
            parA = 1208.22;
            parB = 71.7543;
            parC = 50.0453;
            parD = 36.8673;
            parE = 0.179336;
            parF = -90.0595;
            parG = 0.769378;
            break;
        case 5 :
            parA = 1877.9;
            parB = 47.1401;
            parC = 84.0508;
            parD = 13.722;
            parE = 0.498122;
            parF = -119.185;
            parG = 1.17361;
            break;
        case 6 :
            parA = 1951.68;
            parB = 35.6904;
            parC = 217.758;
            parD = 5.44626;
            parE = 0.684346;
            parF = -320.53;
            parG = 1.78494;
            break;
        case 112 :
            ShiftedRad = Radius[0]+OtherPars[0];
            return fDlmPot(NN_AV18,v18_Coupled3P2,1,1,1,1,1,1,&ShiftedRad);
        default :
            return 0;
        }
        Result = parA*exp(-parB*rad2)+parF*exp(-parG*rad2)+parC*(1.-exp(-parD*rad2))*pow(exp(-parE*rad)/rad,2.);
    }
    //I=1; 3S1
    else if(IsoSpin==1&&Spin==1&&AngMom==0){
        double ShiftedRad = Radius[0]+OtherPars[0];
        switch(DlmPotFlag){
        case 1 :
            parA = 646.405;
            parB = 67.0009;
            parC = 36.7729;
            parD = 29.6041;
            parE = 0.269814;
            parF = -78.1609;
            parG = 0.835816;
            break;
        case 2 :
            parA = 725.683;
            parB = 60.6309;
            parC = 41.4768;
            parD = 24.8009;
            parE = 0.383506;
            parF = -77.5009;
            parG = 0.864818;
            break;
        case 3 :
            parA = 854.558;
            parB = 51.7415;
            parC = 52.2627;
            parD = 17.4994;
            parE = 0.60479;
            parF = -75.5611;
            parG = 0.942424;
            break;
        case 4 :
            parA = 933.869;
            parB = 49.3446;
            parC = 52.3054;
            parD = 16.5628;
            parE = 0.491642;
            parF = -92.3852;
            parG = 1.03321;
            break;
        case 5 :
            parA = 880.126;
            parB = 55.0215;
            parC = 43.9326;
            parD = 21.1401;
            parE = 0.359052;
            parF = -86.3012;
            parG = 0.989076;
            break;
        case 6 :
            parA = 775.005;
            parB = 63.7121;
            parC = 46.2741;
            parD = 22.6468;
            parE = 0.366534;
            parF = -90.0947;
            parG = 1.10268;
            break;
        case 112 :
            return fDlmPot(NN_AV18,v18_Coupled3P2,1,1,1,1,1,1,&ShiftedRad);
        default :
            return 0;
        }
        Result = parA*exp(-parB*rad2)+parF*exp(-parG*rad2)+parC*(1.-exp(-parD*rad2))*pow(exp(-parE*rad)/rad,2.);
    }
    else{
        return 0;
    }

    return Result;
}



struct LatticeValues{
    const unsigned NumPots;
    LatticeValues():NumPots(3){
        PAR_V0 = new double* [NumPots];
        PAR_VS = new double* [NumPots];
        PAR_VT = new double* [NumPots];
        PAR_VST = new double* [NumPots];
        for(unsigned uPot=0; uPot<NumPots; uPot++){
            PAR_V0[uPot] = new double [8];
            PAR_VS[uPot] = new double [8];
            PAR_VT[uPot] = new double [8];
            PAR_VST[uPot] = new double [8];
        }
        for(unsigned uFlag=11; uFlag<=13; uFlag++){
            switch(uFlag){
            case 11:
                PAR_V0[uFlag-11][0] = 832.719;
                PAR_V0[uFlag-11][1] = 0.126567;
                PAR_V0[uFlag-11][2] = 306.315;
                PAR_V0[uFlag-11][3] = 0.262349;
                PAR_V0[uFlag-11][4] = 521.285;
                PAR_V0[uFlag-11][5] = 0.461616;
                PAR_V0[uFlag-11][6] = -80.9157;
                PAR_V0[uFlag-11][7] = 9.41638;

                PAR_VS[uFlag-11][0] = -112.713;
                PAR_VS[uFlag-11][1] = 0.119882;
                PAR_VS[uFlag-11][2] = -60.3916;
                PAR_VS[uFlag-11][3] = 0.219904;
                PAR_VS[uFlag-11][4] = -12.971;
                PAR_VS[uFlag-11][5] = 0.440375;
                PAR_VS[uFlag-11][6] = 0;
                PAR_VS[uFlag-11][7] = 0;

                PAR_VT[uFlag-11][0] = 205.93;
                PAR_VT[uFlag-11][1] = 0.135111;
                PAR_VT[uFlag-11][2] = 93.3825;
                PAR_VT[uFlag-11][3] = 0.278859;
                PAR_VT[uFlag-11][4] = 26.9143;
                PAR_VT[uFlag-11][5] = 0.588477;
                PAR_VT[uFlag-11][6] = 0;
                PAR_VT[uFlag-11][7] = 0;

                PAR_VST[uFlag-11][0] = -79.774;
                PAR_VST[uFlag-11][1] = 0.135531;
                PAR_VST[uFlag-11][2] = -32.564;
                PAR_VST[uFlag-11][3] = 0.275606;
                PAR_VST[uFlag-11][4] = -9.3541;
                PAR_VST[uFlag-11][5] = 0.538997;
                PAR_VST[uFlag-11][6] = -1.75591;
                PAR_VST[uFlag-11][7] = 2.88912;
                break;
            case 12:
                PAR_V0[uFlag-11][0] = 800.944;
                PAR_V0[uFlag-11][1] = 0.125499;
                PAR_V0[uFlag-11][2] = 340.209;
                PAR_V0[uFlag-11][3] = 0.25353;
                PAR_V0[uFlag-11][4] = 528.537;
                PAR_V0[uFlag-11][5] = 0.453636;
                PAR_V0[uFlag-11][6] = -74.9729;
                PAR_V0[uFlag-11][7] = 9.89426;

                PAR_VS[uFlag-11][0] = -124.804;
                PAR_VS[uFlag-11][1] = 0.122483;
                PAR_VS[uFlag-11][2] = -48.2802;
                PAR_VS[uFlag-11][3] = 0.234785;
                PAR_VS[uFlag-11][4] = -12.4604;
                PAR_VS[uFlag-11][5] = 0.450708;
                PAR_VS[uFlag-11][6] = 0;
                PAR_VS[uFlag-11][7] = 0;

                PAR_VT[uFlag-11][0] = 170.527;
                PAR_VT[uFlag-11][1] = 0.121487;
                PAR_VT[uFlag-11][2] = 117.874;
                PAR_VT[uFlag-11][3] = 0.224199;
                PAR_VT[uFlag-11][4] = 42.0892;
                PAR_VT[uFlag-11][5] = 0.510129;
                PAR_VT[uFlag-11][6] = 0;
                PAR_VT[uFlag-11][7] = 0;

                PAR_VST[uFlag-11][0] = -90.3105;
                PAR_VST[uFlag-11][1] = 0.148112;
                PAR_VST[uFlag-11][2] = -26.359;
                PAR_VST[uFlag-11][3] = 0.351332;
                PAR_VST[uFlag-11][4] = -4.72814;
                PAR_VST[uFlag-11][5] = 0.707443;
                PAR_VST[uFlag-11][6] = -1.26899;
                PAR_VST[uFlag-11][7] = 2.55046;
                break;
            case 13:
                PAR_V0[uFlag-11][0] = 807.648;
                PAR_V0[uFlag-11][1] = 0.12642;
                PAR_V0[uFlag-11][2] = 359.403;
                PAR_V0[uFlag-11][3] = 0.255362;
                PAR_V0[uFlag-11][4] = 500.769;
                PAR_V0[uFlag-11][5] = 0.451934;
                PAR_V0[uFlag-11][6] = -67.7968;
                PAR_V0[uFlag-11][7] = 10.0506;

                PAR_VS[uFlag-11][0] = -52.9071;
                PAR_VS[uFlag-11][1] = 0.0916017;
                PAR_VS[uFlag-11][2] = -119.322;
                PAR_VS[uFlag-11][3] = 0.164464;
                PAR_VS[uFlag-11][4] = -22.4783;
                PAR_VS[uFlag-11][5] = 0.408204;
                PAR_VS[uFlag-11][6] = 0;
                PAR_VS[uFlag-11][7] = 0;

                PAR_VT[uFlag-11][0] = 124.788;
                PAR_VT[uFlag-11][1] = 0.112992;
                PAR_VT[uFlag-11][2] = 151.176;
                PAR_VT[uFlag-11][3] = 0.194521;
                PAR_VT[uFlag-11][4] = 48.0954;
                PAR_VT[uFlag-11][5] = 0.490352;
                PAR_VT[uFlag-11][6] = 0;
                PAR_VT[uFlag-11][7] = 0;

                PAR_VST[uFlag-11][0] = -82.1509;
                PAR_VST[uFlag-11][1] = 0.15925;
                PAR_VST[uFlag-11][2] = -26.2495;
                PAR_VST[uFlag-11][3] = 0.358359;
                PAR_VST[uFlag-11][4] = -6.31347;
                PAR_VST[uFlag-11][5] = 0.639725;
                PAR_VST[uFlag-11][6] = -1.32727;
                PAR_VST[uFlag-11][7] = 2.67235;
                break;
            default:break;
            }
        }
    }
    ~LatticeValues(){
        for(unsigned uPot=0; uPot<NumPots; uPot++){
            delete [] PAR_V0[uPot];
            delete [] PAR_VS[uPot];
            delete [] PAR_VT[uPot];
            delete [] PAR_VST[uPot];
        }
        delete [] PAR_V0;
        delete [] PAR_VS;
        delete [] PAR_VT;
        delete [] PAR_VST;
    }

    double** PAR_V0;
    double** PAR_VS;
    double** PAR_VT;
    double** PAR_VST;

    double EvalV(const int& iFlag, const int& Element, const double& Rad){
        if(iFlag>=3 || iFlag<0) return 0;
        double* PAR = Element==0?PAR_V0[iFlag]:Element==1?PAR_VS[iFlag]:Element==2?PAR_VT[iFlag]:Element==3?PAR_VST[iFlag]:NULL;
        if(!PAR) return 0;
        if(!Rad || !PAR[1] || !PAR[3] || !PAR[5]) return 0;
        return  PAR[0]*exp(-pow(Rad/PAR[1],2.))+PAR[2]*exp(-pow(Rad/PAR[3],2.))+PAR[4]*exp(-pow(Rad/PAR[5],2.))+
                PAR[6]*pow(1.-exp(-PAR[7]*Rad*Rad),2.)*exp(-2.*146./197.3269602*Rad)/(Rad*Rad);
    }

    double Eval(const int& DlmPotFlag, const double& IsoSpin, const double& Spin, const double& Rad){
        int iFlag = DlmPotFlag-11;
        if(IsoSpin==0 && Spin==0){
            return EvalV(iFlag,0,Rad)-3.*EvalV(iFlag,1,Rad)-3.*EvalV(iFlag,2,Rad)+9.*EvalV(iFlag,3,Rad);
        }
        else if(IsoSpin==0 && Spin==1){
            return EvalV(iFlag,0,Rad)+EvalV(iFlag,1,Rad)-3.*EvalV(iFlag,2,Rad)-3.*EvalV(iFlag,3,Rad);
        }
        else if(IsoSpin==1 && Spin==0){
            return EvalV(iFlag,0,Rad)-3.*EvalV(iFlag,1,Rad)+EvalV(iFlag,2,Rad)-3.*EvalV(iFlag,3,Rad);
        }
        else if(IsoSpin==1 && Spin==1){
            return EvalV(iFlag,0,Rad)+EvalV(iFlag,1,Rad)+EvalV(iFlag,2,Rad)+EvalV(iFlag,3,Rad);
        }
        else return 0;
    }

};


//updated version from 1th October 2018
double LatticePots_pXi_ver2(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){
    static LatticeValues LVAL;
    if(AngMom) return 0;
    return LVAL.Eval(DlmPotFlag,IsoSpin,Spin,Radius[0]);
}







//! The flags are such that the first two digits are the flags to be used for the I0, the second two for I1
//e.g. flags I0 = 12 and I=1 = 6 would be 1206
double LatticePots_pXi_Avg(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){

    double Result=0;
//LatticePots_pXi(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
    Result += 1./8.*LatticePots_pXi(WhichPot,DlmPotFlag/100,0,t2p1,t2p2,0,AngMom,TotMom,Radius,OtherPars);
    Result += 3./8.*LatticePots_pXi(WhichPot,DlmPotFlag/100,0,t2p1,t2p2,1,AngMom,TotMom,Radius,OtherPars);
    Result += 1./8.*LatticePots_pXi(WhichPot,DlmPotFlag%100,1,t2p1,t2p2,0,AngMom,TotMom,Radius,OtherPars);
    Result += 3./8.*LatticePots_pXi(WhichPot,DlmPotFlag%100,1,t2p1,t2p2,1,AngMom,TotMom,Radius,OtherPars);

    return Result;

}

double LatticePots_pXi_SqrtAvg(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){

    double Result=0;

    Result += sqrt(1./8.)*LatticePots_pXi(WhichPot,DlmPotFlag,0,t2p1,t2p2,0,AngMom,TotMom,Radius,OtherPars);
    Result += sqrt(3./8.)*LatticePots_pXi(WhichPot,DlmPotFlag,0,t2p1,t2p2,1,AngMom,TotMom,Radius,OtherPars);
    Result += sqrt(1./8.)*LatticePots_pXi(WhichPot,DlmPotFlag,1,t2p1,t2p2,0,AngMom,TotMom,Radius,OtherPars);
    Result += sqrt(3./8.)*LatticePots_pXi(WhichPot,DlmPotFlag,1,t2p1,t2p2,1,AngMom,TotMom,Radius,OtherPars);

    return Result;

}

//B1*exp(-B2*r*r)+B3*(1-exp(-b4*r*r))(exp(-2*mpi*r)/(r*r))
double LatticePots_pOmega(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){


    const double mpihal = 146.;
    const double hcbar = 197.3;
    double mpihalfm = mpihal/hcbar;

    double B1;
    double B2;
    double B3;
    double B4;
    double Result;
    const double& rad=Radius[0];
    const double rad2=rad*rad;

    //I=1/2; 5S2
        switch(DlmPotFlag){
        case 11 :
          B1 = -306.5;
          B2 = 73.9;
          B3 = -266.0;
          B4 = 0.78;
            break;
        case 12 :
          B1 = -313.0;
          B2 = 81.7;
          B3 = -252.0;
          B4 = 0.85;
            break;
        case 13 :
          B1 = -316.7;
          B2 = 81.9;
          B3 = -237.0;
          B4 = 0.91;
            break;
        case 14 :
          B1 = -296;
          B2 = 64.0;
          B3 = -272.0;
          B4 = 0.76;
            break;
        //potential I
        case 121 :
            B1 = -248.198;
            B2 = 150.475;
            B3 = -157.83;
            B4 = 2.89924;
            mpihalfm = 1.6;
            break;
        //should be same as 12, it is not!!?? (potential II)
        case 122 :
            B1 = -248.198;
            B2 = 150.475;
            B3 = -157.83;
            B4 = 2.89924;
            mpihalfm = 1.09929;
            break;
        //potential III (most binding)
        case 123 :
            B1 = -248.198;
            B2 = 150.475;
            B3 = -157.83;
            B4 = 2.89924;
            mpihalfm = 0.6;
            break;
        default :
            return 0;
        }
 //f(x) = (x > 0 ? a*exp(-b*x*x)+c*(1-exp(-d*x*x))**1*(exp(-e*x)/x)**1 : a) #  for 1S0 effective
        Result = B1*exp(-B2*rad2)+B3*(1.-exp(-B4*rad2))*exp(-2.0*mpihalfm*rad)/rad2;
        return Result;
}

//p-Omega local potential taken from Hyodo et al. arXiv:1805.04024
//For the moment only the Real part is included (see Fig.9 in the paper), so only the real part out of the coefficinets in table 5-column total
double Tetsuo_pOmega(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){
    const double hbarc = 197.3269602;
    const double coeff[9] = {0.14,-13.76,306.81,-2729.49,11744.48,-26288.42,30558.08,-16730.80,3051.85};
    const double mn[9] = {100./hbarc,200./hbarc,300./hbarc,400./hbarc,500./hbarc,600./hbarc,700./hbarc,800./hbarc,900./hbarc};
    const double pig = 3.1415926535897932384626433832795028841971693993751;
    const double lambda = 1000./hbarc;

    double Result=0;
    const double& rad=Radius[0];

    for(unsigned uPar=0; uPar<9; uPar++){
        Result += hbarc*(((1./(4.*pig*rad))*coeff[uPar]*(pow(pow(lambda,2) /
                        (pow(lambda,2)-pow(mn[uPar],2)),2)))*(exp(-mn[uPar]*rad)-exp(-lambda*rad)*
                        ((pow(lambda,2)-pow(mn[uPar],2))*rad+2*lambda)/(2*lambda)));
        }
    return Result;
}

//flag 0 = p
//flag 1 = n
//flag 2 = p+avgIsospinPotential (i.e. for only a single channel and V=0.5*(VI0+VI1))
//flag 3 = p+avgIsospinPotential (i.e. for only a single channel and V=0.5*(VI0+VI1))
double Tetsuo_pKm(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius){

    const double KI0[11]={-10.833,-2.7962,-0.47980,-0.64480,0.44645,0.089658,-0.23222,0.027650,0.059123,-0.024071,0.0022208};
    const double KI1[11]={-6.2261,-2.1909,-0.37668,-0.14782,2.9791,-0.53283,-3.0760,1.5430,0.64668,-0.59141,0.10746};
    const double sqrtPI = sqrt(3.1415926535897932384626433832795028841971693993751);
    const double hbarc = 197.3269602;

    const double bI0 = 0.38;
    const double bI1 = 0.37;

    const double MassN = DlmPotFlag%2==0?938.272:939.565;
    const double MassKc = 493.677;
    const double RedMass = (MassN*MassKc)/(MassN+MassKc);

    const double& Momentum = Radius[1];
    const double Energy = Momentum*Momentum/(2.*RedMass);
    const double sqrtS = Energy+MassN+MassKc;
    const double EnergyN = (sqrtS*sqrtS-MassKc*MassKc+MassN*MassN)/(2.*sqrtS);
    const double OmegaK = (sqrtS*sqrtS+MassKc*MassKc-MassN*MassN)/(2.*sqrtS);

    double Result=0;

    if(DlmPotFlag==2||DlmPotFlag==3){
        Result = 0.5*(Tetsuo_pKm(WhichPot,DlmPotFlag-2,0,t2p1,t2p2,Spin,AngMom,TotMom,Radius)+Tetsuo_pKm(WhichPot,DlmPotFlag-2,1,t2p1,t2p2,Spin,AngMom,TotMom,Radius));
    }
    else if(IsoSpin==0){
        for(unsigned uPar=0; uPar<11; uPar++){
            Result += KI0[uPar]*pow(Energy*0.01,uPar);
        }
        Result *= ((OmegaK+EnergyN)*MassN*exp(-pow(Radius[0]/bI0,2.)))/(OmegaK*EnergyN*2.*(Energy+MassN+MassKc)*pow(sqrtPI*bI0,3.));
        Result *= hbarc*hbarc;
    }
    else if(IsoSpin==1){
        for(unsigned uPar=0; uPar<11; uPar++){
            Result += KI1[uPar]*pow(Energy*0.01,uPar);
        }
        Result *= ((OmegaK+EnergyN)*MassN*exp(-pow(Radius[0]/bI1,2.)))/(OmegaK*EnergyN*2.*(Energy+MassN+MassKc)*pow(sqrtPI*bI1,3.));
        Result *= hbarc*hbarc;
    }
    else{
        Result = 0;
    }
/*
if(AngMom!=-1){
double TempI0;
double TempI1;
double TempIA;
TempI0 = Tetsuo_pKm(WhichPot,0,0,t2p1,t2p2,Spin,-1,TotMom,Radius);
TempI1 = Tetsuo_pKm(WhichPot,0,1,t2p1,t2p2,Spin,-1,TotMom,Radius);
TempIA = Tetsuo_pKm(WhichPot,2,0,t2p1,t2p2,Spin,-1,TotMom,Radius);
printf("R=%.2f\n",Radius[0]);
printf(" VI0=%.2f\n",TempI0);
printf(" VI1=%.2f\n",TempI1);
printf(" VIA=%.2f\n",TempIA);
}
*/
    return Result;
}

//I took it from Oli, he found it in some paper somewhere...
//the difference to the CRAB potential is only in the 3P1 (triplet) state
double ReidSoftCore(const int& polar, double* Radius){
    double r = Radius[0];
    double pmux,f1,f2,f3,f4,f7,vr;
    vr=0;
    /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
     with Tensor forces neglected */
    if(polar==0){
        /* S=0 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
    }
    else if(polar>0){
        /* S=1 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f2=f1*f1;
        f3=f2*f1;
        f4=f2*f2;
        //f6=f4*f2;

        //ANNALS OF PHYSICS: 50, 411448 (1968)
        //This potential is originally included in CRAB (check), but describes the 3P2 - 3F2 mixing according to Reid paper
        //vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;

        //potential according to Reid:
        vr = 10.463 * ( (1. + 2./pmux + 2./(pmux*pmux))*f1 - (8./pmux + 2./(pmux*pmux))*f4 )/pmux - 135.25 * f2/pmux + 472.81 * f3/pmux;
    }

    return vr;//MeV
}

double ReidSoftCore1S0(double* Radius){
    return ReidSoftCore(0, Radius);
}

double ReidSoftCore3P(double* Radius){
    return ReidSoftCore(1, Radius);
}

//input in fm
double fReidMeVfm(const double& rad,const unsigned& polar){
    double r = rad;//convert in fm
    double pmux,f1,f2,f4,f6,f7,vr;
    /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
     with Tensor forces neglected */
    if(polar==0){
        /* S=0 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
        //vr=0;
    }
    else if(polar>0){
        /* S=1 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f2=f1*f1;
        f4=f2*f2;
        f6=f4*f2;
        vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;
    }

    return vr;//MeV
}

double fReidMeVfm1S0(double* rad){
    return fReidMeVfm(*rad, 0);
}
double fReidMeVfm3P(double* rad){
    return fReidMeVfm(*rad, 1);
}

double fReidDlm(const double& rad,const unsigned short& s,const unsigned short& l,const unsigned short& j){
    const double mu_const = 0.7;
    const double h_const = 10.463;
    double f1,f2,f3,f4,f6,f7,vr;
    const double pmux=rad*mu_const;
    f1=exp(-pmux);
    f2=f1*f1;
    f3=f2*f1;
    f4=f2*f2;
    f6=f4*f2;
    f7=f4*f3;
    vr=0;
    if(s==0 && l==0 && j==0){
        /* 1S0 */
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-h_const*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
        //vr=0;
    }
    else if(s==1 && l==1 && j==0){
        /* 3P0 */
        vr =    -h_const*((1.+4./pmux+4./pmux/pmux)*f1 - (16./pmux+4./pmux/pmux)*f4)/pmux
                + 27.133*f2/pmux - 790.74*f4/pmux + 20662.*f7/pmux;
    }
    else if(s==1 && l==1 && j==1){
        /* 3P1 */
        vr =    h_const*((1.+2./pmux+2./pmux/pmux)*f1 - (8./pmux+2./pmux/pmux)*f4)/pmux
                - 135.25*f2/pmux + 472.81*f3/pmux;
    }
    else if(s==1 && l==1 && j==2){
        /* 3P2 - 3F2 */
        vr =    (h_const*f1/3. - 933.48*f4 + 4152.1*f6)/pmux;
    }
    return vr;//MeV
}

double fReidDlm1S0(double* rad){
    return fReidDlm(*rad,0,0,0);
}
double fReidDlm3P0(double* rad){
    return fReidDlm(*rad,1,1,0);
}
double fReidDlm3P1(double* rad){
    return fReidDlm(*rad,1,1,1);
}
double fReidDlm3P2(double* rad){
    return fReidDlm(*rad,1,1,2);
}
double fReidDlm3P(double* rad){
    return (fReidDlm(*rad,1,1,0)+fReidDlm(*rad,1,1,1)+fReidDlm(*rad,1,1,2))/3.;
}

double fReidVale(const double& rad,const unsigned short& Spin,const unsigned short& AngMom,const unsigned short& TotMom){
    const double mu_const = 0.7;
    const double h_const = 10.463;
    double f1,f2,f3,f4,f6,f7,vc,vt,vls;
    const double pmux=rad*mu_const;
    f1=exp(-pmux);
    f2=f1*f1;
    f3=f2*f1;
    f4=f2*f2;
    f6=f4*f2;
    f7=f4*f3;

    if(Spin==0 && AngMom==0 && TotMom==0){
        /* 1S0 */
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vc=-h_const*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
        vt=0;
        vls=0;
        //vc=0;
    }
    else if(Spin==1 && AngMom==1 && TotMom==0){
        /* 3P0 */
        vc =    -h_const*((1.+4./pmux+4./pmux/pmux)*f1 - (16./pmux+4./pmux/pmux)*f4)/pmux
                + 27.133*f2/pmux - 790.74*f4/pmux + 20662.*f7/pmux;
        vt=0;
        vls=0;
    }
    else if(Spin==1 && AngMom==1 && TotMom==1){
        /* 3P1 */
        vc =    h_const*((1.+2./pmux+2./pmux/pmux)*f1 - (8./pmux+2./pmux/pmux)*f4)/pmux
                - 135.25*f2/pmux + 472.81*f3/pmux;
        vt=0;
        vls=0;
    }
    else if(Spin==1 && AngMom==1 && TotMom==2){
        /* 3P2 - 3F2 */
        vc =    (h_const*f1/3. - 933.48*f4 + 4152.1*f6)/pmux;
        vt=h_const*((1./3.+1./pmux+1./pmux/pmux)*f1-(4./pmux+1./pmux/pmux)*f4)/pmux-34.925*f3/pmux;
        vls=-2074.1*f6/pmux;
    }
    else{
        vc=0; vt=0; vls=0;
    }

    double s12=0;
    double s12m=0;
    //double s12p;
    const double lsm=TotMom-1;
    //const double lsp=-(TotMom+2);
    const int ls=(TotMom*(TotMom+1)-AngMom*(AngMom+1)-Spin*(Spin+1))/2;

    if(Spin==1 && AngMom==TotMom) s12=2;
    else if(AngMom==(TotMom+1)) {s12=-2.*double(TotMom+2)/double(2*TotMom+1);}
    else if(Spin==1 && AngMom==(TotMom-1)){
        s12m=-2.*double(TotMom-1.)/double(2.*TotMom+1.);
        s12=sqrt(double(36.*TotMom*(TotMom+1)))/double(2.*TotMom+1.);
        //s12p=-2.*double(TotMom+2.)/double(2.*TotMom+1.);
    }

    if(Spin==1 && AngMom==(TotMom-1)) return vc+s12m*vt+lsm*vls;
    else return vc+s12*vt+ls*vls;
}


double fReidVale1S0(double* rad){
    return fReidVale(*rad,0,0,0);
}
double fReidVale3P0(double* rad){
    return fReidVale(*rad,1,1,0);
}
double fReidVale3P1(double* rad){
    return fReidVale(*rad,1,1,1);
}
double fReidVale3P2(double* rad){
    return fReidVale(*rad,1,1,2);
}
double fReidVale3P(double* rad){
    return (fReidVale(*rad,1,1,0)+fReidVale(*rad,1,1,1)+fReidVale(*rad,1,1,2))/3.;
}

//N.B. the potentials are never deleted, i.e. they stay there until termination
double fV18potential(const int& V18Pot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius){

    if( (V18Pot<1||V18Pot>24) && (V18Pot<112||V18Pot>114) && (V18Pot<122||V18Pot>124) ){
        return 0;
    }

    #pragma omp critical
    {
    if(!fV18pot){
        NumThreads_DLMPOT = 1;
        NumThreads_DLMPOT = omp_get_num_procs();
        fV18pot = new DLM_StefanoPotentials** [NumThreads_DLMPOT];
        for(unsigned uThread=0; uThread<NumThreads_DLMPOT; uThread++){
            fV18pot[uThread] = new DLM_StefanoPotentials* [30];
            for(unsigned uPot=0; uPot<30; uPot++){
                fV18pot[uThread][uPot] = NULL;
            }
        }
    }
    }


    unsigned StefPotId = V18Pot-1;
    if(V18Pot==112) StefPotId=24;
    else if(V18Pot==113) StefPotId=25;
    else if(V18Pot==114) StefPotId=26;
    else if(V18Pot==122) StefPotId=27;
    else if(V18Pot==123) StefPotId=28;
    else if(V18Pot==124) StefPotId=29;
    unsigned tid = 0;
    tid = omp_get_thread_num();
    if(!fV18pot[tid][StefPotId]){
        fV18pot[tid][StefPotId] = new DLM_StefanoPotentials(V18Pot);
    }

//printf("StefPotId=%i --> V=%f\n",StefPotId,fV18pot[StefPotId]->EvalCATS_v1_0(Radius[0],0));
    //return fV18pot[StefPotId]->Eval_PWprojector_pp(Radius[0],Spin,AngMom,TotMom,DlmPotFlag);
//printf("fV18pot[%u][%u]=%p\n",tid,StefPotId,fV18pot[tid][StefPotId]);
    return fV18pot[tid][StefPotId]->Eval_PWprojector(Radius[0],IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,DlmPotFlag);
    //return 0;
}



////////////////////////////////
//! pLambda

//(this is Oliver's version, not CRAB!)
double UsmaniPotentialOli(const int& Spin, double* Radius)
{

  double r = Radius[0];
  //Values for the potential
  const double vbar = 6.2;

  const double vsigma = 0.25;

  const double wc = 2137;

  double x=r*0.7;
  double vc = wc/(1+exp((r-0.5)/0.2));
  double tpi = (1.0+3.0/x+3.0/(x*x)) * (exp(-x)/x) * pow(1.-exp(-2.*r*r),2.);

  double v = 0.;

  if (Spin == 0) v = vc - (vbar + 0.75*vsigma)*tpi*tpi;//Usmani singlet
  else if (Spin == 1)  v = vc - (vbar - 0.25*vsigma)*tpi*tpi;//Usmani triplet
  else printf ("wrong polarization\n");

  return v;

}


////////////////////////////////////////////


double ppDlmPot(const int& DlmPot, const int& DlmFlag, const int& Spin, const int& AngMom, const int& TotMom, double* Radius){
    return fDlmPot(DlmPot,DlmFlag,1,1,1,Spin,AngMom,TotMom,Radius);
}

//[2] = DlmPot, [3] = DlmFlag
double ppDlmPot1S0(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],0,0,0,Pars);
}
double ppDlmPot3S1(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,0,1,Pars);
}
double ppDlmPot3P0(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,1,0,Pars);
}
double ppDlmPot3P1(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,1,1,Pars);
}
double ppDlmPot3P2(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,1,2,Pars);
}
double ppDlmPot3P(double* Pars){
    return (ppDlmPot3P0(Pars)+ppDlmPot3P1(Pars)+ppDlmPot3P2(Pars))/3.;
}

double pLambdaDlmPot(const int& DlmPot, const int& DlmFlag, const int& Spin, const int& AngMom, const int& TotMom, double* Radius){
    return fDlmPot(DlmPot,DlmFlag,0,0,0,Spin,AngMom,TotMom,Radius);
}
//[2] is Pot Type, 3 is [3] is PotFlag
double pLambdaDlmPot1S0(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],0,0,0,Pars);
}
double pLambdaDlmPot3S1(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,0,1,Pars);
}

//t2p1 - 2xIsospin of particle 1, t2p2 same for particle 2
double fDlmPot(const int& DlmPot, const int& DlmPotFlag,
               const int& IsoSpin, const int& t2p1, const int& t2p2, const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){
    //printf("V=%f\n",fV18potential(9,Spin,AngMom,TotMom,Radius)) ;
    switch(DlmPot){
        case NN_AV18 : return fV18potential(9,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius);
        case NN_ReidV8 : return fV18potential(2,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius);
        case pp_ReidSC : return fReidDlm(Radius[0],Spin,AngMom,TotMom);
        case pp_ReidOli : return ReidSoftCore(Spin,Radius);
        case pp_ReidCrab : return fReidMeVfm(Radius[0],Spin);
        case pp_ReidVale : return fReidVale(Radius[0],Spin,AngMom,TotMom);
        case pL_UsmaniOli : return UsmaniPotentialOli(Spin,Radius);
        case pXim_Lattice : return LatticePots_pXi(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pXim_HALQCD1 : return LatticePots_pXi_ver2(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pXim_LatticeAvg : return LatticePots_pXi_Avg(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pXim_LatticeSqrtAvg : return LatticePots_pXi_SqrtAvg(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pKm_Tetsuo : return Tetsuo_pKm(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius);
        case pOmega_Lattice : return LatticePots_pOmega(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pOmega_Tetsuo : return Tetsuo_pOmega(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        default : return 0;
    }
}

//[0] radius, [1] momentum
//[2] PotentialType, [3] PotentialFlag,
//[4] IsoSpin, [5] ParticleType1 (1 proton, -1 neutron), [6] ParticleType2
//[7] Spin (s), [8] AngMom (l), [9] TotMom (j)
//[10] - optional stuff
double fDlmPot(double* Parameters){
    //printf(" fDlmPot called with %p\n",Parameters);
    return fDlmPot(round(Parameters[2]),round(Parameters[3]),round(Parameters[4]),round(Parameters[5]),
                   round(Parameters[6]),round(Parameters[7]),round(Parameters[8]),round(Parameters[9]),Parameters,&Parameters[10]);
}


void GetDlmPotName(const int& potid, const int& potflag, char* name){
    double Radius=1;
    //fV18potential(1,0,0,0,&Radius);
    fV18potential(1,0,1,1,1,0,0,0,&Radius);
    unsigned tid = 0;
    tid = omp_get_thread_num();
    switch(potid){
        case NN_AV18 :
            fV18pot[tid][8]->PotentialName(9, name);
            break;
        case NN_ReidV8 :
            fV18pot[tid][1]->PotentialName(2, name);
            break;
        case pp_ReidSC :
            strcpy(name,"Castrated Reid SC");
            break;
        case pp_ReidOli :
            strcpy(name,"Reid 3P2-3F2");
            break;
        case pp_ReidCrab :
            strcpy(name,"Reid 3P1");
            break;
        case pp_ReidVale :
            strcpy(name,"Reid Soft-Core");
            break;
        default :
            strcpy(name,"Unknown potential");
            break;
    }
    char Buffer[16];
    sprintf(Buffer, "%i", potflag);
    if(potflag!=0 && strcmp(name,"Unknown potential")){strcat(name, "^{(");strcat(name,Buffer);strcat(name,")}");}
}
