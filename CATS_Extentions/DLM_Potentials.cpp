
#include "DLM_Potentials.h"
#include "DLM_StefanoPotentials.h"
#include "DLM_Histo.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <omp.h>

//int DlmPot=0;
//int DlmPotFlag=0;
DLM_StefanoPotentials*** fV18pot=NULL;
//[lpot][lemp][slj]
DLM_Histo<float>*** fNorfolkPot_pp=NULL;

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

double GaussExpSum(double* Pars){
    return Pars[2]*exp(-pow(Pars[0]/Pars[3],2))+Pars[4]*exp(-pow(Pars[0]/Pars[5],1));
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
            return fDlmPot(NN_AV18,v18_Coupled3P2,1,1,1,1,1,1,&ShiftedRad,0);
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
            return fDlmPot(NN_AV18,v18_Coupled3P2,1,1,1,1,1,1,&ShiftedRad,0);
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

    //if element is 3 => correct potential
    //if element is 4 => the old (wrong) parameterization of the yukawa term (with power of 2)
    double EvalV(const int& iFlag, const int& Element, const double& Rad){
        if(iFlag>=3 || iFlag<-3) return 0;
        double YukawaPowerVST=2.0;
        int ELEMENT = Element;
        if(Element==3) YukawaPowerVST=1.0;//for the 3rd element only
        if(Element==4) ELEMENT=3;//dummy flag to make the 3rd element wrong (as in an older computation)
        double* PAR = ELEMENT==0?PAR_V0[iFlag]:ELEMENT==1?PAR_VS[iFlag]:ELEMENT==2?PAR_VT[iFlag]:ELEMENT==3?PAR_VST[iFlag]:NULL;
        if(!PAR) return 0;
        if(!Rad || !PAR[1] || !PAR[3] || !PAR[5]) return 0;
        return  PAR[0]*exp(-pow(Rad/PAR[1],2.))+PAR[2]*exp(-pow(Rad/PAR[3],2.))+PAR[4]*exp(-pow(Rad/PAR[5],2.))+
                PAR[6]*pow((1.-exp(-PAR[7]*Rad*Rad))*exp(-146./197.3269602*Rad)/Rad,YukawaPowerVST);
    }

    double Eval(const int& DlmPotFlag, const double& IsoSpin, const double& Spin, const double& Rad){
        int iFlag = abs(DlmPotFlag)-11;
        if(IsoSpin==0 && Spin==0){
            return EvalV(iFlag,0,Rad)-3.*EvalV(iFlag,1,Rad)-3.*EvalV(iFlag,2,Rad)+9.*EvalV(iFlag,3+(DlmPotFlag<0),Rad);
        }
        else if(IsoSpin==0 && Spin==1){
            return EvalV(iFlag,0,Rad)+EvalV(iFlag,1,Rad)-3.*EvalV(iFlag,2,Rad)-3.*EvalV(iFlag,3+(DlmPotFlag<0),Rad);
        }
        else if(IsoSpin==1 && Spin==0){
            return EvalV(iFlag,0,Rad)-3.*EvalV(iFlag,1,Rad)+EvalV(iFlag,2,Rad)-3.*EvalV(iFlag,3+(DlmPotFlag<0),Rad);
        }
        else if(IsoSpin==1 && Spin==1){
            return EvalV(iFlag,0,Rad)+EvalV(iFlag,1,Rad)+EvalV(iFlag,2,Rad)+EvalV(iFlag,3,Rad);
        }
        else return 0;
    }
};


struct LatticeValuesPaper{
  std::vector<std::vector<float>> parameter;
  std::vector<std::vector<float>> parameterOLD; 
  LatticeValuesPaper() {
    //Following the most recent hal publication: https://inspirehep.net/literature/1771619
    //Values are summarized: https://docs.google.com/spreadsheets/d/12g3YztBIpL7BlTBpEbNnV_C5vR8K1NiOLOWwRoZb03I/edit?usp=sharing
    //Scheme:
    //          beta_1(0)  , beta_2(1)  , beta_3(2)  , (shared among potentials) 
    //V0:       alpha_1(3) , alpha_2(4) , alpha_3(5) , lmb_2(6) , rho_2(7) ,
    //VSigma:   alpha_1(8) , alpha_2(9) , alpha_3(10),
    //VTau:     alpha_1(11), alpha_2(12), alpha_3(13),
    //VSigmaTau:alpha_1(14), alpha_2(15), alpha_3(16), lmb_1(17), rho_1(18), 
    parameter = {
      //t = 11,  0 - 22 
      {0.1293,0.2574,0.5671,956.5365,551.3278,173.5383,-108.3599,0.6054,-124.1416,-51.0465,-5.4057,192.5351,103.4974,32.1514,-79.912,-37.3872,-7.1677,-1.5995,0.2594},
      {0.1298,0.2597,0.5742,962.4006,551.7967,164.5277,-109.8527,0.6142,-129.9177,-47.4349,-5.7077,195.0046,103.4218,31.0253,-79.8773,-37.0849,-6.9693,-1.5764,0.247},
      {0.1285,0.2548,0.5656,947.132,559.6762,175.2138,-110.3494,0.6105,-124.3225,-51.3004,-6.5205,188.5444,107.7063,31.6191,-80.6019,-37.3876,-6.9862,-1.694,0.2419},
      {0.1281,0.2538,0.5649,941.1004,566.0914,175.7539,-108.3571,0.6071,-122.2735,-52.5149,-6.0959,189.5948,106.987,32.6486,-79.1654,-38.2752,-7.0508,-1.6424,0.2424},
      {0.1286,0.2548,0.5651,944.594,559.7067,177.642,-107.804,0.6024,-124.1264,-51.6175,-5.7903,191.6881,104.8241,32.2946,-79.0914,-37.684,-7.1521,-1.6885,0.2458},
      {0.1284,0.254,0.5633,943.5473,560.7421,178.5771,-108.1242,0.6035,-124.5062,-52.2075,-5.9533,188.337,108.6364,32.1538,-78.8977,-38.0209,-7.499,-1.5837,0.2512},
      {0.1304,0.2626,0.5788,972.7943,538.2606,163.5877,-113.9894,0.6191,-126.1337,-48.3478,-5.4429,197.2591,102.2528,30.7907,-79.5489,-37.6124,-6.5352,-1.6452,0.2495},
      {0.1295,0.2583,0.5686,959.7033,549.2876,171.3685,-109.6581,0.6093,-124.78,-51.3518,-5.2487,191.9771,104.0978,31.5561,-80.4259,-37.0074,-6.9272,-1.655,0.248},
      {0.1299,0.2599,0.5718,966.3629,547.2302,167.0596,-109.7925,0.6115,-125.4904,-50.0556,-5.4705,192.9566,103.6372,30.8222,-79.7209,-37.2283,-6.7919,-1.6256,0.2473},
      {0.1298,0.2599,0.571,963.845,551.4802,166.4173,-107.4425,0.6091,-127.2933,-49.0915,-5.1697,194.1605,102.8876,31.4879,-80.3065,-37.2837,-6.8337,-1.6361,0.2564},
      {0.1297,0.2589,0.5725,961.7472,548.2013,169.1965,-111.6276,0.6117,-126.4412,-50.0892,-5.3419,194.3707,103.273,31.1054,-79.5665,-37.9146,-7.1,-1.5571,0.2501},
      {0.1294,0.258,0.5713,956.4063,553.647,169.6398,-110.686,0.6104,-127.75,-48.3829,-5.6038,191.7904,105.5012,30.9239,-79.3499,-37.6533,-7.0623,-1.621,0.2452},
      {0.1282,0.2531,0.5594,937.9336,558.5323,182.7377,-106.5613,0.5968,-124.384,-51.7904,-5.9424,190.8192,106.4756,32.9703,-78.9275,-39.3348,-7.0515,-1.6618,0.2413},
      {0.1293,0.2581,0.5688,959.1115,552.4117,171.7665,-109.4964,0.6087,-125.5572,-50.7739,-5.5345,194.2351,103.1844,31.8213,-80.0427,-37.6954,-7.064,-1.6312,0.2446},
      {0.1299,0.2604,0.5747,967.2612,548.7219,165.8063,-112.0166,0.6171,-124.6742,-50.1813,-5.5088,193.6907,103.7304,30.5766,-79.4093,-37.235,-6.8329,-1.5997,0.2463},
      {0.13,0.2602,0.5743,965.5491,546.9459,166.4396,-112.4856,0.6153,-127.4257,-49.0366,-5.1588,195.3148,102.6436,31.2775,-79.4646,-37.1564,-6.8706,-1.5455,0.2468},
      {0.1296,0.2583,0.5685,958.0645,548.7949,172.7338,-108.8329,0.6065,-127.4739,-48.9248,-5.8319,193.2001,103.6846,31.5587,-80.1944,-37.1187,-6.7874,-1.6671,0.2508},
      {0.1289,0.2567,0.5649,954.0203,555.0487,174.5725,-109.354,0.6091,-124.4186,-51.1751,-5.5766,191.439,105.4763,32.1902,-78.7879,-38.3427,-6.9149,-1.6105,0.2422},
      {0.1288,0.2559,0.5647,950.7509,556.8509,175.199,-109.3326,0.609,-122.6436,-52.1209,-5.5735,190.6165,105.6,32.3319,-79.081,-37.4545,-7.1825,-1.6177,0.2447},
      {0.1301,0.2607,0.5721,971.8803,543.4366,165.95,-110.589,0.6132,-127.5719,-49.3792,-5.0341,191.841,104.9342,30.7819,-80.6547,-36.5441,-6.7264,-1.6429,0.2514},
      {0.13,0.2604,0.5753,968.628,549.5837,164.7968,-110.4965,0.6137,-126.8373,-49.4348,-5.8063,194.153,101.7749,30.6115,-80.077,-37.0154,-6.7615,-1.6295,0.2495},
      {0.1296,0.2584,0.5688,960.4571,543.9938,174.5126,-110.7417,0.6064,-126.4445,-50.0331,-5.3404,193.1581,105.0561,31.2893,-80.1471,-37.0134,-7.1054,-1.6197,0.2485},
      {0.1292,0.257,0.5694,955.3817,553.8313,172.5263,-109.9073,0.6069,-125.9599,-51.9737,-5.2443,191.3695,104.445,31.6577,-78.6583,-38.4901,-7.2101,-1.5864,0.2814},
      //t = 12, 23 - 45
      {0.1243,0.2399,0.5323,869.2726,626.486,199.3267,-96.6568,0.5932,-121.9414,-53.4535,-7.9636,183.8645,108.8791,36.7804,-84.4872,-31.5858,-12.0533,-1.3531,0.1295},
      {0.1244,0.2412,0.5346,869.3825,633.3597,188.9363,-96.3453,0.6076,-126.2292,-48.3133,-9.0465,187.412,108.4225,37.2313,-85.2586,-32.1331,-11.8374,-1.3908,0.1352},
      {0.1237,0.2377,0.5325,861.7193,638.0694,197.8498,-97.5622,0.6009,-121.3455,-52.6356,-9.4646,182.4489,109.8547,37.2802,-85.8937,-31.6917,-11.668,-1.4942,0.132},
      {0.1229,0.2358,0.5303,845.373,644.7708,203.9443,-95.6018,0.5883,-118.0282,-55.1575,-8.607,182.6999,110.3945,38.0576,-83.4625,-33.7146,-11.7095,-1.3993,0.1333},
      {0.1241,0.2395,0.5337,864.7638,632.5374,197.5935,-95.8526,0.5933,-123.1762,-51.892,-8.9078,184.2065,109.2794,36.6019,-85.7123,-31.4519,-12.0585,-1.4408,0.1314},
      {0.1246,0.2421,0.5316,876.9214,634.9737,185.861,-98.2123,0.6226,-120.6757,-53.7512,-8.4423,184.0762,109.3695,36.4222,-85.3797,-31.8685,-11.8205,-1.3711,0.1318},
      {0.1243,0.2401,0.5403,865.0805,624.295,199.6041,-102.1491,0.5967,-119.3229,-51.7812,-9.1551,184.6361,109.8487,36.5448,-81.3888,-35.1294,-11.0715,-1.4473,0.144},
      {0.1246,0.2418,0.5304,874.5434,631.6375,189.6598,-97.3021,0.615,-121.1861,-52.9023,-7.7991,183.9227,109.2174,36.2288,-86.1726,-30.9641,-11.8553,-1.4055,0.13},
      {0.1244,0.241,0.5333,874.9475,626.0273,196.5466,-96.1963,0.5933,-119.8352,-54.1296,-8.0826,186.0924,107.6127,36.9874,-84.6364,-32.7701,-11.4199,-1.4162,0.1378},
      {0.1246,0.2428,0.5343,877.5143,629.5967,190.1578,-95.0588,0.5998,-123.0794,-51.11,-7.4736,184.7507,107.2605,36.4891,-85.2892,-31.7367,-11.4275,-1.4147,0.1361},
      {0.1242,0.24,0.5361,869.4685,633.0055,194.0039,-99.4875,0.6051,-122.5524,-52.7285,-8.4179,186.8088,108.6526,36.1353,-84.6383,-32.6335,-12.0299,-1.3523,0.1379},
      {0.1253,0.2436,0.5382,884.7147,624.9132,182.6392,-98.8174,0.6174,-122.824,-48.5951,-8.1418,181.1301,110.9874,34.0988,-85.7437,-31.0667,-11.8853,-1.3576,0.1307},
      {0.1227,0.234,0.52,841.7348,644.056,211.4243,-94.5731,0.5921,-118.9633,-56.7354,-8.7035,179.2442,111.4349,39.3407,-82.0135,-35.0272,-11.4835,-1.5067,0.147},
      {0.1245,0.24,0.534,868.8321,624.0716,203.3387,-96.6548,0.5839,-123.5238,-51.0667,-9.1783,186.7109,108.0361,37.6815,-85.7948,-30.5178,-12.0381,-1.4759,0.1337},
      {0.1252,0.2436,0.5409,881.9956,625.6024,185.5533,-100.5601,0.6133,-120.8483,-52.3927,-8.2416,186.2616,108.6618,35.0782,-85.0458,-32.3368,-11.1449,-1.4196,0.1399},
      {0.1247,0.2419,0.534,875.2891,626.7128,189.9732,-99.3135,0.6119,-121.7041,-51.8007,-8.0507,186.0056,107.6007,36.8232,-83.7826,-32.7229,-11.7083,-1.4033,0.1393},
      {0.1244,0.241,0.5288,870.1969,628.7571,195.4659,-93.138,0.5962,-125.3333,-50.8744,-8.6994,184.3621,109.3683,36.7679,-86.0777,-31.4534,-11.1469,-1.5155,0.1356},
      {0.1248,0.2418,0.5319,875.0111,630.6641,190.6844,-98.491,0.6161,-122.8088,-49.8425,-8.3679,187.496,106.5172,37.4021,-85.5914,-32.6351,-11.2962,-1.3305,0.1392},
      {0.1235,0.237,0.5233,858.0896,638.9581,202.7689,-94.1689,0.5972,-120.1618,-53.8632,-8.9756,186.3608,105.337,38.8104,-85.4883,-30.8076,-11.9578,-1.4767,0.1316},
      {0.1254,0.2438,0.5309,891.2166,620.5545,185.0245,-97.3386,0.6149,-122.9415,-52.4436,-7.3945,185.4309,107.4885,36.1271,-85.467,-31.5966,-11.0324,-1.5022,0.1392},
      {0.1246,0.2416,0.5351,878.3103,628.7247,192.5959,-97.6808,0.6018,-122.3331,-51.3059,-8.7188,185.1585,107.08,35.7133,-84.279,-32.7481,-11.3131,-1.4102,0.1396},
      {0.1254,0.2433,0.5332,885.9252,614.7054,192.9299,-97.9995,0.6047,-123.7995,-52.104,-7.3898,188.5747,105.4994,36.3052,-84.8647,-32.3368,-11.2899,-1.4088,0.144},
      {0.1254,0.2435,0.5404,882.8502,623.058,187.5172,-97.8534,0.5994,-124.395,-50.6348,-8.576,186.503,107.7003,35.3144,-85.9443,-31.6568,-11.3147,-1.3851,0.1363},
      //t = 13 46 - 67
      {0.1256,0.2324,0.4998,861.6034,458.9294,369.4892,-83.3086,0.4118,-133.2605,-45.4614,-12.8217,194.0236,91.9739,44.2243,-68.9582,-41.7269,-16.3794,-1.4048,0.2994},
      {0.1218,0.2219,0.5046,807.4034,553.1615,335.3227,-82.3407,0.4365,-122.1551,-58.2538,-11.5397,179.4353,103.1913,44.3769,-64.7202,-47.8337,-16.1288,-1.334,0.2562},
      {0.1228,0.2233,0.5016,818.2958,534.2739,343.9809,-82.0568,0.4303,-119.8854,-56.6616,-11.4019,190.9593,95.316,43.965,-63.1268,-47.2788,-15.4735,-1.3768,0.2521},
      {0.1242,0.2288,0.5033,840.3624,501.7288,349.7359,-82.8958,0.421,-132.3241,-48.909,-13.0284,189.1946,94.9626,43.7524,-66.6964,-45.4421,-16.0307,-1.3946,0.3117},
      {0.1233,0.2254,0.4935,831.9987,487.5463,381.2269,-81.9759,0.4087,-122.3417,-56.3799,-12.5472,186.9325,95.4097,45.4802,-64.6997,-46.1564,-16.8149,-1.369,0.3085},
      {0.1218,0.2188,0.5061,791.6493,557.6506,348.4973,-88.554,0.4416,-121.2703,-57.9451,-12.5107,175.646,105.8907,43.35,-58.7371,-51.0519,-15.5238,-1.3659,0.256},
      {0.1231,0.2247,0.4882,828.4536,465.786,408.1914,-81.8647,0.3936,-124.952,-56.2897,-10.7039,186.6398,96.3516,45.9736,-62.5022,-45.8057,-17.39,-1.3283,0.3679},
      {0.1252,0.2329,0.5018,869.7029,489.0783,341.6121,-82.217,0.4286,-118.8216,-57.3441,-10.4072,196.2683,87.0045,45.5676,-66.1902,-44.2073,-14.9444,-1.4226,0.2681},
      {0.1258,0.2336,0.5022,868.2639,486.1471,340.6827,-80.6898,0.4246,-139.2242,-42.0178,-12.0866,196.6551,85.5578,44.6747,-69.2321,-41.3312,-15.5511,-1.4049,0.2878},
      {0.1237,0.2273,0.5011,841.3998,463.8176,392.496,-87.7526,0.4039,-125.4859,-51.5305,-13.6465,191.601,94.479,44.9427,-68.0537,-43.6687,-16.8054,-1.3306,0.3131},
      {0.1243,0.2277,0.5017,841.7059,485.6206,362.5004,-82.7994,0.4129,-124.2341,-50.074,-12.7079,181.218,98.7631,41.573,-68.6742,-42.9417,-17.4121,-1.304,0.3718},
      {0.1227,0.223,0.4864,823.8422,453.1411,426.8382,-83.1302,0.3857,-122.6792,-56.8454,-13.1594,180.2531,99.0796,46.6662,-61.822,-46.4896,-15.9296,-1.5059,0.2595},
      {0.1219,0.2232,0.5065,813.8276,574.017,317.8289,-83.7579,0.4554,-128.5921,-53.6106,-13.5196,171.0542,112.1573,44.7662,-64.0299,-43.2541,-15.9559,-1.4417,0.2738},
      {0.126,0.2344,0.5093,868.3832,485.1417,332.8766,-84.8049,0.4342,-122.3385,-53.7416,-12.102,194.7484,89.5973,42.9303,-66.2946,-44.0512,-15.1793,-1.3614,0.2877},
      {0.121,0.2177,0.4928,799.085,501.7868,399.4776,-86.2351,0.4071,-121.3964,-57.0124,-13.1324,176.6735,106.6711,45.4701,-58.3644,-50.3742,-18.3112,-1.3031,0.3126},
      {0.1252,0.2337,0.4917,859.5293,462.6938,371.397,-79.2619,0.4045,-127.6561,-51.617,-12.0201,195.3939,88.0511,45.6385,-66.0585,-45.6921,-16.8225,-1.4551,0.4759},
      {0.1258,0.231,0.4991,862.3239,437.7641,394.4049,-87.8782,0.404,-128.1008,-46.9421,-12.2101,193.6237,89.564,45.4565,-68.8651,-42.7661,-15.0458,-1.3315,0.2903},
      {0.1201,0.213,0.4867,757.8726,582.0278,363.0112,-76.6633,0.4231,-119.6642,-56.7003,-13.6544,177.2864,103.6267,48.7892,-60.312,-49.455,-16.9119,-1.4333,0.3043},
      {0.1223,0.2204,0.4902,814.1627,508.2253,370.925,-78.3432,0.4127,-114.8792,-61.787,-11.7645,183.4898,100.6586,47.6946,-59.7193,-51.1181,-18.0322,-1.4734,0.4774},
      {0.1257,0.2317,0.497,862.3302,442.6556,387.1174,-85.5663,0.4042,-132.0333,-48.5274,-12.7652,192.7584,89.9265,44.6653,-65.85,-43.9422,-16.0204,-1.3882,0.3342},
      {0.128,0.2405,0.4999,896.4969,415.9204,376.6059,-85.9244,0.4066,-132.3181,-46.6662,-10.4764,203.886,79.3207,43.7329,-72.6295,-40.1308,-14.7768,-1.4244,0.2629},
      {0.1282,0.2444,0.513,906.2299,435.2305,347.2022,-89.1028,0.4197,-131.5798,-45.768,-11.5168,202.3761,85.1461,39.9807,-77.0847,-35.737,-14.441,-1.2891,0.1796}
    };
    
    //old parameters for QA
    //V0 alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, lambda_2 = d1, rho2 = 1/sqrt(d2) 
    //VS alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, 0, 0
    //VT alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, 0, 0
    //VST alpha_1, beta_1, alpha_2, beta_2, alpha_3, beta_3, lambda_1 = d1, rho2 = 1/sqrt(d2) 
    parameterOLD =
      {
	{ 832.719,0.126567,306.315,0.262349,521.285,0.461616,-80.9157,9.41638,-112.713,0.119882,-60.3916,0.219904,-12.971,0.440375,0,0,205.93,0.135111,93.3825,0.278859,26.9143,0.588477,0,0,-79.774,0.135531,-32.564,0.275606,-9.3541,0.538997,-1.75591,2.88912 },
	{800.944,0.125499,340.209,0.25353,528.537,0.453636,-74.9729,9.89426,-124.804,0.122483,-48.2802,0.234785,-12.4604,0.450708,0,0,170.527,0.121487,117.874,0.224199,42.0892,0.510129,0,0,-90.3105,0.148112,-26.359,0.351332,-4.72814,0.707443,-1.26899,2.55046},
	{807.648,0.12642,359.403,0.255362,500.769,0.451934,-67.7968,10.0506,-52.9071,0.0916017,-119.322,0.164464,-22.4783,0.408204,0,0,124.788,0.112992,151.176,0.194521,48.0954,0.490352,0,0,-82.1509,0.15925,-26.2495,0.358359,-6.31347,0.639725,-1.32727,2.67235}
      };
  }
  
  ~LatticeValuesPaper(){} 

  double EvalVOld(const int& iFlag, const int& Element, const double& Rad){
    if(iFlag > 3 || iFlag < 0) {
      return 0;
    } else  {
      double out = 0;
      if (Element == 0) {
	double alph_1 = parameterOLD[iFlag][0];
	double beta_1 = parameterOLD[iFlag][1];
	double alph_2 = parameterOLD[iFlag][2];
	double beta_2 = parameterOLD[iFlag][3];
	double alph_3 = parameterOLD[iFlag][4];
	double beta_3 = parameterOLD[iFlag][5];
	double lmb_2  = parameterOLD[iFlag][6];
	double rho_2  = parameterOLD[iFlag][7];
	out = EvalYukawa(Rad, 1./sqrt(rho_2));
	out = out*out*lmb_2; 
	out += EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3);
	//std::cout << "Element 0, Radius: " << Rad << " EvalYukawa(Rad, 1./sqrt(rho_2))" << EvalYukawa(Rad, 1./sqrt(rho_2)) << "\n";
      } else if (Element == 1) {
	double alph_1 = parameterOLD[iFlag][8];
	double beta_1 = parameterOLD[iFlag][9];
	double alph_2 = parameterOLD[iFlag][10];
	double beta_2 = parameterOLD[iFlag][11];
	double alph_3 = parameterOLD[iFlag][12];
	double beta_3 = parameterOLD[iFlag][13];
	//14
	//15 
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3); 
      } else if (Element == 2) {
	double alph_1 = parameterOLD[iFlag][16];
	double beta_1 = parameterOLD[iFlag][17];
	double alph_2 = parameterOLD[iFlag][18];
	double beta_2 = parameterOLD[iFlag][19];
	double alph_3 = parameterOLD[iFlag][20];
	double beta_3 = parameterOLD[iFlag][21];
	//22
	//23
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3); 
      } else if (Element == 3) { 
	double alph_1 = parameterOLD[iFlag][24];
	double beta_1 = parameterOLD[iFlag][25];
	double alph_2 = parameterOLD[iFlag][26];
	double beta_2 = parameterOLD[iFlag][27];
	double alph_3 = parameterOLD[iFlag][28];
	double beta_3 = parameterOLD[iFlag][29];
	double lmb_1  = parameterOLD[iFlag][30];
	double rho_1  = parameterOLD[iFlag][31];
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3)+lmb_1*EvalYukawa(Rad, rho_1); 
      }
      return out; 
    } 
  }
  double EvalV(const int& iFlag, const int& Element, const double& Rad){
    if(iFlag > 68 || iFlag < 0) {
      return 0;
    } else  {
      double out = 0;
      double beta_1 = parameter[iFlag][0];
      double beta_2 = parameter[iFlag][1];
      double beta_3 = parameter[iFlag][2];
      if (Element == 0) {//v0
	double alph_1 = parameter[iFlag][3];
	double alph_2 = parameter[iFlag][4];
	double alph_3 = parameter[iFlag][5];
	double lmb_2  = parameter[iFlag][6];
	double rho_2  = parameter[iFlag][7];
	out = EvalYukawa(Rad, rho_2);
	out = out*out*lmb_2; 
	out += EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3); 
      } else if (Element == 1) {//vs
	double alph_1 = parameter[iFlag][8];
	double alph_2 = parameter[iFlag][9];
	double alph_3 = parameter[iFlag][10];
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3); 
      } else if (Element == 2) {//vt
	double alph_1 = parameter[iFlag][11];
	double alph_2 = parameter[iFlag][12];
	double alph_3 = parameter[iFlag][13];
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3); 
      } else if (Element == 3) {//vst
	double alph_1 = parameter[iFlag][14];
	double alph_2 = parameter[iFlag][15];
	double alph_3 = parameter[iFlag][16];
	double lmb_1  = parameter[iFlag][17];
	double rho_1  = parameter[iFlag][18];
	out = EvalExp(Rad,alph_1,beta_1)+EvalExp(Rad,alph_2,beta_2)+EvalExp(Rad,alph_3,beta_3)+lmb_1*EvalYukawa(Rad, rho_1); 
      }
      return out; 
    } 
  }
  
  double EvalExp(const double& Rad, const double& alpha, const double& beta) {
    return alpha*exp(-(Rad/beta)*(Rad/beta)); 
  }
  double EvalYukawa(const double& Rad, const double& rho) {
    double m_pion = 146./197.3269602;
    return (1.-exp(-Rad*Rad/(rho*rho)))*exp(-m_pion*Rad)/Rad; 
  }
  
  double Eval(const int& DlmPotFlag, const double& IsoSpin, const double& Spin, const double& Rad){
    int iFlag = DlmPotFlag-11;
    if(IsoSpin==0 && Spin==0){
      return EvalVOld(iFlag,0,Rad)-3.*EvalVOld(iFlag,1,Rad)-3.*EvalVOld(iFlag,2,Rad)+9.*EvalVOld(iFlag,3,Rad);
    }
    else if(IsoSpin==0 && Spin==1){
      return EvalVOld(iFlag,0,Rad)+EvalVOld(iFlag,1,Rad)-3.*EvalVOld(iFlag,2,Rad)-3.*EvalVOld(iFlag,3,Rad);
    }
    else if(IsoSpin==1 && Spin==0){
      return EvalVOld(iFlag,0,Rad)-3.*EvalVOld(iFlag,1,Rad)+EvalVOld(iFlag,2,Rad)-3.*EvalVOld(iFlag,3,Rad);
    }
    else if(IsoSpin==1 && Spin==1){
      return EvalVOld(iFlag,0,Rad)+EvalVOld(iFlag,1,Rad)+EvalVOld(iFlag,2,Rad)+EvalVOld(iFlag,3,Rad);
    }
    else return 0;
  }
};

//updated version from 1th October 2018
//there was a bug version, where one of the exponential (Yukawa) terms had the wrong power
//to call the broken version, use a negative DlmPotFlag
double LatticePots_pXi_ver2(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){
    static LatticeValues LVAL;
    if(AngMom) return 0;
    return LVAL.Eval(DlmPotFlag,IsoSpin,Spin,Radius[0]);
}

//updated version from june 2020
double LatticePots_pXi_ver3(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars){
    static LatticeValuesPaper LVAL;
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
//OtherPars[0] is a cut-off
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

//------------------------------

double KpProtonEquivalentPotential(double* pars){
    double r = pars[0];
    double& Spin = pars[2];

    double par[7]={0};
    double v = 0;
    /*
      the I=0 and I=1 equivalent potentials are taken from
      Kaon-nucleon scattering amplitudes and S' enhancements from quark Born diagrams
      T. Barnes & E. S. Swanson
      PHYSICAL REVIEW C VOLUME 49, NUMBER 2 FEBRUARY 1994
      https://journals.aps.org/prc/pdf/10.1103/PhysRevC.49.1166
      and interpolated with the sum 2 gaussian with a common normalization factor
    */
    if (Spin == 0){
      //      par[0]    =     0.931149;
      par[0]    =     1;
      par[1]    =     0.612774;
      par[2]    =   0.00379167;
      par[3]    =     0.176426;
      par[4]    =      1.43336;
      par[5]    =   0.00259035;
      par[6]    =    0.0412807;

      v = par[0]*(par[1]/par[2]*exp(-r*r/par[3])+par[4]/par[5]*exp(-r*r/par[6]));

    }
    else if (Spin == 1)  {
      //      par[0]  =     0.844204;
      par[0]  =     1.1;
      par[1]  =      1.42784;
      par[2]  =    0.0029962;
      par[3]  =     0.231457;
      par[4]  =     0.646766;
      par[5]  =   0.00299979;
      par[6]  =    0.0418019;

      v = par[0]*(par[1]/par[2]*exp(-r*r/par[3])+par[4]/par[5]*exp(-r*r/par[6]));

    }
    else printf ("wrong polarization\n");

    return v;
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

double UsmaniPotentialCats(double* Pars){
    double& r = Pars[0];
    double& Spin = Pars[2];
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

//[2] is the amplitude
//[3] is the range
//[4] is the slope
double RepulsiveCore(double* Pars){
    return Pars[2]/(1+exp((Pars[0]-Pars[3])/Pars[4]));
}

void SetUpNorfolk(const char* InputFolder){
    FILE *InFile;
    const unsigned num_lpot = 4;
    const unsigned num_lemp = 3;
    const unsigned num_slj = 5;

    char* cdummy = new char [256];
    float fRadius;
    float fPotential;

    if(fNorfolkPot_pp){
        for(unsigned upot=0; upot<num_lpot; upot++){
            for(unsigned uemp=0; uemp<num_lemp; uemp++){
                delete [] fNorfolkPot_pp[upot][uemp];
            }
            delete [] fNorfolkPot_pp[upot];
        }
        delete [] fNorfolkPot_pp;
        fNorfolkPot_pp = NULL;
    }
//printf("A\n");
    fNorfolkPot_pp = new DLM_Histo<float>** [num_lpot];
    char**** InputFileName = new char*** [num_lpot];
//printf("B\n");
    for(unsigned upot=0; upot<num_lpot; upot++){
//printf("upot=%u\n",upot);
        fNorfolkPot_pp[upot] = new DLM_Histo<float>* [num_lpot];
        InputFileName[upot] = new char** [num_lemp];
        char spot[16];
        switch(upot){
            case 0 : strcpy(spot,"105"); break;
            case 1 : strcpy(spot,"106"); break;
            case 2 : strcpy(spot,"109"); break;
            case 3 : strcpy(spot,"110"); break;
            default : break;
        }
        for(unsigned uemp=0; uemp<num_lemp; uemp++){
//printf(" uemp=%u\n",uemp);
            fNorfolkPot_pp[upot][uemp] = new DLM_Histo<float> [num_slj];
            InputFileName[upot][uemp] = new char* [num_slj];
            char semp[16];
            switch(uemp){
                case 0 : strcpy(semp,"Strong"); break;
                case 1 : strcpy(semp,"CFF"); break;
                case 2 : strcpy(semp,"CFull"); break;
                default : break;
            }
            for(unsigned uslj=0; uslj<num_slj; uslj++){
//printf("  uslj=%u\n",uslj);
                InputFileName[upot][uemp][uslj] = new char [512];
                char sslj[16];
                switch(uslj){
                    case 0 : strcpy(sslj,"1S0"); break;
                    case 1 : strcpy(sslj,"1D2"); break;
                    case 2 : strcpy(sslj,"3P0"); break;
                    case 3 : strcpy(sslj,"3P1"); break;
                    case 4 : strcpy(sslj,"3P2"); break;
                    default : break;
                }
                strcpy(InputFileName[upot][uemp][uslj],InputFolder);
                strcat(InputFileName[upot][uemp][uslj],"/norfolk");
                strcat(InputFileName[upot][uemp][uslj],spot);
                strcat(InputFileName[upot][uemp][uslj],"_");
                strcat(InputFileName[upot][uemp][uslj],semp);
                strcat(InputFileName[upot][uemp][uslj],"_");
                strcat(InputFileName[upot][uemp][uslj],sslj);
                strcat(InputFileName[upot][uemp][uslj],".txt");

                InFile = fopen(InputFileName[upot][uemp][uslj], "r");
//printf("file: %s\n",InputFileName[upot][uemp][uslj]);
                if(!InFile){
                    printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[upot][uemp][uslj]);
                    fNorfolkPot_pp[upot][uemp][uslj].SetUp(1);
                    fNorfolkPot_pp[upot][uemp][uslj].SetUp(0,10,0,100);
                    fNorfolkPot_pp[upot][uemp][uslj].Initialize();
                    continue;
                }
                unsigned NumLines=0;
                while(!feof(InFile)){
                    if(!fgets(cdummy, 255, InFile)) continue;
                    NumLines++;
                }
                fclose(InFile);

                if(NumLines<10){
                    printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName[upot][uemp][uslj]);
                    fNorfolkPot_pp[upot][uemp][uslj].SetUp(1);
                    fNorfolkPot_pp[upot][uemp][uslj].SetUp(0,10,0,100);
                    fNorfolkPot_pp[upot][uemp][uslj].Initialize();
                    continue;
                }

                unsigned NumEntries=NumLines-1;
                double* RadBinCenter = new double [NumEntries];
                double* RadBins = new double [NumEntries+1];
                float* PotVal = new float [NumEntries];

                InFile = fopen(InputFileName[upot][uemp][uslj], "r");

                if(!fgets(cdummy, 255, InFile)){
                    printf("\033[1;33mWARNING!\033[0m Possible BUG in SetUpNorfolk!\n");
                }

                unsigned WhichEntry=0;
                RadBins[0] = 0;
                while(!feof(InFile)){
                    if(!fgets(cdummy, 255, InFile)) continue;
                    sscanf(cdummy, "%f %f",&fRadius,&fPotential);

//printf(" V(%f): %f\n",fRadius,fPotential);
//usleep(125e3);
                    RadBinCenter[WhichEntry] = fRadius;
                    PotVal[WhichEntry] = fPotential;
                    if(WhichEntry>=1){
                        RadBins[WhichEntry] = (RadBinCenter[WhichEntry]+RadBinCenter[WhichEntry-1])*0.5;
                    }
                    WhichEntry++;
                }
                fclose(InFile);
                RadBins[NumEntries] = 2.*RadBinCenter[NumEntries-1]-RadBins[NumEntries-1];

                fNorfolkPot_pp[upot][uemp][uslj].SetUp(1);
                fNorfolkPot_pp[upot][uemp][uslj].SetUp(0,NumEntries,RadBins,RadBinCenter);
                fNorfolkPot_pp[upot][uemp][uslj].Initialize();
                for(unsigned uRad=0; uRad<NumEntries; uRad++){
                    fNorfolkPot_pp[upot][uemp][uslj].SetBinContent(uRad,PotVal[uRad]);
                    //printf("uRad = %u; r=%f; V = %e\n",uRad,fNorfolkPot_pp[upot][uemp][uslj].GetBinCenter(0,uRad),PotVal[uRad]);
                }

                delete [] RadBinCenter;
                delete [] RadBins;
                delete [] PotVal;
            }
        }
    }




    //Clean_SetUpNorfolk:
    for(unsigned upot=0; upot<num_lpot; upot++){
        for(unsigned uemp=0; uemp<num_lemp; uemp++){
            for(unsigned uslj=0; uslj<num_slj; uslj++){
                delete [] InputFileName[upot][uemp][uslj];
            }
            delete [] InputFileName[upot][uemp];
        }
        delete [] InputFileName[upot];
    }
    delete [] InputFileName;
    delete [] cdummy;
}

//[0] radius, [1] momentum
//[2] lpot [3] lemp [4] slj
// for the full nv2s pick lpot between the different models:  lpot=106   RL = 1.2 fm R_S=0.8 fm (up to 125 MeV)
//                             	                              lpot=105   RL = 1.0 fm R_S=0.7 fm (up to 125 MeV)
//                                                            lpot=110   RL = 1.2 fm R_S=0.8 fm (up to 200 MeV)
//                                                            lpot=109   RL = 1.0 fm R_S=0.7 fm (up to 200 MeV)
//lemp:  -1 stong only
//      100 added the effect of the Coulomb strong factor (but WITHOUT the base Coulomb)
//      101 added higher order Coulomb effects (but WITHOUT the base Coulomb)
//slj: 100 = 1S0
//     122 = 1D2
//     310 = 3P0
//     311 = 3P1
//     312 = 3P2
double pp_Norfolk(double* Pars){
    unsigned upot,uemp,uslj;
//printf("0 : %f\n",Pars[0]);
//printf("1 : %f\n",Pars[1]);
//printf("2 : %f\n",Pars[2]);
//printf("3 : %f\n",Pars[3]);
//printf("4 : %f\n",Pars[4]);
    if(Pars[2]==105) upot=0;
    else if(Pars[2]==106) upot=1;
    else if(Pars[2]==109) upot=2;
    else if(Pars[2]==110) upot=3;
    else return 0;

    if(Pars[3]==-1) uemp=0;
    else if(Pars[3]==100) uemp=1;
    else if(Pars[3]==101) uemp=2;
    else return 0;

    if(Pars[4]==100) uslj=0;
    else if(Pars[4]==122) uslj=1;
    else if(Pars[4]==310) uslj=2;
    else if(Pars[4]==311) uslj=3;
    else if(Pars[4]==312) uslj=4;
    else return 0;

    if(!fNorfolkPot_pp) return 0;
    return fNorfolkPot_pp[upot][uemp][uslj].Eval(&Pars[0]);
}

double ppDlmPot(const int& DlmPot, const int& DlmFlag, const int& Spin, const int& AngMom, const int& TotMom, double* Radius){
    return fDlmPot(DlmPot,DlmFlag,1,1,1,Spin,AngMom,TotMom,Radius,0);
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
    return fDlmPot(DlmPot,DlmFlag,0,0,0,Spin,AngMom,TotMom,Radius,0);
}
//[2] is Pot Type, 3 is [3] is PotFlag
double pLambdaDlmPot1S0(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],0,0,0,Pars);
}
double pLambdaDlmPot3S1(double* Pars){
    return ppDlmPot(Pars[2],Pars[3],1,0,1,Pars);
}

//t2p1 - 2xIsospin of particle 1, t2p2 same for particle 2
//in other pars, the 0-th par is assumed to be a cut-off parameter (i.e. V=0 below certain r)
double fDlmPot(const int& DlmPot, const int& DlmPotFlag,
               const int& IsoSpin, const int& t2p1, const int& t2p2, const int& Spin, const int& AngMom, const int& TotMom, double* Radius, const double& CutOff, double* OtherPars){
    //printf("V=%f\n",fV18potential(9,Spin,AngMom,TotMom,Radius)) ;
    if(Radius[0]<CutOff) return 0;
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
        case pXim_HALQCDPaper2020 : return LatticePots_pXi_ver3(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pXim_LatticeAvg : return LatticePots_pXi_Avg(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pXim_LatticeSqrtAvg : return LatticePots_pXi_SqrtAvg(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pKm_Tetsuo : return Tetsuo_pKm(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius);
        case pOmega_Lattice : return LatticePots_pOmega(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        case pOmega_Tetsuo : return Tetsuo_pOmega(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
        //case pXim_ESC16 : return ESC16_pXim(DlmPot,DlmPotFlag,IsoSpin,t2p1,t2p2,Spin,AngMom,TotMom,Radius,OtherPars);
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
                   round(Parameters[6]),round(Parameters[7]),round(Parameters[8]),round(Parameters[9]),Parameters,0,&Parameters[10]);
}
//[0] radius, [1] momentum
//[2] PotentialType, [3] PotentialFlag,
//[4] IsoSpin, [5] ParticleType1 (1 proton, -1 neutron), [6] ParticleType2
//[7] Spin (s), [8] AngMom (l), [9] TotMom (j)
//![10] CutOff, i.e. V(r)=0 below this point
//[11] - optional stuff
double fDlmPotVer2(double* Parameters){
    //printf(" fDlmPot called with %p\n",Parameters);
    return fDlmPot(round(Parameters[2]),round(Parameters[3]),round(Parameters[4]),round(Parameters[5]),
                   round(Parameters[6]),round(Parameters[7]),round(Parameters[8]),round(Parameters[9]),Parameters,Parameters[10],&Parameters[11]);
}

//[2,3,4,5] = B1,2,3,4
//[6] = mpi (pion mass)
//[7] = cutoff
//to use in CATS:
//Parameter 0,1,2,3 are B1,2,3,4 and Parameter 4 is the cutoff
//B1*exp(-B2*r*r)+B3*(1-exp(-b4*r*r))(exp(-2*mpi*r)/(r*r))
double LatticePots(double* Parameters){
    if(Parameters[0]<Parameters[7]) return 0;
    const double mpihalfm = Parameters[6]/197.3;
    const double& B1 = Parameters[2];
    const double& B2 = Parameters[3];
    const double& B3 = Parameters[4];
    const double& B4 = Parameters[5];

    double Result=0;
    const double& rad=Parameters[0];
    const double rad2=rad*rad;
 //f(x) = (x > 0 ? a*exp(-b*x*x)+c*(1-exp(-d*x*x))**1*(exp(-e*x)/x)**1 : a) #  for 1S0 effective
    Result = B1*exp(-B2*rad2)+B3*(1.-exp(-B4*rad2))*exp(-2.0*mpihalfm*rad)/rad2;
    return Result;
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
