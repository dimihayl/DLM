
#include "math.h"

#include "CATS.h"

//this potential is used to model proton-Lambda interaction
double UsmaniPotential(const int& Spin, double* Parameters)
{
  double r = Parameters[0];
  //Values for the potential
  const double vbar = 6.2;

  const double vsigma = 0.25;

  const double wc = Parameters[2];//default value 2137

  double x=r*0.7;
  double vc = wc/(1+exp((r-0.5)/0.2));
  double tpi = (1.0+3.0/x+3.0/(x*x)) * (exp(-x)/x) * pow(1.-exp(-2.*r*r),2.);

  double v = 0.;

  if (Spin == 0) v = vc - (vbar + 0.75*vsigma)*tpi*tpi;//Usmani singlet
  else if (Spin == 1)  v = vc - (vbar - 0.25*vsigma)*tpi*tpi;//Usmani triplet
  else printf ("wrong polarization\n");

  return v;
}

double UsmaniPotential1S0(double* Parameters){
    return UsmaniPotential(0, Parameters);
}

double UsmaniPotential3S1(double* Parameters){
    return UsmaniPotential(1, Parameters);
}

//this potential is used to model proton-proton interaction (only strong!)
double ReidPotential(const double& rad,const double& polar){
    double r = rad;//*197.3269788;//convert in fm
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

double ReidPotential1S0(double* Pars){
    return ReidPotential(Pars[0],0);
}

double ReidPotential3P(double* Pars){
    return ReidPotential(Pars[0],1);
}

double GaussSource(double* Pars){
    const double PI = 3.141592653589793;
    double& Radius = Pars[1];
    double& Size = Pars[3];
    return 4.*PI*Radius*Radius*pow(4.*PI*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}


//In this example we will define a CATS object and use it to compute:
//1) the pp correlation function using the Reid potential and a 1.3 fm Gauss source
//2) the pLambda correlation function using the Usmani potential and a 1.3 fm Gauss source
//3) again pLambda, but we change the source size by 30%
//4) again pLambda, but we change the potential (we make the repulsive core 2% stronger)
void Example1(){

///---------- SOME GENERAL VARIABLES ----------
    unsigned NumMomBins = 25;
    double kMin = 0;
    double kMax = 125;

    double MassProton = 938.272;
    double MassLambda = 1115.683;

///---------- NEXT WE CREATE THE CATS OBJECT ----------
    //ALL UNITS ARE IN fm and MeV
    CATS Kitty;

    //(#NumberOfMomBins, minMom, maxMom)
    Kitty.SetMomBins(NumMomBins,kMin,kMax);

///---------- First let us look at pp correlations --------
    //charge1 x charge2, for pp it is +1
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    //the reduced mass of the particle pair
    Kitty.SetRedMass((MassProton*MassProton)/(MassProton+MassProton));

    //true = use analytic source; false = take the source from a transport model
    Kitty.SetUseAnalyticSource(true);
    //first 3 parameters - dummies, do not set them
    //all other parametres you can pass to your source function
    //in this example we use a Gauss source, which takes only 1 parameter as input (the size)
    double SourcePars[4] = {0,0,0,1.3};
    Kitty.SetAnaSource(GaussSource, SourcePars);

    //the number of different channels to be taken into account.
    Kitty.SetNumChannels(2);
    //the number of partial waves for the 0th channel
    Kitty.SetNumPW(0,2);//(#WhichChannel, #how many partial waves)
    //the number of partial waves for the 1th channel
    Kitty.SetNumPW(1,2);
    //! N.B. the channel number should have the same oddness as the spin
    //important for the wave-function symmetrization of identical particles

    //the weight with which each channel contributes to the final wave function
    //for pp we have two channels -> 1S0 and 3PX, i.e. we treat those as a singlet and as a triplet
    //=> the ratio between them should be 1:3
    Kitty.SetChannelWeight(0, 0.25);//(#WhichChannel, #what weight)
    Kitty.SetChannelWeight(1, 0.75);

    //you can define any potential function you want and pass it any set of parameters you want.
    //the trick is to leave the first 2 parameters as placeholders (used internally by CATS)
    //in this example the potential function does not get any parameters.
    //N.B. the array you pass to CATS should always have a min. size of 2, else you might get memory issues.
    double PotPars[3];
    //the 0,0 means that we set the 0th channel (S=0), the s-wave
    Kitty.SetShortRangePotential(0,0,ReidPotential1S0,PotPars);
    //the 1,1 means that we set the 1st channel (S=1), the p-wave
    Kitty.SetShortRangePotential(1,1,ReidPotential3P,PotPars);

    //this is where the magic happens - we run CATS and all relevant parameters are computed
    Kitty.KillTheCat();

    //GetMomentum gives us the center of the bin, GetCorrFun gives us the corresponding value of C(k)
    printf("1) Result for pp:\n");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("   C(%6.2f)=%.2f  ", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        if( (uBin+1)%5==0 && uBin!=(NumMomBins-1) ) printf("\n");
    }
    printf("\n\n");


///---------- Let us try what happens for pLambda --------
//! One can reuse the same CATS object by changing all settings,
//  we will leave the source the same, i.e. it is not needed to set it up again
//  we want to use 2 channels again (1S0 and 3S1)
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    //the reduced mass of the particle pair
    Kitty.SetRedMass((MassProton*MassLambda)/(MassProton+MassLambda));

    //Note that the for the Usmani potential we use the [2]nd parameter as the strength of its repulsive core
    //the default value of this parameter is 2137 MeV
    PotPars[2] = 2137;
    //the 0,0 means that we set the 0th channel (S=0), the s-wave
    Kitty.SetShortRangePotential(0,0,UsmaniPotential1S0,PotPars);
    //the 1,0 means that we set the 1th channel (S=1), the s-wave
    Kitty.SetShortRangePotential(1,0,UsmaniPotential3S1,PotPars);

    //this is where the magic happens - we run CATS and all relevant parameters are computed
    Kitty.KillTheCat();

    printf("2) Result for pΛ:\n");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("   C(%6.2f)=%.2f  ", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        if( (uBin+1)%5==0 && uBin!=(NumMomBins-1) ) printf("\n");
    }
    printf("\n\n");

///---------- Let us increase source size by 30% --------
    SourcePars[3] *= 1.3;

    //!if you want to customize the output CATS gives you, you can choose between:
    //CATS::nSilent (no output), CATS::nErrors (only errors),
    //CATS::nWarning (errors and warning), CATS::nAll (full output -> the default option)
    Kitty.SetNotifications(CATS::nWarning);

    //!we run CATS again, however do note that if we change the parameters of our source or potential,
    //!CATS cannot know that, and to maximize performance assumes no change has been done!
    //!thus the user has to specify that by using  kSourceChanged, kPotentialChanged or kAllChanged
    Kitty.KillTheCat(CATS::kSourceChanged);

    printf("3) Result for pΛ (+30%% source size):\n");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("   C(%6.2f)=%.2f  ", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        if( (uBin+1)%5==0 && uBin!=(NumMomBins-1) ) printf("\n");
    }
    printf("\n\n");

///---------- Original source value, but we make the repulsive core 2% stronger --------
    SourcePars[3] /= 1.3;
    PotPars[2] *= 1.02;
    //!we run CATS again, however do note that if we change the parameters of our source or potential,
    //!CATS cannot know that, and to maximize performance assumes no change has been done!
    //!thus the user has to specify that by using  kSourceChanged, kPotentialChanged or kAllChanged
    Kitty.KillTheCat(CATS::kAllChanged);

    printf("4) Result for pΛ (+2%% repulsive core):\n");
    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("   C(%6.2f)=%.2f  ", Kitty.GetMomentum(uBin), Kitty.GetCorrFun(uBin));
        if( (uBin+1)%5==0 && uBin!=(NumMomBins-1) ) printf("\n");
    }
    printf("\n\n");
}


int main(int argc, char *argv[]){

    Example1();

    return 0;
}

