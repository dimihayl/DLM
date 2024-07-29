#include "DLM_Source.h"
//#include "DLM_CRAB_PM.h"
#include "DLM_Integration.h"
#include "DLM_Random.h"
#include "DLM_Bessel.h"
#include "DLM_MathFunctions.h"
#include "DLM_Histo.h"
#include "CATSconstants.h"

#include "math.h"
#include <sstream>
#include <iostream>
#include <vector>

//!TEST
//#include <fstream>
//#include <stdio.h>
#include <unistd.h>

double GaussSourceTF1(double* x, double* Pars){
    double& Radius = *x;
    double& Size = Pars[0];
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

double GaussSourceScaledTF1(double* x, double* Pars){
    double& Radius = *x;
    double& Size = Pars[0];
    double& Alpha = Pars[1];
    return Alpha*4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

//2-particle
double GaussSource(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    // printf("Radius = %.2f\n",Radius);
    // printf("Size = %.2f\n",Size);

//printf(" G-> r=%.2f, s=%.2f => %.2f\n",Radius,Size,4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size)));
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

double GaussSourceKstarPol2(double* Pars){
    double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double Size = Pars[3];

    Size += Pars[4]*Momentum*1e-3;
    Size += Pars[5]*Momentum*Momentum*1e-6;
    return 4.*Pi*Radius*Radius*pow(4.*Pi*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}


//this source will not be normalized! However, CATS does a renormalization itself, so it should be fine.
//in case of problems, contact the developers!
double GaussSourceCutOff(double* Pars){
    if(Pars[1]<Pars[4]) return 0;
    return GaussSource(Pars);
}

//same as GaussSource, but we assume that the user calls the function wanting to sample from theta as well.
//since Integral(dTheta) = 2 for a flat theta distribution and the whole Source function needs to be normalized to 1,
//in order to preserve this we should divide the whole function by 2
double GaussSourceTheta(double* Pars){
//printf("HI\n");
//printf("CosTheta=%f\n",Pars[2]);
//if(Pars[2]>0.95) return GaussSource(Pars);
//return 0;
return GaussSource(Pars)*exp(-pow((Pars[2]-1)/(2.*0.4),2.));
    return 0.5*GaussSource(Pars);
}

//double CauchySource(const double& Radius, const double& Size){
//based on single particle emission, which leads to a rather complicated expression, that is approximated with the formula below
double CauchySource(double* Pars){
    double& Radius = Pars[1];
    double& Size = Pars[3];
    return 2.97*Size*sqrt(2)*Radius*Radius/Pi*pow(Radius*Radius+0.125*2.97*2.97*Size*Size,-2.);

}
double CauchySourceTheta(double* Pars){
    return 0.5*CauchySource(Pars);
}

double CauchySource_v2(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 12.*Radius*Size*sqrt(2)/(Pi*Pi*sqrt(Radius*Radius+4.*Size*Size)*(Radius*Radius+6.*Size*Size))*
            atan(sqrt(2.+0.5*pow(Radius/Size,2.))*Radius/Size/sqrt(2));
}

//follows the definition here (Eur. Phys. J. C 36, 67–78 (2004))
//it has the same functional shape as CauchySource, but the parameterization is slightly different,
//i.e. there is going to be a multiplicative offset between the obtained radii
double ExponentialSource(double* Pars){
    double& Radius = Pars[1];
    double& Size = Pars[3];
    return 4.*Radius*Radius*Size*pow(Radius*Radius+Size*Size,-2.);
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

double NormDoubleGaussSource(double* Pars){
  return Pars[6]*DoubleGaussSource(Pars);
}
//[3]*DGS(size1 [0],size2 [1],weight1 [2])
double NormDoubleGaussSourceTF1(double* x, double* Pars){
  double PARS[7];
  PARS[0] = 0;
  PARS[1] = *x;
  PARS[2] = 0;
  for(int i=0; i<4; i++)
    PARS[3+i] = Pars[i];
  return NormDoubleGaussSource(PARS);
}

//Shifted Gaussian will be G(x-S,...), where S is the shift on x axis
//if x<S, the function is zero !!!
//here we have a scaled sum of 3 such Gaussians
//[0]_3: overall normalization
//[1]_4: sigma1
//[2]_5: shift1
//[3]_6: weight1
//[4]_7: sigma2
//[5]_8: shift2
//[6]_9: weight2
//[7]_10: sigma3
//[8]_11: shift3
//weight3 = 1.-weight1-(1-weight1) * weight2 (i.e. limits between 0-1 and weighted such that it cannot go beyond allowed vals)
double NormTripleShiftedGauss(double* Pars){
  static double Pars1[4];
  static double Pars2[4];
  static double Pars3[4];

  Pars1[1] = Pars[1]-Pars[5];;//rad
  Pars1[3] = Pars[4];//size

  Pars2[1] = Pars[1]-Pars[8];//rad
  Pars2[3] = Pars[7];//size

  Pars3[1] = Pars[1]-Pars[11];//rad
  Pars3[3] = Pars[10];//size

  double G1 = Pars1[1]>0?GaussSource(Pars1):0;
  double G2 = Pars2[1]>0?GaussSource(Pars2):0;
  double G3 = Pars3[1]>0?GaussSource(Pars3):0;
  double& W1 = Pars[6];
  double W2 = (1.-W1)*Pars[9];
  double W3 = 1.-W1-W2;
  //double SW = 1;
  if(W1<0 || W1>1 || W2<0 || W2>1 || W3<0 || W3>1){
    static bool ErrorShown = false;
    if(!ErrorShown){
      printf("\033[1;31mERROR:\033[0m (NormTripleShiftedGauss) The weights have unphysical values, possible ERROR in the fit!!!\n");
      printf(" W1, W2, W3 = %.3e %.3e %.3e\n",W1,W2,W3);
      ErrorShown = true;
    }
  }

  return Pars[3]*(W1*G1+W2*G2+W3*G3);
}
double NormTripleShiftedGaussTF1(double* x, double* Pars){
  double PARS[12];
  PARS[0] = 0;
  PARS[1] = *x;
  PARS[2] = 0;
  for(int i=0; i<9; i++)
    PARS[3+i] = Pars[i];
  return NormTripleShiftedGauss(PARS);
}

//[3]_{0} : Number of Gaussians (HAS TO BE FIXED!!!)
//[4]_{1+i*3} : Weight of i-th Gaussian with respect the sum of the weights of all previous ones
//[5]_{2+i*3} : Mean of the i-th Gaussian
//[6]_{3+i*3} : Stdv of the i-th Gaussian
//Table for #Gaussians vs Number of pars:
//            1             4
//            2             7
//            3             10
//etc.
double StupidGaussSum(double* Pars){
  //double& Mom = Pars[0];
  double& Rad = Pars[1];
  if(Rad<0) return 0;
  //double& CosTh = Pars[2];
  unsigned NG = round(Pars[3]);
  double Result=0;
  double WeightSum=0;
  double Weight;
  for(unsigned uG=0; uG<NG; uG++){
    Weight = (1.-WeightSum)*Pars[4+uG*3];
    static bool ErrorShown = false;
    if(!ErrorShown && (Weight<0 || Weight>1)){
      printf("\033[1;31mERROR:\033[0m (StupidGaussSum) The weights have unphysical values, possible ERROR in the fit!!!\n");
      ErrorShown = true;
    }
    WeightSum += Weight;
    double NormDist = pow(2.*Pi,-0.5)/Pars[6+uG*3]*exp(-0.5*pow((Rad-Pars[5+uG*3])/Pars[6+uG*3],2.));
    Result += Weight*NormDist;
    //static int counter=0;
    //if(counter<4 && Rad>1.){
    //  printf("N%u NormDist = %f\n",NG,NormDist);
    //  printf(" W=%f (%f)\n",Weight,Pars[4+uG*3]);
    //  printf(" M=%f\n",Pars[5+uG*3]);
    //  printf(" S=%f\n",Pars[6+uG*3]);
    //  counter++;
    //}
  }
  return Result;
}
double StupidGaussSumTF1(double* x, double* Pars){
  unsigned NG = round(Pars[0]);
  double PARS[NG*3+3+1];
  PARS[0] = 0;
  PARS[1] = *x;
  PARS[2] = 0;
  //printf("N%u StupidGaussSumTF1\n",NG);
  for(unsigned u=0; u<NG*3+1; u++){
    PARS[3+u] = Pars[u];
    //printf("%u %f\n",u,Pars[u]);
    //usleep(1e6);
  }

  return StupidGaussSum(PARS);
}

//sum of many possions, all weighted such that the total weight is still 1
//this is achieved by using the weight parameters as the reletive weight with respect the
//"remaining" weight (i.e. 1 - weight of all previous poissons). As long as all weight pars are within 0-1 this will work
double PoissonSum(double* xVal, double* Pars){
  static KdpPars kdppars;
  for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
    kdppars.mean[uP] = Pars[0+uP*3];
    kdppars.stdv[uP] = Pars[1+uP*3];
    if(uP!=KdpPars::NumDistos-1)
      kdppars.wght[uP] = Pars[2+uP*3];
  }
  return PoissonSum(*xVal,kdppars);
}

double PoissonSum(const double& xVal, const KdpPars& kdppars){
  double Rslt=0;
  double Pssn;
  double Norm;
  double Std2;
  double Mean;
  double RemainingNorm = 1;
  for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
    Mean = kdppars.mean[uP];
    Std2 = kdppars.stdv[uP]*kdppars.stdv[uP];
    if(uP!=KdpPars::NumDistos-1){
      Norm = RemainingNorm*kdppars.wght[uP];
      if(Norm>=RemainingNorm && kdppars.wght[uP]<=1){
        Norm = RemainingNorm;
      }
    }
    else Norm = RemainingNorm;
    if(Norm<0){
      printf("\033[1;33mWARNING!\033[0m PoissonSum a negative (%.3e) norm! kdppars.wght[uP]=%.3e\n",Norm,kdppars.wght[uP]);
    }
    if(Std2) Pssn = DLM_Poisson(xVal*Mean/Std2,Mean*Mean/Std2)*Mean/Std2;
    else Pssn = 0;
    Rslt += Norm*Pssn;
    if(RemainingNorm==Norm) RemainingNorm=0;
    else RemainingNorm -= Norm;
    if(RemainingNorm<0 && RemainingNorm>-1e-6){
      RemainingNorm = 0;
    }
    else if(RemainingNorm<0){
      if(kdppars.wght[uP]>=0 && kdppars.wght[uP]<=1){
        RemainingNorm = 0;
      }
      else{
        printf("\033[1;33mWARNING!\033[0m Could not properly correct the Norm! Check the weights!!! RemainingNorm=%.3e\n",RemainingNorm);
      }
    }
  }

  return Rslt;
}

/*
//These are femto Gaussinas (3D)
//[3]_{0} : Number of Gaussians (HAS TO BE FIXED!!!)
//[4]_{1+i*3} : Weight of i-th Gaussian with respect the sum of the weights of all previous ones
//[5]_{2+i*3} : Mean of the i-th Gaussian
//[6]_{3+i*3} : Stdv of the i-th Gaussian
//Table for #Gaussians vs Number of pars:
//            1             4
//            2             7
//            3             10
//etc.
double StupidShiftedGaussSum(double* Pars){
  //double& Mom = Pars[0];
  double& Rad = Pars[1];
  if(Rad<0) return 0;
  //double& CosTh = Pars[2];
  unsigned NG = round(Pars[3]);
  double Result=0;
  double WeightSum=0;
  double Weight;
  double PAR[4];
  for(unsigned uG=0; uG<NG; uG++){
    Weight = (1.-WeightSum)*Pars[4+uG*3];
    static bool ErrorShown = false;
    if(!ErrorShown && (Weight<0 || Weight>1)){
      printf("\033[1;31mERROR:\033[0m (StupidGaussSum) The weights have unphysical values, possible ERROR in the fit!!!\n");
      ErrorShown = true;
    }
    WeightSum += Weight;
    PAR[0]=Pars[0];
    PAR[1]=Pars[1]-Pars[5];
    PAR[2]=Pars[2];
    PAR[3]=Pars[6];
    Result += Weight*GaussSource(Pars);
  }
  return Result;
}
double StupidShiftedGaussSumTF1(double* x, double* Pars){
  unsigned NG = round(Pars[0]);
  double PARS[NG*2+3+1];
  PARS[0] = 0;
  PARS[1] = *x;
  PARS[2] = 0;
  //printf("N%u StupidGaussSumTF1\n",NG);
  for(unsigned u=0; u<NG*2+1; u++){
    PARS[3+u] = Pars[u];
  }

  return StupidShiftedGaussSum(PARS);
}

*/


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

/*
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
*/
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
    ResoDecayTopology = new RESO_DEC_TOP* [2];
    ResoEmissionAngle = new const DLM_Histo<double>** [2];
    for(unsigned uParticle=0; uParticle<2; uParticle++){
        NumResonances[uParticle]=0;
        ResoWeight[uParticle]=NULL;
        ResoMass[uParticle]=NULL;
        ResoTau[uParticle]=NULL;
        ChildMass0[uParticle]=NULL;
        ChildMass1[uParticle]=NULL;
        SmearResoMomentum[uParticle]=NULL;
        SmearResoMass[uParticle]=NULL;
        ResoDecayTopology[uParticle]=NULL;
        ResoEmissionAngle[uParticle]=NULL;
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
    if(ResoDecayTopology[whichparticle]){delete[]ResoDecayTopology[whichparticle];ResoDecayTopology[whichparticle]=new RESO_DEC_TOP[NumResonances[whichparticle]];}
    if(ResoEmissionAngle[whichparticle]){delete[]ResoEmissionAngle[whichparticle];ResoEmissionAngle[whichparticle]=new const DLM_Histo<double>*[NumResonances[whichparticle]];}
    if(Histo) Histo->SetBinContentAll(1e6);
}
void DLM_CleverMcLevyReso::SetUpReso(const unsigned& whichparticle, const unsigned& whichreso, const double& weight, const double& mass, const double& tau, const double& mass0, const double& mass1, const double& momSmear, const bool& massSmear, const RESO_DEC_TOP& rdt){
    if(whichparticle>=2) {printf("\033[1;33mWARNING:\033[0m You can call SetUpReso only for particle 0 or 1\n"); return;}
    if(whichreso>=NumResonances[whichparticle]) {printf("\033[1;33mWARNING:\033[0m Only %u number of resonances are currently allowed. Change using InitReso\n",NumResonances[whichparticle]); return;}
    if(!ResoWeight[whichparticle]){ResoWeight[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ResoMass[whichparticle]){ResoMass[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ResoTau[whichparticle]){ResoTau[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ChildMass0[whichparticle]){ChildMass0[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!ChildMass1[whichparticle]){ChildMass1[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!SmearResoMomentum[whichparticle]){SmearResoMomentum[whichparticle]=new double[NumResonances[whichparticle]];}
    if(!SmearResoMass[whichparticle]){SmearResoMass[whichparticle]=new bool[NumResonances[whichparticle]];}
    if(!ResoDecayTopology[whichparticle]){ResoDecayTopology[whichparticle]=new RESO_DEC_TOP[NumResonances[whichparticle]];}
    if(ResoWeight[whichparticle][whichreso]==weight&&ResoMass[whichparticle][whichreso]==mass&&ResoTau[whichparticle][whichreso]==tau&&
       ChildMass0[whichparticle][whichreso]==mass0&&ChildMass1[whichparticle][whichreso]==mass1&&
       SmearResoMomentum[whichparticle][whichreso]==momSmear&&SmearResoMass[whichparticle][whichreso]==massSmear&&
       ResoDecayTopology[whichparticle][whichreso]==rdt){return;}
    if(Histo) Histo->SetBinContentAll(1e6);
    ResoWeight[whichparticle][whichreso] = weight;
    ResoMass[whichparticle][whichreso] = mass;
    ResoTau[whichparticle][whichreso] = tau*FmToNu;
    ChildMass0[whichparticle][whichreso] = mass0;
    ChildMass1[whichparticle][whichreso] = mass1;
    SmearResoMomentum[whichparticle][whichreso] = momSmear;
    SmearResoMass[whichparticle][whichreso] = massSmear;
    ResoDecayTopology[whichparticle][whichreso] = rdt;
//printf("ResoWeight=%f\n",ResoWeight[whichparticle][whichreso]);
}
//in which direction is the resonance emitted
//only used if ResoDecayTopology==rdtRandom
void DLM_CleverMcLevyReso::SetUpResoEmission(const unsigned& whichparticle, const unsigned& whichreso, const DLM_Histo<double>* Distr){
    if(whichparticle>=2) {printf("\033[1;33mWARNING:\033[0m You can call SetUpResoEmission only for particle 0 or 1\n"); return;}
    if(whichreso>=NumResonances[whichparticle]) {printf("\033[1;33mWARNING:\033[0m Only %u number of resonances are currently allowed. Change using InitReso\n",NumResonances[whichparticle]); return;}
    if(!ResoEmissionAngle[whichparticle]){ResoEmissionAngle[whichparticle]=new const DLM_Histo<double>*[NumResonances[whichparticle]];}
    ResoEmissionAngle[whichparticle][whichreso] = Distr;
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
//double OLDRAD;
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
//int SmallerHalf=0;
//int BiggerHalf=0;
//double Mean=0;
//int MeanNorm=0;
                //#pragma omp for
                for(unsigned uIter=0; uIter<NumMcIter; uIter++){
                    //vector coordinates for the flight path of the two resonances
                    //double ResoPath0_X,ResoPath0_Y,ResoPath0_Z,ResoPath0_LEN,ResoPath0_THETA,ResoPath0_PHI;
                    //double ResoPath1_X,ResoPath1_Y,ResoPath1_Z,ResoPath1_LEN,ResoPath1_THETA,ResoPath1_PHI;
                    double ResoPathCartesian[2][3];
                    double ResoPathSpherical[2][3];
//if(uIter%10000==0)
//printf("   %u/%u\n",uIter,NumMcIter);
                    //we simulate random 'core' distance RAD
                    if(Type==0) {RAD = RanGen.StableR(3,par_stability,0,par_scale,0);}
                    else if(Type==1) {RAD = RanGen.StableDiffR(3,par_stability,0,par_scale,0);}
                    else {RAD = RanGen.StableNolan(3,par_stability,0,par_scale,0);}
//OLDRAD = RAD;
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
                        ResoPathCartesian[uParticle][0] = 0;
                        ResoPathCartesian[uParticle][1] = 0;
                        ResoPathCartesian[uParticle][2] = 0;
//printf("WhichReso=%i\n",WhichReso);
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
                        //compute the shift of the radius (in terms of length)
                        ResoPathSpherical[uParticle][0] = RanGen.Exponential(TKM)*NuToFm;
//printf(" ResoTau = %f\n",ResoTau[uParticle][WhichReso]*NuToFm);
//if(ResoPathSpherical[uParticle][0]>ResoTau[uParticle][WhichReso]*NuToFm*log(2)) BiggerHalf++;
//else SmallerHalf++;
//Mean+=ResoPathSpherical[uParticle][0];
//MeanNorm++;
//printf(" FlyPath = %f\n",ResoPathSpherical[uParticle][0]);
//ResoPathSpherical[uParticle][0] = 0;
                        switch(ResoDecayTopology[uParticle][WhichReso]){
                        case rdtBackwards :
                            //make sure that the direction is inverted for the second particle
                            ResoPathSpherical[uParticle][1] = uParticle?0:Pi;
                            ResoPathSpherical[uParticle][2] = 0;
                            break;
                        //random
                        case rdtRandom :
                            if(ResoEmissionAngle&&ResoEmissionAngle[uParticle]&&ResoEmissionAngle[uParticle][WhichReso]){
                                //the theta angle we sample from the distribution we gave
                                RanVal = RanGen.Uniform(1e-6,1-1e-6);//0-->1
                                ResoPathSpherical[uParticle][1] = ResoEmissionAngle[uParticle][WhichReso]->Eval(&RanVal);
                                //make sure that the direction is inverted for the second particle
                                if(uParticle) ResoPathSpherical[uParticle][1] = fabs(ResoPathSpherical[uParticle][1]-Pi);

                                //the phi angle is still
                                //this is probably not the most realistic thing, but should be rather conservative.
                                ResoPathSpherical[uParticle][2] = RanGen.Uniform(0,2.*Pi);

//ResoPathSpherical[uParticle][2] = 0;
                            }
                            else
                            {
                                //ResoPathSpherical[uParticle][1] = RanGen.Uniform(0,Pi);//theta
                                //the upper thing was wrong, as actually for theta it is the cosine that is random, now corrected
                                RanVal = acos(RanGen.Uniform(-1.,1.));
                                ResoPathSpherical[uParticle][1] = RanVal;
                                ResoPathSpherical[uParticle][2] = RanGen.Uniform(0,2.*Pi);//phi
                            }
                            break;
                        //rdtRandomBackwards
                        default :
                            RanVal = acos(RanGen.Uniform(-1.,0));
                            ResoPathSpherical[uParticle][1] = RanVal;
                            //make sure that the direction is inverted for the second particle
                            if(uParticle) ResoPathSpherical[uParticle][1] = fabs(ResoPathSpherical[uParticle][1]-Pi);
                            ResoPathSpherical[uParticle][2] = RanGen.Uniform(0,2.*Pi);//phi
                            break;
                        }
                        //RAD += RanGen.Exponential(TKM)*NuToFm;
//OLDRAD += ResoPathSpherical[uParticle][0];
                        ResoPathCartesian[uParticle][0] = ResoPathSpherical[uParticle][0]*sin(ResoPathSpherical[uParticle][1])*cos(ResoPathSpherical[uParticle][2]);
                        ResoPathCartesian[uParticle][1] = ResoPathSpherical[uParticle][0]*sin(ResoPathSpherical[uParticle][1])*sin(ResoPathSpherical[uParticle][2]);
                        ResoPathCartesian[uParticle][2] = ResoPathSpherical[uParticle][0]*cos(ResoPathSpherical[uParticle][1]);
//printf(" up%u: x=%f, y=%f, z=%f\n",uParticle,ResoPathCartesian[uParticle][0],ResoPathCartesian[uParticle][1],ResoPathCartesian[uParticle][2]);
//printf(" eRAD = %f\n",RAD);
                    }

                    //in vectors, |R| = |R0-R1|, with R0 = (coordinates of path0), R1 = (coordinates of path1)
                    //the initial position of R0 is 0,0,0, of R1 is 0,0,r_core
//printf(" RAD = sqrt[(%.3f-%.3f)^2+(%.3f-%.3f)^2+(%.3f-%.3f-%.3f)^2]\n",
//       ResoPathCartesian[0][0],ResoPathCartesian[1][0],
//       ResoPathCartesian[0][1],ResoPathCartesian[1][1],
//       ResoPathCartesian[0][2],ResoPathCartesian[1][2],RAD);
//printf(" RAD=%f --> ",RAD);
                    RAD = sqrt(
                                pow(ResoPathCartesian[0][0]-ResoPathCartesian[1][0],2.)+
                                pow(ResoPathCartesian[0][1]-ResoPathCartesian[1][1],2.)+
                                pow(ResoPathCartesian[0][2]-ResoPathCartesian[1][2]-RAD,2.)
                                );
//printf(" %f (%f)\n",RAD,OLDRAD);


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
//printf("Small/Big = %f\n",double(SmallerHalf)/double(BiggerHalf));
//printf(" Mean = %f\n",double(Mean)/double(MeanNorm));
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







DLM_CleverMcLevyResoTM::DLM_CleverMcLevyResoTM(){
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
    ResoWeight = new double [2];
    BGT_PR = NULL;
    BGT_RP = NULL;
    BGT_RR = NULL;
    MaxBGT_PR = 0;
    MaxBGT_RP = 0;
    MaxBGT_RR = 0;
    NumBGT_PR = 0;
    NumBGT_RP = 0;
    NumBGT_RR = 0;
    NumMcIter = 1000000;
    Normalization = 1;
}
DLM_CleverMcLevyResoTM::~DLM_CleverMcLevyResoTM(){
    if(Histo) {delete Histo;Histo=NULL;}
    delete [] ResoWeight; ResoWeight=NULL;
}
void DLM_CleverMcLevyResoTM::InitStability(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsStability==numPts&&MinStability==minVal&&MaxStability==maxVal) return;
    Reset();
    NumPtsStability=numPts;
    MinStability=minVal;
    MaxStability=maxVal;
}
void DLM_CleverMcLevyResoTM::InitScale(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsScale==numPts&&MinScale==minVal&&MaxScale==maxVal) return;
    Reset();
    NumPtsScale=numPts;
    MinScale=minVal;
    MaxScale=maxVal;
}
void DLM_CleverMcLevyResoTM::InitRad(const unsigned& numPts, const double& minVal, const double& maxVal){
    if(NumPtsRad==numPts&&MinRad==minVal&&MaxRad==maxVal) return;
    Reset();
    NumPtsRad=numPts;
    MinRad=minVal;
    MaxRad=maxVal;
}
void DLM_CleverMcLevyResoTM::InitType(const int& type){
    if(type<0 || type>1) Type = 1;
    Type = type;
}
void DLM_CleverMcLevyResoTM::SetUpReso(const unsigned& whichparticle, const double& weight){
    if(whichparticle>=2) {printf("\033[1;33mWARNING:\033[0m You can call SetUpReso only for particle 0 or 1\n"); return;}
    if(Histo) Histo->SetBinContentAll(1e6);
    ResoWeight[whichparticle] = weight;
}
void DLM_CleverMcLevyResoTM::AddBGT_PR(const float& bgt,const float& a_cp){
    if(NumBGT_PR>=MaxBGT_PR){
        float** tempfloat = NULL;
        //if(BGT_PR){
        tempfloat = new float* [NumBGT_PR];
        for(unsigned uEntry=0; uEntry<NumBGT_PR; uEntry++){
            tempfloat[uEntry] = new float [2];
            for(unsigned uBGT=0; uBGT<2; uBGT++){
                tempfloat[uEntry][uBGT] = BGT_PR[uEntry][uBGT];
            }
            delete [] BGT_PR[uEntry];
        }
        delete [] BGT_PR;
        MaxBGT_PR += 100000;
        BGT_PR = new float* [MaxBGT_PR];
        for(unsigned uEntry=0; uEntry<MaxBGT_PR; uEntry++){
            BGT_PR[uEntry] = new float [2];
            if(uEntry>=NumBGT_PR) continue;
            for(unsigned uBGT=0; uBGT<2; uBGT++){
                BGT_PR[uEntry][uBGT] = tempfloat[uEntry][uBGT];
            }
        }
        //}
        //else{BGT_PR = new float* [MaxBGT_PR];for(unsigned uEntry=0; uEntry<MaxBGT_PR; uEntry++)BGT_PR[uEntry] = new float [2];}
        for(unsigned uEntry=0; uEntry<NumBGT_PR; uEntry++){
          delete [] tempfloat[uEntry];
        }
        delete [] tempfloat;
    }
    BGT_PR[NumBGT_PR][0] = bgt;
    BGT_PR[NumBGT_PR][1] = a_cp;
    NumBGT_PR++;
}
void DLM_CleverMcLevyResoTM::AddBGT_RP(const float& bgt,const float& a_cp){
    if(NumBGT_RP>=MaxBGT_RP){
        float** tempfloat = NULL;
        //if(BGT_RP){
        tempfloat = new float* [NumBGT_RP];
        for(unsigned uEntry=0; uEntry<NumBGT_RP; uEntry++){
            tempfloat[uEntry] = new float [2];
            for(unsigned uBGT=0; uBGT<2; uBGT++){
                tempfloat[uEntry][uBGT] = BGT_RP[uEntry][uBGT];
            }
            delete [] BGT_RP[uEntry];
        }
        delete [] BGT_RP;
        MaxBGT_RP += 100000;
        BGT_RP = new float* [MaxBGT_RP];
        for(unsigned uEntry=0; uEntry<MaxBGT_RP; uEntry++){
            BGT_RP[uEntry] = new float [2];
            if(uEntry>=NumBGT_RP) continue;
            for(unsigned uBGT=0; uBGT<2; uBGT++){
                BGT_RP[uEntry][uBGT] = tempfloat[uEntry][uBGT];
            }
        }
        //}
        //else{BGT_RP = new float* [MaxBGT_RP];for(unsigned uEntry=0; uEntry<MaxBGT_RP; uEntry++)BGT_RP[uEntry] = new float [2];}
        for(unsigned uEntry=0; uEntry<NumBGT_RP; uEntry++){
          delete [] tempfloat[uEntry];
        }
        delete [] tempfloat;
    }
    BGT_RP[NumBGT_RP][0] = bgt;
    BGT_RP[NumBGT_RP][1] = a_cp;
    NumBGT_RP++;
}
void DLM_CleverMcLevyResoTM::AddBGT_RR(const float& bgt0,const float& a_cp0,const float& bgt1,const float& a_cp1,const float& a_p0p1){
    if(NumBGT_RR>=MaxBGT_RR){
        float** tempfloat = NULL;
        //if(BGT_RR){
        tempfloat = new float* [NumBGT_RR];
        for(unsigned uEntry=0; uEntry<NumBGT_RR; uEntry++){
            tempfloat[uEntry] = new float [5];
            for(unsigned uBGT=0; uBGT<5; uBGT++){
                tempfloat[uEntry][uBGT] = BGT_RR[uEntry][uBGT];
            }
            delete [] BGT_RR[uEntry];
        }
        delete [] BGT_RR;
        MaxBGT_RR += 100000;
        BGT_RR = new float* [MaxBGT_RR];
        for(unsigned uEntry=0; uEntry<MaxBGT_RR; uEntry++){
            BGT_RR[uEntry] = new float [5];
            if(uEntry>=NumBGT_RR) continue;
            for(unsigned uBGT=0; uBGT<5; uBGT++){
                BGT_RR[uEntry][uBGT] = tempfloat[uEntry][uBGT];
            }
        }
        //}
        //else{BGT_RR = new float* [MaxBGT_RR];for(unsigned uEntry=0; uEntry<MaxBGT_RR; uEntry++)BGT_RR[uEntry] = new float [5];}
        for(unsigned uEntry=0; uEntry<NumBGT_RR; uEntry++){
          delete [] tempfloat[uEntry];
        }
        delete [] tempfloat;
    }
    BGT_RR[NumBGT_RR][0] = bgt0;
    BGT_RR[NumBGT_RR][1] = a_cp0;
    BGT_RR[NumBGT_RR][2] = bgt1;
    BGT_RR[NumBGT_RR][3] = a_cp1;
    BGT_RR[NumBGT_RR][4] = a_p0p1;
    NumBGT_RR++;
}
void DLM_CleverMcLevyResoTM::InitNumMcIter(const unsigned& numiter){
    if(NumMcIter==numiter) return;
    if(Histo) Histo->SetBinContentAll(1e6);
    NumMcIter=numiter;
}
void DLM_CleverMcLevyResoTM::SetNormalization(const double& norm){
  Normalization = norm;
}
double DLM_CleverMcLevyResoTM::RootEval(double* x, double* Pars){
//printf(" B RootEval\n");
    double PARS[5];
    PARS[1] = *x;
    PARS[3] = Pars[0];
    PARS[4] = Pars[1];
//printf(" E RootEval\n");
    return Normalization*Eval(PARS);
}
double DLM_CleverMcLevyResoTM::RootEvalNorm(double* x, double* Pars){
    double PARS[5];
    PARS[1] = *x;
    PARS[3] = Pars[0];
    PARS[4] = Pars[1];
    return Pars[2]*Eval(PARS);
}
double DLM_CleverMcLevyResoTM::Eval(double* Pars){
//printf("Hello\n");
    if(!Histo) {Init();}
    if(!Histo) {return -1; printf("Possible problem in DLM_CleverMcLevyResoTM::Eval\n");}
    double& Radius = Pars[1];
    double& Scale = Pars[3];
    double& Stability = Pars[4];
    const double RSS[3] = {Radius,Scale,Stability};
    //int RadBin = Histo->GetBin(0,Radius);
    int ScaleBin = Histo->GetBin(1,Scale);
    int StabilityBin = Histo->GetBin(2,Stability);
    if(ScaleBin<0||ScaleBin>=int(Histo->GetNbins())) {printf("\033[1;33mWARNING!\033[0m A bad scaling parameter passed into DLM_CleverMcLevyResoTM::Eval\n"); return 0;}
    if(StabilityBin<0||StabilityBin>=int(Histo->GetNbins())) {printf("\033[1;33mWARNING!\033[0m A bad stability parameter passed into DLM_CleverMcLevyResoTM::Eval\n"); return 0;}
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
                //the beta*gamma*tau correction for each particle
                double BGT[2];BGT[0]=0;BGT[1]=0;
                double CosRcP0=0;
                double CosRcP1=0;
                double CosP0P1=0;

                //#pragma omp for
                for(unsigned uIter=0; uIter<NumMcIter; uIter++){
                    //we simulate random 'core' distance RAD
                    if(Type==0) {RAD = RanGen.StableR(3,par_stability,0,par_scale,0);}
                    else if(Type==1) {RAD = RanGen.StableDiffR(3,par_stability,0,par_scale,0);}
                    else {RAD = RanGen.StableNolan(3,par_stability,0,par_scale,0);}
                    //after that we throw a random dice to decide IF each particle is primary or not

                    bool IsReso[2];
                    IsReso[0] = false;
                    IsReso[1] = false;

                    for(unsigned uParticle=0; uParticle<2; uParticle++){
                        if(!ResoWeight||!ResoWeight[uParticle]) continue;
                        RanVal = RanGen.Uniform(0,1);
                        IsReso[uParticle] = RanVal<ResoWeight[uParticle];
                    }

                    int RanInt;

                    if(IsReso[0]&&IsReso[1]){
                        if(NumBGT_RR){
                            RanInt = RanGen.Integer(0,NumBGT_RR);
                            BGT[0] = BGT_RR[RanInt][0];
                            CosRcP0 = BGT_RR[RanInt][1];
                            BGT[1] = BGT_RR[RanInt][2];
                            CosRcP1 = BGT_RR[RanInt][3];
                            CosP0P1 = BGT_RR[RanInt][4];
                        }
                    }
                    else if(IsReso[0]){
                        if(NumBGT_RP){
                            RanInt = RanGen.Integer(0,NumBGT_RP);
                            BGT[0] = BGT_RP[RanInt][0];
                            CosRcP0 = BGT_RP[RanInt][1];
                            BGT[1] = 0;
                            CosRcP1 = 0;
                            CosP0P1 = 0;
//if(ResoWeight[0]==1) printf("0\n");
                        }
                    }
                    else if(IsReso[1]){
                        if(NumBGT_PR){
                            RanInt = RanGen.Integer(0,NumBGT_PR);
                            BGT[0] = 0;
                            CosRcP0 = 0;
                            BGT[1] = BGT_PR[RanInt][0];
                            CosRcP1 = BGT_PR[RanInt][1];
                            CosP0P1 = 0;
//if(ResoWeight[0]==1) printf("1\n");
                        }
                    }
                    else{
                        BGT[0] = 0;
                        BGT[1] = 0;
                        CosRcP0 = 0;
                        CosRcP1 = 0;
                        CosP0P1 = 0;
//if(ResoWeight[0]==1) printf("C\n");
                    }

//if(ResoWeight[0]==1){
//if(RAD*RAD+BGT[0]*BGT[0]+BGT[1]*BGT[1]
//           -2.*RAD*BGT[0]*CosRcP0+2.*RAD*BGT[1]*CosRcP1-2.*BGT[0]*BGT[1]*CosP0P1<0){
//printf("RAD = %f\n",RAD);
//printf("BGT[0] = %f\n",BGT[0]);
//printf("BGT[1] = %f\n",BGT[1]);
//printf("CosRcP0 = %f\n",CosRcP0);
//printf("CosRcP1 = %f\n",CosRcP1);
//printf("CosP0P1 = %f\n",CosP0P1);
//usleep(100e3);
//printf("\n");
//}

//static unsigned counter=0;
//counter++;
//double rad=RAD;
//static double rad_avg=0;
//rad_avg+=rad;
//printf("r_core = %.1f (%.1f) --> ",RAD,rad_avg/double(counter));
                    //the sign convention is such, that r_core = primary_1 - primary_0

                    RAD = sqrt(RAD*RAD+BGT[0]*BGT[0]+BGT[1]*BGT[1]
                               -2.*RAD*BGT[0]*CosRcP0+2.*RAD*BGT[1]*CosRcP1-2.*BGT[0]*BGT[1]*CosP0P1);
//static double RAD_avg=0;
//RAD_avg+=RAD;

//printf("r_star = %.1f (%.1f)\n",RAD,RAD_avg/double(counter));
//printf(" sqrt(%.1f^2 + %.3f^2 + %.3f^2 - 2*%.1f*%.3f*%.3f + 2*%.1f*%.3f*%.3f - 2*%.3f*%.3f*%.3f)\n",
//rad,BGT[0],BGT[1],rad,BGT[0],CosRcP0,rad,BGT[1],CosRcP1,BGT[0],BGT[1],CosP0P1);
//usleep(25e3);
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
                }
            }
        }
    }
    double RETVAL = Histo->Eval(RSS);
    //if(Normalization*RETVAL<0||Normalization*RETVAL>1){
    //  printf("RETVAL = %f*%f\n",Normalization,RETVAL);
    //  usleep(1000e3);
    //}
    return Normalization*RETVAL;

}
unsigned DLM_CleverMcLevyResoTM::GetNumPars(){
    return 2;
}
void DLM_CleverMcLevyResoTM::Reset(){
    if(Histo) {delete Histo;Histo=NULL;}
}
void DLM_CleverMcLevyResoTM::Init(){
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


DLM_HistoSource::DLM_HistoSource(DLM_Histo<float>* histo):MyOwnHisto(true){
  if(histo->GetDim()>2){
    Histo=NULL;
    printf("\033[1;31mERROR:\033[0m DLM_HistoSource is broken (bad input)\n");
  }
  Histo = histo;
}
DLM_HistoSource::DLM_HistoSource(DLM_Histo<float>& histo):MyOwnHisto(false){
  if(histo.GetDim()>2){
    Histo=NULL;
    printf("\033[1;31mERROR:\033[0m DLM_HistoSource is broken (bad input)\n");
  }
  Histo = new DLM_Histo<float>(histo);
}
DLM_HistoSource::~DLM_HistoSource(){
  if(MyOwnHisto){
    delete Histo;
    Histo = NULL;
  }
}
double DLM_HistoSource::Eval(double* kxc){
  if(Histo->GetDim()==1){
    return Histo->Eval(&kxc[1]);
  }
  else{
    return Histo->Eval(kxc);
  }
}
double DLM_HistoSource::RootEval(double* x, double* kstar){
  double kxc[3];
  kxc[0] = *kstar;
  kxc[1] = *x;
  return Eval(kxc);
}




DLM_CecaSource_v0::DLM_CecaSource_v0(const std::string systype, const std::string anaver, const std::string infolder):
                            SystemType(systype),AnaVersion(anaver),InputFolder(infolder){
  ErrorState = false;
  dlmSource = NULL;
  AnaVersionBase = "";
  ErrorState = !InitHisto();
  //printf("%i\n", ErrorState);
  ErrorState = !LoadSource();
  //printf("%i\n", ErrorState);
}

DLM_CecaSource_v0::~DLM_CecaSource_v0(){
  if(dlmSource){delete dlmSource; dlmSource=NULL;}
}

bool DLM_CecaSource_v0::LoadSource(){
  if(ErrorState) return false;
  if(!dlmSource) return false;
  FILE *InFile;
  int NumPar=0;
  char* cline = new char [512];
  char* cdscr = new char [128];
  char* cval = new char [128];

  const unsigned NumMtBins = dlmSource?dlmSource->GetNbins(0):0;
  KdpPars* SrcPar = new KdpPars [NumMtBins];
  float* MtBinCenter = new float [NumMtBins];
  float* Chi2Ndf = new float [NumMtBins];

  double read_value;
  do{
    std::string InFileName = InputFolder+AnaVersionBase+"."+std::to_string(NumPar)+".ceca.source";
    InFile = fopen(InFileName.c_str(), "r");
    //printf("InFileName = %s\n",InFileName.c_str());
    //0 = the settings
    //1 = kdp
    int ReadMode = 0;
    if(InFile){
      int FitMode;

      float d_x = 0;
      float d_y = 0;
      float d_z = 0;
      float h_x = 0;
      float h_y = 0;
      float h_z = 0;
      float tau = 0;
      float del_proton = 0;
      float del_proton_reso = 0;
      float del_lambda = 0;
      float del_lambda_reso = 0;
      for(unsigned uMt=0; uMt<NumMtBins; uMt++){
        MtBinCenter[uMt] = 0;
        Chi2Ndf[uMt] = 0;
        //d_x[uMt] = 0;
        //d_y[uMt] = 0;
        //d_z[uMt] = 0;
        //h_x[uMt] = 0;
        //h_y[uMt] = 0;
        //h_z[uMt] = 0;
      }

      std::string type="";
      std::string AnaVersion="";
      unsigned ControlMt = 0;
      while(!feof(InFile)){
        if(!fgets(cline, 511, InFile)){
          //printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be properly read (%s)!\n", InputFileName.Data(),cline);
        }
        if(strncmp(cline,"SOURCE:",7)==0){
          ReadMode=1;
          continue;
        }

        if(ReadMode==0){
          sscanf(cline, "%s %s",cdscr,cval);
          if(strcmp(cdscr,"type")==0){
            type = std::string(cval);
          }
          else if(strcmp(cdscr,"AnaVersion")==0){
            AnaVersion = std::string(cval);
          }
          else{
            read_value = stod(cval);
            //      if(strcmp(cdscr,"GLOB_TIMEOUT")==0) {GLOB_TIMEOUT = read_value;}
            //else if(strcmp(cdscr,"multiplicity")==0) {multiplicity = unsigned(read_value);}
            //else if(strcmp(cdscr,"target_yield")==0) {target_yield = unsigned(read_value);}
            //else if(strcmp(cdscr,"femto_region")==0) {femto_region = read_value;}
            if(strcmp(cdscr,"d_x")==0) {d_x = read_value;}
            if(strcmp(cdscr,"d_y")==0) {d_y = read_value;}
            if(strcmp(cdscr,"d_z")==0) {d_z = read_value;}
            if(strcmp(cdscr,"h_x")==0) {h_x = read_value;}
            if(strcmp(cdscr,"h_y")==0) {h_y = read_value;}
            if(strcmp(cdscr,"h_z")==0) {h_z = read_value;}
            //else if(strcmp(cdscr,"h_fct")==0) {h_fct = read_value;}
            if(strcmp(cdscr,"tau")==0) {tau = read_value;}
            //else if(strcmp(cdscr,"tau_fct")==0) {tau_fct = read_value;}
            //else if(strcmp(cdscr,"tau_prp")==0) {tau_prp = bool(read_value);}
            //else if(strcmp(cdscr,"hdr_size")==0) {hdr_size = read_value;}
            //else if(strcmp(cdscr,"hdr_slope")==0) {hdr_slope = read_value;}
            //else if(strcmp(cdscr,"th_kick")==0) {th_kick = read_value;}
            //else if(strcmp(cdscr,"frag_beta")==0) {frag_beta = read_value;}
            //else if(strcmp(cdscr,"fixed_hdr")==0) {fixed_hdr = read_value;}
            //else if(strcmp(cdscr,"momdst_flag")==0) {momdst_flag = int(read_value);}
            //else if(strcmp(cdscr,"reso_flag")==0) {reso_flag = int(read_value);}
            //else if(strcmp(cdscr,"wildcard_flag")==0) {wildcard_flag = int(read_value);}
            //else if(strcmp(cdscr,"m_proton_reso")==0) {m_proton_reso = read_value;}
            //else if(strcmp(cdscr,"tau_proton_reso")==0) {tau_proton_reso = read_value;}
            //else if(strcmp(cdscr,"frac_proton_reso")==0) {frac_proton_reso = read_value;}
            //else if(strcmp(cdscr,"m_lambda_reso")==0) {m_lambda_reso = read_value;}
            //else if(strcmp(cdscr,"tau_lambda_reso")==0) {tau_lambda_reso = read_value;}
            //else if(strcmp(cdscr,"frac_lambda_reso")==0) {frac_lambda_reso = read_value;}
            if(strcmp(cdscr,"del_proton")==0) {del_proton = read_value;}
            if(strcmp(cdscr,"del_proton_reso")==0) {del_proton_reso = read_value;}
            if(strcmp(cdscr,"del_lambda")==0) {del_lambda = read_value;}
            if(strcmp(cdscr,"del_lambda_reso")==0) {del_lambda_reso = read_value;}
          }//strcmp(cdscr,"type")
        }//ReadMode==1
        else{
          if(ControlMt==NumMtBins){
            //printf("\033[1;33mWARNING!\033[0m Possible bad input, as the num of Mt bins does not check out (max = %u)\n",NumMtBins);
            continue;
          }
          else{
            sscanf(cline, "%i %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
            &FitMode,&MtBinCenter[ControlMt],
            &SrcPar[ControlMt].mean[0],&SrcPar[ControlMt].stdv[0],&SrcPar[ControlMt].wght[0],
            &SrcPar[ControlMt].mean[1],&SrcPar[ControlMt].stdv[1],&SrcPar[ControlMt].wght[1],
            &SrcPar[ControlMt].mean[2],&SrcPar[ControlMt].stdv[2],&SrcPar[ControlMt].wght[2],
            &SrcPar[ControlMt].mean[3],&SrcPar[ControlMt].stdv[3],&SrcPar[ControlMt].wght[3],
            &SrcPar[ControlMt].mean[4],&SrcPar[ControlMt].stdv[4],&SrcPar[ControlMt].wght[4],
            &SrcPar[ControlMt].mean[5],&SrcPar[ControlMt].stdv[5],&SrcPar[ControlMt].wght[5],
            &SrcPar[ControlMt].mean[6],&SrcPar[ControlMt].stdv[6],&SrcPar[ControlMt].wght[6],
            &SrcPar[ControlMt].mean[7],&SrcPar[ControlMt].stdv[7],&SrcPar[ControlMt].wght[7],
            &SrcPar[ControlMt].mean[8],&SrcPar[ControlMt].stdv[8],&SrcPar[ControlMt].wght[8],
            &SrcPar[ControlMt].mean[9],&SrcPar[ControlMt].stdv[9],
            &Chi2Ndf[ControlMt]);
            ControlMt++;
          }
        }//ReadMode
      }//while InFile

      if(type!=SystemType){
        if(InFile) fclose(InFile);
        NumPar++;
        continue;
      }
      if(d_x!=d_y || d_y!=d_z){
        ErrorState = true;
        printf("\033[1;31mERROR:\033[0m d_x!=d_y || d_y!=d_z\n");
        break;
      }
      if(h_x!=h_y){
        ErrorState = true;
        printf("\033[1;31mERROR:\033[0m h_x!=h_y\n");
        break;
      }

      for(unsigned uMt=0; uMt<NumMtBins; uMt++){
        unsigned WhichBin[4];
        WhichBin[0] = dlmSource->GetBin(0,MtBinCenter[uMt]);//mT
        if(AnaVersionBase=="Cigar2"||AnaVersionBase=="Jaime1"){
          WhichBin[1] = dlmSource->GetBin(1,d_x);//d
        }
        else if(AnaVersionBase=="JaimeDelay1"){
          if(del_lambda!=del_lambda_reso){
            printf("!!! del_lambda!=del_lambda_reso !!!\n");
          }
          WhichBin[1] = dlmSource->GetBin(1,del_lambda);//d
        }
        WhichBin[2] = dlmSource->GetBin(2,h_x);//ht
        if(AnaVersionBase=="Cigar2")
          WhichBin[3] = dlmSource->GetBin(3,h_z);//hz
        else if(AnaVersionBase=="Jaime1"||AnaVersionBase=="JaimeDelay1")
          WhichBin[3] = dlmSource->GetBin(3,tau);
  //if(fabs(MtBinCenter[uMt]-1110)<1 && fabs(d_x-0.1867)<0.001 && fabs(h_x-4.067)<0.01 && fabs(h_z-10.167)<0.01 ){
  //if(WhichBin[1]>0){
    //printf("file = %s\n",InFileName.c_str());
    //printf("[%u][%u][%u][%u]:\n",WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3]);
    //printf(" {%.2f}{%.2f}{%.2f}{%.2f}{%.2f}{%.2f}\n",MtBinCenter[uMt],d_x,h_x,h_z,tau,del_lambda);
    //for(unsigned uP=0; uP<10; uP++){
    //  printf("  P%u --> %.3f %.3f %.3f\n",uP,SrcPar[uMt].mean[uP],SrcPar[uMt].stdv[uP],SrcPar[uMt].wght[uP]);
    //}
    //usleep(10e3);
  //}


        if(MtBinCenter[uMt]<dlmSource->GetLowEdge(0) || MtBinCenter[uMt]>dlmSource->GetUpEdge(0)){
          printf("\033[1;33mWARNING:\033[0m MtBinCenter outside range\n");
          //if(InFile) fclose(InFile);
          //NumPar++;
          continue;
        }

        if(AnaVersionBase=="Cigar2"||AnaVersionBase=="Jaime1"){
          if(d_x<dlmSource->GetLowEdge(1) || d_x>dlmSource->GetUpEdge(1)){
            printf("\033[1;33mWARNING:\033[0m d_x outside range\n");
            //if(InFile) fclose(InFile);
            //NumPar++;
            continue;
          }
        }
        else if(AnaVersionBase=="JaimeDelay1"){
          if(del_lambda<dlmSource->GetLowEdge(1) || del_lambda>dlmSource->GetUpEdge(1)){
            printf("\033[1;33mWARNING:\033[0m del_lambda outside range\n");
            continue;
          }
        }

        if(h_x<dlmSource->GetLowEdge(2) || h_x>dlmSource->GetUpEdge(2)){
          printf("\033[1;33mWARNING:\033[0m h_x outside range\n");
          //if(InFile) fclose(InFile);
          //NumPar++;
          continue;
        }
        if(AnaVersionBase=="Cigar2"){
          if(h_z<dlmSource->GetLowEdge(3) || h_z>dlmSource->GetUpEdge(3)){
            printf("\033[1;33mWARNING:\033[0m h_z outside range\n");
            //if(InFile) fclose(InFile);
            //NumPar++;
            continue;
          }
        }
        else if(AnaVersionBase=="Jaime1"||AnaVersionBase=="JaimeDelay1"){
          if(tau<dlmSource->GetLowEdge(3) || tau>dlmSource->GetUpEdge(3)){
            printf("\033[1;33mWARNING:\033[0m h_z outside range\n");
            //if(InFile) fclose(InFile);
            //NumPar++;
            continue;
          }
        }


        KdpPars BinCont = dlmSource->GetBinContent(WhichBin);
        if(BinCont!=0){
          printf("\033[1;33mWARNING:\033[0m Attempt to refill a filled bin. This should not happen!\n");
          //if(InFile) fclose(InFile);
          //NumPar++;
          continue;
        }
        dlmSource->SetBinContent(WhichBin,SrcPar[uMt]);
        //printf("mT = %f\n",MtBinCenter[uMt]);
        //printf("d_x = %f\n",d_x);
        //printf("h_x = %f\n",h_x);
        //printf("tau = %f\n",tau);
        //printf("[%u][%u][%u][%u]\n",WhichBin[0],WhichBin[1],WhichBin[2],WhichBin[3]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[0],SrcPar[uMt].stdv[0],SrcPar[uMt].wght[0]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[1],SrcPar[uMt].stdv[1],SrcPar[uMt].wght[1]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[2],SrcPar[uMt].stdv[2],SrcPar[uMt].wght[2]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[3],SrcPar[uMt].stdv[3],SrcPar[uMt].wght[3]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[4],SrcPar[uMt].stdv[4],SrcPar[uMt].wght[4]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[5],SrcPar[uMt].stdv[5],SrcPar[uMt].wght[5]);
        //printf("%f %f %f\n",SrcPar[uMt].mean[6],SrcPar[uMt].stdv[6],SrcPar[uMt].wght[6]);
        //usleep(100e3);
      }//uMt

      NumPar++;
    }//if InFile
    //else{
      //printf("Issue with %s\n",InFileName.c_str());
    //}
    if(InFile) fclose(InFile);
  }
  while(InFile);
  if(InFile) fclose(InFile);

  delete [] cline;
  delete [] cdscr;
  delete [] cval;

  delete [] SrcPar;
  delete [] MtBinCenter;
  delete [] Chi2Ndf;

  return true;
}

bool DLM_CecaSource_v0::InitHisto(){
  if(dlmSource){
    printf("\033[1;31mERROR:\033[0m dlmSource exists when it should not!\n");
    return false;
  }
  //if it starts with that
  if(AnaVersion.rfind("Cigar2_ds",0)==0||AnaVersion.rfind("Jaime1_ds",0)==0||AnaVersion.rfind("JaimeDelay1_dLs",0)==0) { // pos=0 limits the search to the prefix
    if(AnaVersion.rfind("Cigar2_ds",0)==0) AnaVersionBase = "Cigar2";
    else if(AnaVersion.rfind("Jaime1_ds",0)==0) AnaVersionBase = "Jaime1";
    else if(AnaVersion.rfind("JaimeDelay1_dLs",0)==0) AnaVersionBase = "JaimeDelay1";
    vector<string> str_pars;
    std::string str_tmp;
    //used as the LambdaDelay in JaimeDelay1, otherwise the dist
    unsigned NumDist = 0;
    unsigned NumHt = 0;
    //can be used for Tau (Jaime*)
    unsigned NumHz = 0;
    if(SystemType!="pp" && SystemType!="pL"){
      printf("\033[1;31mERROR:\033[0m Unknown SystemType = %s!\n",SystemType.c_str());
      return false;
    }
    //std::string delim = "_";
    std::stringstream ssAnaVersion(AnaVersion);
    while(getline(ssAnaVersion,str_tmp,'_')){
      //printf("str_tmp = %s\n",str_tmp.c_str());
      if(str_tmp.rfind("ds",0)==0){
        //deletes the first two chars
        str_tmp.erase(0,2);
        NumDist = stoi(str_tmp);
      }
      else if(str_tmp.rfind("hts",0)==0){
        str_tmp.erase(0,3);
        NumHt = stoi(str_tmp);
      }
      else if(str_tmp.rfind("hzs",0)==0){
        str_tmp.erase(0,3);
        NumHz = stoi(str_tmp);
      }
      else if(str_tmp.rfind("taus",0)==0){
        str_tmp.erase(0,4);
        NumHz = stoi(str_tmp);
      }
      else if(str_tmp.rfind("dLs",0)==0){
        str_tmp.erase(0,3);
        //printf("%s\n",str_tmp.c_str());
        NumDist = stoi(str_tmp);
        //printf("ND=%u\n",NumDist);
      }
      else if(str_tmp.rfind("Cigar2",0)==0||str_tmp.rfind("Jaime1",0)==0||str_tmp.rfind("JaimeDelay1",0)==0){
        //do nothing
      }
      else{
        printf("\033[1;31mERROR:\033[0m DLM_CecaSource_v0::InitHisto something wrong with the pars of AnaVersion = %s!\n",AnaVersion.c_str());
        return false;
      }
    }//while

    //printf("AnaVersionBase = %s\n",AnaVersionBase.c_str());
    //printf("NumDist = %u\n",NumDist);
    //printf("NumHt = %u\n",NumHt);
    //printf("NumHz = %u\n",NumHz);
    //usleep(1000e3);

    dlmSource = new DLM_Histo<KdpPars>();
    //the 3 pars + mT = 4 dims
    dlmSource->SetUp(4);
    unsigned NumMtBins;
    //all these values are hard coded and currently fixed
    //from the Jaime parameters
    double* BinRange = NULL;
    double* BinCenter = NULL;
    if(SystemType=="pp"){
      NumMtBins = 10;

      BinRange = new double [NumMtBins+1];
      BinRange[0] = 930; //avg  983 ( 985)
      BinRange[1] = 1020;//avg 1054 (1055)
      BinRange[2] = 1080;//avg 1110 (1110)
      BinRange[3] = 1140;//avg 1168 (1170)
      BinRange[4] = 1200;//avg 1228 (1230)
      BinRange[5] = 1260;//avg 1315 (1315)
      BinRange[6] = 1380;//avg 1463 (1460)
      BinRange[7] = 1570;//avg 1681 (1680)
      BinRange[8] = 1840;//avg 1923 (1920)
      BinRange[9] = 2030;//avg 2303 (2300)
      BinRange[10] = 4500;

      BinCenter = new double [NumMtBins];
      BinCenter[0] = 983;
      BinCenter[1] = 1054;
      BinCenter[2] = 1110;
      BinCenter[3] = 1168;
      BinCenter[4] = 1228;
      BinCenter[5] = 1315;
      BinCenter[6] = 1463;
      BinCenter[7] = 1681;
      BinCenter[8] = 1923;
      BinCenter[9] = 2303;

      dlmSource->SetUp(0,NumMtBins,BinRange);
      //
    }
    else if(SystemType=="pL"){
      NumMtBins = 8;

      BinRange = new double [NumMtBins+1];
      BinRange[0] = 1000;//avg 1121 (1120)
      BinRange[1] = 1170;//avg 1210 (1210)
      BinRange[2] = 1250;//avg 1288 (1290)
      BinRange[3] = 1330;//avg 1377 (1380)
      BinRange[4] = 1430;//avg 1536 (1540)
      BinRange[5] = 1680;//avg 1753 (1750)
      BinRange[6] = 1840;//avg 1935 (1935)
      BinRange[7] = 2060;//avg 2334 (2330)
      BinRange[8] = 4800;

      BinCenter = new double [NumMtBins];
      BinCenter[0] = 1121;
      BinCenter[1] = 1210;
      BinCenter[2] = 1288;
      BinCenter[3] = 1377;
      BinCenter[4] = 1536;
      BinCenter[5] = 1753;
      BinCenter[6] = 1935;
      BinCenter[7] = 2334;

      dlmSource->SetUp(0,NumMtBins,BinRange);
    }
    if(!BinRange || !BinCenter){
      printf("\033[1;31mERROR:\033[0m DLM_CecaSource_v0::InitHisto has an impossible bug!\n");
      return false;
    }
    //ranges are fixed for the Cigar2 case
    if(AnaVersionBase=="Cigar2"||AnaVersionBase=="Jaime1"){
      dlmSource->SetUp(1,NumDist,0,1.28);
      dlmSource->SetUp(2,NumHt,0,4.8);
           if(AnaVersionBase=="Cigar2") dlmSource->SetUp(3,NumHz,0,12.0);
      else if(AnaVersionBase=="Jaime1") dlmSource->SetUp(3,NumHz,0,6.0);
    }
    //JaimeDelay1
    else if(AnaVersionBase=="JaimeDelay1"){
      if(SystemType=="pp"){
        if(NumDist!=1){
          printf("JaimeDelay1 pp = wtf\n");
        }
        dlmSource->SetUp(1,NumDist,-0.4,0.1);
      }
      else{
        //printf("NumDist=%u\n",NumDist);
        //usleep(1000e3);
        dlmSource->SetUp(1,NumDist,-0.4,0.1);
      }
      dlmSource->SetUp(2,NumHt,3.45,3.65);
           if(AnaVersionBase=="Cigar2") dlmSource->SetUp(3,NumHz,0,12.0);
      else if(AnaVersionBase=="Jaime1") dlmSource->SetUp(3,NumHz,0,6.0);
      else if(AnaVersionBase=="JaimeDelay1") dlmSource->SetUp(3,NumHz,2.55,2.79);
    }
    dlmSource->Initialize();

    dlmSource->SetBinCenter(0,BinCenter);

    delete [] BinRange;
    delete [] BinCenter;
  }//Cigar2
  else{
    printf("\033[1;31mERROR:\033[0m DLM_CecaSource_v0::InitHisto has an unknown AnaVersion = %s!\n",AnaVersion.c_str());
    return false;
  }
  return true;
}

unsigned DLM_CecaSource_v0::FindMtBin(double Mt){
  if(!dlmSource) return 0;
  if(Mt<0) return 0;
  return dlmSource->GetBin(0,Mt);
}
double DLM_CecaSource_v0::FindMt(unsigned uMt){
  if(!dlmSource) return 0;
  return dlmSource->GetBinCenter(0,uMt);
}

double DLM_CecaSource_v0::Eval(double* kxc){
  return RootEval(&kxc[1],&kxc[3]);
}
//the parameters are: mT, 3 source pars of your choice
//for the Cigar2 these are d,ht,hz (in that order)
//mT == -1 => Gauss
//[4] is the scale factor
double DLM_CecaSource_v0::RootEval(double* x, double* pars){
  //Gaussian source
  if(pars[0]==-1){
    return GaussSourceTF1(x,&pars[1]);
  }

  if(!dlmSource) {
    printf("\033[1;33mWARNING:\033[0m !dlmSource\n");
    return 0;
  }

  for(unsigned uPar=0; uPar<4; uPar++){
    if( pars[uPar]<dlmSource->GetLowEdge(uPar) || pars[uPar]>dlmSource->GetUpEdge(uPar) ){
      printf("\033[1;33mWARNING:\033[0m DLM_CecaSource_v0::RootEval pars[%u] = %.3e is outside the allowed range [%.3e, %.3e]!\n",
      uPar,pars[uPar],dlmSource->GetLowEdge(uPar),dlmSource->GetUpEdge(uPar));
      return 0;
    }
  }

  double xval = (*x)*pars[4];
  KdpPars SrcPars = dlmSource->Eval(pars);
  for(unsigned uP=0; uP<KdpPars::NumDistos-1; uP++){
    if(SrcPars.wght[uP]<0) SrcPars.wght[uP]=0;
    if(SrcPars.wght[uP]>1) SrcPars.wght[uP]=1;
  }
//printf("XXX\n");
//SrcPars.Print();
//printf(" --- (%f) %f\n",xval,PoissonSum(xval,SrcPars));
//usleep(250e3);

  return PoissonSum(xval,SrcPars);
}

double DLM_CecaSource_v0::Low_par(unsigned uP, bool bincenter){
  if(!dlmSource || uP>=dlmSource->GetDim()) return 0;
  if(bincenter) return dlmSource->GetBinCenter(uP,0);
  else return dlmSource->GetLowEdge(uP);
}
double DLM_CecaSource_v0::Up_par(unsigned uP, bool bincenter){
  if(!dlmSource || uP>=dlmSource->GetDim()) return 0;
  if(bincenter) return dlmSource->GetBinCenter(uP,dlmSource->GetNbins(uP)-1);
  else return dlmSource->GetUpEdge(uP);
}

unsigned DLM_CecaSource_v0::GetNbins(unsigned WhichPar){
  return dlmSource->GetNbins(WhichPar);
}
double* DLM_CecaSource_v0::GetBinRange(unsigned WhichPar){
  double* BinRange;
  BinRange = dlmSource->GetBinRange(WhichPar);
  return BinRange;
}
double* DLM_CecaSource_v0::GetBinCenters(unsigned WhichPar){
  double* BinCtr;
  BinCtr = dlmSource->GetBinCenters(WhichPar);
  return BinCtr;
}



DLM_MtKstar_KdpSource::DLM_MtKstar_KdpSource(DLM_Histo<KdpPars>& InputHisto){
    MyOwnCopy = true;
    dlmSource = new DLM_Histo<KdpPars> (InputHisto);
    DIM = dlmSource->GetDim();
}
DLM_MtKstar_KdpSource::DLM_MtKstar_KdpSource(DLM_Histo<KdpPars>* InputHisto){
    dlmSource = NULL;
    MyOwnCopy = false;
    if(!InputHisto){
        printf("\033[1;31mERROR:\033[0m  NULL pointer in the constructor of DLM_MtKstar_KdpSource\n");
        return;
    }
    if(InputHisto->GetDim()<2){
        printf("\033[1;31mERROR:\033[0m  Bad input in the constructor of DLM_MtKstar_KdpSource\n");
        return;
    }
    dlmSource = InputHisto;
    DIM = dlmSource->GetDim();
    printf("%p\n",dlmSource);
}

DLM_MtKstar_KdpSource::~DLM_MtKstar_KdpSource(){
  if(MyOwnCopy && dlmSource){delete dlmSource; dlmSource=NULL;}
}

//[0] = kstar
//[1] = rstar
//[2] = empty
//[3] = mt
//[4...] = extra patameters
double DLM_MtKstar_KdpSource::Eval(double* pars){
    double& kstar = pars[0];
    double& rstar = pars[1];
    double& Mt = pars[3];
    double* EvalAt = new double [DIM];
    EvalAt[0] = Mt;
    EvalAt[1] = kstar;
    for(unsigned uPar=2; uPar<DIM; uPar++){
        EvalAt[uPar] = pars[4+uPar-2];
    }

    //printf("EvalAt: ");
    //for(unsigned uPar=0; uPar<DIM; uPar++){
    //    printf("%.2f ", EvalAt[uPar]);
    //}
    //printf("\n");
    
    //double EvalAt[3];
    //EvalAt[0] = Mt;
    //EvalAt[1] = kstar;
    if(!dlmSource) {printf("ZERO %p\n",dlmSource); return 0;}
    KdpPars SrcPars = dlmSource->Eval(EvalAt);
    double Result = PoissonSum(rstar,SrcPars);
    //if(Result){
        //printf("Result(%.3f) = %.3e\n",rstar,Result);
    //}
    delete [] EvalAt;
    return Result;
}

//pars[0] = mt
//pars[1] = kstar
double DLM_MtKstar_KdpSource::RootEval(double* rstar, double* pars){
    //3 reserved for the variables
    //the rest depend on the DIM, however the kstar counts as var, thus the -1
    double PARS[3+DIM-1];
    PARS[0] = pars[1];//kstar
    PARS[1] = *rstar;
    PARS[2] = 0;
    PARS[3] = pars[0];//mt
    for(unsigned uPar=2; uPar<DIM; uPar++){
        PARS[4+uPar-2] = pars[uPar];
    }
    return Eval(PARS);
}
