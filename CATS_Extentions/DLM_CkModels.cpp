
#include "DLM_CkModels.h"

#include "TMath.h"
#include "TComplex.h"
#include "gsl_sf_dawson.h"

#include "CATSconstants.h"
#include "CATStools.h"

using namespace std;

double Flat_Residual(const double& Momentum, const double* SourcePar, const double* PotPar){
//printf("FR at %f\n", Momentum);
    return 1;
}

double LednickyAsInStar(const double& Momentum, const double& GaussR, const double& ScattLenSin, const double& EffRangeSin,
                        const double& Norm, const double& lambda, const double& ares, const double& RadRes){

    //const double FmToNu=5.067731237e-3;

    //const double Pi(3.141592653589793);

    if(GaussR!=GaussR){
        printf("\033[1;33mWARNING:\033[0m LednickyAsInStar got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR*FmToNu;
    const double sLen1 = ScattLenSin*FmToNu;
    const double eRan1 = EffRangeSin*FmToNu;

    double QMOM = 2*Momentum;

    double F1 = gsl_sf_dawson(QMOM*Radius)/(QMOM*Radius);
    double F2 = (1.-exp(-QMOM*QMOM*Radius*Radius))/(QMOM*Radius);

    complex<double> ScattAmplSin = pow(1./sLen1+0.5*eRan1*Momentum*Momentum-i*Momentum,-1.);

    return Norm*(1.+lambda*(-0.5*exp(-Radius*Radius*QMOM*QMOM)+0.25*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius))+
                     real(ScattAmplSin)*F1/(sqrt(Pi)*Radius)-imag(ScattAmplSin)*F2/Radius/2.)
          +ares*exp(-RadRes*RadRes*QMOM*QMOM)
          );

}

double GeneralLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const double& ScattLenTri, const double& EffRangeTri,
                       const bool& SinOnly, const bool& QS, const bool& InverseScatLen){
    //const double FmToNu=5.067731237e-3;
    //const std::complex<double> i(0,1);
    //const double Pi(3.141592653589793);
//if(ScattLenSin>211){
//printf("ScattLenSin = %e\n",ScattLenSin);
//}
    if(GaussR!=GaussR){
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR*FmToNu;
    const double IsLen1 = InverseScatLen?ScattLenSin/FmToNu:1./(ScattLenSin*FmToNu+1e-64);
    const double eRan1 = EffRangeSin*FmToNu;
    const double IsLen3 = InverseScatLen?ScattLenTri/FmToNu:1./(ScattLenTri*FmToNu+1e-64);
    const double eRan3 = EffRangeTri*FmToNu;

    double F1 = gsl_sf_dawson(2.*Momentum*Radius)/(2.*Momentum*Radius);
    double F2 = (1.-exp(-4.*Momentum*Momentum*Radius*Radius))/(2.*Momentum*Radius);
    complex<double> ScattAmplSin = pow(IsLen1+0.5*eRan1*Momentum*Momentum-i*Momentum,-1.);

    double CkValue = 0.;
    CkValue +=  0.5*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius))+
                2*real(ScattAmplSin)*F1/(sqrt(Pi)*Radius)-imag(ScattAmplSin)*F2/Radius;
    //so far this is the eq. for singlet only, w/o QS

    //if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    //since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if(!SinOnly){
        complex<double> ScattAmplTri = pow(IsLen3+0.5*eRan3*Momentum*Momentum-i*Momentum,-1.);
        CkValue +=  3*(
                    0.5*pow(abs(ScattAmplTri)/Radius,2)*(1-(eRan3)/(2*sqrt(Pi)*Radius))+
                    2*real(ScattAmplTri)*F1/(sqrt(Pi)*Radius)-imag(ScattAmplTri)*F2/Radius);
//printf("CkValue3 = %.3f\n");
        CkValue *= 0.25;
    }
    //if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if(QS){
        CkValue -= exp(-Radius*Radius*4.*Momentum*Momentum);
        CkValue *= 0.5;
    }

    CkValue += 1;

    return CkValue;
}

double GeneralLednicky(const double& Momentum, const double& GaussR,
                       const complex<double>& ScattLenSin, const double& EffRangeSin,
                       const complex<double>& ScattLenTri, const double& EffRangeTri,
                       const bool& SinOnly, const bool& QS, const bool& InverseScatLen){

    if(GaussR!=GaussR){
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR*FmToNu;
    const complex<double> IsLen1 = InverseScatLen?ScattLenSin/FmToNu:1./(ScattLenSin*FmToNu+1e-64);
    const double eRan1 = EffRangeSin*FmToNu;
    const complex<double> IsLen3 = InverseScatLen?ScattLenTri/FmToNu:1./(ScattLenTri*FmToNu+1e-64);
    const double eRan3 = EffRangeTri*FmToNu;

    double F1 = gsl_sf_dawson(2.*Momentum*Radius)/(2.*Momentum*Radius);
    double F2 = (1.-exp(-4.*Momentum*Momentum*Radius*Radius))/(2.*Momentum*Radius);
    complex<double> ScattAmplSin = pow(IsLen1+0.5*eRan1*Momentum*Momentum-i*Momentum,-1.);

    double CkValue = 0.;
    CkValue +=  0.5*pow(abs(ScattAmplSin)/Radius,2)*(1.-(eRan1)/(2*sqrt(Pi)*Radius))+
                2*real(ScattAmplSin)*F1/(sqrt(Pi)*Radius)-imag(ScattAmplSin)*F2/Radius;
    //so far this is the eq. for singlet only, w/o QS

    //if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    //since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if(!SinOnly){
        complex<double> ScattAmplTri = pow(IsLen3+0.5*eRan3*Momentum*Momentum-i*Momentum,-1.);
        CkValue +=  3*(
                    0.5*pow(abs(ScattAmplTri)/Radius,2)*(1-(eRan3)/(2*sqrt(Pi)*Radius))+
                    2*real(ScattAmplTri)*F1/(sqrt(Pi)*Radius)-imag(ScattAmplTri)*F2/Radius);
        CkValue *= 0.25;
    }
    //if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if(QS){
        CkValue -= exp(-Radius*Radius*4.*Momentum*Momentum);
        CkValue *= 0.5;
    }

    CkValue += 1;

    return CkValue;
}

double GeneralCoulombLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const bool& QS, const double& RedMass, const double& Q1Q2){
    //const double FmToNu=5.067731237e-3;
    //const std::complex<double> i(0,1);
    //const double Pi(3.141592653589793);

    if(GaussR!=GaussR){
        printf("\033[1;33mWARNING:\033[0m GeneralCoulombLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Momentum2 = Momentum*Momentum;
    const double Radius = GaussR*FmToNu;
    //const double Radius2 = Radius*Radius;
    const double Rho = Radius*Momentum;
    const double Rho2 = Rho*Rho;
    const double sLen1 = ScattLenSin*FmToNu;
    const double eRan1 = EffRangeSin*FmToNu;

    //double Dawson = gsl_sf_dawson(2.*Rho);
    //double D_2 = (1.-exp(-4.*Rho2));
    double eta = CoulombEta(Momentum,RedMass,Q1Q2);
    double A_c = CoulombPenetrationFactor(eta);

    complex<double> ScattAmplSin=pow(1./sLen1+0.5*eRan1*Momentum2-i*Momentum*A_c-2.*Momentum*eta*CoulombEuler(eta),-1.);

    double CkValue =    0.25*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius)+0.5*pow(A_c-1.,2)*(1.-exp(-4.*Rho2)))
                        +Momentum*real(ScattAmplSin)*gsl_sf_dawson(2.*Rho)/(2.*sqrt(Pi)*Rho2)
                        -Momentum*imag(ScattAmplSin)*(0.25*(1.-exp(-4.*Rho2))/Rho2+(A_c-1.)*cos(Rho)*exp(-Rho2));
    //double CkValue =    0.25*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius))
    //                    +Momentum*real(ScattAmplSin)*gsl_sf_dawson(2.*Rho)/(2.*sqrt(Pi)*Rho2)
    //                    -Momentum*imag(ScattAmplSin)*(0.25*(1.-exp(-4.*Rho2))/Rho2);

    if(QS){
        CkValue -= 0.5*exp(-4.*Rho2);
    }
    else{
        CkValue *= 2.;
        //CkValue += 0.5*exp(-4.*Rho2);
    }

    return A_c*(CkValue+1.);
}
double GeneralCoulombLednicky(const double& Momentum, const double& GaussR,
                       const double& ScattLenSin, const double& EffRangeSin,
                       const double& ScattLenTri, const double& EffRangeTri,
                       const bool& QS, const double& RedMass, const double& Q1Q2){
    return   0.25*GeneralCoulombLednicky(Momentum,GaussR,ScattLenSin,EffRangeSin,QS,RedMass,Q1Q2)
            +0.75*GeneralCoulombLednicky(Momentum,GaussR,ScattLenTri,EffRangeTri,QS,RedMass,Q1Q2);
}

//e.g. ΛΛ
//SourcePar[0] = Radius
//PotPar[0] = a0 for 1S0
//PotPar[1] = Reff for 1S0
//PotPar[2] = a0 for 3S1
//PotPar[3] = Reff for 3S1
double Lednicky_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],0,0,true,true,false);
}
double Lednicky_Identical_Singlet_InvScatLen(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],0,0,true,true,true);
}
double Lednicky_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],0,0,true,false);
}
double Lednicky_Singlet_InvScatLen(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],0,0,true,false,true);
}
double Lednicky_Identical_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],PotPar[2],PotPar[3],false,true);
}
double Lednicky_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],PotPar[2],PotPar[3],false,false);
}

double ComplexLednicky_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen(PotPar[0],PotPar[1]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen,PotPar[2],0,0,true,true,false);
}
double ComplexLednicky_Identical_Singlet_InvScatLen(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen(PotPar[0],PotPar[1]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen,PotPar[2],0,0,true,true,true);
}
double ComplexLednicky_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen(PotPar[0],PotPar[1]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen,PotPar[2],0,0,true,false);
}
double ComplexLednicky_Singlet_InvScatLen(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen(PotPar[0],PotPar[1]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen,PotPar[2],0,0,true,false,true);
}
double ComplexLednicky_Identical_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen1(PotPar[0],PotPar[1]);
    complex<double> ScatLen3(PotPar[3],PotPar[4]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen1,PotPar[2],ScatLen3,PotPar[5],false,true);
}
double ComplexLednicky_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    complex<double> ScatLen1(PotPar[0],PotPar[1]);
    complex<double> ScatLen3(PotPar[3],PotPar[4]);
    return GeneralLednicky(Momentum,SourcePar[0],ScatLen1,PotPar[2],ScatLen3,PotPar[5],false,false);
}

//SourcePar[0] = Radius
//PotPar[0] = a0 for 1S0
//PotPar[1] = Reff for 1S0
//PotPar[2] = RedMass
//PotPar[3] = Q1Q2
double LednickyCoulomb_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    //return CoulombPenetrationFactor(Momentum,PotPar[2],PotPar[3])*Lednicky_Singlet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],false,PotPar[2],PotPar[3]);
}

double LednickyCoulomb_Identical_Singlet(const double& Momentum, const double* SourcePar, const double* PotPar){
    //return CoulombPenetrationFactor(Momentum,PotPar[2],PotPar[3])*Lednicky_Identical_Singlet(Momentum,SourcePar,PotPar);
    //return Lednicky_Identical_Singlet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],true,PotPar[2],PotPar[3]);
}

//SourcePar[0] = Radius
//PotPar[0] = a0 for 1S0
//PotPar[1] = Reff for 1S0
//PotPar[2] = a0 for 3S1
//PotPar[3] = Reff for 3S1
//PotPar[4] = RedMass
//PotPar[5] = Q1Q2
double LednickyCoulomb_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    //return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5])*Lednicky_Triplet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],PotPar[2],PotPar[3],false,PotPar[4],PotPar[5]);
}

double LednickyCoulomb_Identical_Triplet(const double& Momentum, const double* SourcePar, const double* PotPar){
//return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5]);
    //return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5])*Lednicky_Identical_Triplet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],PotPar[2],PotPar[3],true,PotPar[4],PotPar[5]);
}

//e.g. pΛ
//SourcePar[0] = Radius
//PotPar[0] = a0 for 1S0
//PotPar[1] = Reff for 1S0
//PotPar[2] = a0 for 3S1
//PotPar[3] = Reff for 3S1
double Lednicky_SingletTriplet(const double& Momentum, const double* SourcePar, const double* PotPar){
    return GeneralLednicky(Momentum,SourcePar[0],PotPar[0],PotPar[1],PotPar[2],PotPar[3],false,false);
}
/*
//this code is copy-pasted from Oliver's analysis for his thesis
//the only parameter I use is par[0] = Source Size
//The input should be in MeV and fm
 double Lednicky_gauss_Sigma0(const double &Momentum, const double* SourcePar, const double* PotPar){
    //This model tries to calculate the Sigma0 correlation function, one can simply say trying because nothing is known about
    //this kind of interaction

    //Parameter definitions:
    //-----------------------------------------
    //-----------------------------------------
    //             isospin 1/2 (==0)
    //________________________________________
    const complex<double> scatt_length_1s0_0(-1.1,0.);
    const complex<double> scatt_length_3s1_0(-1.1,4.3);
    const complex<double> effrange_1s0_0(-1.5,0.);
    const complex<double> effrange_3s1_0(-2.2,-2.4);

    const double hbarc = 0.197326971844;
    //Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
    const double MomentumGeV = Momentum/1000.;

    //
    //             isospin 3/2 (==1)
    //________________________________________
    const complex<double> scatt_length_1s0_1(2.51,0.);
    const complex<double> scatt_length_3s1_1(-0.73,0.);
    const complex<double> effrange_1s0_1(4.92,0.);
    const complex<double> effrange_3s1_1(-1.22,0.);

    const double mass_proton = 0.938272;
    const double mass_neutron = 0.939565;
    const double mass_sigmaplus = 1.189377;
    const double mass_sigma0 = 1.192642;


    double mu1 = mass_proton*mass_sigma0/(mass_proton+mass_sigma0);
    double mu2 = mass_neutron*mass_sigmaplus/(mass_neutron+mass_sigmaplus);
    //-----------------------------------------
    //-----------------------------------------

    double k1 = MomentumGeV/hbarc;//momentum in the elastic channel
    double k2 = sqrt(2*mu2*(MomentumGeV*MomentumGeV/(2*mu1) + mass_proton + mass_sigma0 - mass_neutron - mass_sigmaplus));
    //sqrt(2*mu2*(MomentumGeV*MomentumGeV/(2*mu1) + mass_proton + mass_sigma0 - mass_neutron - mass_sigmaplus));
    k2 /=hbarc; //momentum in the inelastic channel
    //(k2>k1)

    double RG = SourcePar[0]; //Gaussian radius

    //Define inverse K-matrices for different spin configurations in isopin basis:


    //spin singlet:
    const complex<double> Kmin_1s0_0 = 1./scatt_length_1s0_0 + 0.5 * effrange_1s0_0 * k1*k1;
    const complex<double> Kmin_1s0_1 = 1./scatt_length_1s0_1 + 0.5 * effrange_1s0_1 * k1*k1;
    //spin triplet:
    const complex<double> Kmin_3s1_0 = 1./scatt_length_3s1_0 + 0.5 * effrange_3s1_0 * k1*k1;
    const complex<double> Kmin_3s1_1 = 1./scatt_length_3s1_1 + 0.5 * effrange_3s1_1 * k1*k1;


    //transform to transition basis with help of Clebsch-Gordan coefficients:
    const complex<double> Kmin_11_1s0 = 2./3. * Kmin_1s0_1 + 1./3. * Kmin_1s0_0;
    const complex<double> Kmin_11_3s1 = 2./3. * Kmin_3s1_1 + 1./3. * Kmin_3s1_0;

    const complex<double> Kmin_22_1s0 = 1./3. * Kmin_1s0_1 + 2./3. * Kmin_1s0_0;
    const complex<double> Kmin_22_3s1 = 1./3. * Kmin_3s1_1 + 2./3. * Kmin_3s1_0;

    const complex<double> Kmin_12_1s0 = sqrt(2.)/3. * (Kmin_1s0_1 - Kmin_1s0_0);
    const complex<double> Kmin_12_3s1 = sqrt(2.)/3. * (Kmin_3s1_1 - Kmin_3s1_0);


    //Determinant of scattering amplitude matrix:
    //is this correct, or should we have abs(Kmin_12_1s0)^2 ???
    complex<double> D_1s0 = (Kmin_11_1s0 + complex<double>(0.,-k1))*(Kmin_22_1s0 + complex<double>(0.,-k2)) - Kmin_12_1s0*Kmin_12_1s0;
    complex<double> D_3s1 = (Kmin_11_3s1 + complex<double>(0.,-k1))*(Kmin_22_3s1 + complex<double>(0.,-k2)) - Kmin_12_1s0*Kmin_12_1s0;


    //What is really needed are the scattering amplitudes:

    complex<double> f11_1s0 = (Kmin_22_1s0 + complex<double>(0.,-k2))/D_1s0;
    complex<double> f11_3s1 = (Kmin_22_3s1 + complex<double>(0.,-k2))/D_3s1;

    complex<double> f12_1s0 = -(Kmin_12_1s0)/D_1s0;
    complex<double> f12_3s1 = -(Kmin_12_3s1)/D_3s1;

    double fF1 = gsl_sf_dawson(2*k1*RG)/(2*k1*RG);
    double fF2 = (1.-exp(-pow(2*k1*RG,2)))/(2*k1*RG);

  //elastic part of the correlation function(1 -> 1):
  //-----------------------------------------
    double corr_1s0 = 0.5*pow(abs(f11_1s0),2.)/(RG*RG)* + 2.*real(f11_1s0)/(sqrt(Pi)*RG)*fF1 - imag(f11_1s0)/RG * fF2;
    corr_1s0 *= 0.25;


    double corr_3s1 = 0.5*pow(abs(f11_3s1),2.)/(RG*RG)* + 2.*real(f11_3s1)/(sqrt(Pi)*RG)*fF1 - imag(f11_3s1)/RG * fF2;
    corr_3s1 *= 0.75;
    //-----------------------------------------

    //inelastic part of the correlation function(2 -> 1):
    //-----------------------------------------
    double corr_1s0_inel = 0.5*mu2/mu1*fabs(pow(abs(f12_1s0),2.)/(RG*RG));
    corr_1s0_inel *= 0.25;

    double corr_3s1_inel = 0.5*mu2/mu1*fabs(pow(abs(f12_3s1),2.)/(RG*RG));
    corr_3s1_inel *= 0.75;

    //-----------------------------------------


    //correction term to non spheric distortions on the short range scale
    //Define matrices d_ij = 2Re d(K_ij)/dk²

    complex<double> dKmin_11_1s0_dk = 2./3. * 0.5 * effrange_1s0_1 + 1./3. * 0.5 * effrange_1s0_0;
    complex<double> dKmin_11_3s1_dk = 2./3. * 0.5 * effrange_3s1_1 + 1./3. * 0.5 * effrange_3s1_0;

    complex<double> dKmin_22_1s0_dk = 1./3. * 0.5 * effrange_1s0_1 + 2./3. * 0.5 * effrange_1s0_0;
    complex<double> dKmin_22_3s1_dk = 1./3. * 0.5 * effrange_3s1_1 + 2./3. * 0.5 * effrange_3s1_0;

    complex<double> dKmin_12_1s0_dk = sqrt(2.)/3. * 0.5 * (effrange_1s0_1 - effrange_1s0_0);
    complex<double> dKmin_12_3s1_dk = sqrt(2.)/3. * 0.5 * (effrange_3s1_1 - effrange_3s1_0);

    complex<double> factor1_1s0 = f11_1s0 * conj(f12_1s0);
    complex<double> factor1_3s1 = f11_3s1 * conj(f12_3s1);

    double corr_1s0_distort = pow(abs(f11_1s0),2.)* 2.*real(dKmin_11_1s0_dk) + pow(abs(f12_1s0),2.)* 2.*real(dKmin_22_1s0_dk) + 2.*real(factor1_1s0) * 2.*real(dKmin_12_1s0_dk);
    corr_1s0_distort *= -0.25/(4*sqrt(Pi)*pow(RG,3.));

    double corr_3s1_distort = pow(abs(f11_3s1),2.)* 2.*real(dKmin_11_3s1_dk) + pow(abs(f12_3s1),2.)* 2.*real(dKmin_22_3s1_dk) + 2.* real(factor1_3s1) * 2.*real(dKmin_12_3s1_dk);
    corr_3s1_distort *= -0.75/(4*sqrt(Pi)*pow(RG,3.));

    double corr_fin = 1. + corr_1s0 + corr_3s1 + corr_1s0_inel + corr_3s1_inel + corr_1s0_distort + corr_3s1_distort;

    return corr_fin;
    //if(par[1] == 0.) return 1. + par[2]*(par[3]*corr_fin -1.); //this is with baseline
    //else return par[2]*(par[3]*corr_fin -1.);//this is baseline subtracted
}
*/
//this code is copy-pasted from Oliver's analysis for his thesis
//the only parameter I use is par[0] = Source Size
//The input should be in MeV and fm
 double Lednicky_gauss_Sigma0(const double &Momentum, const double* SourcePar, const double* PotPar){
    //This model tries to calculate the Sigma0 correlation function, one can simply say trying because nothing is known about
    //this kind of interaction

    //Parameter definitions:
    //*****************************************
    //*****************************************
    //             isospin 1/2 (==0)
    //________________________________________
    const TComplex scatt_length_1s0_0(-1.1,0.);
    const TComplex scatt_length_3s1_0(-1.1,4.3);
    const TComplex effrange_1s0_0(-1.5,0.);
    const TComplex effrange_3s1_0(-2.2,-2.4);

    const double hbarc = 0.197326971844;
    //Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
    const double MomentumGeV = Momentum/1000.;

    //
    //             isospin 3/2 (==1)
    //________________________________________
    const TComplex scatt_length_1s0_1(2.51,0.);
    const TComplex scatt_length_3s1_1(-0.73,0.);
    const TComplex effrange_1s0_1(4.92,0.);
    const TComplex effrange_3s1_1(-1.22,0.);

    const Double_t mass_proton = 0.938272;
    const Double_t mass_neutron = 0.939565;
    const Double_t mass_sigmaplus = 1.189377;
    const Double_t mass_sigma0 = 1.192642;


    Double_t mu1 = mass_proton*mass_sigma0/(mass_proton+mass_sigma0);
    Double_t mu2 = mass_neutron*mass_sigmaplus/(mass_neutron+mass_sigmaplus);
    //*****************************************
    //*****************************************

    Double_t k1 = MomentumGeV/hbarc;//momentum in the elastic channel
    Double_t k2 = TMath::Sqrt(2*mu2*(MomentumGeV*MomentumGeV/(2*mu1) + mass_proton + mass_sigma0 - mass_neutron - mass_sigmaplus));
    k2 /=hbarc; //momentum in the inelastic channel
    //(k2>k1)

    Double_t RG = SourcePar[0]; //Gaussian radius

    //Define inverse K-matrices for different spin configurations in isopin basis:


    //spin singlet:
    TComplex Kmin_1s0_0 = 1./scatt_length_1s0_0 + 0.5 * effrange_1s0_0 * k1*k1;
    TComplex Kmin_1s0_1 = 1./scatt_length_1s0_1 + 0.5 * effrange_1s0_1 * k1*k1;
    //spin triplet:
    TComplex Kmin_3s1_0 = 1./scatt_length_3s1_0 + 0.5 * effrange_3s1_0 * k1*k1;
    TComplex Kmin_3s1_1 = 1./scatt_length_3s1_1 + 0.5 * effrange_3s1_1 * k1*k1;


    //transform to transition basis with help of Clebsch-Gordan coefficients:
    TComplex Kmin_11_1s0 = 2./3. * Kmin_1s0_1 + 1./3. * Kmin_1s0_0;
    TComplex Kmin_11_3s1 = 2./3. * Kmin_3s1_1 + 1./3. * Kmin_3s1_0;

    TComplex Kmin_22_1s0 = 1./3. * Kmin_1s0_1 + 2./3. * Kmin_1s0_0;
    TComplex Kmin_22_3s1 = 1./3. * Kmin_3s1_1 + 2./3. * Kmin_3s1_0;

    TComplex Kmin_12_1s0 = TMath::Sqrt(2.)/3. * (Kmin_1s0_1 - Kmin_1s0_0);
    TComplex Kmin_12_3s1 = TMath::Sqrt(2.)/3. * (Kmin_3s1_1 - Kmin_3s1_0);


    //Determinant of scattering amplitude matrix:
    TComplex D_1s0 = (Kmin_11_1s0 + TComplex(0.,-k1))*(Kmin_22_1s0 + TComplex(0.,-k2)) - TComplex::Power(Kmin_12_1s0,2.);
    TComplex D_3s1 = (Kmin_11_3s1 + TComplex(0.,-k1))*(Kmin_22_3s1 + TComplex(0.,-k2)) - TComplex::Power(Kmin_12_3s1,2.);


    //What is really needed are the scattering amplitudes:

    TComplex f11_1s0 = (Kmin_22_1s0 + TComplex(0.,-k2))/D_1s0;
    TComplex f11_3s1 = (Kmin_22_3s1 + TComplex(0.,-k2))/D_3s1;

    TComplex f12_1s0 = -(Kmin_12_1s0)/D_1s0;
    TComplex f12_3s1 = -(Kmin_12_3s1)/D_3s1;

    Double_t fF1 = gsl_sf_dawson(2*k1*RG)/(2*k1*RG);
    Double_t fF2 = (1.-TMath::Exp(-pow(2*k1*RG,2)))/(2*k1*RG);

  //elastic part of the correlation function(1 -> 1):
  //*********************************************
    Double_t corr_1s0 = 0.5*f11_1s0.Rho2()/(RG*RG)* + 2.*f11_1s0.Re()/(TMath::Sqrt(TMath::Pi())*RG)*fF1 - f11_1s0.Im()/RG * fF2;
    corr_1s0 *= 0.25;


    Double_t corr_3s1 = 0.5*f11_3s1.Rho2()/(RG*RG)* + 2.*f11_3s1.Re()/(TMath::Sqrt(TMath::Pi())*RG)*fF1 - f11_3s1.Im()/RG * fF2;
    corr_3s1 *= 0.75;
    //*********************************************

    //inelastic part of the correlation function(2 -> 1):
    //*********************************************
    Double_t corr_1s0_inel = 0.5*mu2/mu1*fabs(f12_1s0.Rho2()/(RG*RG));
    corr_1s0_inel *= 0.25;


    Double_t corr_3s1_inel = 0.5*mu2/mu1*fabs(f12_3s1.Rho2()/(RG*RG));
    corr_3s1_inel *= 0.75;

    //*********************************************


    //correction term to non spheric distortions on the short range scale
    //Define matrices d_ij = 2Re d(K_ij)/dk²

    TComplex dKmin_11_1s0_dk = 2./3. * 0.5 * effrange_1s0_1 + 1./3. * 0.5 * effrange_1s0_0;
    TComplex dKmin_11_3s1_dk = 2./3. * 0.5 * effrange_3s1_1 + 1./3. * 0.5 * effrange_3s1_0;

    TComplex dKmin_22_1s0_dk = 1./3. * 0.5 * effrange_1s0_1 + 2./3. * 0.5 * effrange_1s0_0;
    TComplex dKmin_22_3s1_dk = 1./3. * 0.5 * effrange_3s1_1 + 2./3. * 0.5 * effrange_3s1_0;

    TComplex dKmin_12_1s0_dk = TMath::Sqrt(2.)/3. * 0.5 * (effrange_1s0_1 - effrange_1s0_0);
    TComplex dKmin_12_3s1_dk = TMath::Sqrt(2.)/3. * 0.5 * (effrange_3s1_1 - effrange_3s1_0);

    TComplex factor1_1s0 = f11_1s0 * TComplex::Conjugate(f12_1s0);
    TComplex factor1_3s1 = f11_3s1 * TComplex::Conjugate(f12_3s1);

    Double_t corr_1s0_distort = f11_1s0.Rho2()* 2.*dKmin_11_1s0_dk.Re() + f12_1s0.Rho2()* 2.*dKmin_22_1s0_dk.Re() + 2.*factor1_1s0.Re() * 2.*dKmin_12_1s0_dk.Re();
    corr_1s0_distort *= -0.25/(4*TMath::Sqrt(TMath::Pi())*TMath::Power(RG,3.));

    Double_t corr_3s1_distort = f11_3s1.Rho2()* 2.*dKmin_11_3s1_dk.Re() + f12_3s1.Rho2()* 2.*dKmin_22_3s1_dk.Re() + 2.* factor1_3s1 * 2.*dKmin_12_3s1_dk.Re();
    corr_3s1_distort *= -0.75/(4*TMath::Sqrt(TMath::Pi())*TMath::Power(RG,3.));

    Double_t corr_fin = 1. + corr_1s0 + corr_3s1 + corr_1s0_inel + corr_3s1_inel + corr_1s0_distort + corr_3s1_distort;

    return corr_fin;
    //if(par[1] == 0.) return 1. + par[2]*(par[3]*corr_fin -1.); //this is with baseline
    //else return par[2]*(par[3]*corr_fin -1.);//this is baseline subtracted
}

//this code is copy-pasted from Oliver's analysis for his thesis
double pXi_pheno(const double &Momentum, const double* SourcePar, const double* PotPar){

    const double hbarc = 0.197326971844;
    const double MomentumGeV = Momentum/1000.;

    //phenomenological function developed for fitting p-Xi
    const double& RG = SourcePar[0];

    //very easy definition for a CF that peaks for k->0
    return 1. + exp(-MomentumGeV*RG/hbarc)/(MomentumGeV*RG/hbarc);

}
