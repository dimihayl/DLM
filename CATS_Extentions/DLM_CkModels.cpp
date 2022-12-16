
#include "DLM_CkModels.h"
#include "gsl_sf_dawson.h"
#include "CATSconstants.h"
#include "CATStools.h"
#include "DLM_Integration.h"
#include "DLM_Source.h"

using namespace std;

double (*LEDNI_INT_SFUNCTION)(double *par) = NULL;
double *LEDNI_INT_SFUNCTION_PARS = NULL;
void *LEDNI_INT_SCLASS = NULL;
unsigned LEDNI_INT_SCLASS_NPARS = 0;

void SetLedniIntegral_SourceFunction(double (*AS)(double *), double *par)
{
    LEDNI_INT_SFUNCTION = AS;
    LEDNI_INT_SFUNCTION_PARS = par;
    RemoveLedniIntegral_SourceClass();
}
void SetLedniIntegral_SourceClass(void *context, const unsigned &numparameters)
{
    LEDNI_INT_SCLASS = context;
    LEDNI_INT_SCLASS_NPARS = numparameters;
    RemoveLedniIntegral_SourceFunction();
}
void RemoveLedniIntegral_SourceFunction()
{
    LEDNI_INT_SFUNCTION = NULL;
    LEDNI_INT_SFUNCTION_PARS = NULL;
}
void RemoveLedniIntegral_SourceClass()
{
    LEDNI_INT_SCLASS = NULL;
    LEDNI_INT_SCLASS_NPARS = 0;
}

double Integral1_Function(double *par)
{
    double r = par[1];
    double k = par[0];
    double int1;
    int1 = cos(k * r) * sin(k * r) / (4 * Pi * r * r * k);
    //  printf("r in IntegralFunction1 *hcbar = %.4f \n", r*197.3);
    // printf("GaussSource = %.2f \n", GaussSource(par));
    // printf("multiplic = %.2f \n", int1*GaussSource(par));
    return int1 * GaussSourceCutOff(par);
}

double Integral2_Function(double *par)
{
    double r = par[1];
    double k = par[0];
    double int2;
    int2 = sin(k * r) * sin(k * r) / (4 * Pi * r * r * k);
    return int2 * GaussSourceCutOff(par);
}

double Integral1_FunctionGauss(double *par)
{
    double r = par[1];
    double k = par[0];
    double int1;
    int1 = cos(k * r) * sin(k * r) / (4 * Pi * r * r * k);
    //  printf("r in IntegralFunction1 *hcbar = %.4f \n", r*197.3);
    // printf("GaussSource = %.2f \n", GaussSource(par));
    // printf("multiplic = %.2f \n", int1*GaussSource(par));
    return int1 * GaussSource(par);
}

double Integral2_FunctionGauss(double *par)
{
    double r = par[1];
    double k = par[0];
    double int2;
    int2 = sin(k * r) * sin(k * r) / (4 * Pi * r * r * k);
    return int2 * GaussSource(par);
}

double Flat_Residual(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // printf("FR at %f\n", Momentum);
    return 1;
}

double LednickyAsInStar(const double &Momentum, const double &GaussR, const double &ScattLenSin, const double &EffRangeSin,
                        const double &Norm, const double &lambda, const double &ares, const double &RadRes)
{

    // const double FmToNu=5.067731237e-3;

    // const double Pi(3.141592653589793);

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m LednickyAsInStar got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const double sLen1 = ScattLenSin * FmToNu;
    const double eRan1 = EffRangeSin * FmToNu;

    double QMOM = 2 * Momentum;

    double F1 = gsl_sf_dawson(QMOM * Radius) / (QMOM * Radius);
    double F2 = (1. - exp(-QMOM * QMOM * Radius * Radius)) / (QMOM * Radius);

    complex<double> ScattAmplSin = pow(1. / sLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    return Norm * (1. + lambda * (-0.5 * exp(-Radius * Radius * QMOM * QMOM) + 0.25 * pow(abs(ScattAmplSin) / Radius, 2) * (1 - (eRan1) / (2 * sqrt(Pi) * Radius)) + real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius / 2.) + ares * exp(-RadRes * RadRes * QMOM * QMOM));
}

double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const double &ScattLenSin, const double &EffRangeSin,
                       const double &ScattLenTri, const double &EffRangeTri,
                       const bool &SinOnly, const bool &QS, const bool &InverseScatLen)
{
    // const double FmToNu=5.067731237e-3;
    // const std::complex<double> i(0,1);
    // const double Pi(3.141592653589793);
    // if(ScattLenSin>211){
    // printf("ScattLenSin = %e\n",ScattLenSin);
    // }
    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const double IsLen1 = InverseScatLen ? ScattLenSin / FmToNu : 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;
    const double IsLen3 = InverseScatLen ? ScattLenTri / FmToNu : 1. / (ScattLenTri * FmToNu + 1e-64);
    const double eRan3 = EffRangeTri * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1 - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    // so far this is the eq. for singlet only, w/o QS
    // if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    // since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if (!SinOnly)
    {
        complex<double> ScattAmplTri = pow(IsLen3 + 0.5 * eRan3 * Momentum * Momentum - i * Momentum, -1.);
        CkValue += 3 * (0.5 * pow(abs(ScattAmplTri) / Radius, 2) * (1 - (eRan3) / (2 * sqrt(Pi) * Radius)) +
                        2 * real(ScattAmplTri) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplTri) * F2 / Radius);
        // printf("CkValue3 = %.3f\n");
        CkValue *= 0.25;
    }
    // if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if (QS)
    {
        CkValue -= exp(-Radius * Radius * 4. * Momentum * Momentum);
        CkValue *= 0.5;
    }
    CkValue += 1;

    return CkValue;
}

// Weight1 should be the one corresponding to bose einstein statistics
double GeneralLednicky2channel(const double &Momentum, const double &GaussR,
                               const double &ScattLenSin, const double &EffRangeSin,
                               const double &ScattLenTri, const double &EffRangeTri,
                               const bool &QS, const bool &InverseScatLen,
                               const double &Weight1, const double &Weight2)
{
    // const double FmToNu=5.067731237e-3;
    // const std::complex<double> i(0,1);
    // const double Pi(3.141592653589793);
    // if(ScattLenSin>211){
    // printf("ScattLenSin = %e\n",ScattLenSin);
    // }
    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const double IsLen1 = InverseScatLen ? ScattLenSin / FmToNu : 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;
    const double IsLen3 = InverseScatLen ? ScattLenTri / FmToNu : 1. / (ScattLenTri * FmToNu + 1e-64);
    const double eRan3 = EffRangeTri * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    CkValue += Weight1 * (0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1 - (eRan1) / (2 * sqrt(Pi) * Radius)) +
                          2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius);
    // so far this is the eq. for singlet only, w/o QS
    // if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    // since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    complex<double> ScattAmplTri = pow(IsLen3 + 0.5 * eRan3 * Momentum * Momentum - i * Momentum, -1.);
    CkValue += Weight2 * (0.5 * pow(abs(ScattAmplTri) / Radius, 2) * (1 - (eRan3) / (2 * sqrt(Pi) * Radius)) +
                          2 * real(ScattAmplTri) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplTri) * F2 / Radius);
    // if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if (QS)
    {
        CkValue -= (exp(-Radius * Radius * 4. * Momentum * Momentum)) * (Weight2 - Weight1);
    }
    CkValue += 1;

    return CkValue;
}

double GeneralLednicky(const double &Momentum, const double &GaussR,
                       const complex<double> &ScattLenSin, const double &EffRangeSin,
                       const complex<double> &ScattLenTri, const double &EffRangeTri,
                       const bool &SinOnly, const bool &QS, const bool &InverseScatLen)
{

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;
    const complex<double> IsLen1 = InverseScatLen ? ScattLenSin / FmToNu : 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;
    const complex<double> IsLen3 = InverseScatLen ? ScattLenTri / FmToNu : 1. / (ScattLenTri * FmToNu + 1e-64);
    const double eRan3 = EffRangeTri * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);
    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Momentum * Momentum - i * Momentum, -1.);

    double CkValue = 0.;
    // double term1 = +2*real(ScattAmplSin)*F1/(sqrt(Pi)*Radius);
    // double term2 = -imag(ScattAmplSin)*F2/Radius;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               2 * real(ScattAmplSin) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplSin) * F2 / Radius;
    // so far this is the eq. for singlet only, w/o QS
    //  std::cout << "------- In General Lednicky analytical -------" << std::endl;
    //  printf("term1 = %.4f \n", term1);
    //  printf("term2 = %.4f \n", term2);
    //  std::cout << "----------------------------------------------" << std::endl;

    // if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    // since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if (!SinOnly)
    {
        complex<double> ScattAmplTri = pow(IsLen3 + 0.5 * eRan3 * Momentum * Momentum - i * Momentum, -1.);
        CkValue += 3 * (0.5 * pow(abs(ScattAmplTri) / Radius, 2) * (1 - (eRan3) / (2 * sqrt(Pi) * Radius)) +
                        2 * real(ScattAmplTri) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmplTri) * F2 / Radius);
        CkValue *= 0.25;
    }
    // if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if (QS)
    {
        CkValue -= exp(-Radius * Radius * 4. * Momentum * Momentum);
        CkValue *= 0.5;
    }

    CkValue += 1;

    return CkValue;
}

double GeneralLednickyIntegral(const double &Momentum, const double &GaussR, const double &cutoff,
                               const complex<double> &ScattLenSin, const double &EffRangeSin,
                               const complex<double> &ScattLenTri, const double &EffRangeTri,
                               const bool &SinOnly, const bool &QS, const bool &InverseScatLen)
{
    double params[5];
    params[0] = Momentum;
    params[3] = GaussR * FmToNu;
    params[4] = cutoff * FmToNu;
    // printf("params[0] = %.4f \n", params[0]);
    // printf("params[3] = %.4f \n", params[3]);
    // printf("cutoff = %.4f \n", cutoff);

    double intvalue1;
    double intvalue2;
    double max_int = 20. * FmToNu;
    double min_int = params[4];

    DLM_INT_SetFunction(GaussSourceCutOff, params, 1);
    // DLM_INT_SetFunction(GaussSource, params,1);

    // double NormSource = DLM_INT_aSimpsonWiki(1.e-3,max_int);
    double NormSource = DLM_INT_aSimpsonWiki(min_int, max_int);
    //  double NormSource = DLM_INT_SimpsonWiki(min_int + 1.e-3,max_int,10000);

    printf("NormSource = %.6f \n", NormSource);

    if (NormSource < 1.e-3)
    {
        printf("WARNING: the Normalization of the source is zero!!\n");
        return 1;
    }

    DLM_INT_SetFunction(Integral1_Function, params, 1);
    // DLM_INT_SetFunction(Integral1_FunctionGauss, params,1);
    // intvalue1 = (1./NormSource)*DLM_INT_aSimpsonWiki(1.e-3,max_int);
    intvalue1 = (1. / NormSource) * DLM_INT_aSimpsonWiki(min_int, max_int);
    // intvalue1 = (1./NormSource)*DLM_INT_SimpsonWiki(min_int + 1.e-3,max_int,10000);

    DLM_INT_SetFunction(Integral2_Function, params, 1);
    // DLM_INT_SetFunction(Integral2_FunctionGauss, params,1);

    // intvalue2 = (1./NormSource)*DLM_INT_aSimpsonWiki(1.e-3,max_int);
    intvalue2 = (1. / NormSource) * DLM_INT_aSimpsonWiki(min_int, max_int);
    // intvalue2 = (1./NormSource)*DLM_INT_SimpsonWiki(min_int + 1.e-3,max_int,10000);

    double Mom = params[0];
    double Radius = params[3];

    const complex<double> IsLen1 = InverseScatLen ? ScattLenSin / FmToNu : 1. / (ScattLenSin * FmToNu + 1e-64);
    const double eRan1 = EffRangeSin * FmToNu;
    // const complex<double> IsLen3 = InverseScatLen?ScattLenTri/FmToNu:1./(ScattLenTri*FmToNu+1e-64);
    // const double eRan3 = EffRangeTri*FmToNu;

    complex<double> ScattAmplSin = pow(IsLen1 + 0.5 * eRan1 * Mom * Mom - i * Mom, -1.);

    double CkValue = 0.;
    // double term1 = 8*Pi*real(ScattAmplSin)*intvalue1;
    // double term2 = -8*Pi*imag(ScattAmplSin)*intvalue2;
    CkValue += 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
               8 * Pi * real(ScattAmplSin) * intvalue1 - 8 * Pi * imag(ScattAmplSin) * intvalue2;

    if (!SinOnly)
    {
        // complex<double> ScattAmplTri = pow(IsLen3+0.5*eRan3*Mom*Mom-i*Mom,-1.);
        CkValue += 3 * (0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1. - (eRan1) / (2 * sqrt(Pi) * Radius)) +
                        8 * Pi * real(ScattAmplSin) * intvalue1 - 8 * Pi * imag(ScattAmplSin) * intvalue2);
        CkValue *= 0.25;
    }
    // if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if (QS)
    {
        CkValue -= exp(-Radius * Radius * 4. * Mom * Mom);
        CkValue *= 0.5;
    }

    CkValue += 1;

    return CkValue;
}

double GeneralCoulombLednicky(const double &Momentum, const double &GaussR,
                              const double &ScattLenSin, const double &EffRangeSin,
                              const bool &QS, const double &RedMass, const double &Q1Q2)
{
    // const double FmToNu=5.067731237e-3;
    // const std::complex<double> i(0,1);
    // const double Pi(3.141592653589793);

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralCoulombLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Momentum2 = Momentum * Momentum;
    const double Radius = GaussR * FmToNu;
    // const double Radius2 = Radius*Radius;
    const double Rho = Radius * Momentum;
    const double Rho2 = Rho * Rho;
    const double sLen1 = ScattLenSin * FmToNu;
    const double eRan1 = EffRangeSin * FmToNu;

    // double Dawson = gsl_sf_dawson(2.*Rho);
    // double D_2 = (1.-exp(-4.*Rho2));
    double eta = CoulombEta(Momentum, RedMass, Q1Q2);
    double A_c = CoulombPenetrationFactor(eta);

    complex<double> ScattAmplSin = pow(1. / sLen1 + 0.5 * eRan1 * Momentum2 - i * Momentum * A_c - 2. * Momentum * eta * CoulombEuler(eta), -1.);

    double CkValue = 0.25 * pow(abs(ScattAmplSin) / Radius, 2) * (1 - (eRan1) / (2 * sqrt(Pi) * Radius) + 0.5 * pow(A_c - 1., 2) * (1. - exp(-4. * Rho2))) + Momentum * real(ScattAmplSin) * gsl_sf_dawson(2. * Rho) / (2. * sqrt(Pi) * Rho2) - Momentum * imag(ScattAmplSin) * (0.25 * (1. - exp(-4. * Rho2)) / Rho2 + (A_c - 1.) * cos(Rho) * exp(-Rho2));
    // double CkValue =    0.25*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius))
    //                     +Momentum*real(ScattAmplSin)*gsl_sf_dawson(2.*Rho)/(2.*sqrt(Pi)*Rho2)
    //                     -Momentum*imag(ScattAmplSin)*(0.25*(1.-exp(-4.*Rho2))/Rho2);

    if (QS)
    {
        CkValue -= 0.5 * exp(-4. * Rho2);
    }
    else
    {
        CkValue *= 2.;
        // CkValue += 0.5*exp(-4.*Rho2);
    }

    return A_c * (CkValue + 1.);
}
double GeneralCoulombLednicky(const double &Momentum, const double &GaussR,
                              const double &ScattLenSin, const double &EffRangeSin,
                              const double &ScattLenTri, const double &EffRangeTri,
                              const bool &QS, const double &RedMass, const double &Q1Q2)
{
    return 0.25 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenSin, EffRangeSin, QS, RedMass, Q1Q2) + 0.75 * GeneralCoulombLednicky(Momentum, GaussR, ScattLenTri, EffRangeTri, QS, RedMass, Q1Q2);
}

// General Lednicky strong + Coulomb for spin averaged scattering parameters

double GeneralCoulombLednickyAvg(const double &Momentum, const double &GaussR,
                                 const std::complex<double> &ScattLenSin, const double &EffRangeSin,
                                 const bool &QS, const double &RedMass, const double &Q1Q2)
{
    // const double FmToNu=5.067731237e-3;
    // const std::complex<double> i(0,1);
    // const double Pi(3.141592653589793);

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralCoulombLednicky got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Momentum2 = Momentum * Momentum;
    const double Radius = GaussR * FmToNu;
    // const double Radius2 = Radius*Radius;
    const double Rho = Radius * Momentum;
    const double Rho2 = Rho * Rho;
    const complex<double> sLen1 = ScattLenSin * FmToNu; //(Real, Im)
    const double eRan1 = EffRangeSin * FmToNu;

    // double Dawson = gsl_sf_dawson(2.*Rho);
    // double D_2 = (1.-exp(-4.*Rho2));
    double eta = CoulombEta(Momentum, RedMass, Q1Q2);
    double A_c = CoulombPenetrationFactor(eta);

    complex<double> ScattAmplSin = pow(1. / sLen1 + 0.5 * eRan1 * Momentum2 - i * Momentum * A_c - 2. * Momentum * eta * CoulombEuler(eta), -1.);

    double CkValue = 0.5 * pow(abs(ScattAmplSin) / Radius, 2) * (1 - (eRan1) / (2 * sqrt(Pi) * Radius) + 0.5 * pow(A_c - 1., 2) * (1. - exp(-4. * Rho2))) + Momentum * real(ScattAmplSin) * gsl_sf_dawson(2. * Rho) / (2. * sqrt(Pi) * Rho2) - Momentum * imag(ScattAmplSin) * (0.25 * (1. - exp(-4. * Rho2)) / Rho2 + (A_c - 1.) * cos(Rho) * exp(-Rho2));
    // double CkValue =    0.25*pow(abs(ScattAmplSin)/Radius,2)*(1-(eRan1)/(2*sqrt(Pi)*Radius))
    //                     +Momentum*real(ScattAmplSin)*gsl_sf_dawson(2.*Rho)/(2.*sqrt(Pi)*Rho2)
    //                     -Momentum*imag(ScattAmplSin)*(0.25*(1.-exp(-4.*Rho2))/Rho2);

    if (QS)
    {
        CkValue -= 0.5 * exp(-4. * Rho2);
    }
    else
    {
        CkValue *= 2.;
        // CkValue += 0.5*exp(-4.*Rho2);
    }

    return A_c * (CkValue + 1.);
}

// e.g. ΛΛ
// SourcePar[0] = Radius
// PotPar[0] = a0 for 1S0
// PotPar[1] = Reff for 1S0
// PotPar[2] = a0 for 3S1
// PotPar[3] = Reff for 3S1
double Lednicky_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], 0, 0, true, true, false);
}
double Lednicky_Identical_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], 0, 0, true, true, true);
}
double Lednicky_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], 0, 0, true, false);
}
double Lednicky_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], 0, 0, true, false, true);
}
double Lednicky_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], false, true);
}
double Lednicky_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], false, false);
}

double ComplexLednicky_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, true, false);
}
double ComplexLednicky_Identical_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, true, true);
}
double ComplexLednicky_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false, false);
}
double ComplexLednicky_Singlet_InvScatLen(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false, true);
}
double ComplexLednicky_Singlet_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return SourcePar[3] * (SourcePar[2] * GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false, false) + (1 - SourcePar[2]) * GeneralLednicky(Momentum, SourcePar[1], ScatLen, PotPar[2], 0, 0, true, false, false)) + 1. - SourcePar[3];
}
double ComplexLednicky_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen1(PotPar[0], PotPar[1]);
    complex<double> ScatLen3(PotPar[3], PotPar[4]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen1, PotPar[2], ScatLen3, PotPar[5], false, true);
}
double ComplexLednicky_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen1(PotPar[0], PotPar[1]);
    complex<double> ScatLen3(PotPar[3], PotPar[4]);
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen1, PotPar[2], ScatLen3, PotPar[5], false, false);
}
double Lednicky_2channel(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky2channel(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], false, false, PotPar[4], PotPar[5]);
}

double Lednicky_gauss_pAL_v2(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    // Entries for the GeneralLednicky function:
    //  1. Momentum
    //  2. Radius of Gaussian Source
    //  3. Scattering length f0 Singlet:  Re, Im
    //  4. Effective range (assumed Re)
    //  5. Sc. Len. Triplet
    //  6. Eff. range Triplet
    //  7. TRUE to select singlet only: here the WUT sc. lengths are spin averaged
    //  8. to select Quantum Statistics (false for different fermions)
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false);
}

double Lednicky_gauss_LAL_v2(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    // Entries for the GeneralLednicky function:
    // 1. Momentum
    // 2. Radius of Gaussian Source
    // 3. Scattering length f0 Singlet:  Re, Im
    // 4. Effective range (assumed Re)
    // 5. Sc. Len. Triplet
    // 6. Eff. range Triplet
    // 7. TRUE to select singlet only: here the WUT sc. lengths are spin averaged
    // 8. to select Quantum Statistics (false for different fermions)
    return GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false);
}
double Lednicky_gauss_pAL_Integral(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    // Entries for the GeneralLednicky function:
    //  1. Momentum
    //  2. Radius of Gaussian Source
    //  3. Cutoff in the source
    //  4. Scattering length f0 Singlet:  Re, Im
    //  5. Effective range (assumed Re)
    //  6. Sc. Len. Triplet
    //  7. Eff. range Triplet
    //  8. TRUE to select singlet only: here the WUT sc. lengths are spin averaged
    //  9. to select Quantum Statistics (false for different fermions)
    return GeneralLednickyIntegral(Momentum, SourcePar[0], SourcePar[1], ScatLen, PotPar[2], 0, 0, true, false, false);
}
//
// double Lednicky_LAL_mess(const double& Momentum, const double* SourcePar, const double* PotPar){
//     std::complex<double> ScLenLAL(PotPar[0], PotPar[1]);
//     return Lednicky_gauss_LAL(Momentum, SourcePar, ScLenLAL, PotPar[2]);
// }

// SourcePar[0] = Radius
// PotPar[0] = a0 for 1S0
// PotPar[1] = Reff for 1S0
// PotPar[2] = RedMass
// PotPar[3] = Q1Q2
double LednickyCoulomb_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // return CoulombPenetrationFactor(Momentum,PotPar[2],PotPar[3])*Lednicky_Singlet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], false, PotPar[2], PotPar[3]);
}

double LednickyCoulomb_Identical_Singlet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // return CoulombPenetrationFactor(Momentum,PotPar[2],PotPar[3])*Lednicky_Identical_Singlet(Momentum,SourcePar,PotPar);
    // return Lednicky_Identical_Singlet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], true, PotPar[2], PotPar[3]);
}
// SourcePar[0] = Radius
// PotPar[0] = Real a_avg
// PotPar[1] = Imag a_avg
// PotPar[2] = Reff (assumed Real for the moment!!!!!)
// PotPar[3] = RedMass
// PotPar[4] = Q1Q2

double ComplexLednickyCoulomb_Averaged(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    return GeneralCoulombLednickyAvg(Momentum, SourcePar[0], ScatLen, PotPar[2], false, PotPar[3], PotPar[4]);
}

// Lednicky model for Baryon-Antibaryon analysis with APPROX Coulomb

double Lednicky_gauss_pAp_CATS(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    double radius = 1.188;
    std::complex<double> ScLenpAp = -0.894 + i * 0.88;
    double effrangepAp = 1.0;
    double mprot = 938.;
    double redmass = 0.5 * mprot;
    double c1 = +1.;
    double c2 = -1.;
    double charge = abs(c1 * c2);
    return GeneralCoulombLednickyAvg(Momentum, radius, ScLenpAp, effrangepAp, false, redmass, charge);
}

double Lednicky_gauss_pAp_test(const double &Momentum, const double *SourcePar, const double *PotPar)
{

    double radius = 1.188;
    std::complex<double> ScLenpAp = -0.894 + i * 0.88;
    double effrangepAp = 1.0;
    double mprot = 938.;
    double redmass = 0.5 * mprot;
    double c1 = +1.;
    double c2 = -1.;
    double charge = abs(c1 * c2);

    double eta = CoulombEta(Momentum, redmass, charge);
    double A_c = CoulombPenetrationFactor(eta);
    double Momentum2 = Momentum * Momentum;

    // Need to trasform all the scattering parameters in MeV
    const std::complex<double> sLen = ScLenpAp * FmToNu; //(Real, Im)
    double eRan = effrangepAp * FmToNu;
    // double rho = Momentum*radius;

    std::complex<double> ScattAmplCoul = pow(1. / sLen + 0.5 * eRan * Momentum2 - i * Momentum * A_c - 2. * Momentum * eta * CoulombEuler(eta), -1.);

    double F1 = gsl_sf_dawson(2. * Momentum * radius) / (2. * Momentum * radius);
    double F2 = (1. - std::exp(-4. * Momentum * Momentum * radius * radius)) /
                (2. * Momentum * radius);

    double CkValue = 0.;
    CkValue +=
        0.5 * std::pow(std::abs(ScattAmplCoul) / radius, 2) *
            (1. -
             (eRan) / (2 * std::sqrt(Pi) * radius)) +
        2 * std::real(ScattAmplCoul) * F1 / (std::sqrt(Pi) * radius) -
        std::imag(ScattAmplCoul) * F2 / radius;

    return CkValue + 1;

    // return GeneralCoulombLednickyAvg(Momentum,radius,ScLenpAp,effrangepAp,false,redmass,charge);
}

double Lednicky_gauss_pAL(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    //  double Lednicky_gauss_pAL(const double &Momentum,  const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;

    std::complex<double> ScLenpAL = -1.15 + i * 0.53;
    double effrangepAL = 3.06;
    // double mprot = 938.;
    // double mlam = 1115.;
    // double redmass = (mprot*mlam)/(mprot+mlam);
    // double c1 = +1.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);

    // return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenpAL, effrangepAL, 0, 0, true, false);
}

double Lednicky_gauss_LAL(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // double Lednicky_gauss_LAL(const double &Momentum, const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;

    std::complex<double> ScLenLAL = -0.9 + i * 0.4;
    double effrangeLAL = 2.76;
    // double mlam = 1115.;
    // double redmass = 0.5*mlam;
    // double c1 = +0.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);

    // return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenLAL, effrangeLAL, 0, 0, true, false);
}

double Lednicky_gauss_pAL_varup(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    //  double Lednicky_gauss_pAL(const double &Momentum,  const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;
    std::complex<double> ScLenpAL = -1.10 + i * 0.53;
    double effrangepAL = 3.06;
    // double mprot = 938.;
    // double mlam = 1115.;
    // double redmass = (mprot*mlam)/(mprot+mlam);
    // double c1 = +1.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);
    //  return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenpAL, effrangepAL, 0, 0, true, false);
}

double Lednicky_gauss_pAL_varlow(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    //  double Lednicky_gauss_pAL(const double &Momentum,  const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;
    std::complex<double> ScLenpAL = -1.20 + i * 0.53;
    double effrangepAL = 3.06;
    // double mprot = 938.;
    // double mlam = 1115.;
    // double redmass = (mprot*mlam)/(mprot+mlam);
    // double c1 = +1.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);
    //  return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenpAL, effrangepAL, 0, 0, true, false);
}

double Lednicky_gauss_LAL_varup(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // double Lednicky_gauss_LAL(const double &Momentum, const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;
    std::complex<double> ScLenLAL = -0.86 + i * 0.46;
    double effrangeLAL = 2.48;
    // double mlam = 1115.;
    // double redmass = 0.5*mlam;
    // double c1 = +0.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);
    //  return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenLAL, effrangeLAL, 0, 0, true, false);
}

double Lednicky_gauss_LAL_varlow(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // double Lednicky_gauss_LAL(const double &Momentum, const double* SourcePar, const std::complex<double>& ScattLenSin, const double& EffRangeSin){

    // double radius = 1.188;
    std::complex<double> ScLenLAL = -0.94 + i * 0.34;
    double effrangeLAL = 3.05;
    // double mlam = 1115.;
    // double redmass = 0.5*mlam;
    // double c1 = +0.;
    // double c2 = +0.;
    // double charge = abs(c1*c2);
    //  return GeneralLednicky(Momentum,SourcePar[0],ScattLenSin,EffRangeSin,0,0,true,false);
    return GeneralLednicky(Momentum, SourcePar[0], ScLenLAL, effrangeLAL, 0, 0, true, false);
}

// SourcePar[0] = Radius
// PotPar[0] = a0 for 1S0
// PotPar[1] = Reff for 1S0
// PotPar[2] = a0 for 3S1
// PotPar[3] = Reff for 3S1
// PotPar[4] = RedMass
// PotPar[5] = Q1Q2
double LednickyCoulomb_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5])*Lednicky_Triplet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], false, PotPar[4], PotPar[5]);
}

double LednickyCoulomb_Identical_Triplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5]);
    // return CoulombPenetrationFactor(Momentum,PotPar[4],PotPar[5])*Lednicky_Identical_Triplet(Momentum,SourcePar,PotPar);
    return GeneralCoulombLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], true, PotPar[4], PotPar[5]);
}

// e.g. pΛ
// SourcePar[0] = Radius
// PotPar[0] = a0 for 1S0
// PotPar[1] = Reff for 1S0
// PotPar[2] = a0 for 3S1
// PotPar[3] = Reff for 3S1
double Lednicky_SingletTriplet(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    return GeneralLednicky(Momentum, SourcePar[0], PotPar[0], PotPar[1], PotPar[2], PotPar[3], false, false);
}

// this code is copy-pasted from Oliver's analysis for his thesis
// the only parameter I use is par[0] = Source Size
// The input should be in MeV and fm
double Lednicky_gauss_Sigma0(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    // This model tries to calculate the Sigma0 correlation function, one can simply say trying because nothing is known about
    // this kind of interaction

    // Parameter definitions:
    //*****************************************
    //*****************************************
    //              isospin 1/2 (==0)
    //________________________________________
    const complex<double> scatt_length_1s0_0(-1.1, 0.);
    const complex<double> scatt_length_3s1_0(-1.1, 4.3);
    const complex<double> effrange_1s0_0(-1.5, 0.);
    const complex<double> effrange_3s1_0(-2.2, -2.4);

    const double hbarc = 0.197326971844;
    // Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
    const double MomentumGeV = Momentum / 1000.;

    //
    //             isospin 3/2 (==1)
    //________________________________________
    const complex<double> scatt_length_1s0_1(2.51, 0.);
    const complex<double> scatt_length_3s1_1(-0.73, 0.);
    const complex<double> effrange_1s0_1(4.92, 0.);
    const complex<double> effrange_3s1_1(-1.22, 0.);

    const double mass_proton = 0.938272;
    const double mass_neutron = 0.939565;
    const double mass_sigmaplus = 1.189377;
    const double mass_sigma0 = 1.192642;

    double mu1 = mass_proton * mass_sigma0 / (mass_proton + mass_sigma0);
    double mu2 = mass_neutron * mass_sigmaplus / (mass_neutron + mass_sigmaplus);
    //*****************************************
    //*****************************************

    double k1 = MomentumGeV / hbarc; // momentum in the elastic channel
    double k2 = sqrt(2 * mu2 * (MomentumGeV * MomentumGeV / (2 * mu1) + mass_proton + mass_sigma0 - mass_neutron - mass_sigmaplus));
    k2 /= hbarc; // momentum in the inelastic channel
    //(k2>k1)

    double RG = SourcePar[0]; // Gaussian radius

    // Define inverse K-matrices for different spin configurations in isopin basis:

    // spin singlet:
    complex<double> Kmin_1s0_0 = 1. / scatt_length_1s0_0 + 0.5 * effrange_1s0_0 * k1 * k1;
    complex<double> Kmin_1s0_1 = 1. / scatt_length_1s0_1 + 0.5 * effrange_1s0_1 * k1 * k1;
    // spin triplet:
    complex<double> Kmin_3s1_0 = 1. / scatt_length_3s1_0 + 0.5 * effrange_3s1_0 * k1 * k1;
    complex<double> Kmin_3s1_1 = 1. / scatt_length_3s1_1 + 0.5 * effrange_3s1_1 * k1 * k1;

    // transform to transition basis with help of Clebsch-Gordan coefficients:
    complex<double> Kmin_11_1s0 = 2. / 3. * Kmin_1s0_1 + 1. / 3. * Kmin_1s0_0;
    complex<double> Kmin_11_3s1 = 2. / 3. * Kmin_3s1_1 + 1. / 3. * Kmin_3s1_0;

    complex<double> Kmin_22_1s0 = 1. / 3. * Kmin_1s0_1 + 2. / 3. * Kmin_1s0_0;
    complex<double> Kmin_22_3s1 = 1. / 3. * Kmin_3s1_1 + 2. / 3. * Kmin_3s1_0;

    complex<double> Kmin_12_1s0 = sqrt(2.) / 3. * (Kmin_1s0_1 - Kmin_1s0_0);
    complex<double> Kmin_12_3s1 = sqrt(2.) / 3. * (Kmin_3s1_1 - Kmin_3s1_0);

    // Determinant of scattering amplitude matrix:
    complex<double> D_1s0 = (Kmin_11_1s0 + complex<double>(0., -k1)) * (Kmin_22_1s0 + complex<double>(0., -k2)) - Kmin_12_1s0 * Kmin_12_1s0;
    complex<double> D_3s1 = (Kmin_11_3s1 + complex<double>(0., -k1)) * (Kmin_22_3s1 + complex<double>(0., -k2)) - Kmin_12_3s1 * Kmin_12_3s1;

    // What is really needed are the scattering amplitudes:

    complex<double> f11_1s0 = (Kmin_22_1s0 + complex<double>(0., -k2)) / D_1s0;
    complex<double> f11_3s1 = (Kmin_22_3s1 + complex<double>(0., -k2)) / D_3s1;

    complex<double> f12_1s0 = -(Kmin_12_1s0) / D_1s0;
    complex<double> f12_3s1 = -(Kmin_12_3s1) / D_3s1;

    double fF1 = gsl_sf_dawson(2 * k1 * RG) / (2 * k1 * RG);
    double fF2 = (1. - exp(-pow(2 * k1 * RG, 2))) / (2 * k1 * RG);

    // elastic part of the correlation function(1 -> 1):
    //---------------------------------------------
    double corr_1s0 = 0.5 * pow(abs(f11_1s0), 2.) / (RG * RG) + 2. * real(f11_1s0) / (sqrt(Pi) * RG) * fF1 - imag(f11_1s0) / RG * fF2;
    corr_1s0 *= 0.25;

    double corr_3s1 = 0.5 * pow(abs(f11_3s1), 2.) / (RG * RG) + 2. * real(f11_3s1) / (sqrt(Pi) * RG) * fF1 - imag(f11_3s1) / RG * fF2;
    corr_3s1 *= 0.75;
    //---------------------------------------------

    // inelastic part of the correlation function(2 -> 1):
    //---------------------------------------------
    double corr_1s0_inel = 0.5 * mu2 / mu1 * fabs(pow(abs(f12_1s0), 2.) / (RG * RG));
    corr_1s0_inel *= 0.25;

    double corr_3s1_inel = 0.5 * mu2 / mu1 * fabs(pow(abs(f12_3s1), 2.) / (RG * RG));
    corr_3s1_inel *= 0.75;

    //---------------------------------------------

    // correction term to non spheric distortions on the short range scale
    // Define matrices d_ij = 2Re d(K_ij)/dk²

    complex<double> dKmin_11_1s0_dk = 2. / 3. * 0.5 * effrange_1s0_1 + 1. / 3. * 0.5 * effrange_1s0_0;
    complex<double> dKmin_11_3s1_dk = 2. / 3. * 0.5 * effrange_3s1_1 + 1. / 3. * 0.5 * effrange_3s1_0;

    complex<double> dKmin_22_1s0_dk = 1. / 3. * 0.5 * effrange_1s0_1 + 2. / 3. * 0.5 * effrange_1s0_0;
    complex<double> dKmin_22_3s1_dk = 1. / 3. * 0.5 * effrange_3s1_1 + 2. / 3. * 0.5 * effrange_3s1_0;

    complex<double> dKmin_12_1s0_dk = sqrt(2.) / 3. * 0.5 * (effrange_1s0_1 - effrange_1s0_0);
    complex<double> dKmin_12_3s1_dk = sqrt(2.) / 3. * 0.5 * (effrange_3s1_1 - effrange_3s1_0);

    complex<double> factor1_1s0 = f11_1s0 * conj(f12_1s0);
    complex<double> factor1_3s1 = f11_3s1 * conj(f12_3s1);

    double corr_1s0_distort = pow(abs(f11_1s0), 2.) * 2. * real(dKmin_11_1s0_dk) + pow(abs(f12_1s0), 2.) * 2. * real(dKmin_22_1s0_dk) + 2. * real(factor1_1s0) * 2. * real(dKmin_12_1s0_dk);
    corr_1s0_distort *= -0.25 / (4 * sqrt(Pi) * pow(RG, 3.));

    double corr_3s1_distort = pow(abs(f11_3s1), 2.) * 2. * real(dKmin_11_3s1_dk) + pow(abs(f12_3s1), 2.) * 2. * real(dKmin_22_3s1_dk) + 2. * real(factor1_3s1) * 2. * real(dKmin_12_3s1_dk);
    corr_3s1_distort *= -0.75 / (4 * sqrt(Pi) * pow(RG, 3.));

    double corr_fin = 1. + corr_1s0 + corr_3s1 + corr_1s0_inel + corr_3s1_inel + corr_1s0_distort + corr_3s1_distort;

    return corr_fin;
    // if(par[1] == 0.) return 1. + par[2]*(par[3]*corr_fin -1.); //this is with baseline
    // else return par[2]*(par[3]*corr_fin -1.);//this is baseline subtracted
}

// this code is copy-pasted from Oliver's analysis for his thesis
double pXi_pheno(const double &Momentum, const double *SourcePar, const double *PotPar)
{

    const double hbarc = 0.197326971844;
    const double MomentumGeV = Momentum / 1000.;

    // phenomenological function developed for fitting p-Xi
    const double &RG = SourcePar[0];

    // very easy definition for a CF that peaks for k->0
    return 1. + exp(-MomentumGeV * RG / hbarc) / (MomentumGeV * RG / hbarc);
}

// This implementation of Lednicky takes as input a parametrization of phase
// shift and elasticity coefficient. It uses a pol2 to model the real part of
// the phase shift, and a combination of exp and pol2 for the imaginary part
// PotPar[0] - [2] -> phase shift
// PotPar[3] - [5] -> elasticity coefficient
double LednickySingletScatAmplitude(const double &kStar,
                                    const double *SourcePar,
                                    const double *PotPar)
{
    // Parametrization of the elasticity coefficient
    auto eta = [](const double kStar, const double *PotPar)
    {
        return std::exp(PotPar[3] * kStar) + PotPar[4] * kStar +
               PotPar[5] * kStar * kStar;
    };

    // Parametrization of the phase shift
    auto delta = [](const double kStar, const double *PotPar)
    {
        return PotPar[0] + PotPar[1] * kStar + PotPar[2] * kStar * kStar;
    };

    // Compute the scattering amplitude
    auto ScatteringAmplitude = [&delta, &eta](const double &kStar,
                                              const double *PotPar)
    {
        double etaVal = eta(kStar, PotPar);
        double deltaVal = delta(kStar, PotPar);
        std::complex<double> exponent(0., 2. * deltaVal);
        std::complex<double> exp(std::exp(exponent));
        std::complex<double> nominator(etaVal * exp - 1.);
        std::complex<double> denominator(0, 2. * kStar);
        return nominator / denominator;
    };

    // Compute 1/K for the second derivative
    auto OneOverK = [&delta, &eta](const double &kStar, const double *PotPar)
    {
        const double etaVal = eta(kStar, PotPar);
        const double etaValSq = etaVal * etaVal;
        const double deltaVal = delta(kStar, PotPar);
        std::complex<double> nominator(2 * etaVal * std::sin(2. * deltaVal),
                                       -(1. - etaValSq));
        return nominator * kStar /
               (1 + etaVal * etaVal - 2. * etaValSq * std::cos(2. * deltaVal));
    };

    const double radius = SourcePar[0] * FmToNu;

    // Precision with which the second derivative is conducted
    const double precision = 1.;

    std::complex<double> ScatAmpl = ScatteringAmplitude(kStar, PotPar);

    const double oneOverKSecondDerivative = std::real(
        ((OneOverK(kStar + precision, PotPar) - OneOverK(kStar, PotPar)) -
         (OneOverK(kStar, PotPar) - OneOverK(kStar - precision, PotPar))) /
        (precision * precision));

    double F1 = gsl_sf_dawson(2. * kStar * radius) / (2. * kStar * radius);
    double F2 = (1. - std::exp(-4. * kStar * kStar * radius * radius)) /
                (2. * kStar * radius);

    double CkValue =
        0.5 * std::pow(std::abs(ScatAmpl) / radius, 2) *
            (1. -
             (oneOverKSecondDerivative) / (2 * std::sqrt(Pi) * radius)) +
        2 * std::real(ScatAmpl) * F1 / (std::sqrt(Pi) * radius) -
        std::imag(ScatAmpl) * F2 / radius;

    return CkValue + 1;
}

double Lednicky_CC_pAL(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    //     //This model tries to calculate the CF for 2 channels in pAL
    //     //Elastic scattering parameters taken from WUT measurments in HIC

    //     //Parameter definitions:
    complex<double> ScatLen_1(PotPar[0], PotPar[1]); // elastic
    complex<double> effrange_1(PotPar[2], 0.);
    double RG = SourcePar[0]; // Gaussian radius

    const double hbarc = 0.197326971844;
    //     //Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
    const double MomentumGeV = Momentum / 1000.;

    complex<double> ScatLen_2(PotPar[3], PotPar[4]); // inelastic
    complex<double> effrange_2(PotPar[5], 0.);

    const double mass_proton = 0.938272;
    const double mass_Lambda = 1.115683;
    const double mass_cc1 = 0.1372736; // assuming for now only 2 body πK
    const double mass_cc2 = 0.495644;  // avg mass for π+ and π0 and for K+ and K0

    double mu1 = mass_proton * mass_Lambda / (mass_proton + mass_Lambda);
    double mu2 = mass_cc1 * mass_cc2 / (mass_cc1 + mass_cc2);
    //     //*****************************************
    //     //*****************************************

    double k1 = MomentumGeV / hbarc; // momentum in the elastic channel
    double k2 = sqrt(2 * mu2 * (MomentumGeV * MomentumGeV / (2 * mu1) + mass_proton + mass_Lambda - mass_cc1 - mass_cc2));
    k2 /= hbarc; // momentum in the inelastic channel
    //(k2>k1)

    // Define inverse K-matrices
    complex<double> Kinv_11 = 1. / ScatLen_1 + 0.5 * effrange_1 * k1 * k1;
    complex<double> Kinv_22 = Kinv_11; // we assume the same type of interaction for channel
    complex<double> Kinv_21 = 1. / ScatLen_2 + 0.5 * effrange_2 * k1 * k2;

    // Determinant of scattering amplitude matrix:
    complex<double> Det = (Kinv_11 + complex<double>(0., -k1)) * (Kinv_22 + complex<double>(0., -k2)) - Kinv_21 * Kinv_21;

    // What is really needed are the scattering amplitudes:
    complex<double> f11 = (Kinv_22 + complex<double>(0., -k2)) / Det;
    complex<double> f22 = (Kinv_11 + complex<double>(0., -k1)) / Det;
    complex<double> f21 = -(Kinv_21) / Det;

    // printf("---debug 1 ---\n");
    double fF1 = gsl_sf_dawson(2 * k1 * RG) / (2 * k1 * RG);
    // printf("fF1 = %.6f\n",fF1);
    double fF2 = (1. - exp(-pow(2 * k1 * RG, 2))) / (2 * k1 * RG);

    // elastic part of the correlation function(1 -> 1):
    //---------------------------------------------
    double corr_el = 0.5 * pow(abs(f11), 2.) / (RG * RG) + 2. * real(f11) / (sqrt(Pi) * RG) * fF1 - imag(f11) / RG * fF2;
    //---------------------------------------------

    // inelastic part of the correlation function(2 -> 1):
    //---------------------------------------------
    double corr_inel = 0.5 * mu2 / mu1 * fabs(pow(abs(f21), 2.) / (RG * RG));

    //---------------------------------------------

    // correction term to non spheric distortions on the short range scale
    // Define matrices d_ij = 2Re d(K_ij)/dk²

    complex<double> d11_dk = 0.5 * effrange_1;
    complex<double> d22_dk = 0.5 * effrange_1;
    complex<double> d21_dk = 0.5 * effrange_2;

    complex<double> factor1 = f11 * conj(f21);

    double corr_distort = pow(abs(f11), 2.) * 2. * real(d11_dk) + pow(abs(f21), 2.) * 2. * real(d22_dk) + 2. * real(factor1) * 2. * real(d21_dk);
    corr_distort *= -1. / (4 * sqrt(Pi) * pow(RG, 3.));

    double corr_fin = 1. + corr_el + corr_inel + corr_distort;

    return corr_fin;
    //     //if(par[1] == 0.) return 1. + par[2]*(par[3]*corr_fin -1.); //this is with baseline
    //     //else return par[2]*(par[3]*corr_fin -1.);//this is baseline subtracted
}

double Lednicky_CC_LAL(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    //     //This model tries to calculate the CF for 2 channels in pAL
    //     //Elastic scattering parameters taken from WUT measurments in HIC

    //     //Parameter definitions:
    complex<double> ScatLen_1(PotPar[0], PotPar[1]); // elastic
    complex<double> effrange_1(PotPar[2], 0.);
    double RG = SourcePar[0]; // Gaussian radius

    const double hbarc = 0.197326971844;
    //     //Oli's original code was in GeV and fm though, so we must convert the Momentum to GeV!
    const double MomentumGeV = Momentum / 1000.;

    complex<double> ScatLen_2(PotPar[3], PotPar[4]); // inelastic
    complex<double> effrange_2(PotPar[5], 0.);

    const double mass_Lambda = 1.115683;
    const double mass_cc1 = 0.495644; // assuming for now only 2 body KantiK
    const double mass_cc2 = 0.495644; // avg mass for K+ and K0

    double mu1 = mass_Lambda * mass_Lambda / (mass_Lambda + mass_Lambda);
    double mu2 = mass_cc1 * mass_cc2 / (mass_cc1 + mass_cc2);
    //     //*****************************************
    //     //*****************************************

    double k1 = MomentumGeV / hbarc; // momentum in the elastic channel
    double k2 = sqrt(2 * mu2 * (MomentumGeV * MomentumGeV / (2 * mu1) + mass_Lambda + mass_Lambda - mass_cc1 - mass_cc2));
    k2 /= hbarc; // momentum in the inelastic channel
    //(k2>k1)

    // Define inverse K-matrices
    complex<double> Kinv_11 = 1. / ScatLen_1 + 0.5 * effrange_1 * k1 * k1;
    complex<double> Kinv_22 = Kinv_11; // we assume the same type of interaction for channel
    complex<double> Kinv_21 = 1. / ScatLen_2 + 0.5 * effrange_2 * k1 * k2;

    // Determinant of scattering amplitude matrix:
    complex<double> Det = (Kinv_11 + complex<double>(0., -k1)) * (Kinv_22 + complex<double>(0., -k2)) - Kinv_21 * Kinv_21;

    // What is really needed are the scattering amplitudes:
    complex<double> f11 = (Kinv_22 + complex<double>(0., -k2)) / Det;
    complex<double> f22 = (Kinv_11 + complex<double>(0., -k1)) / Det;
    complex<double> f21 = -(Kinv_21) / Det;

    printf("---debug 1 ---\n");
    double fF1 = gsl_sf_dawson(2. * k1 * RG) / (2. * k1 * RG);
    printf("fF1 = %.2f\n", fF1);
    double fF2 = (1. - exp(-pow(2 * k1 * RG, 2))) / (2 * k1 * RG);

    // double F1 = gsl_sf_dawson(2. * kStar * radius) / (2. * kStar * radius);
    // double F2 = (1. - std::exp(-4. * kStar * kStar * radius * radius)) /
    //             (2. * kStar * radius);
    // elastic part of the correlation function(1 -> 1):
    //---------------------------------------------
    double corr_el = 0.5 * pow(abs(f11), 2.) / (RG * RG) + 2. * real(f11) / (sqrt(Pi) * RG) * fF1 - imag(f11) / RG * fF2;
    //---------------------------------------------

    // inelastic part of the correlation function(2 -> 1):
    //---------------------------------------------
    double corr_inel = 0.5 * mu2 / mu1 * fabs(pow(abs(f21), 2.) / (RG * RG));

    //---------------------------------------------

    // correction term to non spheric distortions on the short range scale
    // Define matrices d_ij = 2Re d(K_ij)/dk²

    complex<double> d11_dk = 0.5 * effrange_1;
    complex<double> d22_dk = 0.5 * effrange_1;
    complex<double> d21_dk = 0.5 * effrange_2;

    complex<double> factor1 = f11 * conj(f21);

    double corr_distort = pow(abs(f11), 2.) * 2. * real(d11_dk) + pow(abs(f21), 2.) * 2. * real(d22_dk) + 2. * real(factor1) * 2. * real(d21_dk);
    corr_distort *= -1. / (4 * sqrt(Pi) * pow(RG, 3.));

    double corr_fin = 1. + corr_el + corr_inel + corr_distort;

    return corr_fin;
    //     //if(par[1] == 0.) return 1. + par[2]*(par[3]*corr_fin -1.); //this is with baseline
    //     //else return par[2]*(par[3]*corr_fin -1.);//this is baseline subtracted
}

/// @brief  Lednicky-Lyuboshits formula with Flatte-like scattering amplitude for the decay into the pair we are measuring (assumed channel 2) based on the Sill distribution (F. Giacosa, Eur. Phys. J. A (2021) 57:336)
/// @param Momentum: momentum of the pair you are measuring
/// @param GaussR: radius of the source
/// @param MassR: mass of the resonance (provide in MeV)
/// @param Gamma1: width of the decay channel 1, lowest thresholds (e.g. πΞ for Ξ(1620))
/// @param Gamma2: width of the decay channel 2, highes thresholds (e.g. ΛK- for Ξ(1620)), typically the pair you are measuring
/// @param En1: energy threshold channel 1, sum of masses of pair 1 (e.g. mπ+mΞ)
/// @param En2: energy threshold channel 2, sum of masses of pair 2 (e.g. mΛ+mK)

double GeneralLednickySill_twochannels(const double &Momentum, const double &GaussR, const double &MassR, const double &Gamma1, const double &Gamma2, const double &m11, const double &m12, const double &m21, const double &m22)
{

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednickySill got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);

    double En1 = m11 + m12;
    double En2 = m21 + m22;
    double RedMass = (m12 * m22) / (m12 + m22);

    /// You are computing the scatt. amplitude for the decay channel = pair you are measuring
    /// the coupling/Gamma above should be always channel 2, the one with the highest threshold
    double En = (Momentum * Momentum) / (2. * RedMass) + En2;

    complex<double> ScattAmpl = (-2 * Gamma2) * (pow(En * En - MassR * MassR + i * Gamma1 * sqrt(En * En - En1 * En1) + i * Gamma2 * sqrt(En * En - En2 * En2), -1.));

    /// correction for small sources, involving
    double DeltaC = pow(abs(ScattAmpl), 2) * (2. + m22 / m12 + m12 / m22) / (2. * sqrt(Pi) * Radius * Radius * Radius * 2 * Gamma2);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmpl) / Radius, 2) +
               2 * real(ScattAmpl) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmpl) * F2 / Radius + DeltaC;

    CkValue += 1;

    return CkValue;
}

double LednickySilltwochannels_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    double MassR = PotPar[0];
    double Gamma1 = PotPar[1];
    double Gamma2 = PotPar[2];
    double m11 = PotPar[3];
    double m12 = PotPar[4];
    double m21 = PotPar[5];
    double m22 = PotPar[6];
    return (SourcePar[3] * (SourcePar[2] * GeneralLednickySill_twochannels(Momentum, SourcePar[0], MassR, Gamma1, Gamma2, m11, m12, m21, m22) + (1 - SourcePar[2]) * GeneralLednickySill_twochannels(Momentum, SourcePar[1], MassR, Gamma1, Gamma2, m11, m12, m21, m22)) + 1. - SourcePar[3]);
}

// ###########################################################################################################################
/// Same as  but here we pass as a second width parameter the tot width, as e.g. measured by Belle
double GeneralLednickySillConstraints_twochannels(const double &Momentum, const double &GaussR, const double &MassR, const double &GammaTilde1, const double &GammaTot, const double &m11, const double &m12, const double &m21, const double &m22)
{

    if (GaussR != GaussR)
    {
        printf("\033[1;33mWARNING:\033[0m GeneralLednickySill got a bad value for the Radius (nan). Returning default value of 1.\n");
        return 1;
    }

    const double Radius = GaussR * FmToNu;

    double F1 = gsl_sf_dawson(2. * Momentum * Radius) / (2. * Momentum * Radius);
    double F2 = (1. - exp(-4. * Momentum * Momentum * Radius * Radius)) / (2. * Momentum * Radius);

    double En1 = m11 + m12;
    double En2 = m21 + m22;
    double RedMass = (m12 * m22) / (m12 + m22);
    double alpha1 = En1 / MassR;
    double alpha2 = En2 / MassR;

    /// You are computing the scatt. amplitude for the decay channel = pair you are measuring
    /// the coupling/Gamma above should be always channel 2, the one with the highest threshold
    double En = (Momentum * Momentum) / (2. * RedMass) + En2;

    double GammaTilde2 = GammaTot / (sqrt(1. - alpha2 * alpha2)) - GammaTilde1 * sqrt(1. - alpha1 * alpha1) / (sqrt(1. - alpha2 * alpha2));

    complex<double> ScattAmpl = (-2 * GammaTilde2) * (pow(En * En - MassR * MassR + i * GammaTilde1 * sqrt(En * En - En1 * En1) + i * GammaTilde2 * sqrt(En * En - En2 * En2), -1.));

    /// correction for small sources, involving
    double DeltaC = pow(abs(ScattAmpl), 2) * (2. + m22 / m12 + m12 / m22) / (2. * sqrt(Pi) * Radius * Radius * Radius * 2 * GammaTilde2);

    double CkValue = 0.;
    CkValue += 0.5 * pow(abs(ScattAmpl) / Radius, 2) +
               2 * real(ScattAmpl) * F1 / (sqrt(Pi) * Radius) - imag(ScattAmpl) * F2 / Radius + DeltaC;

    CkValue += 1;

    return CkValue;
}

double LednickySilltwochannelsConstraints_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar)
{
    double MassR = PotPar[0];
    double GammaTilde1 = PotPar[1];
    double GammaTot = PotPar[2];
    double m11 = PotPar[3];
    double m12 = PotPar[4];
    double m21 = PotPar[5];
    double m22 = PotPar[6];
    return (SourcePar[3] * (SourcePar[2] * GeneralLednickySillConstraints_twochannels(Momentum, SourcePar[0], MassR, GammaTilde1, GammaTot, m11, m12, m21, m22) + (1 - SourcePar[2]) * GeneralLednickySillConstraints_twochannels(Momentum, SourcePar[1], MassR, GammaTilde1, GammaTot, m11, m12, m21, m22)) + 1. - SourcePar[3]);
}
// ###############################################
/// @brief  Lednicky-Lyuboshits formula with Flatte-like scattering amplitude for the decay into the pair we are measuring (assumed channel 2) based on the Sill distribution (F. Giacosa, Eur. Phys. J. A (2021) 57:336), no constraints + LL formula with ere
double LednickySilltwochannelsERE_doublegaussian_lambda(const double &Momentum, const double *SourcePar, const double *PotPar)
{

    /// Ref0 [0]
    /// Imf0 [1]
    /// d0   [2]
    complex<double> ScatLen(PotPar[0], PotPar[1]);
    double cfemto = SourcePar[3] * (SourcePar[2] * GeneralLednicky(Momentum, SourcePar[0], ScatLen, PotPar[2], 0, 0, true, false, false) + (1 - SourcePar[2]) * GeneralLednicky(Momentum, SourcePar[1], ScatLen, PotPar[2], 0, 0, true, false, false)) + 1. - SourcePar[3];

    double MassR = PotPar[3];
    double GammaTilde1 = PotPar[4];
    double GammaTilde2 = PotPar[5];
    double m11 = PotPar[6];
    double m12 = PotPar[7];
    double m21 = PotPar[8];
    double m22 = PotPar[9];

    double csill = (SourcePar[3] * (SourcePar[2] * GeneralLednickySill_twochannels(Momentum, SourcePar[0], MassR, GammaTilde1, GammaTilde2, m11, m12, m21, m22) + (1 - SourcePar[2]) * GeneralLednickySill_twochannels(Momentum, SourcePar[1], MassR, GammaTilde1, GammaTilde2, m11, m12, m21, m22)) + 1. - SourcePar[3]);

    double weight = PotPar[10];

    return (weight * cfemto + (1 - weight) * csill);
}