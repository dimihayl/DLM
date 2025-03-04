#ifndef DLM_Potentials_H
#define DLM_Potentials_H

#include <stdio.h>
#include <complex>
using namespace std;

//which potential to use
//extern int DlmPot;
//pp_ReidSC is the castrated Reid
//pp_ReidOli models S=1 with 3P1
//pp_ReidCrab is the castrated Reid including ONLY 3P2 for the S=1 state
enum DLMPOT {   NN_AV18, NN_ReidV8, pp_ReidSC, pp_ReidOli, pp_ReidCrab, pp_ReidVale,
                pL_UsmaniOli, pXim_Lattice, pXim_HALQCD1, pXim_HALQCDPaper2020, pXim_LatticeAvg, pXim_LatticeSqrtAvg, pKm_Tetsuo, pOmega_Lattice, pOmega_Tetsuo};

//specific flags that are to be passed the potentials. 0 is the default
//extern int DlmPotFlag;
enum V18FLAG { v18_Default, v18_SingleChannelMagic, v18_Coupled3P2 };

void CleanUpV18Pot();

double ZeroPotential(double* Radius);

double SingleGauss(double* Pars);
double SingleGaussDynamic(double* Pars);
double DoubleGaussSum(double* Pars);
//V0*exp(-r^2/β0^2)+V1*exp(-r^2/β1^2)+V2*exp(-r^2/β2^2)
//[0] - r; [1] = k; [2]=V0; [3]=μ0; [4]=V1; [5]=μ1; [6]=V2; [7]=μ2
double TripleGaussSum(double* Pars);
double GaussExpSum(double* Pars);
double UsmaniPotentialCats(double* Pars);
double UsmaniFit(double* Pars);
double UsmaniFitAll(double* Pars);
double UsmaniPotentialLonardoni(double* Pars);
double RepulsiveCore(double* Pars);


double Gaussian(double* Pars);
double Yukawa(double* Pars);
double YukawaRepCore(double* Pars);
double YukawaGaussCore(double* Pars);
double YukawaDimiCore(double* Pars);
double YukawaDimiSmooth(double* Pars);

double Hulthen(double* Pars);
double HulthenSmooth(double* Pars);
double SquareWell(double* Pars);

void SetUpNorfolk(const char* InputFolder);
double pp_Norfolk(double* Pars);

//the AV18 with tons of options to play with
//a repulsive core here is just a Gaussian centered at 0
//P0 [2] = which channel. 0 = 1S0, 1 = 3P0, 2 = 3P1, 3 = 3P2
//P1/2 [3/4] = amplitude/width of a PW
double pp_AV18_Toy(double* Pars);

double LatticePots_pXi(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars);
double LatticePots_pXi_Avg(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars);
double LatticePots_pOmega(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius, double* OtherPars);
struct LatticeValues;
struct LatticeValuesPaper;

double Tetsuo_pKm(const int& WhichPot, const int& DlmPotFlag,
                     const int& IsoSpin, const int& t2p1, const int& t2p2,
                     const int& Spin, const int& AngMom, const int& TotMom, double* Radius);

double KpProtonEquivalentPotential(double* pars);

double ReidSoftCore1S0(double* Radius);
double ReidSoftCore3P(double* Radius);

double fReidMeVfm1S0(double* Radius);
double fReidMeVfm3P(double* Radius);

double fReidDlm1S0(double* Radius);
double fReidDlm3P0(double* Radius);
double fReidDlm3P1(double* Radius);
double fReidDlm3P2(double* Radius);

double fReidVale1S0(double* Radius);
double fReidVale3P0(double* Radius);
double fReidVale3P1(double* Radius);
double fReidVale3P2(double* Radius);

//double fV18potential1S0(double* Radius);
//double fV18potential3P0(double* Radius);
//double fV18potential3P1(double* Radius);
//double fV18potential3P2(double* Radius);

double ppDlmPot(const int& DlmPot, const int& DlmFlag, const int& Spin, const int& AngMom, const int& TotMom, double* Radius);
double ppDlmPot1S0(double* Radius);
double ppDlmPot3S1(double* Radius);
double ppDlmPot3P0(double* Radius);
double ppDlmPot3P1(double* Radius);
double ppDlmPot3P2(double* Radius);
double ppDlmPot3P(double* Radius);

double pLambdaDlmPot1S0(double* Pars);
double pLambdaDlmPot3S1(double* Pars);

double fDlmPot(const int& DlmPot, const int& DlmPotFlag, const int& IsoSpin, const int& t2p1, const int& t2p2, const int& Spin, const int& AngMom, const int& TotMom,
               double* Radius, const double& CutOff, double* OtherPars=NULL);
double fDlmPot(double* Parameters);
double fDlmPotVer2(double* Parameters);
double LatticePots(double* Parameters);

void GetDlmPotName(const int& potid, const int& potflag, char* name);

/// Complex potentials to test
/// 1) used to test iCATS_0: ΛantiK complex potential by Yuki
///     V(r) = V0 * exp(-m^2 r^2), with V0 complex and m = m_ρ
///     V0 = V0R + i V0i = -1547-i*589 MeV
///     [0] -> r; [1] = k; [2]=V0real; [3]=V0imag; [4] = m (m_ρ)
complex<double> ComplexGaussian(double *Pars);

complex<double> ComplexPotLAntiK_mesonexchange(double *Pars);
double RealPotLAntiK_mesonexchange(double *Pars);

#endif
