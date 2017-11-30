#ifndef DLM_Potentials_H
#define DLM_Potentials_H

//which potential to use
//extern int DlmPot;
enum DLMPOT {   NN_AV18, NN_ReidV8, pp_ReidSC, pp_ReidOli, pp_ReidCrab,
                pL_UsmaniOli};

//specific flags that are to be passed the potentials. 0 is the default
//extern int DlmPotFlag;
enum V18FLAG { v18_Default, v18_SingleChannelMagic, v18_Coupled3P2 };

void CleanUpV18Pot();

double ZeroPotential(double* Radius);

double ReidSoftCore1S0(double* Radius);
double ReidSoftCore3P(double* Radius);

double fReidMeVfm1S0(double* Radius);
double fReidMeVfm3P(double* Radius);

double fReidDlm1S0(double* Radius);
double fReidDlm3P0(double* Radius);
double fReidDlm3P1(double* Radius);
double fReidDlm3P2(double* Radius);

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

double fDlmPot(const int& DlmPot, const int& DlmPotFlag, const int& IsoSpin, const int& t2p1, const int& t2p2, const int& Spin, const int& AngMom, const int& TotMom, double* Radius);
double fDlmPot(double* Parameters);

void GetDlmPotName(const int& potid, const int& potflag, char* name);

#endif
