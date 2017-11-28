#ifndef DLM_Potentials_H
#define DLM_Potentials_H

enum V18FLAG { v18_Default, v18_SingleChannelMagic, v18_Coupled3P2 };
//which potential to use
extern int DlmPot;
enum DLMPOT {   pp_AV18, pp_ReidV8, pp_ReidSC, pp_ReidOli, pp_ReidCrab,
                pL_UsmaniOli};

//specific flags that are to be passed the potentials. 0 is the default
extern int DlmPotFlag;

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

double fDlmPot(const int& Spin, const int& AngMom, const int& TotMom, double* Radius);
double fDlmPot(double* Parameters);
double fDlmPot1S0(double* Radius);
double fDlmPot3S1(double* Radius);
double fDlmPot3P0(double* Radius);
double fDlmPot3P1(double* Radius);
double fDlmPot3P2(double* Radius);
double fDlmPot3P(double* Radius);

void GetDlmPotName(const int& potid, const int& potflag, char* name);

#endif
