#ifndef DLM_WFMODEL_H
#define DLM_WFMODEL_H

#include <complex>

class CATS;

using namespace std;

DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
void CleanUpWfHisto(const unsigned short& NumChannels, DLM_Histo<complex<double>>***& Histo);
void CleanUpWfHisto(const CATS& Kitty, DLM_Histo<complex<double>>***& Histo);
#endif
