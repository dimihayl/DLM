#ifndef DLM_WFMODEL_H
#define DLM_WFMODEL_H

#include <complex>

class CATS;
template <class Type> class DLM_Histo;

using namespace std;

DLM_Histo<complex<double>>*** Init_pp_Haidenbauer(const char* CatsFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** Init_pp_Haidenbauer(const char* CatsFolder, CATS* Kitty, const int& TYPE);

DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2019(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2019(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF=600);

DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2023(const char* InputFolder, CATS& Kitty, const char* Singlet, const char* Triplet);
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2023(const char* InputFolder, CATS* Kitty, const char* Singlet, const char* Triplet);

DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pSigmaPlus_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pSigmaPlus_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0, const int& CUTOFF=600);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS& Kitty, const int& TYPE=1);
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS* Kitty, const int& TYPE=1);
DLM_Histo<complex<double>>*** Init_pS0_ESC08(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pS0_ESC08(const char* InputFolder, CATS* Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pS0_ESC16(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pS0_ESC16(const char* InputFolder, CATS* Kitty, const int& TYPE=0);

DLM_Histo<complex<double>>*** Init_pantip_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** Init_LAntiK_Vidana(const char *InputFolder, CATS &Kitty, const int &TYPE);

DLM_Histo<complex<double>>*** Init_pKminus_Kyoto2019(const char* InputFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** Init_pK0_Kyoto2019(const char* InputFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** InitHaidenbauerKaonPlus(const char* InputFolder, CATS& Kitty, const int& TYPE);

DLM_Histo<complex<double>>*** Init_pd_Sebastian(const char* InputFolder, CATS& Kitty, const int& TYPE=0, const int& CUTOFF=400);

DLM_Histo<complex<double>>*** Init_pDminus_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE=0);
DLM_Histo<complex<double>>*** Init_pDminus_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE=0);

/// TYPE selects the different models implemented in Eur.Phys.J.A 56 (2020) 7, 184
/// TYPE = 0 -> LQCDe (CUTOFF=500,CUTOFF=600)
/// TYPE = 1 -> Model A (CV) most attractive
/// TYPE = 2 -> CTNN-d (OKA) bound state
DLM_Histo<complex<double>> ***Init_Lcp_Haidenbauer(const char *InputFolder, CATS &Kitty, const int &TYPE, const int &CUTOFF = 500);
DLM_Histo<complex<double>> ***Init_Lcp_Haidenbauer(const char *InputFolder, CATS *Kitty, const int &TYPE, const int &CUTOFF = 500);

void CleanUpWfHisto(const unsigned short& NumChannels, DLM_Histo<complex<double>>***& Histo);
void CleanUpWfHisto(const CATS& Kitty, DLM_Histo<complex<double>>***& Histo);
#endif
