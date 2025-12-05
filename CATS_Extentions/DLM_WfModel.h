#ifndef DLM_WFMODEL_H
#define DLM_WFMODEL_H

#include <complex>

class CATS;
template <class Type> class DLM_Histo;

using namespace std;

DLM_Histo<complex<double>>*** Init_pp_Haidenbauer(const char* CatsFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** Init_pp_Haidenbauer(const char* CatsFolder, CATS* Kitty, const int& TYPE);

void Init_pp_Epelbaum(const char* CatsFolder, CATS& Kitty, const int& TYPE);
void Init_pp_Epelbaum(const char* CatsFolder, CATS* Kitty, const int& TYPE);

//TYPE = XYZ, where:
//NV2-XY_EMZ:
//X = 1 or 2 = fit of scattering pars up to 125 or 200 MeV in TLab
//Y = 1 or 2 (corresponding to A or B) has different cut off in coordinate space
//Z = 1 or 2 has different treatment of Coulomb, 2 = takes higher order corrections as well
//the default settings, in case we get no TYPE or 0, are 212
//e.g. check here https://arxiv.org/pdf/1412.6446
//Implementation provided by Matthias Goebel, June 2025
void Init_pp_Norfolk(const char* CatsFolder, CATS& Kitty, const int& TYPE);
void Init_pp_Norfolk(const char* CatsFolder, CATS* Kitty, const int& TYPE);

//Implementation provided by Matthias Goebel, June 2025
void Init_pp_AV18_WF(const char* CatsFolder, CATS& Kitty, const int& TYPE);
void Init_pp_AV18_WF(const char* CatsFolder, CATS* Kitty, const int& TYPE);

DLM_Histo<complex<double>>*** Init_pp_Epelbaum_OLD(const char* CatsFolder, CATS& Kitty, const int& TYPE);
DLM_Histo<complex<double>>*** Init_pp_Epelbaum_OLD(const char* CatsFolder, CATS* Kitty, const int& TYPE);


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


DLM_Histo<complex<double>> ***Init_pXi_Haidenbauer2025(const char *InputFolder, CATS &Kitty, const int &TYPE = 0, const int &CUTOFF = 0);
DLM_Histo<complex<double>> ***Init_pXi_Haidenbauer2025(const char *InputFolder, CATS *Kitty, const int &TYPE = 0, const int &CUTOFF = 0);


/// TYPE selects the different models implemented in Eur.Phys.J.A 56 (2020) 7, 184
/// TYPE = 0 -> LQCDe (CUTOFF=500,CUTOFF=600)
/// TYPE = 1 -> Model A (CV) most attractive
/// TYPE = 2 -> CTNN-d (OKA) bound state
DLM_Histo<complex<double>> ***Init_Lcp_Haidenbauer(const char *InputFolder, CATS &Kitty, const int &TYPE, const int &CUTOFF = 500);
DLM_Histo<complex<double>> ***Init_Lcp_Haidenbauer(const char *InputFolder, CATS *Kitty, const int &TYPE, const int &CUTOFF = 500);

void CleanUpWfHisto(const unsigned short& NumChannels, DLM_Histo<complex<double>>***& Histo);
void CleanUpWfHisto(const CATS& Kitty, DLM_Histo<complex<double>>***& Histo);
#endif
