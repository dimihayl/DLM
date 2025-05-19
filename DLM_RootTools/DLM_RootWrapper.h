#ifndef DLM_ROOTWRAPPER_H
#define DLM_ROOTWRAPPER_H

template <class Type> class DLM_Histo;

class TH2F;
class TH1F;
class TGraph;
class CATS;

DLM_Histo<float>* Convert_TH2F_DlmHisto(const TH2F* input);
DLM_Histo<float>* Convert_TH1F_DlmHisto(const TH1F* input);
DLM_Histo<double>* Convert_TH1F_DoubleDlmHisto(const TH1F* input);
TH2F* Convert_DlmHisto_TH2F(const DLM_Histo<float>* input, const char* name);
TH1F* Convert_DlmHisto_TH1F(const DLM_Histo<float>* input, const char* name);
TH1F* GetSourceFromCATS_TH1F(CATS& kitty, const double kstar, const double costheta, const char* name);
TGraph* GetSourceFromCATS_TGraph(CATS& kitty, const double kstar, const double costheta, const char* name);

#endif // DLM_ROOTWRAPPER_H
