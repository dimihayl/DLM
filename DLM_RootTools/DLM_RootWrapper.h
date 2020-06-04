#ifndef DLM_ROOTWRAPPER_H
#define DLM_ROOTWRAPPER_H

#include "DLM_Histo.h"

class TH2F;

DLM_Histo<float>* Convert_TH2F_DlmHisto(const TH2F* input);
TH2F* Convert_DlmHisto_TH2F(const DLM_Histo<float>* input, const char* name);


#endif // DLM_ROOTWRAPPER_H

