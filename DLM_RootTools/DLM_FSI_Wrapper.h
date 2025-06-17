#ifndef DLM_FSI_WRAPPER_H
#define DLM_FSI_WRAPPER_H

#include <complex>

//template <class Type> class DLM_Histo;

class TH2F;
class TH1F;

void Wrap_pp_Epelbaum(const char* InputFileName, const char* OutputFolder);
//used for both Init_pp_Norfolk and Init_pp_AV18_WF
void Wrap_pp_Norfolk(const char* InputFileName, const char* OutputFolder, const char* descriptor);


#endif // DLM_FSI_WRAPPER_H
