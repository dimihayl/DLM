#ifndef DLM_WFMODEL_H
#define DLM_WFMODEL_H

#include "DLM_WfModel.h"
#include <complex>

class CATS;

using namespace std;

DLM_Histo<complex<double>>** Init_pL_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF=600);

#endif
