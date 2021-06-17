#include "TFitResultPtr.h"
#include "Math/MinimizerOptions.h"
//#include "RtypesCore.h"

class TH1;
class TF1;

TFitResultPtr DLM_FitHisto( TH1* h1, TF1* f1, Option_t* option="", Option_t* goption="",
                        const ROOT::Math::MinimizerOptions* moption=NULL,
                        Double_t xxmin=0, Double_t xxmax=0);
