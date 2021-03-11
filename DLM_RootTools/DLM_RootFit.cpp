
#include "DLM_RootFit.h"

#include "TH1.h"
#include "Fit/Fitter.h"
#include "HFitInterface.h"

TFitResultPtr FitHisto( TH1* h1, TF1* f1, Option_t* option, Option_t* goption,
                        const ROOT::Math::MinimizerOptions& moption,
                        Double_t xxmin, Double_t xxmax){

  // implementation of Fit method is in file hist/src/HFitImpl.cxx
  Foption_t fitOption;
  ROOT::Fit::FitOptionsMake(ROOT::Fit::kHistogram,option,fitOption);

  // create range and minimizer options with default values
  ROOT::Fit::DataRange range(xxmin,xxmax);
  ROOT::Math::MinimizerOptions minOption;

  // need to empty the buffer before
  // (t.b.d. do a ML unbinned fit with buffer data)
  h1->BufferEmpty();

  return ROOT::Fit::FitObject(h1, f1 , fitOption , minOption, goption, range);
}
