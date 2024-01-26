
#include "DLM_RootFit.h"

#include "TH1.h"
#include "Fit/Fitter.h"
#include "HFitInterface.h"

#include "TSpline.h"

double DlmTSplineFit(double* xVal, double* pars){
    //[0] = NumKnots
    //[1] = der at 0
    //[2] = der at last
    //[3]... posX
    //[...]... poxY
    const int MAX_KNOTS = 20;
    static double* DLM_FITTER2_ARRAY_SPLINE3_X=NULL;
    static double* DLM_FITTER2_ARRAY_SPLINE3_Y=NULL;
    int NumKnots = TMath::Nint(pars[0]);
    if(NumKnots<2) NumKnots=2;
    if(NumKnots>MAX_KNOTS) NumKnots=MAX_KNOTS;
    if(!DLM_FITTER2_ARRAY_SPLINE3_X) DLM_FITTER2_ARRAY_SPLINE3_X = new double [MAX_KNOTS];
    if(!DLM_FITTER2_ARRAY_SPLINE3_Y) DLM_FITTER2_ARRAY_SPLINE3_Y = new double [MAX_KNOTS];
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        DLM_FITTER2_ARRAY_SPLINE3_X[iKnot] = pars[3+iKnot];
        DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot];
        //fix to the previous one of the value is fixed to 1e6
        if(DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot]==1e6&&iKnot) DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot-1];
    }
    double& derStart = pars[1];
    double& derEnd = pars[2];
    TSpline3 sp3("sp3", DLM_FITTER2_ARRAY_SPLINE3_X, DLM_FITTER2_ARRAY_SPLINE3_Y, NumKnots, "b1e1", derStart, derEnd);
    return sp3.Eval(*xVal);
}

double DlmTSplineFitPositive(double* xVal, double* pars){
    //[0] = NumKnots
    //[1] = der at 0
    //[2] = der at last
    //[3]... posX
    //[...]... poxY
    const int MAX_KNOTS = 20;
    static double* DLM_FITTER2_ARRAY_SPLINE3_X=NULL;
    static double* DLM_FITTER2_ARRAY_SPLINE3_Y=NULL;
    int NumKnots = TMath::Nint(pars[0]);
    if(NumKnots<2) NumKnots=2;
    if(NumKnots>MAX_KNOTS) NumKnots=MAX_KNOTS;
    if(!DLM_FITTER2_ARRAY_SPLINE3_X) DLM_FITTER2_ARRAY_SPLINE3_X = new double [MAX_KNOTS];
    if(!DLM_FITTER2_ARRAY_SPLINE3_Y) DLM_FITTER2_ARRAY_SPLINE3_Y = new double [MAX_KNOTS];
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        DLM_FITTER2_ARRAY_SPLINE3_X[iKnot] = pars[3+iKnot];
        DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot];
        //fix to the previous one of the value is fixed to 1e6
        if(DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot]==1e6&&iKnot) DLM_FITTER2_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot-1];
    }
    double& derStart = pars[1];
    double& derEnd = pars[2];
    TSpline3 sp3("sp3", DLM_FITTER2_ARRAY_SPLINE3_X, DLM_FITTER2_ARRAY_SPLINE3_Y, NumKnots, "b1e1", derStart, derEnd);
    double ReturnValue = sp3.Eval(*xVal);
    if(ReturnValue<0) ReturnValue=0;
    return ReturnValue;
}


TFitResultPtr DLM_FitHisto( TH1* h1, TF1* f1, Option_t* option, Option_t* goption,
                        const ROOT::Math::MinimizerOptions* moption,
                        Double_t xxmin, Double_t xxmax){

  // implementation of Fit method is in file hist/src/HFitImpl.cxx
  Foption_t fitOption;
  ROOT::Fit::FitOptionsMake(ROOT::Fit::EFitObjectType::kHistogram,option,fitOption);

  // create range and minimizer options with default values
  ROOT::Fit::DataRange range(xxmin,xxmax);
  ROOT::Math::MinimizerOptions minOption;

  // need to empty the buffer before
  // (t.b.d. do a ML unbinned fit with buffer data)
  h1->BufferEmpty();

  return ROOT::Fit::FitObject(h1, f1 , fitOption , moption?*moption:minOption, goption, range);
}

//DlmPoisson: x*mu/N, N (i.e. we can modify the shape for a given N)
double DlmPoisson(double x, double mu, double, double sigma, double Norm){
  return TMath::Poisson(x*sigma*sigma/mu, sigma*sigma);
}

//DlmPoisson: x*mu/N, N (i.e. we can modify the shape for a given N)
//par[0] = mu
//par[1] = sigma
//par[2] = Norm
double DlmPoisson(double* xVal, double* pars){
  return DlmPoisson(*xVal,pars[0],pars[1],pars[2]);
}
