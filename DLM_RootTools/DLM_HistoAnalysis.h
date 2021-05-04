

#ifndef DLM_HISTOANALYSIS_H
#define DLM_HISTOANALYSIS_H

#include "TMath.h"
#include "TH1D.h"

class TRandom3;

TH1D * MakeHistoFromTF1(TF1 * f1, double min, double max, unsigned int nbins = 100);
//gets the variance of a function in the interval [min, max]
//saves info about mean and mean2
double FunctionVar(TF1 * f1, double min, double max, double* mean=NULL, double* mean2=NULL, unsigned int steps = 100);

//Gets the central inteval of a functions, examines the function between min and max,
//alfa is the "amount" of integral that is desired, steps is the num of bins for the histogram,
//result is saved in from and to
//returns the median
void GetCentralInterval(TF1 * f1, const double& min, const double& max, double alpha, double& from, double& to, const unsigned int& steps = 100);
double GetCentralInterval(const TH1& h1, const double alpha, double& from, double& to,  bool extrapolate=false);
void DrawCentralInterval(TH1F*& h1, double Median, double LowReach, double UpReach, int NumNewBins,
                         TH1F*& h1main, TH1F*& h1bLow, TH1F*& h1bLow2, TH1F*& h1bUp, TH1F*& h1bUp2, TH1F*& h1median);
bool SameTH1structure(const TH1* h1, const TH1* h2);

//compute the PDF (i.e. errors) numeriacally, given: sum_i Norm_i * Num_i / Denom_i, where Num and Denom are histograms,
//which represent COUNTS! I.e. the uncertainties are Poisson.
class HistoAddRatios{
public:
   HistoAddRatios(const unsigned numterms, const char* OutHistoName, const unsigned numbins, const double ylow, const double yup);
   ~HistoAddRatios();
   //for larger number of counts in the histogram, we assume Gauss
   //void SetGaussLimit(const unsigned gauslim);
   //ignore the uncertaintiy of a certain entry
   void SetIgnoreUncertainty(const unsigned WhichOne, const bool yes=true);
   void SetNormalization(const unsigned WhichOne, const float Norm);
   void SetNormalization(const unsigned WhichOne, const TH1* Norm);
   //N.B. IF the numerator, or denumenator, has uncertainties, and the "ExpectationUncertainty" flag
   //is true, than in the sampling the mean used for the poisson is sampled from a Gaussian
   void SetNumerator(const unsigned WhichOne, const TH1* Num);
   void SetDenominator(const unsigned WhichOne, const TH1* Denom);
   //instead of ratio, we use a constant with a Gaussian uncertainty
   void SetConstant(const unsigned WhichOne, const float value, const float error);
   //we remove the denominator, in which case the numerator is treated as beeing the full ratio
   //the uncertainties are taken from the errors of the bins, and assumed Gaussian!
   void SetRatio(const unsigned WhichOne, const TH1* ratio, const bool MickeyPoisson=false);
   //in this case the ratio is taken from a TH2F, where the projection of eack bin is the pdf from which to
   //sample the uncertainty.
   void SetRatio(const unsigned WhichOne, TH2F* ratio);
   void SetNumIter(const unsigned numiter);
   void SetMinIter(const unsigned numiter);
   //if true, we will also spit out the significance of deviating from zero
   //if the expectation value is non-zero, just plug it in as extra entry with a minus sign
   //N.B. if true, and we have already achieved the required number of entries (MinEntriesForNsigma)
   //to estimate the significance, the iteration WILL terminate even if the NumIter is not yet reached!
   void SetCompareToNullHypothesis(const bool compare=true, const unsigned MinEntriesForNsigma=100);
   void SetRange(const float min, const float max);
   TH2F* GetResult();
   TH2F* CopyResult(TString HistoName);
   double GetTotPval(const bool Fisher=false);
   double GetTotNsig(const bool Fisher=false);
   double GetPval(const unsigned WhichBin);
   double GetNsig(const unsigned WhichBin);
   void SetRandomSeed(const int seed=0);
   void SetExpectationUncertainty(const bool experr);
   void GetCentralInt(const unsigned WhichBin, const double& nsigma, double& Median, double& Low, double& Up);
private:
  bool* IgnoreUncertainty;
  float* fNormalization;
  const TH1** hNormalization;
  const TH1** Numerator;
  const TH1** Denominator;
  TH2F** Ratio;
  float* ConstValue;
  float* ConstError;
  TH2F* hResult;
  const unsigned NumTerms;
  const unsigned NumBinsY;
  const double lowY;
  const double upY;
  //by default it is 10M
  unsigned NumIter;
  unsigned MinIter;
  //unsigned GaussLimit;
  bool CompareToNull;
  unsigned EntriesForNsigma;
  char* hResultName;
  TRandom3* rangen;
  double* pValue;
  double* nSigma;
  unsigned NbinsX;
  float Xmin;
  float Xmax;
  bool ExpectationUncertainty;
  bool* MickeyMousePoisson;
};

#endif
