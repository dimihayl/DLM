

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TRandom3.h"

#include "DLM_HistoAnalysis.h"

using namespace TMath;

bool SameTH1structure(const TH1* h1, const TH1* h2){
  if(h1->GetNbinsX()!=h2->GetNbinsX()) return false;
  for(int uBin=0; uBin<h1->GetNbinsX(); uBin++){
    if(h1->GetXaxis()->GetBinLowEdge(uBin+1)!=h2->GetXaxis()->GetBinLowEdge(uBin+1)) return false;
    if(h1->GetXaxis()->GetBinUpEdge(uBin+1)!=h2->GetXaxis()->GetBinUpEdge(uBin+1)) return false;
  }
  return true;
}

TH1D * MakeHistoFromTF1(TF1 * f1, double min, double max, unsigned int nbins){
    TH1D * h1d = new TH1D("h1d", "h1d", nbins, min, max);

    double x,y;

    for(unsigned int i=1; i<=nbins; i++){
        x = h1d->GetBinCenter(i);
        y = f1->Eval(x);
        h1d->SetBinContent(i, y);
    }

    return h1d;
}

double FunctionVar(TF1 * f1, double min, double max, double* mean, double* mean2, unsigned int steps){

    TH1D * h1d = MakeHistoFromTF1(f1, min, max, steps);

    double Mean = 0;
    double Mean2 = 0;

    Mean = 0;
    Mean2 = 0;

    double hintegral = h1d->Integral();

    for(unsigned int i=1; i<=steps; i++){
        Mean += h1d->GetBinCenter(i)*h1d->GetBinContent(i);
        Mean2  += Power(h1d->GetBinCenter(i), 2.)*h1d->GetBinContent(i);
    }
    Mean /= hintegral;
    Mean2  /= hintegral;

    if(mean && mean2){
        *mean = Mean;
        *mean2 = Mean2;
    }

    delete h1d;
    return Mean2 - Power(Mean, 2.);
}



void GetCentralInterval(TF1 * f1, const double& min, const double& max, double alpha, double& from, double& to, const unsigned int& steps){

    from = min;
    to = max;

    TH1D * h1d = MakeHistoFromTF1(f1, min, max, steps);

    double hintegral = h1d->Integral();
    double fintegral = f1->Integral(min, max);
    double myinterval = fintegral*(0.5-alpha/2.);
    double currentintegral = 0;

    h1d->Scale(fintegral/hintegral);

    for(unsigned int i=1; i<=steps; i++){
        currentintegral += h1d->GetBinContent(i);
        if(currentintegral>=myinterval){
            from = h1d->GetBinCenter(i);
            break;
        }
    }

    currentintegral = 0;
    for(unsigned int i=steps; i>0; i--){
        currentintegral += h1d->GetBinContent(i);
        if(currentintegral>=myinterval){
            to = h1d->GetBinCenter(i);
            break;
        }
    }

    delete h1d;

}



double GetCentralInterval(const TH1& h1, const double alpha, double& from, double& to, bool extrapolate){

    double hintegral = h1.Integral();
    double myinterval = hintegral*(0.5-alpha/2.);
    double currentintegral = 0;

    //Printf("hintegral = %f", hintegral);
    //Printf("myinterval = %f", myinterval);

    //h1->Scale(1/hintegral);
    int nbins = h1.GetNbinsX();

    for(int i=1; i<=nbins; i++){

        currentintegral += h1.GetBinContent(i);
        //double x2 = currentintegral;
        if(currentintegral>=myinterval){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f1 = currentintegral - h1.GetBinContent(i);
                double f2 = currentintegral;
                double fI = myinterval;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                from = (fI-b)/a;
                //Printf("x1 = %f", x1);
                //Printf("x2 = %f", x2);
                //Printf("a = %f", a);
                //Printf("b = %f", b);
                //Printf("myinterval = %f", myinterval);
            }
            else
                from = h1.GetBinLowEdge(i);
            break;
        }
    }
    currentintegral = 0;
    for(unsigned int i=nbins; i>0; i--){
        currentintegral += h1.GetBinContent(i);
        if(currentintegral>=myinterval){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f2 = hintegral-(currentintegral - h1.GetBinContent(i));
                double f1 = hintegral-currentintegral;
                double fI = hintegral-myinterval;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                to = (fI-b)/a;
            }
            else
                to = h1.GetBinLowEdge(i) + h1.GetBinWidth(i);
            break;
        }
    }

    currentintegral = 0;
    for(int i=1; i<=nbins; i++){
        currentintegral += h1.GetBinContent(i);
        if(currentintegral>=0.5*hintegral){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f1 = currentintegral - h1.GetBinContent(i);
                double f2 = currentintegral;
                double fI = 0.5*hintegral;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                return (fI-b)/a;
            }
            return h1.GetBinLowEdge(i);
        }
    }

    return 0;
}

void DrawCentralInterval(TH1F*& h1, double Median, double LowReach, double UpReach, int NumNewBins,
                         TH1F*& h1main, TH1F*& h1bLow, TH1F*& h1bLow2, TH1F*& h1bUp, TH1F*& h1bUp2, TH1F*& h1median){

    //double from = h1->GetBinLowEdge(h1->FindFirstBinAbove(0));
    //double to = h1->GetXaxis()->GetBinUpEdge(h1->FindLastBinAbove(0));
    double from = h1->GetBinLowEdge(1);
    double to = h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());
    int NumBins = h1->GetNbinsX();
    double NewBinRange = (to-from)/double(NumNewBins);
    //from -= NewBinRange;
    //to += NewBinRange;

    Printf("NumBins = %i", NumBins);
    Printf("from = %f", from);
    Printf("to = %f", to);

    h1main = new TH1F(TString::Format("%s_vis", h1->GetName()),
                          TString::Format("%s_vis", h1->GetName()), NumNewBins, from, to);

    int FirstNonZeroBin = h1main->FindBin(LowReach);
    //if(LowReach > h1main->GetBinCenter(FirstNonZeroBin)) FirstNonZeroBin++;
    int LastNonZeroBin = h1main->FindBin(UpReach);
    //if(UpReach < h1main->GetBinCenter(LastNonZeroBin)) LastNonZeroBin--;
    double FNZB_low = h1main->GetBinLowEdge(FirstNonZeroBin);
    double LNZB_up = h1main->GetXaxis()->GetBinUpEdge(LastNonZeroBin);

    //Printf("FirstNonZeroBin = %i; LastNonZeroBin = %i",FirstNonZeroBin, LastNonZeroBin);
    //Printf("%f : %f : %f : %f : %f : %f",from, FNZB_low, LowReach, UpReach, LNZB_up, to);


    h1bLow = new TH1F(TString::Format("%s_vis2l", h1->GetName()),
                          TString::Format("%s_vis2l", h1->GetName()), FirstNonZeroBin-1,
                             from, FNZB_low);
    h1bUp = new TH1F(TString::Format("%s_vis2u", h1->GetName()),
                          TString::Format("%s_vis2u", h1->GetName()), NumNewBins - LastNonZeroBin,
                            LNZB_up, to);
    h1bLow2 = new TH1F(TString::Format("%s_vis2l2", h1->GetName()),
                          TString::Format("%s_vis2l2", h1->GetName()), 1,
                              FNZB_low, LowReach);
    h1bUp2 = new TH1F(TString::Format("%s_vis2u2", h1->GetName()),
                          TString::Format("%s_vis2u2", h1->GetName()), 1,
                            UpReach, LNZB_up);


    h1median = new TH1F(TString::Format("%s_vis3", h1->GetName()),
                          TString::Format("%s_vis3", h1->GetName()), 1, Median-NewBinRange/8., Median+NewBinRange/8.);

    h1main->SetFillColor(h1->GetFillColor());
    h1main->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h1main->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());


    h1bLow->SetFillColor(kGray+1);
    h1bUp->SetFillColor(kGray+1);
    h1bLow2->SetFillColor(kGray+1);
    h1bUp2->SetFillColor(kGray+1);

    h1median->SetFillColor(kBlack);


    for(int i=1; i<=NumBins; i++){
        h1main->Fill(h1->GetBinCenter(i), h1->GetBinContent(i));
        //h1b->Fill(h1->GetBinCenter(i+1), h1->GetBinContent(i+1));
    }

    for(int i=1; i<FirstNonZeroBin; i++)
        h1bLow->SetBinContent(i, h1main->GetBinContent(i));
    h1bLow2->SetBinContent(1, h1main->GetBinContent(FirstNonZeroBin));
    h1bUp2->SetBinContent(1, h1main->GetBinContent(LastNonZeroBin));
    for(int i=LastNonZeroBin+1; i<=NumNewBins; i++)
        h1bUp->SetBinContent(i-LastNonZeroBin, h1main->GetBinContent(i));



    h1main->SetStats(false);
    h1main->SetTitle("");
    h1main->SetLineColor(h1main->GetFillColor());
    h1main->SetLineWidth(1);
    h1main->SetMarkerStyle(0);
    h1main->SetMarkerSize(0);
    h1main->GetYaxis()->SetRange(0, 1.05*h1main->GetMaximum());
    h1main->GetYaxis()->SetRangeUser(0, 1.05*h1main->GetMaximum());
    h1main->GetYaxis()->SetLimits(0, 1.05*h1main->GetMaximum());
    h1main->GetXaxis()->SetLabelSize(0.065);
    h1main->GetXaxis()->CenterTitle();
    h1main->GetXaxis()->SetTitleOffset(1.15);
    h1main->GetXaxis()->SetLabelOffset(0.02);
    h1main->GetXaxis()->SetTitleSize(0.075);
    h1main->GetYaxis()->SetLabelSize(0.065);
    h1main->GetYaxis()->CenterTitle();
    h1main->GetYaxis()->SetTitleOffset(1);
    h1main->GetYaxis()->SetTitleSize(0.075);
    h1main->GetYaxis()->SetNdivisions(505);

    h1bLow->SetStats(false);
    h1bLow->SetTitle("");
    h1bLow->SetLineColor(h1bLow->GetFillColor());
    h1bLow->SetLineWidth(1);
    h1bLow->SetMarkerStyle(0);
    h1bLow->SetMarkerSize(0);

    h1bLow2->SetStats(false);
    h1bLow2->SetTitle("");
    h1bLow2->SetLineColor(h1bLow2->GetFillColor());
    h1bLow2->SetLineWidth(1);
    h1bLow2->SetMarkerStyle(0);
    h1bLow2->SetMarkerSize(0);

    h1bUp->SetStats(false);
    h1bUp->SetTitle("");
    h1bUp->SetLineColor(h1bUp->GetFillColor());
    h1bUp->SetLineWidth(1);
    h1bUp->SetMarkerStyle(0);
    h1bUp->SetMarkerSize(0);

    h1bUp2->SetStats(false);
    h1bUp2->SetTitle("");
    h1bUp2->SetLineColor(h1bUp2->GetFillColor());
    h1bUp2->SetLineWidth(1);
    h1bUp2->SetMarkerStyle(0);
    h1bUp2->SetMarkerSize(0);

    h1median->SetStats(false);
    h1median->SetTitle("");
    h1median->SetLineColor(h1median->GetFillColor());
    h1median->SetLineWidth(1);
    h1median->SetMarkerStyle(0);
    h1median->SetMarkerSize(0);

    //for(int i=FirstNonZeroBin; i<=LastNonZeroBin; i++){
    //    h1b->SetBinContent(i,0);
    //}

    h1median->SetBinContent(1, h1main->GetBinContent(h1main->FindBin(Median)));

}

HistoAddRatios::HistoAddRatios(const unsigned numterms, const char* OutHistoName, const unsigned numbins, const double ylow, const double yup):
NumTerms(numterms),NumBinsY(numbins),lowY(ylow),upY(yup){
  IgnoreUncertainty = NULL;
  MickeyMousePoisson = NULL;
  fNormalization = NULL;
  hNormalization = NULL;
  Numerator = NULL;
  Denominator = NULL;
  Ratio = NULL;
  ConstValue = NULL;
  ConstError = NULL;
  pValue = NULL;
  nSigma = NULL;
  if(NumTerms){
    IgnoreUncertainty = new bool[NumTerms];
    MickeyMousePoisson = new bool[NumTerms];
    fNormalization = new float[NumTerms];
    hNormalization = new const TH1*[NumTerms];
    Numerator = new const TH1* [NumTerms];
    Denominator = new const TH1* [NumTerms];
    Ratio = new TH2F* [NumTerms];
    ConstValue = new float[NumTerms];
    ConstError = new float[NumTerms];
    for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
      IgnoreUncertainty[uTerm] = false;
      MickeyMousePoisson[uTerm] = false;
      fNormalization[uTerm] = 1;
      hNormalization[uTerm] = NULL;
      Numerator[uTerm] = NULL;
      Denominator[uTerm] = NULL;
      Ratio[uTerm] = NULL;
      ConstValue[uTerm] = 0;
      ConstError[uTerm] = 0;
    }
  }
  else{
    printf("\033[1;31mERROR:\033[0m Zero terms in HistoAddRatios\n");
  }
  hResult = NULL;
  NumIter = -1;
  MinIter = 4000;
  //GaussLimit = 32;
  CompareToNull = false;
  EntriesForNsigma = 100;
  hResultName = new char[128];
  strcpy(hResultName,OutHistoName);
  rangen = new TRandom3(11);
  NbinsX = 0;
  Xmin = -1e6;
  Xmax = 1e6;
  ExpectationUncertainty = false;
}
HistoAddRatios::~HistoAddRatios(){
  if(IgnoreUncertainty){delete[]IgnoreUncertainty;IgnoreUncertainty=NULL;}
  if(MickeyMousePoisson){delete[]MickeyMousePoisson;MickeyMousePoisson=NULL;}
  if(fNormalization){delete[]fNormalization;fNormalization=NULL;}
  if(hNormalization){delete[]hNormalization;hNormalization=NULL;}
  if(Numerator){delete[]Numerator;Numerator=NULL;}
  if(Denominator){delete[]Denominator;Denominator=NULL;}
  if(Ratio){delete[]Ratio;Ratio=NULL;}
  if(ConstValue){delete[]ConstValue;ConstValue=NULL;}
  if(ConstError){delete[]ConstError;ConstError=NULL;}
  if(hResult){delete hResult;hResult=NULL;}
  if(hResultName){delete hResultName;hResultName=NULL;}
  if(rangen){delete rangen;rangen=NULL;}
  if(pValue){delete pValue;pValue=NULL;}
  if(nSigma){delete nSigma;nSigma=NULL;}
}
//void HistoAddRatios::SetGaussLimit(const unsigned gauslim){
//  GaussLimit = gauslim;
//}
void HistoAddRatios::SetIgnoreUncertainty(const unsigned WhichOne, const bool yes){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetIgnoreUncertainty, only %u terms are allowed\n",NumTerms);
    return;
  }
  if(IgnoreUncertainty[WhichOne]!=yes&&hResult){delete hResult; hResult=NULL;}
  IgnoreUncertainty[WhichOne] = yes;
}
void HistoAddRatios::SetNormalization(const unsigned WhichOne, const float Norm){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetNormalization, only %u terms are allowed\n",NumTerms);
    return;
  }
  if(fNormalization[WhichOne]!=Norm&&hResult){delete hResult; hResult=NULL;}
  fNormalization[WhichOne] = Norm;
  hNormalization[WhichOne] = NULL;
}
void HistoAddRatios::SetNormalization(const unsigned WhichOne, const TH1* Norm){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetNormalization, only %u terms are allowed\n",NumTerms);
    return;
  }
  //make sure that we have the same binning in all histograms
  for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
    if(hNormalization[uTerm]){
      if(!SameTH1structure(Norm,Numerator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Denominator[uTerm]){
      if(!SameTH1structure(Norm,Denominator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
  }
  if(hResult){delete hResult; hResult=NULL;}
  hNormalization[WhichOne] = Norm;
}
void HistoAddRatios::SetNumerator(const unsigned WhichOne, const TH1* Num){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetNumerator, only %u terms are allowed\n",NumTerms);
    return;
  }
  //make sure that we have the same binning in all histograms
  for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
    if(Numerator[uTerm]){
      if(!SameTH1structure(Num,Numerator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Denominator[uTerm]){
      if(!SameTH1structure(Num,Denominator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Ratio[uTerm]){
      if(!SameTH1structure(Num,Ratio[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
  }
  if(hResult){delete hResult; hResult=NULL;}
  Numerator[WhichOne] = Num;
}
void HistoAddRatios::SetDenominator(const unsigned WhichOne, const TH1* Denom){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetDenominator, only %u terms are allowed\n",NumTerms);
    return;
  }
  //make sure that we have the same binning in all histograms
  for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
    if(Numerator[uTerm]){
      if(!SameTH1structure(Denom,Numerator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Denominator[uTerm]){
      if(!SameTH1structure(Denom,Denominator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Ratio[uTerm]){
      if(!SameTH1structure(Denom,Ratio[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
  }
  if(hResult){delete hResult; hResult=NULL;}
  Denominator[WhichOne] = Denom;
}
void HistoAddRatios::SetConstant(const unsigned WhichOne, const float value, const float error){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetConstant, only %u terms are allowed\n",NumTerms);
    return;
  }
  if(Numerator[WhichOne]){Numerator[WhichOne]=NULL;}
  if(Denominator[WhichOne]){Denominator[WhichOne]=NULL;}
  if(Ratio[WhichOne]){Ratio[WhichOne]=NULL;}
  ConstValue[WhichOne] = value;
  ConstError[WhichOne] = error;
}
void HistoAddRatios::SetRatio(const unsigned WhichOne, const TH1* ratio, const bool MickeyPoisson){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetRatio, only %u terms are allowed\n",NumTerms);
    return;
  }
  if(Denominator[WhichOne]){Denominator[WhichOne]=NULL;}
  SetNumerator(WhichOne,ratio);
  MickeyMousePoisson[WhichOne] = MickeyPoisson;
}
void HistoAddRatios::SetRatio(const unsigned WhichOne, TH2F* ratio){
  if(WhichOne>=NumTerms){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::SetRatio, only %u terms are allowed\n",NumTerms);
    return;
  }
  if(Denominator[WhichOne]){Denominator[WhichOne]=NULL;}
  if(Numerator[WhichOne]){Numerator[WhichOne]=NULL;}
  //make sure that we have the same binning in all histograms
  for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
    if(Numerator[uTerm]){
      if(!SameTH1structure(ratio,Numerator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Denominator[uTerm]){
      if(!SameTH1structure(ratio,Denominator[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
    if(Ratio[uTerm]){
      if(!SameTH1structure(ratio,Ratio[uTerm])){
        printf("\033[1;33mWARNING:\033[0m All histograms in HistoAddRatios should have the same binning\n");
        //return;
      }
      else break;
    }
  }
  Ratio[WhichOne] = ratio;
}
void HistoAddRatios::SetNumIter(const unsigned numiter){
  if(numiter<100){
    printf("\033[1;33mWARNING:\033[0m Very small number of iterations in HistoAddRatios (%u) was set, reconsider?\n",numiter);
  }
  if(NumIter!=numiter&&hResult){delete hResult; hResult=NULL;}
  NumIter = numiter;
}
void HistoAddRatios::SetMinIter(const unsigned numiter){
  if(MinIter!=numiter&&hResult){delete hResult; hResult=NULL;}
  MinIter = numiter;
}

void HistoAddRatios::SetCompareToNullHypothesis(const bool compare, const unsigned MinEntriesForNsigma){
  if(MinEntriesForNsigma<10){
    printf("\033[1;33mWARNING:\033[0m Very small number of entries required for the evaluation of nsigma in HistoAddRatios (%u) was set, reconsider?\n",MinEntriesForNsigma);
  }
  if(CompareToNull!=compare&&hResult){delete hResult; hResult=NULL;}
  if(CompareToNull&&EntriesForNsigma!=MinEntriesForNsigma&&hResult){delete hResult; hResult=NULL;}
  CompareToNull = compare;
  EntriesForNsigma = MinEntriesForNsigma;
}
void HistoAddRatios::SetRange(const float min, const float max){
  Xmin = min;
  Xmax = max;
}
TH2F* HistoAddRatios::GetResult(){
  if(!hResult){
    for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
      if(Numerator[uTerm]||Ratio[uTerm]){
        const TAxis* Xaxis = Numerator[uTerm]?Numerator[uTerm]->GetXaxis():Ratio[uTerm]->GetXaxis();
        NbinsX = Xaxis->GetNbins();
        double* BinningX = new double [NbinsX+1];
        for(unsigned uBin=0; uBin<NbinsX; uBin++){
          BinningX[uBin] = Xaxis->GetBinLowEdge(uBin+1);
        }
        BinningX[NbinsX] = Xaxis->GetBinUpEdge(NbinsX);
        //printf("1\n");
        hResult = new TH2F(hResultName,hResultName,NbinsX,BinningX,NumBinsY,lowY,upY);
        //printf("2\n");
        //hResult = (TH1F*)Numerator[uTerm]->Clone(hResultName);
        delete [] BinningX;
        for(int uBin=0; uBin<=hResult->GetNbinsX(); uBin++){
          hResult->SetBinContent(uBin,0);
          hResult->SetBinError(uBin,0);
        }
        break;
      }
    }
    if(!hResult){
      printf("\033[1;33mWARNING:\033[0m There are no entries set for the computation in HistoAddRatios\n");
      return NULL;
    }
    if(pValue){delete pValue;}
    if(nSigma){delete nSigma;}
    pValue = new double [NbinsX];
    nSigma = new double [NbinsX];
    for(unsigned uBin=0; uBin<NbinsX; uBin++){
      float Xvalue = hResult->GetXaxis()->GetBinCenter(uBin+1);
      if(Xvalue<Xmin) continue;
      if(Xvalue>Xmax) break;
      unsigned AboveZero=0;
      unsigned BelowZero=0;
      unsigned AtZero=0;
      unsigned UPVAL;
      TH1D** RatioProjection = new TH1D* [NumTerms];
      for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
        if(Ratio[uTerm]){
          RatioProjection[uTerm] = Ratio[uTerm]->ProjectionY(TString::Format("RatioProjection_%u",uTerm),uBin+1,uBin+1);
        }
        else{
          RatioProjection[uTerm] = NULL;
        }
      }
      for(unsigned uIter=0; uIter<NumIter; uIter++){
        double numer;
        double denom;
        double norm;
        double sum=0;
        for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
          if(hNormalization[uTerm]) norm = hNormalization[uTerm]->GetBinContent(uBin+1);
          else norm = fNormalization[uTerm];
          //division of yields
          if(Numerator[uTerm]&&Denominator[uTerm]){
            //no errors
            if(IgnoreUncertainty[uTerm]){
              sum += norm*Numerator[uTerm]->GetBinContent(uBin+1)/Denominator[uTerm]->GetBinContent(uBin+1);
            }//with errors
            else{
              double mean_numer = Numerator[uTerm]->GetBinContent(uBin+1);
              //if we want to use expected unertainties for the mean, this is where we do exactly that
              if(Numerator[uTerm]->GetBinError(uBin+1)&&ExpectationUncertainty)
                mean_numer = rangen->Gaus(mean_numer,Numerator[uTerm]->GetBinError(uBin+1));
              if(mean_numer<0) mean_numer = 0;
              //sampling from a Poisson. Below the same for the denom
              numer = rangen->Poisson(mean_numer);
              //printf(" numer = %.1f\n",numer);

              double mean_denom = Denominator[uTerm]->GetBinContent(uBin+1);
              if(Denominator[uTerm]->GetBinError(uBin+1)&&ExpectationUncertainty)
                mean_denom = rangen->Gaus(mean_denom,Denominator[uTerm]->GetBinError(uBin+1));
              if(mean_denom<0) mean_denom = 0;
              denom = rangen->Poisson(mean_denom);
              sum += norm*numer/denom;
            }
          }
          //the case in which we have only set the ratio
          else if(Numerator[uTerm]){
            //no errors
            if(IgnoreUncertainty[uTerm]){
              sum += norm*Numerator[uTerm]->GetBinContent(uBin+1);
            }
            //with errors
            //By Default we assume Gaussian errors, BUT if MickeyMousePoisson is true, we do the following trick:
            //assume that the uncertainty of the ratio is a scaled (by factor A) down Poisson (A*N), where the errors are A*sqrt(N).
            //we can than estimate the effective yield based on the value and error from the bin, and sample the uncertainties from there,
            //before finally rescaling by the corresponding factor A. This way, we avoid problems with negative values and should
            //get a more accurate result in case of investigating corr. functions, where the error is dominated by the numerator (Poisson)
            else{
              if(!MickeyMousePoisson[uTerm]){
                numer = rangen->Gaus(Numerator[uTerm]->GetBinContent(uBin+1),Numerator[uTerm]->GetBinError(uBin+1));
              }
              else{
                double Value = Numerator[uTerm]->GetBinContent(uBin+1);
                double Error = Numerator[uTerm]->GetBinError(uBin+1);
                double EffYield = pow(Value/Error,2.);
                numer = Value/EffYield*rangen->Poisson(EffYield);
              }
              sum += norm*numer;
            }
          }
          else if(Ratio[uTerm]){
            numer = RatioProjection[uTerm]->GetRandom();
            sum += norm*numer;
          }
          //constant
          else{
            if(IgnoreUncertainty[uTerm]){
              sum += norm*ConstValue[uTerm];
            }
            else{
              sum += norm*rangen->Gaus(ConstValue[uTerm],ConstError[uTerm]);
            }
          }
        }
        hResult->Fill(Xvalue,sum);
        pValue[uBin] = 0;
        nSigma[uBin] = 0;
        if(CompareToNull){
          AboveZero += sum>0;
          BelowZero += sum<0;
          AtZero += sum==0;
          if( (AboveZero>=EntriesForNsigma&&BelowZero>=EntriesForNsigma&&uIter+1>=MinIter)||uIter+1==NumIter ){
            if(AboveZero>BelowZero) UPVAL=2.*BelowZero+AtZero;
            else UPVAL=2.*AboveZero+AtZero;
            pValue[uBin] = double(UPVAL)/double(uIter+1);
            //printf(" %u raw pval = %.5f\n",uBin,pValue[uBin]);
            if(pValue[uBin]) nSigma[uBin] = sqrt(2)*TMath::ErfcInverse(pValue[uBin]);
            //printf(" %u raw nsig = %.5f\n",uBin,nSigma[uBin]);
            break;//the uIter
          }
        }
      }//uIter
      for(unsigned uTerm=0; uTerm<NumTerms; uTerm++){
        delete RatioProjection[uTerm];
        RatioProjection[uTerm] = NULL;
      }
      delete [] RatioProjection;
      RatioProjection = NULL;
    }//uBin
  }
  return hResult;
}
TH2F* HistoAddRatios::CopyResult(TString HistoName){
  if(!hResult) GetResult();
  return (TH2F*)hResult->Clone(HistoName);
}
double HistoAddRatios::GetTotPval(const bool Fisher){
  double Chi2=0;
  if(Fisher){
    for(unsigned uBin=0; uBin<NbinsX; uBin++){
      Chi2 -= 2.*log(pValue[uBin]);
    }
    return TMath::Prob(Chi2,2*NbinsX);
  }
  else{
    for(unsigned uBin=0; uBin<NbinsX; uBin++){
      Chi2 += nSigma[uBin]*nSigma[uBin];
    }
    return TMath::Prob(Chi2,NbinsX);
  }
}
double HistoAddRatios::GetTotNsig(const bool Fisher){
  return sqrt(2)*TMath::ErfcInverse(GetTotPval(Fisher));
}
double HistoAddRatios::GetPval(const unsigned WhichBin){
  if(WhichBin>=NbinsX){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::GetPval, only %u bins are allowed\n",NbinsX);
    return 0;
  }
  if(!pValue[WhichBin]&&!nSigma[WhichBin]){
    //return TMath::Erfc(GetNsig(WhichBin)/sqrt(2));
    return -1./double(NumIter);//upper limit
  }
  else return pValue[WhichBin];
}
double HistoAddRatios::GetNsig(const unsigned WhichBin){
  if(WhichBin>=NbinsX){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::GetNsig, only %u bins are allowed\n",NbinsX);
    return 0;
  }
  if(!pValue[WhichBin]&&!nSigma[WhichBin]){
    //TH1D* hProj = hResult->ProjectionY("ProjectionY",WhichBin+1,WhichBin+1);
    //return fabs(hProj->GetMean()/hProj->GetStdDev());
    return -sqrt(2)*TMath::ErfcInverse(1./double(NumIter));
  }
  else return nSigma[WhichBin];
}
void HistoAddRatios::SetRandomSeed(const int seed){
  if(rangen){rangen->SetSeed(seed);}
}
void HistoAddRatios::SetExpectationUncertainty(const bool experr){
  ExpectationUncertainty = experr;
}
void HistoAddRatios::GetCentralInt(const unsigned WhichBin, const double& nsigma, double& Median, double& Low, double& Up){
  if(WhichBin>=NbinsX){
    printf("\033[1;33mWARNING:\033[0m Bad input in HistoAddRatios::GetCentralInt, only %u bins are allowed\n",NbinsX);
    return;
  }
  TH1D* hProj = hResult->ProjectionY("ProjectionY",WhichBin+1,WhichBin+1);
  Median = GetCentralInterval(hProj[0],1.-TMath::Erfc(nsigma/sqrt(2)),Low,Up,true);
}
