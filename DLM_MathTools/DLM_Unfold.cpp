
#include "DLM_Unfold.h"
#include "DLM_CppTools.h"
#include "DLM_RootWrapper.h"
#include "DLM_Histo.h"

#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TSpline.h"
#include "TROOT.h"

#include <iostream>

DLM_Unfold::DLM_Unfold(){
  NumIter = 1;
  Precision = 0.2;
  ReachedPrecision = 0;
  dlmData = NULL;
  dlmData_reb2 = NULL;
  dlmResponse = NULL;
  //dlmResponseReb = NULL;
  hData = NULL;
  hResponse = NULL;
  //15 mins
  MaxTime = 900;
  RangeMinF = 0;
  RangeMaxF = 0;
  RangeMinU = 0;
  RangeMaxU = 0;
  SEEDmin = 0;
  Silent = false;
  DLM_FITTER_ARRAY_SPLINE3_X = NULL;
  DLM_FITTER_ARRAY_SPLINE3_Y = NULL;
  SparseFirst = NULL;
  SparseLast = NULL;
  OutputLevel = 0;
  OutputFileName= "";
  FitFoldOriginal = NULL;
  FitFoldFinal = NULL;
  FitUnfoldFinal = NULL;
  dlmFoldedWorkHorse = NULL;
  dlmUnfoldedWorkHorse = NULL;
  dlmFoldedDimiWorkHorse = NULL;
  dlmUnfoldedDimiWorkHorse = NULL;
  dlmDimiWorkHorse = NULL;
  PcOutput = true;
  BinByBinFolded = NULL;
  BinByBinUnfolded = NULL;
}

DLM_Unfold::~DLM_Unfold(){
  gROOT->cd();
  if(dlmData){delete dlmData; dlmData=NULL;}
  if(dlmData_reb2){delete dlmData_reb2; dlmData_reb2=NULL;}
  if(dlmResponse){delete dlmResponse; dlmResponse=NULL;}
  //if(dlmResponseReb){delete dlmResponseReb; dlmResponseReb=NULL;}
  if(DLM_FITTER_ARRAY_SPLINE3_X){delete DLM_FITTER_ARRAY_SPLINE3_X; DLM_FITTER_ARRAY_SPLINE3_X=NULL;}
  if(DLM_FITTER_ARRAY_SPLINE3_Y){delete DLM_FITTER_ARRAY_SPLINE3_Y; DLM_FITTER_ARRAY_SPLINE3_Y=NULL;}
  if(SparseFirst){delete[]SparseFirst; SparseFirst=NULL;}
  if(SparseLast){delete[]SparseLast; SparseLast=NULL;}
  if(FitFoldOriginal){delete FitFoldOriginal; FitFoldOriginal=NULL;}
  if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}
  if(FitUnfoldFinal){delete FitUnfoldFinal; FitUnfoldFinal=NULL;}
  if(dlmFoldedWorkHorse){delete dlmFoldedWorkHorse; dlmFoldedWorkHorse=NULL;}
  if(dlmUnfoldedWorkHorse){delete dlmUnfoldedWorkHorse; dlmUnfoldedWorkHorse=NULL;}
  if(dlmFoldedDimiWorkHorse){delete dlmFoldedDimiWorkHorse; dlmFoldedDimiWorkHorse=NULL;}
  if(dlmUnfoldedDimiWorkHorse){delete dlmUnfoldedDimiWorkHorse; dlmUnfoldedDimiWorkHorse=NULL;}
  if(dlmDimiWorkHorse){delete dlmDimiWorkHorse; dlmDimiWorkHorse=NULL;}
  if(hData){delete hData; hData=NULL;}
  if(hResponse){delete hResponse; hResponse=NULL;}
  if(BinByBinFolded){delete BinByBinFolded; BinByBinFolded=NULL;}
  if(BinByBinUnfolded){delete BinByBinUnfolded; BinByBinUnfolded=NULL;}
}

void DLM_Unfold::SetData(TH1F* data){
  //if(data.GetDim()!=1){
  //  printf("\033[1;33mWARNING!\033[0m DLM_Unfold got a non 1D data! Not allowed!\n");
  //  return;
  //}
  gROOT->cd();
  if(hData){delete hData; hData=NULL;}
  hData = (TH1F*)data->Clone(TString::Format("DLM_Unfold::hData_%p",this));
  if(dlmData){delete dlmData; dlmData=NULL;}
  //if(dlmResponseReb){delete dlmResponseReb; dlmResponseReb=NULL;}
  if(!data) return;
  dlmData = Convert_TH1F_DlmHisto(hData);
  if(RangeMinF==0&&RangeMaxF==0){
    RangeMinF = dlmData->GetBinLowEdge(0,0);
    RangeMaxF = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
  else if(RangeMinF<dlmData->GetBinLowEdge(0,0)){
    RangeMinF = dlmData->GetBinLowEdge(0,0);
  }
  else if(RangeMaxF>dlmData->GetBinLowEdge(0,dlmData->GetNbins()-1)){
    RangeMaxF = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }

  if(RangeMinU==0&&RangeMaxU==0){
    RangeMinU = dlmData->GetBinLowEdge(0,0);
    RangeMaxU = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
  else if(RangeMinU<dlmData->GetBinLowEdge(0,0)){
    RangeMinU = dlmData->GetBinLowEdge(0,0);
  }
  else if(RangeMaxU>dlmData->GetBinLowEdge(0,dlmData->GetNbins()-1)){
    RangeMaxU = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
}

void DLM_Unfold::SetData(DLM_Histo<float>* data){
  gROOT->cd();
  if(hData){delete hData; hData=NULL;}
  if(dlmData){delete dlmData; dlmData=NULL;}
  if(!data) return;
  dlmData = new DLM_Histo<float>(*data);
  hData = Convert_DlmHisto_TH1F(dlmData,TString::Format("DLM_Unfold::hData_%p",this));
  if(RangeMinF==0&&RangeMaxF==0){
    RangeMinF = dlmData->GetBinLowEdge(0,0);
    RangeMaxF = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
  else if(RangeMinF<dlmData->GetBinLowEdge(0,0)){
    RangeMinF = dlmData->GetBinLowEdge(0,0);
  }
  else if(RangeMaxF>dlmData->GetBinLowEdge(0,dlmData->GetNbins()-1)){
    RangeMaxF = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }

  if(RangeMinU==0&&RangeMaxU==0){
    RangeMinU = dlmData->GetBinLowEdge(0,0);
    RangeMaxU = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
  else if(RangeMinU<dlmData->GetBinLowEdge(0,0)){
    RangeMinU = dlmData->GetBinLowEdge(0,0);
  }
  else if(RangeMaxU>dlmData->GetBinLowEdge(0,dlmData->GetNbins()-1)){
    RangeMaxU = dlmData->GetBinUpEdge(0,dlmData->GetNbins()-1);
  }
}

void DLM_Unfold::SetResponse(TH2F* response){
  //if(!response){dlmResponse = NULL; return;}
  //if(response->GetDim()!=2){
  //  printf("\033[1;33mWARNING!\033[0m DLM_Unfold got a non 2D matrix! Not allowed!\n");
  //  return;
  //}
  gROOT->cd();
  if(hResponse){delete hResponse; hResponse=NULL;}
  hResponse = (TH2F*)response->Clone(TString::Format("DLM_Unfold::hResponse_%p",this));
  if(dlmResponse){delete dlmResponse; dlmResponse=NULL;}
  if(!response) return;
  dlmResponse = Convert_TH2F_DlmHisto(hResponse);

  if(SparseFirst){delete SparseFirst; SparseFirst=NULL;}
  if(SparseLast){delete SparseLast; SparseLast=NULL;}
  SparseFirst = new unsigned [response->GetYaxis()->GetNbins()];
  SparseLast = new unsigned [response->GetYaxis()->GetNbins()];
  for(unsigned uBinY=0; uBinY<response->GetYaxis()->GetNbins(); uBinY++){
    SparseFirst[uBinY]=-1;
    SparseLast[uBinY]=0;
    //bool FirstFound = false;
    for(unsigned uBinX=0; uBinX<response->GetXaxis()->GetNbins(); uBinX++){
      if(response->GetBinContent(uBinX+1,uBinY+1)!=0&&SparseFirst[uBinY]==-1){
        SparseFirst[uBinY]=uBinX;
      }
      if(response->GetBinContent(uBinX+1,uBinY+1)!=0&&SparseLast[uBinY]<=response->GetXaxis()->GetNbins()){
        SparseLast[uBinY]=uBinX;
      }
    }
    //this is done so that we have a single bin left for extrapolation
    if(SparseFirst[uBinY])SparseFirst[uBinY]--;
    if(SparseLast[uBinY]<response->GetXaxis()->GetNbins()-1)SparseLast[uBinY]++;
    //printf("Y %u: %u %u\n",uBinY,SparseFirst[uBinY],SparseLast[uBinY]);
  }
}


void DLM_Unfold::SetResponse(DLM_Histo<float>* response){
  gROOT->cd();
  if(hResponse){delete hResponse; hResponse=NULL;}
  if(dlmResponse){delete dlmResponse; dlmResponse=NULL;}
  if(!response) return;

  dlmResponse = new DLM_Histo<float>(*response);
  hResponse = Convert_DlmHisto_TH2F(dlmResponse,TString::Format("DLM_Unfold::hResponse_%p",this));

  if(SparseFirst){delete SparseFirst; SparseFirst=NULL;}
  if(SparseLast){delete SparseLast; SparseLast=NULL;}
  SparseFirst = new unsigned [hResponse->GetYaxis()->GetNbins()];
  SparseLast = new unsigned [hResponse->GetYaxis()->GetNbins()];
  for(unsigned uBinY=0; uBinY<hResponse->GetYaxis()->GetNbins(); uBinY++){
    SparseFirst[uBinY]=-1;
    SparseLast[uBinY]=0;
    //bool FirstFound = false;
    for(unsigned uBinX=0; uBinX<hResponse->GetXaxis()->GetNbins(); uBinX++){
      if(hResponse->GetBinContent(uBinX+1,uBinY+1)!=0&&SparseFirst[uBinY]==-1){
        SparseFirst[uBinY]=uBinX;
      }
      if(hResponse->GetBinContent(uBinX+1,uBinY+1)!=0&&SparseLast[uBinY]<=hResponse->GetXaxis()->GetNbins()){
        SparseLast[uBinY]=uBinX;
      }
    }
    //this is done so that we have a single bin left for extrapolation
    if(SparseFirst[uBinY])SparseFirst[uBinY]--;
    if(SparseLast[uBinY]<hResponse->GetXaxis()->GetNbins()-1)SparseLast[uBinY]++;
    //printf("Y %u: %u %u\n",uBinY,SparseFirst[uBinY],SparseLast[uBinY]);
  }
}


/*
void DLM_Unfold::SparseTheResponse(){
  unsigned uBin[2];
  if(SparseFirst){delete SparseFirst; SparseFirst=NULL;}
  if(SparseLast){delete SparseLast; SparseLast=NULL;}
  SparseFirst = new unsigned [dlmResponseReb->GetNbins(1)];
  SparseLast = new unsigned [dlmResponseReb->GetNbins(1)];
  for(uBin[1]=0; uBin[1]<dlmResponseReb->GetNbins(1); uBin[1]++){
    SparseFirst[uBin[1]]=-1;
    SparseLast[uBin[1]]=0;
    //bool FirstFound = false;
    for(uBin[0]=0; uBin[0]<dlmResponseReb->GetNbins(0); uBin[0]++){
      if(dlmResponseReb->GetBinContent(uBin)!=0&&SparseFirst[uBin[1]]==-1){
        SparseFirst[uBin[1]]=uBin[0];
      }
      if(dlmResponseReb->GetBinContent(uBin)!=0&&SparseLast[uBin[1]]<=dlmResponseReb->GetNbins(0)){
        SparseLast[uBin[1]]=uBin[0];
      }
    }
    //this is done so that we have a single bin left for extrapolation
    if(SparseFirst[uBin[1]])SparseFirst[uBin[1]]--;
    if(SparseLast[uBin[1]]<<dlmResponseReb->GetNbins(0)-1)SparseLast[uBin[1]]++;
    //printf("Y %u: %u %u\n",uBin[1],SparseFirst[uBin[1]],SparseLast[uBin[1]]);
  }
}
*/
void DLM_Unfold::SetUnfoldPrecision(const double& precision, const double& worst){
  if(precision<=0||precision>1||worst<precision){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold got a weird precision! Ignored!\n");
    return;
  }
  Precision = precision;
  PrecisionWorst = worst;
}

void DLM_Unfold::SetUnfoldBootstrap(const unsigned& numiter, const int& seedmin){
  if(numiter>1000){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold got an overkill (%u) in terms of number of iterations!\n",numiter);
  }
  NumIter = numiter;
  SEEDmin = seedmin;
}

void DLM_Unfold::SetUnfoldSeconds(const double& seconds){
  MaxTime = seconds;
}

void DLM_Unfold::SetUnfoldMinutes(const double& minutes){
  MaxTime = minutes*60.;
}

void DLM_Unfold::SetFoldRange(const double& min, const double& max){
  if(min>max){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold got a weird range! Ignored!\n");
    return;
  }
  RangeMinF = min;
  RangeMaxF = max;
}
void DLM_Unfold::SetUnfoldRange(const double& min, const double& max){
  if(min>max){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold got a weird range! Ignored!\n");
    return;
  }
  RangeMinU = min;
  RangeMaxU = max;
}
void DLM_Unfold::SetSilentMode(const bool& silentmode){
  Silent = silentmode;
}

//void DLM_Unfold::SetOutput(const int& level, const TString& ouputfilename){
//  Level = level;
//  OutputFileName = ouputfilename;
//  if(Level>0){
//
//  }
//  TFile fOutput();
//}



double DLM_Unfold::SplineFitFolded(double* xVal, double* pars){
  static double xLast = *xVal;
  double& xCurrent = *xVal;
  //we test new parameters
  if(xCurrent<=xLast){
    if(!dlmUnfoldedWorkHorse){
      dlmUnfoldedWorkHorse = new DLM_Histo<float>(*dlmData);
    }
    for(unsigned uBin=0; uBin<dlmUnfoldedWorkHorse->GetNbins(); uBin++){
      double BINCEN = dlmUnfoldedWorkHorse->GetBinCenter(0,uBin);
      dlmUnfoldedWorkHorse->SetBinContent(uBin,SplineFit(&BINCEN,pars));
    }
    //if(dlmFoldedWorkHorse){delete dlmFoldedWorkHorse; dlmFoldedWorkHorse=NULL;}
    if(!dlmFoldedWorkHorse) dlmFoldedWorkHorse = new DLM_Histo<float> (*dlmUnfoldedWorkHorse);
    Fold(*dlmUnfoldedWorkHorse,*dlmFoldedWorkHorse);
  }
  return dlmFoldedWorkHorse->Eval(&xCurrent);
}
double DLM_Unfold::SplineFitUnfolded(double* xVal, double* pars){
  static double xLast = *xVal;
  double& xCurrent = *xVal;
  //we test new parameters
  if(xCurrent<=xLast){
    if(!dlmUnfoldedWorkHorse){
      dlmUnfoldedWorkHorse = new DLM_Histo<float>(*dlmData);
    }
    for(unsigned uBin=0; uBin<dlmUnfoldedWorkHorse->GetNbins(); uBin++){
      double BINCEN = dlmUnfoldedWorkHorse->GetBinCenter(0,uBin);
      dlmUnfoldedWorkHorse->SetBinContent(uBin,SplineFit(&BINCEN,pars));
    }
  }
  return dlmUnfoldedWorkHorse->Eval(&xCurrent);
}

double DLM_Unfold::SplineFit(double* xVal, double* pars){
    //[0] = NumKnots
    //[1] = der at 0
    //[2] = der at last
    //[3]... posX
    //[...]... poxY
    const int MAX_KNOTS = 100;
    int NumKnots = TMath::Nint(pars[0]);
    if(NumKnots<2) NumKnots=2;
    if(NumKnots>MAX_KNOTS) NumKnots=MAX_KNOTS;
    if(!DLM_FITTER_ARRAY_SPLINE3_X) DLM_FITTER_ARRAY_SPLINE3_X = new double [MAX_KNOTS];
    if(!DLM_FITTER_ARRAY_SPLINE3_Y) DLM_FITTER_ARRAY_SPLINE3_Y = new double [MAX_KNOTS];
    for(int iKnot=0; iKnot<NumKnots; iKnot++){
        DLM_FITTER_ARRAY_SPLINE3_X[iKnot] = pars[3+iKnot];
        DLM_FITTER_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot];
        //fix to the previous one of the value is fixed to 1e6
        if(DLM_FITTER_ARRAY_SPLINE3_Y[iKnot]==1e6&&iKnot) DLM_FITTER_ARRAY_SPLINE3_Y[iKnot] = pars[3+NumKnots+iKnot-1];
    }
    double& derStart = pars[1];
    double& derEnd = pars[2];
    TSpline3 sp3("sp3", DLM_FITTER_ARRAY_SPLINE3_X, DLM_FITTER_ARRAY_SPLINE3_Y, NumKnots, "b1e1", derStart, derEnd);
    return sp3.Eval(*xVal);
}





double DLM_Unfold::DimiFit(double* xVal, double* pars){
  static double xLast = *xVal;
  double& xCurrent = *xVal;
  if(!dlmDimiWorkHorse) return 0;
  //we test new parameters
  if(xCurrent<=xLast){
    for(unsigned uBin=0; uBin<dlmDimiWorkHorse->GetNbins(); uBin++){
      //if(dlmDimiWorkHorse->GetBinContent(uBin)!=pars[uBin]){
        //printf("at %.2f: %.4f -> %.4f\n",
        //dlmDimiWorkHorse->GetBinCenter(0,uBin),dlmDimiWorkHorse->GetBinContent(uBin),pars[uBin]);
      //}
      dlmDimiWorkHorse->SetBinContent(uBin,pars[uBin]);
    }
  }
  //printf("at %.2f: %.2f\n",*xVal,dlmDimiWorkHorse->Eval(xVal));
  return dlmDimiWorkHorse->Eval(xVal);
}
double DLM_Unfold::DimiFitFolded(double* xVal, double* pars){
  static double xLast = *xVal;
  double& xCurrent = *xVal;
  if(!dlmDimiWorkHorse) return 0;
  //we test new parameters
  if(xCurrent<=xLast){
    for(unsigned uBin=0; uBin<dlmDimiWorkHorse->GetNbins(); uBin++){
      dlmDimiWorkHorse->SetBinContent(uBin,pars[uBin]);
    }
    if(!dlmFoldedDimiWorkHorse) dlmFoldedDimiWorkHorse = new DLM_Histo<float> (*dlmDimiWorkHorse);
    Fold(*dlmDimiWorkHorse,*dlmFoldedDimiWorkHorse);
  }
  //printf("%f vs %f\n",dlmFoldedDimiWorkHorse->Eval(xVal),dlmDimiWorkHorse->Eval(xVal));
  return dlmFoldedDimiWorkHorse->Eval(xVal);
}

//TYPE: 0 - the original (splines)
//TYPE: 1 - DimiStyle, experiment with histo gone wrong
//TYPE: 2 - fit with a function that should be suited in describing ME
DLM_Histo<float>* DLM_Unfold::Unfold(const int& TYPE){
  if(!dlmData||!dlmResponse){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold misses some input to unfold!\n");
    return NULL;
  }

  double YieldScale=1;
  for(unsigned uBin=0; uBin<dlmData->GetNbins(); uBin++){
    if(dlmData->GetBinContent(uBin)){
      double BinCont = dlmData->GetBinContent(uBin);
      double BinErr = dlmData->GetBinError(uBin);
      YieldScale = BinCont/BinErr/BinErr;
      break;
    }
  }


  bool TerminateASAP = false;
  DLM_Timer JobTime;
  TRandom3 rangen(SEEDmin+1);

  const unsigned NumBins = dlmData->GetBin(0,RangeMaxU)-dlmData->GetBin(0,RangeMinU)+1;
  if(NumBins<=1){
    printf("\033[1;33mWARNING!\033[0m Smearing a less than 2 bins is truely the next level! Sorry, cannot do it!\n");
    return NULL;
  }

  if(BinByBinFolded){delete [] BinByBinFolded; BinByBinFolded=NULL;}
  BinByBinFolded = new std::vector<float> [NumBins];
  if(BinByBinUnfolded){delete [] BinByBinUnfolded; BinByBinUnfolded=NULL;}
  BinByBinUnfolded = new std::vector<float> [NumBins];
  DLM_Histo<float>* dlmUnfoldResult = new DLM_Histo<float>(*dlmData);
  dlmUnfoldResult->SetBinContentAll(0);
  dlmUnfoldResult->SetBinErrorAll(0);
  unsigned NumResultEnties = 0;

  double Avg_chi2ndf=0;
  //const unsigned BinDepth = NumBins<5?NumBins:5;
  const unsigned BinDepth = 2;
  //the time (in seconds) beyond which we settle for Okayish_chi2ndf;
  long long PatienceForPerfection;
  //whatever we have after that much time, we take
  long long PatienceMax;
  const unsigned NumBackAndForth = 1024*4;//this is a maximum value
  const unsigned MaxStepsWithoutImprovement1 = 24;//32
  const unsigned MaxStepsWithoutImprovement2 = 12;//16
  const unsigned ErrorDepth = 3;//3

  char* strtime = new char [128];
  char* ctext1 = new char [128];
  char* ctext2 = new char [128];
  char* ctext3 = new char [128];

  //const int NumKnots = NumBins<10?NumBins/2+1:NumBins<20?NumBins/3+1:NumBins/4+1;
  //double* Knot_x = new double [NumKnots];
  //for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
  //  Knot_x[uKnot] = RangeMinU + double(uKnot)*(RangeMaxU-RangeMinU)/double(NumKnots);
  //}
  //const int NumKnots = NumBins/3+1;
  //const int NumKnots = NumBins<10?NumBins/2+1:NumBins<20?NumBins/3+1:NumBins/4+1;
  //double* Knot_x = new double [NumKnots];
  //for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
  //  Knot_x[uKnot] = RangeMinU + double(uKnot)*(RangeMaxU-RangeMinU)/double(NumKnots);
  //}

  //we take: the low edge, the up edge, the first three bin centers, and than each third bin center
  int NumKnots=2;
  for(unsigned uBin=0; uBin<dlmData->GetNbins(); uBin++){
    if(dlmData->GetBinCenter(0,uBin)<RangeMinU||dlmData->GetBinCenter(0,uBin)>RangeMaxU) continue;
    if(uBin<=2){NumKnots++;}
    else if((uBin-2)%3==0){NumKnots++;}
  }
  double* Knot_x = new double [NumKnots];
  Knot_x[0] = RangeMinU;
  NumKnots=1;
  for(unsigned uBin=0; uBin<dlmData->GetNbins(); uBin++){
    if(dlmData->GetBinCenter(0,uBin)<RangeMinU||dlmData->GetBinCenter(0,uBin)>RangeMaxU) continue;
    if(uBin<=2){Knot_x[NumKnots]=dlmData->GetBinCenter(0,uBin); NumKnots++;}
    else if((uBin-2)%3==0){Knot_x[NumKnots]=dlmData->GetBinCenter(0,uBin); NumKnots++;}
  }
  Knot_x[NumKnots]=RangeMaxU;
  NumKnots++;

  if(!Silent) printf("\nUnfolding:\n\n");
  if(!Silent) printf(" NumKnots = %u\n",NumKnots);
  double IterTime = 0;

  double f_chi2;
  int f_ndf;
  double f_chi2ndf;
  double f_pval;
  double f_ns;

  if(TYPE==1){
    if(FitFoldOriginal){delete FitFoldOriginal; FitFoldOriginal=NULL;}
    gROOT->cd();
    const unsigned REB2 = 1;
    if(dlmData_reb2){delete dlmData_reb2; dlmData_reb2=NULL;}
    dlmData_reb2 = new DLM_Histo<float>(*dlmData);
    dlmData_reb2->Rebin(&REB2);
    if(dlmDimiWorkHorse){delete dlmDimiWorkHorse; dlmDimiWorkHorse=NULL;}
    dlmDimiWorkHorse = new DLM_Histo<float>(*dlmData_reb2);

    FitFoldOriginal = new TF1("FitFoldOriginal",this,&DLM_Unfold::DimiFit,RangeMinU,RangeMaxU,dlmData_reb2->GetNbins(),"DLM_Unfold","DimiFit");
    for(unsigned uPar=0; uPar<dlmData_reb2->GetNbins(); uPar++){
      FitFoldOriginal->SetParameter(uPar,dlmData_reb2->GetBinContent(uPar));
      FitFoldOriginal->SetParLimits(uPar,0,2.*dlmData_reb2->GetBinContent(uPar));
      if(dlmData_reb2->GetBinCenter(0,uPar)<RangeMinU||dlmData_reb2->GetBinCenter(0,uPar)>RangeMaxU){
        FitFoldOriginal->FixParameter(uPar,dlmData_reb2->GetBinContent(uPar));
      }
      //printf("p%u = %f\n",uPar,dlmData_reb2->GetBinContent(uPar));
    }
    FitFoldOriginal->SetNpx(256);
    gROOT->cd();
    //printf("hData 0 = %f +/- .3%f\n",hData->GetBinContent(uBin+1),hData->GetBinError(uBin+1));
    hData->Fit(FitFoldOriginal,"S, N, R, M");

    f_chi2 = FitFoldOriginal->GetChisquare();
    f_ndf = FitFoldOriginal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf("Pre-fit (folded):\n");
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}
    gROOT->cd();

    //if(dlmDimiWorkHorse){delete [] dlmDimiWorkHorse; dlmDimiWorkHorse=NULL;}
    //dlmDimiWorkHorse = new DLM_Histo<float>(*dlmData_reb2);
    FitFoldFinal = new TF1("FitFoldFinal",this,&DLM_Unfold::DimiFitFolded,RangeMinU,RangeMaxU,dlmData_reb2->GetNbins(),"DLM_Unfold","DimiFitFolded");

    if(!Silent) printf("Pre-fit:\n");
    f_chi2*=100;
for(unsigned uRefit=0; uRefit<1; uRefit++){
    //for(unsigned uRefit=0; uRefit<11; uRefit++){
      double randomize=1;
      if(uRefit){
        randomize = rangen.Uniform(0.9,1.1);
      }
      for(unsigned uPar=0; uPar<dlmData_reb2->GetNbins(); uPar++){
        //FitFoldFinal->FixParameter(uPar,dlmData_reb2->GetBinContent(uPar));
        FitFoldFinal->SetParameter(uPar,FitFoldOriginal->GetParameter(uPar)*randomize);
        FitFoldFinal->SetParLimits(uPar,0,FitFoldOriginal->GetParameter(uPar)*10.);
        if(dlmData_reb2->GetBinCenter(0,uPar)<RangeMinU||dlmData_reb2->GetBinCenter(0,uPar)>RangeMaxU){
          FitFoldFinal->FixParameter(uPar,FitFoldOriginal->GetParameter(uPar));
        }
        //FitFoldFinal->SetParLimits(uPar,FitFoldOriginal->GetParameter(uPar)*0.1,FitFoldOriginal->GetParameter(uPar)*10.);
      }
      FitFoldFinal->SetNpx(256);
      gROOT->cd();
      hData->Fit(FitFoldFinal,"S, N, R, M");
      f_chi2 = FitFoldFinal->GetChisquare();
      if(f_chi2<FitFoldOriginal->GetChisquare()*100) break;
      if(!Silent) printf(" Refitting (%u, chi2=%.1f)\n",uRefit,f_chi2);
    }
    f_chi2 = FitFoldFinal->GetChisquare();
    f_ndf = FitFoldFinal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitUnfoldFinal){delete FitUnfoldFinal; FitUnfoldFinal=NULL;}
    gROOT->cd();

    //if(dlmDimiWorkHorse){delete [] dlmDimiWorkHorse; dlmDimiWorkHorse=NULL;}
    //dlmDimiWorkHorse = new DLM_Histo<float>(*dlmData_reb2);
    FitUnfoldFinal = new TF1("FitUnfoldFinal",this,&DLM_Unfold::DimiFit,RangeMinU,RangeMaxU,dlmData_reb2->GetNbins(),"DLM_Unfold","DimiFit");
    for(unsigned uPar=0; uPar<dlmData_reb2->GetNbins(); uPar++){
      FitUnfoldFinal->SetParameter(uPar,FitFoldFinal->GetParameter(uPar));
      FitUnfoldFinal->SetParLimits(uPar,0,FitFoldFinal->GetParameter(uPar)*10.);
    }
    FitUnfoldFinal->SetNpx(256);
  }
  else if(TYPE==0){
    if(FitFoldOriginal){delete FitFoldOriginal; FitFoldOriginal=NULL;}
    gROOT->cd();
    FitFoldOriginal = new TF1("FitFoldOriginal",this,&DLM_Unfold::SplineFit,RangeMinU,RangeMaxU,3+NumKnots*2,"DLM_Unfold","SplineFit");
    FitFoldOriginal->FixParameter(0,NumKnots);
    //derivative at the firs knot
//FitFoldOriginal->SetParameter(1,0);
    FitFoldOriginal->FixParameter(1,0);
    FitFoldOriginal->SetParameter(2,0);
    //FitFoldOriginal->SetParLimits(2,-1e-3,1e-3);
    for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
        double HistVal;
        if(uKnot==0) {HistVal = 1e-128;}
        else {HistVal = hData->GetBinContent(hData->FindBin(Knot_x[uKnot]));}
        FitFoldOriginal->FixParameter(3+uKnot,Knot_x[uKnot]);
        FitFoldOriginal->SetParameter(3+NumKnots+uKnot,HistVal*1);
        FitFoldOriginal->SetParLimits(3+NumKnots+uKnot,0,HistVal*2);
        if(uKnot==0){FitFoldOriginal->FixParameter(3+NumKnots+uKnot,HistVal);}
    }
    FitFoldOriginal->SetNpx(256);
    //printf("About to fit hData_pL\n");
    //TH1F* hFIT = (TH1F*)hData->Clone("hFIT");
    //TF1* ftest = new TF1("ftest","[0]+[1]*x",RangeMin,RangeMax);
    gROOT->cd();
    hData->Fit(FitFoldOriginal,"Q, S, N, R, M");

    f_chi2 = FitFoldOriginal->GetChisquare();
    f_ndf = FitFoldOriginal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf("Pre-fit (folded):\n");
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}
    gROOT->cd();
    FitFoldFinal = new TF1("FitFoldFinal",this,&DLM_Unfold::SplineFitFolded,RangeMinU,RangeMaxU,3+NumKnots*2,"DLM_Unfold","SplineFitFolded");

    if(!Silent) printf("Pre-fit:\n");
    f_chi2*=100;
    for(unsigned uRefit=0; uRefit<1; uRefit++){
      double randomize=1;
      if(uRefit){
        randomize = rangen.Uniform(0.9,1.1);
      }

      FitFoldFinal->FixParameter(0,FitFoldOriginal->GetParameter(0));
      FitFoldFinal->FixParameter(1,FitFoldOriginal->GetParameter(1));
      FitFoldFinal->FixParameter(2,FitFoldOriginal->GetParameter(2));
      for(unsigned uKnot=0; uKnot<NumKnots; uKnot++){
          FitFoldFinal->FixParameter(3+uKnot,Knot_x[uKnot]);

          //the edge knot and the first 3 bins (so the first 4 knots) for a ME disto have to have values BELOW the folded solution
          //if(uKnot<4){
          //  FitFoldFinal->SetParameter(3+NumKnots+uKnot,FitFoldOriginal->GetParameter(3+NumKnots+uKnot)*0.1);
          //  FitFoldFinal->SetParLimits(3+NumKnots+uKnot,0,FitFoldOriginal->GetParameter(3+NumKnots+uKnot));
          //}
          //else {
FitFoldFinal->SetParameter(3+NumKnots+uKnot,FitFoldOriginal->GetParameter(3+NumKnots+uKnot));
FitFoldFinal->SetParLimits(3+NumKnots+uKnot,0,FitFoldOriginal->GetParameter(3+NumKnots+uKnot)*10);
            if(uKnot==0){FitFoldFinal->FixParameter(3+NumKnots+uKnot,FitFoldOriginal->GetParameter(3+NumKnots+uKnot));}
          //}
      }
      FitFoldFinal->SetNpx(256);

      //for(unsigned uPar=0; uPar<3+NumKnots*2; uPar++){
      //  FitFoldFinal->SetParameter(uPar,FitFoldOriginal->GetParameter(uPar)*randomize);
      //}
      //FitFoldFinal->SetNpx(256);

      gROOT->cd();
      hData->Fit(FitFoldFinal,"Q, S, N, R, M");
      f_chi2 = FitFoldFinal->GetChisquare();
      if(f_chi2<FitFoldOriginal->GetChisquare()*100) break;
      if(!Silent) printf(" Refitting (%u, chi2=%.1f)\n",uRefit,f_chi2);
    }
    f_ndf = FitFoldFinal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitUnfoldFinal){delete FitUnfoldFinal; FitUnfoldFinal=NULL;}
    gROOT->cd();
    FitUnfoldFinal = new TF1("FitUnfoldFinal",this,&DLM_Unfold::SplineFitUnfolded,RangeMinU,RangeMaxU,3+NumKnots*2,"DLM_Unfold","SplineFitUnfolded");
    for(unsigned uPar=0; uPar<3+NumKnots*2; uPar++){
      FitUnfoldFinal->SetParameter(uPar,FitFoldFinal->GetParameter(uPar));
    }
    FitUnfoldFinal->SetNpx(256);
  }
  else{
    if(FitFoldOriginal){delete FitFoldOriginal; FitFoldOriginal=NULL;}
    gROOT->cd();
    FitFoldOriginal = new TF1("FitFoldOriginal","([0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x)*exp(-pow(x/[4],[5]))",RangeMinU,RangeMaxU);
    FitFoldOriginal->SetParameter(0,0);
    FitFoldOriginal->SetParameter(1,1);
    FitFoldOriginal->SetParameter(2,0);
    FitFoldOriginal->SetParameter(3,0);
    //derivative at the firs knot
    FitFoldOriginal->SetParameter(4,800);
    FitFoldOriginal->SetParameter(5,1.5);
    FitFoldOriginal->SetNpx(256);
    gROOT->cd();
    hData->Fit(FitFoldOriginal,"Q, S, N, R, M");

    f_chi2 = FitFoldOriginal->GetChisquare();
    f_ndf = FitFoldOriginal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf("Pre-fit (folded):\n");
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}
    gROOT->cd();
    FitFoldFinal = new TF1("FitFoldFinal","([0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x)*exp(-pow(x/[4],[5]))",RangeMinU,RangeMaxU);

    if(!Silent) printf("Pre-fit:\n");
    f_chi2*=100;
    for(unsigned uRefit=0; uRefit<11; uRefit++){
      for(unsigned uPar=0; uPar<6; uPar++){
        FitFoldFinal->SetParameter(uPar,FitFoldOriginal->GetParameter(uPar));
      }
      FitFoldFinal->SetNpx(256);
      gROOT->cd();
      hData->Fit(FitFoldFinal,"Q, S, N, R, M");
      f_chi2 = FitFoldFinal->GetChisquare();
      if(f_chi2<FitFoldOriginal->GetChisquare()*10) break;
      if(!Silent) printf(" Refitting (%u, chi2=%.1f)\n",uRefit,f_chi2);
    }
    f_ndf = FitFoldFinal->GetNumberFitPoints();
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));

    if(FitUnfoldFinal){delete FitUnfoldFinal; FitUnfoldFinal=NULL;}
    gROOT->cd();
    FitUnfoldFinal = new TF1("FitUnfoldFinal","([0]*x+[1]*x*x+[2]*x*x*x+[3]*x*x*x*x)*exp(-pow(x/[4],[5]))",RangeMinU,RangeMaxU);
    for(unsigned uPar=0; uPar<6; uPar++){
      FitUnfoldFinal->SetParameter(uPar,FitFoldFinal->GetParameter(uPar));
    }
    FitUnfoldFinal->SetNpx(256);
  }

//return NULL;
  for(unsigned uSeed=SEEDmin; uSeed<SEEDmin+NumIter; uSeed++){
    if(!Silent) printf(" uSeed = %u (%u)\n",uSeed,SEEDmin+NumIter);
/*
    DLM_Histo<float> Chi2Ndf_Unfold;// = new TH1F("Chi2Ndf_Unfold","Chi2Ndf_Unfold",128,0,Unacceptable_chi2ndf);
    DLM_Histo<float> CPU_Unfold;// = new TH1F("CPU_Unfold","CPU_Unfold",128,0,4.*double(JobTimeLimit)/double(NumIter)/60.);
    DLM_Histo<float> Chi2Ndf_CPU_Unfold;// = new TH2F("Chi2Ndf_CPU_Unfold","Chi2Ndf_CPU_Unfold",
                                        //128,0,Unacceptable_chi2ndf,
                                        //128,0,4.*double(JobTimeLimit)/double(NumIter)/60.);
    Chi2Ndf_Unfold.SetUp(1);
    Chi2Ndf_Unfold.SetUp(0,128,0,PrecisionWorst);
    Chi2Ndf_Unfold.Initialize();
*/

    bool TerminateIter = false;
    DLM_Timer ComputingTime;
    long long TotalTime = JobTime.Stop()/1000000.;
//printf("TotalTime = %lli",TotalTime);
    double TimeRemaining = double(MaxTime-TotalTime);
    double IterRemaining = double(SEEDmin+NumIter-uSeed);
    double AvgTimePerSeedRemaining = TimeRemaining/IterRemaining;
    double NumIterDone = double(uSeed-SEEDmin);//completed iterations (not counting the current)
    double NumIterRemain = double(SEEDmin+NumIter-uSeed-1);//remaining iteration (not counting the current)
    double AvgIterTime = NumIterDone?IterTime/NumIterDone:0;
    //a bonus that we add to the remaining average time, based on the difference between the
    //remaining avg. and the current average time. If we are bad on time, this bonus becomes negative
    //and serves as a regulator to bring us back on schedule
    double BufferTime = AvgTimePerSeedRemaining-AvgIterTime;
    //the maximum allowed time is: on the first iter, the average
    //afterwards, it represents the remaining time, minus the expected time needed for the remaining iterations,
    //based on the average run time so far. This is in case we do good on time
    double MaxIterTime = NumIterDone?TimeRemaining-AvgIterTime*NumIterRemain:AvgTimePerSeedRemaining;
    //in case we are doing bad on average, we try to give some extra bonus on the current runs,
    //by increasing by up to 33% the allowed run time
    double BonusFactor = Avg_chi2ndf/Precision;
    if(BonusFactor<1) BonusFactor=1;
    if(BonusFactor>1.33) BonusFactor=1.33;
    //if we do bad on time, the upper definition can lead to very very small max time.
    //if too small, we set the max time to its 'minimal' max time, which is the average time remaining per iter
    if(MaxIterTime<AvgTimePerSeedRemaining*BonusFactor) MaxIterTime = AvgTimePerSeedRemaining*BonusFactor;

    //double MaxTime = NumIterDone?TimeRemaining/(NumIterDone+1.):AvgTimePerSeedRemaining;
    double AllowedTime = AvgTimePerSeedRemaining+BufferTime;
    AllowedTime *= BonusFactor;
    if(AllowedTime>MaxIterTime) AllowedTime=MaxIterTime;
    PatienceMax = (long long)(AllowedTime);
    PatienceForPerfection = (long long)(AvgTimePerSeedRemaining);
    ShowTime(PatienceMax,strtime,2,true,5);
    if(!Silent) printf("  Time limit: %s\n",strtime);

    double OriginalValue;
    //double OriginalValueMin;
    //double OriginalValueMax;
    double Error;
    double RandomValue;
    double EffectiveCounts;//assuming E/N = sqrt(N)/N = 1/sqrt(N) = relErr; N = 1/relErr^2
    //double RandomMean;
    double xAxis;

    DLM_Histo<float> dlmData_boot = (*dlmData);
    for(unsigned uBin=0; uBin<dlmData->GetNbins(); uBin++){
      xAxis = dlmData->GetBinCenter(0,uBin);
      if(xAxis<RangeMinU) continue;
      if(xAxis>RangeMaxU) break;
      OriginalValue = dlmData->GetBinContent(uBin);
      Error = dlmData->GetBinError(uBin);
      EffectiveCounts = pow(Error/OriginalValue,-2.);
      if(uSeed){
        //RandomValue = rangen.Gaus(OriginalValue,Error);
        RandomValue = rangen.Poisson(EffectiveCounts)*OriginalValue/EffectiveCounts;
      }
      else RandomValue = dlmData->GetBinContent(uBin);
      dlmData_boot.SetBinContent(uBin,RandomValue);
      dlmData_boot.SetBinError(uBin,Error);
    }

    //TH1F* hData_Unfolded = (TH1F*)hData_boot->Clone("hData_Unfolded");
    DLM_Histo<float> dlmUnfolded(*dlmData);
    DLM_Histo<float> dlmFolded(*dlmData);
    double BestChi2=1e6;
    double BestChi2_Data=1e6;
    double BestChi2_Fit=1e6;
    double BestChi2_BinData=1e6;
    double BestNDF=0;
    DLM_Histo<float> dlmSampleUnfolded(*dlmData);
    DLM_Histo<float> dlmBestUnfolded(*dlmData);
    DLM_Histo<float> dlmBestFolded(*dlmData);

    unsigned StepsWithoutImprovement=0;
    //unsigned NumSuccess=0;
    //unsigned NumAttempts=0;
    //double SuccessRate=1;

    unsigned* NumSuccess = new unsigned[dlmSampleUnfolded.GetNbins()];
    unsigned* NumAttempts = new unsigned[dlmSampleUnfolded.GetNbins()];
    double* SuccessRate = new double[dlmSampleUnfolded.GetNbins()];
    double* Impact = new double[dlmSampleUnfolded.GetNbins()];
    const unsigned HistoryLength = 16;
    std::vector<float>* BinByBinChi2 = new std::vector<float>[dlmSampleUnfolded.GetNbins()];
    double MinSuccessRate=0.25;
    unsigned MinAttempts=2;

    //bootstrap the sampling histogram ones, to avoid a bias in the chi2 determination (where the original solution is used as a benchmark)
    //here we should sample from the unfolded soltion, in fact we use the fit from upstairs
    for(unsigned uBin=0; uBin<dlmSampleUnfolded.GetNbins(); uBin++){
      xAxis = dlmData->GetBinCenter(0,uBin);
      if(xAxis<RangeMinU) continue;
      if(xAxis>RangeMaxU) break;
      OriginalValue = FitUnfoldFinal->Eval(xAxis);//dlmData->GetBinContent(uBin);
      Error = dlmData->GetBinError(uBin+1);
      //RandomValue = rangen.Gaus(OriginalValue,Error);
      //what the hell is this??
      RandomValue=OriginalValue;
      dlmSampleUnfolded.SetBinContent(uBin,RandomValue);
if(RandomValue<0){
printf("uBin=%u; OriginalValue=%.3e; RandomValue=%.3e\n",uBin,OriginalValue,RandomValue);
}

      NumSuccess[uBin]=0;
      NumAttempts[uBin]=0;
      SuccessRate[uBin]=1;
      Impact[uBin]=0;
    }

    //double* ChiBinByBin = new double [];

    for(unsigned uBF=0; uBF<NumBackAndForth; uBF++){
      sprintf(ctext1,"   uBF = %u (%u)",uBF,NumBackAndForth);
      strcpy(ctext2,"");
      strcpy(ctext3,"");
      std::cout << ctext1 << std::flush;
      bool REFRESHED_INFO = false;

      //MinAttempts = 10 + uBF*4;
      //MinSuccessRate = 1./(10 + uBF*16);
      for(unsigned uBin=0; uBin<dlmSampleUnfolded.GetNbins(); uBin++){
        NumSuccess[uBin]=0;
        NumAttempts[uBin]=0;
        SuccessRate[uBin]=1;
        Impact[uBin]=0;
      }

      for(unsigned uBin=0; uBin<2*dlmSampleUnfolded.GetNbins(); uBin++){
        //we iterate twice, from below and from above.
        //The second iteration is for 'polishing',
        //which is performed for fewer iterations, but higher depth
        bool Polishing = uBin>=dlmSampleUnfolded.GetNbins();
        //keep track of time, in case we exceed our limit, terminate asap
        TotalTime = JobTime.Stop()/1000000.;
        if(TotalTime>MaxTime){
          TerminateASAP = true;
          break;
        }
        unsigned WhichBin = (!Polishing)?uBin:2*dlmSampleUnfolded.GetNbins()-uBin-1;
        double HighestImpact = -1000;
        //printf("\nIMPACT:\n");
        for(unsigned ub=0; ub<dlmSampleUnfolded.GetNbins(); ub++){
          xAxis = dlmData->GetBinCenter(0,ub);
          if(xAxis<RangeMinU) continue;
          if(xAxis>RangeMaxU) break;
          if(Impact[ub]>HighestImpact){HighestImpact=Impact[ub];}
        //  printf("Impact[%u]=%.3e\n",ub,Impact[ub]);
        //  usleep(50e3);
        }
        //printf("HighestImpact = %f\n",HighestImpact);
        double ProbToCompute=exp(-(HighestImpact-Impact[WhichBin])/(fabs(HighestImpact)+1e-16));
        //if(ProbToCompute<1) printf("PTC = %.2f%%, HI=%.3e\n",ProbToCompute*100.,HighestImpact);
        if(rangen.Uniform()>ProbToCompute) continue;
        if(!Polishing&&WhichBin>dlmSampleUnfolded.GetNbins()-BinDepth) continue;
        //if(ProbToCompute<1) printf(" -> go\n");
        if(Polishing&&WhichBin<BinDepth-1) continue;
        unsigned MaxStepsWithoutImprovement = (!Polishing)?MaxStepsWithoutImprovement1:MaxStepsWithoutImprovement2;
        //if we are in a tough spot, increase the MaxSteps
        MaxStepsWithoutImprovement *= (uBF/2+1);
        //we monitor how often we get an improvement. If often, we keep on digging, as it will
        //help reduce the chi2 from regions that can be easily improved

        //if(HighestImpact==Impact[WhichBin]) printf("\n HIB = %u\n",WhichBin);

        for(unsigned uED=1*unsigned(Polishing); uED<ErrorDepth+1*unsigned(Polishing); uED++){
          double Scale = pow(0.5,double(uED));
          StepsWithoutImprovement=0;
          //NumSuccess = 0;
          //NumAttempts = 0;
          //SuccessRate = 1;
          unsigned uBin_From =(!Polishing)?WhichBin:WhichBin-BinDepth+1;
          unsigned uBin_To = (!Polishing)?WhichBin+BinDepth:WhichBin+1;
          //while(StepsWithoutImprovement<MaxStepsWithoutImprovement || SuccessRate>MinimumSuccessRate){
          while(SuccessRate[uBin_From]>MinSuccessRate||NumAttempts[uBin_From]<MinAttempts){

            //printf("%u: SR=%.2f/%.2f; NA=%u/%u\n",uBin_From,SuccessRate[uBin_From]*100.,MinSuccessRate*100.,NumAttempts[uBin_From],MinAttempts);
            dlmUnfolded = dlmSampleUnfolded;
            //we iterate over the desired bin range and bootstrap a bit
            for(unsigned uBin2=uBin_From; uBin2<uBin_To; uBin2++){
              OriginalValue = dlmSampleUnfolded.GetBinContent(WhichBin);
              //Error = dlmData_boot.GetBinError(WhichBin);
              //if(NumAttempts[uBin_From]%2==1) Error*=2;
              //else if(NumAttempts[uBin_From]%2==2) Error*=4;
              //EffectiveCounts = pow(Error*Scale/OriginalValue,-2.);
              EffectiveCounts = OriginalValue*YieldScale;
              //if(NumAttempts[uBin_From]%2==1) EffectiveCounts/=2.;
              //if(NumAttempts[uBin_From]%2==2) EffectiveCounts/=4.;

              /*
              //if(StepsWithoutImprovement%2) RandomValue = rangen.Uniform(OriginalValue-0.5*Error*Scale,OriginalValue+0.5*Error*Scale);
              if(NumAttempts[uBin_From]%3==0) RandomValue = rangen.Uniform(OriginalValue-0.5*Error*Scale,OriginalValue+0.5*Error*Scale);
              //half of the time we allow for a larger spread
              //else RandomValue = rangen.Gaus(OriginalValue,Error*Scale);
              else if(NumAttempts[uBin_From]%3==1){
                RandomValue = rangen.Uniform(OriginalValue-3*Error*Scale,OriginalValue+3*Error*Scale);
              }
              else{
                RandomValue = rangen.Poisson(EffectiveCounts);
                RandomValue*=OriginalValue;
                RandomValue/=EffectiveCounts;
              }
              */
              RandomValue = rangen.Poisson(EffectiveCounts);
              RandomValue*=OriginalValue;
              RandomValue/=EffectiveCounts;

              RandomValue = fabs(RandomValue);
              dlmUnfolded.SetBinContent(WhichBin,RandomValue);
if(RandomValue<0){
printf("WhichBin=%u; OriginalValue=%.3e; RandomValue=%.3e\n",WhichBin,OriginalValue,RandomValue);
}
            }
            //Ck_Unfolded->Update(true);
            //CkDec_Unfolded->Update(true);
            Fold(dlmUnfolded,dlmFolded);
            double Chi2=0;
            double Chi2_Data=0;
            double Chi2_Fit=0;
            double Momentum;
            double CkValCorrected;
            double CkValData;
            double Error;
            double Chi2_BinData;
            double Chi2_BinFit;
            double Chi2_BinDataMax=0;
            double Chi2_BinFitMax=0;
            double Ndf=0;

            for(unsigned uBin2=0; uBin2<dlmData_boot.GetNbins(); uBin2++){
              xAxis = dlmData_boot.GetBinCenter(0,uBin2);
              if(xAxis<RangeMinU) continue;
              if(xAxis>RangeMaxU) break;
              CkValCorrected = dlmFolded.GetBinContent(uBin2);
              CkValData = dlmData_boot.GetBinContent(uBin2);
              Error = dlmData_boot.GetBinError(uBin2);
              //folded - current solution vs data
              Chi2_BinData = pow((CkValData-CkValCorrected)/(Error),2.);
              //unfolded - current solution vs fit
              Chi2_BinFit = pow((dlmUnfolded.GetBinContent(uBin2)-FitUnfoldFinal->Eval(xAxis))/(Error),2.);
              Chi2 += Chi2_BinData;
              Chi2 += Chi2_BinFit;
              Chi2_Data += Chi2_BinData;
              Chi2_Fit += Chi2_BinFit;
              if(uBin2>=uBin_From&&uBin2<uBin_To){
                  if(Chi2_BinData>Chi2_BinDataMax) Chi2_BinDataMax=Chi2_BinData;
                  if(Chi2_BinFit>Chi2_BinFitMax) Chi2_BinFitMax=Chi2_BinFit;
              }
              Ndf++;
            }

            //if positive -> good
            double Improvement = BestChi2_Data-Chi2_Data;
            if(BinByBinChi2[WhichBin].size()>=HistoryLength){BinByBinChi2[WhichBin].erase(BinByBinChi2[WhichBin].begin());}
            BinByBinChi2[WhichBin].push_back(Improvement);
            Impact[WhichBin]=0;
            for(unsigned uHist=0; uHist<BinByBinChi2[WhichBin].size(); uHist++){
              Impact[WhichBin]+=BinByBinChi2[WhichBin].at(uHist);
            }
            //the average Improvement over the last few iterations (bin-by-bin)
            Impact[WhichBin]/=double(BinByBinChi2[WhichBin].size());

            //the idea: on the first iteration, make the solution close to the expectation (find the desired global min)
            //          on the next iterations, loosen up the condition to explore local minima
            if(Chi2_Data<BestChi2_Data&&Chi2_BinFitMax<=double(uBF/16+1)*Chi2_BinDataMax){
            //if(Chi2_BinDataMax<BestChi2_BinData){
                sprintf(ctext3," --> %.1f",Chi2_Data);
                if(!Silent&&PcOutput){
                    printf("\r\033[K%s%s%s",ctext1,ctext2,ctext3);
                    std::cout << std::flush;
                }
                if(!REFRESHED_INFO){
                    sprintf(ctext2," χ2 = %.1f",Chi2_Data);
                    if(!Silent&&!PcOutput) std::cout << ctext2 << std::flush;
                    REFRESHED_INFO = true;
                }

                BestChi2=Chi2;
                BestChi2_Data=Chi2_Data;
                BestChi2_BinData=Chi2_BinDataMax;

                BestChi2_Fit=Chi2_Fit;
                BestNDF=Ndf;
                dlmBestFolded = dlmFolded;
                dlmBestUnfolded = dlmUnfolded;
                dlmSampleUnfolded = dlmUnfolded;
                StepsWithoutImprovement=0;
                if(BestChi2_Data/BestNDF<Precision){
                    TerminateIter = true;
                    break;
                }
                long long ExeTime = ComputingTime.Stop()/1000000.;
                if(ExeTime>PatienceMax){
                    TerminateIter = true;
                    break;
                }
                NumSuccess[uBin_From]++;
            }
            else{
                StepsWithoutImprovement++;
            }
            NumAttempts[uBin_From]++;
            SuccessRate[uBin_From] = double(NumSuccess[uBin_From])/double(NumAttempts[uBin_From]);
          }//while
          if(TerminateIter) break;
        }//uED
        if(TerminateIter) break;
      }//uBin

      long long ExeTime = ComputingTime.Stop()/1000000.;
      f_chi2 = BestChi2_Data;
      f_ndf = dlmData->GetBin(0,RangeMaxU)-dlmData->GetBin(0,RangeMinU)+1;
      f_chi2ndf = f_chi2/double(f_ndf);
      f_pval = TMath::Prob(f_chi2,f_ndf);
      f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);

      if(!REFRESHED_INFO){
          sprintf(ctext2," χ2 = %.1f",f_chi2);
          std::cout << ctext2 << std::flush;
          REFRESHED_INFO = true;
      }
      if(!Silent&&PcOutput){printf("\r\033[K%s%s",ctext1,ctext2);}
      if(!Silent)printf(" --> %.1f\n",f_chi2);

      //terminate ASAP
      if(f_chi2ndf<Precision){
          ShowTime(ExeTime,strtime,2,true,5);
          if(!Silent)printf("  Very good solution found at %s\n",strtime);
          break;
      }
      //terminate if we are slow. If we are fast, try to get to a very good solution
      if(f_chi2ndf<Avg_chi2ndf&&double(ExeTime)>AvgTimePerSeedRemaining&&double(ExeTime)>AvgIterTime){
          ShowTime(ExeTime,strtime,2,true,5);
          if(!Silent)printf("  Good solution found at %s\n",strtime);
          break;
      }
      //terminate ASAP if we have exceeded the time limit
      if(ExeTime>PatienceMax){
          ShowTime(ExeTime,strtime,2,true,5);
          if(!Silent)printf("  Over the time limit at %s\n",strtime);
          break;
      }

      TotalTime = JobTime.Stop()/1000000.;
      if(TotalTime>MaxTime){
          printf(" TT = %lli; MT = %f\n",TotalTime,MaxTime);
          ShowTime(TotalTime,strtime,2,true,5);
          if(!Silent)printf("  Time limit exceeded. Terminating at %s\n",strtime);
          TerminateASAP = true;
          break;
      }

    }//uBF

    delete [] NumAttempts;
    delete [] NumSuccess;
    delete [] SuccessRate;
    delete [] Impact;
    delete [] BinByBinChi2;

    long long ExeTime = ComputingTime.Stop()/1000000.;
    Avg_chi2ndf = (Avg_chi2ndf*NumIterDone+f_chi2ndf)/(NumIterDone+1);
    IterTime += ExeTime;
    //Chi2Ndf_Unfold->Fill(f_chi2ndf);
    //CPU_Unfold->Fill(double(ExeTime)/60.);
    //Chi2Ndf_CPU_Unfold->Fill(f_chi2ndf,double(ExeTime)/60.);

    if(!Silent)printf("  Final fit quality:\n");
    if(!Silent)printf("   chi2/ndf = %.2f/%i (%.2f)\n",f_chi2,f_ndf,Avg_chi2ndf*f_ndf);
    //printf("   pval = %.4f\n",f_pval);
    //printf("   nsigma = %.3f\n",f_ns);

    bool AcceptSolution = (f_chi2ndf<PrecisionWorst);
    if(!AcceptSolution){
        if(!Silent)printf(" REJECTED\n\n");
        continue;
    }
    else{
        if(!Silent)printf(" ACCEPTED\n\n");
    }

    if(AcceptSolution){
      unsigned BinCounter=0;
      double Mean=0;
      double Mean2=0;
      for(unsigned uBin=0; uBin<dlmData_boot.GetNbins(); uBin++){
        xAxis = dlmData_boot.GetBinCenter(0,uBin);
        if(xAxis<RangeMinU) continue;
        if(xAxis>RangeMaxU) break;
        BinByBinFolded[BinCounter].push_back(dlmBestFolded.GetBinContent(uBin));
        BinByBinUnfolded[BinCounter].push_back(dlmBestUnfolded.GetBinContent(uBin));

        if(!NumResultEnties){
          Mean = dlmBestUnfolded.GetBinContent(uBin);
//if(uBin<5)
//printf("1. Mean = %f\n",Mean);
          Mean2 = dlmBestUnfolded.GetBinContent(uBin)*dlmBestUnfolded.GetBinContent(uBin);
        }
        else{
          Mean = dlmUnfoldResult->GetBinContent(uBin);//dlmUnfoldResult
//if(uBin<5)
//printf("2. Mean = %f\n",Mean);
          Mean += dlmBestUnfolded.GetBinContent(uBin)/double(NumResultEnties);
//if(uBin<5)
//printf("3. Mean = %f\n",Mean);
          Mean *= double(NumResultEnties)/double(NumResultEnties+1);

          Mean2 = dlmUnfoldResult->GetBinError(uBin);//dlmUnfoldResult
          Mean2 += dlmBestUnfolded.GetBinContent(uBin)/double(NumResultEnties);
          Mean2 *= double(NumResultEnties)/double(NumResultEnties+1);
        }
        dlmUnfoldResult->SetBinContent(uBin,Mean);
//if(uBin<5)
//printf("4. Mean = %f\n\n",Mean);
        dlmUnfoldResult->SetBinError(uBin,Mean);

        BinCounter++;
      }
      NumResultEnties++;
    }
    if(TerminateASAP) break;
  }//uSeed

  for(unsigned uBin=0; uBin<dlmUnfoldResult->GetNbins(); uBin++){
    double xAxis = dlmUnfoldResult->GetBinCenter(0,uBin);
    if(xAxis<RangeMinU) continue;
    if(xAxis>RangeMaxU) break;
    dlmUnfoldResult->SetBinError(uBin,sqrt(dlmUnfoldResult->GetBinError(uBin)-dlmUnfoldResult->GetBinContent(uBin)*dlmUnfoldResult->GetBinContent(uBin)));
  }

  delete [] strtime;
  delete [] ctext1;
  delete [] ctext2;
  delete [] ctext3;
  delete [] Knot_x;
//  delete fit_original;
  return dlmUnfoldResult;
}



DLM_Histo<float>* DLM_Unfold::UnfoldDimi(){
  if(!dlmData||!dlmResponse){
    printf("\033[1;33mWARNING!\033[0m DLM_Unfold misses some input to unfold!\n");
    return NULL;
  }

  DLM_Timer JobTime;
  TRandom3 rangen(SEEDmin+1);

  const unsigned NumBins = dlmData->GetBin(0,RangeMaxU)-dlmData->GetBin(0,RangeMinU)+1;
  if(NumBins<=1){
    printf("\033[1;33mWARNING!\033[0m Smearing a less than 2 bins is truely the next level! Sorry, cannot do it!\n");
    return NULL;
  }

  if(BinByBinFolded){delete [] BinByBinFolded; BinByBinFolded=NULL;}
  BinByBinFolded = new std::vector<float> [NumBins];
  if(BinByBinUnfolded){delete [] BinByBinUnfolded; BinByBinUnfolded=NULL;}
  BinByBinUnfolded = new std::vector<float> [NumBins];
  DLM_Histo<float>* dlmUnfoldResult = new DLM_Histo<float>(*dlmData);
  dlmUnfoldResult->SetBinContentAll(0);
  dlmUnfoldResult->SetBinErrorAll(0);
  unsigned NumResultEnties = 0;

  char* strtime = new char [128];
  char* ctext1 = new char [128];

  if(!Silent) printf("\nUnfolding:\n\n");
  double IterTime = 0;

  double f_chi2;
  int f_ndf;
  double f_chi2ndf;
  double f_pval;
  double f_ns;

  if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}

  for(unsigned uReb=1; uReb>=1; uReb/=2){
//printf("uReb %u\n",uReb);
    if(dlmData_reb2){delete dlmData_reb2; dlmData_reb2=NULL;}
    dlmData_reb2 = new DLM_Histo<float>(*dlmData);
    dlmData_reb2->Rebin(&uReb);
    if(dlmDimiWorkHorse){delete dlmDimiWorkHorse; dlmDimiWorkHorse=NULL;}
    dlmDimiWorkHorse = new DLM_Histo<float>(*dlmData_reb2);
    gROOT->cd();
    unsigned NumPars = dlmData_reb2->GetBin(0,RangeMaxU)-dlmData_reb2->GetBin(0,RangeMinU);
    TF1* FitTempFold = new TF1(TString::Format("FitFoldFinal_%u",uReb),this,&DLM_Unfold::DimiFitFolded,RangeMinU,RangeMaxU,dlmData_reb2->GetNbins(),"DLM_Unfold","DimiFitFolded");
    double xVal;
    for(unsigned uPar=0; uPar<dlmData_reb2->GetNbins(); uPar++){
      xVal = dlmData_reb2->GetBinCenter(0,uPar);
      //FitTempFold->SetParameter(uPar,FitFoldFinal?FitFoldFinal->Eval(xVal):dlmData_reb2->GetBinContent(uPar));
      //FitTempFold->SetParLimits(uPar,0,1.5*FitTempFold->GetParameter(uPar));
      FitTempFold->SetParameter(uPar,dlmData_reb2->GetBinContent(uPar));
      double UpLimit;
      if(FitFoldFinal) UpLimit=FitFoldFinal->Eval(xVal)*1.5;
      else UpLimit=dlmData_reb2->GetBinContent(uPar)*1.5;
      FitTempFold->SetParLimits(uPar,0,UpLimit);
      //printf(" p%u: %f --> %f\n",uPar,FitTempFold->GetParameter(uPar),1.5*FitTempFold->GetParameter(uPar));
    }
    FitTempFold->SetNpx(256);
    if(FitFoldFinal){delete FitFoldFinal; FitFoldFinal=NULL;}
    FitFoldFinal = FitTempFold;

    gROOT->cd();
    //printf("hData 0 = %.3f +/- %.3f\n",hData->GetBinContent(1),hData->GetBinError(1));
    //printf("hData 1 = %.3f +/- %.3f\n",hData->GetBinContent(2),hData->GetBinError(2));
    hData->Fit(FitFoldFinal,"S, N, R, M");
//printf("fit done");
    f_chi2 = FitFoldFinal->GetChisquare();
    f_ndf = FitFoldFinal->GetNumberFitPoints();
//printf("f_chi2 %f; f_ndf %i\n",f_chi2,f_ndf);
    f_chi2ndf = f_chi2/double(f_ndf);
    f_pval = TMath::Prob(f_chi2,f_ndf);
    f_ns = sqrt(2)*TMath::ErfcInverse(f_pval+1e-32);
    if(!Silent) printf("Fit (rebin %u):\n",uReb);
    if(!Silent) printf(" chi2/ndf = %.1f/%i = %.2f\n",f_chi2,f_ndf,f_chi2/double(f_ndf));
//printf("gROOT\n");
    gROOT->cd();
    if(FitUnfoldFinal){delete FitUnfoldFinal; FitUnfoldFinal=NULL;}
//printf("new\n");
    FitUnfoldFinal = new TF1(TString::Format("FitUnfoldFinal_%u",uReb),this,&DLM_Unfold::DimiFit,RangeMinU,RangeMaxU,dlmData_reb2->GetNbins(),"DLM_Unfold","DimiFit");
//printf("uPar\n");
    for(unsigned uPar=0; uPar<dlmData_reb2->GetNbins(); uPar++){
      xVal = dlmData_reb2->GetBinCenter(0,uPar);
      FitUnfoldFinal->FixParameter(uPar,FitFoldFinal->GetParameter(uPar));
    }
    FitUnfoldFinal->SetNpx(256);

  }//uReb
printf("outside\n");
  for(unsigned uBin=0; uBin<dlmUnfoldResult->GetNbins(); uBin++){
    dlmUnfoldResult->SetBinContent(uBin,FitUnfoldFinal->Eval(dlmUnfoldResult->GetBinCenter(0,uBin)));
  }
printf("deleting\n");
  delete [] strtime;
  delete [] ctext1;
  return dlmUnfoldResult;
printf("ending\n");
}



DLM_Histo<float>* DLM_Unfold::Fold(){
  if(!dlmData) return NULL;
  DLM_Histo<float>* dlmFolded = new DLM_Histo<float>(*dlmData);
  Fold(*dlmData,*dlmFolded);
  return dlmFolded;
}

//N.B. I think I rely on having the same binning between response and data, not to mention that
//I pick only a single value from the response matrix, not the ingegral... too bad
void DLM_Unfold::Fold(const DLM_Histo<float>& DATA, DLM_Histo<float>& RESULT){
  if(!dlmResponse){
    RESULT=DATA;
    return;
  }
  /*
  if(!dlmResponseReb){
    double* BinRange = DATA.GetBinRange(0);
    dlmResponseReb = new DLM_Histo<float>();
    dlmResponseReb->SetUp(2);
    dlmResponseReb->SetUp(0,DATA.GetNbins(),BinRange);
    dlmResponseReb->SetUp(1,DATA.GetNbins(),BinRange);
    dlmResponseReb->Initialize();
    //printf("About to rebin\n");
    dlmResponse->Rebin(*dlmResponseReb);
    //printf("About to sparse\n");
    SparseTheResponse();
    delete [] BinRange;
  }
  */
  const unsigned NumBins = DATA.GetNbins();
  unsigned uBin[2];
  double Momentum[2];
  //DLM_Histo<float>* dlmFolded = new DLM_Histo<float>(DATA);
  for(uBin[1]=0; uBin[1]<NumBins; uBin[1]++){
    //printf(" [%u ... %u ... %u]\n",0,uBin[1],NumBins);
    Momentum[1] = DATA.GetBinCenter(0,uBin[1]);
    if(Momentum[1]>RangeMaxF||Momentum[1]<RangeMinF) continue;
    double Integral=0;
    double Normalization=0;
    //[low/up][sDim]
    double BinRange[2][2];
    //double RespMatrInt;
    //double RespMatrNorm;
    BinRange[0][1] = DATA.GetBinLowEdge(0,uBin[1]);
    BinRange[1][1] = DATA.GetBinUpEdge(0,uBin[1]);
    for(uBin[0]=SparseFirst[uBin[1]]; uBin[0]<=SparseLast[uBin[1]]; uBin[0]++){
    //for(uBin[0]=0; uBin[0]<NumBins; uBin[0]++){
      //printf(" X: [%u ... %u ... %u]\n",SparseFirst[uBin[1]],uBin[0],SparseLast[uBin[1]]);
      Momentum[0] = DATA.GetBinCenter(0,uBin[0]);
      if(Momentum[0]>RangeMaxF||Momentum[0]<RangeMinF) continue;
      BinRange[0][0] = DATA.GetBinLowEdge(0,uBin[0]);
      BinRange[1][0] = DATA.GetBinUpEdge(0,uBin[0]);
      if(!DATA.GetBinContent(uBin)) continue;
      //RespMatrInt = dlmResponse->Integral(BinRange[0],BinRange[1],true);
      //printf(" RM = %f vs %f\n",dlmResponse->Eval(Momentum),RespMatrInt);
      Integral += DATA.GetBinContent(uBin)*dlmResponse->Eval(Momentum);
      //Integral += DATA.GetBinContent(uBin)*RespMatrInt;
      Normalization += dlmResponse->Eval(Momentum);
    }
    Integral /= Normalization;
    RESULT.SetBinContent(uBin[1],Integral);
  }
  RESULT.SetBinErrorAll(0);
}
