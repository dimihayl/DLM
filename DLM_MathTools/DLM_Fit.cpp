#include "DLM_Fit.h"
#include <iostream>
#include "omp.h"
#include <unistd.h>
#include <thread>

#include "DLM_Histo.h"
#include "DLM_Random.h"

DLM_Fit::DLM_Fit():MaxThreads(std::thread::hardware_concurrency()?std::thread::hardware_concurrency():1){

  FitFnct = NULL;
  NumPars = 0;
  NumData = 0;
  chi2 = 0;
  NumDataPts = 0;
  DataLow = NULL;
  DataUp = NULL;
  Data = new std::vector<DLM_Histo<float>*> ();
  Model = new std::vector<DLM_Histo<float>*> ();
  Solution = new std::vector<DLM_FitSolution> ();
  BestSols = new std::vector<DLM_FitSolution> ();
  NumBestSols = 0;
  NumWildCards = 8;
  NumThreads = MaxThreads;
  RanGen = new DLM_Random* [MaxThreads];
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    RanGen[uTh] = new DLM_Random(uTh+1);
  }
}


DLM_Fit::~DLM_Fit(){
  if(DataLow){delete[]DataLow;DataLow=NULL;}
  if(DataUp){delete[]DataUp;DataUp=NULL;}
  if(Model){
    for(DLM_Histo<float>* mdl : *Model){
      if(mdl){delete mdl; mdl=NULL;}
    }
    delete Model; Model=NULL;
  }
  if(Data){
    delete Data; Data=NULL;
  }
  if(Solution){
    delete Solution; Solution=NULL;
  }
  if(BestSols){
    delete BestSols; BestSols=NULL;
  }
  if(RanGen){
    for(unsigned uTh=0; uTh<MaxThreads; uTh++){
      if(RanGen[uTh]){
        delete RanGen[uTh];
        RanGen[uTh] = NULL;
      }
    }
    delete [] RanGen;
    RanGen = NULL;
  }
}



bool DLM_Fit::SetUp(const unsigned& numdata, const unsigned& numpars){

  if(numdata==0){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetUp) At least 1 data set is needed!\n");
    return false;
  }
  if(Data->size()!=numdata){
    Data->resize(numdata);
    for(DLM_Histo<float>* mdl : *Model){
      if(mdl){delete mdl; mdl=NULL;}
    }
    Model->resize(numdata);
    if(DataLow){delete[]DataLow;DataLow=NULL;}
    if(DataUp){delete[]DataUp;DataUp=NULL;}
    DataLow = new std::vector<float> [numdata];
    DataUp = new std::vector<float> [numdata];
  }
  for(DLM_Histo<float>* mdl : *Model){
    mdl=NULL;
  }
//??? do some check on numpars
  NumData = numdata;
  NumPars = numpars;

  return true;

}

bool DLM_Fit::SetData(const unsigned& WhichSet, DLM_Histo<float>& data){

  if(WhichSet>=NumData){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetData) Only %u are set up!\n",NumData);
    return false;
  }
  Data->at(WhichSet) = &data;
  if(Model->at(WhichSet)){delete Model->at(WhichSet); Model->at(WhichSet)=NULL;}
  Model->at(WhichSet) = new DLM_Histo<float> (data);

  DataLow[WhichSet].clear();
  for(unsigned uDim=0; uDim<data.GetDim(); uDim++){
    DataLow[WhichSet].push_back(data.GetLowEdge(uDim));
  }

  DataUp[WhichSet].clear();
  for(unsigned uDim=0; uDim<data.GetDim(); uDim++){
    DataUp[WhichSet].push_back(data.GetUpEdge(uDim));
  }

  return true;

}

void DLM_Fit::SetFitFnct(void (*fitfnc)(const std::vector<float>&, std::vector<DLM_Histo<float>*>&)){
  if(FitFnct==fitfnc) return;
  FitFnct = fitfnc;
}

std::vector<DLM_Histo<float>*> DLM_Fit::Eval(const std::vector<float>& pars){

  if(!FitFnct){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::Eval) The fit function is not set up!\n");
    return *Model;
  }
  //the values of the Model should be set here
  FitFnct(pars,*Model);
  //evaluate the chi2
  chi2 = 0;
  NumDataPts = 0;
  for(unsigned uData=0; uData<Data->size(); uData++){
    if(Data->at(uData)==NULL){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::Eval) The data is not set up!\n");
      return *Model;
    }
    if(Model->at(uData)==NULL){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::Eval) The model is not set up!\n");
      return *Model;
    }
    double* axisVal = new double [Data->at(uData)->GetDim()];
    unsigned* binID = new unsigned [Data->at(uData)->GetDim()];
    for(unsigned uBin=0; uBin<Model->at(uData)->GetNbins(); uBin++){
      Model->at(uData)->GetBinCoordinates(uBin,binID);
      bool GoodBin = true;
      for(unsigned uDim=0; uDim<Model->at(uData)->GetDim(); uDim++){
        if( Model->at(uData)->GetBinCenter(uDim,binID[uDim])<DataLow[uData].at(uDim)||
            Model->at(uData)->GetBinCenter(uDim,binID[uDim])>DataUp[uData].at(uDim)){
          GoodBin = false;
          break;
        }
      }
      if(GoodBin){
        float DataVal, DataErr, ModVal, ModErr;
        DataVal = Data->at(uData)->GetBinContent(uBin);
        DataErr = Data->at(uData)->GetBinError(uBin);
        ModVal = Model->at(uData)->GetBinContent(uBin);
        ModErr = Model->at(uData)->GetBinError(uBin);
        //printf("DV %.3f DE %.3f MV %.3f ME %.3f\n",DataVal, DataErr, ModVal, ModErr);
        chi2 += pow((DataVal-ModVal),2.)/(DataErr*DataErr+ModErr*ModErr);
        NumDataPts++;
      }
    }//uBin
    delete [] axisVal;
  }
  return *Model;

}

float DLM_Fit::Chi2() const{
  return chi2;
}
unsigned DLM_Fit::Npts() const{
  return NumDataPts;
}

void DLM_Fit::SetNumBestSols(const unsigned& num){
  if(num<1){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetNumBestSols) NumBestSols has to be at least 1\n");
      return;
  }
  NumBestSols = num;
  BestSols->clear();
  if(NumWildCards>num){
    NumWildCards = num;
  }
}

void DLM_Fit::SetNumWildCards(const unsigned& num){
  if(num<1||num>BestSols->size()){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetNumWildCards) NumWildCards has to be in [1, NumBestSols]\n");
      return;
  }
  NumWildCards = num;
}

void DLM_Fit::WalkAround(){
  OrderSolutions();

  //add the new solutions
  //order the solutions
}
void DLM_Fit::OrderSolutions(){

}

void DLM_Fit::SetNumThreads(const unsigned& num){
  if(num>MaxThreads||!num){
    NumThreads = MaxThreads;
  }
  else{
    NumThreads = num;
  }
}
void DLM_Fit::SetSeed(const unsigned& thread, const unsigned& seed){
  if(thread>=MaxThreads) return;
  RanGen[thread]->SetSeed(seed);
}
