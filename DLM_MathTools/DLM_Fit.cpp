#include "DLM_Fit.h"
#include <iostream>
#include "omp.h"
#include <unistd.h>
#include <thread>
#include <algorithm>

#include "DLM_Histo.h"
#include "DLM_Random.h"
#include "DLM_CppTools.h"
#include "DLM_MathFunctions.h"
#include "DLM_Sort.h"

DLM_FitSolution::DLM_FitSolution(){
  Par = new std::vector<float> ();
  Chi2=-1; ID=0;
}
DLM_FitSolution::DLM_FitSolution(const DLM_FitSolution& other):DLM_FitSolution(){
  operator=(other);
}
DLM_FitSolution::~DLM_FitSolution(){
  delete Par;
}

DLM_FitSolution DLM_FitSolution::operator-(const DLM_FitSolution& other) const{
  DLM_FitSolution result;
  if(Par->size()!=other.Par->size()) return result;
  for(unsigned uPar=0; uPar<Par->size(); uPar++){
    result.Par->push_back(Par->at(uPar)-other.Par->at(uPar));
  }
  return result;
}
DLM_FitSolution DLM_FitSolution::operator+(const DLM_FitSolution& other) const{
  DLM_FitSolution result;
  if(Par->size()!=other.Par->size()) return result;
  for(unsigned uPar=0; uPar<Par->size(); uPar++){
    result.Par->push_back(Par->at(uPar)+other.Par->at(uPar));
  }
  return result;
}
void DLM_FitSolution::Add(const DLM_FitSolution& other, const float& scale){
  if(Par->size()!=other.Par->size()) return;
  for(unsigned uPar=0; uPar<Par->size(); uPar++){
    Par->at(uPar) = Par->at(uPar)+scale*other.Par->at(uPar);
  }
}
void DLM_FitSolution::operator+=(const DLM_FitSolution& other){
  Add(other,1);
}
void DLM_FitSolution::operator-=(const DLM_FitSolution& other){
  Add(other,-1);
}
void DLM_FitSolution::operator*=(const float& value){
  for(unsigned uPar=0; uPar<Par->size(); uPar++){
    Par->at(uPar) = Par->at(uPar)*value;
  }
}
DLM_FitSolution& DLM_FitSolution::operator=(const DLM_FitSolution& other) {
  Par->clear();
  for(unsigned uPar=0; uPar<other.Par->size(); uPar++){
    Par->push_back(other.Par->at(uPar));
  }
  Chi2=other.Chi2;
  ID=other.ID;

  return *this;
}

DLM_Fit::DLM_Fit(const unsigned num_threads):MaxThreads(std::thread::hardware_concurrency()?num_threads?num_threads<=std::thread::hardware_concurrency()?num_threads:std::thread::hardware_concurrency():std::thread::hardware_concurrency():1){
  InterProcess = false;
  FitFnct = NULL;
  NumPars = 0;
  NumData = 0;
  chi2 = new float [MaxThreads];
  NumDataPts = new unsigned [MaxThreads];
  DataLow = NULL;
  DataUp = NULL;
  Data = new std::vector<DLM_Histo<float>*> ();
  //Model = new std::vector<DLM_Histo<float>*> ();
  Model = new std::vector<DLM_Histo<float>*>* [MaxThreads];
  Solution = new std::vector<DLM_FitSolution> ();
  //BestSols = new std::vector<DLM_FitSolution> ();
  Par = new std::vector<float> ();
  ParL = new std::vector<float> ();
  ParU = new std::vector<float> ();
  ParOL = new std::vector<float> ();
  ParOU = new std::vector<float> ();
  NumBestSols = 32;
  NumWildCards = 8;
  NumThreads = MaxThreads;
  RanGen = new DLM_Random* [MaxThreads];
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    RanGen[uTh] = new DLM_Random(uTh+1);
    Model[uTh] = new std::vector<DLM_Histo<float>*> ();
    chi2[uTh] = -1;
    NumDataPts[uTh] = 0;
  }
}


DLM_Fit::~DLM_Fit(){
  if(chi2){delete[]chi2;chi2=NULL;}
  if(NumDataPts){delete[]NumDataPts;NumDataPts=NULL;}
  if(DataLow){delete[]DataLow;DataLow=NULL;}
  if(DataUp){delete[]DataUp;DataUp=NULL;}
  if(Par){delete Par;Par=NULL;}
  if(ParL){delete ParL;ParL=NULL;}
  if(ParU){delete ParU;ParU=NULL;}
  if(ParOL){delete ParOL;ParOL=NULL;}
  if(ParOU){delete ParOU;ParOU=NULL;}
  //if(Model){
  //  for(DLM_Histo<float>* mdl : *Model){
  //    if(mdl){delete mdl; mdl=NULL;}
  //  }
  //  delete Model; Model=NULL;
  //}
  if(Data){
    delete Data; Data=NULL;
  }
  if(Solution){
    delete Solution; Solution=NULL;
  }
  //if(BestSols){
  //  delete BestSols; BestSols=NULL;
  //}
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
  if(Model){
    for(unsigned uTh=0; uTh<MaxThreads; uTh++){
      if(Model[uTh]){
        delete Model[uTh];
        Model[uTh] = NULL;
      }
    }
    delete [] Model;
    Model = NULL;
  }
}

void DLM_Fit::SetInterProcess(const bool& yesno){
  InterProcess = yesno;
}

void DLM_Fit::SetParLimits(const unsigned& WhichPar, const float& min, const float& max){
  if(WhichPar>=NumPars) return;
  if(min>max){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetParLimits) min>max, check your code!\n");
    return;
  }
  ParL->at(WhichPar) = min;
  ParU->at(WhichPar) = max;
  if(min==max){
    Par->at(WhichPar) = min;
  }
  else if(Par->at(WhichPar)<min || Par->at(WhichPar)>max){
    Par->at(WhichPar) = (min+max)*0.5;
  }
}
void DLM_Fit::SetParOrdMag(const unsigned& WhichPar, const int& min, const int& max){
  if(WhichPar>=NumPars) return;
  if(min>max){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetParOrdMag) min>max, check your code!\n");
    return;
  }
  if(max<-38 || min>38){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetParOrdMag) The min/max allowed order of magnitude is +/-38\n");
    return;
  }

  if(min<-38){
    ParOL->at(WhichPar) = -127;
  }
  else{
    double min_val = pow(10.,min);
    int min_exp;
    frexp(min_val,&min_exp);
    ParOL->at(WhichPar) = min_exp<-127?-127:min_exp;
    //printf("min_exp = %i\n",min_exp);
  }

  if(max>38){
    ParOU->at(WhichPar) = 127;
  }
  else{
    double max_val = pow(10.,max+1);
    int max_exp;
    frexp(max_val,&max_exp);
    ParOU->at(WhichPar) = max_exp>127?127:max_exp;
    //printf("max_exp = %i\n",max_exp);
  }

}
void DLM_Fit::SetParMinMag(const unsigned& WhichPar, const int& min){
  SetParOrdMag(WhichPar, min, 127);
}
void DLM_Fit::SetParMaxMag(const unsigned& WhichPar, const int& max){
  SetParOrdMag(WhichPar, -127, max);
}


bool DLM_Fit::SetUp(const unsigned& numdata, const unsigned& numpars){

  if(numdata==0){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetUp) At least 1 data set is needed!\n");
    return false;
  }
  if(Data->size()!=numdata){
    Data->resize(numdata);
    for(unsigned uTh=0; uTh<MaxThreads; uTh++){
      for(DLM_Histo<float>* mdl : *Model[uTh]){
        if(mdl){delete mdl; mdl=NULL;}
      }
      Model[uTh]->resize(numdata);
    }

    if(DataLow){delete[]DataLow;DataLow=NULL;}
    if(DataUp){delete[]DataUp;DataUp=NULL;}
    DataLow = new std::vector<float> [numdata];
    DataUp = new std::vector<float> [numdata];
  }

  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    for(DLM_Histo<float>* mdl : *Model[uTh]){
      mdl=NULL;
    }
  }

//??? do some check on numpars
  NumData = numdata;
  NumPars = numpars;
  Par->resize(numpars);
  ParL->resize(numpars);
  ParU->resize(numpars);
  ParOL->resize(numpars);
  ParOU->resize(numpars);
  for(float& par : *Par){
    par = 0;
  }
  for(float& lim : *ParL){
    lim = -pow(2.,127.);
  }
  for(float& lim : *ParU){
    lim = pow(2.,127.);
  }
  for(float& lim : *ParOL){
    lim = -127;
  }
  for(float& lim : *ParOU){
    lim = 127;
  }

  return true;
}

bool DLM_Fit::SetData(const unsigned& WhichSet, DLM_Histo<float>& data){

  if(WhichSet>=NumData){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetData) Only %u are set up!\n",NumData);
    return false;
  }
  Data->at(WhichSet) = &data;
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    if(Model[uTh]->at(WhichSet)){delete Model[uTh]->at(WhichSet); Model[uTh]->at(WhichSet)=NULL;}
    Model[uTh]->at(WhichSet) = new DLM_Histo<float> (data);
  }


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

bool DLM_Fit::ProbePars(const std::vector<float>* pars, const unsigned& ThId){
  if(!pars){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::ProbePars) The paremeters are not set!\n");
    return false;
  }
  if(!FitFnct){
    printf("\033[1;31mERROR:\033[0m (DLM_Fit::ProbePars) The fit function is not set up!\n");
    return false;
  }
  //the values of the Model should be set here
//printf("FitFnct...\n");
//printf("size = %u\n",pars->size());
  FitFnct(*pars,*Model[ThId]);
//printf("Okay...\n");
  //evaluate the chi2
  chi2[ThId] = 0;
  NumDataPts[ThId] = 0;
  for(unsigned uData=0; uData<Data->size(); uData++){
    if(Data->at(uData)==NULL){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::Eval) The data is not set up!\n");
      return false;
    }
    if(Model[ThId]->at(uData)==NULL){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::Eval) The model is not set up!\n");
      return false;
    }
//printf("true...\n");
    double* axisVal = new double [Data->at(uData)->GetDim()];
    unsigned* binID = new unsigned [Data->at(uData)->GetDim()];
    for(unsigned uBin=0; uBin<Model[ThId]->at(uData)->GetNbins(); uBin++){
      Model[ThId]->at(uData)->GetBinCoordinates(uBin,binID);
      bool GoodBin = true;
      for(unsigned uDim=0; uDim<Model[ThId]->at(uData)->GetDim(); uDim++){
        //a bin outside of the fit range is bad
        if( Model[ThId]->at(uData)->GetBinCenter(uDim,binID[uDim])<DataLow[uData].at(uDim)||
            Model[ThId]->at(uData)->GetBinCenter(uDim,binID[uDim])>DataUp[uData].at(uDim)){
          GoodBin = false;
          break;
        }
        //a bin that will result in zero error is bad (inf chi2), e.g. empty bins with model uncertainty of 0
        if(Data->at(uData)->GetBinError(uBin)==0 && Model[ThId]->at(uData)->GetBinError(uBin)==0){
          GoodBin = false;
          break;
        }
      }
      if(GoodBin){
        float DataVal, DataErr, ModVal, ModErr;
        DataVal = Data->at(uData)->GetBinContent(uBin);
        DataErr = Data->at(uData)->GetBinError(uBin);
        ModVal = Model[ThId]->at(uData)->GetBinContent(uBin);
        ModErr = Model[ThId]->at(uData)->GetBinError(uBin);
        //printf("DV %.3f DE %.3f MV %.3f ME %.3f\n",DataVal, DataErr, ModVal, ModErr);
        chi2[ThId] += pow((DataVal-ModVal),2.)/(DataErr*DataErr+ModErr*ModErr);
        //chi2[ThId] += pow((DataVal-ModVal),2.)/fabs(ModVal);
        NumDataPts[ThId]++;
      }
    }//uBin
    delete [] axisVal;
  }
  //printf("probe DONE\n");
  return true;
}

std::vector<DLM_Histo<float>*> DLM_Fit::Eval(const std::vector<float>& pars){
  unsigned ThId = omp_get_thread_num();
  if(!ProbePars(&pars,ThId)){
    std::vector<DLM_Histo<float>*> dummy;
    return dummy;
  }
  else{
    ThIdSol = ThId;
    //printf("ThIdSol = %u\n",ThIdSol);
    return *Model[ThId];
  }
/*
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
*/

}

float DLM_Fit::Chi2() const{
  return chi2[ThIdSol];
}
unsigned DLM_Fit::Npts() const{
  return NumDataPts[ThIdSol];
}

void DLM_Fit::SetNumBestSols(const unsigned& num){
  if(num<2){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetNumBestSols) NumBestSols has to be at least 2\n");
      return;
  }
  NumBestSols = num;
  //BestSols->clear();
  if(NumWildCards>num){
    NumWildCards = num;
  }
}

void DLM_Fit::SetNumWildCards(const unsigned& num){
  if(num<1||num>NumBestSols){
      printf("\033[1;31mERROR:\033[0m (DLM_Fit::SetNumWildCards) NumWildCards has to be in [1, NumBestSols]\n");
      return;
  }
  NumWildCards = num;
}

//we will walk along "lines" in the parameter space, connecting two solutions (0 and i).
//the new "solutions" are not evaluated yet, recongnized by their Chi2 = -1
//returns true if we converged
bool DLM_Fit::PrepareForWalk(){
  unsigned ThId = omp_get_thread_num();

  //completely random
  if(Solution->size()<NumBestSols){
    if(Solution->size()){
      printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) Solution->size() shows a bug, contact the developers!\n");
      return false;
    }

    float Rnd;
    float RndMnt;
    int RndExp;
    int RndSign;

    float* MaxAbsVal_Pos = new float [NumPars];
    float* MinAbsVal_Pos = new float [NumPars];
    float* MaxAbsVal_Neg = new float [NumPars];
    float* MinAbsVal_Neg = new float [NumPars];

    int* MinExp_Pos = new int [NumPars];
    int* MaxExp_Pos = new int [NumPars];
    int* MinExp_Neg = new int [NumPars];
    int* MaxExp_Neg = new int [NumPars];

    //mantissa [0.5,1.0]
    //the min/max value is evaluated ASSUMING we have the minimum/maximum order in place
    //non-trivial: if we generate a number that has the max exp, than the mantissa should be 0.5 - MaxManVal
    //             if we generate a number that has the min exp, than the mantissa should be MinManVal - 1.0
    float* MaxManVal_Pos = new float [NumPars];
    float* MinManVal_Pos = new float [NumPars];
    float* MaxManVal_Neg = new float [NumPars];
    float* MinManVal_Neg = new float [NumPars];

    for(unsigned uPar=0; uPar<NumPars; uPar++){
      if(ParIsFixed(uPar)) continue;

      //set the fabs of maximum allowed non-negative number (neg value => no positive values)
      if(ParU->at(uPar)>=0){
        MaxAbsVal_Pos[uPar] = fabs(ParU->at(uPar));
      }
      else{
        MaxAbsVal_Pos[uPar] = -1;
      }

      //set the fabs of minimum allowed non-negative number
      if(ParL->at(uPar)>=0){
        MinAbsVal_Pos[uPar] = fabs(ParL->at(uPar));
      }
      else{
        MinAbsVal_Pos[uPar] = 0;
      }

      //set the fabs of maximum allowed non-positive number (neg value => no negative values)
      if(ParL->at(uPar)<=0){
        MaxAbsVal_Neg[uPar] = fabs(ParL->at(uPar));
      }
      else{
        MaxAbsVal_Neg[uPar] = -1;
      }

      //set the fabs of minimum allowed non-positive number
      if(ParU->at(uPar)<=0){
        MinAbsVal_Neg[uPar] = fabs(ParU->at(uPar));
      }
      else{
        MinAbsVal_Neg[uPar] = 0;
      }

      MinManVal_Pos[uPar] = frexp(MinAbsVal_Pos[uPar],&MinExp_Pos[uPar]);
      MaxManVal_Pos[uPar] = frexp(MaxAbsVal_Pos[uPar],&MaxExp_Pos[uPar]);
      MinManVal_Neg[uPar] = frexp(MinAbsVal_Neg[uPar],&MinExp_Neg[uPar]);
      MaxManVal_Neg[uPar] = frexp(MaxAbsVal_Neg[uPar],&MaxExp_Neg[uPar]);

      //limiting value, better reduce to mantissa 1 and an exponent lower
      if(MaxManVal_Pos[uPar]==0.5){
        MaxManVal_Pos[uPar]=1.0;
        MaxExp_Pos[uPar]--;
      }
      if(MaxManVal_Neg[uPar]==0.5){
        MaxManVal_Neg[uPar]=1.0;
        MaxExp_Neg[uPar]--;
      }

      //if zero is containted, the mantissa has to be set to min value of 0.5, max of 1.0
      //while the exponent is set to minimum value (-127)
      if(MinManVal_Pos[uPar]==0){
        MinManVal_Pos[uPar] = 0.5;
        MinExp_Pos[uPar] = -127;
      }
      if(MinManVal_Neg[uPar]==0){
        MinManVal_Neg[uPar] = 0.5;
        MinExp_Neg[uPar] = -127;
      }


      //this is the case in which the upper limit of the order of magnitude is inbetween the parameter limits
      //in case a signle ord of mag is studied, do nothing, otherwise set the max ord mag to the corresponding value
      //and make sure the manitssa is reset to reach the maximum value (which will still be within the set limits)
      //N.B. in any other case nothing is done. 1) if the max order of mag is larger than the limit, it is obvious why.
      //2) if the max ord mag is lower than the actual parameter limits, is impossible. Thus if set so -> ignored (no warnings presently)
      if(ParOU->at(uPar)>=MinExp_Pos[uPar] && ParOU->at(uPar)<MaxExp_Pos[uPar]){
          MaxManVal_Pos[uPar]=1.0;
          MaxExp_Pos[uPar] = ParOU->at(uPar);
      }
      //the procedure is similar for the min ord mag, but inverted logic
      if(ParOL->at(uPar)<=MaxExp_Pos[uPar] && ParOL->at(uPar)>MinExp_Pos[uPar]){
          MinManVal_Pos[uPar]=0.5;
          MinExp_Pos[uPar] = ParOL->at(uPar);
      }

      //same is repeated for the neg values
      if(ParOU->at(uPar)>=MinExp_Neg[uPar] && ParOU->at(uPar)<MaxExp_Neg[uPar]){
          MaxManVal_Neg[uPar]=1.0;
          MaxExp_Neg[uPar] = ParOU->at(uPar);
      }
      //the procedure is similar for the min ord mag, but inverted logic
      if(ParOL->at(uPar)<=MaxExp_Neg[uPar] && ParOL->at(uPar)>MinExp_Neg[uPar]){
          MinManVal_Neg[uPar]=0.5;
          MinExp_Neg[uPar] = ParOL->at(uPar);
      }


//printf("uPar %u: [%.4e, %.4e]\n",uPar,ParL->at(uPar),ParU->at(uPar));
//printf("   (+): exp in [%i, %i]; mantissa in [%.4f, %.4f]\n",MinExp_Pos[uPar],MaxExp_Pos[uPar],MinManVal_Pos[uPar],MaxManVal_Pos[uPar]);
//printf("   (-): exp in [%i, %i]; mantissa in [%.4f, %.4f]\n",MinExp_Neg[uPar],MaxExp_Neg[uPar],MinManVal_Neg[uPar],MaxManVal_Neg[uPar]);
    }

    //setting NumBestSols initial random solutions
    for(unsigned uBS=0; uBS<NumBestSols*NumBestSols; uBS++){
      Solution->push_back(DLM_FitSolution());
      for(unsigned uPar=0; uPar<NumPars; uPar++){
        if( ParIsFixed(uPar) ){
          //the set value
          Solution->back().Par->push_back(Par->at(uPar));
          //printf(" -> F %f\n",Solution->back().Par->at(uPar));
        }
        else{
          if(ParL->at(uPar)>=0){
            RndSign = 1;
          }
          else if(ParU->at(uPar)<=0){
            RndSign = -1;
          }
          else{
            //how many orders of magnitudes are covered on the positive/negative side
            int PosRange = MaxExp_Pos[uPar]-MinExp_Pos[uPar]+1;
            int NegRange = MaxExp_Neg[uPar]-MinExp_Neg[uPar]+1;
            //printf("+ : - = %i:%i\n",PosRange,NegRange);
            RndSign = RanGen[ThId]->Integer(PosRange+NegRange)<PosRange?1:-1;
          }

          //non-trivial: if we generate a number that has the max exp, than the mantissa should be 0.5 - MaxManVal
          //             if we generate a number that has the min exp, than the mantissa should be MinManVal - 1.0
          if(RndSign==1){
            //this is the added maximum range of the manitssa of each exponent (each mantissa is 0.5 to 1)
            float MntRange = float(MaxExp_Pos[uPar]-MinExp_Pos[uPar]+1)*0.5;
            //the first and last Exp may have a reduced range, leading to an offset
            float MntRangeMin = MinManVal_Pos[uPar] - 0.5;
            float MntRangeMax = MntRange - (1.-MaxManVal_Pos[uPar]);
            //single manitssa range is 0.5, we x2 to get it effectively 1 for ease
            float fRndExp = RanGen[ThId]->Uniform(MntRangeMin,MntRangeMax)*2;
            RndExp = MinExp_Pos[uPar]+floor(fRndExp);
            //RndExp = RanGen[ThId]->Int(MinExp_Pos[uPar],MaxExp_Pos[uPar]);

            if(RndExp==MaxExp_Pos[uPar]&&RndExp!=MinExp_Pos[uPar]){
              RndMnt = RanGen[ThId]->Uniform(0.5,MaxManVal_Pos[uPar]);
            }
            else if(RndExp==MinExp_Pos[uPar]&&RndExp!=MaxExp_Pos[uPar]){
              RndMnt = RanGen[ThId]->Uniform(MinManVal_Pos[uPar],1.0);
            }
            else if(RndExp==MaxExp_Pos[uPar]&&RndExp==MinExp_Pos[uPar]){
              RndMnt = RanGen[ThId]->Uniform(MinManVal_Pos[uPar],MaxManVal_Pos[uPar]);
            }
            else{
              RndMnt = RanGen[ThId]->Uniform(0.5,1.0);
            }
          }
          else{
            //this is the added maximum range of the manitssa of each exponent (each mantissa is 0.5 to 1)
            float MntRange = float(MaxExp_Neg[uPar]-MinExp_Neg[uPar]+1)*0.5;
            //the first and last Exp may have a reduced range, leading to an offset
            float MntRangeMin = MinManVal_Neg[uPar] - 0.5;
            float MntRangeMax = MntRange - (1.-MaxManVal_Neg[uPar]);
            //single manitssa range is 0.5, we x2 to get it effectively 1 for ease
            float fRndExp = RanGen[ThId]->Uniform(MntRangeMin,MntRangeMax)*2;
            RndExp = MinExp_Neg[uPar]+floor(fRndExp);
            //RndExp = RanGen[ThId]->Int(MinExp_Neg[uPar],MaxExp_Neg[uPar]);
            if(RndExp==MaxExp_Neg[uPar]&&RndExp!=MinExp_Neg[uPar]){
              RndMnt = -RanGen[ThId]->Uniform(0.5,MaxManVal_Neg[uPar]);
            }
            else if(RndExp==MinExp_Neg[uPar]&&RndExp!=MaxExp_Neg[uPar]){
              RndMnt = -RanGen[ThId]->Uniform(MinManVal_Neg[uPar],1.0);
            }
            else if(RndExp==MaxExp_Neg[uPar]&&RndExp==MinExp_Neg[uPar]){
              RndMnt = -RanGen[ThId]->Uniform(MinManVal_Neg[uPar],MaxManVal_Neg[uPar]);
            }
            else{
              RndMnt = -RanGen[ThId]->Uniform(0.5,1.0);
            }
          }

          Rnd = RndMnt*pow(2.,RndExp);
          Solution->back().Par->push_back(Rnd);
        }
      }//uPar
    }//uBS

    delete [] MaxAbsVal_Pos;
    delete [] MinAbsVal_Pos;
    delete [] MaxAbsVal_Neg;
    delete [] MinAbsVal_Neg;
    delete [] MinExp_Pos;
    delete [] MaxExp_Pos;
    delete [] MinExp_Neg;
    delete [] MaxExp_Neg;
    delete [] MaxManVal_Pos;
    delete [] MinManVal_Pos;
    delete [] MaxManVal_Neg;
    delete [] MinManVal_Neg;

  }//the random walk done initially

  //this is the case we already have done the initial random sample of parameters
  //i.e. this is now the guided random walk
  else{
    //order the solutions according to their chi2 value
    OrderSolutions();
    if(Solution->back().Chi2 < 0){
      printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) Chi2<0 shows a bug, contact the developers!\n");
      return false;
    }

    //check for convergence
    //the DeltaChi2 corresponding to 1sigma with respect to our best solution
    double DeltaChi2 = GetDeltaChi2(1.,Npar(true));
    //if(Solution->size()>=256){
    //  if(Solution->at(Solution->size()-256).Chi2 - Solution->back().Chi2 <= DeltaChi2) return true;
    //}
    //if( Solution->at(Solution->size()-NumBestSols*NumBestSols).Chi2 - Solution->back().Chi2 <= DeltaChi2) return true;
    //if( Solution->back().Chi2 == Solution->at(Solution->size()-2).Chi2) return true;
//we can make different convergence strategies.
//1. Precision: 1-5, makes the number above fixed to 10^p, demands that we have at least that num of entries
//2. Chi2 (or nsigma): until you reach that number. In case you cant, revert back to 1
//3. Hit: until you manage to go through the errors of all bins. If you fail: revert back to 1

    unsigned BestSolIdUp = Solution->size()-1;
    unsigned BestSolIdLow = Solution->size()-NumBestSols;

//for(DLM_FitSolution& sol : *Solution){
//  printf("SOL chi2 = %e; ps=%u\n",sol.Chi2,sol.Par->size());
//}

    //the best chi2 from all solutions so far
    //double BestChi2 = Solution->at(0).Chi2;

  //printf("DeltaChi2 = %f\n",DeltaChi2);

    //const unsigned NumFreePars = Npar(true);
    //A target is the ID of a solution that we want to use as the 0-th solution for the walk
    //initally, all of the best (NumBestSols) number of solutions are potential targets
    //than, we randomly select only a fraction (NumWildCards) of them. The best solution is always selected.
    std::vector<unsigned> Targets;
    //this is vector that helps us select our targets. At first, it contains all (NumBestSols) solutions
    //and than each time we select a new real target, we push this ID into the Target and erase it from
    //the PotentialTargets, thus by the next random sample we wont duplicate
    std::vector<unsigned> PotentialTargets;
    for(unsigned uBS=BestSolIdLow; uBS<=BestSolIdUp; uBS++){
      PotentialTargets.push_back(uBS);
//printf("PT chi2 = %e\n", Solution->at(uBS).Chi2);
    }
    if(PotentialTargets.size() != NumBestSols){
      printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) PotentialTargets.size() != NumBestSols shows a bug, contact the developers!\n");
      return false;
    }

    for(unsigned uWild=0; uWild<NumWildCards; uWild++){
      if(uWild==0){
        Targets.push_back(PotentialTargets.back());//0
        //PotentialTargets.erase(PotentialTargets.begin());//begin
        PotentialTargets.pop_back();
      }
      else{
        unsigned target = RanGen[ThId]->Integer(PotentialTargets.size());
        Targets.push_back(PotentialTargets.at(target));
        PotentialTargets.erase(PotentialTargets.begin()+target);
      }
    }
//for(unsigned target : Targets){
//printf("TRG %u chi2 = %e\n",target,Solution->at(target).Chi2);
//}
  //printf("Targets in sight\n");
    for(unsigned target : Targets){
      //printf("T%u sz %lu =? %lu\n",target,Solution->at(target).Par->size(),Solution->at(target).Par->size());
      for(unsigned uBS=BestSolIdLow; uBS<=BestSolIdUp; uBS++){
        //we must have two differen elements
        if(uBS==target) continue;
        //we have to check if this combi was already used before.
        //if uBS is a target that was already checked, this is indeed the case
        if(uBS>target && ElementInVector(Targets,uBS)) continue;

        //Solution->at(target) and Solution->at(uBS) represents the two sets of parameters, along the line of which
        //the parameters for the next soltion are evaluated in GetSuitableParameters

        if(Solution->at(uBS).Chi2==Solution->at(target).Chi2){
          printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) sol2.Chi2==sol1.Chi2 shows a bug, contact the developers!\n");
          //return true;
        }

        unsigned SolA = Solution->at(target).Chi2<Solution->at(uBS).Chi2 ? target : uBS;
        unsigned SolB = Solution->at(target).Chi2>Solution->at(uBS).Chi2 ? target : uBS;

        if(Solution->at(SolB).Chi2<Solution->at(SolA).Chi2){
          printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) solB.Chi2<solA.Chi2 shows a bug, contact the developers!\n");
          return false;
        }

        float modif = DeltaChi2/(Solution->at(SolB).Chi2-Solution->at(SolA).Chi2);

        for(unsigned iter=0; iter<1; iter++){
          /*
          //the random magnitude to be used
          //the sign will tell us in which direction to go
          float RndMod = RanGen[ThId]->Gauss(0,2*modif);
          //float RndMod = RanGen[ThId]->Gauss(2.*modif,1*modif);
          //this is the limiting case of getting back to the worst of the two solutions
          //it is pointless to explore this region
          //if(RndMod<=-2*modif) RndMod = -RndMod;
          //float RndMod = RanGen[ThId]->Uniform(-2.0*modif,2*modif);
          //for(int which_way=-1; which_way<=1; which_way+=2){
            Solution->push_back(DLM_FitSolution(Solution->at(SolA)-Solution->at(SolB)));
            //Solution->back() *= ( float(which_way)*RndMod );
            //Solution->push_back(DLM_FitSolution(Solution->at(SolA)));
            Solution->back() *= ( RndMod );
            //Solution->back() *= fabs( RndMod );
            Solution->back().Add(Solution->at(SolA),2.);
            Solution->back().Add(Solution->at(SolB),-1.);
            */
            Solution->push_back(DLM_FitSolution(Solution->at(SolA)-Solution->at(SolB)));
            float RndMod=RanGen[ThId]->Gauss(modif,0.5*modif);
            //going back from the worst of the two solutions, pointless
            //we resample such as to end up by +/-(a-b) away from vector a
            if(RndMod<0){
              RndMod = RanGen[ThId]->Uniform(0,2);
            }

            Solution->back() *= ( RndMod );
            Solution->back().Add(Solution->at(SolB),1.);

            Solution->back().Chi2 = -1;

            ///!!!! NB AT THE END MAKE SURE TO SET THE FIXED PARS EXPLICITELY TO THEIR VALUE
            //NEEDED TO AVOID NUMERICAL FLUCTUATIONS AROUND THE MUST VALUE !!!
            //AS QA VERIFY THAT ACTUALLY THEY ARE AT LEAST CLOSE TO THE MUST VALUE AFTER THIS PROCEDURE ENDS
            for(unsigned uPar=0; uPar<NumPars; uPar++){
              if(ParIsFixed(uPar) && Solution->back().Par->at(uPar)!=Solution->at(SolA).Par->at(uPar)){
                printf("\033[1;33mWARNING:\033[0m Suspicious behaviour DLM_Fit::PrepareForWalk, contact the developers! Report: a fixed par was changed!\n");
                Solution->back().Par->at(uPar) = Solution->at(SolA).Par->at(uPar);
              }
            }
          //}//which_way
        }//iter
      }//uBS
    }//target
  }
  return false;
}

void DLM_Fit::WanderAround(){
  static unsigned COUNTER = 0;
  unsigned CurCntr = 0;
  int SolZero=-1;
  int SolFinal=-1;
  for( int iSol=Solution->size()-1; iSol>=0; iSol--){
    if(Solution->at(iSol).Chi2>=0) break;
    if(SolFinal<0) SolFinal=iSol;
    SolZero=iSol;
  }
  if(SolFinal==-1){
    printf("\033[1;33mWARNING:\033[0m DLM_Fit::WanderAround, does not have any parameters to explore!\n");
    return;
  }

//printf("SolZero = %i / %i %i\n",SolZero,SolFinal,Solution->size());
  if(InterProcess){
    printf("InterProcess NOT DONE YET\n");
  }
  //in case of local evaluation
  else{
    //#pragma omp parallel for
    for( int iSol=SolZero; iSol<=SolFinal; iSol++){
      unsigned ThId = omp_get_thread_num();
      //printf("Probe %i (%u)\n",iSol,ThId);
      //usleep(2000e3);
      ProbePars(Solution->at(iSol).Par,ThId);
      //printf(" done\n");
      //usleep(2000e3);
      Solution->at(iSol).Chi2 = chi2[ThId];
      //printf(" going out\n");
      COUNTER++;
      CurCntr++;
    }
  }
  printf("Goodbye %u <> %u\n",CurCntr,COUNTER);
}

//descending order, i.e. last element is best chi2
void DLM_Fit::OrderSolutions(){
  std::sort(Solution->rbegin(), Solution->rend() );
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

bool DLM_Fit::ParIsFixed(const unsigned& uPar) const{
  //if both of these numbers are zero
  return (ParL->at(uPar)==0 && ParU->at(uPar)==0);
}
bool DLM_Fit::ParIsFree(const unsigned& uPar) const{
  return !ParIsFixed(uPar);
}

unsigned DLM_Fit::Npar(const bool& free_pars) const{
  if(!free_pars) return NumPars;
  unsigned num_pars = NumPars;
  for(unsigned uPar=0; uPar<NumPars; uPar++){
    if( ParIsFixed(uPar) ) num_pars--;
  }
  return num_pars;
}

unsigned DLM_Fit::ThIdBestSol(){
  unsigned thid=0;
  float best_chi2 = 1e37;
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    if(chi2[uTh]<best_chi2 && chi2[uTh]>=0){
      best_chi2 = chi2[uTh];
      thid = uTh;
    }
  }
  ThIdSol = thid;
  return thid;
}




bool DLM_Fit::DEBUG_PrepareForWalk(){
  return PrepareForWalk();
  for(unsigned uSol=0; uSol<Solution->size(); uSol++){
    //printf("sol%u: ",uSol);
    for(unsigned uPar=0; uPar<Solution->at(uSol).Par->size(); uPar++){
      //printf("p%u=%.3e ",uPar,Solution->at(uSol).Par->at(uPar));
    }
    //printf("\n");
  }
}













  DLM_SA_Fit::DLM_SA_Fit(const unsigned& ndata, const unsigned& npars):N_Data(ndata),N_Pars(npars),min_n_grid_pts(30),max_n_grid_pts(268435456),max_n_solutions(8){
    Par = NULL;
    ParMin = NULL;
    ParMax = NULL;
    ParSol = NULL;
    ParErr = NULL;
    input_data = NULL;
    theory_model = NULL;
    fit_min = NULL;
    fit_max = NULL;
    chi2_grid = NULL;
    rangen = new DLM_Random (0);
    par_id = NULL;
    function_type = -1;

    if(N_Data==0||N_Data>128){
      printf("\033[1;31mERROR:\033[0m DLM_SA_Fit needs to have at least 1 data set, and no more than 128\n");
      return;
    }
    if(N_Pars>64){
      printf("\033[1;31mERROR:\033[0m DLM_SA_Fit allows a maximum of 64 fit parameters\n");
      return;
    }    

    Par = new std::vector<float> (N_Pars, 0);
    ParMin = new std::vector<float> (N_Pars, 0);
    ParMax = new std::vector<float> (N_Pars, 0);
    ParSteps = NULL;
    if(N_Pars==1) ParSteps = new std::vector<unsigned> (N_Pars, 20);
    else ParSteps = new std::vector<unsigned> (N_Pars, 5);
    
    par_id = new std::vector<int>* [N_Data];
    input_data = new DLM_Histo<float>* [N_Data];
    theory_model = new SA_FitFunction [N_Data];
    fit_min = new std::vector<float> (N_Data, 0);
    fit_max = new std::vector<float> (N_Data, 0);
    for(unsigned uData=0; uData<N_Data; uData++){
      input_data[uData] = NULL;
      theory_model[uData] = 0;
      par_id[uData] = NULL;
    }

    chi2_grid = new DLM_Histo<float>* [max_n_solutions];
    ParSol = new std::vector<float>* [max_n_solutions];
    ParErr = new std::vector<float>* [max_n_solutions];
    for(unsigned isol=0; isol<max_n_solutions; isol++){
      chi2_grid[isol] = NULL;
      ParSol[isol] = new std::vector<float>(N_Pars, 0);
      ParErr[isol] = new std::vector<float>(N_Pars, 0);
    }

    epsChi2 = 0.25;
    zero_bin_err = 1e-37;
    num_refinements = 20;
    refinment_factor = 4;
    currentChi2 = -1;
    currentBestChi2 = 1e37;
    num_data_pts = 0;
    num_free_fit_pars = 0;
    num_solutions = 0;


    which_bin = new unsigned [N_Pars];
    which_bin_best = new unsigned [N_Pars];
    tot_best_bin = new std::vector<unsigned> (max_n_solutions,0);

  }
  DLM_SA_Fit::~DLM_SA_Fit(){
    if(Par){
      delete Par;
      Par = NULL;
    }
    if(ParMin){
      delete ParMin;
      ParMin = NULL;
    }    
    if(ParMax){
      delete ParMax;
      ParMax = NULL;
    }
    if(ParSteps){
      delete ParSteps;
      ParSteps = NULL;
    }
    if(par_id){
      delete par_id;
      par_id = NULL;
    }
    if(input_data){
      //for(unsigned uData=0; uData<N_Data; uData++){
      //  if(input_data[uData]){
      //  }
      //}
      delete [] input_data;
      input_data = NULL;
    }
    if(theory_model){
      delete [] theory_model;
      theory_model = NULL;
    }
      if(fit_min){
      delete fit_min;
      fit_min = NULL;
    }
    if(fit_max){
      delete fit_max;
      fit_max = NULL;
    }
    if(chi2_grid){
      del_chi2_grid();
    }
    if(ParSol){
      for(unsigned isol=0; isol<max_n_solutions; isol++){
        delete ParSol[isol];
        ParSol[isol] = NULL;
      }
      delete [] ParSol;
      ParSol = NULL;
    }
    if(ParErr){
      for(unsigned isol=0; isol<max_n_solutions; isol++){
        delete ParErr[isol];
        ParErr[isol] = NULL;
      }
      delete [] ParErr;
      ParErr = NULL;
    }
    if(rangen){
      delete rangen;
      rangen = NULL;
    }
    if(which_bin){
      delete [] which_bin;
      which_bin = NULL;
    }
    if(which_bin_best){
      delete [] which_bin_best;
      which_bin_best = NULL;
    }    
    if(tot_best_bin){
      delete tot_best_bin;
      tot_best_bin = NULL;
    }
  }

  void DLM_SA_Fit::del_chi2_grid(){
    if(chi2_grid){
      for(unsigned isol=0; isol<num_solutions; isol++){
        if(chi2_grid[isol]){
          delete chi2_grid[isol];
          chi2_grid[isol] = NULL;
        }
      }
      delete [] chi2_grid;
      chi2_grid = NULL;
      num_solutions = 0;
    }
  }

  void DLM_SA_Fit::SetParameter(const unsigned ipar, const float value){
    if(ipar>=N_Pars || !Par){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParameter given index out of scope\n");
      return;
    }
    Par->at(ipar) = value;
  }
  void DLM_SA_Fit::SetParLimits(const unsigned ipar, const float value_min, const float value_max){
    if(ipar>=N_Pars || !ParMin || !ParMax){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParLimits given index out of scope\n");
      return;
    }
    if(value_min>value_max){
      ParMin->at(ipar) = value_max;
      ParMax->at(ipar) = value_min;
    }
    else{
      ParMin->at(ipar) = value_min;
      ParMax->at(ipar) = value_max;
    }
  }
  void DLM_SA_Fit::FixParameter(const unsigned ipar, const float value){
    SetParLimits(ipar, value, value);
  }


  //void SetParLogScale(const bool yesno);
  //how many evaluation step we will do +/- the default.
  //by default this is 5. The total number of iterations is 1 + par_st*2
  void DLM_SA_Fit::SetParSteps(const unsigned ipar, const unsigned par_st){
    if(ipar>=N_Pars || !ParSteps){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParSteps given index out of scope\n");
      return;
    }    
    if(par_st<=1 || par_st>200){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParSteps can take numbers between [2,200]\n");
      return;
    }
    ParSteps->at(ipar) = par_st;
  }
  //while we fit, is the fitter allowed to make the binning finer to evaluate the parameter 
  //to better precision. A refinment halves the current range. By default this is 7
  void DLM_SA_Fit::SetParRefinements(const unsigned ref_num, const float ref_fac){     
    if(ref_num>100){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParRefinements can take numbers between [0,100]\n");
      return;
    }
    if(ref_fac<1.2 || ref_fac>10){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetParRefinements can take refinement factos in [1.2,10]\n");
      return;
    }    
    num_refinements = ref_num;
    refinment_factor = ref_fac;
  }
  //default 0.25. This parameter gives us a cut-off at which we stop refining the parameter space, given in units 
  //of chi2. I.e. if neighbouring grid pts do not result in larger chi2 differences, we stop the whole fit
  //lower value will give more precision, but will increase computational time.
  void DLM_SA_Fit::SetEpsChi2(const double eps_chi2){
    //if(eps_chi2<0.1 || eps_chi2>1){
    //  printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetEpsChi2 can take numbers between [0.1,1]\n");
    //  return;
    //}
    epsChi2 = eps_chi2;    
  }
  //fit range
  void DLM_SA_Fit::SetDataRange(const unsigned idata, const float value_min, const float value_max){
    if(idata>=N_Data || !fit_min || !fit_max){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetDataRange given index out of scope\n");
      return;
    }
    if(value_min>value_max){
      fit_min->at(idata) = value_max;
      fit_max->at(idata) = value_min;
    }
    else{
      fit_min->at(idata) = value_min;
      fit_max->at(idata) = value_max;      
    }
  }
  void DLM_SA_Fit::SetData(const unsigned idata, DLM_Histo<float>& hdata){
    if(idata>=N_Data || !input_data){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetData given index out of scope\n");
      return;
    }
    if(hdata.GetDim()!=1){
      printf("\033[1;31mERROR:\033[0m Bad input in DLM_SA_Fit::SetData, the histogram has to be 1D\n");
      return;
    }
    input_data[idata] = &hdata;
  }

  void DLM_SA_Fit::SetZeroBinError(const double& zbe){
    zero_bin_err = zbe;
  }

  //the parameter_ids should be a vector, containing the ID of the parameter (as defined within this object)
  //which should be used by the function fit_fun. E.g. parameter_ids = {0,5,6} would imply that fit_fun needs 3 parameters,
  //end they correspond ot the 0-th, 5-th and 6-th parameters set by SetParameter(const unsigned ipar, const float value)
  void DLM_SA_Fit::SetFitFunction(const unsigned idata, SA_FitFunction fit_fun, std::vector<int>& parameter_ids){
    if(idata>=N_Data || !theory_model){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetFitFunction given index out of scope\n");
      if(theory_model){
        theory_model[idata] = 0;
      }
      function_type = -1;
      return;
    }
    
    
    if(par_id[idata]){
      delete par_id[idata];
      par_id[idata] = NULL;
    }
    par_id[idata] = new std::vector<int> (parameter_ids);
    for(unsigned impar=0; impar<parameter_ids.size(); impar++){
      if(par_id[idata]->at(impar)<0 || impar>=N_Pars){
        printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetFitFunction given index out of scope\n");
        theory_model[idata] = 0;
        function_type = -1;
        return;
      }
    }
    
    
    theory_model[idata] = fit_fun;
    function_type = 1;
  }

  //void DLM_SA_Fit::SetFitFunction(SA_FitFunction fit_fun){
  //  if(!theory_model){
  //    printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::SetFitFunction given index out of scope\n");
  //    function_type = -1;
  //    return;
  //  }
  //  theory_model[0] = fit_fun;
  //  function_type = 0;
  //}

  float DLM_SA_Fit::GetParameter(const unsigned& ipar){
    if(ipar>=N_Pars || !Par){
      return 0;
    }
    else{
      //return Par->at(ipar);sol
      return ParSol[0]->at(ipar);
    }
  }
  float DLM_SA_Fit::GetParError(const unsigned& ipar){
    if(ipar>=N_Pars || !Par){
      return 0;
    }
    //sol
    return ParErr[0]->at(ipar);
//TBD
  }
  void DLM_SA_Fit::GetParLimits(unsigned ipar, float& parmin, float& parmax) const{
    if(ipar>=N_Pars || !ParMin || !ParMax){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::GetParLimits given index out of scope\n");
      parmin = 0;
      parmax = 0;
      return;
    }
    parmin = ParMin->at(ipar);
    parmax = ParMax->at(ipar);
  }
  //bool DLM_SA_Fit::GetParLogScale(){

  //}
  unsigned DLM_SA_Fit::GetParSteps(const unsigned ipar){
    if(ipar>=N_Pars || !ParSteps){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::GetParSteps given index out of scope\n");
      return 0;
    }
    return ParSteps->at(ipar);
  }
  unsigned DLM_SA_Fit::GetParRefinements(){
    return num_refinements;
  }
  double DLM_SA_Fit::GetEpsChi2(){
    return epsChi2;
  }
  void DLM_SA_Fit::GetDataRange(unsigned idata, float& datamin, float& datamax) const{
    if(idata>=N_Data || !fit_min || !fit_max){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::GetDataRange given index out of scope\n");
      datamin = 0;
      datamax = 0;
      return;
    }
    datamin = fit_min->at(idata);
    datamax = fit_max->at(idata);
  }

  bool DLM_SA_Fit::GoBabyGo(){
    if(!Par || !ParMin || !ParMax || !ParSteps || !input_data || !fit_min || !fit_max || !theory_model){
      printf("\033[1;31mERROR:\033[0m DLM_SA_Fit::GoBabyGo problematic set up\n");
      return false;
    }

    for(unsigned uData=0; uData<N_Data; uData++){
      if(input_data[uData] == NULL || function_type==-1){
        printf("\033[1;31mERROR:\033[0m DLM_SA_Fit::GoBabyGo problematic set up\n");
        return false;        
      }
    }

    std::vector<double> original_range(N_Pars);
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      original_range.at(ipar) = ParMax->at(ipar)-ParMin->at(ipar);
    }

    //check the number of grid pts
    unsigned TotNumGridPts = 1;
    std::vector<bool> fixed_par(N_Pars);
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      fixed_par.at(ipar) = ParMax->at(ipar)==ParMin->at(ipar);
      if(fixed_par.at(ipar)==true) ParSteps->at(ipar) = 0;
      else TotNumGridPts *= (2*ParSteps->at(ipar) + 1);  
    }
    if(TotNumGridPts>max_n_grid_pts){
      printf("\033[1;31mERROR:\033[0m DLM_SA_Fit::GoBabyGo found way too many grid pts. Cannot proceed. Reduce the points using SetParSteps\n");
      return false;
    }
    if(TotNumGridPts<min_n_grid_pts){
      printf("\033[1;31mERROR:\033[0m DLM_SA_Fit::GoBabyGo found way too few grid pts (%u). Cannot proceed. Increase the poins using SetParSteps to at least %u\n",TotNumGridPts,min_n_grid_pts);
      return false;
    }   

    num_free_fit_pars = 0;
    
    std::vector<int> free_fit_pars;
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      if(fixed_par.at(ipar)==false){
        num_free_fit_pars++;
        free_fit_pars.push_back(ipar);
      }
    }
    if(num_free_fit_pars==0){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::GoBabyGo has zero fit parameters, nothing was done\n");
      return false;      
    }

    //create the actual chi2 grid
    //it is kept as a 1D vector, where the indexing can be evaluated
    del_chi2_grid();
    chi2_grid = new DLM_Histo<float>* [max_n_solutions];
    std::vector<std::vector<double>> refinement_penalty(max_n_solutions);

    for(unsigned isol=0; isol<max_n_solutions; isol++){
      chi2_grid[isol] = NULL;
      refinement_penalty[isol].resize(N_Pars);
      for(unsigned ipar=0; ipar<N_Pars; ipar++){
        refinement_penalty[isol].at(ipar) = 1;
      }
    }
    std::vector<bool> converged_sol(max_n_solutions, false);

    //aim at evaluating (check in that order)
    //a minimum of 1% of the grid points
    //a minimum of 30 points
    //a maximum of 300 points
    //a maximum of TotNumGridPts
    //dont allow more than half of TotNumGridPts, as otherwise you will be stick in random picking
    unsigned num_pts_to_explore = TotNumGridPts/100;
    if(num_pts_to_explore<30) num_pts_to_explore = 30;
    if(num_pts_to_explore>300) num_pts_to_explore = 300;
    if(num_pts_to_explore>TotNumGridPts/2) num_pts_to_explore = TotNumGridPts;
    num_solutions = 1;
    
    std::vector<std::vector<double>> ParLimMin(max_n_solutions);
    std::vector<std::vector<double>> ParLimMax(max_n_solutions);
    std::vector<std::vector<bool>> freeze_parameter(max_n_solutions);
    for(unsigned isol=0; isol<max_n_solutions; isol++){
      ParLimMin[isol].resize(N_Pars);
      ParLimMax[isol].resize(N_Pars);
      freeze_parameter[isol].resize(N_Pars);
      for(unsigned ipar=0; ipar<N_Pars; ipar++){
        ParLimMin[isol].at(ipar) = ParMin->at(ipar);
        ParLimMax[isol].at(ipar) = ParMax->at(ipar);
        freeze_parameter[isol].at(ipar) = false;
      }
    }

    double* par_values = new double [N_Pars];
    //<= as we need to do at least one iteration
    for(unsigned iref=0; iref<=num_refinements; iref++){
printf("iref = %u/%u\n",iref,num_refinements);
      int add_new_solutions = 0;
      for(unsigned isol=0; isol<num_solutions; isol++){
        //here we only keep the best solutions, with a fixed maximal length
        //the 0-th element is always the worst best solution, the last element is the best current solution
        std::vector<unsigned> explored_points;
        std::vector<float> explored_points_chi2;
        const unsigned max_kept_elements = min_n_grid_pts/2;

        if(converged_sol.at(isol)==true) continue;

        if(chi2_grid[isol]){
          delete chi2_grid[isol];
          chi2_grid[isol] = NULL;
        }
        chi2_grid[isol] = new DLM_Histo<float>();
        chi2_grid[isol]->SetUp(N_Pars);
        for(unsigned ipar=0; ipar<N_Pars; ipar++){
          if(fixed_par.at(ipar)==false){
            chi2_grid[isol]->SetUp(ipar, 2*ParSteps->at(ipar)+1, ParLimMin[isol].at(ipar),ParLimMax[isol].at(ipar));
          }
          else{
            double small_eps = (ParLimMax[isol].at(ipar) + ParLimMin[isol].at(ipar))*0.5;
            small_eps *= 1e-4;
            small_eps = fabs(small_eps);
            if( small_eps<1e-32 ) small_eps = 1e-32;
            chi2_grid[isol]->SetUp(ipar, 1, ParLimMin[isol].at(ipar)-small_eps,ParLimMax[isol].at(ipar)+small_eps);
          }
        }
        chi2_grid[isol]->Initialize();
        chi2_grid[isol]->SetBinContentAll(-1);


        double worst_chi2 = 0;
        //the last point is special, where in the case we did not explore all points, we will look at all 
        //neighbours around the best current solution
        std::vector<unsigned> best_neighbours_id(0);
        bool best_neighbours_init = false;
        for(unsigned irnd=0; irnd<=num_pts_to_explore; irnd++){
          unsigned irnd_next = irnd;
          //we have the random explore only if we do not need to check all points
          if(num_pts_to_explore!=TotNumGridPts && irnd!=num_pts_to_explore){
            //at the start we only explore new points, later on we reduce the probability
            float prob_new_explore = 1. - float(irnd)/float(num_pts_to_explore);            
            //we explore a new random point
            if(rangen->Uniform()<prob_new_explore || explored_points.size()==0){
              do irnd_next = rangen->Integer(TotNumGridPts);
              while(chi2_grid[isol]->GetBinContent(irnd_next)!=-1); 
            }//new point
            //we explore around an old (good) solution
            else{
              //this is an explored point, we do now a random walk to find an inexplored one
              irnd_next = explored_points.at(rangen->Integer(explored_points.size()));
              
              unsigned rnd_dim;
              int rnd_dir;
              
              while(chi2_grid[isol]->GetBinContent(irnd_next)!=-1){
                chi2_grid[isol]->GetBinCoordinates(irnd_next, which_bin);
                
                rnd_dim = 0;
                //rnd_dim = rangen->Integer(N_Pars);
                if(num_free_fit_pars>1){
                  rnd_dim = free_fit_pars.at(rangen->Integer(num_free_fit_pars));
                }
                

                //if this is the first/last bin, only one way to go
                     if(which_bin[rnd_dim]==0) rnd_dir = 1;
                else if(which_bin[rnd_dim]==chi2_grid[isol]->GetNbins(rnd_dim)-1) rnd_dir = -1;
                else rnd_dir = -1 + 2*rangen->Integer(2);//+/- 1
                which_bin[rnd_dim] += rnd_dir;
                irnd_next = chi2_grid[isol]->GetTotBin(which_bin);
              }
            }//random walk
          }//explore
          if(irnd==num_pts_to_explore){
            //if we computed it all, no need to redo the neighbours
            if(num_pts_to_explore==TotNumGridPts){continue;}
            if(!best_neighbours_init){
              best_neighbours_id.clear();
              best_neighbours_id = get_all_neighbours(isol, explored_points.back());
//printf("get neighb\n");
              best_neighbours_init = true;
            }
            if(best_neighbours_id.size()==0) continue;
            else{
              irnd_next = best_neighbours_id.back();
              best_neighbours_id.pop_back();
              irnd--;
            }
          }
          if(chi2_grid[isol]->GetBinContent(irnd_next)!=-1){
            continue;
          }
          chi2_grid[isol]->GetBinAxisCenters(irnd_next, par_values);

          currentChi2 = EvalChisquare(par_values);
          chi2_grid[isol]->SetBinContent(irnd_next, currentChi2);
          if(currentChi2<currentBestChi2){
            currentBestChi2 = currentChi2;
            chi2_grid[isol]->GetBinCoordinates(irnd_next, which_bin_best);
            tot_best_bin->at(isol) = irnd_next;
            best_neighbours_init = false;
          }
          if(currentChi2>worst_chi2){
            worst_chi2 = currentChi2;
          }
          //the next bit is AI, it is simply making sure that we keep 
          //explored_points_chi2 as a vector storing the best solutions we have so far, in descending chi2 order
          //explored_points keeps the corresponding bin_ids
          // Ensure the vectors are not empty before checking the front element
          //the last conditions makes sure we dont do this in a all-point scan
          if(explored_points_chi2.empty()){
            explored_points_chi2.push_back(currentChi2);
            explored_points.push_back(irnd_next);
          }
          else{
              //add new element
              if(currentChi2 < explored_points_chi2[0] || explored_points_chi2.size()<max_kept_elements){
                if(currentChi2<explored_points_chi2.back()) best_neighbours_init = false;

                // 1. Pop the worst (largest) chi2 and its corresponding point from the front
                if(explored_points_chi2.size()==max_kept_elements && currentChi2 < explored_points_chi2[0]){
                  explored_points_chi2.erase(explored_points_chi2.begin());
                  explored_points.erase(explored_points.begin());
                }

                // 2. Use binary search (O(log M)) to find where currentChi2 belongs.
                // std::lower_bound with std::greater finds the first element LESS THAN OR EQUAL TO currentChi2.
                auto it_chi2 = std::lower_bound(
                    explored_points_chi2.begin(), 
                    explored_points_chi2.end(), 
                    currentChi2, 
                    std::greater<float>()
                );

                // 3. Calculate the numerical index distance from the beginning
                size_t insert_index = std::distance(explored_points_chi2.begin(), it_chi2);

                // 4. Insert both values at the exact same relative position to maintain parallel alignment
                explored_points_chi2.insert(explored_points_chi2.begin() + insert_index, currentChi2);
                explored_points.insert(explored_points.begin() + insert_index, irnd_next);
              }
          }
//printf("%i chi2 = %f (%f)\n",irnd, chi2_grid[isol]->GetBinContent(explored_points.back()), currentBestChi2);
        }//irnd



        //printf("b chi2 = %f (%f)\n",chi2_grid[isol]->GetBinContent(explored_points.back()), currentBestChi2);
        //chi2_grid[isol]->GetBinCoordinates(explored_points.back(), which_bin_best);


        if(explored_points.size()){
          
          //check if we are done with this solution, i.e. all neighbours of the best solution are close in chi2 space
          //unsigned* which_bin = new unsigned [N_Pars];
          //unsigned* which_bin_best = new unsigned [N_Pars];
          chi2_grid[isol]->GetBinCoordinates(explored_points.back(), which_bin);
          chi2_grid[isol]->GetBinCoordinates(explored_points.back(), which_bin_best);
          tot_best_bin->at(isol) = explored_points.back();
          bool close_chi2 = true;
//printf("which_bin_best");

          for(unsigned ipar=0; ipar<N_Pars; ipar++){
//printf(" %f",which_bin[ipar]);    
            freeze_parameter[isol].at(ipar) = true;
            if(which_bin_best[ipar]>0){
              which_bin[ipar] = which_bin_best[ipar]-1;
//if(chi2_grid[isol]->GetBinContent(which_bin)==-1) {printf("NEGATIVE -\n"); usleep(1000e3);}
//if(chi2_grid[isol]->GetBinContent(which_bin_best)==-1) printf("NEGATIVE B -\n");
//printf(" chi2 %f vs %f\n",chi2_grid[isol]->GetBinContent(which_bin),chi2_grid[isol]->GetBinContent(which_bin_best));
              if(chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best) > epsChi2){
printf("Xdelta%i %f\n",ipar,chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best));
//usleep(1000e3);
                close_chi2 = false;
                freeze_parameter[isol].at(ipar) = false;
              }
              else{
printf("delta %i %f\n",ipar,chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best));
              }
              which_bin[ipar] = which_bin_best[ipar];
            }
            if(which_bin_best[ipar]<2*ParSteps->at(ipar)){
              which_bin[ipar] = which_bin_best[ipar]+1;
//if(chi2_grid[isol]->GetBinContent(which_bin)==-1) {printf("NEGATIVE +\n"); usleep(1000e3);}     
//if(chi2_grid[isol]->GetBinContent(which_bin_best)==-1) printf("NEGATIVE B +\n");      
              if(chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best) > epsChi2){
printf("Ydelta%i %f\n",ipar,chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best));
                close_chi2 = false;
                freeze_parameter[isol].at(ipar) = false;
              }
              else{
printf("deltaa %i %f\n",ipar,chi2_grid[isol]->GetBinContent(which_bin)-chi2_grid[isol]->GetBinContent(which_bin_best));
              }              
              which_bin[ipar] = which_bin_best[ipar];
            }
          }//ipar
//printf("\n");      
//for(unsigned ipar=0; ipar<N_Pars; ipar++) freeze_parameter[isol].at(ipar) = false;


//close_chi2 = false;

          //we are done here
          if(close_chi2){
            converged_sol.at(isol) = true;
//printf("DONE AT %i\n",iref);
            continue;
          }
            
        }







//IT NEVER FINDS ANY ALTERNATIVES
//OKAY, MADE IT WORK, BUT IT DOES NOT WORK
//TIME CONSTRAINTS SAY: LEAVE IT FOR NOW
        //find suitable solutions
        //the best is always taken, but we do allow to split into up to 4
        //creteria to follow up a solution: delta_chi2 (wrt best) < delta_chi2 (wrt worst) * 0.1 (max 5 sigma)
        //and no overalpping points i.e. at least one dimension where the solutions are 3 steps away
        /*
        chi2_grid[isol]->GetBinCoordinates(explored_points.back(), which_bin_best);
        double not_so_bad_delta_chi2 = worst_chi2*0.1;
        if(not_so_bad_delta_chi2>25) not_so_bad_delta_chi2 = 25;
        //dont check the best solution, only the rest
        for(unsigned ipoint=0; ipoint<max_kept_elements-1; ipoint++){
          bool useful_solution = false;
          //printf("AS: %f %f\n",explored_points_chi2.at(ipoint),explored_points_chi2.back());
          if(explored_points_chi2.at(ipoint) - explored_points_chi2.back() < not_so_bad_delta_chi2){
            useful_solution = true;
            chi2_grid[isol]->GetBinCoordinates(explored_points.at(ipoint), which_bin);
            for(unsigned ipar=0; ipar<N_Pars; ipar++){
              //printf("ABS %i\n", abs(int(which_bin_best[ipar])-int(which_bin[ipar])));
              if( abs(int(which_bin_best[ipar])-int(which_bin[ipar]))<1 ){
                useful_solution = false;
                break;
              }
            }
          }
          if(useful_solution && num_solutions+add_new_solutions<max_n_solutions){
            unsigned new_sol_id = num_solutions+add_new_solutions;
            for(unsigned ipar=0; ipar<N_Pars; ipar++){
              ParLimMin[new_sol_id].at(ipar) = chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]) - original_range.at(ipar)*0.5*pow(1./refinment_factor, iref+1)*refinement_penalty[isol].at(ipar);
              if(ParLimMin[new_sol_id].at(ipar) < ParMin->at(ipar)) ParLimMin[new_sol_id].at(ipar) = ParMin->at(ipar);
              ParLimMax[new_sol_id].at(ipar) = chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]) + original_range.at(ipar)*0.5*pow(1./refinment_factor, iref+1)*refinement_penalty[isol].at(ipar);
              if(ParLimMax[new_sol_id].at(ipar) > ParMax->at(ipar)) ParLimMax[new_sol_id].at(ipar) = ParMax->at(ipar);
            }
printf("alternative solution found\n");
usleep(1000e3);
            add_new_solutions++;
          }
        }//ipoint
*/

        double dist = 0;
        double cent_pos = 0;
        double half_len = 0;
        double dist_chi2 = 0;
        chi2_grid[isol]->GetBinCoordinates(explored_points.back(), which_bin_best);
        for(unsigned ipar=0; ipar<N_Pars; ipar++){
          which_bin[ipar] = which_bin_best[ipar];
        }
        for(unsigned ipar=0; ipar<N_Pars; ipar++){
          if(freeze_parameter[isol].at(ipar)==true) continue;
          //we check where the best solution is located. if it is on the edge of the bin, 
          //we reduce the amount by which we shrink the bin at the next step
          cent_pos = (ParLimMax[isol].at(ipar) + ParLimMin[isol].at(ipar))*0.5;
          half_len = (ParLimMax[isol].at(ipar) - ParLimMin[isol].at(ipar))*0.5;
          dist = fabs(cent_pos - chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar])) / half_len;

          if(dist>1){
            dist = 1;
            printf("\033[1;33mWARNING:\033[0m dist>1, points to a bug, contact the developers\n");
          }

          //we do now something similar, but looking at the chi2. If it becomes lower compared to epsChi2,
          //we also try to keep the bin size. This will help with the error estimation
          dist_chi2 = 0;

          if(which_bin_best[ipar]>0){
            which_bin[ipar] = which_bin_best[ipar]-1;
            if(chi2_grid[isol]->GetBinContent(which_bin) - chi2_grid[isol]->GetBinContent(which_bin_best)>dist_chi2){
              dist_chi2 = chi2_grid[isol]->GetBinContent(which_bin) - chi2_grid[isol]->GetBinContent(which_bin_best);
            }
            which_bin[ipar] = which_bin_best[ipar];
          }
          if(which_bin_best[ipar]<2*ParSteps->at(ipar)){
            which_bin[ipar] = which_bin_best[ipar]+1;
            if(chi2_grid[isol]->GetBinContent(which_bin) - chi2_grid[isol]->GetBinContent(which_bin_best)>dist_chi2){
              dist_chi2 = chi2_grid[isol]->GetBinContent(which_bin) - chi2_grid[isol]->GetBinContent(which_bin_best);
            }
            which_bin[ipar] =  which_bin_best[ipar];
          }
  //printf("dist_chi2 %i = %f\n",ipar,dist_chi2);
          dist_chi2 /= epsChi2;
          dist_chi2 = exp(-pow(dist_chi2,2.));

          if(dist_chi2>1){
            dist_chi2 = 1;
            printf("\033[1;33mWARNING:\033[0m dist_chi2>1, points to a bug, contact the developers\n");
          }
//printf("%f %f, ipar%i, dch2 %f\n",dist,dist_chi2,ipar,chi2_grid[isol]->GetBinContent(which_bin) - chi2_grid[isol]->GetBinContent(which_bin_best));
          if(dist_chi2>dist){
            dist = dist_chi2;
          }
//printf(" dist_chi2 = %f\n",dist_chi2);
          refinement_penalty[isol].at(ipar) *= 1+(refinment_factor-1)*dist;

          ParLimMin[isol].at(ipar) = chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]) - original_range.at(ipar)*0.5*pow(1./refinment_factor, iref+1)*refinement_penalty[isol].at(ipar);
          if(ParLimMin[isol].at(ipar) < ParMin->at(ipar)) ParLimMin[isol].at(ipar) = ParMin->at(ipar);
          ParLimMax[isol].at(ipar) = chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]) + original_range.at(ipar)*0.5*pow(1./refinment_factor, iref+1)*refinement_penalty[isol].at(ipar);
          if(ParLimMax[isol].at(ipar) > ParMax->at(ipar)) ParLimMax[isol].at(ipar) = ParMax->at(ipar);

        }


      }//isol
      num_solutions += add_new_solutions;
printf("currentBestChi2 = %f\n",currentBestChi2);
    }//iref
    delete [] par_values;

    printf("best chi2 = %f\n",currentBestChi2);
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      which_bin[ipar] = which_bin_best[ipar];
    }
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      double par_val = chi2_grid[0]->GetBinCenter(ipar, which_bin_best[ipar]);
      //double chi2_val = chi2_grid[0]->GetBinContent(ipar, which_bin_best[ipar]);
      //printf("pv %f %f %f\n", chi2_grid[0]->GetLowEdge(ipar), par_val, chi2_grid[0]->GetUpEdge(ipar));
      which_bin[ipar] = 0;
      double chi2_limL = chi2_grid[0]->GetBinContent(which_bin);
      double chi2_limB = chi2_grid[0]->GetBinContent(which_bin_best);
      which_bin[ipar] = ParSteps->at(ipar)-1;
      double chi2_limU = chi2_grid[0]->GetBinContent(which_bin);
      which_bin[ipar] = which_bin_best[ipar];
      //printf("chi2 %f %f %f\n", chi2_limL, chi2_limB, chi2_limU);
      which_bin[ipar] = which_bin_best[ipar]-1;
      chi2_limL = chi2_grid[0]->GetBinContent(which_bin);
      which_bin[ipar] = which_bin_best[ipar]+1;
      chi2_limU = chi2_grid[0]->GetBinContent(which_bin);
      which_bin[ipar] = which_bin_best[ipar];
      //printf("chi2 %f %f %f\n", chi2_limL, chi2_limB, chi2_limU);
    }


    //fine-tune the values using extrapolation
    const unsigned fine_pts = 100;
    double true_min_chi2 = 1e37;


    
    for(unsigned isol=0; isol<num_solutions; isol++){

      //calculate even more points
      std::vector<unsigned> the_neighborhood = get_the_neighborhood(isol, tot_best_bin->at(isol));
      double* par_values = new double [N_Pars];
      for(unsigned ipoint=0; ipoint<the_neighborhood.size(); ipoint++){
        if(chi2_grid[isol]->GetBinContent(the_neighborhood.at(ipoint))==-1){
          chi2_grid[isol]->GetBinAxisCenters(the_neighborhood.at(ipoint),par_values);
          chi2_grid[isol]->SetBinContent(the_neighborhood.at(ipoint), EvalChisquare(par_values));
        }
      }
      delete [] par_values;




// Get the exact coordinates and chi2 value of the minimum
    chi2_grid[isol]->GetBinCoordinates(tot_best_bin->at(isol), which_bin_best);
    double y0 = chi2_grid[isol]->GetBinContent(tot_best_bin->at(isol));

    // Sync our working bin array to the best bin position
    for (unsigned ipar = 0; ipar < N_Pars; ipar++) {
        which_bin[ipar] = which_bin_best[ipar];
    }

    std::vector<double> errors(N_Pars);
    std::vector<std::vector<double>> H(N_Pars, std::vector<double>(N_Pars, 0.0));
    // 2. Compute the second derivatives (Hessian Matrix elements)
    for (unsigned ipar = 0; ipar < N_Pars; ipar++) {
        double h_i = chi2_grid[isol]->GetBinSize(ipar, which_bin[ipar]);

        // --- Diagonal Elements H[i][i] ---
        which_bin[ipar] = which_bin_best[ipar] + 1;
        double y_plus = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

        which_bin[ipar] = which_bin_best[ipar] - 1;
        double y_minus = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

        which_bin[ipar] = which_bin_best[ipar]; // Reset coordinate i

        H[ipar][ipar] = (y_plus - 2.0 * y0 + y_minus) / (h_i * h_i);

        // --- Off-Diagonal Elements H[i][j] ---
        for (unsigned jpar = ipar + 1; jpar < N_Pars; jpar++) {
            double h_j = chi2_grid[isol]->GetBinSize(jpar, which_bin[jpar]);

            // Bin shifts for the 4 diagonal corners around the minimum
            which_bin[ipar] = which_bin_best[ipar] + 1; which_bin[jpar] = which_bin_best[jpar] + 1;
            double y_pp = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

            which_bin[jpar] = which_bin_best[jpar] - 1;
            double y_pm = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

            which_bin[ipar] = which_bin_best[ipar] - 1; which_bin[jpar] = which_bin_best[jpar] + 1;
            double y_mp = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

            which_bin[jpar] = which_bin_best[jpar] - 1;
            double y_mm = chi2_grid[isol]->GetBinContent(chi2_grid[isol]->GetTotBin(which_bin));

            // Reset coordinates back to the minimum
            which_bin[ipar] = which_bin_best[ipar];
            which_bin[jpar] = which_bin_best[jpar];

            double mixed_derivative = (y_pp - y_pm - y_mp + y_mm) / (4.0 * h_i * h_j);
            H[ipar][jpar] = mixed_derivative;
            H[jpar][ipar] = mixed_derivative; // Symmetric matrix
        }
    }

// 3. Set up an Augmented Matrix [H | I] for Gauss-Jordan Elimination
    // Size is N_Pars rows by 2*N_Pars columns
    std::vector<std::vector<double>> aug(N_Pars, std::vector<double>(2 * N_Pars, 0.0));
    for (unsigned ipar = 0; ipar < N_Pars; ipar++) {
        // FIX: Loop through ALL columns (jpar = 0) to copy the full Hessian matrix
        for (unsigned jpar = 0; jpar < N_Pars; jpar++) {
            aug[ipar][jpar] = H[ipar][jpar];
        }
        aug[ipar][N_Pars + ipar] = 1.0; // Identity matrix on the right side
    }

    // 4. Perform Gauss-Jordan Elimination to invert the matrix
    for (unsigned ipar = 0; ipar < N_Pars; ipar++) {
        // Partial pivoting: find the largest pivot element for numerical stability
        unsigned pivot_row = ipar;
        for (unsigned jpar = ipar + 1; jpar < N_Pars; jpar++) {
            if (std::abs(aug[jpar][ipar]) > std::abs(aug[pivot_row][ipar])) {
                pivot_row = jpar;
            }
        }

        // If the pivot is zero (or dangerously close), the matrix is singular
        if (std::abs(aug[pivot_row][ipar]) < 1e-12) {
            // FIX: Fill with -1.0 and return early to stop division by zero
            for (unsigned k = 0; k < N_Pars; k++) errors[k] = -1.0;
            //return errors; 
        }

        // Swap the current row with the pivot row
        if (pivot_row != ipar) {
            std::swap(aug[ipar], aug[pivot_row]);
        }

        // Scale the pivot row so the leading coefficient is 1.0
        double pivot_val = aug[ipar][ipar];
        for (unsigned jpar = 0; jpar < 2 * N_Pars; jpar++) {
            aug[ipar][jpar] /= pivot_val;
        }

        // Eliminate the current column from all other rows
        for (unsigned jpar = 0; jpar < N_Pars; jpar++) {
            if (jpar != ipar) {
                double factor = aug[jpar][ipar];
                for (unsigned kpar = 0; kpar < 2 * N_Pars; kpar++) {
                    aug[jpar][kpar] -= factor * aug[ipar][kpar];
                }
            }
        }
    }

// 5. Extract the true marginalized 1-sigma errors from the inverted matrix side
    for (unsigned ipar = 0; ipar < N_Pars; ipar++) {
        // FIX: Multiply by 2.0 because Covariance = 2 * Inverse(Hessian)
        double variance = 2.0 * aug[ipar][N_Pars + ipar]; 
        
        if (variance > 0.0) {
            errors[ipar] = std::sqrt(variance);
        } else {
            errors[ipar] = -1.0; 
        }
        printf("%i val %e err %e\n", ipar, chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]), errors[ipar]);
    }



    //numerical error based on what is on the grid (isol)
    std::vector<unsigned> best_neighbours_id = get_the_neighborhood(0, tot_best_bin->at(0));
    std::vector<double> largest_error = errors;
    par_values = new double [N_Pars];
    
    for(unsigned ineigh=0; ineigh<best_neighbours_id.size(); ineigh++){
      chi2_grid[0]->GetBinAxisCenters(best_neighbours_id.at(ineigh), par_values);
      for(unsigned ipar=0; ipar<N_Pars; ipar++){

        double delta_chi2 = chi2_grid[0]->GetBinContent(best_neighbours_id.at(ineigh)) - currentBestChi2;
        double par_diff = fabs(par_values[ipar] - chi2_grid[0]->GetBinCenter(ipar, which_bin_best[ipar]));
        //comes from definition of chi2
        double par_err = par_diff/sqrt(delta_chi2);
        
        //if(par_err > largest_error.at(ipar)){
          largest_error.at(ipar) = par_err;
        //}
      }
    }
    delete [] par_values;

    //
    printf("grid estimate of the error:\n");
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      printf(" %i val %e err %e\n", ipar, chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]), largest_error.at(ipar));
    }

    //sol missing
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      ParSol[isol]->at(ipar) = chi2_grid[isol]->GetBinCenter(ipar, which_bin_best[ipar]);
      ParErr[isol]->at(ipar) = largest_error.at(ipar);
    }
    


      /*
      DLM_Histo<float> dlm_fine_chi2;
      dlm_fine_chi2.SetUp(N_Pars);
      for(unsigned ipar=0; ipar<N_Pars; ipar++){
        double low_edge = chi2_grid[isol]->GetBinLowEdge(ipar,which_bin_best[ipar]);
        double up_edge = chi2_grid[isol]->GetBinUpEdge(ipar,which_bin_best[ipar]);
        dlm_fine_chi2.SetUp(ipar, fine_pts, low_edge, up_edge);
      }
      dlm_fine_chi2.Initialize();

      double* par_val = new double [N_Pars];

      for(unsigned uBin=0; uBin<dlm_fine_chi2.GetNbins(); uBin++){
        dlm_fine_chi2.GetBinAxisCenters(uBin, par_val);
        dlm_fine_chi2.SetBinContent(uBin, chi2_grid[isol]->Eval(par_val));
        if(dlm_fine_chi2.GetBinContent(uBin)<true_min_chi2){
          true_min_chi2 = dlm_fine_chi2.GetBinContent(uBin);
        }
      }
      
      delete [] par_val;
      */
    }
    //printf("true_min_chi2 = %f\n",true_min_chi2);

    return true;
  }

  float DLM_SA_Fit::EvalChisquare(double* pars){
    if(!input_data) return -1;
    if(!theory_model || function_type==-1) return -1;
    double chi2_val = 0;
    double data_val, theory_val, err_val, var_val;
    num_data_pts = 0;
    for(unsigned idata=0; idata<N_Data; idata++){
      if(!input_data[idata]) return -1;
      double* model_pars = new double [par_id[idata]->size()];
      for(unsigned impar=0; impar<par_id[idata]->size(); impar++){
        model_pars[impar] = pars[par_id[idata]->at(impar)];
      }
      for(unsigned uBin=0; uBin<input_data[idata]->GetNbins(); uBin++){
        var_val = input_data[idata]->GetBinCenter(0, uBin);
        data_val = input_data[idata]->GetBinContent(uBin);
        err_val = input_data[idata]->GetBinError(uBin);
        if(var_val < fit_min->at(idata) || var_val > fit_max->at(idata)){
          continue;
        }
        if(data_val==0 && err_val==0){
          err_val = zero_bin_err;
        }
        else if(err_val==0){
          err_val = 1e-37;
        }
        if(err_val==0) continue;
        num_data_pts++;
        if(function_type==0){
          theory_val = theory_model[0](&var_val, model_pars);
        }
        else{
          theory_val = theory_model[idata](&var_val, model_pars);
        }
        //if(uBin==0){
        //  printf("%f %f %f %f\n",var_val,data_val,err_val,theory_val);
        //}
        
        chi2_val += pow((data_val - theory_val)/err_val,2.);
      }
      delete [] model_pars;
    }
    return float(chi2_val);
  }

  float DLM_SA_Fit::GetChisquare(){
    return currentBestChi2;
  }
  int DLM_SA_Fit::GetNDF(){
    return int(num_data_pts) - int(num_free_fit_pars);
  }
  unsigned DLM_SA_Fit::GetNumberFitPoints(){
    return num_data_pts;
  }
  unsigned DLM_SA_Fit::GetNumberFreeParameters(){
    return num_free_fit_pars;
  }
  float DLM_SA_Fit::GetProb(){
//TBD
  }

  void DLM_SA_Fit::SetRandomSeed(unsigned rnd_seed){
    if(!rangen){
      rangen = new DLM_Random(rnd_seed);
    }
    else{
      rangen->SetSeed(rnd_seed);
    }
  }

  //neighbours shidted by 1 unit in 1 of the par space
  std::vector<unsigned> DLM_SA_Fit::get_all_neighbours(unsigned isol, unsigned tot_bin_id){
    std::vector<unsigned> coords(0);
    chi2_grid[isol]->GetBinCoordinates(tot_bin_id, which_bin_best);
    chi2_grid[isol]->GetBinCoordinates(tot_bin_id, which_bin);
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      if(which_bin_best[ipar]>0){
        which_bin[ipar] = which_bin_best[ipar]-1;
        coords.push_back(chi2_grid[isol]->GetTotBin(which_bin));
        which_bin[ipar] = which_bin_best[ipar];
      }
      if(which_bin_best[ipar]<chi2_grid[isol]->GetNbins(ipar)-1){
        which_bin[ipar] = which_bin_best[ipar]+1;
        coords.push_back(chi2_grid[isol]->GetTotBin(which_bin));
        which_bin[ipar] = which_bin_best[ipar];
      }
    }
    return coords;
  }


  //neighbours shifted by 1 unit in ANY of the par space
  std::vector<unsigned> DLM_SA_Fit::get_the_neighborhood(unsigned isol, unsigned tot_bin_id){
    std::vector<unsigned> coords(0);
    chi2_grid[isol]->GetBinCoordinates(tot_bin_id, which_bin_best);
    //chi2_grid[isol]->GetBinCoordinates(tot_bin_id, which_bin);
    for(unsigned uBin=0; uBin<chi2_grid[isol]->GetNbins(); uBin++){
      chi2_grid[isol]->GetBinCoordinates(uBin, which_bin);
      bool in_the_neighborhood = true;
      for(unsigned ipar=0; ipar<N_Pars; ipar++){
        if( abs(int(which_bin[ipar])-int(which_bin_best[ipar]))>1 ){
          in_the_neighborhood = false;
          break;
        }
      }
      if(in_the_neighborhood){
        coords.push_back(chi2_grid[isol]->GetTotBin(which_bin));
      }
    }
    return coords;
  }


/*
 * @brief Maps a flat 1D linear index into an N-dimensional vector of bin indices.
 * * Convention: The last dimension changes the fastest in the 1D grid layout.
 * * @param linear_index The 1D global index on your grid.
 * @return std::vector<unsigned> N-dimensional coordinate indices, where element i is in [0, N_i).
 */
 /*
std::vector<unsigned> DLM_SA_Fit::unflatten_index(unsigned linear_index){
    
    std::vector<unsigned> null_vector(0);
    if(!ParSteps){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::unflatten_index cannot be called before intializing the parameters\n");
      
      return null_vector;
    }  
    std::vector<unsigned> coords(N_Pars);

    unsigned TotNumGridPts=1;
    for(unsigned ipar=0; ipar<N_Pars; ipar++){
      TotNumGridPts *= (2*ParSteps->at(ipar)+1);
    }
    if(!TotNumGridPts){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::unflatten_index says the parameter space contains a dimension with zero parameters\n");
      return null_vector;
    }
    if(linear_index>=TotNumGridPts){
      printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::unflatten_index has a linear_index that is out of scope\n");
      return null_vector;
    }    

    unsigned current = linear_index;
    
    // Loop backwards from the last dimension to the first.
    // Using 'i' from dimensions down to 1 to safely handle underflows in loops.
    for (unsigned ipar = N_Pars; ipar > 0; ipar--) {
        unsigned dim_idx = ipar - 1;
        
        // The remainder gives the coordinate for the current dimension
        coords[dim_idx] = current % ParSteps->at(dim_idx);
        
        // The quotient is passed to the next slower-changing dimension
        current /= ParSteps->at(dim_idx);
    }
    
    return coords;
    
}

unsigned DLM_SA_Fit::flatten_index(const std::vector<unsigned>& coords) {
    
    if(!ParSteps) {
        printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::flatten_index cannot be called before initializing the parameters\n");
        return 0;
    }  

    // Check if the input vector size matches the number of parameters
    if(coords.size() != N_Pars) {
        printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::flatten_index dimension mismatch. Expected %u, got %lu\n", 
               N_Pars, coords.size());
        return 0;
    }

    unsigned linear_index = 0;

    // Loop forward through the dimensions (from slowest changing to fastest changing).
    // Math behind this: Horner's method evaluates the index without explicit nested products.
    // Index = (...((coord[0] * N_1 + coord[1]) * N_2 + coord[2]) * ...) * N_(n-1) + coord[n-1]
    for(unsigned ipar = 0; ipar < N_Pars; ipar++) {
        
        // Validation: Ensure the requested index is strictly within the allocated steps for this parameter
        if(coords.at(ipar) >= ParSteps->at(ipar)) {
            printf("\033[1;33mWARNING:\033[0m DLM_SA_Fit::flatten_index coordinate for parameter %u is out of scope (%u >= %u)\n", 
                   ipar, coords.at(ipar), ParSteps->at(ipar));
            return 0;
        }

        linear_index = linear_index * ParSteps->at(ipar) + coords.at(ipar);
    }
    
    return linear_index;
}
*/