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
