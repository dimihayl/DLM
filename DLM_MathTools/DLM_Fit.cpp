#include "DLM_Fit.h"
#include <iostream>
#include "omp.h"
#include <unistd.h>
#include <thread>

#include "DLM_Histo.h"
#include "DLM_Random.h"
#include "DLM_CppTools.h"
#include "DLM_MathFunctions.h"


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
  Par = new std::vector<float> ();
  ParL = new std::vector<float> ();
  ParU = new std::vector<float> ();
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
  if(Par){delete[]Par;Par=NULL;}
  if(ParL){delete[]ParL;ParL=NULL;}
  if(ParU){delete[]ParU;ParU=NULL;}
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
  Par->resize(numpars);
  ParL->resize(numpars);
  ParU->resize(numpars);
  for(float& par : *Par){
    par = 0;
  }
  for(float& lim : *ParL){
    lim = -1e64;
  }
  for(float& lim : *ParU){
    lim = 1e64;
  }

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

//we will walk along "lines" in the parameter space, connecting two solutions (0 and i).
void DLM_Fit::PrepareForWalk(){
  unsigned ThId = omp_get_thread_num();

  //completely random
  if(Solution->size()<NumBestSols){
    if(Solution->size()){
      printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) Solution->size() shows a bug, contact the developers!\n");
      return;
    }

    //setting NumBestSols initial random solutions
    for(unsigned uBS=0; uBS<NumBestSols; uBS++){
      Solution->push_back(DLM_FitSolution());
      for(unsigned uPar=0; uPar<NumPars; uPar++){
        if( ParIsFixed(uPar) ){
          //the set value
          Solution->back().Par->push_back(Par->at(uPar));
        }
        else{
          //a uniform random value within the limits of the parameter
          Solution->back().Par->push_back(RanGen[ThId]->Uniform(ParL->at(uPar),ParU->at(uPar)));
        }
      }
    }

    return;
  }

  //order the solutions according to their chi2 value
  OrderSolutions();
  //the best chi2 from all solutions so far
  //double BestChi2 = Solution->at(0).Chi2;
  //the DeltaChi2 corresponding to 1sigma with respect to our best solution
  double DeltaChi2 = GetDeltaChi2(1.,Npar(true));
  //const unsigned NumFreePars = Npar(true);
  //A target is the ID of a solution that we want to use as the 0-th solution for the walk
  //initally, all of the best (NumBestSols) number of solutions are potential targets
  //than, we randomly select only a fraction (NumWildCards) of them. The best solution is always selected.
  std::vector<unsigned> Targets;
  //this is vector that helps us select our targets. At first, it contains all (NumBestSols) solutions
  //and than each time we select a new real target, we push this ID into the Target and erase it from
  //the PotentialTargets, thus by the next random sample we wont duplicate
  std::vector<unsigned> PotentialTargets;
  for(unsigned uBS=0; uBS<NumBestSols; uBS++){
    PotentialTargets.push_back(uBS);
  }

  for(unsigned uWild=0; uWild<NumWildCards; uWild++){
    if(uWild==0){
      Targets.push_back(0);
      PotentialTargets.erase(PotentialTargets.begin());
    }
    else{
      unsigned target = RanGen[ThId]->Integer(PotentialTargets.size());
      Targets.push_back(PotentialTargets.at(target));
      PotentialTargets.erase(PotentialTargets.begin()+target);
    }
  }

  for(unsigned target : Targets){
    DLM_FitSolution& sol1 = Solution->at(target);
    for(unsigned uBS=0; uBS<NumBestSols; uBS++){
      //we must have two differen elements
      if(uBS==target) continue;
      //we have to check if this combi was already used before.
      //if uBS is a target that was already checked, this is indeed the case
      if(uBS<target && ElementInVector(Targets,uBS)) continue;
      DLM_FitSolution& sol2 = Solution->at(uBS);
      //sol1 and sol2 represents the two sets of parameters, along the line of which
      //the parameters for the next soltion are evaluated in GetSuitableParameters

      float L_i = 0;
      for(unsigned uPar=0; uPar<NumPars; uPar++){
        L_i += pow(sol1.Par->at(uPar)-sol2.Par->at(uPar),2.);
      }
      L_i = sqrt(L_i);

      if(sol2.Chi2==sol1.Chi2){
        printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) sol2.Chi2==sol1.Chi2 shows a bug, contact the developers!\n");
        return;
      }
      float L_tot = DeltaChi2/(sol2.Chi2-sol1.Chi2)*L_i;
      if(L_tot<0){
        printf("\033[1;31mFATAL ERROR:\033[0m (DLM_Fit::PrepareForWalk) L_tot<0 shows a bug, contact the developers!\n");
        return;
      }

      //the random magnitude to be used
      //the sign will tell us in which direction to go
      float RndMag = RanGen[ThId]->Gauss(0,L_tot);

      //Solution->push_back(sol1.GetSuitableParameters(sol2,NumFreePars,RanGen[ThId]));

      ///!!!! NB AT THE END MAKE SURE TO SET THE FIXED PARS EXPLICITELY TO THEIR VALUE
      //NEEDED TO AVOID NUMERICAL FLUCTUATIONS AROUND THE MUST VALUE !!!
      //AS QA VERIFY THAT ACTUALLY THEY ARE AT LEAST CLOSE TO THE MUST VALUE AFTER THIS PROCEDURE ENDS
    }
  }
}

void DLM_Fit::WanderAround(){

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
