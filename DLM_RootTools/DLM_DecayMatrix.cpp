
#include "DLM_DecayMatrix.h"

#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH2F.h"


DLM_DecayMatrix::DLM_DecayMatrix(){
    //UseGeV = false;
    FileName = new char [128];
    HistoName = new char [128];
    strcpy(FileName,"./DLM_DecayMatrix.root");
    strcpy(HistoName,"./DLM_DecayMatrix.root");
    NumBins = 1000;
    kMin = 0;
    kMax = 1000;
    Mass_Daughter1 = NULL;
    Mass_Daughter2 = NULL;
    SetNumDaughters1(2);
    SetDaughterMass1(0,0);
    SetDaughterMass1(1,0);
    SetNumDaughters2(2);
    SetDaughterMass2(0,0);
    SetDaughterMass2(1,0);
    Mass1 = 134.976;
    Mass2 = 134.976;
    MomMean = 0;
    MomSpread = 350;
}

DLM_DecayMatrix::~DLM_DecayMatrix(){
    if(Mass_Daughter1) delete [] Mass_Daughter1;
    if(Mass_Daughter2) delete [] Mass_Daughter2;
}

double DLM_DecayMatrix::relKcalc(const TLorentzVector& track1, const TLorentzVector& track2){

  TLorentzVector trackSum, track1_cms, track2_cms;
  trackSum = track1 + track2;

  double beta = trackSum.Beta();
  double beta_x = beta*cos(trackSum.Phi())*sin(trackSum.Theta());
  double beta_y = beta*sin(trackSum.Phi())*sin(trackSum.Theta());
  double beta_z = beta*cos(trackSum.Theta());

  track1_cms = track1;
  track2_cms = track2;

  track1_cms.Boost(-beta_x,-beta_y,-beta_z);
  track2_cms.Boost(-beta_x,-beta_y,-beta_z);

  TLorentzVector track_relK;

  track_relK = track1_cms - track2_cms;
  double relK = 0.5*track_relK.P();


  return relK;
}

void DLM_DecayMatrix::SetFileName(const char* name){
    strcpy(FileName,name);
}
void DLM_DecayMatrix::SetHistoName(const char* name){
    strcpy(HistoName,name);
}
void DLM_DecayMatrix::SetBins(const unsigned& numbins, const unsigned& kmin, const unsigned& kmax){
    NumBins = numbins;
    kMin = kmin;
    kMax = kmax;
}
void DLM_DecayMatrix::SetNumDaughters1(const unsigned& numdaughters){
    if(NumDaughters1==numdaughters) return;
    if(Mass_Daughter1) delete [] Mass_Daughter1;
    NumDaughters1 = numdaughters;
    Mass_Daughter1 = new double [NumDaughters1];
}
void DLM_DecayMatrix::SetDaughterMass1(const unsigned& daughter, const double& mass){
    if(daughter>=NumDaughters1){
        printf("\033[1;31mERROR:\033[0m DLM_DecayMatrix: daughter>=NumDaughters1");
        return;
    }
    Mass_Daughter1[daughter] = mass;
}
void DLM_DecayMatrix::SetMotherMass1(const double& mass){
    Mass1 = mass;
}
void DLM_DecayMatrix::SetNumDaughters2(const unsigned& numdaughters){
    if(NumDaughters2==numdaughters) return;
    if(Mass_Daughter2) delete [] Mass_Daughter2;
    NumDaughters2 = numdaughters;
    Mass_Daughter2 = new double [NumDaughters2];
}
void DLM_DecayMatrix::SetDaughterMass2(const unsigned& daughter, const double& mass){
    if(daughter>=NumDaughters2){
        printf("\033[1;31mERROR:\033[0m DLM_DecayMatrix: daughter>=NumDaughters2");
        return;
    }
    Mass_Daughter2[daughter] = mass;
}
void DLM_DecayMatrix::SetMotherMass2(const double& mass){
    Mass2 = mass;
}
void DLM_DecayMatrix::SetMeanMomentum(const double& mean){
    MomMean = mean;
}
void DLM_DecayMatrix::SetMomentumSpread(const double& spread){
    MomSpread = spread;
}
//void DLM_DecayMatrix::SetUnitsMeV(){
//    UseGeV = false;
//}
//void DLM_DecayMatrix::SetUnitsGeV(){
//    UseGeV = true;
//}

void DLM_DecayMatrix::Run(const int& SEED, const unsigned& NumIter)
{
    TRandom3 rangen(SEED);

    TFile *outputfile = new TFile(FileName,"update");
    if(!outputfile) outputfile = new TFile(FileName,"recreate");

    TH2F* hRes = new TH2F(HistoName,"x: relk mothers, y: relk daughters",NumBins,kMin,kMax,NumBins,kMin,kMax);

    TLorentzVector tvec_m1, tvec_m2;
    TLorentzVector* tvec_d1;
    TLorentzVector* tvec_d2;
    double relK_m;
    double relK_d;
    TGenPhaseSpace decay1;
    TGenPhaseSpace decay2;

    for(unsigned uIter=0;uIter<NumIter;uIter++){
        tvec_m1.SetXYZM(rangen.Gaus(MomMean,MomSpread),rangen.Gaus(MomMean,MomSpread),rangen.Gaus(MomMean,MomSpread),Mass1);
        tvec_m2.SetXYZM(rangen.Gaus(MomMean,MomSpread),rangen.Gaus(MomMean,MomSpread),rangen.Gaus(MomMean,MomSpread),Mass2);
        relK_m = relKcalc(tvec_m1,tvec_m2);
        if(NumDaughters1>1){
            decay1.SetDecay(tvec_m1, NumDaughters1, Mass_Daughter1);
            decay1.Generate();
            tvec_d1 = decay1.GetDecay(0);
        }
        else{
            tvec_d1 = new TLorentzVector(tvec_m1);
        }

        if(NumDaughters2>1){
            decay2.SetDecay(tvec_m2, NumDaughters2, Mass_Daughter2);
            decay2.Generate();
            tvec_d2 = decay2.GetDecay(0);
        }
        else{
            tvec_d2 = new TLorentzVector(tvec_m2);
        }
        relK_d = relKcalc(*tvec_d1,*tvec_d2);
        hRes->Fill(relK_m,relK_d);

        if(NumDaughters1<=1){
            delete tvec_d1;
        }
        if(NumDaughters2<=1){
            delete tvec_d2;
        }
    }

    outputfile->cd();

    hRes->Write("",TObject::kOverwrite);

    delete hRes;
    delete outputfile;

}
