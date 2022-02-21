
#include "TREPNI.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "DLM_CppTools.h"
#include "DLM_Histo.h"
#include "DLM_Random.h"
#include "DLM_MathFunctions.h"
#include "CATStools.h"
#include "CATSconstants.h"

TreParticle::TreParticle(TREPNI& database):Database(database){
  TreName = new char [Database.Len_PrtclName];
  Mass = new float [3];
  Width = new float [3];
  Abundance = new float [3];
  Radius = new float [3];
  RadSlope = new float [3];
  DelayTau = new float [3];
  for(char i=0; i<3; i++){
    Mass[i]=0;Width[i]=0;Abundance[i]=0;
    Radius[i]=0;RadSlope[i]=0;
    DelayTau[i]=0;
  }
  NumDecays = 0;
  Decay = NULL;
  //CurrentDecay = NULL;
  PtEtaPhi = NULL;
  PxPyPz_Width = 0;
  Acceptance_pT[0] = 0;
  Acceptance_pT[1] = 1e128;
  Acceptance_Eta[0] = -1e128;
  Acceptance_Eta[1] = 1e128;
  Acceptance_CosTh[0] = -1;
  Acceptance_CosTh[1] = 1;
  Acceptance_Phi[0] = 0;
  Acceptance_Phi[1] = 2.*Pi;
}

TreParticle::~TreParticle(){
  delete [] TreName; TreName=NULL;
  delete [] Mass; Mass=NULL;
  delete [] Width; Width=NULL;
  delete [] Abundance; Abundance=NULL;
  delete [] Radius; Radius=NULL;
  delete [] RadSlope; RadSlope=NULL;
  delete [] DelayTau; DelayTau=NULL;
  if(Decay){
    for(unsigned char uDec=0; uDec<NumDecays; uDec++){
      if(Decay[uDec]){delete Decay[uDec];Decay[uDec]=NULL;}
    }
    delete[]Decay;
    Decay=NULL;
  }
  if(PtEtaPhi){
    delete PtEtaPhi;
    PtEtaPhi = NULL;
  }
}

void TreParticle::SetPtEtaPhi(const DLM_Histo<float>& pdf){
  if(pdf.GetDim()>3){
    static bool ShowMessage=true;
    if(Database.PrintLevel>=1 && ShowMessage){
      printf("\033[1;31mERROR:\033[0m (TreParticle::SetPtEtaPhi) The momentum distribution must have 1,2 or 3 dimensions\n");
    }
    return;
  }
  if(PtEtaPhi){delete PtEtaPhi;}
  PtEtaPhi = new DLM_Histo<float>(pdf);
}
void TreParticle::SetPtEtaPhi(const float& width){
  if(PtEtaPhi){delete PtEtaPhi; PtEtaPhi = NULL;}
  PxPyPz_Width = width;
}

void TreParticle::SamplePxPyPz(double* axisValues, DLM_Random* RanGen, const bool& UnderOverFlow) const{
  if(!RanGen) RanGen = Database.RanGen;
  if(PtEtaPhi){
    PtEtaPhi->Sample(axisValues,UnderOverFlow,RanGen);
    double pt = axisValues[0];
    double eta = axisValues[1];
    double phi = axisValues[2];
    double& px = axisValues[0];
    double& py = axisValues[1];
    double& pz = axisValues[2];
    double cos_th,sin_th,cotg_th,ptot;

    bool Auto_eta = true;
    bool Auto_phi = true;
    if(PtEtaPhi->GetDim()>1)
      Auto_eta = false;
    if(PtEtaPhi->GetDim()>2)
      Auto_phi = false;
    //PtEtaPhi->Sample(axisValues,true,RanGen[ThId]);

    if(Auto_phi){
      phi = RanGen->Uniform(Acceptance_Phi[0],Acceptance_Phi[1]);
    }
    if(Auto_eta){
      cos_th = RanGen->Uniform(Acceptance_CosTh[0],Acceptance_CosTh[1]);
      //sin always positive here
      sin_th = sqrt(1.-cos_th*cos_th);
      cotg_th = cos_th/sin_th;
      //verified that this relation is true
      //eta = -0.5*log((1.-cos_th)/(1.+cos_th));
    }
    else{
      sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
      cotg_th = (1.-exp(-2.*eta))/(2.*exp(-eta));
      cos_th = (1.-exp(-2.*eta))/(1.+exp(-2.*eta));
    }
    pz = pt*cotg_th;
    ptot = sqrt(pt*pt+pz*pz);
    px = ptot*cos(phi)*sin_th;
    py = ptot*sin(phi)*sin_th;
//pT = sqrt(px2+py2) = ptot*|sin_th|
//theta, phi, pT = pz*tg_th. We can sample pT from a Gauss of mean mu and sigma = sigma0 * tg_th
  }
  else{
    double& px = axisValues[0];
    double& py = axisValues[1];
    double& pz = axisValues[2];

    double cos_th = RanGen->Uniform(Acceptance_CosTh[0],Acceptance_CosTh[1]);
    double sin_th = sqrt(1.-cos_th*cos_th)+1e-128;
    double tan_th = sin_th/cos_th;
    double phi = RanGen->Uniform(Acceptance_Phi[0],Acceptance_Phi[1]);
    double ptot,pt;
    do{
      ptot = sqrt(pow(RanGen->Gauss(0,PxPyPz_Width),2.)+pow(RanGen->Gauss(0,PxPyPz_Width),2.)+pow(RanGen->Gauss(0,PxPyPz_Width),2.));
      pt = ptot*sin_th;
    }
    while(pt<Acceptance_pT[0]||pt>Acceptance_pT[1]);
    px = ptot*cos(phi)*sin_th;
    py = ptot*sin(phi)*sin_th;
    pz = pt/tan_th;
  }
}

void TreParticle::SetAcceptance_pT(const float& min, const float& max){
  if(min>=max||max<=0){
    printf("\033[1;33mWARNING:\033[0m (TreParticle::SetAcceptance_pT) Zero acceptance set\n");
    Acceptance_pT[0] = 0;
    Acceptance_pT[1] = 0;
  }
  Acceptance_pT[0] = min;
  Acceptance_pT[1] = max;
}

void TreParticle::SetAcceptance_Eta(const float& min, const float& max){
  if(min>=max){
    printf("\033[1;33mWARNING:\033[0m (TreParticle::SetAcceptance_Eta) Zero acceptance set\n");
    Acceptance_Eta[0] = 0;
    Acceptance_Eta[1] = 0;
    Acceptance_CosTh[0] = 0;
    Acceptance_CosTh[1] = 0;
  }
  Acceptance_Eta[0] = min;
  Acceptance_Eta[1] = max;
  Acceptance_CosTh[0] = (1.-exp(-2.*min)/(1.+exp(-2.*min)));
  Acceptance_CosTh[1] = (1.-exp(-2.*max)/(1.+exp(-2.*max)));
}

void TreParticle::SetAcceptance_Phi(const float& min, const float& max){
  if(min>=max||max<=0){
    printf("\033[1;33mWARNING:\033[0m (TreParticle::SetAcceptance_Phi) Zero acceptance set\n");
    Acceptance_Phi[0] = 0;
    Acceptance_Phi[1] = 0;
  }
  Acceptance_Phi[0] = min<0?0:min;
  Acceptance_Phi[1] = max>2.*Pi?2.*Pi:max;
}

double TreParticle::AcceptanceMin_pT() const{
  return Acceptance_pT[0];
}
double TreParticle::AcceptanceMax_pT() const{
  return Acceptance_pT[1];
}
double TreParticle::AcceptanceMin_Eta() const{
  return Acceptance_Eta[0];
}
double TreParticle::AcceptanceMax_Eta() const{
  return Acceptance_Eta[1];
}
double TreParticle::AcceptanceMin_CosTh() const{
  return Acceptance_CosTh[0];
}
double TreParticle::AcceptanceMax_CosTh() const{
  return Acceptance_CosTh[1];
}
double TreParticle::AcceptanceMin_Phi() const{
  return Acceptance_Phi[0];
}
double TreParticle::AcceptanceMax_Phi() const{
  return Acceptance_Phi[1];
}


DLM_Histo<float>* TreParticle::GetPtEtaPhi() const{
  return PtEtaPhi;
}

void TreParticle::FillPxPyPz(const float& xval, const float& yval, const float& zval){

}

void TreParticle::FillPtEtaPhi(const float& pt, const float& eta, const float& phi){

}

void TreParticle::FillPtEtaPhi(CatsLorentzVector& cats_vector){

}


std::string TreParticle::GetName() const{
  return TreName;
}

float TreParticle::GetMass() const{
  return Mass[1];
}

float TreParticle::GetMassLow() const{
  return Mass[0];
}

float TreParticle::GetMassUp() const{
  return Mass[2];
}

float TreParticle::GetWidth() const{
  return Width[1];
}

float TreParticle::GetWidthLow() const{
  return Width[0];
}

float TreParticle::GetWidthUp() const{
  return Width[2];
}

float TreParticle::GetAbundance() const{
  //printf("    GA %f\n",Abundance[1]);
  return Abundance[1];
}

float TreParticle::GetAbundanceLow() const{
  return Abundance[0];
}

float TreParticle::GetAbundanceUp() const{
  return Abundance[2];
}

float TreParticle::GetRadius() const{
  //printf("    GA %f\n",Abundance[1]);
  return Radius[1];
}

float TreParticle::GetRadiusLow() const{
  return Radius[0];
}

float TreParticle::GetRadiusUp() const{
  return Radius[2];
}

float TreParticle::GetRadiusSlope() const{
  //printf("    GA %f\n",Abundance[1]);
  return RadSlope[1];
}

float TreParticle::GetRadiusSlopeLow() const{
  return RadSlope[0];
}

float TreParticle::GetRadiusSlopeUp() const{
  return RadSlope[2];
}

float TreParticle::GetDelayTau() const{
  return DelayTau[1];
}
float TreParticle::GetDelayTauLow() const{
  return DelayTau[0];
}
float TreParticle::GetDelayTauUp() const{
  return DelayTau[2];
}

unsigned char TreParticle::GetNumDecays() const{
  return NumDecays;
}

const TREPNI* TreParticle::GetDatabase() const{
  return &Database;
}

void TreParticle::AppendInBinary(std::ofstream& file){

}

void TreParticle::LoadFromBinary(std::ifstream& file){

}

void TreParticle::SetName(const char* name){
  if(strlen(name)>=Database.Len_PrtclName){
    static bool ShowMessage=true;
    printf("\033[1;33mWARNING:\033[0m (TreParticle::SetName) The name of the particle is capped at %u characters\n",
    unsigned(Database.Len_PrtclName));
    if(Database.SingleError) ShowMessage=false;
    return;
  }
  strcpy(TreName,name);
}

void TreParticle::SetMass(const float& mass){
  Mass[1] = mass;
  if(Mass[0]==Mass[2] && Mass[0]==0){
    Mass[0] = Mass[1];
    Mass[2] = Mass[1];
  }
}

void TreParticle::SetMassLimit(const float& mass_low, const float& mass_up){
  Mass[0] = mass_low;
  Mass[2] = mass_up;
  if(Mass[0]==Mass[2] && Mass[0]==0){
    Mass[0] = Mass[1];
    Mass[2] = Mass[1];
  }
}

void TreParticle::SetWidth(const float& width){
  Width[1] = width;
  if(Width[0]==Width[2] && Width[0]==0){
    Width[0] = Width[1];
    Width[2] = Width[1];
  }
}

void TreParticle::SetWidthLimit(const float& width_low, const float& width_up){
  Width[0] = width_low;
  Width[2] = width_up;
  if(Width[0]==Width[2] && Width[0]==0){
    Width[0] = Width[1];
    Width[2] = Width[1];
  }
}

void TreParticle::SetAbundance(const float& abundance){
  Abundance[1] = abundance;
  if(Abundance[0]==Abundance[2] && Abundance[0]==0){
    Abundance[0] = Abundance[1];
    Abundance[2] = Abundance[1];
  }
}

void TreParticle::SetAbundanceLimit(const float& abundance_low, const float& abundance_up){
  Abundance[0] = abundance_low;
  Abundance[2] = abundance_up;
  if(Abundance[0]==Abundance[2] && Abundance[0]==0){
    Abundance[0] = Abundance[1];
    Abundance[2] = Abundance[1];
  }
}

void TreParticle::SetRadius(const float& rad){
  //QA!!!
  Radius[1] = rad;
  if(Radius[0]==Radius[2] && Radius[0]==0){
    Radius[0] = Radius[1];
    Radius[2] = Radius[1];
  }
}
void TreParticle::SetRadiusSlope(const float& slope){
  //QA!!!
  RadSlope[1] = slope;
  if(RadSlope[0]==RadSlope[2] && RadSlope[0]==0){
    RadSlope[0] = RadSlope[1];
    RadSlope[2] = RadSlope[1];
  }
}
void TreParticle::SetDelayTau(const float& delay){
  //QA!!!
  DelayTau[1] = delay;
  if(DelayTau[0]==DelayTau[2] && DelayTau[0]==0){
    DelayTau[0] = DelayTau[1];
    DelayTau[2] = DelayTau[1];
  }
}

TreChain* TreParticle::NewDecay(){
  ResizeArray(Decay,NumDecays,NumDecays+1);
  Decay[NumDecays] = new TreChain(*this);
  return Decay[NumDecays++];
}

TreChain* TreParticle::GetDecay(const unsigned char& whichone) const{
  if(!Decay) return NULL;
  //CurrentDecay = Decay[whichone];
  return Decay[whichone];
}
TreChain* TreParticle::GetRandomDecay(DLM_Random* RanGen) const{
  if(!RanGen) RanGen = Database.RanGen;
  float rnd = RanGen->Uniform(0,100);
  float cum = 0;
  for(unsigned char uDec=0; uDec<NumDecays; uDec++){
    cum += Decay[uDec]->GetBranching();
    if(cum>=rnd) return Decay[uDec];
  }
  return NULL;
}
//TreChain* TreParticle::GetCurrentDecay(){
//  return CurrentDecay;
//}

void TreParticle::Print(){
  printf("--- Particle information ---\n");
  printf(" Name : %s\n",TreName);
  printf(" Mass : %e +(%e) -(%e)\n", Mass[1], Mass[2]-Mass[1], Mass[1]-Mass[0]);
  printf(" Width: %e +(%e) -(%e)\n", Mass[1], Width[2]-Width[1], Width[1]-Width[0]);
  printf(" Abund: %e +(%e) -(%e)\n", Abundance[1], Abundance[2]-Abundance[1], Abundance[1]-Abundance[0]);
  printf(" #decs: %u\n", NumDecays);
  for(unsigned char uDec=0; uDec<NumDecays; uDec++){
    printf("  -> ");
    for(unsigned char uDaugh=0; uDaugh<Decay[uDec]->NumDaughters; uDaugh++){
      if(uDaugh) printf(" ");
      printf("%s", Decay[uDec]->Daughter[uDaugh]->TreName);
    }
    printf("\n");
  }
}

TreChain::TreChain(TreParticle& mother):Mother(mother){
  NumDaughters = 0;
  Daughter = NULL;
  Branching = new float [3];
  Branching[0]=0;
  Branching[1]=0;
  Branching[2]=0;
  //DaughterMasses = NULL;
}
TreChain::~TreChain(){
  if(Daughter){delete[]Daughter;Daughter=NULL;NumDaughters=0;}
  //if(DaughterMasses){delete[]DaughterMasses;DaughterMasses=NULL;}
  delete [] Branching;
}

void TreChain::AddDaughter(const TreParticle& daughter){
  ResizeArray(Daughter,NumDaughters,NumDaughters+1);
  Daughter[NumDaughters] = &daughter;
  NumDaughters++;
  //if(DaughterMasses){delete[]DaughterMasses;DaughterMasses=NULL;}
}

std::string TreChain::GetName() const{
  std::string name;
  name = Mother.GetName();
  name += " -> ";
  for(unsigned char uDaugh=0; uDaugh<NumDaughters; uDaugh++){
    if(uDaugh) name += " ";
    name += Daughter[uDaugh]->GetName();
  }
  return name;
}

void TreChain::SetBranching(const float& br){
  Branching[1] = br;
  if(Branching[0]==Branching[2] && Branching[0]==0){
    Branching[0] = Branching[1];
    Branching[2] = Branching[1];
  }
}

void TreChain::SetBranchingLimit(const float& br_low, const float& br_up){
  Branching[0] = br_low;
  Branching[2] = br_up;
  if(Branching[0]==Branching[2] && Branching[0]==0){
    Branching[0] = Branching[1];
    Branching[2] = Branching[1];
  }
}

float TreChain::GetBranching() const{
  return Branching[1];
}

unsigned char TreChain::GetNumDaughters() const{
  return NumDaughters;
}
const TreParticle* TreChain::GetMother() const{
  return &Mother;
}
const TreParticle* TreChain::GetDaughter(const unsigned char& whichone) const{
  if(whichone>=NumDaughters) return NULL;
  return Daughter[whichone];
}
const double* TreChain::GetDaughterMasses() const{
  //if(!DaughterMasses){
    double* DaughterMasses = new double [NumDaughters];
    for(unsigned char uDaugh=0; uDaugh<NumDaughters; uDaugh++){
      DaughterMasses[uDaugh] = Daughter[uDaugh]->GetMass();
    }
  //}
  return DaughterMasses;
}

//void TreChain::SetDaughters(const unsigned char& numdaughters, const TreParticle* daughter){
//  if(Daughters){delete[]Daughters;Daughters=NULL;NumDaughters=0;}
//  Daughters = new TreParticle* [numdaughters];
//  NumDaughters = numdaughters;
//  for(unsigned char uch=0; uch<NumDaughters; uch++){
//    Daughters[uch] = daughter[uch];
//  }
//}

TREPNI::TREPNI(const unsigned short& version):Len_DtbsName(32),Len_PrtclName(24),
Version(version),MaxMemSteps(1024),NumFunctions(64){

  Particle = NULL;
  DatabaseName = new char [Len_DtbsName];
  NumParticles = 0;
  MaxParticles = 0;
  PrintLevel = 2;
  SingleError = true;
  ErrorOccured = new int[NumFunctions];
  for(short us=0; us<NumFunctions; us++) ErrorOccured[us] = 0;
  //QA_passed = false;
  TotAbundance = 0;
  RanGen = new DLM_Random(0);
}

TREPNI::~TREPNI(){
  delete [] DatabaseName; DatabaseName=NULL;
  delete [] ErrorOccured; ErrorOccured=NULL;
  if(Particle){
    for(unsigned uPart=0; uPart<MaxParticles; uPart++){
      if(Particle[uPart]){
        delete Particle[uPart];
        Particle[uPart] = NULL;
      }
    }
    delete [] Particle;
    Particle = NULL;
  }
  delete RanGen; RanGen=NULL;
}

bool TREPNI::QA(const int& type) const{
  bool qa = true;
  if(type<Full || type>Abundance){return false;}
  if(type==Full || type==Name){
    qa *= QA_Name();
  }
  if(type==Full || type==Daughters){
    qa *= QA_Daughters();
  }
  if(type==Full || type==Mass){
    qa *= QA_Mass();
  }
  if(type==Full || type==Width){
    qa *= QA_Width();
  }
  if(type==Full || type==BR){
    qa *= QA_BR();
  }
  if(type==Full || type==Abundance){
    qa *= QA_Abundance();
  }
  //QA_passed = qa;
  return qa;
}

bool TREPNI::QA_Name() const{
  bool AllesGut = true;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    for(unsigned uPart2=uPart+1; uPart2<NumParticles; uPart2++){
      if(strcmp(Particle[uPart]->TreName, Particle[uPart2]->TreName)==0){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) Multiple instances of particle '%s'\n",Particle[uPart]->TreName);
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
    }
    if(strcmp(Particle[uPart]->TreName,"")==0){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The name of the particle cannot be blank\n");
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    if(strncmp(Particle[uPart]->TreName,"anti_",5)==0){
      static bool ShowMessage=true;
      //printf("\033[1;35m%s %i\033[0m\n",name,int(ShowMessage));
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) Keyword 'anti_' is designated for anti-particles, "
        "which are auto-generated. Please only define particles!\n");
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    for(unsigned uChar=0; uChar<strlen(Particle[uPart]->TreName); uChar++){
      if( strncmp(&Particle[uPart]->TreName[uChar]," ",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],",",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],".",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],";",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],"\"",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],"'",1)==0 ||
          strncmp(&Particle[uPart]->TreName[uChar],"\n",1)==0
        ){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The special characters , . ; empty spaces "
          "new lines or quation marks are not allowed within the naming convention.\n");
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
    }
  }
  return AllesGut;
}

bool TREPNI::QA_Daughters() const{
  bool AllesGut = true;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    for(unsigned char uDec=0; uDec<Particle[uPart]->NumDecays; uDec++){
      if(Particle[uPart]->Decay[uDec]->NumDaughters < 2){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The particle '%s' has a decay channel with less than 2 daughters\n",
          Particle[uPart]->TreName);
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
    }
  }
  return AllesGut;
}

bool TREPNI::QA_Mass() const{
  float MotherMass;
  float MotherLow;
  float MotherUp;
  float DaughtersMass;
  bool AllesGut = true;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    MotherMass = Particle[uPart]->Mass[1];
    MotherLow = Particle[uPart]->Mass[0];
    MotherUp = Particle[uPart]->Mass[2];
    for(unsigned char uDec=0; uDec<Particle[uPart]->NumDecays; uDec++){
      DaughtersMass = 0;
      for(unsigned char uDaugh=0; uDaugh<Particle[uPart]->Decay[uDec]->NumDaughters; uDaugh++){
        DaughtersMass += Particle[uPart]->Decay[uDec]->Daughter[uDaugh]->Mass[1];
      }
      if(MotherMass<DaughtersMass){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The decay %s is impossible (Mass checksum)\n",
          Particle[uPart]->Decay[uDec]->GetName().c_str());
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
    }
    if(MotherMass<0){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) %s has negative mass\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    if(MotherLow>MotherUp || MotherLow<0){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) Invalid mass range for %s\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    if(MotherMass<MotherLow || MotherMass>MotherUp){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The %s mass is outside of the allowed range\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
  }
  return AllesGut;
}

bool TREPNI::QA_Width() const{
  float MotherWidth;
  float MotherLow;
  float MotherUp;
  bool AllesGut = true;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    MotherWidth = Particle[uPart]->Width[1];
    MotherLow = Particle[uPart]->Width[0];
    MotherUp = Particle[uPart]->Width[2];
    if(MotherWidth<0){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) %s has negative decay width\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    if(MotherLow>MotherUp || MotherLow<0){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) Invalid decay width range for %s\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
    if(MotherWidth<MotherLow || MotherWidth>MotherUp){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The %s decay width is outside of the allowed range\n",
        Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
  }
  return AllesGut;
}

bool TREPNI::QA_BR() const{
  bool AllesGut = true;
  const float TotalBR = 100;
  float MinBR = 0;
  float MaxBR = 0;
  float LowBR;
  float UpBR;
  float LowValue;
  float MeanValue;
  float UpValue;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    if(Particle[uPart]->NumDecays==0) continue;
    MinBR = 0;
    MaxBR = 0;
    for(unsigned char uDec=0; uDec<Particle[uPart]->NumDecays; uDec++){
      LowValue = Particle[uPart]->Decay[uDec]->Branching[0];
      MeanValue = Particle[uPart]->Decay[uDec]->Branching[1];
      UpValue = Particle[uPart]->Decay[uDec]->Branching[2];
      if(MeanValue<0){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) %s has negative BR\n",
          Particle[uPart]->Decay[uDec]->GetName().c_str());
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
      if(LowValue>UpValue || LowValue<0){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) Invalid BR range for %s\n",
          Particle[uPart]->Decay[uDec]->GetName().c_str());
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
      if(MeanValue<LowValue || MeanValue>UpValue){
        static bool ShowMessage=true;
        if(PrintLevel>=1 && ShowMessage){
          printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The %s BR is outside of the allowed range\n",
          Particle[uPart]->Decay[uDec]->GetName().c_str());
          if(SingleError) ShowMessage=false;
        }
        AllesGut = false;
      }
      MinBR += LowValue;
      MaxBR += UpValue;
    }
    LowBR = MinBR + 0.16*(MaxBR-MinBR);
    UpBR = MaxBR - 0.16*(MaxBR-MinBR);
    if(LowBR>TotalBR || UpBR<TotalBR){
      static bool ShowMessage=true;
      if(PrintLevel>=1 && ShowMessage){
        printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The BRs of '%s' do not sum up to 100%%\n",Particle[uPart]->TreName);
        if(SingleError) ShowMessage=false;
      }
      AllesGut = false;
    }
  }
  return AllesGut;
}

bool TREPNI::QA_Abundance() const{
  if(TotAbundance<=0){return true;}
  bool AllesGut = true;
  float MinAbund = 0;
  float MaxAbund = 0;
  float LowAbund;
  float UpAbund;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    MinAbund += Particle[uPart]->Abundance[0];
    MaxAbund += Particle[uPart]->Abundance[2];
  }
  LowAbund = MinAbund + 0.16*(MaxAbund-MinAbund);
  UpAbund = MaxAbund - 0.16*(MaxAbund-MinAbund);
  if(LowAbund>TotAbundance || UpAbund<TotAbundance){
    static bool ShowMessage=true;
    if(PrintLevel>=1 && ShowMessage){
      printf("\033[1;31mERROR:\033[0m (TREPNI::QA) The total yield of all particles is outside the allowed limit\n");
      if(SingleError) ShowMessage=false;
    }
    AllesGut = false;
  }
  return AllesGut;
}

void TREPNI::SetTotalYield(const float& totyield){
  TotAbundance = totyield;
}

float TREPNI::GetYield() const{
  if(TotAbundance>0) return TotAbundance;
  float yield = 0;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    //printf(" gy uPart=%u %p\n",uPart,Particle[uPart]);
    yield += Particle[uPart]->GetAbundance();
    //printf("  %f\n",yield);
  }
  return yield;
}

TreParticle* TREPNI::NewParticle(const char* name){
  ResizeArray(Particle,NumParticles,NumParticles+1);
  Particle[NumParticles] = new TreParticle(*this);
  if(name){Particle[NumParticles]->SetName(name);}
  return Particle[NumParticles++];
}

TreParticle* TREPNI::GetParticle(const unsigned& whichone) const{
  return Particle[whichone];
}

TreParticle* TREPNI::GetParticle(const char* name) const{
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    if(strcmp(Particle[uPart]->TreName,name)==0){
      return Particle[uPart];
    }
  }
  return NULL;
}
TreParticle* TREPNI::GetParticle(const std::string& name) const{
  return GetParticle(name.c_str());
}

//for the sampling, some node structure for log performance would be nice
TreParticle* TREPNI::GetRandomParticle(DLM_Random* rangen) const{
//printf("GetRandomParticle %p\n",Particle);
  if(!rangen) rangen = RanGen;
  const float Yield = GetYield();
//printf(" %f\n",Yield);
  float RndYield = rangen->Uniform(0,Yield);
  float Yield_Last = 0;
  float Yield_New = 0;
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
//printf("%u\n",uPart);
    Yield_New += Particle[uPart]->GetAbundance();
    if(Yield_Last<=RndYield && RndYield<=Yield_New){
      return Particle[uPart];
    }
  }
  return NULL;
}

unsigned TREPNI::GetNumParticles() const{
  return NumParticles;
}

void TREPNI::Randomize(){

}

void TREPNI::RandomizeMass(){

}

void TREPNI::RandomizeWidth(){

}

void TREPNI::RandomizeAbundance(){

}

void TREPNI::RandomizeBR(){

}


void TREPNI::SetSeed(const unsigned& seed){
  RanGen->SetSeed(seed);
}

void TREPNI::SetPrintLevel(const char& lvl, const bool& single){
  PrintLevel = lvl;
  if(PrintLevel<0) PrintLevel=0;
  if(PrintLevel>3) PrintLevel=3;
  SingleError = single;
}


/*
//FuntionID = 0
TREPNI::TREPNI(const unsigned short& version):Len_DtbsName(32),Len_PrtclName(24),
Version(version),MaxMemSteps(1024),NumFunctions(64),MaxDecayCh(16),MaxDaughters(8){
  DatabaseName = new char [Len_DtbsName];
  NumParticles = 0;
  MaxParticles = 0;
  PrintLevel = 2;
  SingleError = true;
  TrepName = NULL;
  Mass = NULL;
  Gamma = NULL;
  Nch = NULL;
  Ndaughter = NULL;
  Branching = NULL;
  DaughterID = NULL;
  ErrorOccured = new int[NumFunctions];
  for(short us=0; us<NumFunctions; us++) ErrorOccured[us] = 0;
}

//FuntionID = 1
TREPNI::~TREPNI(){
  delete [] DatabaseName;
  if(Mass){delete[]Mass;Mass=NULL;}
  if(Gamma){delete[]Gamma;Gamma=NULL;}
  if(Nch){delete[]Nch;Nch=NULL;}
  if(TrepName){
    for(unsigned uPart=0; uPart<MaxParticles; uPart++)
      if(TrepName[uPart]){delete[]TrepName[uPart];TrepName[uPart]=NULL;}
    delete[]TrepName;TrepName=NULL;
  }
  if(Ndaughter){
    for(unsigned uPart=0; uPart<MaxParticles; uPart++)
      if(Ndaughter[uPart]){delete[]Ndaughter[uPart];Ndaughter[uPart]=NULL;}
    delete[]Ndaughter;Ndaughter=NULL;
  }
  if(Branching){
    for(unsigned uPart=0; uPart<MaxParticles; uPart++)
      if(Branching[uPart]){delete[]Branching[uPart];Branching[uPart]=NULL;}
    delete[]Branching;Branching=NULL;
  }
  if(DaughterID){
    for(unsigned uPart=0; uPart<MaxParticles; uPart++){
      if(DaughterID[uPart]){
        for(char uDch=0; uDch<MaxDecayCh; uDch++){
          if(DaughterID[uPart][uDch])
          {delete[]DaughterID[uPart][uDch];DaughterID[uPart][uDch]=NULL;}
        }
        delete[]DaughterID[uPart];DaughterID[uPart]=NULL;
      }
    }
    delete[]DaughterID;DaughterID=NULL;
  }

  delete [] ErrorOccured;
}

void TREPNI::MemoryManager(const bool& destroy){
  if(NumParticles<MaxParticles){return;}

  unsigned NewSlots=0;
  unsigned TotSlots=0;
  if(!destroy){
    NewSlots = MaxParticles;
    if(NewSlots==0)NewSlots=1;
    if(NewSlots>MaxMemSteps) NewSlots=MaxMemSteps;
    TotSlots = MaxParticles+NewSlots;
  }

  ResizeArray(TrepName,MaxParticles,TotSlots);
  ResizeArray(Mass,3*MaxParticles,3*TotSlots);
  ResizeArray(Gamma,3*MaxParticles,3*TotSlots);
  ResizeArray(Nch,MaxParticles,TotSlots);
  ResizeArray(Ndaughter,MaxParticles,TotSlots);
  ResizeArray(Branching,MaxParticles,TotSlots);
  ResizeArray(DaughterID,MaxParticles,TotSlots);
  for(unsigned uPart=MaxParticles; uPart<TotSlots; uPart++){
    TrepName[uPart] = new char [24];
    Ndaughter[uPart] = new int [MaxDecayCh];
    Branching[uPart] = new float [3*MaxDecayCh];
    DaughterID[uPart] = new int* [MaxDecayCh];
    for(unsigned uDch=0; uDch<MaxDecayCh; uDch++){
      DaughterID[uPart][uDch] = new int [MaxDaughters];
    }
  }
  MaxParticles = TotSlots;
}

std::string TREPNI::GetParticleName(const int& id) const{
  int ErrorID=0;
  int LineID = fabs(id)-1;
  if(id==0){
    if(PrintLevel>=1){
      ErrorID=0;
      if( (ErrorOccured[getparticlename]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::GetParticleName) The particle numbering (ID) starts from 1\n");
        //occupy the bit in case we dont want to repeat this error
        if(SingleError){
          ErrorOccured[getparticlename] ^= (1 << ErrorID);
        }
      }
    }
    return "";
  }

  if(LineID>=NumParticles&&NumParticles){
    if(PrintLevel>=2){
      //if the corresponding bit (0 in this example)
      //is zero, than we go ahead and print out the error
      ErrorID=1;
      if( (ErrorOccured[getparticlename]&(1 << ErrorID))==0 ){
        printf("\033[1;33mWARNING:\033[0m (TREPNI::GetParticleName) There are only %u number of particles defined\n",NumParticles);
        //occupy the bit in case we dont want to repeat this error
        if(SingleError){
          ErrorOccured[getparticlename] ^= (1 << ErrorID);
        }
      }
    }
  }
  else if(LineID>=MaxParticles){
    if(PrintLevel>=1){
      ErrorID=2;
      if( (ErrorOccured[getparticlename]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::GetParticleName) There are only %u number of particles defined and %u of maximum slots\n",NumParticles,MaxParticles);
        //occupy the bit in case we dont want to repeat this error
        if(SingleError){
          ErrorOccured[getparticlename] ^= (1 << ErrorID);
        }
      }
    }
    return "";
  }
  std::string str;
  str = TrepName[LineID];
  if(id<0) str.insert(0,"anti_");
  return str;
}

int TREPNI::GetParticleId(const char* name) const{
  char* search_name = new char [29];
  int Particle = 1;
  if(strncmp(name,"anti_",5)==0){
    strcpy(search_name,&name[5]);
    Particle = -1;
  }
  else{
    strcpy(search_name,name);
  }

  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    if(strcmp(search_name,TrepName[uPart])==0){
      delete [] search_name;
      return Particle*int(uPart+1);
    }
  }

  if(PrintLevel>=2){
    int ErrorID=0;
    if( (ErrorOccured[getparticleid]&(1 << ErrorID))==0 ){
      printf("\033[1;33mWARNING:\033[0m (TREPNI::GetParticleId) The particle %s does not exist\n",name);
      //occupy the bit in case we dont want to repeat this error
      if(SingleError){
        ErrorOccured[getparticleid] ^= (1 << ErrorID);
      }
    }
  }
  delete [] search_name;
  return 0;
}


void TREPNI::SetPrintLevel(const char& lvl, const bool& single){
  PrintLevel = lvl;
  if(PrintLevel<0) PrintLevel=0;
  if(PrintLevel>3) PrintLevel=3;
  SingleError = single;
}

void TREPNI::SetParticle(const char* name, const double& mass_min, const double& mass_max,
                          const double& gamma_min, const double& gamma_max){

  int ErrorID;

  if(strcmp(name,"")==0){
    if(PrintLevel>=1){
      ErrorID=0;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The name of the particle cannot be black\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  if(strncmp(name,"anti_",5)==0){
    if(PrintLevel>=1){
      ErrorID=1;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) Keyword 'anti_' is designated for anti-particles, "
        "which are auto-generated. Please only define particles!\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  for(unsigned uChar=0; uChar<strlen(name); uChar++){
    if( strncmp(&name[uChar]," ",1)==0 ||
        strncmp(&name[uChar],",",1)==0 ||
        strncmp(&name[uChar],".",1)==0 ||
        strncmp(&name[uChar],";",1)==0 ||
        strncmp(&name[uChar],"\"",1)==0 ||
        strncmp(&name[uChar],"'",1)==0 ||
        strncmp(&name[uChar],"\n",1)==0
      ){
      if(PrintLevel>=1){
        ErrorID=2;
        if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
          printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The special characters , . ; empty spaces "
          "new lines or quation marks are not allowed within the naming convention.\n");
          if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
        }
      }
      return;
    }
  }
  if(mass_min<0||mass_min!=mass_min){
    if(PrintLevel>=1){
      ErrorID=3;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The mass of the particle is either negative or n/a\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  if(mass_max<0||mass_max!=mass_max){
    if(PrintLevel>=1){
      ErrorID=3;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The mass of the particle is either negative or n/a\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  if(mass_max<mass_min){
    if(PrintLevel>=1){
      ErrorID=4;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The lower limit of the mass is larger than the upper limit\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }

  if(gamma_min<0||gamma_min!=gamma_min){
    if(PrintLevel>=1){
      ErrorID=5;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The width of the particle is either negative or n/a\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  if(gamma_max<0||gamma_max!=gamma_max){
    if(PrintLevel>=1){
      ErrorID=5;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The width of the particle is either negative or n/a\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }
  if(gamma_max<gamma_min){
    if(PrintLevel>=1){
      ErrorID=6;
      if( (ErrorOccured[setparticle]&(1 << ErrorID))==0 ){
        printf("\033[1;31mERROR:\033[0m (TREPNI::SetParticle) The lower limit of the width is larger than the upper limit\n");
        if(SingleError) ErrorOccured[setparticle] ^= (1 << ErrorID);
      }
    }
    return;
  }

  //stop here if you find a particle with that name (rewrite it, no memory update)
  //btw, if we have already defined decay channels, we might get into a conflict by
  //changing the mass to an unrealistic value. Best implement a QA function
  //to be able to run whenever you save to a file
  for(unsigned uPart=0; uPart<NumParticles; uPart++){
    if(strcmp(name,TrepName[uPart])==0){
      Mass[3*uPart] = mass_min;
      Mass[3*uPart+1] = (mass_min+mass_max)*0.5;
      Mass[3*uPart+2] = mass_max;

      Gamma[3*uPart] = gamma_min;
      Gamma[3*uPart+1] = (gamma_min+gamma_max)*0.5;
      Gamma[3*uPart+2] = gamma_max;
      return;
    }
  }

  int id = NumParticles++;
  MemoryManager();

  strcpy(TrepName[id],name);

  Mass[3*id] = mass_min;
  Mass[3*id+1] = (mass_min+mass_max)*0.5;
  Mass[3*id+2] = mass_max;

  Gamma[3*id] = gamma_min;
  Gamma[3*id+1] = (gamma_min+gamma_max)*0.5;
  Gamma[3*id+2] = gamma_max;
}

*/
