
#include "CATSconstants.h"
#include "CECA.h"
#include "TREPNI.h"
#include "CATStools.h"
#include "DLM_Random.h"
#include "DLM_Histo.h"
#include "DLM_Source.h"
#include "DLM_CppTools.h"
#include "DLM_MathFunctions.h"

#include "omp.h"
#include <unistd.h>
#include <thread>


CecaParticle::CecaParticle(){
  cats = new CatsParticle();
  //printf("new %p %p\n",cats,this);
  trepni = NULL;
  mother = NULL;
  decay = NULL;
  Origin = 0;
}
CecaParticle::CecaParticle(const CecaParticle &other){
  //printf("copy\n");
  CecaParticle();
  *this=other;
}
CecaParticle::~CecaParticle(){
  //printf("del %p %p\n",cats,this);
  if(cats) {delete cats; cats = NULL;}
  if(mother) {delete mother; mother = NULL;}
}
const TreParticle* CecaParticle::Trepni() const{
  return trepni;
}
CatsParticle* CecaParticle::Cats() const{
  return cats;
}
const TreChain* CecaParticle::Decay() const{
  return decay;
}
CatsParticle* CecaParticle::Mother() const{
  return mother;
}
void CecaParticle::SetTrepni(const TreParticle& prt_tre){
  trepni = &prt_tre;
}
void CecaParticle::SetCats(const CatsParticle& prt_cats){
  *cats = prt_cats;
}
void CecaParticle::SetDecay(const TreChain& prt_dec){
  decay = &prt_dec;
}
void CecaParticle::SetTrepni(const TreParticle* prt_tre){
  //printf("SetTrepni(const TreParticle* prt_tre)\n");
  trepni = prt_tre;
}
void CecaParticle::SetCats(const CatsParticle* prt_cats){
  if(!cats) cats = new CatsParticle();
  *cats = *prt_cats;
}
void CecaParticle::SetDecay(const TreChain* prt_dec){
  decay = prt_dec;
}
void CecaParticle::SetMother(const CatsParticle* mama){
  if(!mother) mother = new CatsParticle();
  *mother = *mama;
}
void CecaParticle::RandomDecay(DLM_Random* RanGen){
  if(trepni){
    decay = trepni->GetRandomDecay(RanGen);
  }
}
void CecaParticle::SetOrigin(const char& origin){
  Origin = origin;
}
bool CecaParticle::IsUseful() const{
  if(Origin==1||Origin==2) return true;
  return false;
}//is of the required type
bool CecaParticle::IsUsefulPrimordial() const{
  return Origin==1;
}//+primordial
bool CecaParticle::IsUsefulProduct() const{
  return Origin==2;
}//+decay product
bool CecaParticle::WithinAcceptance() const{
  if(cats->GetPt()<trepni->AcceptanceMin_pT()) return false;
  if(cats->GetPt()>trepni->AcceptanceMax_pT()) return false;
  if(cats->GetPseudoRap()<trepni->AcceptanceMin_Eta()) return false;
  if(cats->GetPseudoRap()>trepni->AcceptanceMax_Eta()) return false;
  if(cats->GetPphi()<trepni->AcceptanceMin_Phi()) return false;
  if(cats->GetPphi()>trepni->AcceptanceMax_Phi()) return false;
  return true;
}
CecaParticle& CecaParticle::operator=(const CecaParticle& other){
  //printf("=\n");
  trepni = other.trepni;
  if(!cats&&other.cats) cats = new CatsParticle();
  if(other.cats) *cats = *other.cats;
  decay = other.decay;
  if(!mother&&other.mother) mother = new CatsParticle();
  if(other.mother) *mother = *other.mother;
  Origin = other.Origin;
  return *this;
}


//! nothing done on the errors (single, level etc), do it afterwards
CECA::CECA(const TREPNI& database,const std::vector<std::string>& list_of_particles):
//CECA::CECA(const TREPNI& database):
  Database(database),MaxThreads(std::thread::hardware_concurrency()?std::thread::hardware_concurrency():1){
  Displacement = new float [6];
  DisplacementAlpha = new float [6];
  Hadronization = new float [6];
  HadronizationAlpha = new float [6];
  for(int i=0; i<6; i++){
    Displacement[i]=0;
    DisplacementAlpha[i]=2;
    Hadronization[i]=0;
    HadronizationAlpha[i]=2;
  }
  HadrFluct = 0;
  //Hadronization = 0;
  //HadronizationAlpha = 2;
  Tau = 0;
  TauEbe = 0;
  TauFluctuation = 0;
  ProperTau = true;
  FixedHadr = true;
  FragmentBeta = 0;
  EqualFsiTau = true;
  ThermalKick = 0;
  PropagateMother = false;
  SDIM = 2;
  TargetYield = 100000;
  AchievedYield = 0;
  FemtoLimit = 200;
  UpperLimit = 300;
  EMULT = 0;
  SrcCnv = 1;
  DebugMode = false;
  exp_file_name = "";
  exp_file_flag = 0;
  //CLV.clear();

  ///////////////////////////////////////////////
  Ghetto_NumMtBins = 24;
  Ghetto_MtMin = 0;
  Ghetto_MtMax = 4096;
  Ghetto_MtBins = NULL;

  Ghetto_NumMomBins = 256;
  Ghetto_MomMin = 0;
  Ghetto_MomMax = 4096;

  Ghetto_NumRadBins = 1024;
  Ghetto_RadMin = 0;
  Ghetto_RadMax = 128;

  Ghetto_rstar = NULL;
  Ghetto_rcore = NULL;
  GhettOld_rstar = NULL;
  Old_rstar = NULL;
  Old_rcore = NULL;
  Old_source = NULL;
  SetUp_RSM = NULL;
  SetUp_RSM_UNI = NULL;
  SetUp_RSM_BB = NULL;
  Buffer_RSM = NULL;
  Old_CosRcP1 = NULL;
  Old_CosRcP2 = NULL;
  Old_CosP1P2 = NULL;
  Old_RcP1 = NULL;
  Old_RcP2 = NULL;
  Old_P1P2 = NULL;
  Ghetto_kstar = NULL;
  Ghetto_kstar_rstar = NULL;
  Ghetto_kstar_rstar_PP = NULL;
  Ghetto_kstar_rstar_PR = NULL;
  Ghetto_kstar_rstar_RP = NULL;
  Ghetto_kstar_rstar_RR = NULL;
  Ghetto_kstar_rstar_mT = NULL;
  Ghetto_kstar_rcore_mT = NULL;
  Ghetto_kstar_reso_mT = NULL;
  Ghetto_mT_rstar = NULL;
  GhettoFemto_rstar = NULL;
  GhettoFemto_rcore = NULL;
  Ghetto_mT_costh = NULL;
  GhettoSP_pT_th = NULL;
  GhettoSP_pT_1 = NULL;
  GhettoSP_pT_2 = NULL;
  GhettoSPr_X = NULL;
  GhettoSPr_Y = NULL;
  GhettoSPr_Z = NULL;
  GhettoSPr_Rho = NULL;
  GhettoSPr_R = NULL;
  GhettoSP_X = NULL;
  GhettoSP_Y = NULL;
  GhettoSP_Z = NULL;
  GhettoSP_Rho = NULL;
  GhettoSP_R = NULL;
  Ghetto_RP_AngleRcP1 = NULL;
  Ghetto_PR_AngleRcP2 = NULL;
  Ghetto_RR_AngleRcP1 = NULL;
  Ghetto_RR_AngleRcP2 = NULL;
  Ghetto_RR_AngleP1P2 = NULL;
  Ghetto_PP_AngleRcP1 = NULL;
  Ghetto_PP_AngleRcP2 = NULL;
  Ghetto_PP_AngleP1P2 = NULL;
  Ghetto_ScatteringAngle = NULL;
  GhettoFemto_mT_rstar = NULL;
  GhettoFemto_mT_rcore = NULL;
  GhettoFemto_mT_kstar = NULL;
  GhettoFemtoPrimordial_mT_kstar = NULL;
  GhettoFemto_pT1_pT2 = NULL;
  GhettoFemto_pT1_div_pT = NULL;
  Ghetto_mT_mTwrong = NULL;
  GhettoFemto_mT_mTwrong = NULL;
  GhettoPrimReso[0]=0;
  GhettoPrimReso[1]=0;
  GhettoPrimReso[2]=0;
  GhettoPrimReso[3]=0;
  GhettoFemtoPrimReso[0]=0;
  GhettoFemtoPrimReso[1]=0;
  GhettoFemtoPrimReso[2]=0;
  GhettoFemtoPrimReso[3]=0;

  ThreadClock = new DLM_Timer [MaxThreads];
  //30 seconds as a default timeout
  Timeout = 30*10000000;
  GlobalTimeout = -1;
  RanGen = new DLM_Random* [MaxThreads];
  for(unsigned uTh=0; uTh<MaxThreads; uTh++){
    RanGen[uTh] = new DLM_Random(uTh+1);
  }
  //NumSystVars = 1;
  ListOfParticles = list_of_particles;
  for(std::string& particle : ListOfParticles){
    if(!Database.GetParticle(particle)){
      printf("\033[1;31mERROR:\033[0m (CECA::CECA) The particle '%s' is not in the database\n",particle.c_str());
    }
  }

  GhettoInit();
}

CECA::~CECA(){
  delete [] Displacement;
  delete [] DisplacementAlpha;
  delete [] Hadronization;
  delete [] HadronizationAlpha;
  //CLV.clear();

  ///////////////////////////////////////////////
  if(Ghetto_rstar){delete Ghetto_rstar; Ghetto_rstar=NULL;}
  if(Ghetto_rcore){delete Ghetto_rcore; Ghetto_rcore=NULL;}
  if(GhettOld_rstar){delete GhettOld_rstar; GhettOld_rstar=NULL;}
  if(Old_rstar){delete Old_rstar; Old_rstar=NULL;}
  if(Old_rcore){delete Old_rcore; Old_rcore=NULL;}
  if(Old_CosRcP1){delete Old_CosRcP1; Old_CosRcP1=NULL;}
  if(Old_CosRcP2){delete Old_CosRcP2; Old_CosRcP2=NULL;}
  if(Old_CosP1P2){delete Old_CosP1P2; Old_CosP1P2=NULL;}
  if(Old_RcP1){delete Old_RcP1; Old_RcP1=NULL;}
  if(Old_RcP2){delete Old_RcP2; Old_RcP2=NULL;}
  if(Old_P1P2){delete Old_P1P2; Old_P1P2=NULL;}
  if(Old_source){delete Old_source; Old_source=NULL;}
  if(Ghetto_kstar){delete Ghetto_kstar; Ghetto_kstar=NULL;}
  if(Ghetto_kstar_rstar){delete Ghetto_kstar_rstar; Ghetto_kstar_rstar=NULL;}
  if(Ghetto_kstar_rstar_PP){delete Ghetto_kstar_rstar_PP; Ghetto_kstar_rstar_PP=NULL;}
  if(Ghetto_kstar_rstar_PR){delete Ghetto_kstar_rstar_PR; Ghetto_kstar_rstar_PR=NULL;}
  if(Ghetto_kstar_rstar_RP){delete Ghetto_kstar_rstar_RP; Ghetto_kstar_rstar_RP=NULL;}
  if(Ghetto_kstar_rstar_RR){delete Ghetto_kstar_rstar_RR; Ghetto_kstar_rstar_RR=NULL;}
  if(Ghetto_kstar_rstar_mT){delete Ghetto_kstar_rstar_mT; Ghetto_kstar_rstar_mT=NULL;}
  if(Ghetto_kstar_rcore_mT){delete Ghetto_kstar_rcore_mT; Ghetto_kstar_rcore_mT=NULL;}
  if(Ghetto_kstar_reso_mT){delete Ghetto_kstar_reso_mT; Ghetto_kstar_reso_mT=NULL;}
  if(Ghetto_mT_rstar){delete Ghetto_mT_rstar; Ghetto_mT_rstar=NULL;}
  if(GhettoFemto_rstar){delete GhettoFemto_rstar; GhettoFemto_rstar=NULL;}
  if(GhettoFemto_rcore){delete GhettoFemto_rcore; GhettoFemto_rcore=NULL;}
  if(Ghetto_mT_costh){delete Ghetto_mT_costh; Ghetto_mT_costh=NULL;}
  if(GhettoSP_pT_th){delete GhettoSP_pT_th; GhettoSP_pT_th=NULL;}
  if(GhettoSP_pT_1){delete GhettoSP_pT_1; GhettoSP_pT_1=NULL;}
  if(GhettoSP_pT_2){delete GhettoSP_pT_2; GhettoSP_pT_2=NULL;}
  if(GhettoSPr_X){delete GhettoSPr_X; GhettoSPr_X=NULL;}
  if(GhettoSPr_Y){delete GhettoSPr_Y; GhettoSPr_Y=NULL;}
  if(GhettoSPr_Z){delete GhettoSPr_Z; GhettoSPr_Z=NULL;}
  if(GhettoSPr_Rho){delete GhettoSPr_Rho; GhettoSPr_Rho=NULL;}
  if(GhettoSPr_R){delete GhettoSPr_R; GhettoSPr_R=NULL;}
  if(GhettoSP_X){delete GhettoSP_X; GhettoSP_X=NULL;}
  if(GhettoSP_Y){delete GhettoSP_Y; GhettoSP_Y=NULL;}
  if(GhettoSP_Z){delete GhettoSP_Z; GhettoSP_Z=NULL;}
  if(GhettoSP_Rho){delete GhettoSP_Rho; GhettoSP_Rho=NULL;}
  if(GhettoSP_R){delete GhettoSP_R; GhettoSP_R=NULL;}
  if(Ghetto_RP_AngleRcP1){delete Ghetto_RP_AngleRcP1; Ghetto_RP_AngleRcP1=NULL;}
  if(Ghetto_PR_AngleRcP2){delete Ghetto_PR_AngleRcP2; Ghetto_PR_AngleRcP2=NULL;}
  if(Ghetto_RR_AngleRcP1){delete Ghetto_RR_AngleRcP1; Ghetto_RR_AngleRcP1=NULL;}
  if(Ghetto_RR_AngleRcP2){delete Ghetto_RR_AngleRcP2; Ghetto_RR_AngleRcP2=NULL;}
  if(Ghetto_RR_AngleP1P2){delete Ghetto_RR_AngleP1P2; Ghetto_RR_AngleP1P2=NULL;}
  if(Ghetto_ScatteringAngle){delete Ghetto_ScatteringAngle; Ghetto_ScatteringAngle=NULL;}
  if(Ghetto_PP_AngleRcP1){delete Ghetto_PP_AngleRcP1; Ghetto_PP_AngleRcP1=NULL;}
  if(Ghetto_PP_AngleRcP2){delete Ghetto_PP_AngleRcP2; Ghetto_PP_AngleRcP2=NULL;}
  if(Ghetto_PP_AngleP1P2){delete Ghetto_PP_AngleP1P2; Ghetto_PP_AngleP1P2=NULL;}
  if(GhettoFemto_mT_rstar){delete GhettoFemto_mT_rstar; GhettoFemto_mT_rstar=NULL;}
  if(GhettoFemto_mT_rcore){delete GhettoFemto_mT_rcore; GhettoFemto_mT_rcore=NULL;}
  if(GhettoFemto_mT_kstar){delete GhettoFemto_mT_kstar; GhettoFemto_mT_kstar=NULL;}
  if(GhettoFemtoPrimordial_mT_kstar){delete GhettoFemtoPrimordial_mT_kstar; GhettoFemtoPrimordial_mT_kstar=NULL;}
  if(GhettoFemto_pT1_pT2){delete GhettoFemto_pT1_pT2; GhettoFemto_pT1_pT2=NULL;}
  if(GhettoFemto_pT1_div_pT){delete GhettoFemto_pT1_div_pT; GhettoFemto_pT1_div_pT=NULL;}
  if(Ghetto_mT_mTwrong){delete Ghetto_mT_mTwrong; Ghetto_mT_mTwrong=NULL;}
  if(GhettoFemto_mT_mTwrong){delete GhettoFemto_mT_mTwrong; GhettoFemto_mT_mTwrong=NULL;}
  if(Ghetto_MtBins){delete [] Ghetto_MtBins; Ghetto_MtBins=NULL;}
  if(ThreadClock){delete [] ThreadClock; ThreadClock=NULL;}
  if(RanGen){
    for(unsigned uTh=0; uTh<MaxThreads; uTh++){
      delete RanGen[uTh]; RanGen[uTh]=NULL;
    }
    delete [] RanGen; RanGen=NULL;
  }
}

void CECA::SetDisplacementX(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[0] = fabs(width);
  DisplacementAlpha[0] = levy;
}
float CECA::GetDisplacementX() const{
  return Displacement[0];
}

void CECA::SetFixedHadr(const bool& yesno){
  FixedHadr = yesno;
}
void CECA::SetFragmentBeta(const float& fragbeta){
  if(fragbeta<0||fragbeta>1){
    printf("ERROR SetFragmentBeta\n");
    return;
  }
  FragmentBeta = fragbeta;
}


void CECA::SetDisplacementY(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[1] = fabs(width);
  DisplacementAlpha[1] = levy;
}
float CECA::GetDisplacementY() const{
  return Displacement[1];
}

void CECA::SetDisplacementZ(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[2] = fabs(width);
  DisplacementAlpha[2] = levy;
}
float CECA::GetDisplacementZ() const{
  return Displacement[2];
}

//identical X,Y
void CECA::SetDisplacementT(const float& width, const float& levy){
  SetDisplacementX(width,levy);
  SetDisplacementY(width,levy);
}
float CECA::GetDisplacementT() const{
  return sqrt(Displacement[0]*Displacement[0]+Displacement[1]*Displacement[1])/sqrt(2.);
}

//identical X,Y,Z
void CECA::SetDisplacement(const float& width, const float& levy){
  SetDisplacementX(width,levy);
  SetDisplacementY(width,levy);
  SetDisplacementZ(width,levy);
}
float CECA::GetDisplacement() const{
  return sqrt(Displacement[0]*Displacement[0]+Displacement[1]*Displacement[1]+Displacement[2]*Displacement[2])/sqrt(3.);
}

void CECA::SetExportPairs(int flag, std::string file_name){
  exp_file_flag = flag;
  exp_file_name = file_name;
}
//event-by-event fluctuations of the parameters
void CECA::SetDisplacementEbeX(const float& fwidth, const float& flevy){
  Displacement[3] = fwidth;
  DisplacementAlpha[3] = flevy;
}
void CECA::SetDisplacementEbeY(const float& fwidth, const float& flevy){
  Displacement[4] = fwidth;
  DisplacementAlpha[4] = flevy;
}
void CECA::SetDisplacementEbeZ(const float& fwidth, const float& flevy){
  Displacement[5] = fwidth;
  DisplacementAlpha[5] = flevy;
}
void CECA::SetDisplacementEbeT(const float& fwidth, const float& flevy){
  SetDisplacementEbeX(fwidth,flevy);
  SetDisplacementEbeY(fwidth,flevy);
}
void CECA::SetDisplacementEbe(const float& fwidth, const float& flevy){
  SetDisplacementEbeX(fwidth,flevy);
  SetDisplacementEbeY(fwidth,flevy);
  SetDisplacementEbeZ(fwidth,flevy);
}

void CECA::SetHadronizationEbeX(const float& fwidth, const float& flevy){
  Hadronization[3] = fwidth;
  HadronizationAlpha[3] = flevy;
}
void CECA::SetHadronizationEbeY(const float& fwidth, const float& flevy){
  Hadronization[4] = fwidth;
  HadronizationAlpha[4] = flevy;
}
void CECA::SetHadronizationEbeZ(const float& fwidth, const float& flevy){
  Hadronization[5] = fwidth;
  HadronizationAlpha[5] = flevy;
}
void CECA::SetHadronizationEbeT(const float& fwidth, const float& flevy){
  SetHadronizationEbeX(fwidth,flevy);
  SetHadronizationEbeY(fwidth,flevy);  
}
void CECA::SetHadronizationEbe(const float& fwidth, const float& flevy){
  SetHadronizationEbeX(fwidth,flevy);
  SetHadronizationEbeY(fwidth,flevy);
  SetHadronizationEbeZ(fwidth,flevy);
}
void CECA::SetTauEbe(const float& fwidth){
  TauEbe = fwidth;
}


void CECA::SetHadronizationX(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[0] = width;
  HadronizationAlpha[0] = levy;
}
float CECA::GetHadronizationX() const{
  return Hadronization[0];
}

void CECA::SetHadronizationY(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[1] = width;
  HadronizationAlpha[1] = levy;
}
float CECA::GetHadronizationY() const{
  return Hadronization[1];
}

void CECA::SetHadronizationZ(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[2] = width;
  HadronizationAlpha[2] = levy;
}
float CECA::GetHadronizationZ() const{
  return Hadronization[2];
}

//identical X,Y
void CECA::SetHadronizationT(const float& width, const float& levy){
  SetHadronizationX(width,levy);
  SetHadronizationY(width,levy);
}
float CECA::GetHadronizationT() const{
  return sqrt(Hadronization[0]*Hadronization[0]+Hadronization[1]*Hadronization[1])/sqrt(2.);
}

//identical X,Y,Z
void CECA::SetHadronization(const float& width, const float& levy){
  SetHadronizationX(width,levy);
  SetHadronizationY(width,levy);
  SetHadronizationZ(width,levy);
  //if(levy<1||levy>2){
  //  printf("ERROR levy\n");
  //  return;
  //}
  //Hadronization = fabs(width);
  //HadronizationAlpha = levy;
}
float CECA::GetHadronization() const{
  return sqrt(Hadronization[0]*Hadronization[0]+Hadronization[1]*Hadronization[1]+Hadronization[2]*Hadronization[2])/sqrt(3.);
}

void CECA::SetHadrFluctuation(const float& fluct){
  HadrFluct = fluct;
}
float CECA::GetHadrFluctuation() const{
  return HadrFluct;
}
void CECA::SetPropagateMother(const bool& yesno){
  PropagateMother = yesno;
}

void CECA::SetTau(const float& tau, const bool& proper){
  if(tau<0){
    printf("ERROR tau\n");
    return;
  }
  Tau = tau;
  ProperTau = proper;
}
void CECA::SetTauFluct(const float& taufluct){
  if(taufluct>1){
    printf("ERROR taufluct\n");
    return;
  }
  TauFluctuation = taufluct;
}
float CECA::GetTau() const{
  return Tau;
}

void CECA::SetThermalKick(const float& kick){
  ThermalKick = kick;
}

void CECA::SetSourceDim(const unsigned char& sdim){
  if(sdim<2||sdim>=16){
    printf("ERROR sdim\n");
    return;
  }
  SDIM = sdim;
}

void CECA::SetTargetStatistics(const unsigned long long& yield){
  TargetYield = yield;
}

unsigned long long CECA::GetStatistics(){
  return AchievedYield;
}

void CECA::SetFemtoRegion(const float& femto, const float& info){
  if(femto<0||info<0){
    printf("ERROR femto\n");
    return;
  }
  FemtoLimit = femto;
  UpperLimit = info;
}

void CECA::SetEventMult(const unsigned short& emult){
  EMULT = emult;
}

//void CECA::SetSystVars(const unsigned& howmany){
//  if(!howmany) NumSystVars = 1;
//  else NumSystVars = howmany;
//}

void CECA::SetSourceConvention(const char& srccnv){
  if(srccnv<0||srccnv>20){
    printf("ERROR srccnv\n");
    return;
  }
  SrcCnv = srccnv;
}
void CECA::SetDebugMode(const bool& debugmode){
  DebugMode = debugmode;
}
void CECA::SetThreadTimeout(const unsigned& seconds){
  Timeout = seconds?seconds:1;//a minimum of 1 second
}
void CECA::SetGlobalTimeout(const unsigned& seconds){
  GlobalTimeout = seconds?seconds:-1;
}
void CECA::SetSeed(const unsigned& thread, const unsigned& seed){
  if(thread>=MaxThreads) return;
  RanGen[thread]->SetSeed(seed);
}
void CECA::EqualizeFsiTime(const bool& yesno){
  EqualFsiTau = yesno;
}

//returns the number of generated multiplets
unsigned CECA::GoSingleCore(const unsigned& ThId){
  ThreadClock[ThId].Start();
  unsigned ExeTime;
  unsigned NumMultiplets = 0;

  unsigned DebugCounter=0;
  if(DebugMode){
    //#pragma omp critical
    //{
    //  printf("Single core event generator (%u)\n",ThId);
    //  printf(" -> NumMult = %i; ExeTime = %u\n",NumMultiplets,ExeTime);
    //}
  }


  do{
    NumMultiplets += GenerateEvent(ThId);
    ExeTime = unsigned(ThreadClock[ThId].Stop()/(long long)(1000000));
    //printf("ExeTime = %u\n",ExeTime);
    DebugCounter++;
  }
  while(ExeTime<Timeout);
  if(DebugMode){
    //#pragma omp critical
    //{
    //  printf(" -> The exectution finished after %u calls, generating %u multiplets\n",DebugCounter,NumMultiplets);
    //  printf("\n");
    //}
  }

  return NumMultiplets;
}

void CECA::OptimizeThreadCount(){

}

void CECA::SaveBuffer(){

}

void CECA::GoBabyGo(const unsigned& num_threads){
  GhettoInit();
  if(!Database.QA()){
      return;
  }
  bool DynamicThreads;
  if(num_threads==0){
    NumThreads = std::thread::hardware_concurrency();
    if(NumThreads>MaxThreads){
      printf("\033[1;33mWARNING:\033[0m omp_get_num_threads>hardware_concurrency (%u>%u),"
      "likely a bug (not affecting physics output), please contact the developers.\n",NumThreads,MaxThreads);
      NumThreads = MaxThreads;
    }
    DynamicThreads = true;
  }
  else{
    NumThreads = num_threads;
    if(NumThreads>MaxThreads){
      printf("\033[1;33mWARNING:\033[0m CECA::GoBabyGo says you demand more threads (%u) than available on your system (%u).\n",NumThreads,MaxThreads);
      printf("  If you are sure everything is properly set up, please inform the developers.\n");
      NumThreads = MaxThreads;
    }
    DynamicThreads = false;
  }

  if(DebugMode){
    printf("Running GoBabyGo\n");
    printf(" Detected threads: %u\n",NumThreads);
  }

  if(exp_file_name!=""){
    FILE *file_ptr;
    // Open the file in write mode
    file_ptr = fopen(exp_file_name.c_str(), "w");
    if (file_ptr == NULL) {
      printf("Error opening file: %s\n", exp_file_name.c_str());
      exp_file_flag = 0;
    }
    else{
      //only rd, hT, tau
      if(exp_file_flag==2){
        fprintf(file_ptr, " kstar\t\trstar\t\tmT\t\trd\t\th\t\ttau\n");
      }
      else{
        fprintf(file_ptr, " kstar\t\trstar\t\tmT\t\trx\t\try\t\trz\t\thx\t\thy\t\thz\t\ttau\n");
      }
      
      fclose(file_ptr);
    }
  }
  else{
    exp_file_flag = 0;
  }

  unsigned* BufferYield = new unsigned [NumThreads];
  unsigned ExeTime=0;
  DLM_Timer GlobalClock;
  GlobalClock.Start();
  AchievedYield = 0;
  //const unsigned TargetPerSyst = TargetYield / NumSystVars;
  //unsigned AchievedSystYield = 0;
  //we iterate until we have our target yield
  //this refers to the global yield, but also the minimum yiled per systematics
  //the latter is done to minimize the bias within the yield of each syst. variation
  //while(AchievedYield<TargetYield || AchievedSystYield<TargetPerSyst){
  while(AchievedYield<TargetYield && ExeTime<GlobalTimeout){
    //we run each thread for a maximum of some preset amount of time
    #pragma omp parallel for
    for(unsigned uThread=0; uThread<NumThreads; uThread++){
      BufferYield[uThread] = GoSingleCore(uThread);
    }
    for(unsigned uThread=0; uThread<NumThreads; uThread++){
      AchievedYield += BufferYield[uThread];
      //AchievedSystYield += BufferYield[uThread];
    }
    ExeTime = unsigned(GlobalClock.Stop()/(long long)(1000000));
    if(DebugMode){
      char* strtime = new char [128];
      ShowTime((long long)(ExeTime),strtime,2,true,5);
      if(TargetYield<10*1000) printf(" Achieved/Target Yield = %llu / %llu @ %s\n",AchievedYield,TargetYield,strtime);
      else if(TargetYield<10*1000*1000) printf(" Achieved/Target Yield = %llu / %llu k @ %s\n",AchievedYield/1000,TargetYield/1000,strtime);
      else printf(" Achieved/Target Yield = %llu / %llu M @ %s\n",AchievedYield/1000000,TargetYield/1000000,strtime);

      delete [] strtime;
    }
    //if(AchievedSystYield>=TargetPerSyst){
    //  //Database.Randomize();
    //  AchievedSystYield = 0;
    //}
    //after the timeout, we optimize the thread count and, in case
    //TargetYield is not achieved yet, than we continue
    if(DynamicThreads) OptimizeThreadCount();
    //saves permanantly the output that was in the buffer of each thread
    SaveBuffer();
  }

  delete [] BufferYield;
}

bool CECA::ParticleInList(const std::string& name) const{
  for(std::string str : ListOfParticles){
    if(str==name)
      {return true;}
  }
  return false;
}
bool CECA::ParticleInList(const TreParticle* prt) const{
  return ParticleInList(prt->GetName());
}

unsigned CECA::GenerateEventTEMP(){
  unsigned ThId = omp_get_thread_num();
  std::vector<CecaParticle*> Primordial;
  Primordial.clear();
  for(int i=0; i<4; i++){
    for(unsigned short uMult=0; uMult<EMULT; uMult++){
      Primordial.push_back(new CecaParticle());
      //delete Primordial.back();
    }
  }

  for(CecaParticle* particle : Primordial){
    delete particle;
  }

  return 0;
}


//generates all particles, propagates and decays into the particles of interest
//and lastly builds up the
unsigned CECA::GenerateEvent(const unsigned& ThId){
    //unsigned ThId = omp_get_thread_num();

    //the event-by-event fluctuations, in absolute value
    float disp_ebe[3];
    float hadr_ebe[3];
    float tau_ebe;
    for(int xyz=0; xyz<3; xyz++){
      if(Displacement[3+xyz]){
        disp_ebe[xyz] = RanGen[ThId]->Gauss(0,Displacement[xyz]*Displacement[xyz+3]);
      }
      else{
        disp_ebe[xyz] = 0;
      }
      if(Hadronization[3+xyz]){
        hadr_ebe[xyz] = RanGen[ThId]->Gauss(0,Hadronization[xyz]*Hadronization[xyz+3]);
      }
      else{
        hadr_ebe[xyz] = 0;
      }
    }
    if(TauEbe){
      tau_ebe = RanGen[ThId]->Gauss(0,Tau*TauEbe);
    }
    else{
      tau_ebe = 0;
    }
    
    //there was some issue using the objects
    //a silly workaround: I will only keep track of pointers,
    //however we will need to call delete for each object whenever required
    std::vector<CecaParticle*> Primordial;
    std::vector<CecaParticle*> Primary;

    //while we create an event containing all final state speciies of interest
    //N.B. so far no sanity check if this is even possible, i.e. an infinite loop is more than possible here!!!
    while(true){
      for(CecaParticle* particle : Primordial) {delete particle;}
      Primordial.clear();
      for(unsigned short uMult=0; uMult<EMULT; uMult++){
        //--- SELECT A PRIMORDIAL ---//
        Primordial.push_back(new CecaParticle());

        TreParticle* tre = Database.GetRandomParticle(RanGen[ThId]);
        Primordial.back()->SetTrepni(tre);

        //check if we need this guy, i.e. is it, or one of its daughters,
        //a particle that we would like to study.
        //N.B. if it is already the particle of interest, it will NOT be decaying
        bool FSP_is_primordial = false;
        FSP_is_primordial = ParticleInList(Primordial.back()->Trepni());
        bool FSP_is_product = false;

        if(Primordial.back()->Trepni()->GetNumDecays()){
          //it should not happen, that our final state particle can decay
          if(FSP_is_primordial){
            printf("\033[1;33mWARNING:\033[0m (CECA::GenerateEvent) The final state particle '%s' has a decay channel. Ignoring!\n",Primordial.back()->Trepni()->GetName().c_str());
          }
          else{
            //we select a decay chain of this particle
            Primordial.back()->RandomDecay(RanGen[ThId]);
            for(unsigned uDaught=0; uDaught<Primordial.back()->Decay()->GetNumDaughters(); uDaught++){
              FSP_is_product = ParticleInList(Primordial.back()->Decay()->GetDaughter(uDaught));
              if(FSP_is_product) break;
            }
          }
        }
        else{
          Primordial.back()->SetDecay(NULL);
        }
        bool Useless_particle = (!FSP_is_primordial&&!FSP_is_product);

        //the particle will NOT be used
        if(Useless_particle){
          delete Primordial.back();
          Primordial.pop_back();
        }
        //here we can either get a primordial of interest, or a resonance of interest
        else{
          if(FSP_is_primordial) Primordial.back()->SetOrigin(1);
          else Primordial.back()->SetOrigin(-1);
        }

        //we end the uMult here, and will continue only after
        //we have selected all of our primordials with which to work

      }//uMult<EMULT

      //at this point, we need to check if the primordials, and their decay chains,
      //are providing the multiplets we need. If yes, we proceed with the
      //momentum generation and propagation

      //create a copy of the list of particles, and each time you find one
      //you delete it from the list. If this list is empty at the end, it means
      //we are okay, i.e. found all particles of interest
      std::vector<std::string> particle_list = ListOfParticles;

      for(int iNeeded=particle_list.size()-1; iNeeded>=0; iNeeded--){
        std::string Needed = particle_list.at(iNeeded).c_str();
        //printf("On the lookout for %s\n",Needed.c_str());
        for(const CecaParticle* particle : Primordial){
          //printf("Investigating the primoridal %s\n",particle->Trepni()->GetName().c_str());
          bool FoundIt=false;
          if(particle->Trepni()->GetName()==Needed){
            //printf(" Found a primordial %s\n",particle_list.back().c_str());
            particle_list.pop_back();
            FoundIt = true;
          }
          else if(particle->Decay()){
            unsigned char ndaugh = particle->Decay()->GetNumDaughters();
            for(unsigned char uDaugh=0; uDaugh<ndaugh; uDaugh++){
              if(particle->Decay()->GetDaughter(uDaugh)->GetName()==particle_list.at(iNeeded)){
                //printf(" Found a primary %s\n",particle_list.back().c_str());
                particle_list.pop_back();
                FoundIt = true;
                break;
              }
            }
          }
          if(FoundIt) break;
        }//particle
      }//iNeed

      //we found all particles
      if(particle_list.size()==0)
        break;

    }//the inifinite while loop


    for(CecaParticle* primordial : Primordial){
      //--- SAMPLE THE MOMENTUM ---//
      double axisValues[3];
      double& pT = axisValues[0];
      double& eta = axisValues[1];
      double& phi = axisValues[2];
      double px,py,pz,ptot,sin_th,cos_th,cotg_th,cos_phi,sin_phi;

      primordial->Trepni()->SamplePxPyPz(axisValues,RanGen[ThId],true);
      px = axisValues[0];
      py = axisValues[1];
      pz = axisValues[2];

      ptot = sqrt(px*px+py*py+pz*pz);
      cos_th = pz/ptot;
      sin_th = sqrt(1.-cos_th*cos_th);

      cos_phi = px/(ptot*sin_th);
      sin_phi = py/(ptot*sin_th);

      primordial->Cats()->SetMXYZ(primordial->Trepni()->GetMass(),px,py,pz);
      primordial->Cats()->SetWidth(primordial->Trepni()->GetWidth());
      primordial->Cats()->SetDecayRanGen(RanGen[ThId]);

      //--- EMISSION ---//
      //--- PROPAGATE BASED ON THE PROPERTIES OF THE CORE SOURCE ---//
      double rd[3],beta[3],rtot[3],rh[3],mom[3],mom0[3];
//printf("0R = %.2f\n",sqrt(rtot[0]*rtot[0]+rtot[1]*rtot[1]+rtot[2]*rtot[2]));
      double energy;
      double tau = Tau;
      //the model where we assume the source is an ellipsoid around the displacement point,
      //and that direction of velocity is what determines the "crossing point" of the particle with the emission source
      double rh_len=-1;
      int ResampleCount = 0;
      const int ResampleLimit = 1000;

      double RejectProb = 0;
      double dr_tot;
      //save the original momentum (used when resampling)
      for(int xyz=0; xyz<3; xyz++) mom0[xyz] = primordial->Cats()->GetP(xyz);

      while(true){
        for(int xyz=0; xyz<3; xyz++) primordial->Cats()->SetMomXYZ(mom0[0],mom0[1],mom0[2]);
        if(ResampleCount==0){
          for(int xyz=0; xyz<3; xyz++){
            rd[xyz] = RanGen[ThId]->Gauss(0,(Displacement[xyz]+disp_ebe[xyz]));
            //rh is the extension of the ellipsoid in each direction
            rh[xyz] = RanGen[ThId]->Gauss((Hadronization[xyz]+hadr_ebe[xyz]),(Hadronization[xyz]+hadr_ebe[xyz])*HadrFluct);
          }
          //this comes from the definition of an ellipsoid, where each term inside is x,y,z coordinate evaluated
          //from the angles and spacial extension of the ellipsoid
          rh_len = sqrt(pow(rh[0]*sin_th*cos_phi,2.)+pow(rh[1]*sin_th*sin_phi,2.)+pow(rh[2]*cos_th,2.));
          if(TauFluctuation<0) tau = RanGen[ThId]->Exp(fabs(Tau+tau_ebe));
          else tau = (Tau+tau_ebe)+RanGen[ThId]->Gauss(0,(Tau+tau_ebe)*TauFluctuation);
          if(ProperTau) tau *= primordial->Cats()->Gamma();
        }
        //at resampling, we do not change momentum or hadronization, or time component, we just shift a little bit
        //the spacial vector, so that we have a higher chance not to be on top of another particle.
        //the resampling is a fluctuation around the original values, scaled by the size of the hadron and
        //the probability with which it has been rejected (i.e. resampling is largest for particles that really are on top of one another)
        else{
          for(int xyz=0; xyz<3; xyz++){
            rd[xyz] += RanGen[ThId]->Gauss(0,RejectProb*primordial->Trepni()->GetRadius());
          }
        }
/*
double MASS = primordial->Cats()->Mag();
double FragCorr;
if(MASS<400) FragCorr = 1;
else if(MASS<800) FragCorr = 0.8;
else if(MASS<1600) FragCorr = 0.6;
else FragCorr = 0.4;
FragCorr = 1;
*/
        //add the displacement and the beta*tau components
        energy=primordial->Cats()->Mag2();
        for(int xyz=0; xyz<3; xyz++){
          //thermal kick
          //N.B. because of it, we need to reevaluate beta mom etc of the particle and
          //we cannot use the primordial->Cats() !!!
          mom[xyz] = RanGen[ThId]->Gauss(mom0[xyz],ThermalKick);
          energy += mom[xyz]*mom[xyz];
        }
        energy = sqrt(energy);
        double mom_tot=0;

        for(int xyz=0; xyz<3; xyz++) beta[xyz] = mom[xyz]/energy;
        double beta_tot = sqrt(beta[0]*beta[0]+beta[1]*beta[1]+beta[2]*beta[2]);
        for(int xyz=0; xyz<3; xyz++){
          rtot[xyz] = rd[xyz];
          if(FragmentBeta){
            if(beta_tot) rtot[xyz] += beta[xyz]*tau*(FragmentBeta/beta_tot);
          }
          else rtot[xyz] += beta[xyz]*tau;
          mom_tot += mom[xyz]*mom[xyz];
        }
//printf("gR = %.2f\n",sqrt(rd[0]*rd[0]+rd[1]*rd[1]+rd[2]*rd[2]));
//printf("tR = %.2f\n",sqrt(rtot[0]*rtot[0]+rtot[1]*rtot[1]+rtot[2]*rtot[2]));
        mom_tot = sqrt(mom_tot);
        //we need to add the hadronization part separately, as we demand it to have
        //the same direction as the velocity, i.e. we need beta first
        for(int xyz=0; xyz<3; xyz++){
          if(FragmentBeta){
            if(beta_tot){
              rtot[xyz]+=beta[xyz]*(FragmentBeta/beta_tot)*rh_len;
            }
          }
          else{
            rtot[xyz]+=beta[xyz]/beta_tot*rh_len;
          }
        }
//printf("hR = %.2f\n",sqrt(rtot[0]*rtot[0]+rtot[1]*rtot[1]+rtot[2]*rtot[2]));
        //in the last step, particles with delayed time of formation are set to be produced with
        //a time offset. This offset is concidered to be given as proper time, and the particle
        //is simply propagated in a straight line
        if(primordial->Trepni()->GetDelayTau()){
          double gamma = energy/primordial->Cats()->Mag();
          //the time traveled evaluated in LAB
          double dtau = gamma*primordial->Trepni()->GetDelayTau();
          //printf("delayed by %f\n",dtau);
          tau += dtau;
          for(int xyz=0; xyz<3; xyz++){
            rtot[xyz] += beta[xyz]*dtau;
          }
        }
//printf("dR = %.2f\n",sqrt(rtot[0]*rtot[0]+rtot[1]*rtot[1]+rtot[2]*rtot[2]));
        //we update tau based on the time it took the partile to travel rh_len
        if(!FixedHadr){
          if(FragmentBeta) tau += rh_len*(FragmentBeta/beta_tot);
          else tau += rh_len;
        }


        //the final position is saved. The time corresponds to the time elapsed
        //in the laboratory
        primordial->Cats()->SetTXYZ(tau,rtot[0],rtot[1],rtot[2]);
        primordial->Cats()->SetMomXYZ(mom[0],mom[1],mom[2]);

//printf(" R = %.2f\n",sqrt(rtot[0]*rtot[0]+rtot[1]*rtot[1]+rtot[2]*rtot[2]));
//usleep(100e3);

        double p_tot,p_x,p_y,p_z;
        double dr_x,dr_y,dr_z;
        double LorentzWeight[2];
        double LorentzSize[2];
        double Size[2];
        double Slope[2];
        CecaParticle* prim[2];
        prim[0] = primordial;
        //the probability to accept this position sampling
        RejectProb = 0;
        //iterate over all particles, to see if they overalap. If need, resample
        for(CecaParticle* primordial2 : Primordial){
          prim[1] = primordial2;
          //this break statement makes sure we only concider 12 combinations, not 21
          //moreover, this here is needed as we are inside the primoridal loop, and the entries
          //after the current primordial are still empty.
          if(prim[0]==prim[1]) break;

          //this is needed to avoid having RejectProb==1 for the first time we iterate
          if(!RejectProb) RejectProb=1;
          for(int ip=0; ip<2; ip++){
            //we need the effective distance between the two particles, that is modulated based on the lorentz contaction
            //to do that, we find the projection of the unity vector of the primordial onto the radius vector connecting
            //the two particles. The length of the projection is the weight with which the contracted radius is taken
            Size[ip] = prim[ip]->Trepni()->GetRadius();
            Slope[ip] = prim[ip]->Trepni()->GetRadiusSlope();
            LorentzSize[ip] = Size[ip]/prim[ip]->Cats()->Gamma();
            double ls2 = LorentzSize[ip]*LorentzSize[ip];
            double sz2 = Size[ip]*Size[ip];

            if(ls2&&sz2){
              p_tot = prim[ip]->Cats()->GetP();
              p_x = prim[ip]->Cats()->GetPx();
              p_y = prim[ip]->Cats()->GetPy();
              p_z = prim[ip]->Cats()->GetPz();

              dr_x = prim[ip]->Cats()->GetX()-prim[!ip]->Cats()->GetX();
              dr_y = prim[ip]->Cats()->GetY()-prim[!ip]->Cats()->GetY();
              dr_z = prim[ip]->Cats()->GetZ()-prim[!ip]->Cats()->GetZ();
              dr_tot = sqrt(dr_x*dr_x+dr_y*dr_y+dr_z*dr_z);

              //this the projection of the unity vector. The denum. takes care of this unity.
              LorentzWeight[ip] = fabs(p_x*dr_x+p_y*dr_y+p_z*dr_z)/(p_tot*dr_tot);
              double lw2 = LorentzWeight[ip]*LorentzWeight[ip];
              double EffectiveSize = sqrt(lw2*ls2+(1.-lw2)*sz2);
              double FD = exp((dr_tot-EffectiveSize)/(EffectiveSize*Slope[ip]))+1.;
              RejectProb /= FD?FD:1;
            }
            else RejectProb=0;

          }//ip
        }//iter over primoridal2

        if(RanGen[ThId]->Uniform(0,1)>=RejectProb){
          break;
        }

        if(ResampleCount>=ResampleLimit){
          printf("\033[1;33mWARNING:\033[0m CECA::GenerateEvent says that it cannot separate the particles at the set requirement.\n");
          printf("   RejectProb = %e\n",RejectProb);
          printf("   To solve the issue, verify there is enought displacement and the particle radius is not too large.\n");
          break;
        }
        ResampleCount++;
      }//while(rh_len<0)

      //--- DECAY OF RESONANCES + PROPAGATION OF THE MOTHERS ---/
      //if the width is zero, the Decay function returns the daughters
      bool IsResonance = true;
      IsResonance = !ParticleInList(primordial->Trepni());

      bool PrimoridalUsed = false;

      if(IsResonance){
        CatsParticle* daughters =
          primordial->Cats()->Decay(
          primordial->Decay()->GetDaughterMasses(),
          PropagateMother);

        for(unsigned char nd=0; nd<primordial->Decay()->GetNumDaughters(); nd++){
          if(ParticleInList(primordial->Decay()->GetDaughter(nd))){
            Primary.push_back(new CecaParticle());
            Primary.back()->SetTrepni(primordial->Decay()->GetDaughter(nd));
            Primary.back()->SetDecay(NULL);
            Primary.back()->SetCats(daughters[nd]);
            Primary.back()->SetOrigin(2);
            Primary.back()->SetMother(primordial->Cats());

            if(!Primary.back()->WithinAcceptance()){
              delete Primary.back();
              Primary.pop_back();
            }
            else{
              PrimoridalUsed = true;
              GhettoSP_pT_th->Fill(Primary.back()->Cats()->GetPt(),Primary.back()->Cats()->GetTheta());
              if(Primary.back()->Trepni()->GetName()==ListOfParticles.at(0)){
                GhettoSP_pT_1->Fill(Primary.back()->Cats()->GetPt());
              }
              if(Primary.back()->Trepni()->GetName()==ListOfParticles.at(1)){
                GhettoSP_pT_2->Fill(Primary.back()->Cats()->GetPt());
              }
            }
          }
        }
        delete [] daughters;
      }
      else{
        Primary.push_back(new CecaParticle());
        *Primary.back() = *primordial;
        if(!Primary.back()->WithinAcceptance()){
          delete Primary.back();
          Primary.pop_back();
        }
        else{
          PrimoridalUsed = true;
          GhettoSP_pT_th->Fill(Primary.back()->Cats()->GetPt(),Primary.back()->Cats()->GetTheta());
          if(Primary.back()->Trepni()->GetName()==ListOfParticles.at(0)){
            GhettoSP_pT_1->Fill(Primary.back()->Cats()->GetPt());
          }
          if(Primary.back()->Trepni()->GetName()==ListOfParticles.at(1)){
            GhettoSP_pT_2->Fill(Primary.back()->Cats()->GetPt());
          }
        }
      }

      if(PrimoridalUsed){
        double crd_x,crd_y,crd_z;
        crd_x = primordial->Cats()->GetX();
        crd_y = primordial->Cats()->GetY();
        crd_z = primordial->Cats()->GetZ();
        GhettoSPr_X->Fill(crd_x);
        GhettoSPr_Y->Fill(crd_y);
        GhettoSPr_Z->Fill(crd_z);
        GhettoSPr_Rho->Fill(sqrt(crd_x*crd_x+crd_y*crd_y));
        GhettoSPr_R->Fill(sqrt(crd_x*crd_x+crd_y*crd_y+crd_z*crd_z));
      }

    }//iteration over all primordials


    //BUILD THE MULTIPLETS, EVALUATE THEIR NUMBER AND RETURN THE CORRECT VALUE
    //the array position of each particle, which is to be used to build the multiplet
    //the length SDIM represents the number of particles in each multiplet
    //e.g. 5 particles, SDIM=3 has to build all permutations: 012,013,014,023,024,034,123,124,134,234
    std::vector<std::vector<unsigned>> Permutations = BinomialPermutations(Primary.size(),SDIM);
    //the pid is a single permutation, e.g. {0,1,2}
    unsigned FemtoPermutations = 0;
    for(std::vector<unsigned> pid : Permutations){

/////////////////////////// TESTO TO INCREASE STAT, MAYBE WRONG!!!
      /*
      CecaParticle* prt_tmp = new CecaParticle[SDIM];
      unsigned char ud=0;
      for(unsigned ID : pid){
          prt_tmp[ud] = *Primary.at(ID);
          ud++;
      }
      CatsLorentzVector cm_rel_tmp = *prt_tmp[1].Cats()-*prt_tmp[0].Cats();
      delete [] prt_tmp;
      if(cm_rel_tmp.GetP()>FemtoLimit*4) continue;
      */
/////////////////////////////


//GHETTO: make multiplets and simply drop the output for the source as a function of rstar, no kstar, no shit
//this so that you can show something next FemTUM
      //these are all particles we need to include in a multiplet
      CatsLorentzVector boost_v;
      CecaParticle* prt_cm = new CecaParticle[SDIM];
      CecaParticle* prt_lab = new CecaParticle[SDIM];
      unsigned char ud=0;
      std::vector<float> cos_th;
      for(unsigned ID : pid){
        boost_v = boost_v+*(Primary.at(ID)->Cats());//
        prt_cm[ud] = *Primary.at(ID);
        prt_lab[ud] = *Primary.at(ID); 
        //prt_cm[ud].SetCats(Primary.at(ID)->Cats());
        //is_promordial[ud] = Primary.at(ID)->IsUsefulPrimordial();
        if(Primary.at(ID)->IsUseful()==false){printf("How did this happen!?!?!\n");}
        cos_th.push_back(cos(prt_cm[ud].Cats()->GetTheta()));
        ud++;
      }

      //the starting time of the interaction
      //given by the last particle to form
      double fsi_tau=-1e64;

      for(unsigned char ud=0; ud<SDIM; ud++){
        prt_cm[ud].Cats()->Boost(boost_v);
        if(prt_cm[ud].Mother()) prt_cm[ud].Mother()->Boost(boost_v);
        if(prt_cm[ud].Cats()->GetT()>fsi_tau){
          fsi_tau = prt_cm[ud].Cats()->GetT();
        }
      }

      //unify the time of all particles. I.e. in their rest frame, tau should be the same
      //if not, the particles that are formed earlier are propagated along a straight line until
      //the time of formation of the last particle is reached
      if(EqualFsiTau){
        for(unsigned char ud=0; ud<SDIM; ud++){
          double dtau = fsi_tau-prt_cm[ud].Cats()->GetT();
          if(dtau<0){
            printf("How did this dtau happen!?!?!\n");
            continue;
          }
          if(fabs(dtau/fsi_tau)>1e-12){
            //we propagate the particle in THIS frame of reference by dtau
            prt_cm[ud].Cats()->Propagate(dtau,false);
          }
        }
      }

      CatsLorentzVector cm_sumQA;
      for(unsigned char ud=0; ud<SDIM; ud++){
        cm_sumQA = cm_sumQA+*prt_cm[ud].Cats();
      }


      for(CecaParticle* primary : Primary){
        double crd_x,crd_y,crd_z;
        crd_x = primary->Cats()->GetX();
        crd_y = primary->Cats()->GetY();
        crd_z = primary->Cats()->GetZ();
        GhettoSP_X->Fill(crd_x);
        GhettoSP_Y->Fill(crd_y);
        GhettoSP_Z->Fill(crd_z);
        GhettoSP_Rho->Fill(sqrt(crd_x*crd_x+crd_y*crd_y));
        GhettoSP_R->Fill(sqrt(crd_x*crd_x+crd_y*crd_y+crd_z*crd_z));

      }
//NEXT_STEPS
//the tau correction, based on largest tau, and than build up R and Q, plot R for Q<FemtoLimit.
//test for two particles first!!!
//so plot rstar for femto pairs and see how it looks
//also plot the angles relevant for epos comparison

#pragma omp critical
{
//GHETTO, works for pairs only
if(SDIM==2){
CatsLorentzVector cm_rel = *prt_cm[1].Cats()-*prt_cm[0].Cats();
CatsLorentzVector cm_core;

if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
  cm_core = *prt_cm[1].Cats()-*prt_cm[0].Cats();
}
if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulProduct()){
  cm_core = *prt_cm[1].Mother()-*prt_cm[0].Cats();
}
if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulPrimordial()){
  cm_core = *prt_cm[1].Cats()-*prt_cm[0].Mother();
}
if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulProduct()){
  cm_core = *prt_cm[1].Mother()-*prt_cm[0].Mother();
}

double kstar = 0.5*cm_rel.GetP();
double rstar = cm_rel.GetR();
double rcore = cm_core.GetR();
double klab = 0.5*boost_v.GetP();
double kT = 0.5*boost_v.GetPt();
double m1 = prt_cm[0].Cats()->GetMass();
double m2 = prt_cm[1].Cats()->GetMass();
double mavg = (m1+m2)*0.5;
double mT = 0.5*boost_v.GetMt();
double mT_wrong = sqrt(kT*kT+mavg*mavg);


double AngleP1P2=0;
double AngleRcP1=0;
double AngleRcP2=0;
double BGT,BGT1,BGT2;
if(kstar<FemtoLimit){

  if(exp_file_flag){
    FILE *file_ptr;
    // Open the file in append mode
    file_ptr = fopen(exp_file_name.c_str(), "a");
    if(exp_file_flag==2){
      fprintf(file_ptr, "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", kstar, rstar, mT, 
              sqrt(pow(Displacement[0]+disp_ebe[0],2.)+pow(Displacement[1]+disp_ebe[1],2.)+pow(Displacement[2]+disp_ebe[2],2.))/sqrt(3.),
              sqrt(pow(Hadronization[0]+hadr_ebe[0],2.)+pow(Hadronization[1]+hadr_ebe[1],2.))/sqrt(2.),
              Tau+tau_ebe);
    }
    else{
      fprintf(file_ptr, "%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", kstar, rstar, mT, 
              Displacement[0], Displacement[1], Displacement[2],
              Hadronization[0], Hadronization[1], Hadronization[2],
              Tau);
    }

    fclose(file_ptr);
  }

  if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
    double cosine = cm_core.GetCosAngleRP(prt_cm[0].Cats());
    AngleRcP1 = acos(cosine);
    cosine = cm_core.GetCosAngleRP(prt_cm[1].Cats());
    AngleRcP2 = acos(cosine);
    cosine = prt_cm[0].Cats()->GetCosAngleP(prt_cm[1].Cats());
    if(cosine>1 || cosine<-1) cosine = round(cosine);
    AngleP1P2 = acos(cosine);
    if( prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)&&
        prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1)){
          Ghetto_PP_AngleRcP1->Fill(AngleRcP1);
          Ghetto_PP_AngleRcP2->Fill(AngleRcP2);
          Ghetto_PP_AngleP1P2->Fill(AngleP1P2);
    }
    if( prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1)&&
        prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0)){
          Ghetto_PP_AngleRcP1->Fill(Pi-AngleRcP2);
          Ghetto_PP_AngleRcP2->Fill(Pi-AngleRcP1);
          Ghetto_PP_AngleP1P2->Fill(AngleP1P2);
    }
  }
  else if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulPrimordial()){
    double cosine = cm_core.GetCosAngleRP(prt_cm[0].Mother());
    AngleRcP1 = acos(cosine);
    if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)){
      Ghetto_RP_AngleRcP1->Fill(AngleRcP1);
      if(SetUp_RSM){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM->AddBGT_RP(BGT,cos(AngleRcP1));
      }
      if(SetUp_RSM_UNI){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM_UNI->AddBGT_RP(BGT,RanGen[ThId]->Uniform(-1,1));
      }
      if(SetUp_RSM_BB){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM_BB->AddBGT_RP(BGT,-1.);
      }
      if(Buffer_RSM){
        Buffer_RSM->push_back(new float[11]);
        Buffer_RSM->back()[0]=110;
        Buffer_RSM->back()[1]=kstar*2;
        Buffer_RSM->back()[2]=prt_cm[0].Mother()->GetP();
        Buffer_RSM->back()[3]=0;
        Buffer_RSM->back()[4]=prt_cm[0].Mother()->GetMass();
        Buffer_RSM->back()[5]=0;
        Buffer_RSM->back()[6]=hbarc/prt_cm[0].Mother()->GetWidth();
        Buffer_RSM->back()[7]=0;
        Buffer_RSM->back()[8]=AngleRcP1;
        Buffer_RSM->back()[9]=0;
        Buffer_RSM->back()[10]=0;
      }
    }
    if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1)){
      Ghetto_PR_AngleRcP2->Fill(Pi-AngleRcP1);
      if(SetUp_RSM){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM->AddBGT_PR(BGT,cos(Pi-AngleRcP1));
      }
      if(SetUp_RSM_UNI){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM_UNI->AddBGT_PR(BGT,RanGen[ThId]->Uniform(-1,1));
      }
      if(SetUp_RSM_BB){
        BGT = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
        SetUp_RSM_BB->AddBGT_PR(BGT,1.);
      }
      if(Buffer_RSM){
        Buffer_RSM->push_back(new float[11]);
        Buffer_RSM->back()[0]=101;
        Buffer_RSM->back()[1]=kstar*2;
        Buffer_RSM->back()[2]=0;
        Buffer_RSM->back()[3]=prt_cm[0].Mother()->GetP();
        Buffer_RSM->back()[4]=0;
        Buffer_RSM->back()[5]=prt_cm[0].Mother()->GetMass();
        Buffer_RSM->back()[6]=0;
        Buffer_RSM->back()[7]=hbarc/prt_cm[0].Mother()->GetWidth();
        Buffer_RSM->back()[8]=0;
        Buffer_RSM->back()[9]=Pi-AngleRcP1;
        Buffer_RSM->back()[10]=0;
      }
    }
  }
  else if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulProduct()){
    double cosine = cm_core.GetCosAngleRP(prt_cm[1].Mother());
    AngleRcP2 = acos(cosine);
    if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1)){
      Ghetto_PR_AngleRcP2->Fill(AngleRcP2);
      if(SetUp_RSM){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM->AddBGT_PR(BGT,cos(AngleRcP2));
      }
      if(SetUp_RSM_UNI){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM_UNI->AddBGT_PR(BGT,RanGen[ThId]->Uniform(-1,1));
      }
      if(SetUp_RSM_BB){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM_BB->AddBGT_PR(BGT,1.);
      }
      if(Buffer_RSM){
        Buffer_RSM->push_back(new float[11]);
        Buffer_RSM->back()[0]=101;
        Buffer_RSM->back()[1]=kstar*2;
        Buffer_RSM->back()[2]=0;
        Buffer_RSM->back()[3]=prt_cm[1].Mother()->GetP();
        Buffer_RSM->back()[4]=0;
        Buffer_RSM->back()[5]=prt_cm[1].Mother()->GetMass();
        Buffer_RSM->back()[6]=0;
        Buffer_RSM->back()[7]=hbarc/prt_cm[1].Mother()->GetWidth();
        Buffer_RSM->back()[8]=0;
        Buffer_RSM->back()[9]=AngleRcP2;
        Buffer_RSM->back()[10]=0;
      }
    }
    if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0)){
      Ghetto_RP_AngleRcP1->Fill(Pi-AngleRcP2);
      if(SetUp_RSM){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM->AddBGT_RP(BGT,cos(Pi-AngleRcP2));
      }
      if(SetUp_RSM_UNI){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM_UNI->AddBGT_RP(BGT,RanGen[ThId]->Uniform(-1,1));
      }
      if(SetUp_RSM_BB){
        BGT = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
        SetUp_RSM_BB->AddBGT_RP(BGT,-1.);
      }
      if(Buffer_RSM){
        Buffer_RSM->push_back(new float[11]);
        Buffer_RSM->back()[0]=110;
        Buffer_RSM->back()[1]=kstar*2;
        Buffer_RSM->back()[2]=prt_cm[1].Mother()->GetP();
        Buffer_RSM->back()[3]=0;
        Buffer_RSM->back()[4]=prt_cm[1].Mother()->GetMass();
        Buffer_RSM->back()[5]=0;
        Buffer_RSM->back()[6]=hbarc/prt_cm[1].Mother()->GetWidth();
        Buffer_RSM->back()[7]=0;
        Buffer_RSM->back()[8]=Pi-AngleRcP2;
        Buffer_RSM->back()[9]=0;
        Buffer_RSM->back()[10]=0;
      }
    }
  }
  else if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulProduct()){
    double cosine = cm_core.GetCosAngleRP(prt_cm[0].Mother());
    AngleRcP1 = acos(cosine);

    cosine = cm_core.GetCosAngleRP(prt_cm[1].Mother());
    AngleRcP2 = acos(cosine);

    cosine = prt_cm[0].Mother()->GetCosAngleP(prt_cm[1].Mother());
    AngleP1P2 = acos(cosine);
    if( prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)&&
        prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1)){
          Ghetto_RR_AngleRcP1->Fill(AngleRcP1);
          Ghetto_RR_AngleRcP2->Fill(AngleRcP2);
          Ghetto_RR_AngleP1P2->Fill(AngleP1P2);
          if(SetUp_RSM){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            SetUp_RSM->AddBGT_RR(BGT1,cos(AngleRcP1),BGT2,cos(AngleRcP2),cos(AngleP1P2));
          }
          if(SetUp_RSM_UNI){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            CatsLorentzVector clvp1;
            CatsLorentzVector clvp2;
            CatsLorentzVector clvrc;
            clvrc.SetTXYZ(0,1,0,0);
            clvp1.SetTXYZ(0,RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1));
            clvp2.SetTXYZ(0,RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1));
            SetUp_RSM->AddBGT_RR(BGT1,clvrc.GetCosAngleR(clvp1),BGT2,clvrc.GetCosAngleR(clvp2),clvp1.GetCosAngleR(clvp2));
          }
          if(SetUp_RSM_BB){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            SetUp_RSM->AddBGT_RR(BGT1,-1.0,BGT2,1.0,-1.0);
          }
          if(Buffer_RSM){
            Buffer_RSM->push_back(new float[11]);
            Buffer_RSM->back()[0]=111;
            Buffer_RSM->back()[1]=kstar*2;
            Buffer_RSM->back()[2]=prt_cm[0].Mother()->GetP();
            Buffer_RSM->back()[3]=prt_cm[1].Mother()->GetP();
            Buffer_RSM->back()[4]=prt_cm[0].Mother()->GetMass();
            Buffer_RSM->back()[5]=prt_cm[1].Mother()->GetMass();
            Buffer_RSM->back()[6]=hbarc/prt_cm[0].Mother()->GetWidth();
            Buffer_RSM->back()[7]=hbarc/prt_cm[1].Mother()->GetWidth();
            Buffer_RSM->back()[8]=AngleRcP1;
            Buffer_RSM->back()[9]=AngleRcP2;
            Buffer_RSM->back()[10]=AngleP1P2;
          }
    }
    if( prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1)&&
        prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0)){
          Ghetto_RR_AngleRcP1->Fill(Pi-AngleRcP2);
          Ghetto_RR_AngleRcP2->Fill(Pi-AngleRcP1);
          //printf("Pi-RcP2 is core to %s\n",prt_cm[0].Trepni()->GetName().c_str());
          //usleep(250e3);
          Ghetto_RR_AngleP1P2->Fill(AngleP1P2);
          if(SetUp_RSM){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            SetUp_RSM->AddBGT_RR(BGT1,cos(Pi-AngleRcP2),BGT2,cos(Pi-AngleRcP1),cos(AngleP1P2));
          }
          if(SetUp_RSM_UNI){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            CatsLorentzVector clvp1;
            CatsLorentzVector clvp2;
            CatsLorentzVector clvrc;
            clvrc.SetTXYZ(0,1,0,0);
            clvp1.SetTXYZ(0,RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1));
            clvp2.SetTXYZ(0,RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1),RanGen[ThId]->Gauss(0,1));
            SetUp_RSM->AddBGT_RR(BGT1,clvrc.GetCosAngleR(clvp1),BGT2,clvrc.GetCosAngleR(clvp2),clvp1.GetCosAngleR(clvp2));
          }
          if(SetUp_RSM_BB){
            BGT1 = RanGen[ThId]->Exponential(prt_cm[0].Mother()->GetWidth()*prt_cm[0].Mother()->GetMass()/prt_cm[0].Mother()->GetP())*hbarc;
            BGT2 = RanGen[ThId]->Exponential(prt_cm[1].Mother()->GetWidth()*prt_cm[1].Mother()->GetMass()/prt_cm[1].Mother()->GetP())*hbarc;
            SetUp_RSM->AddBGT_RR(BGT1,-1.0,BGT2,1.0,-1.0);
          }

          if(Buffer_RSM){
            Buffer_RSM->push_back(new float[11]);
            Buffer_RSM->back()[0]=111;
            Buffer_RSM->back()[1]=kstar*2;
            Buffer_RSM->back()[2]=prt_cm[1].Mother()->GetP();
            Buffer_RSM->back()[3]=prt_cm[0].Mother()->GetP();
            Buffer_RSM->back()[4]=prt_cm[1].Mother()->GetMass();
            Buffer_RSM->back()[5]=prt_cm[0].Mother()->GetMass();
            Buffer_RSM->back()[6]=hbarc/prt_cm[1].Mother()->GetWidth();
            Buffer_RSM->back()[7]=hbarc/prt_cm[0].Mother()->GetWidth();
            Buffer_RSM->back()[8]=Pi-AngleRcP2;
            Buffer_RSM->back()[9]=Pi-AngleRcP1;
            Buffer_RSM->back()[10]=AngleP1P2;
          }
    }
  }
Ghetto_ScatteringAngle->Fill(cm_rel.GetScatAngle());
GhettoFemto_mT_mTwrong->Fill(mT,mT_wrong);
//the heavier particle is on x
if(prt_lab[0].Cats()->GetMass()>prt_lab[1].Cats()->GetMass()){
  GhettoFemto_pT1_pT2->Fill(prt_lab[0].Cats()->GetPt(),prt_lab[1].Cats()->GetPt());
  GhettoFemto_pT1_div_pT->Fill(prt_lab[0].Cats()->GetPt()/boost_v.GetPt());
}
else{
  GhettoFemto_pT1_pT2->Fill(prt_lab[1].Cats()->GetPt(),prt_lab[0].Cats()->GetPt());
  GhettoFemto_pT1_div_pT->Fill(prt_lab[1].Cats()->GetPt()/boost_v.GetPt());
}
//printf("%f, %f\n",prt_cm[0].Cats()->GetPt()/prt_cm[1].Cats()->GetPt(),prt_cm[0].Cats()->GetPz(),prt_lab[0].Cats()->GetPt()/prt_lab[1].Cats()->GetPt());


FemtoPermutations++;
}//femto particles

//if(kstar>150 && kstar<250){
  Ghetto_mT_mTwrong->Fill(mT,mT_wrong);
//}


  Ghetto_rstar->Fill(rstar);
  Ghetto_kstar->Fill(kstar);
  Ghetto_kstar_rstar->Fill(kstar,rstar);
  //printf("===============\n");
  //printf("IsPrim: %i %i\n",prt_cm[0].IsUsefulPrimordial(),prt_cm[1].IsUsefulPrimordial());
  //printf("Particles (list): %s %s\n",ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());
  //printf("Particles (slct): %s %s\n",prt_cm[0].Trepni()->GetName().c_str(),prt_cm[1].Trepni()->GetName().c_str());
  //printf("Categorized as ");
  if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
    Ghetto_kstar_rstar_PP->Fill(kstar,rstar);
    //printf("PP\n");
  }
  else if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulPrimordial()){
    if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)){
      Ghetto_kstar_rstar_RP->Fill(kstar,rstar);//Ghetto_RP_AngleRcP1->Fill(AngleRcP1);
      //printf("RP\n");
    }
    else if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1)){
      Ghetto_kstar_rstar_PR->Fill(kstar,rstar);
      //printf("PR\n");
    }
  }
  else if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulProduct()){
    if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1)){
      Ghetto_kstar_rstar_PR->Fill(kstar,rstar);//Ghetto_RP_AngleRcP1->Fill(AngleRcP1);
      //printf("PR\n");
    }
    else if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0)){
      Ghetto_kstar_rstar_RP->Fill(kstar,rstar);
      //printf("RP\n");
    }
  }
  else{
    Ghetto_kstar_rstar_RR->Fill(kstar,rstar);
    //printf("RR\n");
  }
//usleep(200e3);


  Ghetto_kstar_rstar_mT->Fill(kstar,rstar,mT);
  if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
    Ghetto_kstar_rcore_mT->Fill(kstar,rcore,mT);
    //printf("%.0f %.2f %.0f\n",kstar,rcore,mT);
  }

  Ghetto_mT_rstar->Fill(mT,rstar);
  for(float ct : cos_th){
    Ghetto_mT_costh->Fill(mT,ct);
  }
  GhettoFemto_mT_kstar->Fill(mT,kstar);
  if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
    GhettoFemtoPrimordial_mT_kstar->Fill(mT,kstar);
  }
  
  //GhettoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
  //GhettoPrimReso[1] += (prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
  //GhettoPrimReso[2] += (!prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
  //GhettoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
  GhettoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
  GhettoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
  if(prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial()){
         if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0))GhettoPrimReso[1]++;
    else if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1))GhettoPrimReso[2]++;
    else printf("BAD BAD BAD 1\n");
  }
  if(prt_cm[1].IsUsefulPrimordial() && !prt_cm[0].IsUsefulPrimordial()){
         if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1))GhettoPrimReso[2]++;
    else if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0))GhettoPrimReso[1]++;
    else printf("BAD BAD BAD 2\n");
  }
  //if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)&&prt_cm[0].Trepni()->GetName()!=ListOfParticles.at(1)){
  //  GhettoPrimReso[1]
  //}

  if(kstar<FemtoLimit){
    if(mT>=Ghetto_MinMt && mT<=Ghetto_MaxMt){
      GhettoFemto_rstar->Fill(rstar);
    }
    GhettoFemto_mT_rstar->Fill(mT,rstar);
    //N.B. without this if statement, we get a different rcore, as it depends on the mass of the particles
    //i.e. the core is not the same for all species, even if the single particle emission IS the same
    if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){
      if(mT>=Ghetto_MinMt && mT<=Ghetto_MaxMt){
        GhettoFemto_rcore->Fill(rcore);
      }
      GhettoFemto_mT_rcore->Fill(mT,rcore);
      Ghetto_kstar_reso_mT->Fill(kstar,0,mT);
    }
    //GhettoFemtoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
    //GhettoFemtoPrimReso[1] += (prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
    //GhettoFemtoPrimReso[2] += (!prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
    //GhettoFemtoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
    GhettoFemtoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
    GhettoFemtoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
    if(prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulProduct()){
       if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(0)){
         GhettoFemtoPrimReso[1]++;
         Ghetto_kstar_reso_mT->Fill(kstar,1,mT);
       }
      else if(prt_cm[0].Trepni()->GetName()==ListOfParticles.at(1)){
        GhettoFemtoPrimReso[2]++;
        Ghetto_kstar_reso_mT->Fill(kstar,2,mT);
      }
      else printf("BAD BAD BAD 1\n");

    }
    if(prt_cm[1].IsUsefulPrimordial() && prt_cm[0].IsUsefulProduct()){
      if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(1)){
        GhettoFemtoPrimReso[2]++;
        Ghetto_kstar_reso_mT->Fill(kstar,2,mT);
      }
      else if(prt_cm[1].Trepni()->GetName()==ListOfParticles.at(0)){
        GhettoFemtoPrimReso[1]++;
        Ghetto_kstar_reso_mT->Fill(kstar,1,mT);
      }
      else printf("BAD BAD BAD 2\n");
    }
    if(prt_cm[1].IsUsefulProduct() && prt_cm[0].IsUsefulProduct()){
      Ghetto_kstar_reso_mT->Fill(kstar,3,mT);
    }
  }
//}
}

/////////////////////////////// RAFA
  else if(SDIM==3){

  }
///////////////////////////////

//printf("OUT OF THE GHETTO\n");
///////////////////////////
}
      delete [] prt_cm;
      delete [] prt_lab;
    }//permutations over possible multiplets



  for(CecaParticle* particle : Primordial){
    delete particle;
  }
  for(CecaParticle* particle : Primary){
    delete particle;
  }
  //return Permutations.size();
  return FemtoPermutations;
}



//random pick two particles, in terms of species, either proton or Reso
//random pick their momentum (fixed and the same)
//1)  propagate and decay them in lab
//    evaluate their rstar value in the end, save it in a histo
//    evaluate their r_core (distance between mothers) and save it in a histo
//forget about 2) and 3) for now, this can be made after seeing the results from 1)
//2)  with the same particles, do the classical approach:
//    decay them asap and go to the CM frame. Inside, fix one particle and its daughter,
//    but shift the second one (and daughter) such that the relative distance between the
//    two mothers is r_core (random sampled).
//3)  By trial, find a correspondence between the r_core and r_SP to get the best
//    coincidence between the two methods.
void CECA::GhettoTest1(const unsigned NumPairs, const float r_SP, const float p_SP){
  //const float M_proton = (938+1116)*0.5;
  //const float M_reso = (1362+1462)*0.5;
  //const float tau_reso = (1.65+4.69)*0.5;

  //const float M_proton = 938;
  //const float M_reso = 1362;
  ////const float tau_reso = 1.65;
  //const float tau_reso = 5;
  //const float M_pi = 140;

  const float M_proton = 140;
  const float M_reso = 1124;
  const float tau_reso = 1.5;
  const float M_pi = 938;

  const float F_prim = 0.36;
  //const float F_prim = 1.0;
  //const float F_prim = 0.0;

//some issue with the random numbers, EVEN if they have their own class....
//bummer
  const unsigned NumThr = 7;
  const double SP_rew = 1.0;
  const double PK_rew = 0.0;

  const unsigned NumRadBins = 1024;
  const float RadMin = 0;
  const float RadMax = 16;

  const unsigned NumMomBins = 64;
  const float MomMin = 0;
  const float MomMax = 4096;
//printf("h1 1\n");
  if(Old_source) delete Old_source;
  Old_source = new DLM_CleverMcLevyResoTM();
  Old_source->InitStability(1,2-1e-6,2+1e-6);
  Old_source->InitScale(38,0.15,2.0);
  Old_source->InitRad(257*2,0,64);
  Old_source->InitType(2);
  Old_source->SetUpReso(0,1.-F_prim);
  Old_source->SetUpReso(1,1.-F_prim);
//printf("h1 2\n");
  omp_set_dynamic(0);
  omp_set_num_threads(NumThr);
  DLM_Random** RanGen = new DLM_Random* [NumThr];
  for(unsigned uThr=0; uThr<NumThr; uThr++){
    RanGen[uThr] = new DLM_Random(uThr+1);
  }

  //each histogram here is 2dim, first dim is kstar (coarse binning)
  //DLM_Histo<float>* dlm_r = new DLM_Histo<float> [NumThr+1];
  DLM_Histo<float>* dlm_rstar = new DLM_Histo<float> [NumThr];
  DLM_Histo<float>* dlm_rcore = new DLM_Histo<float> [NumThr];

  //DLM_Histo<float>* old_rstar = new DLM_Histo<float> [NumThr];
  //DLM_Histo<float>* old_rcore = new DLM_Histo<float> [NumThr];

  for(unsigned uThr=0; uThr<NumThr; uThr++){
    //dlm_r[uThr].SetUp(2);
    //dlm_r[uThr].SetUp(0,NumMomBins,MomMin,MomMax);
    //dlm_r[uThr].SetUp(1,NumRadBins,RadMin,RadMax);
    //dlm_r[uThr].Initialize();

    dlm_rstar[uThr].SetUp(2);
    dlm_rstar[uThr].SetUp(0,NumMomBins,MomMin,MomMax);
    dlm_rstar[uThr].SetUp(1,NumRadBins,RadMin,RadMax);
    dlm_rstar[uThr].Initialize();

    dlm_rcore[uThr].SetUp(2);
    dlm_rcore[uThr].SetUp(0,NumMomBins,MomMin,MomMax);
    dlm_rcore[uThr].SetUp(1,NumRadBins,RadMin,RadMax);
    dlm_rcore[uThr].Initialize();

    //old_rstar[uThr].SetUp(2);
    //old_rstar[uThr].SetUp(0,NumMomBins,MomMin,MomMax);
    //old_rstar[uThr].SetUp(1,NumRadBins,RadMin,RadMax);
    //old_rstar[uThr].Initialize();

    //old_rcore[uThr].SetUp(2);
    //old_rcore[uThr].SetUp(0,NumMomBins,MomMin,MomMax);
    //old_rcore[uThr].SetUp(1,NumRadBins,RadMin,RadMax);
    //old_rcore[uThr].Initialize();
  }
  if(Ghetto_rstar) delete Ghetto_rstar;
  Ghetto_rstar = new DLM_Histo<float>(dlm_rstar[0]);
  if(Ghetto_rcore) delete Ghetto_rcore;
  Ghetto_rcore = new DLM_Histo<float>(dlm_rcore[0]);
  if(GhettOld_rstar) delete GhettOld_rstar;

  if(GhettOld_rstar) delete GhettOld_rstar;
  GhettOld_rstar = new DLM_Histo<float>();
  GhettOld_rstar->SetUp(1);
  GhettOld_rstar->SetUp(0,NumRadBins,RadMin,RadMax);
  GhettOld_rstar->Initialize();
//printf("h1 3\n");
  if(Old_rstar) delete Old_rstar;
  //Old_rstar = new DLM_Histo<float>(dlm_rstar[0]);
  Old_rstar = new DLM_Histo<float>();
  Old_rstar->SetUp(1);
  Old_rstar->SetUp(0,NumRadBins,RadMin,RadMax);
  Old_rstar->Initialize();
  if(Old_rcore) delete Old_rcore;
  //Old_rcore = new DLM_Histo<float>(dlm_rcore[0]);
  Old_rcore = new DLM_Histo<float>();
  Old_rcore->SetUp(1);
  Old_rcore->SetUp(0,NumRadBins,RadMin,RadMax);
  Old_rcore->Initialize();

  if(Old_CosRcP1) delete Old_CosRcP1;
  //Old_rcore = new DLM_Histo<float>(dlm_rcore[0]);
  Old_CosRcP1 = new DLM_Histo<float>();
  Old_CosRcP1->SetUp(1);
  Old_CosRcP1->SetUp(0,NumRadBins,-1,1);
  Old_CosRcP1->Initialize();

  if(Old_CosRcP2) delete Old_CosRcP2;
  //Old_CosRcP2 = new DLM_Histo<float>(dlm_rcore[0]);
  Old_CosRcP2 = new DLM_Histo<float>();
  Old_CosRcP2->SetUp(1);
  Old_CosRcP2->SetUp(0,NumRadBins,-1,1);
  Old_CosRcP2->Initialize();

  if(Old_CosP1P2) delete Old_CosP1P2;
  //Old_CosP1P2 = new DLM_Histo<float>(dlm_rcore[0]);
  Old_CosP1P2 = new DLM_Histo<float>();
  Old_CosP1P2->SetUp(1);
  Old_CosP1P2->SetUp(0,NumRadBins,-1,1);
  Old_CosP1P2->Initialize();

  if(Old_RcP1) delete Old_RcP1;
  //Old_rcore = new DLM_Histo<float>(dlm_rcore[0]);
  Old_RcP1 = new DLM_Histo<float>();
  Old_RcP1->SetUp(1);
  Old_RcP1->SetUp(0,NumRadBins,0,3.1416);
  Old_RcP1->Initialize();

  if(Old_RcP2) delete Old_RcP2;
  //Old_RcP2 = new DLM_Histo<float>(dlm_rcore[0]);
  Old_RcP2 = new DLM_Histo<float>();
  Old_RcP2->SetUp(1);
  Old_RcP2->SetUp(0,NumRadBins,0,3.1416);
  Old_RcP2->Initialize();

  if(Old_P1P2) delete Old_P1P2;
  //Old_P1P2 = new DLM_Histo<float>(dlm_rcore[0]);
  Old_P1P2 = new DLM_Histo<float>();
  Old_P1P2->SetUp(1);
  Old_P1P2->SetUp(0,NumRadBins,0,3.1416);
  Old_P1P2->Initialize();

//printf("h1 4\n");
  #pragma omp parallel for
  for(unsigned uPair=0; uPair<NumPairs; uPair++){
//printf("uPair = %u\n",uPair);
    unsigned ThId = omp_get_thread_num();
//if(omp_get_num_threads()==1) ThId = uPair%7;
//ThId = uPair%7;
    float r_x,r_y,r_z,p_x,p_y,p_z,rad,mom,mass,kstar,rstar,rcore;
    bool prim1,prim2;
    CatsParticle* Daughter = NULL;
    //for now, we will be setting the initial time (in LAB) to zero

    //#pragma omp critical
    {
    r_x = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    r_y = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    r_z = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    p_x = RanGen[ThId]->Gauss(0,p_SP);
    p_y = RanGen[ThId]->Gauss(0,p_SP);
    p_z = RanGen[ThId]->Gauss(0,p_SP);
    r_x += p_x * r_SP / p_SP * PK_rew ;
    r_y += p_y * r_SP / p_SP * PK_rew;
    r_z += p_z * r_SP / p_SP * PK_rew;
    prim1 = (RanGen[ThId]->Uniform(0,1)<F_prim);
    }
/*
    //#pragma omp critical
    {
    //if(uPair<64){
    if(ThId==3){
      static int cntr=0;
      if(cntr<16){
        printf("%u(%u): %.3f %.3f %.3f %.0f %.0f %.0f %i\n",
        uPair,ThId,r_x,r_y,r_z,p_x,p_y,p_z,prim1);
      }
      cntr++;
    }
    }
*/
//prim = false;
    //rad = sqrt(r_x*r_x+r_y*r_y+r_z*r_z);
    rad = 0;
    mom = sqrt(p_x*p_x+p_y*p_y+p_z*p_z);
    mass = prim1?M_proton:M_reso;
    CatsParticle Particle1;
    Particle1.SetMass(mass);
    Particle1.SetWidth( 1./(tau_reso*FmToNu) );
    Particle1.Set(rad,r_x,r_y,r_z,sqrt(mass*mass+mom*mom),p_x,p_y,p_z);
    Particle1.SetDecayRanGen(RanGen[ThId]);
//Particle1.Print();
    if(!prim1) Daughter = Particle1.Decay(M_proton,M_pi,false);
    CatsParticle FsiParticle1 = Daughter?Daughter[0]:Particle1;
    if(Daughter) {delete[]Daughter; Daughter=NULL;}
    //#pragma omp critical
    {
    r_x = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    r_y = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    r_z = RanGen[ThId]->Gauss(0,r_SP*SP_rew);
    p_x = RanGen[ThId]->Gauss(0,p_SP);
    p_y = RanGen[ThId]->Gauss(0,p_SP);
    p_z = RanGen[ThId]->Gauss(0,p_SP);
    r_x += p_x * r_SP / p_SP * PK_rew ;
    r_y += p_y * r_SP / p_SP * PK_rew;
    r_z += p_z * r_SP / p_SP * PK_rew;
    prim2 = (RanGen[ThId]->Uniform(0,1)<F_prim);
    }
//prim = false;
    //rad = sqrt(r_x*r_x+r_y*r_y+r_z*r_z);
    rad = 0;
    mom = sqrt(p_x*p_x+p_y*p_y+p_z*p_z);
    mass = prim2?M_proton:M_reso;
    CatsParticle Particle2;
    Particle2.SetMass(mass);
    Particle2.SetWidth( 1./(tau_reso*FmToNu) );
    Particle2.Set(rad,r_x,r_y,r_z,sqrt(mass*mass+mom*mom),p_x,p_y,p_z);
    Particle2.SetDecayRanGen(RanGen[ThId]);
//Particle2.Print();
    if(!prim2) Daughter = Particle2.Decay(M_proton,M_pi,false);
    CatsParticle FsiParticle2 = Daughter?Daughter[0]:Particle2;
    if(Daughter) {delete[]Daughter; Daughter=NULL;}
/*
    if(!prim1&&false){
      CatsLorentzVector DiffBefore = Particle2-Particle1;
      printf("-------------\n");
      printf("%i%i\n",!prim1,!prim2);
      //printf("P1  = %f\n",Particle1.GetP());
      //printf("D   = %f\n",DiffBefore.GetP());
      printf("P1BB: "); Particle1.Print();
      printf("P2BB: "); Particle2.Print();
      printf(" B1BB %f\n",Particle1.Beta());
      printf(" B2BB %f\n",Particle2.Beta());
    }
*/
    //printf("FSI_1: "); FsiParticle1.Print();
    //printf("FSI_2: "); FsiParticle2.Print();
    //printf("MTR_1: "); Particle1.Print();
    //printf(" %.2f %.2f %.2f\n",FsiParticle1.GetX()/Particle1.GetX(),FsiParticle1.GetY()/Particle1.GetY(),FsiParticle1.GetZ()/Particle1.GetZ());
    //printf("MTR_2: "); Particle2.Print();
/*
    if(!prim1&&prim2&&kstar<FemtoLimit){
      printf("-----------------------\n");
      printf("%i%i\n",!prim1,!prim2);
      printf(" rx ry rz: %5.2f %5.2f %5.2f\n",Particle1.GetX(),Particle1.GetY(),Particle1.GetZ());
      printf(" px py pz: %5.0f %5.0f %5.0f\n",Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
      printf(" dx dy dz: %5.2f %5.2f %5.2f\n",
                Particle1.GetX()+Particle1.GetPx()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm,
                Particle1.GetY()+Particle1.GetPy()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm,
                Particle1.GetZ()+Particle1.GetPz()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm);
      printf(" fx fy fz: %5.2f %5.2f %5.2f\n",FsiParticle1.GetX(),FsiParticle1.GetY(),FsiParticle1.GetZ());
      printf("BOOOOOOOOOOOST\n");
    }
*/


    CatsLorentzVector BoostVector = FsiParticle1+FsiParticle2;
    FsiParticle1.Boost(BoostVector);
    FsiParticle2.Boost(BoostVector);
    Particle1.Boost(BoostVector);
    Particle2.Boost(BoostVector);

    kstar = sqrt( pow(FsiParticle1.GetPx()-FsiParticle2.GetPx(),2.)+
                  pow(FsiParticle1.GetPy()-FsiParticle2.GetPy(),2.)+
                  pow(FsiParticle1.GetPz()-FsiParticle2.GetPz(),2.));
    rstar = sqrt( pow(FsiParticle1.GetX()-FsiParticle2.GetX(),2.)+
                  pow(FsiParticle1.GetY()-FsiParticle2.GetY(),2.)+
                  pow(FsiParticle1.GetZ()-FsiParticle2.GetZ(),2.));
    rcore = sqrt( pow(Particle1.GetX()-Particle2.GetX(),2.)+
                  pow(Particle1.GetY()-Particle2.GetY(),2.)+
                  pow(Particle1.GetZ()-Particle2.GetZ(),2.));

    double core_x = Particle2.GetX()-Particle1.GetX();
    double core_y = Particle2.GetY()-Particle1.GetY();
    double core_z = Particle2.GetZ()-Particle1.GetZ();
    double core_r = sqrt(core_x*core_x+core_y*core_y+core_z*core_z);
    //funny thing, but for both old and new method to agree I need the same Tau as in CATStools Decay...
    //so it will only really work if we fix the width
    double BGT_1 = Particle1.Beta()*Particle1.Gamma()/Particle1.GetWidth()*NuToFm;
    double BGT_2 = Particle2.Beta()*Particle2.Gamma()/Particle2.GetWidth()*NuToFm;
    //double BGT_1 = Particle1.Beta()*Particle1.Gamma()*RanGen[ThId]->Exponential(Particle1.GetWidth())*NuToFm;
    //double BGT_2 = Particle2.Beta()*Particle2.Gamma()*RanGen[ThId]->Exponential(Particle2.GetWidth())*NuToFm;
    double CosRcP1 = (core_x*Particle1.GetPx()+core_y*Particle1.GetPy()+core_z*Particle1.GetPz())
                    /(core_r*Particle1.GetP());
    double CosRcP2 = (core_x*Particle2.GetPx()+core_y*Particle2.GetPy()+core_z*Particle2.GetPz())
                    /(core_r*Particle2.GetP());
    double CosP1P2 = (Particle1.GetPx()*Particle2.GetPx()+Particle1.GetPy()*Particle2.GetPy()+Particle1.GetPz()*Particle2.GetPz())
                    /(Particle1.GetP()*Particle2.GetP());
    if(CosRcP1<-1) CosRcP1=-1;
    if(CosRcP1>1) CosRcP1=1;
    if(CosRcP2<-1) CosRcP2=-1;
    if(CosRcP2>1) CosRcP2=1;
    if(CosP1P2<-1) CosP1P2=-1;
    if(CosP1P2>1) CosP1P2=1;
/*
    if(!prim1&&prim2&&kstar<FemtoLimit){
      printf(" rx ry rz: %5.2f %5.2f %5.2f\n",Particle1.GetX(),Particle1.GetY(),Particle1.GetZ());
      printf(" px py pz: %5.0f %5.0f %5.0f\n",Particle1.GetPx(),Particle1.GetPy(),Particle1.GetPz());
      printf(" dx dy dz: %5.2f %5.2f %5.2f\n",
                Particle1.GetX()+Particle1.GetPx()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm,
                Particle1.GetY()+Particle1.GetPy()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm,
                Particle1.GetZ()+Particle1.GetPz()/Particle1.GetWidth()/Particle1.GetMass()*NuToFm);
      printf(" fx fy fz: %5.2f %5.2f %5.2f\n",FsiParticle1.GetX(),FsiParticle1.GetY(),FsiParticle1.GetZ());
      printf("SECOND:\n");
      printf(" rx ry rz: %5.2f %5.2f %5.2f\n",Particle2.GetX(),Particle2.GetY(),Particle2.GetZ());
      printf(" px py pz: %5.0f %5.0f %5.0f\n",Particle2.GetPx(),Particle2.GetPy(),Particle2.GetPz());
      printf(" fx fy fz: %5.2f %5.2f %5.2f\n",FsiParticle2.GetX(),FsiParticle2.GetY(),FsiParticle2.GetZ());
      printf(" dr = %.3f (%.3f)\n",
        sqrt( pow(FsiParticle1.GetX()-FsiParticle2.GetX(),2.)+
              pow(FsiParticle1.GetY()-FsiParticle2.GetY(),2.)+
              pow(FsiParticle1.GetZ()-FsiParticle2.GetZ(),2.)),
              rstar
            );
      printf(" CosRcP1 = %.3f\n",CosRcP1);
      printf(" RSTAR = %.3f\n",sqrt(rcore*rcore+BGT_1*BGT_1-2.*rcore*BGT_1*CosRcP1));
      usleep(100e3);
    }
*/

    //printf("kstar = %0f; rstar = %.3f; rcore = %.3f\n",kstar,rstar,rcore);

/*
if(!prim1&&false){
  CatsLorentzVector DiffAfter = Particle2-Particle1;
  //printf(" BGT_1 = %f vs %f\n",BGT_1,Particle1.GetP()/Particle1.Mag()/Particle1.GetWidth()*NuToFm);
  //printf(" P1 = %f\n",Particle1.GetP());
  //printf(" D  = %f\n",DiffAfter.GetP());
  //printf("BGT_2 = %f\n",BGT_2);
  //printf("\n");
  printf("P1AB: "); Particle1.Print();
  printf("P2AB: "); Particle2.Print();
  printf(" B1AB %f\n",Particle1.Beta());
  printf(" B2AB %f\n",Particle2.Beta());
  usleep(200e3);
}
*/



    //CosP1P2 = RanGen[ThId]->Uniform(-1,1);

    if(kstar<FemtoLimit){
      #pragma omp critical
      {
        double OldSoureValue;
      if(prim1&&prim2){
        OldSoureValue = core_r;
      }
      else if(prim1&&!prim2){
        Old_source->AddBGT_PR(BGT_2,CosRcP2);
        //Old_source->AddBGT_PR(BGT_1,-1);
        //Old_source->AddBGT_PR(0,0.5);
        OldSoureValue = sqrt(core_r*core_r+BGT_2*BGT_2+2.*core_r*BGT_2*CosRcP2);
      }
      else if(!prim1&&prim2){
        Old_source->AddBGT_RP(BGT_1,CosRcP1);
        //Old_source->AddBGT_RP(BGT_1,1);
        //Old_source->AddBGT_RP(0,0.5);
        OldSoureValue = sqrt(core_r*core_r+BGT_1*BGT_1-2.*core_r*BGT_1*CosRcP1);
      }
      else{
        Old_source->AddBGT_RR(BGT_1,CosRcP1,BGT_2,CosRcP2,CosP1P2);
        //Old_source->AddBGT_RR(BGT_1,-1,BGT_2,1,-1);
//printf("%f %f %f %f %f\n", BGT_1,CosRcP1,BGT_2,CosRcP2,CosP1P2);
//printf("%f %f %f %f %f\n", BGT_1,acos(CosRcP1)*180./3.14159,BGT_2,acos(CosRcP2)*180./3.14159,acos(CosP1P2)*180./3.14159);
//printf("%.0f: %.2f (%.2f) -> %.2f (%.2f)\n",kstar,core_r,rcore,
//sqrt(core_r*core_r+BGT_1*BGT_1+BGT_2*BGT_2-2.*core_r*BGT_1*CosRcP1+
//2.*core_r*BGT_2*CosRcP2-2.*core_r*BGT_1*BGT_2*CosP1P2),rstar);
        OldSoureValue = sqrt(core_r*core_r+BGT_1*BGT_1+BGT_2*BGT_2-2.*core_r*BGT_1*CosRcP1+
        2.*core_r*BGT_2*CosRcP2-2.*BGT_1*BGT_2*CosP1P2);
//usleep(200e3);
        //Old_source->AddBGT_RR(BGT_1,0.5,BGT_1,0.5,1);
      }
      Old_CosRcP1->AddAt(&CosRcP1);
      Old_CosRcP2->AddAt(&CosRcP2);
      if(!prim1&&!prim2) {Old_CosP1P2->AddAt(&CosP1P2);}

      GhettOld_rstar->AddAt(&OldSoureValue);

      //if(!prim1&&prim2){
      //if(rstar>6){
      //  printf("%i%i ",!prim1,!prim2);
      //  printf("%f vs %f\n",rstar,OldSoureValue);
      //  usleep(200e3);
      //}



      double Angle;
      Angle = acos(CosRcP1); Old_RcP1->AddAt(&Angle);
      Angle = acos(CosRcP2); Old_RcP2->AddAt(&Angle);
      if(!prim1&&!prim2) {Angle = acos(CosP1P2); Old_P1P2->AddAt(&Angle);}
      }
    }

    double axis_val[2];

    axis_val[0] = kstar;
    axis_val[1] = rstar;
    dlm_rstar[ThId].AddAt(axis_val);

    axis_val[0] = kstar;
    axis_val[1] = rcore;
    dlm_rcore[ThId].AddAt(axis_val);

  }

  Old_source->InitNumMcIter(1000000);
  for(unsigned uBin=0; uBin<NumRadBins; uBin++){
    double rad = Old_rcore->GetBinCenter(0,uBin);
    double pars[2];
    pars[0] = r_SP;
    pars[1] = 2.;
    Old_rstar->SetBinContent(uBin,Old_source->RootEval(&rad,pars));
    Old_rcore->SetBinContent(uBin,GaussSourceTF1(&rad,pars));
  }

  for(unsigned uThr=0; uThr<NumThr; uThr++){
    *Ghetto_rstar += dlm_rstar[uThr];
    *Ghetto_rcore += dlm_rcore[uThr];
  }

  Ghetto_rstar->ComputeError();
  Ghetto_rcore->ComputeError();

  for(unsigned uThr=0; uThr<NumThr; uThr++){
    delete RanGen[uThr];
  }
  delete [] RanGen;

  delete [] dlm_rstar;
  delete [] dlm_rcore;
}

void CECA::GhettoInit(){

  if(Ghetto_RP_AngleRcP1) delete Ghetto_RP_AngleRcP1;
  Ghetto_RP_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RP_AngleRcP1->SetUp(1);
  Ghetto_RP_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RP_AngleRcP1->Initialize();

  if(Ghetto_RP_AngleRcP1) delete Ghetto_RP_AngleRcP1;
  Ghetto_RP_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RP_AngleRcP1->SetUp(1);
  Ghetto_RP_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RP_AngleRcP1->Initialize();

  if(Ghetto_PR_AngleRcP2) delete Ghetto_PR_AngleRcP2;
  Ghetto_PR_AngleRcP2 = new DLM_Histo<float>();
  Ghetto_PR_AngleRcP2->SetUp(1);
  Ghetto_PR_AngleRcP2->SetUp(0,128,0,Pi);
  Ghetto_PR_AngleRcP2->Initialize();

  if(Ghetto_RR_AngleRcP1) delete Ghetto_RR_AngleRcP1;
  Ghetto_RR_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RR_AngleRcP1->SetUp(1);
  Ghetto_RR_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleRcP1->Initialize();

  if(Ghetto_RR_AngleRcP2) delete Ghetto_RR_AngleRcP2;
  Ghetto_RR_AngleRcP2 = new DLM_Histo<float>();
  Ghetto_RR_AngleRcP2->SetUp(1);
  Ghetto_RR_AngleRcP2->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleRcP2->Initialize();

  if(Ghetto_RR_AngleP1P2) delete Ghetto_RR_AngleP1P2;
  Ghetto_RR_AngleP1P2 = new DLM_Histo<float>();
  Ghetto_RR_AngleP1P2->SetUp(1);
  Ghetto_RR_AngleP1P2->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleP1P2->Initialize();

  if(Ghetto_ScatteringAngle) delete Ghetto_ScatteringAngle;
  Ghetto_ScatteringAngle = new DLM_Histo<float>();
  Ghetto_ScatteringAngle->SetUp(1);
  Ghetto_ScatteringAngle->SetUp(0,128,0,Pi);
  Ghetto_ScatteringAngle->Initialize();

  if(Ghetto_PP_AngleRcP1) delete Ghetto_PP_AngleRcP1;
  Ghetto_PP_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_PP_AngleRcP1->SetUp(1);
  Ghetto_PP_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_PP_AngleRcP1->Initialize();

  if(Ghetto_PP_AngleRcP2) delete Ghetto_PP_AngleRcP2;
  Ghetto_PP_AngleRcP2 = new DLM_Histo<float>();
  Ghetto_PP_AngleRcP2->SetUp(1);
  Ghetto_PP_AngleRcP2->SetUp(0,128,0,Pi);
  Ghetto_PP_AngleRcP2->Initialize();

  if(Ghetto_PP_AngleP1P2) delete Ghetto_PP_AngleP1P2;
  Ghetto_PP_AngleP1P2 = new DLM_Histo<float>();
  Ghetto_PP_AngleP1P2->SetUp(1);
  Ghetto_PP_AngleP1P2->SetUp(0,128,0,Pi);
  Ghetto_PP_AngleP1P2->Initialize();


  if(GhettoSP_pT_th) delete GhettoSP_pT_th;
  GhettoSP_pT_th = new DLM_Histo<float>();
  GhettoSP_pT_th->SetUp(2);
  GhettoSP_pT_th->SetUp(0,64,0,4096);
  GhettoSP_pT_th->SetUp(1,64,-0.1,3.2);
  GhettoSP_pT_th->Initialize();

  if(GhettoSP_pT_1) delete GhettoSP_pT_1;
  GhettoSP_pT_1 = new DLM_Histo<float>();
  GhettoSP_pT_1->SetUp(1);
  GhettoSP_pT_1->SetUp(0,512,0,16384);
  GhettoSP_pT_1->Initialize();

  if(GhettoSP_pT_2) delete GhettoSP_pT_2;
  GhettoSP_pT_2 = new DLM_Histo<float>();
  GhettoSP_pT_2->SetUp(1);
  GhettoSP_pT_2->SetUp(0,512,0,16384);
  GhettoSP_pT_2->Initialize();

  if(GhettoSPr_X) delete GhettoSPr_X;
  GhettoSPr_X = new DLM_Histo<float>();
  GhettoSPr_X->SetUp(1);
  GhettoSPr_X->SetUp(0,256,-24,24);
  GhettoSPr_X->Initialize();

  if(GhettoSPr_Y) delete GhettoSPr_Y;
  GhettoSPr_Y = new DLM_Histo<float>();
  GhettoSPr_Y->SetUp(1);
  GhettoSPr_Y->SetUp(0,256,-24,24);
  GhettoSPr_Y->Initialize();

  if(GhettoSPr_Z) delete GhettoSPr_Z;
  GhettoSPr_Z = new DLM_Histo<float>();
  GhettoSPr_Z->SetUp(1);
  GhettoSPr_Z->SetUp(0,256,-24,24);
  GhettoSPr_Z->Initialize();

  if(GhettoSPr_Rho) delete GhettoSPr_Rho;
  GhettoSPr_Rho = new DLM_Histo<float>();
  GhettoSPr_Rho->SetUp(1);
  GhettoSPr_Rho->SetUp(0,256,0,32);
  GhettoSPr_Rho->Initialize();

  if(GhettoSPr_R) delete GhettoSPr_R;
  GhettoSPr_R = new DLM_Histo<float>();
  GhettoSPr_R->SetUp(1);
  GhettoSPr_R->SetUp(0,256,0,32);
  GhettoSPr_R->Initialize();

  if(GhettoSP_X) delete GhettoSP_X;
  GhettoSP_X = new DLM_Histo<float>();
  GhettoSP_X->SetUp(1);
  GhettoSP_X->SetUp(0,256,-24,24);
  GhettoSP_X->Initialize();

  if(GhettoSP_Y) delete GhettoSP_Y;
  GhettoSP_Y = new DLM_Histo<float>();
  GhettoSP_Y->SetUp(1);
  GhettoSP_Y->SetUp(0,256,-24,24);
  GhettoSP_Y->Initialize();

  if(GhettoSP_Z) delete GhettoSP_Z;
  GhettoSP_Z = new DLM_Histo<float>();
  GhettoSP_Z->SetUp(1);
  GhettoSP_Z->SetUp(0,256,-24,24);
  GhettoSP_Z->Initialize();

  if(GhettoSP_Rho) delete GhettoSP_Rho;
  GhettoSP_Rho = new DLM_Histo<float>();
  GhettoSP_Rho->SetUp(1);
  GhettoSP_Rho->SetUp(0,256,0,32);
  GhettoSP_Rho->Initialize();

  if(GhettoSP_R) delete GhettoSP_R;
  GhettoSP_R = new DLM_Histo<float>();
  GhettoSP_R->SetUp(1);
  GhettoSP_R->SetUp(0,256,0,32);
  GhettoSP_R->Initialize();

  if(Ghetto_rstar) delete Ghetto_rstar;
  Ghetto_rstar = new DLM_Histo<float>();
  Ghetto_rstar->SetUp(1);
  Ghetto_rstar->SetUp(0,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_rstar->Initialize();

  if(Ghetto_rcore) delete Ghetto_rcore;
  Ghetto_rcore = new DLM_Histo<float>();
  Ghetto_rcore->SetUp(1);
  Ghetto_rcore->SetUp(0,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_rcore->Initialize();

  if(Ghetto_kstar) delete Ghetto_kstar;
  Ghetto_kstar = new DLM_Histo<float>();
  Ghetto_kstar->SetUp(1);
  Ghetto_kstar->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar->Initialize();

  if(Ghetto_kstar_rstar) delete Ghetto_kstar_rstar;
  Ghetto_kstar_rstar = new DLM_Histo<float>();
  Ghetto_kstar_rstar->SetUp(2);
  Ghetto_kstar_rstar->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_kstar_rstar->Initialize();

  if(Ghetto_kstar_rstar_PP) delete Ghetto_kstar_rstar_PP;
  Ghetto_kstar_rstar_PP = new DLM_Histo<float>();
  Ghetto_kstar_rstar_PP->SetUp(2);
  Ghetto_kstar_rstar_PP->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar_PP->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_kstar_rstar_PP->Initialize();

  if(Ghetto_kstar_rstar_PR) delete Ghetto_kstar_rstar_PR;
  Ghetto_kstar_rstar_PR = new DLM_Histo<float>();
  Ghetto_kstar_rstar_PR->SetUp(2);
  Ghetto_kstar_rstar_PR->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar_PR->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_kstar_rstar_PR->Initialize();

  if(Ghetto_kstar_rstar_RP) delete Ghetto_kstar_rstar_RP;
  Ghetto_kstar_rstar_RP = new DLM_Histo<float>();
  Ghetto_kstar_rstar_RP->SetUp(2);
  Ghetto_kstar_rstar_RP->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar_RP->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_kstar_rstar_RP->Initialize();

  if(Ghetto_kstar_rstar_RR) delete Ghetto_kstar_rstar_RR;
  Ghetto_kstar_rstar_RR = new DLM_Histo<float>();
  Ghetto_kstar_rstar_RR->SetUp(2);
  Ghetto_kstar_rstar_RR->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar_RR->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_kstar_rstar_RR->Initialize();

  if(Ghetto_kstar_rstar_mT) delete Ghetto_kstar_rstar_mT;
  Ghetto_kstar_rstar_mT = new DLM_Histo<float>();
  Ghetto_kstar_rstar_mT->SetUp(3);
  Ghetto_kstar_rstar_mT->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rstar_mT->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  if(Ghetto_MtBins){
    Ghetto_kstar_rstar_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_kstar_rstar_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_kstar_rstar_mT->Initialize();

  if(Ghetto_kstar_rcore_mT) delete Ghetto_kstar_rcore_mT;
  Ghetto_kstar_rcore_mT = new DLM_Histo<float>();
  Ghetto_kstar_rcore_mT->SetUp(3);
  Ghetto_kstar_rcore_mT->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_rcore_mT->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  if(Ghetto_MtBins){
    Ghetto_kstar_rcore_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_kstar_rcore_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_kstar_rcore_mT->Initialize();

  if(Ghetto_kstar_reso_mT) delete Ghetto_kstar_reso_mT;
  Ghetto_kstar_reso_mT = new DLM_Histo<float>();
  Ghetto_kstar_reso_mT->SetUp(3);
  Ghetto_kstar_reso_mT->SetUp(0,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  Ghetto_kstar_reso_mT->SetUp(1,4,-0.5,3.5);
  if(Ghetto_MtBins){
    Ghetto_kstar_reso_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_kstar_reso_mT->SetUp(2,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_kstar_reso_mT->Initialize();


  if(Ghetto_mT_rstar) delete Ghetto_mT_rstar;
  Ghetto_mT_rstar = new DLM_Histo<float>();
  Ghetto_mT_rstar->SetUp(2);
  if(Ghetto_MtBins){
    Ghetto_mT_rstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_mT_rstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_mT_rstar->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  Ghetto_mT_rstar->Initialize();


  if(GhettoFemto_rstar) delete GhettoFemto_rstar;
  GhettoFemto_rstar = new DLM_Histo<float>();
  GhettoFemto_rstar->SetUp(1);
  GhettoFemto_rstar->SetUp(0,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  GhettoFemto_rstar->Initialize();

  if(GhettoFemto_rcore) delete GhettoFemto_rcore;
  GhettoFemto_rcore = new DLM_Histo<float>();
  GhettoFemto_rcore->SetUp(1);
  GhettoFemto_rcore->SetUp(0,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  GhettoFemto_rcore->Initialize();

  if(Ghetto_mT_costh) delete Ghetto_mT_costh;
  Ghetto_mT_costh = new DLM_Histo<float>();
  Ghetto_mT_costh->SetUp(2);
  if(Ghetto_MtBins){
    Ghetto_mT_costh->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_mT_costh->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_mT_costh->SetUp(1,128,-1,1);
  Ghetto_mT_costh->Initialize();


  if(GhettoFemto_mT_rstar) delete GhettoFemto_mT_rstar;
  GhettoFemto_mT_rstar = new DLM_Histo<float>();
  GhettoFemto_mT_rstar->SetUp(2);
  if(Ghetto_MtBins){
    GhettoFemto_mT_rstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    GhettoFemto_mT_rstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  GhettoFemto_mT_rstar->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  GhettoFemto_mT_rstar->Initialize();

  if(GhettoFemto_mT_rcore) delete GhettoFemto_mT_rcore;
  GhettoFemto_mT_rcore = new DLM_Histo<float>();
  GhettoFemto_mT_rcore->SetUp(2);
  if(Ghetto_MtBins){
    GhettoFemto_mT_rcore->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    GhettoFemto_mT_rcore->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  GhettoFemto_mT_rcore->SetUp(1,Ghetto_NumRadBins,Ghetto_RadMin,Ghetto_RadMax);
  GhettoFemto_mT_rcore->Initialize();

  if(GhettoFemto_mT_kstar) delete GhettoFemto_mT_kstar;
  GhettoFemto_mT_kstar = new DLM_Histo<float>();
  GhettoFemto_mT_kstar->SetUp(2);
  if(Ghetto_MtBins){
    GhettoFemto_mT_kstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    GhettoFemto_mT_kstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }

  GhettoFemto_mT_kstar->SetUp(1,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  GhettoFemto_mT_kstar->Initialize();


  if(GhettoFemtoPrimordial_mT_kstar) delete GhettoFemtoPrimordial_mT_kstar;
  GhettoFemtoPrimordial_mT_kstar = new DLM_Histo<float>();
  GhettoFemtoPrimordial_mT_kstar->SetUp(2);
  if(Ghetto_MtBins){
    GhettoFemtoPrimordial_mT_kstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    GhettoFemtoPrimordial_mT_kstar->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }

  GhettoFemtoPrimordial_mT_kstar->SetUp(1,Ghetto_NumMomBins,Ghetto_MomMin,Ghetto_MomMax);
  GhettoFemtoPrimordial_mT_kstar->Initialize();





  if(GhettoFemto_pT1_pT2) delete GhettoFemto_pT1_pT2;
  GhettoFemto_pT1_pT2 = new DLM_Histo<float>();
  GhettoFemto_pT1_pT2->SetUp(2);
  GhettoFemto_pT1_pT2->SetUp(0,1024,0,4096);
  GhettoFemto_pT1_pT2->SetUp(1,1024,0,4096);
  GhettoFemto_pT1_pT2->Initialize();

  if(GhettoFemto_pT1_div_pT) delete GhettoFemto_pT1_div_pT;
  GhettoFemto_pT1_div_pT = new DLM_Histo<float>();
  GhettoFemto_pT1_div_pT->SetUp(1);
  GhettoFemto_pT1_div_pT->SetUp(0,1024,0,1);
  GhettoFemto_pT1_div_pT->Initialize();


  if(Ghetto_mT_mTwrong) delete Ghetto_mT_mTwrong;
  Ghetto_mT_mTwrong = new DLM_Histo<float>();
  Ghetto_mT_mTwrong->SetUp(2);
  if(Ghetto_MtBins){
    Ghetto_mT_mTwrong->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
    Ghetto_mT_mTwrong->SetUp(1,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    Ghetto_mT_mTwrong->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
    Ghetto_mT_mTwrong->SetUp(1,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  Ghetto_mT_mTwrong->Initialize();


  if(GhettoFemto_mT_mTwrong) delete GhettoFemto_mT_mTwrong;
  GhettoFemto_mT_mTwrong = new DLM_Histo<float>();
  GhettoFemto_mT_mTwrong->SetUp(2);
  if(Ghetto_MtBins){
    GhettoFemto_mT_mTwrong->SetUp(0,Ghetto_NumMtBins,Ghetto_MtBins);
    GhettoFemto_mT_mTwrong->SetUp(1,Ghetto_NumMtBins,Ghetto_MtBins);
  }
  else{
    GhettoFemto_mT_mTwrong->SetUp(0,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
    GhettoFemto_mT_mTwrong->SetUp(1,Ghetto_NumMtBins,Ghetto_MtMin,Ghetto_MtMax);
  }
  GhettoFemto_mT_mTwrong->Initialize();

}
