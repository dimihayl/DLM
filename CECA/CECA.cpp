
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
void CecaParticle::RandomDecay(){
  if(trepni){
    decay = trepni->GetDecay();
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

  Displacement = new float [3];
  DisplacementAlpha = new float [3];
  Hadronization = new float [3];
  HadronizationAlpha = new float [3];
  for(int i=0; i<3; i++){
    Displacement[i]=0;
    DisplacementAlpha[i]=2;
    Hadronization[i]=0;
    HadronizationAlpha[i]=2;
  }
  HadrFluct = 0;
  //Hadronization = 0;
  //HadronizationAlpha = 2;
  Tau = 0;
  TauFluctuation = 0;
  ProperTau = true;
  ThermalKick = 0;
  SDIM = 2;
  TargetYield = 100000;
  AchievedYield = 0;
  FemtoLimit = 200;
  UpperLimit = 300;
  EMULT = 0;
  SrcCnv = 1;
  DebugMode = false;
  //CLV.clear();

  ///////////////////////////////////////////////
  Ghetto_rstar = NULL;
  Ghetto_rcore = NULL;
  GhettOld_rstar = NULL;
  Old_rstar = NULL;
  Old_rcore = NULL;
  Old_source = NULL;
  Old_CosRcP1 = NULL;
  Old_CosRcP2 = NULL;
  Old_CosP1P2 = NULL;
  Old_RcP1 = NULL;
  Old_RcP2 = NULL;
  Old_P1P2 = NULL;
  Ghetto_kstar = NULL;
  Ghetto_kstar_rstar = NULL;
  Ghetto_mT_rstar = NULL;
  GhettoFemto_rstar = NULL;
  Ghetto_mT_costh = NULL;
  GhettoSP_pT_th = NULL;
  Ghetto_RP_AngleRcP1 = NULL;
  Ghetto_PR_AngleRcP2 = NULL;
  Ghetto_RR_AngleRcP1 = NULL;
  Ghetto_RR_AngleRcP2 = NULL;
  Ghetto_RR_AngleP1P2 = NULL;
  GhettoFemto_mT_rstar = NULL;
  GhettoFemto_mT_kstar = NULL;
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
  if(Ghetto_mT_rstar){delete Ghetto_mT_rstar; Ghetto_mT_rstar=NULL;}
  if(GhettoFemto_rstar){delete GhettoFemto_rstar; GhettoFemto_rstar=NULL;}
  if(Ghetto_mT_costh){delete Ghetto_mT_costh; Ghetto_mT_costh=NULL;}
  if(GhettoSP_pT_th){delete GhettoSP_pT_th; GhettoSP_pT_th=NULL;}
  if(Ghetto_RP_AngleRcP1){delete Ghetto_RP_AngleRcP1; Ghetto_RP_AngleRcP1=NULL;}
  if(Ghetto_PR_AngleRcP2){delete Ghetto_PR_AngleRcP2; Ghetto_PR_AngleRcP2=NULL;}
  if(Ghetto_RR_AngleRcP1){delete Ghetto_RR_AngleRcP1; Ghetto_RR_AngleRcP1=NULL;}
  if(Ghetto_RR_AngleRcP2){delete Ghetto_RR_AngleRcP2; Ghetto_RR_AngleRcP2=NULL;}
  if(Ghetto_RR_AngleP1P2){delete Ghetto_RR_AngleP1P2; Ghetto_RR_AngleP1P2=NULL;}
  if(GhettoFemto_mT_rstar){delete GhettoFemto_mT_rstar; GhettoFemto_mT_rstar=NULL;}
  if(GhettoFemto_mT_kstar){delete GhettoFemto_mT_kstar; GhettoFemto_mT_kstar=NULL;}
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

void CECA::SetDisplacementY(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[1] = fabs(width);
  DisplacementAlpha[1] = levy;
}

void CECA::SetDisplacementZ(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[2] = fabs(width);
  DisplacementAlpha[2] = levy;
}

//identical X,Y
void CECA::SetDisplacementT(const float& width, const float& levy){
  SetDisplacementX(width,levy);
  SetDisplacementY(width,levy);
}

//identical X,Y,Z
void CECA::SetDisplacement(const float& width, const float& levy){
  SetDisplacementX(width,levy);
  SetDisplacementY(width,levy);
  SetDisplacementZ(width,levy);
}

void CECA::SetHadronizationX(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[0] = width;
  HadronizationAlpha[0] = levy;
}

void CECA::SetHadronizationY(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[1] = width;
  HadronizationAlpha[1] = levy;
}

void CECA::SetHadronizationZ(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization[2] = width;
  HadronizationAlpha[2] = levy;
}

//identical X,Y
void CECA::SetHadronizationT(const float& width, const float& levy){
  SetHadronizationX(width,levy);
  SetHadronizationY(width,levy);
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

void CECA::SetHadrFluctuation(const float& fluct){
  HadrFluct = fluct;
}

void CECA::SetTau(const float& tau, const float& fluct, const bool& proper){
  if(tau<0||fluct<0){
    printf("ERROR tau\n");
    return;
  }
  Tau = tau;
  TauFluctuation = fluct;
  ProperTau = proper;
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

void CECA::SetTargetStatistics(const unsigned& yield){
  TargetYield = yield;
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
//void CECA::SetSeed(const unsigned& seed){
//  RanGen->SetSeed(seed);
//}

//returns the number of generated multiplets
unsigned CECA::GoSingleCore(const unsigned& ThId){
  ThreadClock[ThId].Start();
  unsigned ExeTime;
  unsigned NumMultiplets = 0;
  do{
    if(DebugMode){
      //printf("Single core event generator\n");
      //printf("-> NumMult = %i; ExeTime = %u\n",NumMultiplets,ExeTime);
    }
    NumMultiplets += GenerateEvent();
    ExeTime = unsigned(ThreadClock[ThId].Stop()/(long long)(1000000));
    //if(DebugMode) printf("\n");
  }
  while(ExeTime<Timeout);
  return NumMultiplets;
}

void CECA::OptimizeThreadCount(){

}

void CECA::SaveBuffer(){

}

void CECA::GoBabyGo(const unsigned& num_threads){
  if(!Database.QA()){
      return;
  }
  bool DynamicThreads;
  if(num_threads==0){
    NumThreads = omp_get_num_threads();
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

  //CLV.clear();
  //for(unsigned uThread=0; uThread<NumThreads; uThread++){
  //  std::vector<CatsLorentzVector> vtemp;
  //  CLV.push_back(vtemp);
  //}

  if(DebugMode){
    printf("Running GoBabyGo\n");
    printf(" Detected threads: %u\n",NumThreads);
  }

  unsigned* BufferYield = new unsigned [NumThreads];
  AchievedYield = 0;
  //const unsigned TargetPerSyst = TargetYield / NumSystVars;
  //unsigned AchievedSystYield = 0;
  //we iterate until we have our target yield
  //this refers to the global yield, but also the minimum yiled per systematics
  //the latter is done to minimize the bias within the yield of each syst. variation
  //while(AchievedYield<TargetYield || AchievedSystYield<TargetPerSyst){
  while(AchievedYield<TargetYield){
    //we run each thread for a maximum of some preset amount of time
    #pragma omp parallel for
    for(unsigned uThread=0; uThread<NumThreads; uThread++){
      BufferYield[uThread] = GoSingleCore(uThread);
    }
    for(unsigned uThread=0; uThread<NumThreads; uThread++){
      AchievedYield += BufferYield[uThread];
      //AchievedSystYield += BufferYield[uThread];
    }
    if(DebugMode){
      printf(" Achieved/Target Yield = %u / %u\n",AchievedYield,TargetYield);
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

//generates all particles, propagates and decays into the particles of interest
//and lastly builds up the
unsigned CECA::GenerateEvent(){
//printf("\n<<<Welcome to the generator>>>\n");
    unsigned ThId = omp_get_thread_num();

    //there was some issue using the objects
    //a silly workaround: I will only keep track of pointers,
    //however we will need to call delete for each object whenever required
    std::vector<CecaParticle*> Primordial;
    std::vector<CecaParticle*> Primary;

    //while we create an event containing all final state speciies of interest
    //N.B. so far no sanity check if this is even possible, i.e. an infinite loop is more than possible here!!!
    while(true){
      //printf("Lets give it a shot\n");
      Primordial.clear();
      //printf("A fresh start\n");
      for(unsigned short uMult=0; uMult<EMULT; uMult++){
        //--- SELECT A PRIMORDIAL ---//
        Primordial.push_back(new CecaParticle());
        TreParticle* tre = Database.GetRandomParticle();
        Primordial.back()->SetTrepni(tre);
//std::cout<< Primordial.back()->Trepni()->GetName() << std::endl;

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
            Primordial.back()->RandomDecay();
            //for(std::string& str : ListOfParticles){
              for(unsigned uDaught=0; uDaught<Primordial.back()->Decay()->GetNumDaughters(); uDaught++){
                FSP_is_product = ParticleInList(Primordial.back()->Decay()->GetDaughter(uDaught));
                if(FSP_is_product) break;
              }
              //if(FSP_is_product) break;
            //}
          }
        }
        else{
          Primordial.back()->SetDecay(NULL);
        }
        bool Useless_particle = (!FSP_is_primordial&&!FSP_is_product);
//printf("UP=%i\n",Useless_particle);
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
//usleep(2000e3);
      /*
      //for(const TreParticle* particle : Primordial){
      printf("Едно две три ... %lu %lu\n",particle_list.size(),ListOfParticles.size());
      for(CecaParticle* particle : Primordial){
        for(int iNeeded=particle_list.size()-1; iNeeded>=0; iNeeded--){
          printf("On the lookout for %s\n",particle_list.at(iNeeded).c_str());
          //check if primordial
          if(particle->Trepni()->GetName()==particle_list.at(iNeeded)){
            printf("Found a primordial %s\n",particle->Trepni()->GetName().c_str());
            particle_list.erase(particle_list.begin()+iNeeded);
            //iNeeded--; if(iNeeded<0) break;
            printf(" - removed it from the list\n");
          }
          //check if decay product
          else if(particle->Decay()){
            printf("Found a decay %s\n",particle->Trepni()->GetName().c_str());
            unsigned char ndaugh = particle->Decay()->GetNumDaughters();
            for(unsigned char uDaugh=0; uDaugh<ndaugh; uDaugh++){
              std::string pname = particle->Decay()->GetDaughter(uDaugh)->GetName();
              if(pname==particle_list.at(iNeeded)){
                printf("Found a primary %s\n",pname.c_str());
                particle_list.erase(particle_list.begin()+iNeeded);
                printf(" - removed it from the list\n");
                if(!particle_list.size()) break;
              }
              else{
                printf("A rejected %s\n",pname.c_str());
              }
            }
          }
        }
      }
      */

      //printf("ела ме настигни ... \n");
      //we found all particles
      if(particle_list.size()==0)
        break;
    }//the inifinite while loop
//printf("Out of the infinite loop\n");
    for(CecaParticle* primordial : Primordial){
      //--- SAMPLE THE MOMENTUM ---//
      double axisValues[3];
      double& pT = axisValues[0];
      double& eta = axisValues[1];
      double& phi = axisValues[2];
      double px,py,pz,ptot,sin_th,cos_th,cotg_th,cos_phi,sin_phi;

      //if we sample from a predefined PDF
      if(primordial->Trepni()->GetMomPDF()){
        bool Auto_eta = true;
        bool Auto_phi = true;
        //Auto_pT = false;
        if(primordial->Trepni()->GetMomPDF()->GetDim()>1)
          Auto_eta = false;
        if(primordial->Trepni()->GetMomPDF()->GetDim()>2)
          Auto_phi = false;
        primordial->Trepni()->GetMomPDF()->Sample(axisValues,true);

        if(Auto_phi){
          phi = RanGen[ThId]->Uniform(0,2.*Pi);
        }
        if(Auto_eta){
          cos_th = RanGen[ThId]->Uniform(-1.,1.);
          //sin always positive here
          sin_th = sqrt(1.-cos_th*cos_th);
          cotg_th = cos_th/sin_th;
          //verified that this relation is true
          //eta = -0.5*log((1.-cos_th)/(1.+cos_th));
        }
        else{
          sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
          cotg_th = (1.-exp(-2.*eta))/(2.*exp(-eta));
          cos_th = (1.-exp(-2.*eta)/(1.+exp(-2.*eta)));
        }
        pz = pT*cotg_th;
        ptot = sqrt(pT*pT+pz*pz);
        px = ptot*cos(phi)*sin_th;
        py = ptot*sin(phi)*sin_th;

      }
      //automatic sampling of all components, following Gaussian x,y,z of certain width
      //if nothing is set, the width is assumed zero by default, so no motion at all
      else{
        primordial->Trepni()->SampleMomXYZ(axisValues,true);
        px = axisValues[0];
        py = axisValues[1];
        pz = axisValues[2];
        //pT = sqrt(px*px+py*py+pz*pz);
        ptot = sqrt(px*px+py*py+pz*pz);
        cos_th = pz/ptot;
        sin_th = sqrt(1.-cos_th*cos_th);
      }
      cos_phi = px/(ptot*sin_th);
      sin_phi = py/(ptot*sin_th);

      primordial->Cats()->SetMXYZ(primordial->Trepni()->GetMass(),px,py,pz);
      primordial->Cats()->SetWidth(primordial->Trepni()->GetWidth());
      primordial->Cats()->SetDecayRanGen(RanGen[ThId]);

      //printf("\n");
      //primordial->Cats()->Print();

      //--- EMISSION ---//
      //--- PROPAGATE BASED ON THE PROPERTIES OF THE CORE SOURCE ---//
      double rd[3],beta[3],rtot[3],rh[3],mom[3];
      double energy;
      double tau = Tau;
      //double rh[3];
      //standard/original procedure
      //double rh_len = RanGen[ThId]->StableR(3,HadronizationAlpha[0],0,Hadronization[0]);
      //the model where we assume the source is an ellipsoid around the displacement point,
      //and that direction of velocity is what determines the "crossing point" of the particle with the emission source
      double rh_len=-1;

      int ResampleCount = 1000;
      while(true){
        for(int xyz=0; xyz<3; xyz++){
          //rh[xyz] = RanGen[ThId]->Gauss(0,Hadronization[xyz]);
          //rh[xyz] = Hadronization[xyz];
          //rh[xyz] = RanGen[ThId]->Gauss(Hadronization[xyz],Displacement[xyz]);
          //rh[xyz] = RanGen[ThId]->Gauss(Hadronization[xyz],0);
          //rh[xyz] = Hadronization[xyz];

          rd[xyz] = RanGen[ThId]->Gauss(0,Displacement[xyz]);
          rh[xyz] = RanGen[ThId]->Gauss(Hadronization[xyz],Hadronization[xyz]*HadrFluct);
        }
        //this comes from the definition of an ellipsoid
        rh_len = sqrt(pow(rh[0]*sin_th*cos_phi,2.)+pow(rh[1]*sin_th*sin_phi,2.)+pow(rh[2]*cos_th,2.));
        //rh_len = RanGen[ThId]->Gauss(rh_len,Displacement[2]*rh_len);
        tau += RanGen[ThId]->Gauss(0,TauFluctuation);

        if(ProperTau) tau *= primordial->Cats()->Gamma();
        //add the displacement and the beta*tau components
        energy=primordial->Cats()->Mag2();
        for(int xyz=0; xyz<3; xyz++){
          //thermal kick
          mom[xyz] = RanGen[ThId]->Gauss(primordial->Cats()->GetP(xyz),ThermalKick);
          energy += mom[xyz]*mom[xyz];
        }
        energy = sqrt(energy);
        for(int xyz=0; xyz<3; xyz++){
          //beta[xyz] = primordial->Cats()->Beta(xyz);
          //mom[xyz] = primordial->Cats()->GetP(xyz);
          //beta[xyz] = RanGen[ThId]->Gauss(beta[xyz],0.075);
          beta[xyz] = mom[xyz]/energy;
          //printf("beta %f --> %f\n",primordial->Cats()->Beta(xyz),beta[xyz]);
          //usleep(200e3);

          //rd[xyz]  = RanGen[ThId]->Stable(DisplacementAlpha[xyz],0,Displacement[xyz]);
          rtot[xyz] = rd[xyz]+beta[xyz]*tau;
          //rtot[xyz] = beta[xyz]*tau;
          //rh[xyz] = RanGen[ThId]->Stable(HadronizationAlpha[xyz],0,Hadronization[xyz]);
        }
        //we need to add the hadronization part separately, as we demand it to have
        //the same direction as the velocity, i.e. we need beta first
        double beta_tot = sqrt(beta[0]*beta[0]+beta[1]*beta[1]+beta[2]*beta[2]);
//if(beta_tot>1)beta_tot=1;
        for(int xyz=0; xyz<3; xyz++){
          if(beta_tot) rtot[xyz]+=beta[xyz]/beta_tot*rh_len;
        }

        //the final position is saved. The time corresponds to the time elapsed
        //in the laboratory
        primordial->Cats()->SetTXYZ(tau,rtot[0],rtot[1],rtot[2]);

        //if(ResampleCount<1000){
        //  printf("  Resampled\n");
        //  printf("   t,x,y,z (l) = %4.2f %4.2f %4.2f %4.2f (%4.2f)\n",tau,rtot[0],rtot[1],rtot[2],rh_len);
        //  usleep(1000e3);
        //}

        double p_tot,p_x,p_y,p_z;
        double dr_tot,dr_x,dr_y,dr_z;
        double LorentzWeight[2];
        double LorentzSize[2];
        double Size[2];
        double Slope[2];
        CecaParticle* prim[2];
        prim[0] = primordial;
        //the probability to accept this position sampling
        double RejectProb = 0;
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
            //if(LorentzWeight[ip]<0 || LorentzWeight[ip]>1){
              //printf(" BIG BUG IN CECA::GenerateEvent --> LorentzWeight = %e\n",LorentzWeight[ip]);
            //}
            double EffectiveSize = sqrt(lw2*ls2+(1.-lw2)*sz2);
            double FD = exp((dr_tot-EffectiveSize)/(EffectiveSize*Slope[ip]))+1.;
            RejectProb /= FD?FD:1;
            //if(dr_tot<0.5){
            //  printf("   p_tot = %5.0f; dr_tot = %5.2f; rp = %6.4f; lw = %5.3f; sz = %4.2f; esz = %4.2f\n",p_tot,dr_tot,RejectProb,LorentzWeight[ip],Size[ip],EffectiveSize);
            //  usleep(500e3);
            //}
          }//ip
        }//iter over primoridal2
//RejectProb=0;
//printf("RP = %f\n",RejectProb);
        if(RanGen[ThId]->Uniform(0,1)>=RejectProb){
          break;
        }
        else{
          //printf(" REJETED\n");
          //printf("  t,x,y,z (l) = %4.2f %4.2f %4.2f %4.2f (%4.2f)\n",tau,rtot[0],rtot[1],rtot[2],rh_len);
          //usleep(1000e3);
        }
        if(ResampleCount<0){
          printf("\033[1;33mWARNING:\033[0m CECA::GenerateEvent says that it cannot separate the particles at the set requirement.\n");
          printf("   To solve the issue, verify there is enought displacement and the particle radius is not too large.\n");
          break;
        }
        ResampleCount--;
      }//while(rh_len<0)

      //--- DECAY OF RESONANCES + PROPAGATION OF THE MOTHERS ---/
      //if the width is zero, the Decay function returns the daughters
      bool IsResonance = true;
      IsResonance = !ParticleInList(primordial->Trepni());
      //printf("IR %i\n",IsResonance);
//usleep(500e3);
      if(IsResonance){
        //primordial->Cats()->Print();
        CatsParticle* daughters =
          primordial->Cats()->Decay(
          primordial->Decay()->GetNumDaughters(),
          primordial->Decay()->GetDaughterMasses(),
          true);
        //primordial->Cats()->Print();
        //printf(" is reso\n");
        for(unsigned char nd=0; nd<primordial->Decay()->GetNumDaughters(); nd++){
          if(ParticleInList(primordial->Decay()->GetDaughter(nd))){
            //printf(" daughter found\n");
            //clv_primary.push_back(daughters[nd]);
            //Primary.push_back(primordial.Decay()->GetDaughter(nd));
            Primary.push_back(new CecaParticle());
            Primary.back()->SetTrepni(primordial->Decay()->GetDaughter(nd));
            Primary.back()->SetDecay(NULL);
            Primary.back()->SetCats(daughters[nd]);
            Primary.back()->SetOrigin(2);
            Primary.back()->SetMother(primordial->Cats());
            //Primary.back()->Cats()->Print();
            //printf("\n\n");
            //usleep(500e3);
          }
        }
        delete [] daughters;
      }
      else{
        Primary.push_back(new CecaParticle());
        //printf("crash?\n");
        *Primary.back() = *primordial;
      }
      //printf("IR END\n");

      if(!GhettoSP_pT_th){
        GhettoSP_pT_th = new DLM_Histo<float>();
        GhettoSP_pT_th->SetUp(2);
        GhettoSP_pT_th->SetUp(0,64,0,4096);
        GhettoSP_pT_th->SetUp(1,64,-0.1,3.2);
        GhettoSP_pT_th->Initialize();
      }
      GhettoSP_pT_th->Fill(Primary.back()->Cats()->GetPt(),Primary.back()->Cats()->GetTheta());
    }//iteration over all primordials
//printf("iter over primordials over\n");
    //BUILD THE MULTIPLETS, EVALUATE THEIR NUMBER AND RETURN THE CORRECT VALUE

    //the array position of each particle, which is to be used to build the multiplet
    //the length SDIM represents the number of particles in each multiplet
    //e.g. 5 particles, SDIM=3 has to build all permutations: 012,013,014,023,024,034,123,124,134,234
    std::vector<std::vector<unsigned>> Permutations = BinomialPermutations(Primary.size(),SDIM);
    //printf("PERM %lu %lu\n",Primary.size(),Permutations.size());
    //the pid is a single permutation, e.g. {0,1,2}
    for(std::vector<unsigned> pid : Permutations){
//GHETTO: make multiplets and simply drop the output for the source as a function of rstar, no kstar, no shit
//this so that you can show something next FemTUM
      //these are all particles we need to include in a multiplet
      //printf("\n");
      CatsLorentzVector boost_v;

      //CatsLorentzVector* prt_cm = new CatsLorentzVector[SDIM];
      CecaParticle* prt_cm = new CecaParticle[SDIM];

      //bool* is_promordial = new bool[SDIM];
      //CatsLorentzVector* mother_cm = new CatsLorentzVector[SDIM];

      unsigned char ud=0;
      std::vector<float> cos_th;
      for(unsigned ID : pid){
        boost_v = boost_v+*(Primary.at(ID)->Cats());//
        prt_cm[ud] = *Primary.at(ID);
        //prt_cm[ud].SetCats(Primary.at(ID)->Cats());
        //is_promordial[ud] = Primary.at(ID)->IsUsefulPrimordial();
        if(Primary.at(ID)->IsUseful()==false){printf("How did this happen!?!?!\n");}
        cos_th.push_back(cos(prt_cm[ud].Cats()->GetTheta()));
        //prt_cm[ud].Print();
        ud++;
      }
      //printf("-boost-\n");boost_v.Print();printf("-------\n");
      for(unsigned char ud=0; ud<SDIM; ud++){
        prt_cm[ud].Cats()->Boost(boost_v);
        if(prt_cm[ud].Mother()) prt_cm[ud].Mother()->Boost(boost_v);
        //prt_cm[ud].Print();
      }

      CatsLorentzVector cm_sumQA;
      for(unsigned char ud=0; ud<SDIM; ud++){
        cm_sumQA = cm_sumQA+*prt_cm[ud].Cats();
      }
      //printf("-boosted sum-\n");cm_sumQA.Print();printf("-------\n");
      //usleep(1000e3);
/*
NEXT_STEPS
the tau correction, based on largest tau, and than build up R and Q, plot R for Q<FemtoLimit.
test for two particles first!!!
so plot rstar for femto pairs and see how it looks
also plot the angles relevant for epos comparison
*/

//GHETTO, works for pairs only
if(SDIM==2){
//printf("INSIDE THE GHETTO\n");
CatsLorentzVector cm_rel = *prt_cm[1].Cats()-*prt_cm[0].Cats();
double drx,dry,drz;
drx = cm_rel.GetX();
dry = cm_rel.GetY();
drz = cm_rel.GetZ();

double kstar = cm_rel.GetP();
double rstar = cm_rel.GetR();
double mT = 0.5*boost_v.GetMt();

double AngleP1P2=0;
double AngleRcP1=0;
double AngleRcP2=0;


if(!Ghetto_RP_AngleRcP1){
  Ghetto_RP_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RP_AngleRcP1->SetUp(1);
  Ghetto_RP_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RP_AngleRcP1->Initialize();
}
if(!Ghetto_RP_AngleRcP1){
  Ghetto_RP_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RP_AngleRcP1->SetUp(1);
  Ghetto_RP_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RP_AngleRcP1->Initialize();
}
if(!Ghetto_PR_AngleRcP2){
  Ghetto_PR_AngleRcP2 = new DLM_Histo<float>();
  Ghetto_PR_AngleRcP2->SetUp(1);
  Ghetto_PR_AngleRcP2->SetUp(0,128,0,Pi);
  Ghetto_PR_AngleRcP2->Initialize();
}
if(!Ghetto_RR_AngleRcP1){
  Ghetto_RR_AngleRcP1 = new DLM_Histo<float>();
  Ghetto_RR_AngleRcP1->SetUp(1);
  Ghetto_RR_AngleRcP1->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleRcP1->Initialize();
}
if(!Ghetto_RR_AngleRcP2){
  Ghetto_RR_AngleRcP2 = new DLM_Histo<float>();
  Ghetto_RR_AngleRcP2->SetUp(1);
  Ghetto_RR_AngleRcP2->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleRcP2->Initialize();
}
if(!Ghetto_RR_AngleP1P2){
  Ghetto_RR_AngleP1P2 = new DLM_Histo<float>();
  Ghetto_RR_AngleP1P2->SetUp(1);
  Ghetto_RR_AngleP1P2->SetUp(0,128,0,Pi);
  Ghetto_RR_AngleP1P2->Initialize();
}

if(kstar<300){

  if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulPrimordial()){

  }
  else if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulPrimordial()){
    double cosine = cm_rel.GetX()*prt_cm[0].Mother()->GetPx()+
                    cm_rel.GetY()*prt_cm[0].Mother()->GetPy()+
                    cm_rel.GetZ()*prt_cm[0].Mother()->GetPz();
    cosine /= (cm_rel.GetR()*prt_cm[0].Mother()->GetP());
    AngleRcP1 = acos(cosine);
    Ghetto_RP_AngleRcP1->Fill(AngleRcP1);
  }
  else if(prt_cm[0].IsUsefulPrimordial()&&prt_cm[1].IsUsefulProduct()){
    double cosine = cm_rel.GetX()*prt_cm[1].Mother()->GetPx()+
                    cm_rel.GetY()*prt_cm[1].Mother()->GetPy()+
                    cm_rel.GetZ()*prt_cm[1].Mother()->GetPz();
    cosine /= (cm_rel.GetR()*prt_cm[1].Mother()->GetP());
    AngleRcP2 = acos(cosine);
    Ghetto_PR_AngleRcP2->Fill(AngleRcP2);
  }
  else if(prt_cm[0].IsUsefulProduct()&&prt_cm[1].IsUsefulProduct()){
    double cosine = cm_rel.GetX()*prt_cm[0].Mother()->GetPx()+
                    cm_rel.GetY()*prt_cm[0].Mother()->GetPy()+
                    cm_rel.GetZ()*prt_cm[0].Mother()->GetPz();
    cosine /= (cm_rel.GetR()*prt_cm[0].Mother()->GetP());
    AngleRcP1 = acos(cosine);

    cosine =        cm_rel.GetX()*prt_cm[1].Mother()->GetPx()+
                    cm_rel.GetY()*prt_cm[1].Mother()->GetPy()+
                    cm_rel.GetZ()*prt_cm[1].Mother()->GetPz();
    cosine /= (cm_rel.GetR()*prt_cm[1].Mother()->GetP());
    AngleRcP2 = acos(cosine);

    cosine =  prt_cm[0].Mother()->GetPx()*prt_cm[1].Mother()->GetPx()+
              prt_cm[0].Mother()->GetPy()*prt_cm[1].Mother()->GetPy()+
              prt_cm[0].Mother()->GetPz()*prt_cm[1].Mother()->GetPz();
    cosine /= (prt_cm[0].Mother()->GetP()*prt_cm[1].Mother()->GetP());//
    AngleP1P2 = acos(cosine);

    Ghetto_RR_AngleRcP1->Fill(AngleRcP1);
    Ghetto_RR_AngleRcP2->Fill(AngleRcP2);
    Ghetto_RR_AngleP1P2->Fill(AngleP1P2);
  }

}


const double NumMomBins = 32;
const double MomMin = 0;
const double MomMax = 2048;

const double NumRadBins = 128;
const double RadMin = 0;
const double RadMax = 8;

const double NumMtBins = 32;
const double MtMin = 0;
const double MtMax = 4096;

if(!Ghetto_rstar){
  Ghetto_rstar = new DLM_Histo<float>();
  Ghetto_rstar->SetUp(1);
  Ghetto_rstar->SetUp(0,NumRadBins,RadMin,RadMax);
  Ghetto_rstar->Initialize();
}
if(!Ghetto_kstar){
  Ghetto_kstar = new DLM_Histo<float>();
  Ghetto_kstar->SetUp(1);
  Ghetto_kstar->SetUp(0,NumMomBins,MomMin,MomMax);
  Ghetto_kstar->Initialize();
}
if(!Ghetto_kstar_rstar){
  Ghetto_kstar_rstar = new DLM_Histo<float>();
  Ghetto_kstar_rstar->SetUp(2);
  Ghetto_kstar_rstar->SetUp(0,NumMomBins,MomMin,MomMax);
  Ghetto_kstar_rstar->SetUp(1,NumRadBins,RadMin,RadMax);
  Ghetto_kstar_rstar->Initialize();
}

if(!Ghetto_mT_rstar){
  Ghetto_mT_rstar = new DLM_Histo<float>();
  Ghetto_mT_rstar->SetUp(2);
  Ghetto_mT_rstar->SetUp(0,NumMtBins,MtMin,MtMax);
  Ghetto_mT_rstar->SetUp(1,NumRadBins,RadMin,RadMax);
  Ghetto_mT_rstar->Initialize();
}

if(!GhettoFemto_rstar){
  GhettoFemto_rstar = new DLM_Histo<float>();
  GhettoFemto_rstar->SetUp(1);
  GhettoFemto_rstar->SetUp(0,NumRadBins,RadMin,RadMax);
  GhettoFemto_rstar->Initialize();
}

if(!Ghetto_mT_costh){
  Ghetto_mT_costh = new DLM_Histo<float>();
  Ghetto_mT_costh->SetUp(2);
  Ghetto_mT_costh->SetUp(0,NumMtBins,MtMin,MtMax);
  Ghetto_mT_costh->SetUp(1,128,-1,1);
  Ghetto_mT_costh->Initialize();

}

if(!GhettoFemto_mT_rstar){
  GhettoFemto_mT_rstar = new DLM_Histo<float>();
  GhettoFemto_mT_rstar->SetUp(2);
  GhettoFemto_mT_rstar->SetUp(0,NumMtBins,MtMin,MtMax);
  GhettoFemto_mT_rstar->SetUp(1,NumRadBins,RadMin,RadMax);
  GhettoFemto_mT_rstar->Initialize();
}
if(!GhettoFemto_mT_kstar){
  GhettoFemto_mT_kstar = new DLM_Histo<float>();
  GhettoFemto_mT_kstar->SetUp(2);
  GhettoFemto_mT_kstar->SetUp(0,NumMtBins,MtMin,MtMax);
  GhettoFemto_mT_kstar->SetUp(1,NumMomBins,MomMin,MomMax);
  GhettoFemto_mT_kstar->Initialize();
}


//if(kstar<200){
  //printf("\nA pair of\n");
  //printf(" kstar = %.2f\n", kstar);
  //printf(" rstar = %.2f\n", rstar);
  //printf(" mT = %.2f\n", mT);
  //printf("-------------------\n");
  //usleep(500e3);
  Ghetto_rstar->Fill(rstar);
  Ghetto_kstar->Fill(kstar);
  Ghetto_kstar_rstar->Fill(kstar,rstar);
  Ghetto_mT_rstar->Fill(mT,rstar);
  for(float ct : cos_th){
    Ghetto_mT_costh->Fill(mT,ct);
  }
  GhettoFemto_mT_kstar->Fill(mT,kstar);
//printf("%i%i\n",is_promordial[0],is_promordial[1]);
//usleep(200e3);
  GhettoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
  GhettoPrimReso[1] += (prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
  GhettoPrimReso[2] += (!prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
  GhettoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());

  if(kstar<150){
    GhettoFemto_rstar->Fill(rstar);
    GhettoFemto_mT_rstar->Fill(mT,rstar);
    //printf(" mT = %.2f\n", mT);
    GhettoFemtoPrimReso[0] += (prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
    GhettoFemtoPrimReso[1] += (prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
    GhettoFemtoPrimReso[2] += (!prt_cm[0].IsUsefulPrimordial() && prt_cm[1].IsUsefulPrimordial());
    GhettoFemtoPrimReso[3] += (!prt_cm[0].IsUsefulPrimordial() && !prt_cm[1].IsUsefulPrimordial());
  }
//}
}
//printf("OUT OF THE GHETTO\n");
///////////////////////////


      delete [] prt_cm;
      //delete [] is_promordial;
      //Primary
      //CatsMultiplet Multiplet(*RanGen[ThId],SDIM,false);
//printf("ALL BAD MEMORIES ARE GONE\n");
    }//permutations over possible multiplets
//printf("We want to delete shit\n");
  for(CecaParticle* particle : Primordial){
    delete particle;
  }
//printf("We did so for the primordials\n");
  for(CecaParticle* particle : Primary){
    delete particle;
  }
//printf("We did so for the primaries\n");
//printf("We will return %lu\n",Permutations.size());
  return Permutations.size();
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
