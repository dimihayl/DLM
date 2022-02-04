
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
  decay = NULL;
}
CecaParticle::CecaParticle(const CecaParticle &other){
  //printf("copy\n");
  CecaParticle();
  *this=other;
}
CecaParticle::~CecaParticle(){
  //printf("del %p %p\n",cats,this);
  if(cats) {delete cats; cats = NULL;}
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
  *cats = *prt_cats;
}
void CecaParticle::SetDecay(const TreChain* prt_dec){
  decay = prt_dec;
}
void CecaParticle::RandomDecay(){
  if(trepni){
    decay = trepni->GetDecay();
  }
}
CecaParticle& CecaParticle::operator=(const CecaParticle& other){
  //printf("=\n");
  trepni = other.trepni;
  if(!cats) cats = new CatsParticle();
  *cats = *other.cats;
  decay = other.decay;
  return *this;
}


//! nothing done on the errors (single, level etc), do it afterwards
CECA::CECA(const TREPNI& database,const std::vector<std::string>& list_of_particles):
//CECA::CECA(const TREPNI& database):
  Database(database),MaxThreads(std::thread::hardware_concurrency()?std::thread::hardware_concurrency():1){

  Displacement = new float [3];
  DisplacementAlpha = new float [3];
  //Hadronization = new float [3];
  //HadronizationAlpha = new float [3];
  Hadronization = 0;
  HadronizationAlpha = 2;
  Tau = 0;
  ProperTau = true;
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
  //delete [] Hadronization;
  //delete [] HadronizationAlpha;
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
/*
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
*/
//identical X,Y,Z
void CECA::SetHadronization(const float& width, const float& levy){
  //SetHadronizationX(width,levy);
  //SetHadronizationY(width,levy);
  //SetHadronizationZ(width,levy);
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Hadronization = fabs(width);
  HadronizationAlpha = levy;
}

void CECA::SetTau(const float& tau, const bool& proper){
  if(tau<0){
    printf("ERROR tau\n");
    return;
  }
  Tau = tau;
  ProperTau = proper;
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

    unsigned ThId = omp_get_thread_num();

    //there was some issue using the objects
    //a silly workaround: I will only keep track of pointers,
    //however we will need to call delete for each object whenever required
    std::vector<CecaParticle*> Primordial;
    std::vector<CecaParticle*> Primary;

    //while we create an event containing all final state speciies of interest
    //N.B. so far no sanity check if this is even possible, i.e. an infinite loop is more than possible here!!!
    while(true){
      Primordial.clear();
      for(unsigned short uMult=0; uMult<EMULT; uMult++){
        //--- SELECT A PRIMORDIAL ---//
        Primordial.push_back(new CecaParticle());
        TreParticle* tre = Database.GetRandomParticle();
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
            Primordial.back()->RandomDecay();
            for(std::string& str : ListOfParticles){
              for(unsigned uDaught=0; uDaught<Primordial.back()->Decay()->GetNumDaughters(); uDaught++){
                FSP_is_product = ParticleInList(Primordial.back()->Decay()->GetDaughter(uDaught));
              }
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
      //for(const TreParticle* particle : Primordial){
      for(CecaParticle* particle : Primordial){
        for(int iNeeded=0; iNeeded<particle_list.size(); iNeeded++){
          //check if primary
          if(particle->Trepni()->GetName()==particle_list.at(iNeeded)){
            particle_list.erase(particle_list.begin()+iNeeded);
          }
          //check if decay product
          else if(particle->Decay()){
            unsigned char ndaugh = particle->Decay()->GetNumDaughters();
            for(unsigned char uDaugh=0; uDaugh<ndaugh; uDaugh++){
              std::string pname = particle->Decay()->GetDaughter(uDaugh)->GetName();
              if(pname==particle_list.at(iNeeded)){
                particle_list.erase(particle_list.begin()+iNeeded);
                if(!particle_list.size()) break;
              }
            }
          }
        }
      }
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
      double px,py,pz,ptot,sin_th,cos_th,cotg_th;

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
      }

      primordial->Cats()->SetMXYZ(primordial->Trepni()->GetMass(),px,py,pz);
      primordial->Cats()->SetWidth(primordial->Trepni()->GetWidth());
      primordial->Cats()->SetDecayRanGen(RanGen[ThId]);

      //--- EMISSION ---//
      //--- PROPAGATE BASED ON THE PROPERTIES OF THE CORE SOURCE ---//
      double rd[3], beta[3], rtot[3];
      double tau = Tau;
      double rh = RanGen[ThId]->StableR(3,HadronizationAlpha,0,Hadronization);
      if(ProperTau) tau *= primordial->Cats()->Gamma();
      //add the displacement and the beta*tau components
      for(int xyz=0; xyz<3; xyz++){
        beta[xyz] = primordial->Cats()->Beta(xyz);
        rd[xyz]  = RanGen[ThId]->Stable(DisplacementAlpha[xyz],0,Displacement[xyz]);
        rtot[xyz] = rd[xyz]+beta[xyz]*tau;
      }
      //we need to add the hadronization part separately, as we demand it to have
      //the same direction as the velocity, i.e. we need beta first
      double beta_tot = sqrt(beta[0]*beta[0]+beta[1]*beta[1]+beta[2]*beta[2]);
      for(int xyz=0; xyz<3; xyz++){
        if(beta_tot) rtot[xyz]+=beta[xyz]/beta_tot*rh;
      }
      //the final position is saved. The time corresponds to the time elapsed
      //in the laboratory
      primordial->Cats()->SetTXYZ(tau,rtot[0],rtot[1],rtot[2]);

      //--- DECAY OF RESONANCES + PROPAGATION OF THE MOTHERS ---/
      //if the width is zero, the Decay function returns the daughters
      bool IsResonance = true;
      IsResonance = !ParticleInList(primordial->Trepni());

      if(IsResonance){
        CatsParticle* daughters =
          primordial->Cats()->Decay(
          primordial->Decay()->GetNumDaughters(),
          primordial->Decay()->GetDaughterMasses(),
          true);
        for(unsigned char nd=0; nd<primordial->Decay()->GetNumDaughters(); nd++){
          if(ParticleInList(primordial->Decay()->GetDaughter(nd))){
            //clv_primary.push_back(daughters[nd]);
            //Primary.push_back(primordial.Decay()->GetDaughter(nd));
            Primary.push_back(new CecaParticle());
            Primary.back()->SetTrepni(primordial->Decay()->GetDaughter(nd));
            Primary.back()->SetDecay(NULL);
            Primary.back()->SetCats(daughters[nd]);
          }
        }
        delete [] daughters;
      }
      else{
        Primary.push_back(new CecaParticle());
        *Primary.back() = *primordial;
      }
    }//iteration over all primordials

    //BUILD THE MULTIPLETS, EVALUATE THEIR NUMBER AND RETURN THE CORRECT VALUE

    //the array position of each particle, which is to be used to build the multiplet
    //the length SDIM represents the number of particles in each multiplet
    //e.g. 5 particles, SDIM=3 has to build all permutations: 012,013,014,023,024,034,123,124,134,234
    std::vector<std::vector<unsigned>> Permutations = BinomialPermutations(Primary.size(),SDIM);
    //the pid is a single permutation, e.g. {0,1,2}
    for(std::vector<unsigned> pid : Permutations){
//GHETTO: make multiplets and simply drop the output for the source as a function of rstar, no kstar, no shit
//this so that you can show something next FemTUM
      //these are all particles we need to include in a multiplet
      printf("\n");
      CatsLorentzVector boost_v;

      CatsLorentzVector* prt_cm = new CatsLorentzVector[SDIM];

      unsigned char ud=0;
      for(unsigned ID : pid){
        boost_v = boost_v+*(Primary.at(ID)->Cats());//
        prt_cm[ud] = *Primary.at(ID)->Cats();
        prt_cm[ud].Print();
        ud++;
      }
      printf("-boost-\n");boost_v.Print();printf("-------\n");
      for(unsigned char ud=0; ud<SDIM; ud++){
        prt_cm[ud].Boost(boost_v);
        prt_cm[ud].Print();
      }

      CatsLorentzVector cm_sumQA;
      for(unsigned char ud=0; ud<SDIM; ud++){
        cm_sumQA = cm_sumQA+prt_cm[ud];
      }
      printf("-boosted sum-\n");cm_sumQA.Print();printf("-------\n");
      usleep(1000e3);
/*
NEXT_STEPS
the tau correction, based on largest tau, and than build up R and Q, plot R for Q<FemtoLimit.
test for two particles first!!!
so plot rstar for femto pairs and see how it looks
also plot the angles relevant for epos comparison
*/
      delete [] prt_cm;

      //Primary
      //CatsMultiplet Multiplet(*RanGen[ThId],SDIM,false);

    }

    return 0;

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
