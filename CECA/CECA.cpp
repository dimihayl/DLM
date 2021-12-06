
#include "CATSconstants.h"
#include "CECA.h"
#include "TREPNI.h"
#include "CATStools.h"
#include "DLM_Random.h"
#include "DLM_Histo.h"
#include "DLM_Source.h"
#include "DLM_CppTools.h"

#include "omp.h"
#include <unistd.h>
#include <thread>

//! nothing done on the errors (single, level etc), do it afterwards
CECA::CECA(const TREPNI& database):Database(database),MaxThreads(std::thread::hardware_concurrency()?std::thread::hardware_concurrency():1){
  Displacement = new float [3];
  DisplacementAlpha = new float [3];
  Hadronization = new float [3];
  HadronizationAlpha = new float [3];
  Tau = 0;
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
  if(ThreadClock){delete [] ThreadClock; ThreadClock=NULL;}
}

void CECA::SetDisplacementX(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[0] = width;
  DisplacementAlpha[0] = levy;
}

void CECA::SetDisplacementY(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[1] = width;
  DisplacementAlpha[1] = levy;
}

void CECA::SetDisplacementZ(const float& width, const float& levy){
  if(levy<1||levy>2){
    printf("ERROR levy\n");
    return;
  }
  Displacement[2] = width;
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
}

void CECA::SetTau(const float& tau){
  if(tau<0){
    printf("ERROR tau\n");
    return;
  }
  Tau = tau;
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

void CECA::GoSingleCore(const unsigned& ThId){
  ThreadClock[ThId].Start();
  unsigned ExeTime;
  do{
    GenerateEvent();

    ExeTime = unsigned(ThreadClock[ThId].Stop()/(long long)(1000000));
  }
  while(ExeTime<Timeout);
}

void CECA::OptimizeThreadCount(){

}

void CECA::SaveBuffer(){

}

void CECA::GoBabyGo(const unsigned& num_threads){
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

  AchievedYield = 0;
  //we iterate until we have our target yield
  while(AchievedYield<TargetYield){
    //we run each thread for a maximum of some preset amount of time
    #pragma omp parallel for
    for(unsigned uThread=0; uThread<NumThreads; uThread++){
      GoSingleCore(uThread);
    }
    //after the timeout, we optimize the thread count and, in case
    //TargetYield is not achieved yet, than we continue
    if(DynamicThreads) OptimizeThreadCount();
    //saves permanantly the output that was in the buffer of each thread
    SaveBuffer();
  }

}

//generates all particles, propagates and decays into the particles of interest
//and lastly builds up the
void CECA::GenerateEvent(){
  double axis_values[3];
  for(unsigned short uMult=0; uMult<EMULT; uMult++){
    //Get a random particle specie
    //TreParticle* Particle = Database.GetRandomParticle();
    //TreChain* Particle->GetRandomDecay();


    //TEMP
    //Test with pp interaction, where you only have 3 particles in the database:
    //proton, Reso, pion
    //hardcode all the fractions here, just try to simulate the propagation and
    //final 2-particle source to see if it works
    //TreParticle* Particle = Database.GetParticle();
    //TreChain* Particle->GetRandomDecay();

  }
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
