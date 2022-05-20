//Common Emission in CAts

#ifndef CECA_H
#define CECA_H

#include <vector>

template <class Type> class DLM_Histo;
class CatsLorentzVector;
class CatsParticle;
class TREPNI;
class TreParticle;
class TreChain;
class DLM_CleverMcLevyResoTM;
class DLM_Timer;
class DLM_Random;

class CecaParticle{
public:
  CecaParticle();
  CecaParticle(const CecaParticle &other);
  ~CecaParticle();
  const TreParticle* Trepni() const;
  CatsParticle* Cats() const;
  CatsParticle* Mother() const;

  const TreChain* Decay() const;
  void SetTrepni(const TreParticle& prt_tre);
  void SetCats(const CatsParticle& prt_cats);
  void SetDecay(const TreChain& prt_dec);
  void RandomDecay(DLM_Random* RanGen=NULL);
  void SetTrepni(const TreParticle* prt_tre);
  void SetCats(const CatsParticle* prt_cats);
  void SetDecay(const TreChain* prt_dec);
  //0 - useless, 1 - primordial, 2 - decay product
  void SetOrigin(const char& origin);
  void SetMother(const CatsParticle* mama);
  bool IsUseful() const;//is of the required type
  bool IsUsefulPrimordial() const;//+primordial
  bool IsUsefulProduct() const;//+decay product
  bool WithinAcceptance() const;
  CecaParticle& operator=(const CecaParticle& other);
private:
  //n.b. the CatsParticle is owned!!!
  const TreParticle* trepni;
  CatsParticle* cats;
  CatsParticle* mother;
  const TreChain* decay;
  //1 - primordial of the required type
  //2 - decay product of the required type
  //-1 - primordial of non-required type
  //-2 - decay product of non-required type
  char Origin;
};

//we can create one CECA only with a specific database, that is not allowed to be modified
class CECA{
public:
  //list_of_particles contains the names of the particles separated by commas and/or spaces
  CECA(const TREPNI& database, const std::vector<std::string>& list_of_particles);
  //CECA(const TREPNI& database);
  ~CECA();

  //the width (X), alpha factor (Levy) related to the displacement source
  //the displ. is relative to the "perfect" collision point
  void SetDisplacementX(const float& width, const float& levy=2);//done
  float GetDisplacementX() const;//done
  void SetDisplacementY(const float& width, const float& levy=2);//done
  float GetDisplacementY() const;//done
  void SetDisplacementZ(const float& width, const float& levy=2);//done
  float GetDisplacementZ() const;//done
  //identical X,Y
  void SetDisplacementT(const float& width, const float& levy=2);//done
  float GetDisplacementT() const;//done
  //identical X,Y,Z
  void SetDisplacement(const float& width, const float& levy=2);//done
  float GetDisplacement() const;//done

  void SetHadronizationX(const float& width, const float& levy=2);//done
  float GetHadronizationX() const;//done
  void SetHadronizationY(const float& width, const float& levy=2);//done
  float GetHadronizationY() const;//done
  void SetHadronizationZ(const float& width, const float& levy=2);//done
  float GetHadronizationZ() const;//done
  //identical X,Y
  void SetHadronizationT(const float& width, const float& levy=2);//done
  float GetHadronizationT() const;//done
  //identical X,Y,Z
  void SetHadronization(const float& width, const float& levy=2);//done
  float GetHadronization() const;//done
  void SetHadrFluctuation(const float& fluct);//done no QA
  float GetHadrFluctuation() const;//done
  //this should be false, however during testing it was found out that the
  //EPOS model is doing exactly that (propagating the primordial resonances to the decay point)
  //this makes the angles that we extract to equal those between rstar and XX, not rcore and XX
  //This function was introduced only with the intent of testing these effects!!!
  void SetPropagateMother(const bool& yesno);//done

  //if proper == true, it means the time is defined in the rest frame
  //of the particle (i.e. property of the particle)
  //if false, the time is let to run in the LAB, i.e. we treat this parameter
  //as a property of the system itself
  void SetTau(const float& tau, const float& fluct=0, const bool& proper=true);//done
  float GetTau() const;
  void SetThermalKick(const float& kick);//done w/o qa

  //all source up to DIM will be evaluated
//THIS SHOULD BE FIX BY THE PARTICLE LIST
  void SetSourceDim(const unsigned char& sdim);//done
  //the yield of multipletes in the highest dimension
  //e.g. if DIM is set to 3, this will be the target for 3-body yield
  //the target is evaluated for the femto region
  void SetTargetStatistics(const unsigned& yield);//done
  //femto is the k* for the target statistics.
  //info the upper limit until which information is still saved
  //N.B. the corresponding Q3 etc is automatically evaluated, here its only k*!
  void SetFemtoRegion(const float& femto, const float& info=0);//done
  //higher number should make for less CPU time
  //questions about accuracy though, so far no effect seen
  //still, if Mult = Dim is the configuration that will for sure
  //give the correct Dim-body source. Usefull parameter to QA
  //if zero, than by default this is set to Dim*Dim
  void SetEventMult(const unsigned short& emult=0);//done

  //the number of desired systematic variations
  //dont go overboard, as generating a new variation within the database
  //takes some time. The idea is that we equally split our TargetYield
  //amoung each (random) systematic variation, done automatically in GoBabyGo
  //void SetSystVars(const unsigned& howmany);//done

  //0 = single particle (SP)
  //1 = nolan (NL)
  //2 = gaussian pair (GP)
  //the default is 1, to keep is on same foot with the two-body analysis
  //this is very confusing, but in short:
  //The convention is related to the source size (R).
  //SP means that the width of the SINGLE particle source distribution is R in each (X,Y,Z) direction
  //GP is a convention where the pair source extends by R in each X,Y,Z
  //Nolan is somewhere in-between, and has been historically used for femto...
  //the relations amoung those are the following:
  //R_SP = 2/sqrt(alpha)*R_NL = sqrt(alpha)R_GP
  //alpha is the Levy factor
  void SetSourceConvention(const char& srccnv);//done
  void SetDebugMode(const bool& debugmode);//done
  //the time that each thread will run on its own, before communicating with the master class.
  //This makes so that the AchievedYield will only be checked ones in this time period, i.e.
  //the final TargetYield can be exceeded. Longer timeout maximizes single thread throughput,
  //however the overall CPU throughput is optimized in steps of timeout, hence it is not good to go crazy long.
  //also, the GoBabyGo will run for a minimum time set by the timeout.
  //Unless there is a very good reason, do not change the default value!
  //The time is given in seconds (integer!)
  void SetThreadTimeout(const unsigned& seconds);//done
  void SetSeed(const unsigned& thread, const unsigned& seed);//done
  //void SetThreadTimeout(const unsigned& seconds);
  //true by default. Within a multiplet it propagates the particles until
  //they all have the same tau component
  void EqualizeFsiTime(const bool& yesno);//done

  //if num_threads=0, we will dynamically adjust based on efficiency
//check if EMULT>=SDIM on start!!!
  void GoBabyGo(const unsigned& num_threads=0);

  void GhettoTest1(const unsigned NumPairs, const float r_SP, const float p_SP);
  //CatsParticle* GhettoDecay(CatsParticle& particle, const float mass1, const float mass2);
  DLM_Histo<float>* Ghetto_rstar;
  DLM_Histo<float>* Ghetto_rcore;
  DLM_Histo<float>* GhettOld_rstar;
  DLM_Histo<float>* Old_rstar;
  DLM_Histo<float>* Old_rcore;
  DLM_Histo<float>* Old_CosRcP1;
  DLM_Histo<float>* Old_CosRcP2;
  DLM_Histo<float>* Old_CosP1P2;
  DLM_Histo<float>* Old_RcP1;
  DLM_Histo<float>* Old_RcP2;
  DLM_Histo<float>* Old_P1P2;

  DLM_Histo<float>* Ghetto_kstar;
  DLM_Histo<float>* Ghetto_kstar_rstar;
  DLM_Histo<float>* Ghetto_mT_rstar;
  DLM_Histo<float>* GhettoFemto_rstar;
  DLM_Histo<float>* GhettoFemto_rcore;
  DLM_Histo<float>* GhettoFemto_mT_rstar;
  DLM_Histo<float>* GhettoFemto_mT_rcore;
  DLM_Histo<float>* GhettoFemto_mT_kstar;
  DLM_Histo<float>* Ghetto_mT_costh;
  DLM_Histo<float>* GhettoSP_pT_th;
  DLM_Histo<float>* GhettoSP_pT_1;
  DLM_Histo<float>* GhettoSP_pT_2;

  DLM_Histo<float>* Ghetto_PP_AngleRcP1;
  DLM_Histo<float>* Ghetto_PP_AngleRcP2;
  DLM_Histo<float>* Ghetto_PP_AngleP1P2;

  DLM_Histo<float>* Ghetto_RP_AngleRcP1;
  DLM_Histo<float>* Ghetto_PR_AngleRcP2;

  DLM_Histo<float>* Ghetto_RR_AngleRcP1;
  DLM_Histo<float>* Ghetto_RR_AngleRcP2;
  DLM_Histo<float>* Ghetto_RR_AngleP1P2;

  DLM_Histo<float>* Ghetto_ScatteringAngle;

  //pp / pr / rp / rr
  unsigned GhettoFemtoPrimReso[4];
  unsigned GhettoPrimReso[4];
  //unsigned GhettoSpPrim[2];


  DLM_CleverMcLevyResoTM* Old_source;
  DLM_CleverMcLevyResoTM* SetUp_RSM;
  DLM_CleverMcLevyResoTM* SetUp_RSM_UNI;
  DLM_CleverMcLevyResoTM* SetUp_RSM_BB;
  std::vector<float*>* Buffer_RSM;

double GHETTO_ResoAbundance[2];
double GHETTO_ResoMass[2];
double GHETTO_ResoTau[2];
double GHETTO_DaughterMass[2][2];

bool GHETTO_EVENT=false;
private:
  const TREPNI& Database;
  std::vector<std::string> ListOfParticles;

  short SourceConvention;
  //XYZ
  float* Displacement;
  float* DisplacementAlpha;
  float* Hadronization;
  float* HadronizationAlpha;
  float HadrFluct;
  float Tau;
  float TauFluctuation;
  bool ProperTau;
  bool EqualFsiTau;
  float ThermalKick;
  bool PropagateMother;
//THIS SHOULD BE FIX BY THE PARTICLE LIST
  unsigned char SDIM;
//bug prone: if this is smaller then SDIM, we are up for trouble!
  unsigned short EMULT;
  unsigned TargetYield;
  unsigned AchievedYield;
  char SrcCnv;
  bool DebugMode;
  DLM_Random** RanGen;

  //the k* (MeV) below which a FemtoPair is concidered such
  //N.B. for more particles QN = sqrt(N*kstarlimit)
  float FemtoLimit;
  //the upper limit refers to the femto limit, and information about
  //all particle pairs is still saved, e.g. for systematics
  float UpperLimit;
  //the number of desired systematic variations
  //unsigned NumSystVars;

  //how many particles to simulate
  //it would be nice to have some basic auto-sampling to estimate how many
  //pairs are needed to reach a certain amount of pairs, triplets, or whatever
  unsigned NumParticles;//this is not used so far...

  //CatsMultiplet

  bool ParticleInList(const std::string& name) const;
  bool ParticleInList(const TreParticle* prt) const;
  unsigned GenerateEvent();
unsigned GenerateEventTEMP();
  //returns the number of generated multiplets
  unsigned GoSingleCore(const unsigned& ThId);
  void OptimizeThreadCount();
  void SaveBuffer();


  //variables to be used for the multi-threading
  //set up at GoBabyGo
  unsigned NumThreads;
  const unsigned MaxThreads;
  //all stuff used by each thread
  //[thread][mult]
  //inside we save, per thread, all particles that we generated
  //std::vector<std::vector<CatsLorentzVector>> CLV;

  std::vector<CatsLorentzVector>* ParticleBuffer;


  DLM_Timer* ThreadClock;
  //in seconds
  unsigned Timeout;
  //

  void GhettoInit();

};



#endif // CECA_H
