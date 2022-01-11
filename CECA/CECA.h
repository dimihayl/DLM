//Common Emission in CAts

#ifndef CECA_H
#define CECA_H

#include <vector>

template <class Type> class DLM_Histo;
class CatsLorentzVector;
class CatsParticle;
class TREPNI;
class DLM_CleverMcLevyResoTM;
class DLM_Timer;

//we can create one CECA only with a specific database, that is not allowed to be modified
class CECA{
public:
  CECA(const TREPNI& database);
  ~CECA();

  //the width (X), alpha factor (Levy) related to the displacement source
  //the displ. is relative to the "perfect" collision point
  void SetDisplacementX(const float& width, const float& levy=2);//done
  void SetDisplacementY(const float& width, const float& levy=2);//done
  void SetDisplacementZ(const float& width, const float& levy=2);//done
  //identical X,Y
  void SetDisplacementT(const float& width, const float& levy=2);//done
  //identical X,Y,Z
  void SetDisplacement(const float& width, const float& levy=2);//done

  void SetHadronizationX(const float& width, const float& levy=2);//done
  void SetHadronizationY(const float& width, const float& levy=2);//done
  void SetHadronizationZ(const float& width, const float& levy=2);//done
  //identical X,Y
  void SetHadronizationT(const float& width, const float& levy=2);//done
  //identical X,Y,Z
  void SetHadronization(const float& width, const float& levy=2);//done

  void SetTau(const float& tau);//done

  //all source up to DIM will be evaluated
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
  void SetSystVars(const unsigned& howmany);//done

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
  void SetThreadTimeout(const unsigned& seconds);
  //void SetThreadTimeout(const unsigned& seconds);

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

  DLM_CleverMcLevyResoTM* Old_source;

private:
  const TREPNI& Database;

  short SourceConvention;
  //XYZ
  float* Displacement;
  float* DisplacementAlpha;
  float* Hadronization;
  float* HadronizationAlpha;
  float Tau;
  unsigned char SDIM;
//bug prone: if this is smaller then SDIM, we are up for trouble!
  unsigned short EMULT;
  unsigned TargetYield;
  unsigned AchievedYield;
  char SrcCnv;
  bool DebugMode;

  //the k* (MeV) below which a FemtoPair is concidered such
  //N.B. for more particles QN = sqrt(N*kstarlimit)
  float FemtoLimit;
  //the upper limit refers to the femto limit, and information about
  //all particle pairs is still saved, e.g. for systematics
  float UpperLimit;
  //the number of desired systematic variations
  unsigned NumSystVars;

  //how many particles to simulate
  //it would be nice to have some basic auto-sampling to estimate how many
  //pairs are needed to reach a certain amount of pairs, triplets, or whatever
  unsigned NumParticles;//this is not used so far...

  void GenerateEvent();
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



};

#endif // CECA_H
