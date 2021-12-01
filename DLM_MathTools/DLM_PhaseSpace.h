#ifndef DLM_PHASESPACE_H
#define DLM_PHASESPACE_H

//#include <vector>

template <class Type> class DLM_Histo;
class DLM_Random;


class DLM_DecayGen{
public:
  DLM_DecayGen(const unsigned& seed=0);
  ~DLM_DecayGen();

  //generate a single decay
  void GenerateEvent();
  //generate the kinematic distribution of the decay
  void GenerateDist(const unsigned& numdecays);
  //my first implementation, that needed to be set up within few hours
  //hopefully will be obsolete soon, in favour of the above functions
  void DirtyDist(const unsigned& numdecays);

  //some decays are saved as a seed for good solutions. It is needed
  //to make the generation of new events faster.
  void SetBufferDepth(const unsigned& depth);

  //N.B. These are setter functions, cannot be called after an event has
  //already been generated
  void SetMother(const double& mass, const double& px=0, const double& py=0, const double& pz=0);
  //the code will not run for few than 2 daughters
  //at the moment the maximum # daughters is only 4, due to CPU issues
  //if time allows, I will work on improving this for the future
  void AddDaughter(const double& mass);// ADD QA !!!
  void SetEpsilon(const double& epsilon, const unsigned& steps);// ADD QA !!!
  //silly events are those that do not respect energy conservation
  //by default they are not allowed
  //still, one can include them in the analysis, if they deviate by less
  //than nthr of Epsilon, i.e. something like nsigma
  void SetSillyThreshold(const double& nthr);

private:
  const unsigned SEED;
  const unsigned MaxNumDaughters;
  const unsigned MinNumDaughters;
  unsigned BufferDepth;
  unsigned NumDaughters;
  double MotherMass;
  double* DaughterMass;
  double EpsilonStart;
  unsigned EpsSteps;
  double SillyThreshold;
  std::vector DaughterMass;
  DLM_Random* RanGen;
  //DLM_Histo<double>* dResult;

};

void TwoBodyDecay(){

}

#endif
