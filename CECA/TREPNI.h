//The Resonance Emission Particle Nomenclature In CECA

#ifndef TREPNI_H
#define TREPNI_H

#include <iostream>

class TREPNI;
class TreChain;
template <class Type> class DLM_Histo;
class CatsLorentzVector;
class DLM_Random;

//the idea of this class is to only be used by TREPNI. Thus, there are NO
//set functions available for the outside world. Else we might end up f***ing up
//the whole TREPNI database

//so, TREPNI should own a list of particles,
//each particle should own a list of decay chains,
//and each decay chain should contain a list of linked daughter particles
class TreParticle{
friend class TREPNI;
public:
  std::string GetName() const;//done
  float GetMass() const;//done
  float GetMassLow() const;//done
  float GetMassUp() const;//done
  float GetWidth() const;//done
  float GetWidthLow() const;//done
  float GetWidthUp() const;//done
  float GetAbundance() const;//done
  float GetAbundanceLow() const;//done
  float GetAbundanceUp() const;//done
  unsigned char GetNumDecays() const;//done
  const TREPNI* GetDatabase() const;//done

  //these we can make private, ones the testing is done
  TreParticle(TREPNI& database);//done
  ~TreParticle();//done

  //N.B. the pdf is copied into TREPNI, which sometimes might take a while,
  //so do not over-use this function. The decision to copy the pdf and not point to it
  //is related to the fact that ultimately these particles will be saved into a file,
  //so we better own all the information straight away.
  void SetMomPDF(const DLM_Histo<float>& pdf);//done, the others below might not be needed
  void FillMomXYZ(const float& xval, const float& yval, const float& zval);
  void FillMomPtEtaPhi(const float& pt, const float& eta, const float& phi);
  void FillMomPDF(CatsLorentzVector& cats_vector);
  const DLM_Histo<float>* GetMomPDF() const;

  //the QA of these setters is done from TREPNI
  void SetName(const char* name);//done
  void SetMass(const float& mass);//done
  void SetMassLimit(const float& mass_low, const float& mass_up);//done
  void SetWidth(const float& width);//done
  void SetWidthLimit(const float& width_low, const float& width_up);//done
  void SetAbundance(const float& abundance);//done W/O QA !!!
  void SetAbundanceLimit(const float& abundance_low, const float& abundance_up);//done W/O QA !!!
  TreChain* NewDecay();//done W/O QA !!!
  TreChain* GetDecay(const unsigned char& whichone);//done W/O QA !!!
  void Print();

  //randomize all
  void Randomize();
  void RandomizeMass();
  void RandomizeWidth();
  void RandomizeAbundance();
  void RandomizeBR();
  ////////////////////////////////////////////////////////

private:

  const TREPNI& Database;
  //it is a 3D disto, given in terms of pT eta and phi
  //by convention, the dimensions are [0] = pT, [1] = eta, [2] = phi
  //the DIM can be reduced, i.e. 1D or 2D histos will be accepted, where if that is the case
  //for 2D we will only have pT,eta (phi assumed flat)
  //while for 1D we will only have pT (cos theta and phi both assumed flat)
  const DLM_Histo<float>* MomPDF;
  char* TreName;
  float* Mass;
  float* Width;
  float* Abundance;
  unsigned char NumDecays;
  //important: the decay chains will be owned by the particle
  TreChain** Decay;

  //in both cases, we assume that the streamer is set to the correct position
  void AppendInBinary(std::ofstream& file);
  void LoadFromBinary(std::ifstream& file);

};

//a decay chain
//just as for TreParticle, only TREPNI can be setting up this class,
//no outside user interaction will be provided!
class TreChain{
friend class TREPNI;
friend class TreParticle;
public:
  //check if the mass/widths are set up fine to make possible for this decay
  //bool QA() const;

//these we can make private, ones the testing is done
//void SetDaughters(const unsigned char& numdaughters, const TreParticle* daughter);
void AddDaughter(const TreParticle& daughter);//done W/O QA !!!
std::string GetName();//done W/O QA !!!
void SetBranching(const float& br);//done W/O QA !!!
void SetBranchingLimit(const float& br_low, const float& br_up);//done W/O QA !!!
void RandomizeBR();
////////////////////////////////////////////////////////

private:
  TreChain(TreParticle& mother);//done
  ~TreChain();//done
  unsigned char NumDaughters;
  const TreParticle& Mother;
  //in percent
  float* Branching;
  //important: the daughters are only addresses, and NOT owned by the decay chain
  const TreParticle** Daughter;

};

class TREPNI{
friend class TreParticle;
public:
  TREPNI(const unsigned short& version);
  //for the future, when you LOAD from file
  //TREPNI(char* InputFile);
  ~TREPNI();

//for the QA, perhaps introduce the required mother!!!
//any decay chain should end with the mother!!!
  enum QA_type { Full, Name, Daughters, Mass, Width, BR, Abundance };
  bool QA(const int& type=0);//done
  void SetTotalYield(const float& totyield);//done
  //n.b. if the TotalYield is not fixed, the return value is dynamically
  //evaluated by summing up all yields of all particles
  float GetYield() const;//done
  TreParticle* NewParticle(const char* name=NULL);//done W/O QA !!!
  TreParticle* GetParticle(const unsigned& whichone) const;//done W/O QA !!!
  TreParticle* GetParticle(const char* name) const;//done
  //based on abundance
  TreParticle* GetRandomParticle() const;
  unsigned GetNumParticles();//done

  //this randomizes all particles for:
  //all properties, or just mass,width,abundance and their BRs
  void Randomize();
  void RandomizeMass();
  void RandomizeWidth();
  void RandomizeAbundance();
  void RandomizeBR();
  void SetSeed(const unsigned& seed);//done
  //0-3, 3 is the most
  void SetPrintLevel(const char& lvl, const bool& single=true);

private:
  //an enum containing an ID of each function. To keep track, each time you
  //define a new function, add it here
  enum FunName {  constructor_0, destructor_0, getparticlename, getparticleid, setprintlevel,
                  setparticle, deleteparticle, setparticlemass, setpartiglegamma, memorymanager };
  //max length of the database name (string)
  const unsigned Len_DtbsName;
  //max length of the particle name (string)
  const unsigned Len_PrtclName;
  //version of this database
  const unsigned short Version;
  //the allocated mem will be doubled each time it is not sufficient
  //however, we have a limit on this doubling procedure, e.g. dont create more than 1024 new slots
  //the length corresponds to the number of particles within
  const unsigned MaxMemSteps;
  //number of functions for which we can save error information (see ErrorOccured)
  const short NumFunctions;
  //this will be important when creating the binary file
  //if false, you are not allowed to create one!
  //i.e. the data base needs to set up properly
  //it should also work in the other direction, if it is false when reading from a binary,
  //e.g. created externaly via web interface, you will not be allowed to load the code!
  bool QA_passed;

  DLM_Random* RanGen;

  //if this paramter is >0, than the QA on the Abundance is enforced to
  //make sure the total Abundance sums up to this number. Otherwise the
  //Abundance is treated as a paramter with arbitrary units and the total Abundance
  //is NOT forced to be a specific number
  float TotAbundance;

  //0 = no output
  //1 = only errors
  //2 = err + warning
  //3 = err + warning + hint
  char PrintLevel;
  //if true, each error is displayed only one time
  bool SingleError;
  //here we will use a bit mask. Each bit will be treated as a bool for
  //if an error occured within a certain function. I.e. per function we
  //will have a miximum of 32 error messages
  int* ErrorOccured;

  char* DatabaseName;


  //void SetParticle();

  //no two particles with the same name
  bool QA_Name();//done
  //checks if each decay channel has at least two daughters
  bool QA_Daughters();//done
  //the MEAN mass of the mother should be LARGER than the
  //sum of the MEAN masses of all daughters
  bool QA_Mass();//done
  bool QA_Width();//done
  //the sum of all BRs (all decays per particle) sould be 100%, thus we damand
  //that 100% should be within 68% central interval of the uncertainties,
  //i.e. we have high probability of sampling meaningful BRs
  bool QA_BR();//done
  //the sum of the abundancies (of all particles) should sum up to TotAbundance, thus we demand
  //that 100% should be within 68% central interval of the uncertainties,
  //i.e. we have high probability of sampling meaningful abundancies
  //NOT relevent if TotAbundance<=0
  bool QA_Abundance();//done

  unsigned NumParticles;
  //for the mem allocation
  unsigned MaxParticles;
  //important: the particles are owned by the database
  TreParticle** Particle;
//make a checksum function to add all abundancies

};


class TREPNI_old{
public:
  TREPNI_old(const unsigned short& version);
  //for the future, when you LOAD from file
  //TREPNI(char* InputFile);
  ~TREPNI_old();

private:



/*
  std::string GetParticleName(const int& id) const;//done
  //the counting starts from 1, minus sign is anti particle
  int GetParticleId(const char* name) const;//done

  //0-3, 3 is the most
  void SetPrintLevel(const char& lvl, const bool& single=true);//done
  //if the particle name exists, it rewrites it. Otherwise it makes a new entry
  void SetParticle(const char* name, const double& mass_min, const double& mass_max,
                            const double& gamma_min, const double& gamma_max);
  void DeleteParticle(const char* name);
  void SetParticleMass(const char* name, const double& mass);
  void SetParticleGamma(const char* name, const double& gamma);



private:
  //an enum containing an ID of each function. To keep track, each time you
  //define a new function, add it here
  enum FunName {  constructor_0, destructor_0, getparticlename, getparticleid, setprintlevel,
                  setparticle, deleteparticle, setparticlemass, setpartiglegamma, memorymanager};
  //max length of the database name (string)
  const unsigned Len_DtbsName;
  //max length of the particle name (string)
  const unsigned Len_PrtclName;
  //version of this database
  const unsigned short Version;
  //the allocated mem will be doubled each time it is not sufficient
  //however, we have a limit on this doubling procedure, e.g. dont create more than 1024 new slots
  //the length corresponds to the number of particles within
  const unsigned MaxMemSteps;
  //number of functions for which we can save error information (see ErrorOccured)
  const short NumFunctions;
  const unsigned char MaxDecayCh;
  const unsigned char MaxDaughters;
  //0 = no output
  //1 = only errors
  //2 = err + warning
  //3 = err + warning + hint
  char PrintLevel;
  //if true, each error is displayed only one time
  bool SingleError;
  //here we will use a bit mask. Each bit will be treated as a bool for
  //if an error occured within a certain function. I.e. per function we
  //will have a miximum of 32 error messages
  int* ErrorOccured;

  char* DatabaseName;
  //unsigned char MaxDecayCh;

  unsigned NumParticles;
  //for the mem allocation
  unsigned MaxParticles;
  //properties of each particle
  //N.B. NEVER call these objects directly, but use the Get functions
  //Why? Because of special conventions, to keep it all close by in memory
  //the Mass and Gamma have 3*NumParticles elements, where each two neighboring
  //elements correspond to the minimum/currently used/maximum value of M or G
  //this comes into play for the uncertainties
  //int* TrepnID;
  //name of the particle
  char** TrepName;
  float* Mass;
  float* Gamma;
  int* Nch;

  TreParticle* Particle;


  //for each decay channel of each particle:
  //[Particle][DecayChannel]
  int** Ndaughter;
  //this also has 3 times the elements, to accomodate min/current/max
  float** Branching;
  //[Particle][DecayChannel][Daughter]
  int*** DaughterID;

  void MemoryManager(const bool& destroy=false);
  */
};

#endif // TREPNI_H
