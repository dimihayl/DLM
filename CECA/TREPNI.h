//The Resonance Emission Particle Nomenclature In CECA

#ifndef TREPNI_H
#define TREPNI_H

#include <iostream>

class TREPNI;

class TreParticle{
friend class TREPNI;
private:
  TreParticle(TREPNI& database);
  ~TreParticle();
  const TREPNI& Database;
  char* TreName;
  float* Mass;
  float* Width;
  //in both cases, we assume that the streamer is set to the correct position
  void AppendInBinary(std::ofstream& file);
  void LoadFromBinary(std::ifstream& file);
};

class TreChain{
private:
  TreChain(TreParticle& mother);
  ~TreChain();
  const TreParticle& Mother;

};

class TREPNI{
public:
  TREPNI(const unsigned short& version);
  //for the future, when you LOAD from file
  //TREPNI(char* InputFile);
  ~TREPNI();


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
};

#endif // TREPNI_H
