
#include "TREPNI.h"

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "DLM_CppTools.h"

TreParticle::TreParticle(TREPNI& database):Database(database){
  
}





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
