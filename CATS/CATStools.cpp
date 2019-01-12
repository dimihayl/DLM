#include "CATStools.h"

#include <iostream>
#include <math.h>

#include "gsl_sf_gamma.h"

#include "DLM_CppTools.h"
#include "CATSconstants.h"
#include "CATS.h"

//!Needed only for testing (contains usleep)
//#include <unistd.h>

using namespace std;

CatsLorentzVector::CatsLorentzVector(){
    FourSpace[0]=0;
    FourSpace[1]=0;
    FourSpace[2]=0;
    FourSpace[3]=0;
    FourMomentum[0]=0;
    FourMomentum[0]=0;
    FourMomentum[0]=0;
    FourMomentum[0]=0;
    Length=0;
    Length2=0;
    TotMom=0;
    TotMom2=0;
    Magnitude=0;
    Magnitude2=0;
    gamma=0;
    betaX=0;
    betaY=0;
    betaZ=0;
    beta=0;
}
CatsLorentzVector::~CatsLorentzVector(){

}

double CatsLorentzVector::GetT() const{
    return FourSpace[0];
}
double CatsLorentzVector::GetX() const{
    return FourSpace[1];
}
double CatsLorentzVector::GetY() const{
    return FourSpace[2];
}
double CatsLorentzVector::GetZ() const{
    return FourSpace[3];
}
double CatsLorentzVector::GetE() const{
    return FourMomentum[0];
}
double CatsLorentzVector::GetPx() const{
    return FourMomentum[1];
}
double CatsLorentzVector::GetPy() const{
    return FourMomentum[2];
}
double CatsLorentzVector::GetPz() const{
    return FourMomentum[3];
}
void CatsLorentzVector::Set(const double& tCrd, const double& xCrd, const double& yCrd, const double& zCrd,
             const double& engy, const double& xMom, const double& yMom, const double& zMom){
    FourSpace[0]=tCrd;
    FourSpace[1]=xCrd;
    FourSpace[2]=yCrd;
    FourSpace[3]=zCrd;
    FourMomentum[0]=engy;
    FourMomentum[1]=xMom;
    FourMomentum[2]=yMom;
    FourMomentum[3]=zMom;
    ComputeBetaGamma();
}
void CatsLorentzVector::RotateMomPhi(const double& angle){
    double XNEW = FourMomentum[1]*cos(angle)-FourMomentum[2]*sin(angle);
    double YNEW = FourMomentum[2]*cos(angle)+FourMomentum[1]*sin(angle);
    FourMomentum[1] =  XNEW;
    FourMomentum[2] =  YNEW;
    ComputeBetaGamma();
}

void CatsLorentzVector::RenormSpacialCoordinates(const double& Renorm){
    FourSpace[0]*=Renorm;
    FourSpace[1]*=Renorm;
    FourSpace[2]*=Renorm;
    FourSpace[3]*=Renorm;
}

CatsLorentzVector const CatsLorentzVector::operator+(const CatsLorentzVector& other){
    CatsLorentzVector RESULT;
    RESULT.FourSpace[0]=FourSpace[0]+other.FourSpace[0];
    RESULT.FourSpace[1]=FourSpace[1]+other.FourSpace[1];
    RESULT.FourSpace[2]=FourSpace[2]+other.FourSpace[2];
    RESULT.FourSpace[3]=FourSpace[3]+other.FourSpace[3];

    RESULT.FourMomentum[0]=FourMomentum[0]+other.FourMomentum[0];
    RESULT.FourMomentum[1]=FourMomentum[1]+other.FourMomentum[1];
    RESULT.FourMomentum[2]=FourMomentum[2]+other.FourMomentum[2];
    RESULT.FourMomentum[3]=FourMomentum[3]+other.FourMomentum[3];

    RESULT.ComputeBetaGamma();

    return RESULT;
}

CatsLorentzVector const CatsLorentzVector::operator-(const CatsLorentzVector& other){
    CatsLorentzVector RESULT;
    RESULT.FourSpace[0]=FourSpace[0]-other.FourSpace[0];
    RESULT.FourSpace[1]=FourSpace[1]-other.FourSpace[1];
    RESULT.FourSpace[2]=FourSpace[2]-other.FourSpace[2];
    RESULT.FourSpace[3]=FourSpace[3]-other.FourSpace[3];

    RESULT.FourMomentum[0]=FourMomentum[0]-other.FourMomentum[0];
    RESULT.FourMomentum[1]=FourMomentum[1]-other.FourMomentum[1];
    RESULT.FourMomentum[2]=FourMomentum[2]-other.FourMomentum[2];
    RESULT.FourMomentum[3]=FourMomentum[3]-other.FourMomentum[3];

    RESULT.ComputeBetaGamma();

    return RESULT;
}

void CatsLorentzVector::operator=(const CatsLorentzVector& other){
    FourSpace[0]=other.FourSpace[0];
    FourSpace[1]=other.FourSpace[1];
    FourSpace[2]=other.FourSpace[2];
    FourSpace[3]=other.FourSpace[3];

    FourMomentum[0]=other.FourMomentum[0];
    FourMomentum[1]=other.FourMomentum[1];
    FourMomentum[2]=other.FourMomentum[2];
    FourMomentum[3]=other.FourMomentum[3];

    Length=other.Length;
    Length2=other.Length2;
    TotMom=other.TotMom;
    TotMom2=other.TotMom2;
    Magnitude=other.Magnitude;
    Magnitude2=other.Magnitude2;
    gamma=other.gamma;
    betaX=other.betaX;
    betaY=other.betaY;
    betaZ=other.betaZ;
    beta=other.beta;
}

void CatsLorentzVector::Boost(const CatsLorentzVector& boostVec){
    Boost(boostVec, FourSpace, FourSpace);
    Boost(boostVec, FourMomentum, FourMomentum);
    ComputeBetaGamma();
}
void CatsLorentzVector::Boost(const CatsLorentzVector& boostVec, const double* InVec, double* OutVec){
    double GammaMomBeta = boostVec.gamma*(InVec[1]*boostVec.betaX + InVec[2]*boostVec.betaY + InVec[3]*boostVec.betaZ);
    double GammaVec0 = boostVec.gamma*InVec[0];
    double GammaDevGammaPlusOne = boostVec.gamma/(boostVec.gamma+1);

    OutVec[1] += boostVec.betaX*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[2] += boostVec.betaY*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[3] += boostVec.betaZ*(GammaDevGammaPlusOne*GammaMomBeta - GammaVec0);
    OutVec[0] = GammaVec0 - GammaMomBeta;
}
void CatsLorentzVector::ComputeBetaGamma(){
    Length2 = FourSpace[1]*FourSpace[1]+FourSpace[2]*FourSpace[2]+FourSpace[3]*FourSpace[3];
    Length = sqrt(Length2);
    TotMom2 = FourMomentum[1]*FourMomentum[1]+FourMomentum[2]*FourMomentum[2]+FourMomentum[3]*FourMomentum[3];
    TotMom = sqrt(TotMom2);
    Magnitude2 = FourMomentum[0]*FourMomentum[0]-TotMom2;
    Magnitude = sqrt(Magnitude2);
    betaX = FourMomentum[1]/FourMomentum[0];
    betaY = FourMomentum[2]/FourMomentum[0];
    betaZ = FourMomentum[3]/FourMomentum[0];
    beta = TotMom/FourMomentum[0];
    gamma = FourMomentum[0]/Magnitude;
}

CatsLorentzVector CatsLorentzVector::GetBoost(const CatsLorentzVector& boostVec){
    CatsLorentzVector RESULT=(*this);
    Boost(boostVec, FourSpace, RESULT.FourSpace);
    Boost(boostVec, FourMomentum, RESULT.FourMomentum);
    RESULT.ComputeBetaGamma();
    return RESULT;
}

double CatsLorentzVector::GetR() const{
    return Length;
}
double CatsLorentzVector::GetR2() const{
    return Length2;
}
double CatsLorentzVector::GetP() const{
    return TotMom;
}
double CatsLorentzVector::GetP2() const{
    return TotMom2;
}
double CatsLorentzVector::GetPt() const{
    return sqrt(FourMomentum[1]*FourMomentum[1]+FourMomentum[2]*FourMomentum[2]);
}
double CatsLorentzVector::Mag() const{
    return Magnitude;
}
double CatsLorentzVector::Mag2() const{
    return Magnitude2;
}
double CatsLorentzVector::GetPseudoRap() const{
    return 0.5*log((TotMom+FourMomentum[3])/(TotMom-FourMomentum[3]));
}
double CatsLorentzVector::GetRapidity() const{
    return 0.5*log((FourMomentum[0]+FourMomentum[3])/(FourMomentum[0]-FourMomentum[3]));
}
double CatsLorentzVector::GetScatAngle() const{
    return acos(GetCosScatAngle());
}
double CatsLorentzVector::GetCosScatAngle() const{
    return (FourSpace[1]*FourMomentum[1]+FourSpace[2]*FourMomentum[2]+FourSpace[3]*FourMomentum[3])/(GetR()*GetP());
}

CatsParticle::CatsParticle(){
    Pid=0;
    Mass=0;
}
CatsParticle::~CatsParticle(){

}
void CatsParticle::ReadFromOscarFile(FILE *InFile){
    int ParticleNr;
    if(
        !fscanf(InFile,"%i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",
        &ParticleNr,&Pid,
        &FourMomentum[1],&FourMomentum[2],&FourMomentum[3],&FourMomentum[0],
        &Mass,
        &FourSpace[1],&FourSpace[2],&FourSpace[3],&FourSpace[0])){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading the OscarFile!\n");
    }
    ComputeBetaGamma();
}
void CatsParticle::SetPid(const int& pid){
    Pid=pid;
}
void CatsParticle::SetMass(const int& mass){
    Mass=mass;
}
int CatsParticle::GetPid() const{
    return Pid;
}
double CatsParticle::GetMass() const{
    return Mass;
}

void CatsParticle::operator=(const CatsParticle& other){
    CatsLorentzVector::operator = (other);
    Pid = other.Pid;
    Mass = other.Mass;
}

CatsParticlePair::CatsParticlePair(){

}
CatsParticlePair::~CatsParticlePair(){

}

void CatsParticlePair::SetPair(const CatsParticle& particle1, const CatsParticle& particle2, const bool& TauCorrection, const bool& BOOST){
    Particle1=particle1;
    Particle2=particle2;
    ParticleSum=Particle1+Particle2;
    CatsLorentzVector::operator = (Particle1-Particle2);

    if(BOOST){
        Particle1.Boost(ParticleSum);
        Particle2.Boost(ParticleSum);
        Boost(ParticleSum);
    }

    if(TauCorrection){
        CatsParticle* FirstParticle;
        double deltaT;
        if(Particle1.GetT()<Particle2.GetT()){
            deltaT = -FourSpace[0];
            FirstParticle = &Particle1;
        }
        else{
            deltaT = FourSpace[0];
            FirstParticle = &Particle2;
        }

        FirstParticle->FourSpace[0] += deltaT;
        FirstParticle->FourSpace[1] += FirstParticle->betaX*deltaT;
        FirstParticle->FourSpace[2] += FirstParticle->betaY*deltaT;
        FirstParticle->FourSpace[3] += FirstParticle->betaZ*deltaT;

        ParticleSum=Particle1+Particle2;
        CatsLorentzVector::operator = (Particle1-Particle2);

    }
}
const CatsParticle& CatsParticlePair::GetParticle(const int& WhichParticle) const{
    switch(WhichParticle){
        case 0 : return Particle1;
        case 1 : return Particle2;
        default : printf("ERROR in CatsParticlePair::GetParticle: TRYING TO ACCESS A NON-EXISTENT ADDRESS!\n"); return Particle1;
    }
}
const CatsLorentzVector& CatsParticlePair::GetSum() const{
    return ParticleSum;
}

CatsEvent::CatsEvent(const int& pid1, const int& pid2):Pid1(pid1),Pid2(pid2){
    BufferSize = 64;
    ParticleType1 = new CatsParticle [BufferSize];
    ParticleType2 = new CatsParticle [BufferSize];
    ParticlePair = NULL;
    NumParticles1 = 0;
    NumParticles2 = 0;
    NumPairs = 0;
}
CatsEvent::~CatsEvent(){
    delete [] ParticleType1;
    delete [] ParticleType2;
    if(ParticlePair) {delete[]ParticlePair; ParticlePair=NULL;}
}
void CatsEvent::Reset(){
    NumParticles1 = 0;
    NumParticles2 = 0;
    NumPairs = 0;
}
void CatsEvent::AddParticle(const CatsParticle& Particle){
    if(Particle.GetPid()==Pid1){
        if(NumParticles1==BufferSize){
            BufferSize *= 2;
            CatsParticle* Temp = new CatsParticle[BufferSize];
            for(unsigned uPart=0; uPart<NumParticles1; uPart++){
                Temp[uPart] = ParticleType1[uPart];
            }
            delete [] ParticleType1;
            ParticleType1 = Temp;
        }
        ParticleType1[NumParticles1] = Particle;
        NumParticles1++;
    }
    else if(Particle.GetPid()==Pid2){
        if(NumParticles2==BufferSize){
            BufferSize *= 2;
            CatsParticle* Temp = new CatsParticle[BufferSize];
            for(unsigned uPart=0; uPart<NumParticles2; uPart++){
                Temp[uPart] = ParticleType2[uPart];
            }
            delete [] ParticleType2;
            ParticleType2 = Temp;
        }
        ParticleType2[NumParticles2] = Particle;
        NumParticles2++;
    }
    else{
        return;
    }
}
void CatsEvent::ComputeParticlePairs(const bool& TauCorrection, const bool& BOOST){
    if(Pid1==Pid2){
        if(NumParticles2){
            printf("WARNING in CatsEvent::ComputeParticlePairs: Potential bug detected! Please contact the developers!\n");
        }
        NumPairs = ( (NumParticles1-1)*NumParticles1 )/2;
    }
    else{
        NumPairs = NumParticles1*NumParticles2;
    }
    if(ParticlePair) {delete [] ParticlePair; ParticlePair=NULL;}
    if(!NumPairs) return;

    ParticlePair = new CatsParticlePair [NumPairs];

    CatsParticle* PartType1 = ParticleType1;
    CatsParticle* PartType2 = Pid1==Pid2?ParticleType1:ParticleType2;
    unsigned& NumPart1 = NumParticles1;
    unsigned& NumPart2 = Pid1==Pid2?NumParticles1:NumParticles2;
    unsigned uPair=0;

    for(unsigned uPart1=0; uPart1<NumPart1; uPart1++){
        for(unsigned uPart2=(Pid1==Pid2?uPart1+1:0); uPart2<NumPart2; uPart2++){
            ParticlePair[uPair].SetPair(PartType1[uPart1],PartType2[uPart2],TauCorrection,BOOST);
            uPair++;
        }
    }

    if(uPair!=NumPairs){
        printf("WARNING in CatsEvent::ComputeParticlePairs: Please contact the developers regarding uPair!=NumPairs\n");
    }

}
unsigned CatsEvent::GetNumPairs() const{
    return NumPairs;
}
CatsParticlePair& CatsEvent::GetParticlePair(const unsigned& WhichPair) const{
    if(WhichPair>=NumPairs){
        printf("WARNING in CatsEvent::GetParticlePair: Pointing to a non-existing pair!\n");
    }
    return ParticlePair[WhichPair];
}
unsigned CatsEvent::GetNumParticles1() const{
    return NumParticles1;
}
unsigned CatsEvent::GetNumParticles2() const{
    return Pid1==Pid2?NumParticles1:NumParticles2;
}
const CatsParticle& CatsEvent::GetParticleType1(const unsigned& WhichPart) const{
    if(WhichPart>=NumParticles1){
        printf("WARNING in CatsEvent::GetParticleType1: Pointing to a non-existing pair!\n");
    }
    return ParticleType1[WhichPart];
}
const CatsParticle& CatsEvent::GetParticleType2(const unsigned& WhichPart) const{
    if(WhichPart>=GetNumParticles2()){
        printf("WARNING in CatsEvent::GetParticleType2: Pointing to a non-existing pair!\n");
    }
    return Pid1==Pid2?ParticleType1[WhichPart]:ParticleType2[WhichPart];
}
bool CatsEvent::GetSameType() const{
    return Pid1==Pid2;
}

CatsDataBuffer::CatsDataBuffer(const unsigned& bsize, const int& pid1, const int& pid2):NumEvents(bsize){
    NumSePairs=0;
    NumMePairs=0;
    TotalNumPairs=0;
    DataEvent = new const CatsEvent* [NumEvents];
    for(unsigned uEve=0; uEve<NumEvents; uEve++){
        DataEvent[uEve] = NULL;
    }
    MixedParticlePair=NULL;
    PointerToPair=NULL;
}
CatsDataBuffer::~CatsDataBuffer(){
    delete [] DataEvent;
    DataEvent = NULL;

    if(MixedParticlePair){
        delete [] MixedParticlePair;
        MixedParticlePair = NULL;
    }

    if(PointerToPair){
        delete [] PointerToPair;
        PointerToPair = NULL;
    }
}
void CatsDataBuffer::SetEvent(const unsigned& WhichEvent, const CatsEvent& Event){
    if(WhichEvent>=NumEvents) return;
    DataEvent[WhichEvent] = &Event;
}
unsigned CatsDataBuffer::GetNumPairsSameEvent() const{
    return NumSePairs;
}
unsigned CatsDataBuffer::GetNumPairsMixedEvent() const{
    return NumMePairs;
}
unsigned CatsDataBuffer::GetNumPairs() const{
    return TotalNumPairs;
}
double CatsDataBuffer::GetAvgNumPairs() const{
    double Average = 0;
    for(unsigned uEve=0; uEve<NumEvents; uEve++){
        Average += DataEvent[uEve]->GetNumPairs();
    }
    return Average/double(NumEvents);
}
const CatsParticlePair* CatsDataBuffer::GetPair(const unsigned& WhichPair) const{
    if(WhichPair>=TotalNumPairs) return NULL;
    return PointerToPair[WhichPair];
}
const CatsParticlePair* CatsDataBuffer::GetSePair(const unsigned& WhichPair) const{
    if(WhichPair>=NumSePairs) return NULL;
    return PointerToPair[WhichPair];
}
const CatsParticlePair* CatsDataBuffer::GetMePair(const unsigned& WhichPair) const{
    if(WhichPair>=NumMePairs) return NULL;
    return PointerToPair[NumSePairs+WhichPair];
}
void CatsDataBuffer::GoBabyGo(const bool& TauCorrection, const bool& BOOST){
    NumSePairs=0;
    NumMePairs=0;
    TotalNumPairs=0;

    for(unsigned uEve=0; uEve<NumEvents; uEve++){
        if(DataEvent[uEve]==NULL) continue;
        NumSePairs+=DataEvent[uEve]->GetNumPairs();
    }

    for(unsigned uEve1=0; uEve1<NumEvents; uEve1++){
        if(DataEvent[uEve1]==NULL) continue;
        for(unsigned uEve2=uEve1+1; uEve2<NumEvents; uEve2++){
            if(DataEvent[uEve2]==NULL) continue;
            NumMePairs += DataEvent[uEve1]->GetNumParticles1()*DataEvent[uEve2]->GetNumParticles2();
            //if we have just one type of particle, than ParticleType1 and ParticleType2 are the same
            //and we will double-count unless we continue at this point!
            if(DataEvent[uEve1]->GetSameType()) continue;
            NumMePairs += DataEvent[uEve1]->GetNumParticles2()*DataEvent[uEve2]->GetNumParticles1();
        }
    }
    if(PointerToPair) {delete[]PointerToPair; PointerToPair=NULL;}
    PointerToPair = new const CatsParticlePair* [NumSePairs+NumMePairs];

    unsigned uSePair=0;
    for(unsigned uEve=0; uEve<NumEvents; uEve++){
        if(DataEvent[uEve]==NULL) continue;
        for(unsigned uPair=0; uPair<DataEvent[uEve]->GetNumPairs(); uPair++){
            PointerToPair[TotalNumPairs] = &DataEvent[uEve]->GetParticlePair(uPair);
            uSePair++;
            TotalNumPairs++;
        }
    }
    if(uSePair!=NumSePairs){
        printf("WARNING in CatsDataBuffer::GoBabyGo: Please contact the developers regarding uSePair!=NumSePairs\n");
    }
    if(MixedParticlePair) delete [] MixedParticlePair;
    MixedParticlePair = new CatsParticlePair [NumMePairs];
    unsigned uMePair=0;
    for(unsigned uEve1=0; uEve1<NumEvents; uEve1++){
        if(DataEvent[uEve1]==NULL) continue;
        for(unsigned uEve2=uEve1+1; uEve2<NumEvents; uEve2++){
            if(DataEvent[uEve2]==NULL) continue;
            for(unsigned uPair1=0; uPair1<DataEvent[uEve1]->GetNumParticles1(); uPair1++){
                for(unsigned uPair2=0; uPair2<DataEvent[uEve2]->GetNumParticles2(); uPair2++){
                    MixedParticlePair[uMePair].SetPair(
                            DataEvent[uEve1]->GetParticleType1(uPair1),
                            DataEvent[uEve2]->GetParticleType2(uPair2),
                                                       TauCorrection,BOOST);
                    PointerToPair[TotalNumPairs] = &MixedParticlePair[uMePair];
                    uMePair++;
                    TotalNumPairs++;
                }
            }
            //if we have just one type of particle, than ParticleType1 and ParticleType2 are the same
            //and we will double-count unless we continue at this point!
            if(DataEvent[uEve1]->GetSameType()) continue;
            for(unsigned uPair1=0; uPair1<DataEvent[uEve2]->GetNumParticles1(); uPair1++){
                for(unsigned uPair2=0; uPair2<DataEvent[uEve1]->GetNumParticles2(); uPair2++){
                    MixedParticlePair[uMePair].SetPair(
                            DataEvent[uEve2]->GetParticleType1(uPair1),
                            DataEvent[uEve1]->GetParticleType2(uPair2),
                                                       TauCorrection,BOOST);
                    PointerToPair[TotalNumPairs] = &MixedParticlePair[uMePair];
                    uMePair++;
                    TotalNumPairs++;
                }
            }
        }
    }
    if(uMePair!=NumMePairs){
        printf("WARNING in CatsDataBuffer::GoBabyGo: Please contact the developers regarding uMePair!=NumMePairs\n");
    }
}

CATSnode::CATSnode(CATSelder* elder, const short& depth, const unsigned& firstid, const unsigned& lastid, double* mean, double* len,
                   const CATSnode* TemplateNode):
                   Elder(elder),Depth(depth),FirstID(firstid),LastID(lastid){
    MeanVal = NULL;
    IntLen = NULL;
    if(elder!=this) StandardNodeInit(mean, len, TemplateNode);
}

CATSnode::~CATSnode(){
    if(MeanVal) {delete [] MeanVal; MeanVal=NULL;}
    if(IntLen) {delete [] IntLen; IntLen=NULL;}

    if(child){
        for(unsigned uSub=0; uSub<Elder->NumSubNodes; uSub++){
            delete child[uSub];
        }
        delete [] child;
        child = NULL;
    }
}

unsigned CATSnode::GetNumOfBoxes(){
    return LastID-FirstID+1;
}
unsigned CATSnode::GetNumOfEl(){
    if(Elder->SourceContext) return 0;
    if(Elder==this) return Elder->NumOfEl;
    return round(SourceValue*double(Elder->NumOfEl));
}
double CATSnode::GetGridSize(){
    return GridSize;
}

void CATSnode::Update(){
    Update(false);
}

void CATSnode::Update(const bool& ThisNodeOnly){
    if(Elder->SourceContext){
        for(short sDim=0; sDim<Elder->Dim; sDim++){
            Elder->SourcePars[1+sDim] = MeanVal[sDim];
        }
        SourceValue = Elder->SourceFunction(Elder->SourceContext)*GridSize;
    }
    else if(Elder->GridBoxId){
        unsigned first, last;
        first = Elder->FindFirstParticleWithID(FirstID);
        last = Elder->FindLastParticleWithID(LastID);

        //the normal situation
        if(first<Elder->NumOfEl && last<Elder->NumOfEl){
            SourceValue = double(last-first+1)/double(Elder->NumOfEl);
        }
        //the case where we should include all particles
        //(both boxes are outside the range, first on the low side and last on the upper side)
        else if(first==Elder->NumOfEl && last==first+1){
            SourceValue = 1;
        }
        //the case where the first box is outside range on the low side
        else if(first==Elder->NumOfEl && last<Elder->NumOfEl){
            SourceValue = double(last+1)/double(Elder->NumOfEl);
        }
        //the case where the last box is outside range on the up side
        else if(first<Elder->NumOfEl && last>Elder->NumOfEl){
            SourceValue = double(Elder->NumOfEl-first)/double(Elder->NumOfEl);
        }
        else{
            SourceValue = 0;
        }
    }
    else{
       SourceValue = 0;
    }
    if(child && !ThisNodeOnly){
        for(unsigned uSub=0; uSub<Elder->NumSubNodes; uSub++) child[uSub]->Update();
    }
}

void CATSnode::StandardNodeInit(double* mean, double* len, const CATSnode* TemplateNode){
    SourceValue = 0;
    MeanVal = new double [Elder->Dim];
    IntLen = new double [Elder->Dim];
    GridSize=1;
    for(short sDim=0; sDim<Elder->Dim; sDim++){
        MeanVal[sDim] = mean[sDim];
        IntLen[sDim] = len[sDim];
        GridSize *= IntLen[sDim];
    }

    Update(true);

    child = NULL;
    bool HasChildren;

    //if we use a template node, we just look if that one has children
    if(TemplateNode){
        HasChildren = TemplateNode->child;
    }
    //this is the standard condition
    else{
        HasChildren = (Depth<Elder->MinDepth || (Depth<Elder->MaxDepth && SourceValue>Elder->Epsilon &&
                                                 GetNumOfEl()>=Elder->MinEntries));
    }

    if( HasChildren ){
        child = new CATSnode* [Elder->NumSubNodes];
        //we want to divide our total interval in two for each parameter on the grid.
        //in order to keep track in which "quadrant" we are, we introduce a very simple counter WhichPart for each
        //of the parameters, that can only take values 0 or 1. Each time WhichPart[x] is increased to 2, than it is set to zero
        //and WhichPart[x+1] is increased, i.e. we continue to iterate over the next parameter.
        char WhichPart[Elder->Dim];
        double ChildMean[Elder->Dim];
        double ChildLen[Elder->Dim];
        unsigned ChildNumBoxes = (LastID-FirstID+1)/Elder->NumSubNodes;
        unsigned ChildFirstID = FirstID;
        for(short sDim=0; sDim<Elder->Dim; sDim++){
            WhichPart[sDim] = 0;
            ChildMean[sDim] = MeanVal[sDim]-IntLen[sDim]*0.25;
            ChildLen[sDim] = IntLen[sDim]*0.5;
        }
        for(unsigned uSub=0; uSub<Elder->NumSubNodes; uSub++){
            child[uSub] = new CATSnode(Elder,Depth+1,ChildFirstID,ChildFirstID+ChildNumBoxes-1,ChildMean,ChildLen,
                                       TemplateNode?TemplateNode->child[uSub]:NULL);
            ChildFirstID += ChildNumBoxes;
            for(short sDim=0; sDim<Elder->Dim; sDim++){
                WhichPart[sDim] = (WhichPart[sDim]+1)%2;
                ChildMean[sDim] = MeanVal[sDim]-IntLen[sDim]*0.25+0.5*IntLen[sDim]*WhichPart[sDim];
                if(WhichPart[sDim]) break;
            }
        }
    }
    else{
        Elder->AddEndNode(this);
    }
}

//! see what happens if epsilon==0
CATSelder::CATSelder(const short& dim, const short& mindep, const short& maxdep, const double& epsilon, double* mean, double* len,
                     void* context, double* Pars, int64_t* gbid, const unsigned& numel):
    CATSnode(this, 0, 0, uipow(2,dim*maxdep)-1, mean, len),
    Dim(dim),MinDepth(mindep),MaxDepth(maxdep),Epsilon(epsilon),NumSubNodes(uipow(2,Dim)),NumOfEl(numel){
    BaseConstructor(mean, len, context, Pars, gbid, numel, NULL);
}

CATSelder::CATSelder(const CATSelder* TemplateElder, void* context, double* Pars, int64_t* gbid, const unsigned& numel):
    CATSnode(this, 0, 0, uipow(2,TemplateElder->Dim*TemplateElder->MaxDepth)-1, TemplateElder->MeanVal, TemplateElder->IntLen),
    Dim(TemplateElder->Dim),MinDepth(TemplateElder->MinDepth),MaxDepth(TemplateElder->MaxDepth),
    Epsilon(TemplateElder->Epsilon),NumSubNodes(uipow(2,Dim)),NumOfEl(numel){

    BaseConstructor(TemplateElder->MeanVal, TemplateElder->IntLen, context, Pars, gbid, numel, TemplateElder);

}

void CATSelder::BaseConstructor(double* mean, double* len, void* context, double* Pars, int64_t* gbid, const unsigned& numel,
                                const CATSelder* TemplateElder){

    if(TemplateElder){
        MaxNumEndNodes = TemplateElder->NumEndNodes;
        NumEndNodes=0;
    }
    else{
        MaxNumEndNodes = unsigned(1./Epsilon)*NumSubNodes;
        NumEndNodes=0;
    }

    SourceRenormError = 0;
    EndNode = new CATSnode* [MaxNumEndNodes];

    SourceContext = context;
    SourcePars = Pars;
    GridBoxId = gbid;

    //the min number of entries required in a node
    //I believe this should be zero in case we work with an Ana Source
    if(Elder->SourceContext) MinEntries=0;
    else MinEntries = Elder->NumSubNodes*4;

    if( (!SourceContext && !SourcePars && !GridBoxId) ||
        (SourceContext && SourcePars && GridBoxId) ||
        (!SourceContext && SourcePars) || (SourceContext && !SourcePars) ||
        (!GridBoxId && NumOfEl) || (GridBoxId && !NumOfEl)
       ){
        printf("WARNING! CATSelder say that the input to the constructor makes no sense!\n");
        SourceContext = NULL;
        SourcePars = NULL;
        GridBoxId = NULL;
    }
    StandardNodeInit(mean, len, TemplateElder);
    if(TemplateElder && NumEndNodes!=TemplateElder->NumEndNodes){
        printf("WARNING! A potential huge bug in CATSelder::BaseConstructor, please contact the developers!\n");
    }

    double Integral=0;
    for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
        Integral += EndNode[uNode]->SourceValue;
    }

    SourceRenormError = Integral?1./Integral:1e64;
    for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
        EndNode[uNode]->SourceValue *= SourceRenormError;
    }

    if(SourceRenormError<1) SourceRenormError = SourceRenormError?1./SourceRenormError:1e64;
    SourceRenormError -= 1;

    double EffectiveEpsilon = double(MinEntries)/double(NumOfEl);
    if(Epsilon>EffectiveEpsilon) EffectiveEpsilon=Epsilon;

    //if(SourceRenormError>Epsilon*sqrt(double(NumEndNodes))){
    if(SourceRenormError>EffectiveEpsilon*double(NumEndNodes)){
        printf("WARNING: CATSelder says that SourceRenormError=%.6f, which seems odd.\n", SourceRenormError);
        printf("         Either the source function is not properly normalized or there is a bug in CATS!\n");
        printf("         In case its the letter, please contact the developers!\n");
    }
}

CATSelder::~CATSelder(){
    delete [] EndNode;
    EndNode = NULL;
}

short CATSelder::GetMaxDepth(){
    return MaxDepth;
}

unsigned CATSelder::GetNumEndNodes(){
    return NumEndNodes;
}

void CATSelder::GetParValues(const unsigned& WhichNode, double* values){
    if(WhichNode>=NumEndNodes) return;
    for(short sDim=0; sDim<Dim; sDim++){
        values[sDim] = EndNode[WhichNode]->MeanVal[sDim];
    }
}
double CATSelder::GetParValue(const unsigned& WhichNode, const short& WhichPar){
    if(WhichNode>=NumEndNodes || WhichPar<0 || WhichPar>=Dim) return 0;
    return EndNode[WhichNode]->MeanVal[WhichPar];
}

double CATSelder::GetGridValue(const unsigned& WhichNode, const bool& Normalized){
    if(WhichNode>=NumEndNodes) return 0;
    return EndNode[WhichNode]->SourceValue/(Normalized?EndNode[WhichNode]->GridSize:1);
}

double CATSelder::GetGridError(const unsigned& WhichNode, const bool& Normalized){
    if(WhichNode>=NumEndNodes) return 0;
    if(SourceContext) return SourceRenormError/(Normalized?EndNode[WhichNode]->GridSize:1);
    else return pow(double(EndNode[WhichNode]->GetNumOfEl()),-0.5)/(Normalized?EndNode[WhichNode]->GridSize:1);
}
void CATSelder::GetGridAxis(const unsigned& WhichNode, double* Axis){
    for(short sDim=0; sDim<Dim; sDim++){
        Axis[sDim]=0;
    }
    if(WhichNode>=NumEndNodes) return;
    for(short sDim=0; sDim<Dim; sDim++){
        Axis[sDim]=EndNode[WhichNode]->MeanVal[sDim];
    }
}

//btw, if the range is outside the limits, the return value will be equal
//to the NumberOfBoxes. Used somewhere else this might lead to potential segmentation faults, so
//make sure to take care of that!
unsigned CATSelder::GetBoxId(double* particle){
    double ParentMean[Dim];
    double ParentLen[Dim];
    for(short sDim=0; sDim<Dim; sDim++){
        ParentMean[sDim] = MeanVal[sDim];
        ParentLen[sDim] = IntLen[sDim];
    }
    short ChildDepth=0;
    unsigned ChildLastID = LastID;
    unsigned ChildFirstID = FirstID;
    double ChildMean[Dim];
    double ChildLen[Dim];
    unsigned ChildNumBoxes;

    //we want to divide our total interval in two for each parameter on the grid.
    //in order to keep track in which "quadrant" we are, we introduce a very simple counter WhichPart for each
    //of the parameters, that can only take values 0 or 1. Each time WhichPart[x] is increased to 2, than it is set to zero
    //and WhichPart[x+1] is increased, i.e. we continue to iterate over the next parameter.
    char WhichPart[Dim];
    bool ThisBox;

    while( ChildDepth<MaxDepth ){
        ChildNumBoxes = (ChildLastID-ChildFirstID+1)/NumSubNodes;
        //ChildFirstID = ParentFirstID;
        //this is the first node
        for(short sDim=0; sDim<Dim; sDim++){
            WhichPart[sDim] = 0;
            ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25;
            ChildLen[sDim] = ParentLen[sDim]*0.5;

        }
        ChildDepth++;
        for(unsigned uSub=0; uSub<NumSubNodes; uSub++){
            ChildLastID = ChildFirstID+ChildNumBoxes-1;
            ThisBox = true;
            //see if the particles in located in one of the following boxes
            for(short sDim=0; sDim<Dim; sDim++){
                ThisBox *= (particle[sDim]>=ChildMean[sDim]-0.5*ChildLen[sDim] &&
                            particle[sDim]<=ChildMean[sDim]+0.5*ChildLen[sDim]);
            }
            //when we find the correct bin, we break out of the loop
            if(ThisBox){
                break;
            }
            ChildFirstID += ChildNumBoxes;
            for(short sDim=0; sDim<Dim; sDim++){
                WhichPart[sDim] = (WhichPart[sDim]+1)%2;
                ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25+0.5*ParentLen[sDim]*WhichPart[sDim];
                if(WhichPart[sDim]) break;
            }
        }
        //update the values for the parent
        for(short sDim=0; sDim<Dim; sDim++){
            ParentMean[sDim] = ChildMean[sDim];
            ParentLen[sDim] = ChildLen[sDim];
        }
    }
    return ChildFirstID;
}

unsigned CATSelder::FindFirstParticleWithID(const unsigned& gbid){
    if(!GridBoxId) return 0;
    if(NumOfEl<=1) return 0;
    if(GridBoxId[0]==gbid) return 0;

    unsigned Position=NumOfEl/2;
    unsigned Mod=4;
    unsigned Step;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(gbid<GridBoxId[0]) return NumOfEl;
    if(gbid>GridBoxId[NumOfEl-1]) return NumOfEl+1;

    while(true){
        if(GridBoxId[Position]<gbid && GridBoxId[Position+1]>=gbid){
            return Position+1;
        }
        else if(gbid<=GridBoxId[Position]){
            Step = NumOfEl/Mod;
            Position -= Step?Step:1;
            Mod *= 2;
        }
        else{
            Step = NumOfEl/Mod;
            Position += Step?Step:1;
            Mod *= 2;
        }
    }
}

unsigned CATSelder::FindLastParticleWithID(const unsigned& gbid){
    if(!GridBoxId) return 0;
    if(NumOfEl<=1) return 0;
    if(GridBoxId[NumOfEl-1]==gbid) return NumOfEl-1;

    unsigned Position=NumOfEl/2;
    unsigned Mod=4;
    unsigned Step;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(gbid<GridBoxId[0]) return NumOfEl;
    if(gbid>GridBoxId[NumOfEl-1]) return NumOfEl+1;

    while(true){
        if(GridBoxId[Position]>gbid && GridBoxId[Position-1]<=gbid){
            return Position-1;
        }
        else if(gbid<GridBoxId[Position]){
            Step = NumOfEl/Mod;
            Position -= Step?Step:1;
            Mod *= 2;
        }
        else{
            Step = NumOfEl/Mod;
            Position += Step?Step:1;
            Mod *= 2;
        }
    }
}

void CATSelder::AddEndNode(CATSnode* node){
    if(!node) return;
    if(MaxNumEndNodes==NumEndNodes){
        MaxNumEndNodes *= 2;
        CATSnode** TempNode = new CATSnode* [MaxNumEndNodes];
        for(unsigned uNode=0; uNode<NumEndNodes; uNode++){
            TempNode[uNode] = EndNode[uNode];
        }
        delete [] EndNode;
        EndNode = TempNode;
    }
    EndNode[NumEndNodes] = node;
    NumEndNodes++;
}

double CATSelder::SourceFunction(void* context){
    return static_cast<CATS*>(context)->EvaluateTheSource(SourcePars);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double CoulombEta(const double& Momentum, const double& RedMass, const double& Q1Q2){
    if(!Momentum) return 0;
    //return 0.5*AlphaFS*RedMass*Q1Q2/Momentum;//I think that one is wrong, verify again
    return AlphaFS*RedMass*Q1Q2/Momentum;
}

//h function, as defined in Lednicky 1981 paper (Yad.Fiz. 35 (1981) 1316-1330)
//the x^2 is replaced with 1/eta^2
double CoulombEuler(const double& eta){
    if(!eta) return 0;
    double RESULT = 0;
    double ADD;
    const double eta2 = eta*eta;
    for(double dIter=1; dIter<=11; dIter++){
        ADD = 1./(dIter*(dIter*dIter+eta2));
        RESULT += ADD;
        if(fabs(ADD/RESULT)<1e-7) break;
    }
    RESULT *= eta2;
    RESULT -= log(eta2)+EulerConst;
    return RESULT;
}

//Momentum = k (the Ac function)
double CoulombPenetrationFactor(const double& eta){
    //if Q1Q2 is zero, than we have no correction.
    return eta?2.*Pi*eta/(exp(2.*Pi*eta)-1.):1;
}

complex<double> GamowCorrection(const double& Momentum, const double& RedMass, const double& Q1Q2){
    //if Q1Q2 is zero, than we have no correction.
    double eta = CoulombEta(Momentum,RedMass,Q1Q2);
    gsl_sf_result lnr;
    gsl_sf_result arg;
    //this functions returns in lnr the ln(abs(Gamma)), while arg is the arg(Gamma)
    //this is so, as numerically ln(Gamma) is evaluated, and ln(Gamma) = lnr+i*arg
    gsl_sf_lngamma_complex_e(1.,eta,&lnr,&arg);
    return eta?exp(i*arg.val)*sqrt(CoulombPenetrationFactor(eta)):1;
}
