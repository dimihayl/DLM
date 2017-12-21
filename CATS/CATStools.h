
#ifndef CATSTOOLS_H
#define CATSTOOLS_H

#include <stdio.h>
#include <stdint.h>

class CatsParticle;

class CatsLorentzVector{
friend class CatsParticlePair;
public:
    CatsLorentzVector();
    ~CatsLorentzVector();

    void Boost(const CatsLorentzVector& boostVec);
    CatsLorentzVector GetBoost(const CatsLorentzVector& boostVec);

    double GetR() const;
    double GetR2() const;

    double GetP() const;
    double GetP2() const;

    double Mag() const;
    double Mag2() const;

    double GetT() const;
    double GetX() const;
    double GetY() const;
    double GetZ() const;
    double GetE() const;
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    void Set(const double& tCrd, const double& xCrd, const double& yCrd, const double& zCrd,
             const double& engy, const double& xMom, const double& yMom, const double& zMom);
    void RenormSpacialCoordinates(const double& Renorm);
    CatsLorentzVector const operator+(const CatsLorentzVector& other);
    CatsLorentzVector const operator-(const CatsLorentzVector& other);
    void operator=(const CatsLorentzVector& other);
protected:
    double FourSpace[4];
    double FourMomentum[4];

    double Length;
    double Length2;
    double TotMom;
    double TotMom2;
    double Magnitude;
    double Magnitude2;
    double gamma;
    double betaX;
    double betaY;
    double betaZ;
    double beta;

    void Boost(const CatsLorentzVector& boostVec, const double* InVec, double* OutVec);
    void ComputeBetaGamma();
};

class CatsParticle:public CatsLorentzVector{
public:
    CatsParticle();
    ~CatsParticle();
    void ReadFromOscarFile(FILE *InFile);
    void SetPid(const int& pid);
    void SetMass(const int& mass);
    int GetPid() const;
    double GetMass() const;
    void operator=(const CatsParticle& other);
protected:
    int Pid;
    double Mass;
};

//contains all info about the particles in their CM system.
//the object itself has the coordinates of the difference of the two particles.
//ParticleSum is NOT transformed in CM, but is rather given in LAB!
class CatsParticlePair:public CatsLorentzVector{
public:
    CatsParticlePair();
    ~CatsParticlePair();

    void SetPair(const CatsParticle& particle1, const CatsParticle& particle2, const bool& TauCorrection=false);
    const CatsParticle& GetParticle(const int& WhichParticle) const;
    const CatsLorentzVector& GetSum() const;
protected:
    CatsParticle Particle1;
    CatsParticle Particle2;
    CatsLorentzVector ParticleSum;
};

class CatsEvent{
public:
    CatsEvent(const int& pid1, const int& pid2);
    ~CatsEvent();
    void Reset();
    void AddParticle(const CatsParticle& Particle);
    void ComputeParticlePairs(const bool& TauCorrection=false);
    unsigned GetNumPairs() const;
    CatsParticlePair& GetParticlePair(const unsigned& WhichPair) const;
    unsigned GetNumParticles1() const;
    unsigned GetNumParticles2() const;
    const CatsParticle& GetParticleType1(const unsigned& WhichPart) const;
    const CatsParticle& GetParticleType2(const unsigned& WhichPart) const;
    bool GetSameType() const;
private:
    CatsParticle* ParticleType1;
    CatsParticle* ParticleType2;
    CatsParticlePair* ParticlePair;

    const int Pid1;
    const int Pid2;

    unsigned NumParticles1;
    unsigned NumParticles2;
    unsigned NumPairs;

    unsigned BufferSize;
};

//at the moment I do not check if the events loaded are of the same PID type (which should be the case!)
//either be careful with that or add some check about it!
class CatsDataBuffer{
public:
    CatsDataBuffer(const unsigned& bsize, const int& pid1, const int& pid2);
    ~CatsDataBuffer();
    void SetEvent(const unsigned& WhichEvent, const CatsEvent& Event);
    unsigned GetNumPairsSameEvent() const;
    unsigned GetNumPairsMixedEvent() const;
    unsigned GetNumPairs() const;
    const CatsParticlePair* GetPair(const unsigned& WhichPair) const;
    const CatsParticlePair* GetSePair(const unsigned& WhichPair) const;
    const CatsParticlePair* GetMePair(const unsigned& WhichPair) const;
    void GoBabyGo(const bool& TauCorrection=false);
private:
    const unsigned NumEvents;
    unsigned NumSePairs;
    unsigned NumMePairs;
    unsigned TotalNumPairs;
    const CatsEvent** DataEvent;
    CatsParticlePair* MixedParticlePair;
    const CatsParticlePair** PointerToPair;
};

class CATSelder;

class CATSnode{
friend class CATSelder;
public:
    CATSnode(CATSelder* elder, const short& depth, const unsigned& firstid, const unsigned& lastid,
             double* mean, double* len, const CATSnode* TemplateNode=NULL);
    ~CATSnode();
    unsigned GetNumOfBoxes();
    unsigned GetNumOfEl();
    void Update();
protected:
    CATSelder* Elder;
    const short Depth;
    //id of the boxes at max depth
    const unsigned FirstID;
    const unsigned LastID;
    double SourceValue;
    double* MeanVal;
    double* IntLen;
    double GridSize;
    CATSnode** child;

    void Update(const bool& ThisNodeOnly);
    void StandardNodeInit(double* mean, double* len, const CATSnode* TemplateNode=NULL);
};



//! THE IDEA FOR TOMORROW:
//бате махни kitty от конструктура и сложи пойнтър към GridBoxId и double (*AnalyticSource)(double*).
//съответно от CATS винаги викай конструктура с един от двата пойнтъра NULL. Този който е зададен ще
//бъде използван от старейшината за да си смята източника. За CATS: запомни, че имаш мултиплисити бинс само
//когато имаш комбиниране на събития! Т.е. недей да създаваш излишен на брой старейшини, а само толкова
//колкото са ти необходими
class CATSelder:public CATSnode{
friend class CATSnode;
public:
    CATSelder(const short& dim, const short& mindep, const short& maxdep, const double& epsilon,
              //double* mean, double* len, double (CATS::*sfun)(const double*, const double&));
              double* mean, double* len, double (*AS)(double*), double* Pars, int64_t* gbid, const unsigned& numel);
    CATSelder(const CATSelder* TemplateElder,
              double (*AS)(double*), double* Pars, int64_t* gbid, const unsigned& numel);
    void BaseConstructor(double* mean, double* len, double (*AS)(double*), double* Pars, int64_t* gbid, const unsigned& numel,
                         const CATSelder* TemplateElder);
    ~CATSelder();

    short GetMaxDepth();
    unsigned GetNumEndNodes();
    void GetParValues(const unsigned& WhichNode, double* values);
    double GetParValue(const unsigned& WhichNode, const short& WhichPar);
    double GetGridValue(const unsigned& WhichNode);
    double GetGridError(const unsigned& WhichNode);

    unsigned GetBoxId(double* particle);
    unsigned FindFirstParticleWithID(const unsigned& gbid);
    unsigned FindLastParticleWithID(const unsigned& gbid);

    void AddEndNode(CATSnode* node);

protected:
    const short Dim;
    const short MinDepth;
    const short MaxDepth;
    const double Epsilon;
    const unsigned NumSubNodes;
    unsigned NumEndNodes;
    unsigned MaxNumEndNodes;
    double SourceRenormError;
    unsigned MinEntries;

    CATSnode** EndNode;

    //pars and grid-size
    double (*SourceFunction)(double*);
    double* SourcePars;
    int64_t* GridBoxId;
    const unsigned NumOfEl;
    //CATS* Kitty;
};

template <class Type> class CATShisto{
public:
    CATShisto(const unsigned& numbin, const Type* bins):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        for(unsigned uBin=0; uBin<=NumBins; uBin++){
            BinRange[uBin] = bins[uBin];
        }
        //for(unsigned uBin=0; uBin<NumBins; uBin++){
            //BinValue[uBin] = GetBinCenter(uBin);
        //}
    }
    CATShisto(const unsigned& numbin, const Type& xmin, const Type& xmax):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        if(NumBins==1){
            BinRange[0] = xmin; BinRange[1] = xmax;
            //BinValue[0] = (xmin+xmax)*0.5;
        }
        else{
            Type BinWidth = (xmax-xmin)/Type(NumBins);
//printf("NumBins=%u; BinWidth=%f;\n",NumBins,BinWidth);
            for(unsigned uBin=0; uBin<=NumBins; uBin++){
                BinRange[uBin] = xmin + Type(uBin)*BinWidth;
//printf(" uBin=%u -> BR=%f\n",uBin, BinRange[uBin]);
            }
            //for(unsigned uBin=0; uBin<NumBins; uBin++){
            //    BinValue[uBin] = GetBinCenter(uBin);
            //}
        }

    }
    CATShisto(const CATShisto& other):NumBins(other.NumBins){
        BinRange = NULL;
        BinValue = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        operator=(other);
    }
    ~CATShisto(){
        if(BinRange) {delete [] BinRange; BinRange=NULL;}
        if(BinValue) {delete [] BinValue; BinValue=NULL;}
    }

    void SetBinContent(const unsigned& WhichBin, const double& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]=Val;
    }
    void Add(const unsigned& WhichBin, const double& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]+=Val;
    }
    void SetBinAt(const double& xVal, const double& Val){
        unsigned WhichBin = GetBin(xVal);
        SetBinContent(WhichBin, Val);
    }
    void AddAt(const double& xVal, const double& Val){
        unsigned WhichBin = GetBin(xVal);
        Add(WhichBin, Val);
    }

    unsigned GetBin(const double& xVal) const{
        if(NumBins<=1) return 0;
        unsigned WhichBin=(NumBins+1)/2;
        unsigned BinMod=4;
        unsigned BinStep;
        //makes sure that the xVal is in BinRange. If not, the returned value is either
        //NumBins or NumBins+1, depending on if we have an underflow or overflow
        if(xVal<BinRange[0]) return NumBins;
        if(xVal>BinRange[NumBins]) return NumBins+1;
        while(true){
            if(BinRange[WhichBin]<=xVal && BinRange[WhichBin+1]>=xVal){
                return WhichBin;
            }
    //!check the logic of this <=, I made it so, since else we crashed on limit values, before that it was <
            else if(xVal<=BinRange[WhichBin]){
                BinStep = (NumBins+1)/BinMod;
                WhichBin -= BinStep?BinStep:1;
                BinMod *= 2;
            }
            else{
                BinStep = (NumBins+1)/BinMod;
                WhichBin += BinStep?BinStep:1;
                BinMod *= 2;
            }
        }
    }
    Type GetBinCenter(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return 0.5*(BinRange[WhichBin]+BinRange[WhichBin+1]);
    }
    Type GetBinContent(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return BinValue[WhichBin];
    }
    Type GetBinLowEdge(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return BinRange[WhichBin];
    }
    Type GetBinUpEdge(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return BinRange[WhichBin+1];
    }
    Type Eval(const double& xVal) const{
        if(xVal<BinRange[0] || xVal>BinRange[NumBins]) return 0;
        if(NumBins==1) return BinValue[0];
        unsigned WhichBin = GetBin(xVal);

        Type Value[3];
        Value[0] = WhichBin?GetBinCenter(WhichBin-1):-1;
        Value[1] = GetBinCenter(WhichBin);
        Value[2] = WhichBin<(NumBins-1)?GetBinCenter(WhichBin+1):-1;

        Type* InterpolRange;
        const Type* FunRange;

        if(Value[0]==-1){
            InterpolRange = &Value[1];
            FunRange = &BinValue[WhichBin];
        }
        else if(Value[2]==-1){
            InterpolRange = &Value[0];
            FunRange = &BinValue[WhichBin-1];
        }
        else if(xVal<Value[1]){
            InterpolRange = &Value[0];
            FunRange = &BinValue[WhichBin-1];
        }
        else if(Value[1]<xVal){
            InterpolRange = &Value[1];
            FunRange = &BinValue[WhichBin];
        }
        else{//Value[1]==xVal
            return BinValue[WhichBin];
        }

        if(InterpolRange[1]-InterpolRange[0]){
            return (FunRange[1]*(xVal-InterpolRange[0])-
                    FunRange[0]*(xVal-InterpolRange[1]))/
                    (InterpolRange[1]-InterpolRange[0]);
        }
        else{
            printf("WTF!!! at bin %u\n",WhichBin);
            return 1e6;
        }
    }

    unsigned GetNbins() const{
        return NumBins;
    }

    bool Copy(const CATShisto& other, const Type& Scale){
        if(NumBins!=other.NumBins) return false;
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            BinRange[uBin] = other.BinRange[uBin];
            BinValue[uBin] = other.BinValue[uBin]*Scale;
        }
        BinRange[NumBins] = other.BinRange[NumBins];
/*
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            if(other.BinValue[uBin]>0.95452 && other.BinValue[uBin]<0.95453){
                printf("I found it: the copy is %f\n",BinValue[uBin]);
            }
            printf("%u: OBV=%f; BV=%f; Scale=%f\n", uBin, other.BinValue[uBin], BinValue[uBin], Scale);
        }
*/
        return true;
    }
    bool Add(const CATShisto& other, const Type& Scale){
        if(NumBins!=other.NumBins) return false;
        for(unsigned uBin=0; uBin<=NumBins; uBin++) if(BinRange[uBin]!=other.BinRange[uBin]) return false;
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] += other.BinValue[uBin]*Scale;
        return true;
    }
    void Add(const Type& Value){
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] += Value;
    }

    bool operator=(const CATShisto& other){
        return Copy(other, 1);
    }
    bool operator+=(const CATShisto& other){
        return Add(other, 1);
    }

protected:

    CATShisto(const unsigned& numbin):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
    }

    const unsigned NumBins;
    Type* BinRange;
    Type* BinValue;

};

#endif // CATSTOOLS_H
