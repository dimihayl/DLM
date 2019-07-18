
#ifndef CATSTOOLS_H
#define CATSTOOLS_H

#include <stdio.h>
#include <stdint.h>
#include <complex>

class CatsParticle;

class CATSparameters{
friend class CATSelder;
public:
    enum Type {tSource=3,tPotential=2};
    CATSparameters(const unsigned type, const unsigned numpar, const bool threadsafe);
    CATSparameters(const CATSparameters& other);
    ~CATSparameters();
    double* GetParameters() const;
    //the parameter is set the same for all threads!
    void SetParameter(const unsigned& WhichPar, const double& Value, const bool& CurrentThread=false);
    void SetParameters(const double* pars, const bool& CurrentThread=false);
    //the variable is set only for the current thread!
    void SetVariable(const unsigned& WhichVar, const double& Value, const bool& CurrentThread);
    double GetParameter(const unsigned& WhichPar) const;
    double GetVariable(const unsigned& WhichVar) const;
    unsigned GetNumPars() const;
    unsigned GetTotNumPars() const;
    bool operator ==(const CATSparameters &other) const;
protected:
    //the number of dummy parameters (variables)
    const unsigned NumVars;
    //the number of actual parameters
    const unsigned NumPars;
    //dummy+actual parameters
    const unsigned TotNumPars;
    const bool ThreadSafe;
    const unsigned NumThreads;
    double** Parameter;
};

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
    double GetPt() const;

    double Mag() const;
    double Mag2() const;

    double GetPseudoRap() const;
    double GetRapidity() const;
    //returns the angle between the spacial and momentum vector
    double GetScatAngle() const;
    double GetCosScatAngle() const;

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
    //rotates the Momentum vector in Phi
    void RotateMomPhi(const double& angle);
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
    void operator=(const CatsLorentzVector& other);
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

    void SetPair(const CatsParticle& particle1, const CatsParticle& particle2, const bool& TauCorrection=false, const bool& BOOST=true);
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
    void ComputeParticlePairs(const bool& TauCorrection=false, const bool& BOOST=true);
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
    double GetAvgNumPairs() const;
    const CatsParticlePair* GetPair(const unsigned& WhichPair) const;
    const CatsParticlePair* GetSePair(const unsigned& WhichPair) const;
    const CatsParticlePair* GetMePair(const unsigned& WhichPair) const;
    void GoBabyGo(const bool& TauCorrection=false, const bool& BOOST=true);
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
    double GetGridSize();
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
              double* mean, double* len, void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel);
    CATSelder(const CATSelder* TemplateElder,
              void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel);
    void BaseConstructor(double* mean, double* len, void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel,
                         const CATSelder* TemplateElder);
    ~CATSelder();

    short GetMaxDepth();
    unsigned GetNumEndNodes();
    void GetParValues(const unsigned& WhichNode, double* values);
    double GetParValue(const unsigned& WhichNode, const short& WhichPar);
    double GetGridValue(const unsigned& WhichNode, const bool& Normalized=false);
    double GetGridError(const unsigned& WhichNode, const bool& Normalized=false);
    void GetGridAxis(const unsigned& WhichNode, double* Axis);

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
    void* SourceContext;
    CATSparameters* SourcePars;
    int64_t* GridBoxId;
    const unsigned NumOfEl;

    double SourceFunction(void* context);

    //CATS* Kitty;
};

double CatsSourceForwarder(void* context, double* Pars);
class CatsSource{
public:
    virtual ~CatsSource();
    virtual double Eval(double* Pars);
    virtual void SetParameter(const unsigned& WhichPar, const double& Value);
    virtual unsigned GetNumPars();
    double Eval(const double& Momentum, const double Radius, const double& Angle);
private:
    double PARS[3];
};

double CoulombEta(const double& Momentum, const double& RedMass, const double& Q1Q2);
double CoulombEuler(const double& eta);
double CoulombPenetrationFactor(const double& eta);
std::complex<double> GamowCorrection(const double& Momentum, const double& RedMass, const double& Q1Q2);

//pLab to pCm, Mass2 is the mass of the particle at rest
double pLab_pCm(const double& pLab, const double& Mass1, const double& Mass2);
//tLab to kCm, Mass2 is the mass of the particle at rest
double tLab_pCm(const double& tLab, const double& Mass1, const double& Mass2);
double pCm_pLab(const double& pCm, const double& Mass1, const double& Mass2);
double pCm_tLab(const double& pCm, const double& Mass1, const double& Mass2);

#endif // CATSTOOLS_H
