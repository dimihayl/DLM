
#ifndef CATSTOOLS_H
#define CATSTOOLS_H

#include <stdio.h>
#include <stdint.h>
#include <complex>
#include <vector>

class CatsParticle;
class DLM_Random;

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
    void BoostBack(const CatsLorentzVector& boostVec);
    //CatsLorentzVector GetBoost(const CatsLorentzVector& boostVec);

    double GetR() const;
    double GetR2() const;

    double GetP() const;
    double GetP2() const;
    double GetPt() const;
    double GetMt() const;

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
    double GetR(const int& xyz) const;
    double GetPhi() const;
    double GetTheta() const;
    double GetE() const;
    double GetP(const int& xyz) const;
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    double GetPtheta() const;
    double GetPphi() const;
    double Gamma() const;
    double Beta() const;
    double Beta(const int& xyz) const;
    double BetaX() const;
    double BetaY() const;
    double BetaZ() const;


    void Set(const double& tCrd, const double& xCrd, const double& yCrd, const double& zCrd,
             const double& engy, const double& xMom, const double& yMom, const double& zMom);
    //sets the momentum components, keeping the MASS constant (reevaluates energy)
    void SetMomXYZ(const double& xMom, const double& yMom, const double& zMom);
    void SetMXYZ(const double& mass, const double& xMom, const double& yMom, const double& zMom);
    void SetMPtEtaPhi(const double& mass, const double& pt, const double& eta, const double& phi);
    void SetTXYZ(const double& tCrd, const double& xCrd, const double& yCrd, const double& zCrd);

    //propagates the particle for a time tau, moving along a straight line
    //by default we consider that the time given is the proper time, meaning that we
    //have gamma*tau for the propagation time
    //if the proper time is false, we simply use tau
    void Propagate(const double& tau, const bool& proper_time=true);

    //rotates the Momentum vector in Phi
    void RotateMomPhi(const double& angle);
    void RenormSpacialCoordinates(const double& Renorm);
    void Print();
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
    void BoostBack(const CatsLorentzVector& boostVec, const double* InVec, double* OutVec);
    void ComputeBetaGamma();
};

class CatsParticle:public CatsLorentzVector{
public:
    CatsParticle();
    ~CatsParticle();
    void ReadFromOscarFile(FILE *InFile);
    void SetPid(const int& pid);
    void SetMass(const double& mass);
    void SetZeroMassEpsilon(const double& eps);
    void SetWidth(const double& width);
    int GetPid() const;
    double GetMass() const;
    double GetWidth() const;
    //CatsParticle*& Decay(const unsigned& Nbody, const double* mass);
    //propagate==true means that the coordinates of the mother will be shifted
    //N.B. this is always done for the daughters
    //to the point at which the decay occurs. The lifetime is assumed 1/Width. If width is zero, no propagation is done
    //N.B. if also only works if we have a SetDecayRanGen activated!!!
    //N.B. for the propagation to work, we need the correct unit conversion between space and momentum
    //by default we assume that we work with MeV and fm, i.e. we have an hbarc conversion factor of c.a. 197 MeV*fm
    //you can change this conversion by Set_hbarc, i.e. setting it to 1 assumes we work in natural units
    CatsParticle* Decay(const double& mass1, const double& mass2, const bool& propagate=true);
    //for the Nbody decay. NOT DONE YET, WORKS ONLY WITH 2-BODY
    CatsParticle* Decay(const std::vector<double> masses, const bool& propagate=true);
    CatsParticle* DecaySimple(const std::vector<double> masses, const bool& propagate=true);
    void SetDecayRanGen(DLM_Random* rangen);
    void SetDecayRanGen(DLM_Random& rangen);
    void Set_hbarc(const double& HBARC);
    void operator=(const CatsParticle& other);
    void operator=(const CatsLorentzVector& other);
protected:
    int Pid;
    //the mass will be concidered zero if below some limit
    //(for numerical stability)
    double ZeroMass;
    double Width;
    DLM_Random* RanGen;
    double hbarc;
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
  //what might be possible to make more efficient, make pointers here
  //at first glance there is no need to copy the particles, but make sure this is so!
  //on second thought... at some point there is a boost performed, and perhaps it is not wise
  //to boost the original particles. Thus the copy is needed.
    CatsParticle Particle1;
    CatsParticle Particle2;
    CatsLorentzVector ParticleSum;
};

class CatsEvent{
public:
    CatsEvent(const int& pid1, const int& pid2);
    ~CatsEvent();
    void Reset();
    //the particle is copied!
    void AddParticle(const CatsParticle& Particle);
    void ComputeParticlePairs(const bool& TauCorrection=false, const bool& BOOST=true);
    unsigned GetNumPairs() const;
    CatsParticlePair& GetParticlePair(const unsigned& WhichPair) const;
    unsigned GetNumParticles1() const;
    unsigned GetNumParticles2() const;
    int GetPidParticle1() const;
    int GetPidParticle2() const;
    const CatsParticle& GetParticleType1(const unsigned& WhichPart) const;
    const CatsParticle& GetParticleType2(const unsigned& WhichPart) const;
    bool GetSameType() const;
    void SetRandomSeed(const unsigned& SEED);
    //randomizes all particles in the event, according to their phi distribution
    //distribution = "Gauss", "Uniform"
    //option =
    //      "" : each particle smeared independently
    //      "EnergyConservation" : each second particle is rotated such, that E,p conservation is restored
    void RandomizeMomentumPhi(const double& smearValue, const char* option = "", const char* distribution = "Uniform");
private:
    CatsParticle* ParticleType1;
    CatsParticle* ParticleType2;
    CatsParticlePair* ParticlePair;
    DLM_Random* RanGen;

    const int Pid1;
    const int Pid2;

    unsigned NumParticles1;
    unsigned NumParticles2;
    unsigned NumPairs;

    unsigned BufferSize1;
    unsigned BufferSize2;
};

//an object used to study a specific N-body problem, e.g. 2,3,4 whatever body problem,
//and this object corresponds to one set of N-bodies.
class CatsMultiplet{
public:
  enum RefFrame { LAB, CM };
  //copy_particles==true means that each time a particle is set, it will be copied to this class
  //else you will point to the original object, and it is up to the user to keep track
  //of the memory. N.B. the latter can be a big deal in memory saving and reduced CPU time due to copy,
  //which may be useful in large scale simulations.
  CatsMultiplet(DLM_Random& ran_gen, const unsigned& nbody, const bool& copy_particles);
  //void SetParticle(const unsigned& which_one, CatsParticle& particle, );
  ~CatsMultiplet();

  //KEEP THE TAU CORRECTION AS AN OPTION!!!
  void Compute(const bool& TauCorrection=true);
  CatsParticle* GetParticle(const int& ref_frame, const unsigned& which_one) const;

//SO HERE: TRY TO SOMEHOW PUT THE GETQ (4-mom) DEFINITION + THE OPTION OF KSTAR (3-vec) FOR TWO_BODY CASE
//FOR THE GETR, LEAVE ROOM OPEN FOR OTHER OPTIONS LATER ON, FOR NOW DO THE STANDARD 2-body 3-vec
//AND FOR MULTIPLETS... WELL PERHAPS JUST FOR FUN DO THE 4-VECTOR STUFF, WHY THE HELL NOT
  //the relative distance between two particles
  double GetRstar(const int& ref_frame, const unsigned& part1, const unsigned& part2) const;
  //based on 4-vectors
  double GetRinv(const int& ref_frame, const unsigned& part1, const unsigned& part2) const;
  //same as above but for all particles (R_N similar to Q_N)
  double GetRinv(const int& ref_frame);
  //relative momentum between two particles defined based on 3-vectors
  double GetQstar(const int& ref_frame, const unsigned& part1, const unsigned& part2) const;
  //1/2 of Qstar
  double GetKstar(const int& ref_frame, const unsigned& part1, const unsigned& part2) const;
  //the relative momentum between two particles, based on 4-momenta
  //N.B. for two body this is approx 2*kstar !!!
  double GetQinv(const int& ref_frame, const unsigned& part1, const unsigned& part2) const;
  //for QN definition: https://www.annualreviews.org/doi/pdf/10.1146/annurev.nucl.55.090704.151533
  //aslo info here: https://arxiv.org/pdf/1502.02121.pdf
  double GetQinv(const int& ref_frame) const;
private:
  DLM_Random& RanGen;
  const unsigned Nbody;
  const bool CopyParticles;
  //these particles can be pointers or objects
  CatsParticle** Particle;
  //these guys are always owned by the CatsMultiplet
  CatsParticle** ParticleCm;

};

//at the moment I do not check if the events loaded are of the same PID type (which should be the case!)
//either be careful with that or add some check about it!
class CatsDataBuffer{
public:
    CatsDataBuffer(const unsigned& bsize, const int& pid1, const int& pid2);
    ~CatsDataBuffer();
    //the event is pointed to
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
              double* mean, double* len, void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel, const bool& renorm=true);
    CATSelder(const CATSelder* TemplateElder,
              void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel, const bool& renorm=true);
    void BaseConstructor(double* mean, double* len, void* context, CATSparameters* Pars, int64_t* gbid, const unsigned& numel,
                         const CATSelder* TemplateElder, const bool& renorm=true);
    ~CATSelder();

    short GetMaxDepth();
    unsigned GetNumEndNodes();
    void GetParValues(const unsigned& WhichNode, double* values);
    double GetParValue(const unsigned& WhichNode, const short& WhichPar);
    double GetGridValue(const unsigned& WhichNode, const bool& Normalized=false);
    double GetGridError(const unsigned& WhichNode, const bool& Normalized=false);
    void GetGridAxis(const unsigned& WhichNode, double* Axis);
    void Renormalize();

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
