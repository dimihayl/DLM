//! Product:         Correlation Analysis Tools using the Schrödinger equation (CATS)
//! Current Version: 2.2 (16 October 2017)
//! Copyright:       Dimitar Lubomirov Mihaylov (Technical University of Munich)
//! Support:         dimitar.mihaylov(at)mytum.de
//! Documentation:   a full documentation is not available yet

#ifndef CATS_H
#define CATS_H

#include <stdint.h>

#include "DLM_MergeSort.h"
#include "CATStools.h"

//this typedef it used to save the potentials for the
//different channels as an array of function pointers.
typedef double (*CatsPotential)(double*);

//assumptions: the potential is radial-symmetric.
//internally only Gaussian natural units (in MeV !!!) are used,
//i.e. c=ħ=kB=4πε0=1

//!the input is assumed to be in [MeV] for M,P,E
//!for the radii: [fm]

class CATS;
class CATScontainerOLD;
class CATSelder;
class CATSnode;

class CATS{
//friend class CATSelder;
friend class CATSnode;
public:
    CATS();
    ~CATS();

    //!Sets and gets
    void SetRedMass(const double& redMass);
    double GetRedMass();

    void SetPdgId(const int& id1, const int& id2);
    void GetPdgId(int& id1, int& id2);

    //If the number of polarizations is changed, all previous input about the
    //polarization themselves is lost (i.e. NumPW is reset!)
    void SetNumChannels(const unsigned short& numCh);
    unsigned short GetNumChannels();

    //!N.B. here usCh plays the role of the spin quantum number
    void SetNumPW(const unsigned short& usCh, const unsigned short& numPW);
    unsigned short GetNumPW(const unsigned short& usCh);

    void SetQ1Q2(const int& q1q2);
    int GetQ1Q2();

    unsigned GetNumMomBins();
    unsigned GetNumIpBins();
    unsigned GetNumPairs();

    //N.B. the size of mombins should be NumMomBins+1, where each element represents the low-edge of
    //the corresponding bin, and the one extra element is the upper edge of the last bin.
    void SetMomBins(const unsigned& nummombins, const double* mombins);
    void SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom);

    void SetIpBins(const unsigned& numBbins, const double* imppar);
    void SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar);

    void SetChannelWeight(const unsigned short& usCh, const double& weight);
    double GetChannelWeight(const unsigned short& usCh);

    void SetStartRad(const double& srad);
    double GetStartRad();

    void SetEpsilonProp(const double& epsp);
    double GetEpsilonProp();

    void SetEpsilonConv(const double& epsc);
    double GetEpsilonConv();

    //in fm
    void SetMaxRad(const double& maxrad);
    double GetMaxRad();

    void SetMaxRho(const double& maxrho);
    double GetMaxRho();

    void SetExcludeFailedBins(const bool& efb);
    bool GetExcludeFailedBins();

    void SetMaxGridDepth(const short& mgd);
    short GetMaxGridDepth();
    void SetSourceMinValOnGrid(const double& smvg);
    double GetSourceMinValOnGrid();
    void SetMaxNumGridPart(const unsigned& mngp);
    unsigned GetMaxNumGridPart();

    void SetGridMinDepth(const short& val);
    short GetGridMinDepth();
    void SetGridMaxDepth(const short& val);
    short GetGridManDepth();
    void SetGridEpsilon(const double& val);
    double GetGridEpsilon();

    void SetUseAnalyticSource(const bool& val);
    bool GetUseAnalyticSource();

    void SetThetaDependentSource(const bool& val);
    bool GetThetaDependentSource();

    void SetTransportRenorm(const double& val);
    double GetTransportRenorm();

    //the input values should be non-negative
    //please set this condition BEFORE you load the source, else CATS will not save the TotalMomentum at all
    //The values should be in MeV
    //if noReload is true, than CATS will not reload the data from the file. The user should
    void SetTotPairMomCut(const double& minval, const double& maxval);
    void GetTotPairMomCut(double& minval, double& maxval);
    void RemoveTotPairMomCut();

    void SetNotifications(const short& notify);

    void SetMaxPairsPerBin(unsigned mpp);
    unsigned GetMaxPairsPerBin();

    void SetMaxPairsToRead(unsigned mpp);
    unsigned GetMaxPairsToRead();

    //void SetMaxPairsToLoad(unsigned mpp);
    //unsigned GetMaxPairsToLoad();

    //void SetEventMixing(const bool& mix);
    //bool GetEventMixing();
    void SetMixingDepth(const unsigned short& mix);
    unsigned short GetMixingDepth();

    void SetBufferEventMix(const unsigned& bem);
    unsigned GetBufferEventMix();

    void SetTauCorrection(const bool& tc);
    bool GetTauCorrection();

    void SetInputFileName(const char* fname);
    void GetInputFileName(char* fname);

    unsigned GetNumPairsPerBin(const unsigned& uMomBin, const unsigned& uIpBin);
    unsigned GetNumPairsPerBin(const unsigned& uMomBin);

    void GetPairInfo(const unsigned& uWhichPair,
                     double& RelMom, double& RelPos, double& RelCosTh, double& TotMom);
    void GetPairInfo(const unsigned& uWhichPair, double* Output);

    unsigned GetLoadedPairs(const unsigned& WhichMomBin, const unsigned& WhichIpBin);
    unsigned GetRelativeMomentum(const unsigned& WhichParticle);
    unsigned GetRelativePosition(const unsigned& WhichParticle);
    unsigned GetRelativeCosTheta(const unsigned& WhichParticle);
    unsigned GetTotalPairMomentum(const unsigned& WhichParticle);

//    void SetLogFileName(const char* fname);
//    void GetLogFileName(char* fname);

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFun(const unsigned& WhichMomBin);
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFun(const unsigned& WhichMomBin, double& Momentum);
    //evaluates Ck at this point based on interpolation
    double EvalCorrFun(const double& Momentum);

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFunErr(const unsigned& WhichMomBin);
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFunErr(const unsigned& WhichMomBin, double& Momentum);
    //evaluates Ck at this point based on interpolation
    double EvalCorrFunErr(const double& Momentum);

    //The same, but for a specific impact parameter bin
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin);
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin, double& Momentum, double& ImpPar);

    double GetPhaseShift(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW);
    double EvalPhaseShift(const double& Momentum, const unsigned short& usCh, const unsigned short& usPW);

    unsigned GetNumRadialWFpts(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW);
    double GetRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const unsigned& WhichRadBin);
    double EvalRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const double& Radius,
                                  const bool& DevideByR=true);

    //The momentum in the WhichMomBin-th bin
    double GetMomentum(const unsigned& WhichMomBin);
    double GetMomBinLowEdge(const unsigned& WhichMomBin);
    double GetMomBinUpEdge(const unsigned& WhichMomBin);

    unsigned GetMomBin(const double& Momentum);
    unsigned GetIpBin(const double& bVal);

    //convert fm to 1/MeV
    const double& FmNu();
    //convert 1/MeV to fm
    const double& NuFm();

    //[0]-[2] reserved for CATS ([0] is Momentum, [1] is Radius (fm), [2] is CosTheta)
    //n.b. for the time being CATS assumes radial symmetric potential, thus [2] is actually never used,
    //i.e. please always use only radial symmetric potential
    double*** PotPar;

    void RemoveShortRangePotential();
    void RemoveShortRangePotential(const unsigned& usCh, const unsigned& usPW);
    //Pars[0] should be the radius, Pars[1] should be the momentum.
    //!Pars should be an array with at least 2 elements
    void SetShortRangePotential(const unsigned& usCh, const unsigned& usPW,
                           double (*pot)(double* Pars), double* Pars);


    void RemoveAnaSource();
    //input vars: [0] should always be the momentum (MeV), [1] the radius (fm) and [2] 'cosθ'
    double* AnaSourcePar;
    void SetAnaSource(double (*AS)(double*), double* Pars);

    //!------------------------------------------------

    //!Running the analysis
    //CATS tries to be clever and runs many of the function only if it is needed and they were not computed before.
    //However it keeps no track if the input file, the parameters of the analytic source or of the potential
    //have changed. By default it assumes that they have not.
    //!Thus run CATS with KillTheCat(KillOptions) if you have changed anything!
    void KillTheCat(const int& Options=kNothingChanged);
    void ComputeTheRadialWaveFunction();
    //!------------------------------------------------

    enum KillOptions { kNothingChanged, kSourceChanged, kPotentialChanged, kAllChanged };
    enum NotificationOptions { nSilent, nError, nWarning, nAll };

private:

    enum PrevNext { kNext=1, kPrevious=-1 };

    //!Variables needed as an input for the physics calculation
    double RedMass;
    int pdgID[2];

    //Number of channels
    unsigned short NumCh;
    //Number of partial waves for each polarization
    unsigned short* NumPW;
    bool IdenticalParticles;
    //!CATS will only work if there is either a short range potential (falls off faster than 1/r in the asymptotic region)
    //!or in case there is a Coulomb-like potential (Q1Q2/r). This parameter declares which case is valid for the current calculation.
    //charge x charge of the two particles, for the Coulomb potential it should be != 0
    int Q1Q2;
    //!------------------------------------------------

    //!Variables needed as settings

    //Maximum number of pairs to be analyzed per bin (momentum<->ImpPar)
    unsigned MaxPairsPerBin;
    //the maximum pairs to be read from the input file
    unsigned MaxPairsToRead;
    //the maximum # pairs to load from the file
    //unsigned MaxPairsToLoad;

    //if true, CATS will treat events with the same impact parameter as a single event.
    //This is sometimes a necessity in order to gain more statistics, however the implementation of event-mixing
    //assumes that 'b' and 'k' are not strongly correlated in the range if interest. In case of a small correlation the user
    //is advised to make further investigation of possible systematical errors.
    //By default this option is switched off!
    //bool EventMixing;

    //in case of EventMixing, this variable specifies what is the max. number of particles pairs to mix.
    //increasing this number will improve the error, but the computational cost goes up.
    unsigned short MixingDepth;
    bool TauCorrection;
    bool UseAnalyticSource;
    bool ThetaDependentSource;
    //multiply the source of the transport model by some coefficient
    double TransportRenorm;

    double MinTotPairMom;
    double MaxTotPairMom;
    double LoadedMinTotPairMom;
    double LoadedMaxTotPairMom;
    bool UseTotMomCut;
    //0 => no notifications of any kind
    //1 => show only errors
    //2 => show errors and warnings
    //3 => full output
    short Notifications;
    //!------------------------------------------------

    //!Variables needed as an input for the numerical calculation
    //N.B. the size of mombins should be NumMomBins+1, where each element represents the low-edge of
    //the corresponding bin, and the one extra element is the upper edge of the last bin.
    unsigned NumMomBins;
    double* MomBin;
    double* ChannelWeight;

    unsigned NumIpBins;
    double* IpBin;

    //the very first radius to be computed. Related with the shape of the potential. As a rule of thumb,
    //this value should be at least an order of magnitude smaller than the smallest desirable resolution of the numerical method.
    //by default this is set to 0.005 fm. This parameter relates to the initial step size as well as the minimal step-size
    //that the solver will be allowed to use.
    double StartRad;
    //practically scales linearly with the step size. A perfect value would be around 1e-7, which would balance perfectly
    //between numerical and machine precision. It is recommended to keep this value between 1e-4 and 1e-9. The default value is 5e-6
    double EpsilonProp;
    //determines the criteria for a convergence. It is the threshold value for the relative difference between
    //the propagating function with or without potential. Similarly as for EpsilonProp it is assumed that the perfect value
    //should be around 1e-7. By default CATS uses EpsilonConv = 5e-6
    double EpsilonConv;

    //break the numerical computation if the algorithm has failed to converge to the asymptotic solution up to
    //a particular Radius value. By default MaxRad==32 fm
    //N.B. this parameter has an influence on the data that is loaded for the source. I.e. changing it will influence not only
    //the Schr. solver but also the source.
    double MaxRad;

    //the same, but this time a condition for rho. Both conditions are useful and needed, since at high momenta the
    //convergence region is much more determined by the radius, but at low momenta the rho coefficient may be much
    //more important in the presence of a short-range potential. The default value is 16
    //N.B. this parameter influences ONLY the Schroedinger solver (unlike MaxRad)
    double MaxRho;

    bool ExcludeFailedConvergence;

    //5 = default value
    short GridMinDepth;
    //0 = default value (14 for 1D grid, 10 for 2D grid)
    short GridMaxDepth;
    //0 = default value (1/1024 for 1D grid, 1/8192 for 2D grid), max value is 0.125
    double GridEpsilon;

    char* InputFileName;
//    char* LogFileName;
    //total number of selected pairs
    unsigned NumPairs;
    //percentage of selected same event pairs with a specific impact parameter
    double* WeightIp;
    //N.B. at the moment the error is overestimated, due to the assumption that all WeightIps are independent.
    //this should be a fairly small effect though, thus it is probably okay to leave the code as it is.
    double* WeightIpError;

    //in bins of momentum
    unsigned* LoadedPairsPerMomBin;
    //in bins of momentum/ImpactParameter
    unsigned** LoadedPairsPerBin;
    //all particle pairs
    double* RelativeMomentum;
    double* RelativePosition;
    double* RelativeCosTheta;
    double* TotalPairMomentum;
    unsigned* PairMomBin;
    unsigned* PairIpBin;
    int64_t* GridBoxId;

    //for all particle pairs
    CATSelder* BaseSourceGrid;
    //in bins of momentum only
    CATSelder** kSourceGrid;
    double** WaveFunction2;
    //in bins of momentum/ImpactParameter
    CATSelder*** kbSourceGrid;

    //CATScontainerOLD*** ParticleContainer1D;
    //CATScontainerOLD*** ParticleContainer2D;

    bool LoadedData;//i.e. the data-file was read
    bool SourceGridReady;//i.e. the Particle container is set up properly
    bool ComputedWaveFunction;
    bool ComputedCorrFunction;

    //!INFO ABOUT THE ABOVE 3 VARIABLES
    //one should be mindful that at large relative momenta (k above 200 MeV) the solution converges at higher rho values.
    //this means that, especially for a Coulomb potential, that one can be in a situation where the result does not converge
    //within the limits set. Thus the above 3 parameters should be carefully adjusted.

    //after the wave-function is computed the result is saved (as a function of rho) in a equidistant grid
    //double FinalRhoGrid;
    //!------------------------------------------------
    //!Stuff needed as an input for the numerical calculation

    //!THE INPUT FOR THE POTENTIAL IS ASSUMED TO BE IN [fm]
    //!THE OUTPUT SHOULD BE IN [MeV]
    CatsPotential** ShortRangePotential;

    double CoulombPotential(const double& Radius);

    //input vars: [0] should always be the momentum, [1] the radius and [2] 'cosθ'
    double (*AnalyticSource)(double*);
    //!------------------------------------------------

    //!Any other variables or functions used at runtime internally by CATS

    double CurrentRhoStep;
    //!------------------------------------------------

    //!Constants
    const double Pi;
    //fine-structure constant (==e^2 in Gaussian units)
    const double AlphaFS;
    const double RevSqrt2;
    //convert fm into natural units (1/MeV)
    const double FmToNu;
    const double NuToFm;
    //!------------------------------------------------

    //!Functions used internally by CATS
    //the differential equation for the Schroedinger equation
    void PropagatingFunction(double& Basic, double& Full,
                               const double& Radius, const double& Momentum,
                               const unsigned short& AzQN, const unsigned short& Pol);

    void ComputeWaveFunction();
    void ComputeTotWaveFunction(const bool& ReallocateTotWaveFun);
    short LoadData(const unsigned short& NumBlankHeaderLines=3);
    unsigned LoadDataBuffer(const unsigned& WhichIpBin, CatsDataBuffer* KittyBuffer);
    void FoldSourceAndWF();
    void SortAllData();
    void SetUpSourceGrid();

    float ProgressCompute;
    float ProgressLoad;
    float ProgressFold;

    //plane partial wave as a solution from the gsl libraries
    double PlanePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);
    //coulomb partial wave as a solution from the gsl libraries
    double CoulombPartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);
    //radial/coulomb partial wave as a solution from the gsl libraries
    double ReferencePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW);

    double AsymptoticRatio(const double& Radius, const double& Momentum, const unsigned short& usPW);

    //a numerical root-finder. Very fast and accurate for well-behaved (near to linear) functions
    double NewtonRapson(double (CATS::*Function)(const double&, const double&, const unsigned short&),
                        const double& EpsilonX, const unsigned short& usPW, const double& Momentum,
                          const double&  xMin, const double&  xMax, const double& fValShift);

    //void ResortData(double* input, DLM_MergeSort <int64_t, unsigned>& Sorter);
    //void ResortData(unsigned* input, DLM_MergeSort <int64_t, unsigned>& Sorter);
    template <class Type> void ResortData(Type* input, DLM_MergeSort <int64_t, unsigned>& Sorter){
        unsigned NumOfEl = Sorter.GetNumOfEl();
        Type* Temp;
        Temp = new Type[NumOfEl];
        for(unsigned uEl=0; uEl<NumOfEl; uEl++){
            Temp[uEl] = input[Sorter.GetKey()[uEl]];
        }
        for(unsigned uEl=0; uEl<NumOfEl; uEl++){
            input[uEl] = Temp[uEl];
        }
        delete [] Temp;
    }
    unsigned GetBoxId(double* particle);

    //evaluates the solution to the radial equation based on the numerical result and the computed phaseshift.
    //I.e. if Radius is within the computed range, we extrapolate based on the result. If Radius is outside
    //the computed range we use the shifted reference wave. If DivideByR==true, computed is R = u/r.
    //N.B. The result would differ from EvalWaveFunctionU/Radius due to the extrapolation done.
    double EvalWaveFunctionU(const unsigned& uMomBin, const double& Radius,
                             const unsigned short& usCh, const unsigned short& usPW, const bool& DivideByR);
    double EffectiveFunction(const unsigned& uMomBin, const double& Radius, const unsigned short& usCh);
    double EffectiveFunction(const unsigned& uMomBin, const double& Radius);

    double EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta, const unsigned short& usCh);
    double EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta);

    //computes the momentum bin corresponding to Momentum
    unsigned GetBin(const double& Value, const double* Range, const unsigned& NumBins);
    //computes the radius bin corresponding to Radius (for a certain l and S)
    unsigned GetRadBin(const double& Radius, const unsigned& uMomBin,
                       const unsigned short& usCh, const unsigned short& usPW);


    //double SourceFunction(const double* Pars, const double& GridSize);

    //delete all variables that depend only on the number of momentum bins, b-bins and MaxPairs
    void DelMomIpMp();
    //delete all variables that depend only on the number of b-bins
    void DelIp();
    //delete all variables that depend only on the number of channels
    void DelCh();
    //delete all variables that depend on the number of momentum bins, number of channels and number of partial waves
    void DelMomChPw();
    //delete all variables that depend only on the number of momentum bins
    void DelMom();
    //delete all variables that depend only on the number of momentum and b-bins
    void DelMomIp();

    //delete all variables that depend on the number of momentum bins
    void DelAllMom();
    //delete all variables that depend on the number of b-bins
    void DelAllIp();
    //delete all variables that depend on the number of channels
    void DelAllCh();

    void DelAll();

    //!------------------------------------------------

    //!Variables used to save the output
    //limits the amount of memory used to save information about the wave function.
    //unsigned MaxWaveFunBins;
    //the WaveFunRad show the value of r at which the WaveFunction is evaluated
    unsigned*** SavedWaveFunBins;//in bins of mom/pol/pw
    //double*** RadStepWF;//in bins of mom/pol/pw
    double*** PhaseShift;//in bins of mom/pol/pw, saved only until the end of each k-iteration
    double**** WaveFunRad;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    double**** WaveFunctionU;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    bool* MomBinConverged;//bins of mom, marked as true in case the num. comp. failed and this bin should not be used

    //in bins of momentum/ImpactParameter
    double** kbCorrFun;
    double** kbCorrFunErr;

    //in bins of momentum
    double* kCorrFun;
    double* kCorrFunErr;


};



#endif // CATS_H
