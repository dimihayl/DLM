//! Product:         Correlation Analysis Tools using the Schrödinger equation (CATS)
//! Current Version: 2.19 (16th October 2018)
//! Copyright:       Dimitar Lubomirov Mihaylov (Technical University of Munich)
//! Support:         dimitar.mihaylov(at)mytum.de
//! Documentation:   a full documentation is not available yet

//! In case you use CATS for your analysis please cite:
//! D. L. Mihaylov, V. M. Sarti, O. W. Arnold, L. Fabbietti, B. Hohlweger and A. M. Mathis, Eur.Phys.J. C78 (2018) no.5, 394
/*BibTeX:
@article{Mihaylov:2018rva,
      author         = "Mihaylov, D. L. and Mantovani Sarti, V. and Arnold, O. W.
                        and Fabbietti, L. and Hohlweger, B. and Mathis, A. M.",
      title          = "{A femtoscopic Correlation Analysis Tool using the
                        Schr\"odinger equation (CATS)}",
      journal        = "Eur. Phys. J.",
      volume         = "C78",
      year           = "2018",
      number         = "5",
      pages          = "394",
      doi            = "10.1140/epjc/s10052-018-5859-0",
      eprint         = "1802.08481",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      SLACcitation   = "%%CITATION = ARXIV:1802.08481;%%"
}
*/

#ifndef CATS_H
#define CATS_H

#include <stdint.h>
#include <complex>

#include "DLM_Sort.h"
#include "CATStools.h"
//#include "DLM_Histo.h"

using namespace std;

//this typedef it used to save the potentials for the
//different channels as an array of function pointers.
typedef double (*CatsPotential)(double*);

//assumptions: the potential is radial-symmetric.
//internally only Gaussian natural units (in MeV !!!) are used,
//i.e. c=ħ=kB=4πε0=1

//!the input is assumed to be in [MeV] for M,P,E
//!for the radii: [fm]

class CATS;
class CATSelder;
class CATSnode;

class CATS{
friend class CATSelder;
friend class CATSnode;
public:
    CATS();
    ~CATS();

    //!Sets and gets
    void SetRedMass(const double& redMass);
    double GetRedMass() const;

    void SetPdgId(const int& id1, const int& id2);
    void GetPdgId(int& id1, int& id2) const;

    //0 = turned off; 1 = turned on; else = automatically set based on the PdgID
    void SetQuantumStatistics(short qs);
    short GetQuantumStatistics() const;

    //If the number of polarizations is changed, all previous input about the
    //polarization themselves is lost (i.e. NumPW is reset!)
    void SetNumChannels(const unsigned short& numCh);
    unsigned short GetNumChannels() const;

    void SetNumPW(const unsigned short& usCh, const unsigned short& numPW);
    unsigned short GetNumPW(const unsigned short& usCh) const;

    void SetSpin(const unsigned short& usCh, const unsigned short& spin);
    unsigned short GetSpin(const unsigned short& spin) const;

    void SetQ1Q2(const int& q1q2);
    int GetQ1Q2() const;
    void SetGamow(const bool& gamow);
    bool GetGamow() const;

    unsigned GetNumMomBins() const;
    unsigned GetNumIpBins() const;
    unsigned GetNumPairs() const;

    //N.B. the size of mombins should be NumMomBins+1, where each element represents the low-edge of
    //the corresponding bin, and the one extra element is the upper edge of the last bin.
    void SetMomBins(const unsigned& nummombins, const double* mombins, const double* bincenter=NULL);
    void SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom);

    void SetIpBins(const unsigned& numBbins, const double* imppar);
    void SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar);

    void SetChannelWeight(const unsigned short& usCh, const double& weight);
    double GetChannelWeight(const unsigned short& usCh) const;

    void SetStartRad(const double& srad);
    double GetStartRad() const;

    void SetEpsilonProp(const double& epsp);
    double GetEpsilonProp() const;

    void SetEpsilonConv(const double& epsc);
    double GetEpsilonConv() const;

    //in fm
    void SetMaxRad(const double& maxrad);
    double GetMaxRad() const;

    void SetMaxPw(const unsigned short& maxpw);
    unsigned short GetMaxPw() const;

    //if true, the total wave function is computed ONLY by adding the available numerical partial waves
    //this was implemented with the idea to be used whenever the WF for a coupled channel is available
    void SetOnlyNumericalPw(const unsigned short& usCh, const bool& val);
    bool GetOnlyNumericalPw(const unsigned short& usCh) const;

    void SetMaxNumThreads(const unsigned short& maxnumthreads);
    unsigned short GetMaxNumThreads() const;

    void SetMaxRho(const double& maxrho);
    double GetMaxRho() const;

    void SetExcludeFailedBins(const bool& efb);
    bool GetExcludeFailedBins() const;

    void SetMaxGridDepth(const short& mgd);
    short GetMaxGridDepth() const;
    void SetSourceMinValOnGrid(const double& smvg);
    double GetSourceMinValOnGrid() const;
    void SetMaxNumGridPart(const unsigned& mngp);
    unsigned GetMaxNumGridPart() const;

    void SetGridMinDepth(const short& val);
    short GetGridMinDepth() const;
    void SetGridMaxDepth(const short& val);
    short GetGridManDepth() const;
    void SetGridEpsilon(const double& val);
    double GetGridEpsilon() const;

    void SetMomentumDependentSource(const bool& val);

    void SetUseAnalyticSource(const bool& val);
    bool GetUseAnalyticSource() const;

    void SetThetaDependentSource(const bool& val);
    bool GetThetaDependentSource() const;

    void SetTransportRenorm(const double& val);
    double GetTransportRenorm() const;

    void SetPoorManRenorm(const double& val);
    double GetPoorManRenorm() const;

    void SetSourceMinRange(const double& val);
    void SetSourceMaxRange(const double& val);
    double GetSourceMinRange() const;
    double GetSourceMaxRange() const;
    //if true, the source is assumed to be normalized over its full domain.
    //any deviations from unity are attributed to long-range correlations which are assumed
    //to contribute with a flat correlation, which is added automatically
    void SetNormalizedSource(const bool& val=true);
    bool GetNormalizedSource() const;

    //the input values should be non-negative
    //please set this condition BEFORE you load the source, else CATS will not save the TotalMomentum at all
    //The values should be in MeV
    //if noReload is true, than CATS will not reload the data from the file. The user should
    void SetTotPairMomCut(const double& minval, const double& maxval);
    void GetTotPairMomCut(double& minval, double& maxval) const;
    void RemoveTotPairMomCut();

    void SetNotifications(const short& notify);
    short GetNotifications() const;

    void SetMaxPairsPerBin(unsigned mpp);
    unsigned GetMaxPairsPerBin() const;

    void SetMaxPairsToRead(unsigned mpp);
    unsigned GetMaxPairsToRead() const;

    //void SetMaxPairsToLoad(unsigned mpp);
    //unsigned GetMaxPairsToLoad();

    void SetMixingDepth(const unsigned short& mix);
    unsigned short GetMixingDepth() const;

    void SetBufferEventMix(const unsigned& bem);
    unsigned GetBufferEventMix() const;

    void SetTauCorrection(const bool& tc);
    bool GetTauCorrection() const;

    void SetInputFileName(const char* fname);
    void GetInputFileName(char* fname) const;

    unsigned GetNumPairsPerBin(const unsigned& uMomBin, const unsigned& uIpBin) const;
    unsigned GetNumPairsPerBin(const unsigned& uMomBin) const;

    void GetPairInfo(const unsigned& uWhichPair,
                     double& RelMom, double& RelPos, double& RelCosTh, double& TotMom) const;
    void GetPairInfo(const unsigned& uWhichPair, double* Output) const;

    unsigned GetLoadedPairs(const unsigned& WhichMomBin, const unsigned& WhichIpBin) const;
    unsigned GetRelativeMomentum(const unsigned& WhichParticle) const;
    unsigned GetRelativePosition(const unsigned& WhichParticle) const;
    unsigned GetRelativeCosTheta(const unsigned& WhichParticle) const;
    unsigned GetTotalPairMomentum(const unsigned& WhichParticle) const;

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFun(const unsigned& WhichMomBin) const;
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFun(const unsigned& WhichMomBin, double& Momentum) const;
    //evaluates Ck at this point based on interpolation
    double EvalCorrFun(const double& Momentum) const;

    //returns the CorrFun evaluated for some MomBin
    double GetCorrFunErr(const unsigned& WhichMomBin) const;
    //the same, but the value of the corresponding Momentum as saved in Momentum
    double GetCorrFunErr(const unsigned& WhichMomBin, double& Momentum) const;
    //evaluates Ck at this point based on interpolation
    double EvalCorrFunErr(const double& Momentum) const;

    //The same, but for a specific impact parameter bin
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin) const;
    double GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin, double& Momentum, double& ImpPar) const;

    double GetPhaseShift(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW) const;
    float EvalPhaseShift(const double& Momentum, const unsigned short& usCh, const unsigned short& usPW) const;

    unsigned GetNumRadialWFpts(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW) const;
    complex<double> GetRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const unsigned& WhichRadBin) const;
    complex<double> EvalRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const double& Radius,
                                  const bool& DivideByR=true) const;
    double EvalWaveFun2(const unsigned& uMomBin, const double& Radius, const double& CosTheta, const unsigned short& usCh);
    double EvalWaveFun2(const unsigned& uMomBin, const double& Radius, const unsigned short& usCh);

    complex<double> EvalAsymptoticRadialWF(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const double& Radius,
                                  const bool& DivideByR=true);
    complex<double> EvalReferenceRadialWF(const unsigned& WhichMomBin,const unsigned short& usPW, const double& Radius, const bool& DivideByR=true);

    //The momentum in the WhichMomBin-th bin
    double GetMomentum(const unsigned& WhichMomBin) const;
    double GetMomBinLowEdge(const unsigned& WhichMomBin) const;
    double GetMomBinUpEdge(const unsigned& WhichMomBin) const;
    //use with care, the memory allocation is not monitored by cats, i.e. the pointer you
    //get from here should be deleted later in your code!
    double* CopyMomBin();

    unsigned GetMomBin(const double& Momentum) const;
    unsigned GetIpBin(const double& bVal) const;
    //computes the momentum bin corresponding to Momentum
    unsigned GetBin(const double& Value, const double* Range, const unsigned& NumBins) const;
    //computes the radius bin corresponding to Radius (for a certain l and S)
    unsigned GetRadBin(const double& Radius, const unsigned& uMomBin,
                       const unsigned short& usCh, const unsigned short& usPW) const;
    double EvaluateTheSource(CATSparameters* Pars) const;
    double EvaluateTheSource(const double& Momentum, const double& Radius, const double& CosTheta) const;
    unsigned GetNumSourcePars() const;
    double EvaluateThePotential(const unsigned short& usCh, const unsigned short& usPW, const double& Momentum, const double& Radius) const;
    double EvaluateCoulombPotential(const double& Radius) const;
    unsigned GetNumPotPars(const unsigned short& usCh, const unsigned short& usPW) const;
    CATSelder* GetTheElder(const double& Momentum);
    //convert fm to 1/MeV
    const double& FmNu() const;
    //convert 1/MeV to fm
    const double& NuFm() const;


    void RemoveShortRangePotential();
    void RemoveShortRangePotential(const unsigned& usCh, const unsigned& usPW);
    //Pars[0] should be the radius, Pars[1] should be the momentum.
    //!Pars should be an array with at least 2 elements
    void SetShortRangePotential(const unsigned& usCh, const unsigned& usPW, double (*pot)(double* Pars), CATSparameters& Pars);
    //void SetShortRangePotential(const unsigned& usCh, const unsigned& usPW, double (*pot)(double* Pars), double* Pars);
    //set the value of the WhichPar-th parameter of the potential corresponding to the usCh,usPW
    //N.B. WhichPar counts from zero, i.e. CATS sets the value of PotPar[usCh][usPW][3+WhichPar]. Since CATS
    //has no information of the length of this array, it is the responsibility of the user to make source there is
    //no segmentation violation!!!
    void SetShortRangePotential(const unsigned& usCh, const unsigned& usPW, const unsigned& WhichPar, const double& Value);

    void RemoveAnaSource();


    void SetAnaSource(double (*AS)(double*), CATSparameters& Pars);
    //void SetAnaSource(double (*AS)(double*), double* Pars);
    //this definition uses a C++ trick, allowing to pass the member function of any object as an input argument.
    //!This class needs to have a function Eval(double* Pars), where Pars[0,1,2] are the momentum,radius and cosθ
    //Imagine you have a class called "MYCLASS" and it has the Eval function. To use this source from an object MYCLASS OBJ you just need to do:
    //1) define a "forwarder" function double FORWARDER(void* context, double* Pars){return static_cast<MYCLASS*>(context)->Eval(Pars);}
    //2) define a source function as double SOURCE(double (*fptr)(void*, double*), void* context, double* Pars){return fptr(context,Pars);}
    //3) pass to your CATS object the source by calling .SetAnaSource(FORWARDER,OBJ);
    void SetAnaSource(double (*FS)(void*, double*), void* context, const unsigned& numparameters=0);
    //void SetAnaSource(const CatsSource& SOURCE);
/*
///////////////////////////////////////////////////////////////////////////////////////////
double test_ad5_function(double (*fptr)(void*, const double*), void* context, const double* Pars){
    return fptr(context,Pars);
}
double test_ad5_forward(void* context, const double* Pars){
    return static_cast<af_class1*>(context)->afc1_7(Pars);
}

//!so here I can use a member function of an object as an argument
//the way it is done is by introducing a forwarder, which basically tells us where to look for the function of the object
//i.e. we have an object afc1 that is our "context". Then we define a forwarder (test_ad4_forward) that takes the context as an argument
//and calls some member function. The return value of the forwarder should be the same as of the member function.
//To execute all of that we need to define a function that takes as arguments both the forwarder and the context.
void test_ad5(){
    af_class1 afc1;
    double Pars[3] = {1,2,3};
    double result = test_ad5_function(&test_ad5_forward,&afc1,Pars);
    printf("result = %.4f\n",result);
}
///////////////////////////////////////////////////////////////////////////////////////////
*/
    //set the value of the WhichPar-th parameter of the source
    //N.B. WhichPar counts from zero, i.e. CATS sets the value of AnaSourcePar[2+WhichPar].
    //It is the responsibility of the user to respect the size of the array and avoid segmentation fault!
    //if(SmallChange==true) => force CATS to reuse the computing grid
    //This will  gain performance, might decrease the accuracy though, please use only when fine-tuning the source
    void SetAnaSource(const unsigned& WhichPar, const double& Value, const bool& SmallChange=false);
    double GetAnaSourcePar(const unsigned& WhichPar) const;

    //void SetPotPar(const unsigned& WhichPar, const double& Value);
    double GetPotPar(const unsigned& usCh, const unsigned& usPW, const unsigned& WhichPar) const;

    //The input should be in u_l = r*R_l
    void SetExternalWaveFunction(const unsigned& usCh, const unsigned& usPW, DLM_Histo<complex<double>>& histWF, DLM_Histo<complex<double>>& histPS);
    void SetExternalWaveFunction(DLM_Histo<complex<double>>*** ExternalWF);
    void RemoveExternalWaveFunction(const unsigned& usCh, const unsigned& usPW);
    //!------------------------------------------------

    //!Running the analysis
    //CATS tries to be clever and runs many of the function only if it is needed and they were not computed before.
    //However it keeps no track if the input file, the parameters of the analytic source or of the potential
    //have changed. By default it assumes that they have not.
    //!Thus run CATS with KillTheCat(KillOptions) if you have changed anything!
    void KillTheCat(const int& Options=kNothingChanged);
    //true = C(k) is computed; false = C(k) needs to be reevaluated
    bool CkStatus();
    bool PotentialStatus();
    bool SourceStatus();
    void ComputeTheRadialWaveFunction();
    //!------------------------------------------------
    enum KillOptions { kNothingChanged, kSourceChanged, kPotentialChanged, kAllChanged };
    enum NotificationOptions { nSilent, nError, nWarning, nAll };
//DLM_Histo<double> SourceHistoTemp;
protected:

    enum PrevNext { kNext=1, kPrevious=-1 };

    //!Variables needed as an input for the physics calculation
    double RedMass;
    int pdgID[2];
    short QuantumStatistics;

    //Number of channels
    unsigned short NumCh;
    //total spin of each channel
    unsigned short* Spin;
    //Number of partial waves for each polarization
    unsigned short* NumPW;
    bool IdenticalParticles;
    //!CATS will only work if there is either a short range potential (falls off faster than 1/r in the asymptotic region)
    //!or in case there is a Coulomb-like potential (Q1Q2/r). This parameter declares which case is valid for the current calculation.
    //charge x charge of the two particles, for the Coulomb potential it should be != 0
    int Q1Q2;
    //if true it applies a Gamow correction based on Q1Q2
    bool Gamow;
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
    //By default this option is switched off, i.e. MixingDepth==1
    //this variable specifies what is the max. number of events to mix
    unsigned short MixingDepth;
    bool TauCorrection;
    bool UseAnalyticSource;
    bool ThetaDependentSource;
    //multiply the source of the transport model by some coefficient
    double TransportRenorm;
    //an additional renormalization done not when loading the source, but applied simply to the 'r' in CM frame.
    //the advantage of using this option is computational performance when fitting.
    //This condition is only applied to sources from transport models
    double PoorManRenorm;
    //minimum/maximum allowed radius of the source
    //the result is being renormalized such that the integral of the source is unchanged
    double SourceMinRad;
    double SourceMaxRad;
    bool NormalizedSource;

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
    double* MomBinCenter;
    double* ChannelWeight;

    unsigned NumIpBins;
    double* IpBin;

    //this guy should only be modified in ComputeTotWaveFunction, else the WaveFunction2 memory management will fail
    unsigned NumGridPts;

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
    //the max. 'l' to be computed by CATS
    unsigned short MaxPw;
    bool* OnlyNumPw;
    unsigned short MaxNumThreads;

    bool ExcludeFailedConvergence;
    bool MomDepSource;

    //5 = default value
    short GridMinDepth;
    //0 = default value (14 for 1D grid, 10 for 2D grid)
    short GridMaxDepth;
    //0 = default value (1/1024 for 1D grid, 1/8192 for 2D grid), max value is 0.125
    double GridEpsilon;

    char* InputFileName;
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
//double** WaveFunction2;
    //in bins of momentum/ImpactParameter
    CATSelder*** kbSourceGrid;

    bool LoadedData;//i.e. the data-file was read
    bool SourceGridReady;//i.e. the Particle container is set up properly
    bool SourceUpdated;
    bool ComputedWaveFunction;
    bool ComputedCorrFunction;
    bool GamowCorrected;

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

    double CoulombPotential(const double& Radius) const;

    //input vars: [0] should always be the momentum, [1] the radius and [2] 'cosθ'
    double (*AnalyticSource)(double*);
    double (*ForwardedSource)(void*, double*);
    void* SourceContext;
    //CatsSource* MemberSource;


    //[0]-[2] reserved for CATS ([0] is Momentum, [1] is Radius (fm), [2] is CosTheta)
    //n.b. for the time being CATS assumes radial symmetric potential, thus [2] is actually never used,
    //i.e. please always use only radial symmetric potential
    //[usCh][usPW][...]
    CATSparameters*** PotPar;
    //double*** PotParArray;

    //input vars: [0] should always be the momentum (MeV), [1] the radius (fm) and [2] 'cosθ'
    CATSparameters* AnaSourcePar;
    //double* AnaSourceParArray;
    //CATSparameters* ForwardedSourcePar;

    //!------------------------------------------------

    //!Any other variables or functions used at runtime internally by CATS

    double CurrentRhoStep;
    //!------------------------------------------------

    //!Constants
    const unsigned short NumPotPars;
    const unsigned short NumSourcePars;
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
    void UpdateSourceGrid();

    float ProgressCompute;
    float ProgressLoad;
    float ProgressFold;

    //plane partial wave as a solution from the gsl libraries
    double PlanePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW) const;
    //coulomb partial wave as a solution from the gsl libraries
    double CoulombPartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const;
    //radial/coulomb partial wave as a solution from the gsl libraries
    double ReferencePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const;

    double AsymptoticRatio(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const;

    //a numerical root-finder. Very fast and accurate for well-behaved (near to linear) functions
    double NewtonRapson(double (CATS::*Function)(const double&, const double&, const unsigned short&, const int&) const,
                        const double& EpsilonX, const unsigned short& usPW, const double& Momentum, const int& q1q2,
                          const double&  xMin, const double&  xMax, const double& fValShift) const;

    template <class Type> void ResortData(Type* input, DLM_Sort <int64_t, unsigned>& Sorter);
    unsigned GetBoxId(double* particle);

    //evaluates the solution to the radial equation based on the numerical result and the computed phaseshift.
    //I.e. if Radius is within the computed range, we extrapolate based on the result. If Radius is outside
    //the computed range we use the shifted reference wave. If DivideByR==true, computed is R = u/r.
    //N.B. The result would differ from EvalWaveFunctionU/Radius due to the extrapolation done.
    complex<double> EvalWaveFunctionU(const unsigned& uMomBin, const double& Radius,
                             const unsigned short& usCh, const unsigned short& usPW, const bool& DivideByR, const bool& Asymptotic=false) const;
    double EffectiveFunction(const unsigned& uMomBin, const double& Radius, const unsigned short& usCh);
    double EffectiveFunction(const unsigned& uMomBin, const double& Radius);

    double EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta, const unsigned short& usCh);
    double EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta);

    template <class Type> Type GetBinCenter(const Type* Bins, const unsigned& WhichBin) const;
    template <class Type> Type EvalBinnedFun(const double& xVal, const unsigned& NumBins, const double* Bins, const double* BinCent, const Type* Function) const;

    //double SourceFunction(const double* Pars, const double& GridSize);
    void UpdateCPF();
    void UpdateExternalPhaseShifts(const unsigned& usCh, const unsigned& usPW);

    //delete all variables that depend only on the number of momentum bins, b-bins and MaxPairs
    void DelMomIpMp();
    //delete all variables that depend only on the number of b-bins
    void DelIp();
    //delete all variables that depend only on the number of channels
    void DelCh();
    //delete all variables that depend on the number of momentum bins, number of channels and number of partial waves
    void DelMomChPw();
    //delete all variables that depend on the number of momentum bins and number of channels
    void DelMomCh();
    //delete all variables that depend only on the number of momentum bins
    void DelMom();
    //delete all variables that depend only on the number of momentum and b-bins
    void DelMomIp();
    //delete the potential parameters associated with a particular PW
    void DelPotPw(const unsigned short& usCh, const unsigned short& usPW);
    //delete the potential parameters associated with a particular channel
    void DelPotCh(const unsigned short& usCh);

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
    float*** PhaseShiftF;//in bins of pol/pw/mom, saved only until the end of each k-iteration
    double**** WaveFunRad;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    complex<double>**** WaveFunctionU;//in bins of mom/pol/pw/rad, saved only until the end of each k-iteration
    bool* MomBinConverged;//bins of mom, marked as true in case the num. comp. failed and this bin should not be used

    //in bins of momentum, channel, GridPoints
    double*** WaveFunction2;

    //in bins of momentum/ImpactParameter
    double** kbCorrFun;
    double** kbCorrFunErr;

    //in bins of momentum
    double* kCorrFun;
    double* kCorrFunErr;

    //!further input variables
    ////bool*** UseExternalWF;//in bins of mom/pol/pw
    //const complex<double>**** ExternalWF;//in bins of mom/pol/pw (reserved mem) / rad (provided by the user). If ExternalWF[x][y][z]=NULL => Do not use ext. wf.
    //unsigned*** NumExtWfRadBins;//in bins of mom/pol/pw
    //const double**** ExtWfRadBins;//in bins of mom/pol/pw (reserved mem) / rad (provided by the user -> BinRanges).
    //[usCh][usPW]
    DLM_Histo<complex<double>>*** ExternalWF;
    DLM_Histo<complex<double>>*** ExternalPS;

    //these are used as buffers when it comes to computing the Reference Partial Waves and the Legendre Polynomials
    //in particular, when we loop over all PWs twice, we actually evaluate the same functions multiple times => save them in an array to save CPU time
    complex<double>* RefPartWave;
    complex<double>* SolvedPartWave;
    double* LegPol;
    //the gamow correction factors (Coulomb penetration factor) pre-computed for all momentum bins
    complex<double>* CPF;
int DEBUG;
};

#endif // CATS_H
