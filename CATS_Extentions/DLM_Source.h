
#ifndef DLM_SOURCE_H
#define DLM_SOURCE_H

#include "DLM_Histo.h"
#include "CATStools.h"

class DLM_Random;

double GaussSource(double* Pars);
double GaussSourceCutOff(double* Pars);
double GaussSourceTF1(double* x, double* Pars);
double GaussSourceScaledTF1(double* x, double* Pars);

double GaussSourceTheta(double* Pars);
double CauchySource(double* Pars);
double CauchySourceTheta(double* Pars);
double CauchySource_v2(double* Pars);
double ExponentialSource(double* Pars);

double DoubleGaussSource(double* Pars);
double GaussCauchySource(double* Pars);
//double LevyIntegral1D(double* Pars);
double LevySource3D_2particle(double* Pars);
double LevySource3D_single(double* Pars);
double LevySource3D(double* Pars);
//double LevySource_A(double* Pars);
//a monte-carlo out-side-long Gaussian source. Works very slowly!
double GaussOSL_MC(double* Pars);

double Gauss_Exp_Approx(double* Pars);
double Gauss_Exp(double* Pars);
double GaussExpSimple(double* Pars);
double GaussExpTotSimple(double* Pars);
double GaussExpTotIdenticalSimple(double* Pars);
double GaussExpTotIdenticalSimple_2body(double* Pars);
double GaussExpTotSimple_2body(double* Pars);

/*
class MS_Gauss:public CatsSource{
public:
    //MS_Gauss(){}
    //~MS_Gauss(){}
    double Size;
    void SetParameter(const unsigned& WhichPar, const double& Value);
    double Eval(double* Pars);
};
*/

//A source taking into account resonances and mT scaling. The Simple part is that resonances are back to back
//and in case of two resonances, the t*p/m are just added up. Also we use the approximate relation for small t*p/m
class MS_GaussExp_mT_Simple:public CatsSource{
public:
    MS_GaussExp_mT_Simple();
    ~MS_GaussExp_mT_Simple();

    void SetNum_mT(const unsigned& nmt);
    void SetMean_mT(const unsigned& umt, const double& mmt);
    void SetWeight_mT(const unsigned& umt, const double& wmt);
    void SetLinear_mT(const double& lin);
    void SetSlope_mT(const double& slope);
    void SetCustomFunction(const unsigned& umt, const double& value);
    void RemoveCustomFunction();
    void SetMass(const unsigned short& particle, const double& mass);
    void SetMassR(const unsigned short& particle, const double& mass);
    void SetMassD(const unsigned short& particle, const double& mass);
    void SetTau(const unsigned short& particle, const double& tau);
    void SetResonanceWeight(const unsigned short& particle, const double& weight);

    void SetParameter(const unsigned& WhichPar, const double& Value);
    double Eval(double* Pars);
    double EvalROOT(double* x, double* Pars);
    unsigned GetNumPars();
private:
    unsigned Num_mT;//number of mT bins (common for the particle pair)
    //[mT]
    double* Mean_mT;//the mean mT in each bin
    double* Weight_mT;//the weight of each mT bin
    double Linear_mT;//a+b*mT = r, linear is a, slope is b
    double Slope_mT;
    //it this is not NULL, we use the values here instead of a linear function
    double* FunctionValue;
    //[particle 1,2]
    double* Mass;//mass of the main particles we investigate
    double* MassR;//mass of the resonance from which each particles comes from
    double* MassD;//mass of the second daughter in the decay of the resonance
    double* Tau;//the mean lifetime of the resonance
    double* Weight_R;//the amount of resonances

    double* Parameters;
};

class DLM_StableDistribution:public CatsSource{
public:
    DLM_StableDistribution(const unsigned& numgridpts=512*2);
    ~DLM_StableDistribution();
    void SetStability(const double& val);
    void SetLocation(const double& val);
    void SetScale(const double& val);
    void SetSkewness(const double& val);
    void SetNumIter(const unsigned& val);
    void SetParameter(const unsigned& WhichPar, const double& Value);
    double Eval(double* Pars);
    unsigned GetNumPars();
private:
    const unsigned NumGridPts;
    double Stability;
    double Skewness;
    double Scale;
    double Location;
    unsigned NumIter;
    bool Generated;
    void Generate(const double& stability, const double& location, const double& scale, const double& skewness);
    DLM_Histo1D<double>* Histo;
    DLM_Random* RanGen;
};

class DLM_CleverLevy:public CatsSource{
public:
    DLM_CleverLevy();
    ~DLM_CleverLevy();
    double Eval(double* Pars);
    double RootEval(double* x, double* Pars);
    void InitStability(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitScale(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitRad(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitType(const int& type);
    unsigned GetNumPars();
private:
    //0 = single particle
    //1 = pair
    //2 = Nolan notation
    int Type;
    unsigned NumPtsStability;
    unsigned NumPtsScale;
    unsigned NumPtsRad;
    double MinStability;
    double MaxStability;
    double MinScale;
    double MaxScale;
    double MinRad;
    double MaxRad;
    DLM_Histo<double>* Histo;
    void Reset();
    void Init();
};

class DLM_CleverMcLevyReso:public CatsSource{
public:
    DLM_CleverMcLevyReso();
    ~DLM_CleverMcLevyReso();
    double Eval(double* Pars);
    double RootEval(double* x, double* Pars);

    void InitStability(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitScale(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitRad(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitType(const int& type);
    //this should be called individually for both of the particles in the pair
    void InitReso(const unsigned& whichparticle, const unsigned& numreso);
    //momSmear is in percent
    //the smearing is ONLY applied to the original resonance
    //additional information on the decay topology of the resonance
    //rdtBackwards = emission only of back-to-back particles (angle = 180 deg)
    //rdtRandom = uniform emission of particles (0-180 deg)
    //rdtRandomBackwards = uniform emission of particles, but particles are only allowed to move away from one another (90-180)
    enum RESO_DEC_TOP { rdtBackwards, rdtRandom, rdtRandomBackwards };
    void SetUpReso(const unsigned& whichparticle, const unsigned& whichreso, const double& weight, const double& mass, const double& tau, const double& mass0, const double& mass1,
                   const double& momSmear=0, const bool& massSmear=false, const RESO_DEC_TOP& rdt=rdtBackwards);
    void SetUpResoEmission(const unsigned& whichparticle, const unsigned& whichreso, const DLM_Histo<double>* Distr);
    void InitNumMcIter(const unsigned& numiter);
    unsigned GetNumPars();
private:
    //0 = single particle
    //1 = pair
    //2 = Nolan notation
    int Type;
    unsigned NumPtsStability;
    unsigned NumPtsScale;
    unsigned NumPtsRad;
    double MinStability;
    double MaxStability;
    double MinScale;
    double MaxScale;
    double MinRad;
    double MaxRad;

    unsigned* NumResonances;
    double** ResoWeight;
    double** ResoMass;
    double** ResoTau;
    //ChildMass0 is the 'primary' particle of interest
    //we assume that the decay is two-body and that the final k* is zero. Base on that we can get an estimate for the
    //momentum of the resonance, which we will use
    double** ChildMass0;
    double** ChildMass1;
    double** SmearResoMomentum;
    bool** SmearResoMass;
    RESO_DEC_TOP** ResoDecayTopology;
    //the angle between r and k of the resonance
    //this should be the cumulative distribution of the probability of emission at a certain angle (in RAD between 0 and Pi),
    //where the axis are inverted, i.e. y(x)
    const DLM_Histo<double>*** ResoEmissionAngle;
    //the number of MC iterations with which each source is to be initialized
    unsigned NumMcIter;
    DLM_Histo<double>* Histo;
    void Reset();
    void Init();
};

//with input from transport model to evaluate modification in r
class DLM_CleverMcLevyResoTM:public CatsSource{
public:
    DLM_CleverMcLevyResoTM();
    ~DLM_CleverMcLevyResoTM();
    double Eval(double* Pars);
    double RootEval(double* x, double* Pars);
    double RootEvalNorm(double* x, double* Pars);

    void InitStability(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitScale(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitRad(const unsigned& numPts, const double& minVal, const double& maxVal);
    void InitType(const int& type);
    void SetUpReso(const unsigned& whichparticle, const double& weight);
    //EVERYTHING BELOW IS IN THE CM OF THE DAUGHTERS
    //bgt = beta*gamma*tau
    //a_cp = angle between the r_core and the momentum of the resonance
    void AddBGT_PR(const float& bgt,const float& a_cp);
    void AddBGT_RP(const float& bgt,const float& a_cp);
    //the same quantities for both resonances, as well as the angle between the momenta of the two resonances.
    void AddBGT_RR(const float& bgt0,const float& a_cp0,const float& bgt1,const float& a_cp1,const float& a_p0p1);
    void InitNumMcIter(const unsigned& numiter);
    unsigned GetNumPars();
    void SetNormalization(const double& norm);
private:
    //0 = single particle
    //1 = pair
    //2 = Nolan notation
    int Type;
    unsigned NumPtsStability;
    unsigned NumPtsScale;
    unsigned NumPtsRad;
    double MinStability;
    double MaxStability;
    double MinScale;
    double MaxScale;
    double MinRad;
    double MaxRad;
    double Normalization;

    double* ResoWeight;
    //const DLM_Histo<double>* ResoEmissionPR;
    //const DLM_Histo<double>* ResoEmissionRP;
    //const DLM_Histo<double>* ResoEmissionRR;

    //a list of possible beta*gamma*tau
    //for primordial(P), reso (R)
    //[#entry number][0] = bgt0
    //[#entry number][1] = bgt1
    //[#entry number][2] = a_cp0
    //[#entry number][3] = a_cp1
    //[#entry number][4] = a_p0p1
    //N.B. in case one particle is primordial and one reso, than the BGT of the reso is clear,
    //i.e. beta*gamma*tau of the reso itself. An approximation that one could consider: set beta*gamma of the prim is lower
    //then that of the reso (due to the lower momentum), hence we neglect this term and set it to 0
    unsigned MaxBGT_PR;
    unsigned NumBGT_PR;
    float** BGT_PR;
    //for particle0 reso, particle1 primordial
    unsigned MaxBGT_RP;
    unsigned NumBGT_RP;
    float** BGT_RP;
    //for particle0 reso, particle1 reso
    unsigned MaxBGT_RR;
    unsigned NumBGT_RR;
    float** BGT_RR;

    //the number of MC iterations with which each source is to be initialized
    unsigned NumMcIter;
    DLM_Histo<double>* Histo;
    void Reset();
    void Init();
};


/*
//use DLM_CleverMcLevyReso as the baseline for a differential analysis
class DLM_CleverMcLevyReso_Diff:public DLM_CleverMcLevyReso{
public:
    DLM_CleverMcLevyReso_Diff();
    ~DLM_CleverMcLevyReso_Diff();
    enum DiffType { dtConst, dtLiniar, dtSquare, dtExp };

    double Eval(double* Pars);
    double RootEval(double* x, double* Pars);
    void InitNumDiffBins(const unsigned& numbins);//!
    //The type of computation:
    //0 is the case where the fit is done independently for each bin
    //else we use two digit number, the first describing the scale parameter, the second the stability, where:
    // 1 is linear scaling; 2 is pol2;
    void InitDiffType(const int& type);//!
    unsigned GetNumPars();
private:
    //0 = single particle
    //1 = pair
    //2 = Nolan notation
    int Type;
    unsigned NumPtsStability;
    unsigned NumPtsScale;
    unsigned NumPtsRad;
    double MinStability;
    double MaxStability;
    double MinScale;
    double MaxScale;
    double MinRad;
    double MaxRad;
    unsigned* NumResonances;
    double** ResoWeight;
    double** ResoMass;
    double** ResoTau;
    //ChildMass0 is the 'primary' particle of interest
    //we assume that the decay is two-body and that the final k* is zero. Base on that we can get an estimate for the
    //momentum of the resonance, which we will use
    double** ChildMass0;
    double** ChildMass1;
    //the number of MC iterations with which each source is to be initialized
    unsigned NumMcIter;
    DLM_Histo<double>* Histo;
    void Reset();
    void Init();
}
*/
#endif
