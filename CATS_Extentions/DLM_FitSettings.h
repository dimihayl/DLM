#ifndef DLM_FITSETTINGS_H
#define DLM_FITSETTINGS_H

#include "TString.h"
#include <complex>

class CATSparameters;

using namespace std;

//forward declaration
template<typename Type> class CATShisto;
class DLM_Ck;
class CATS;
class TH1F;
class TH2F;

//the settings specific for each particle
class FemtoParticle{
public:
    FemtoParticle(const TString& name, const unsigned& numfeed, const unsigned& nummisid);
    ~FemtoParticle();
    void SetFeed(const unsigned& WhichOne, FemtoParticle& Particle);
    void SetMisid(const unsigned& WhichOne, FemtoParticle& Particle);
    bool NameIsAvailable(const TString& name) const;
    TString GetName() const;
    //only within this object
    unsigned GetNumFeed() const;
    //returns the total number of contributors to this correlation function. This includes the main fraction,
    //all feed-downs (iteratively), and the feed objects related to the misid branch of the mother.
    //Nevertheless the misid of the higher branches (children) is not taken into account.
    //the reason is that in an exp. we can never measure the misid branch of a feed-down channel.
    unsigned GetNumContrib() const;
    unsigned GetNumMisid() const;

    int GetFeedID(const TString& name) const;
    //void SetMass(const double& mass);
    //double GetMass();
private:
    const TString Name;
    //const double Mass;
    const unsigned NumFeedChannels;
    const unsigned NumMisidChannels;
    FemtoParticle** FeedParticle;
    FemtoParticle** MisidParticle;

    //there just for protecting outside call of GetTotNumFeed(true/false)
    unsigned GetNumContrib(const bool& IncludeMisid) const;
};

class FemtoExperiment{
public:
    FemtoExperiment(const TString& name, const unsigned& numpart);
    ~FemtoExperiment();

    void SetParticle(const unsigned& WhichOne, const FemtoParticle& particle);
    void SetPurity(const unsigned& WhichOne, const double& pur);
    void SetPurity(const TString& WhichOne, const double& pur);

    void SetFraction(const unsigned& WhichMother, const unsigned& WhichFeed, const double& feed);
    void SetFraction(const TString& WhichMother, const TString& WhichFeed, const double& feed);

    int FindParticle(const TString& WhichOne) const;

    TString GetName() const;

private:
    const TString Name;
    //number of primary particles which are to be analyzed
    const unsigned NumParticles;
    const FemtoParticle** Particle;
    double* Purity;
    double** Fraction;


};

class FemtoPair{
public:
    FemtoPair(const FemtoParticle& part1, const FemtoParticle& part2);
    ~FemtoPair();

    //check if this is our pair, ordering of the names does not matter
    bool IsIt(const TString& Name1, const TString& Name2) const;

    //if hfeed=NULL => flat contribution!
    void SetFeedDownMatrix(const TString& Mother1, const TString& Mother2,
                           const TString& Child1, const TString& Child2,
                           TH2F* hfeed);

    void SetFeedDownMatrix(const TString& Mother1, const TString& Mother2,
                           TH2F* hfeed);

    const FemtoParticle* GetParticle(const TString& NAME) const;
    const FemtoParticle* GetParticle(const unsigned& WhichOne) const;

private:
    const FemtoParticle& Particle1;
    const FemtoParticle& Particle2;
    //it is actually num feed down contributions + the main contribution,
    //so NumFeed+1 in a sense
    const unsigned NumFeed1;
    const unsigned NumFeed2;
    TH2F*** hFeedDown;

};

class FemtoExpPair{
public:
    FemtoExpPair(const FemtoPair& ppair, const FemtoExperiment& experiment);
    ~FemtoExpPair();

    void SetFemtoRegion(const double& MIN, const double& MAX);
    void SetBlRegion(const double& MIN, const double& MAX);
    //void AddCorrelationModel(DLM_Ck& CorrMod);
    void SetData(TH1F* data);
    void SetResolutionMatrix(TH2F* hres);
    void SetLargeMomCk(const int& TYPE); enum cktype { ckLinear, ckFlat, ckUnity };
    void SetBlType(const int& TYPE); enum bltype { blPol0, blPol1, blPol2, blSeparatePol0, blSeparatePol1, blSeparatePol2};

    //complementary to SetCorrelationModel
    void SetStandardInteraction(const TString& Inter, TH1F* TemplateData=NULL);
    TString GetStandardInteraction();

    const FemtoExperiment* GetExperiment() const;
    const FemtoParticle* GetParticle(const unsigned& WhichOne) const;

private:
    const FemtoPair& PartPair;
    const FemtoExperiment& Experiment;

    double FemtoRegionLow;
    double FemtoRegionUp;
    //how to model the C(k) in the BL region: choose between a linear function, 1, and a constant
    int LargeMomCk;
    double BlRegionLow;
    double BlRegionUp;
    //how to model the baseline. Choose between pol0,1,2, as well as the same options but performing separate fitting of the BL
    int BlType;

    CATS* Kitty;
    CATSparameters* CatSourcePars;
    //double** CatPotPars;
    CATSparameters** CatPotPars;
    complex<double>**** WaveFunctionU; double*** PhaseShifts; double* RadBins; unsigned NumRadBins;
    DLM_Ck* CorrelationModel;
    TString StandardInter;

    TH2F* hResolution;

    TH1F* hData;

    //if != NULL => add the errors to the data quadratically
    //TH1F* hSystErr;

    //double* FemtoRegion;
    //unsigned NumFemtoRegion;//number of variations to be done
    //double* BaselineRegion;
    //unsigned NumBaselineRegion;
    //DLM_Ck* CorrelationModel;
    //unsigned NumCorrelationModel;

    //TH1F* hResolution;
    //TH1F** hFeedDown;
    //TH1F** hData;
    //unsigned NumData;
    //TH1F* hSystErr;
    //int SystErrType
};


class FemtoGlobalFit{
public:
    FemtoGlobalFit();
    ~FemtoGlobalFit();



private:

};

#endif
