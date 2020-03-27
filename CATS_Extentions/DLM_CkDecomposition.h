#ifndef DLM_CKDECOMPOSITION_H
#define DLM_CKDECOMPOSITION_H

#include "CATS.h"
#include "CATStools.h"
#include "DLM_ResponseMatrix.h"

#include "TH2F.h"

class DLM_Ck : public DLM_Histo<double>{

public:
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double* bins);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double& minMom, const double& maxMom);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
           const unsigned& numbin, const double* bins, double (*CorrFun)(const double&, const double*, const double*));
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
           const unsigned& numbin, const double& minMom, const double& maxMom, double (*CorrFun)(const double&, const double*, const double*));
    ~DLM_Ck();

    //change UpToData when you do that
    void SetSourcePar(const unsigned& WhichPar, const double& Value);
    double GetSourcePar(const unsigned& WhichPar);
    unsigned GetNumSourcePar();
    void SetPotPar(const unsigned& WhichPar, const double& Value);
    double GetPotPar(const unsigned& WhichPar);
    //double GetPotPar(const unsigned& WhichPar);
    unsigned GetNumPotPar();
    //the momentum value after which the correlation will be interpolated
    //N.B. if we go outside the maximum momentum, the same is assumed
    //the interpolation is done so, as the C(k) is linearly propagated between the final computed point and C(kc) = 1.
    //If kc<maximum momentum C(k>maximum momentum) = C(maximum momentum)
    //If kc<0, than outside the maximum momentum C(k) = |kc|
    void SetCutOff(const double& Momentum=1e6, const double& kc=-1);
    double GetCutOff() const;
    double GetCutOff_kc() const;

    bool Status();
    void Update(const bool& FORCE=false);
    double Eval(const double& Momentum);

private:
    const unsigned NumSourcePar;
    const unsigned NumPotPar;
    //double* MomBinCopy;
    CATS* Kitty;
    double (*CkFunction)(const double&, const double*, const double*);
    double* SourcePar;
    bool SourceUpToDate;
    double* PotPar;
    bool PotUpToDate;
    double CutOff;
    double CutOff_kc;

    void DefaultConstructor();
};


class DLM_CkDecomposition{

public:
    enum CkDecType { cFeedDown, cFake };

    DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const TH2F* hSigmaMatrix, const bool& InvertAxis=false);
//DLM_CkDecomposition():ERROR_STATE(0),NumChildren(0){}
    ~DLM_CkDecomposition();

    //! all matrices that are passed on should have the same definitions of x,y axis, else the whole thing will crash!
    //void AddFeeddown(const unsigned& WhichCk, const double& fraction, DLM_CkDecomposition* child,
    //                 TH2F* hResidualMatrix, const bool& InvertedAxis=false);
    //void AddImpurity(const unsigned& WhichCk, const double& fraction, DLM_CkDecomposition* child);
    void AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomposition* child=NULL,
                         const TH2F* hResidualMatrix=NULL, const bool& InvertedAxis=false);

    //full Ck
    double EvalCk(const double& Momentum);
    //only the part related to the main contribution (normalized to unity at large k)
    double EvalMain(const double& Momentum);
    double EvalSmearedMain(const double& Momentum);
    //excluding the missid part (normalized to unity at large k)
    double EvalMainFeed(const double& Momentum);
    double EvalSmearedMainFeed(const double& Momentum);
    //the effect of a particular child on Ck, normalized to unity at large k
    double EvalSmearedFeed(const unsigned& WhichChild,const double& Momentum);
    
    //evaluate the fully smeared correlation, returning lambda*(Ck - 1)
    double EvalSignal(const double& Momentum);
    double EvalSignalSmeared(const double& Momentum);
    double EvalSignalMain(const double& Momentum);
    double EvalSignalSmearedMain(const double& Momentum);
    double EvalSignalChild(const unsigned& WhichChild,const double& Momentum);
    double EvalSignalSmearedChild(const unsigned& WhichChild,const double& Momentum);
    
    //if(Renorm) -> spit out the C(k) normalized to 1
/*
    //Only the main correlation
    double EvalMainCk(const double& Momentum, const bool& Renorm=true);
    //Only the feed-down part
    double EvalFeedCk(const double& Momentum, const bool& Renorm=true);
    //Only the fake part
    double EvalFakeCk(const double& Momentum, const bool& Renorm=true);
    //Only a specific contribution
    double EvalSpecificCk(const unsigned& WhichCk, const double& Momentum, const bool& Renorm=true);
    //The main+feed down part
    //double EvalMainFeedCk(const double& Momentum, const bool& Renorm=true);
*/

    unsigned GetNumChildren();
    DLM_CkDecomposition* GetChild(const unsigned& WhichChild);
    DLM_CkDecomposition* GetContribution(const char* name);
    DLM_Histo<double>* GetChildContribution(const unsigned& WhichChild, const bool& WithLambda=false);
    DLM_Histo<double>* GetChildContribution(const char* name, const bool& WithLambda=false);
    DLM_Ck* GetCk();

    const TH2F* GetResolutionMatrix();

    void GetName(char* name);

    bool Status();
    //run with flag true before fitting (i.e. after you finish the full set up)
    void Update(const bool& FORCE_FULL_UPDATE=true, const bool& UpdateDecomp=true);

int DEBUGFLAG;

private:

    const bool ERROR_STATE;
    const unsigned NumChildren;
    char* Name;
    DLM_CkDecomposition** Child;
    double LambdaMain;
    //LambdaMain + fraction of feed-down
    double MuPar;
    double* LambdaPar;
    int* Type;
    //the main C(k)
    DLM_Ck* CkMain;
    //the smeared main C(k)
    DLM_Histo<double>* CkMainSmeared;
    //the C(k) containing the main and feed-down contributions
    DLM_Histo<double>* CkMainFeed;
    //the smeared Main+FeedDown C(k)
    DLM_Histo<double>* CkSmearedMainFeed;
    //the residual C(k) (corrected with the residual matrix) stemming from the CkMainFeed of the children
    DLM_Histo<double>** CkChildMainFeed;
    DLM_Histo<double>** CkSmearedChildMainFeed;
    
    //the signal (lambda*(Ck - 1)) for each contribution. All corrections are applied!
    DLM_Histo<double>* SignalMain;
    DLM_Histo<double>* SignalSmearedMain;
    DLM_Histo<double>** SignalChild;
    DLM_Histo<double>** SignalSmearedChild;
    
    bool SignalsUpdated;

    //N.B. I assume that the MomResolution is only taken into account if you look into that particular C(k) or the feed-down contributions.
    //For the children which are of type fake, the Sigma matrix is taken from their definition.
    DLM_ResponseMatrix* RM_MomResolution;
    DLM_ResponseMatrix** RM_Child;
    //this is the smearing matrix of the PARENT, but with the binning of the child
    //need for plotting
    DLM_ResponseMatrix** SM_Child;
    //status of the main C(k)
    bool CurrentStatus;
    //status of the main C(k) of the child
    bool* CurrentStatusChild;
    bool DecompositionStatus;

    bool UniqueName(const char* name);
    //bool CheckStatus();

    //void SmearOLD(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared);
    void Smear(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared);
};


#endif
