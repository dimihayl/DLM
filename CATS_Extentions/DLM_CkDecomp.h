#ifndef DLM_CKDECOMP_H
#define DLM_CKDECOMP_H

#include "CATS.h"
#include "CATStools.h"
#include "DLM_ResponseMatrix.h"

class DLM_Ck;

class DLM_CkDecomp{

public:
    enum CkDecType { cFeedDown, cFake };

    DLM_CkDecomp(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const DLM_Histo<float>* hSigmaMatrix, const bool& InvertAxis=false);
//DLM_CkDecomp():ERROR_STATE(0),NumChildren(0){}
    ~DLM_CkDecomp();

    //! all matrices that are passed on should have the same definitions of x,y axis, else the whole thing will crash!
    //void AddFeeddown(const unsigned& WhichCk, const double& fraction, DLM_CkDecomp* child,
    //                 TH2F* hResidualMatrix, const bool& InvertedAxis=false);
    //void AddImpurity(const unsigned& WhichCk, const double& fraction, DLM_CkDecomp* child);
    void AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomp* child=NULL,
                         const DLM_Histo<float>* hResidualMatrix=NULL, const bool& InvertedAxis=false);
    void AddPhaseSpace(const unsigned& WhichCk, const DLM_Histo<float>* hPhaseSpace);
    //this is for the momentum smearing
    void AddPhaseSpace(const DLM_Histo<float>* hPhaseSpace);
/*
    //if true, the data is unfolded, and later on all Eval functions applied
    //if the data has not been unfolded yet, it is performed with the default settings
    void UseUnfoldedData(const bool& unfold){}
    //performs the actual unfolding (only for the momentum resolution)
    //nboot = number of iterations for the bootstrap
    void UnfoldData(const unsigned& nboot=100);
*/
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

    double GetLambdaMain();
    double GetLambdaChild(const unsigned& WhichChild);

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
    DLM_CkDecomp* GetChild(const unsigned& WhichChild);
    DLM_CkDecomp* GetContribution(const char* name);
    DLM_Histo<double>* GetChildContribution(const unsigned& WhichChild, const bool& WithLambda=false);
    DLM_Histo<double>* GetChildContribution(const char* name, const bool& WithLambda=false);
    DLM_Ck* GetCk();

    const DLM_Histo<float>* GetResolutionMatrix();

    void GetName(char* name);

    bool Status();
    //run with flag true before fitting (i.e. after you finish the full set up)
    void Update(const bool& FORCE_FULL_UPDATE=true, const bool& UpdateDecomp=true);

int DEBUGFLAG;

protected:

    DLM_CkDecomp(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction);
    void StandardSetup(const char* name);
    void AddTheMatrix(const DLM_Histo<float>* hSigmaMatrix=NULL, const bool& InvertAxis=false);

    const bool ERROR_STATE;
    const unsigned NumChildren;
    char* Name;
    DLM_CkDecomp** Child;
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
    DLM_Histo<float>* PS_Main;
    DLM_Histo<float>** PS_Child;

    //status of the main C(k)
    bool CurrentStatus;
    //status of the main C(k) of the child
    bool* CurrentStatusChild;
    bool DecompositionStatus;

    bool UniqueName(const char* name);
    //bool CheckStatus();

    //void SmearOLD(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared);
    void Smear(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared, DLM_Histo<float>* PhaseSpace=NULL);
};


#endif
