#ifndef DLM_FITTERS_H
#define DLM_FITTERS_H

//#include "CATS.h"
//#include "DLM_SmearedCats.h"

#include "TString.h"

//#include "DLM_CkDecomposition.h"

//this typedef it used to save the potentials for the
//different channels as an array of function pointers.
//typedef double (*CatsAnaSource)(double*);

class DLM_CkDecomposition;
class TGraph;
class TH1F;
class TF1;

//this typedef it used to save the functions for the
//different baselines as an array of function pointers.
typedef double (*FitterFunction)(double*,double*);

double DLM_FITTER2_FUNCTION_POL(double* xVal, double* pars);
double DLM_FITTER2_FUNCTION_SPLINE3(double* xVal, double* pars);

//!FOR TEST!
//#include "TCanvas.h"
//#include "TFile.h"

class DLM_Fitter1{

public:
    //(a+b*x+c*x*x)*C(k)
    enum fFitPar { p_a, p_b, p_c, p_3, p_4, p_ab_0, p_ab_1, p_ab_2, p_ab_3, p_ab_4, p_ab_5, p_ab_6, p_Cl, p_kc,
    p_sor0, p_sor1, p_sor2, p_sor3, p_sor4, p_sor5, p_pot0, p_pot1, p_pot2, p_pot3, p_pot4, p_pot5, p_spline  };
    enum fFitRange { kmin, kf, kl, kmax };

    DLM_Fitter1(const unsigned& maxnumsyst=16);
    ~DLM_Fitter1();

    //sets up a certain system that we want to fit
    void SetSystem(const unsigned& WhichSyst, const TH1F& histo, const double& FromMeV ,
                   DLM_CkDecomposition& decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX);
    int GetSystem(const TString& System);
    TString GetSystem(const unsigned& WhichSyst);
    bool ChangeCkModel(const unsigned& WhichSyst, DLM_CkDecomposition& decomp);
    //in case there are multiple systems sharing a source this should be set here
    void AddSameSource(const TString& System, const TString& EqualTo, const int& numpars);
    void RemoveSameSource(const TString& System);
    void AddSamePotential(const TString& System, const TString& EqualTo, const int& numpars);
    void RemoveSamePotential(const TString& System);
    void AddSameParameter(const TString& System, const unsigned& WhichPar, const TString& ParentSystem, const unsigned& ParentPar);
    void RemoveSameParameter(const TString& System, const unsigned& WhichPar);
    //returns not just the parent, but parent of the parent of the .....
    int GetBaseParameter(const int& WhichSyst, const int& WhichPar);
    int GetBaseParameter(const int& WhichSyst, const int& WhichPar, int& ParentSystem, int& ParentPar);
    int GetBaseParameter(const TString& System, const int& WhichPar);
    int GetBaseParameter(const TString& System, const int& WhichPar, TString& ParentSystem, int& ParentPar);

    //0 is a polynomial of type a+bx+...
    //1 is a polynomial of type a*(1+bx+...)
    void SetMultBaselineType(const int& TYPE);
    void SetAddBaselineType(const int& TYPE);

    //
    //void SetSameSourceChildren(const TString& WhichSyst, const TString& EqualTo, const int& NumPars);

    void SetParameter(const unsigned& WhichSyst, const unsigned& WhichPar, const double& Value, const double& ValueDown=-1e64, const double& ValueUp=1e64);
    void SetParameter(const TString& WhichSyst, const unsigned& WhichPar, const double& Value, const double& ValueDown=-1e64, const double& ValueUp=1e64);
    void FixParameter(const unsigned& WhichSyst, const unsigned& WhichPar, const double& Value);
    void FixParameter(const TString& WhichSyst, const unsigned& WhichPar, const double& Value);
    //fix the parameter to the current value
    void FixParameter(const unsigned& WhichSyst, const unsigned& WhichPar);
    void FixParameter(const TString& WhichSyst, const unsigned& WhichPar);
    double GetParameter(const unsigned& WhichSyst, const unsigned& WhichPar);
    double GetParError(const unsigned& WhichSyst, const unsigned& WhichPar);
    double GetParameter(const TString& WhichSyst, const unsigned& WhichPar);
    double GetParError(const TString& WhichSyst, const unsigned& WhichPar);
    double GetChi2();
    int GetNdf();
    double GetChi2Ndf();
    double GetPval();
    //double Eval(const unsigned& WhichSyst, const double& Momentum);
    void GetFitGraph(const unsigned& WhichSyst, TGraph& OutGraph, const bool& DataBinning=true);
    void GetMultBaselineGraph(const unsigned& WhichSyst, TGraph& OutGraph, const bool& DataBinning=true);
    void GetCkDecompGraph(const unsigned& WhichSyst, TGraph& OutGraph);
    void GetCkTheoryGraph(const unsigned& WhichSyst, TGraph& OutGraph);

    void SetOutputDir(const TString& outdirname);
    void SetSeparateBL(const unsigned& WhichSyst, const bool& yesno);
    //if true, there is not really a baseline region, in a sense that Ck is evaluated fully in the baseline region as well.
    //this was build in case a region of the fit needs to be cut-off, e.g. a peak
    //N.B. if set to yes, it overrides the SeparateBL!
    void SetFullCkForBaseline(const unsigned& WhichSyst, const bool& yesno);

    void RemoveNegativeCk(const bool& yesno);
    bool CheckNegativeCk();

    const unsigned GetNumParPerSyst(){return NumPar;}

    TF1* GetFit();
    TF1* GetBaselineFit(const unsigned& WhichSyst) const;
    const TH1F* GetGlobalHisto() const;


    //set up the global histogram and perform the fit
    void GoBabyGo(const bool& show_fit_info=false);

//void TEST1(const unsigned& WhichSyst, TH1F* histo, const double& FromMeV ,
//                   DLM_CkDecomposition* decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX);

    //double TESTEVAL(double xVal, double* PARS) {return EvalGlobal(&xVal, PARS);}


private:
    const unsigned MaxNumSyst;
    const unsigned NumPar;
    //const unsigned NumPotPar;
    const unsigned NumRangePar;
    bool RemoveNegCk;
    int TypeMultBl;
    int TypeAddBl;

    TString OutputDirName;

    const TH1F** HistoOriginal;
    TH1F** HistoToFit;

    //[WhichSyst]
    DLM_CkDecomposition** SystemToFit;
    unsigned NumSourceSystems;
    //an array containing all main Ck's + all secondaries that have their source fixed to one of the mains
    DLM_CkDecomposition** SourceSystems;
    //the ID in the main array of the Ck to which the source of a system is fixed
    int* ParentSource;
    int* ParentParameter;

    unsigned NumPotentialSystems;
    //an array containing all main Ck's + all secondaries that have their source fixed to one of the mains
    DLM_CkDecomposition** PotentialSystems;
    //the ID in the main array of the Ck to which the source of a system is fixed
    int* ParentPotential;

    //[WhichSyst][0-3] 0-3: kMin,kFemto,kLinear,kMax
    double** FitRange;

    //[Entry][WhichSystem/TemplateSystem]
    TString** SameSourceMap;
    int* NumSameSourcePar;
    unsigned NumSourceMapEntries;

    TString** SamePotentialMap;
    int* NumSamePotentialPar;
    unsigned NumPotentialMapEntries;

    //for the future. Might be useful for UrFAT
    //TString** SameInterMap;
    //int* NumSameInterMap;

    double** ParValue;
    double** ParDownLimit;
    double** ParUpLimit;
    bool** FixPar;
    bool* SeparateBaseLineFit;
    bool* FullCkForBaseline;

    TH1F* HistoGlobal;
    TF1* FitGlobal;
    TF1** FitBL;

    double* GlobalToMomentum;
    unsigned* NumBinsSyst;
    unsigned* CumulativeNumBinsSyst;

    bool ShowFitInfo;

    //unsigned SourceAnchoredTo(const);

    double EvalGlobal(double* xVal, double* Pars);
    int GetBaseParameter(const int& WhichSyst, const int& WhichPar, int& ParentSystem, int& ParentPar, const int& StartSystem, const int& StartPar);
    int GetBaseParameter(const TString& System, const int& WhichPar, TString& ParentSystem, int& ParentPar, const TString& StartSystem, const int& StartPar);

};





class DLM_Fitter2{

public:
    //AddBl(k)*C(k)+MultBl(k)
    enum fFitPar { p_Ck, p_sor, p_pot, p_AddBl, p_MultBl };
    enum fFitRange { kmin, kf, kl, kmax };
    //Default functions for the baseline, information:
    //f_pol:        Info = Order of the polynomial
    //              [0] = Order of the polynomial (N) (SHOULD BE ALWAYS FIXED!!!)
    //              [1...N+1] = the parameter for the 0...N order
    //f_spline3:    Info = Number of knots in the spline
    //              [0] = Number of knots in the spline (N)
    //              [1...N] = the x-values of the knots
    //              [N+1...2N] = the y-value of the knots
    //              [2N+1] = the value of the derivative at the first knot
    //              [2N+2] = the value of the derivative at the last knot
    enum fBlFunction { f_pol, f_spline3, f_na=-1 };
    enum fBlType { b_add, b_mult };

    DLM_Fitter2(const int& maxnumsyst=16);
    ~DLM_Fitter2();

    int GetSystemID(const TString& WhichSyst);

    //set up the baseline function, by passing a pointed
    void SetUpAddBl(const TString& WhichSyst, double (*FUNCTION)(double*,double*),const int& numpars);
    //setup the baseline function for the predefined type. Info is specific for each type and my have different meaning
    void SetUpAddBl(const TString& WhichSyst, const fBlFunction& TYPE, const int& Info);
    void SetUpMultBl(const TString& WhichSyst, double (*FUNCTION)(double*,double*),const int& numpars);
    void SetUpMultBl(const TString& WhichSyst, const fBlFunction& TYPE, const int& Info);

    void SetPolParX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& Value);
    void FixPolParX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& Value);
    void SetLimitPolParX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& down, const double& up);
    void SetPolParY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& Value);
    void FixPolParY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& Value);
    void SetLimitPolParY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichPar, const double& down, const double& up);

    void SetSplineKnotX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& Value);
    void FixSplineKnotX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& Value);
    void SetLimitSplineKnotX(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& down, const double& up);
    void SetSplineKnotY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& Value);
    void FixSplineKnotY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& Value);
    void SetLimitSplineKnotY(const TString& WhichSyst, const fBlType& WhichBl, const int& WhichKnot, const double& down, const double& up);

private:

    const int MaxNumSyst;

    //Number of parameters per system
    int* NumPars;
    int* NumCkPars;
    int* NumSorPars;
    int* NumPotPars;
    int* NumAddBlPars;
    int* NumMultBlPars;

    //the position of the first parameter of the corresponding system, within the lists below
    int* StartPars;
    int* StartCkPars;
    int* StartSorPars;
    int* StartPotPars;
    int* StartAddBlPars;
    int* StartMultBlPars;

    //all parameters. To get one associated with a particular system, use the StartPars above
    int CurrentNumPars;
    int MaxNumPars;
    double* Pars;
    double* ParsErr;
    double* ParsDownLim;
    double* ParsUpLim;
    /*
    double* Pars_Ck;
    double* ParsErr_Ck;
    double* ParsDownLim_Ck;
    double* ParsUpLim_Ck;
    double* Pars_Sor;
    double* ParsErr_Sor;
    double* ParsDownLim_Sor;
    double* ParsUpLim_Sor;
    double* Pars_Pot;
    double* ParsErr_Pot;
    double* ParsDownLim_Pot;
    double* ParsUpLim_Pot;
    double* Pars_AddBl;
    double* ParsErr_AddBl;
    double* ParsDownLim_AddBl;
    double* ParsUpLim_AddBl;
    double* Pars_MultBl;
    double* ParsErr_MultBl;
    double* ParsDownLim_MultBl;
    double* ParsUpLim_MultBl;
    */

    //double (*Function_AddBl)(double* xVal, double* pars);
    //double (*Function_MultBl)(double* xVal, double* pars);

    FitterFunction* Function_AddBl;
    FitterFunction* Function_MultBl;
    fBlFunction* FunctionType_AddBl;
    fBlFunction* FunctionType_MultBl;

    DLM_CkDecomposition** SystemToFit;

    void IncreaseMaxNumPars();
};


/*
CATS* FIT_CATS1; CATS* FIT_CATS2; CATS* FIT_CATS3; CATS* FIT_CATS4;

//!Info for the global fits
//to fit globally, we assume that there is a single histogram to fit. This histogram contains all bins of the pp and pL
//data, where the first set of bins are those corresponding to pp. The x-val is in arbitrary units and takes values from 0 to 1.
//the BorderValue (the [0]th parameter) is the value which separates pp from pL bins
//in order to convert the x-axis to momentum, we need to know the minimum and maximum momentum for pp and pL,
//those are set in the [1] and [2] fit parameters

//The fitting function: N*(a+b*k)*C(k,r)
//pars:
//[0] is the BorderValue (below which we fit pp) --> FIX
//[1] is kMin pp, [2] is kMin_pL --> FIX
//[2] is kMax pp, [3] is kMax_pL --> FIX
//[4] a pp, [5] b pp, [6] Norm pp, [7] Radius pp
//[8] a pL, [9] b pL, [10] Norm pL, [11] Radius pL
double DLM_GlobalFit1_pp_pL(double* mom, double* par){
    double Momentum;
    double& histoX = *mom;
    double& Border = par[0];
    double& kMin_pp = par[1];
    double& kMin_
    double& kMax_pp = par[1];
    double& kMax_pL = par[2];
    //pp
    if(histoX<Border){
        Momentum = histoX/Border*kMax_pp;
        FIT_KITTY_pp->SetAnaSource(0,par[6]);
        //needed for the residual
        FIT_KITTY_pL->SetAnaSource(0,par[10]==-1?par[6]:par[10]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 0, par[11]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 1, par[12]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 2, par[13]);
        if(!FIT_KITTY_pL->CkStatus()){
            FIT_KITTY_pL->KillTheCat();
            FIT_SMEARED_KITTY_pL->Correct(false);
        }
        if(!FIT_KITTY_pp->CkStatus()){
            FIT_KITTY_pp->KillTheCat();
            FIT_SMEARED_KITTY_pp->Correct(false);
        }
        return par[5]*(par[3]+par[4]*Momentum)*FIT_SMEARED_KITTY_pp->EvalCorrectedCk(Momentum);
    }
    //pL
    else{
        Momentum = (histoX-Border)/(1.-Border)*kMax_pL;
        FIT_KITTY_pL->SetAnaSource(0,par[10]==-1?par[6]:par[10]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 0, par[11]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 1, par[12]);
        FIT_KITTY_pL->SetShortRangePotential(0, 0, 2, par[13]);
        if(!FIT_KITTY_pL->CkStatus()){
            FIT_KITTY_pL->KillTheCat();
            FIT_SMEARED_KITTY_pL->Correct(false);
        }
        return par[9]*(par[7]+par[8]*Momentum)*FIT_SMEARED_KITTY_pL->EvalCorrectedCk(Momentum);
    }

    return 0;
}
*/

#endif
