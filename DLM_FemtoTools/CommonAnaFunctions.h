#ifndef COMMONANAFUNCTIONS_H
#define COMMONANAFUNCTIONS_H

#include <iostream>

template <class Type> class DLM_Histo;

class TString;
class TH1F;
class TH2F;
class TF1;
class CATS;
class DLM_CleverLevy;
class DLM_CleverMcLevyReso;
class DLM_CleverMcLevyResoTM;
class DLM_Ck;

//class MS_GaussExp_mT_Simple;

const double Mass_pi0 = 134.9766;
const double Mass_pic = 139.57018;
const double Mass_Kch = 493.677;
const double Mass_K0 = 497.648;
const double Mass_phi = 1019.461;
const double Mass_p = 938.272;
const double Mass_n = 939.565;
const double Mass_L = 1115.683;
const double Mass_Xim = 1321.7;
const double Mass_Xi0 = 1314.86;
const double Mass_Xim1530 = 1535;
const double MassOmega = 1672.45;
const double Mass_d = 1875.613;
const double Mass_Dch = 1869.62;
const double Mass_D0 = 1864.84;
const double Mass_Dch_star = 2010;
const double Mass_D0_star = 2007;


class DLM_CommonAnaFunctions{

public:

    DLM_CommonAnaFunctions();
    ~DLM_CommonAnaFunctions();

    //! ALWAY CALL THESE FUNCTIONS AFTER YOU HAVE DEFINED THE MOMENTUM BINS!
    //SOURCE:
    //"Gauss"
    //"Cauchy"
    //"Levy_Nolan"
    //"Levy_Single"
    //"Levy_Diff"
    //"CleverLevy_Nolan"
    //"CleverLevy_Single"
    //"CleverLevy_Diff"
    //"GaussExpTotSimple_2body" (the first version)
    //"McLevyNolan_Reso" (the Monte-Carlo version without mT scaling)
    //"EPOS"
    //"EPOSrescaled" -> starts with a basis rescaling of 1.5
    //"Levy_mT_Reso" (the MC version created for pLambda analysis)
    //POT:
    //  "AV18"
    void SetUpCats_pp(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);
    void SetUpCats_ppic(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=200);
    void SetUpCats_pipi(CATS& Kitty, const TString& SOURCE, const int& SourceVar);
    void SetUpCats_pipi_broken(CATS& Kitty, const TString& SOURCE, const int& SourceVar);
    //POT:
    //  "LO"
    //  "NLO"
    //  "NLO_Coupled_S"
    //  "Usmani"
    //no potential variations at the moment
    //the source variation is at the moment only relevant for the McReso sources
    //  SourceVar is considered 3 digit, the first two are for the momentum smearing, the second for the mass smearing
    //   the first two:
    //    0 = no smear; else the smearing in percent (up to 99%)
    //   the second:
    //    0 = no smear; 1 = smear according to the life-time
    void SetUpCats_pL(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);
    void SetUpCats_pS0(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);
    //POT:
    //  "pXim_Lattice" (the first version)
    //  "pXim_HALQCD1" (the second version, THE ONE TO USE)
    void SetUpCats_pXim(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);
    void SetUpCats_pXi0(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);
    void SetUpCats_pOmegam(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar=0, const int& SourceVar=0);

    DLM_Ck* SetUpLednicky_pL(const unsigned& NumMomBins, const double* MomBins,  const TString& POT);

    void SetUpBinning_pp(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion, const int& MomBinVar=0, const int& FitRegVar=0);
    void SetUpBinning_pL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion, const int& MomBinVar, const int& FitRegVar);

    //DataSamples: SystemEnergy_Trigger_Version
    //the version is there to mark the different versions based on our own analysis, it can be some short description
    //Versions:
    //  Run2paper: as used for all of the first Run2 papers (LL, pXim etc)
    //DataSamples:
    //  pp13TeV_MB_Run2paper
    //  pp13TeV_HM_March19
    //  pPb5TeV_Run2paper
    //  pPb5TeV_CPR_Mar19 (with close pair rejection, as in end of March 2019)
    //The Variation flag is there for the systematics, refer to the functions themselves for more information
    void GetPurities_p(const TString& DataSample, const int& Variation, double* Purities);
    void GetPurities_L(const TString& DataSample, const int& Variation, double* Purities);
    void GetPurities_Xim(const TString& DataSample, const int& Variation, double* Purities);
    void GetFractions_p(const TString& DataSample, const int& Variation, double* Fractions);
    void GetFractions_L(const TString& DataSample, const int& Variation, double* Fractions);
    void GetFractions_Xim(const TString& DataSample, const int& Variation, double* Fractions);
    //primary, pL->pp, XX->pp, pp fake
    void SetUpLambdaPars_pp(const TString& DataSample, const int& Variation_p, double* lambda_pars);
    //primary, pS0->pL, pXim->pL, XX->pL, pp fake
    void SetUpLambdaPars_pL(const TString& DataSample, const int& Variation_p, const int& Variation_L, double* lambda_pars);
    void SetUpLambdaPars_pXim(const TString& DataSample, const int& Variation_p, const int& Variation_Xim, double* lambda_pars);

    void SetCatsFilesFolder(const TString& folder);

    TH2F* GetResolutionMatrix(const TString& DataSample,const TString&& System);
    TH2F* GetResidualMatrix(const TString&& FinalSystem, const TString& InitialSystem);
    TH1F* GetAliceExpCorrFun(const TString& DataSample,const TString& System,const TString& CutVar,const int& iReb, const bool& AddSyst=false,const int mTbin=-1);

    DLM_CleverMcLevyReso* GetCleverMcLevyReso_pp();
    DLM_CleverMcLevyReso* GetCleverMcLevyReso_pL();
    DLM_CleverMcLevyReso* GetCleverMcLevyReso_pXim();
    DLM_CleverMcLevyReso* GetCleverMcLevyReso_pOmegam();

    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_pp();
    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_pL();
    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_pXim();
    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_pOmegam();
    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_pipi();
    DLM_CleverMcLevyResoTM* GetCleverMcLevyResoTM_ppic();


    DLM_CleverMcLevyResoTM* GaussCoreRsm_LK(const int& SourceVar);
    DLM_CleverMcLevyResoTM* GaussCoreRsm_pK(const int& SourceVar);
    DLM_CleverMcLevyResoTM* GaussCoreRsm_pp(const int& SourceVar);

private:
    void Clean_CommonAnaFunctions();
    //MS_GaussExp_mT_Simple* Simple_Reso;
    DLM_CleverLevy* CleverLevy;
    DLM_CleverMcLevyReso* CleverMcLevyReso;
    DLM_CleverMcLevyResoTM* CleverMcLevyResoTM;
    const unsigned NumCleverLevyObjects;
    TString* CatsFilesFolder;
};

DLM_Histo<double>* ConvertThetaAngleHisto(const TString& FileName, const TString& HistoName, const double kMin, const double kMax);

void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot);
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, CATS* Kitty);
void RootFile_DlmSource(const TString& RootFileName, const TString& GraphName, CATS* Kitty, const unsigned& NumBins, const double& rMin, const double& rMax, const double& lambda=1, const bool& FourPi=true);
double Get_reff(TH1F* hsource, const float lambda=1, const float CEI=0.9);
double Get_reff_TF1(TH1F* hsource, TF1* fsource, const float lambda=1, const float CEI=0.9);
double GetRcore(DLM_CleverMcLevyResoTM& MagicSource, const double& reff);
double GetReff(DLM_CleverMcLevyResoTM& MagicSource, const double& rcore);

//quick way to set up the pp or pL source with default settings (default meaning RUN2 paper)
//flag X0 is the RSM with the correct sign, flag X1 is the RSM with wrong sign (original)
//flag 0X is the EPOS angulars, flag 1X is flat angles
void BasicSetUp_MS(DLM_CleverMcLevyResoTM& MagicSource, double frac1, double frac2);
void SetUp_RSMflat_pp(DLM_CleverMcLevyResoTM& MagicSource);
void SetUp_RSM_pp(DLM_CleverMcLevyResoTM& MagicSource, const TString InputFolder, const int flag=0);
void SetUp_RSMflat_pL(DLM_CleverMcLevyResoTM& MagicSource);
void SetUp_RSM_pL(DLM_CleverMcLevyResoTM& MagicSource, const TString InputFolder, const int flag=0);
void SetUp_RSMflat_pXi(DLM_CleverMcLevyResoTM& MagicSource);
void SetUp_RSM_pXi(DLM_CleverMcLevyResoTM& MagicSource, const TString InputFolder, const int flag=0);
void SetUp_RSMflat_pOmega(DLM_CleverMcLevyResoTM& MagicSource);
void SetUp_RSM_pOmega(DLM_CleverMcLevyResoTM& MagicSource, const TString InputFolder, const int flag=0);
/*
class DLM_Analyzer{

public:

    DLM_Analyzer();
    ~DLM_Analyzer();


private:
    TH1F* hData;
    DLM_Fitter1* fitter;
    TString System;

};
*/

#endif
