
#include "Worksheet.h"
#include "ExtendedCk.h"
#include "Basics.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "DLM_Ck.h"
#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"

CATS* KITTY_CATS_WORKFIT_PL;

double CATS_WORKFIT_PL(double* x, double* pars){
    double& Norm = pars[0];
    double& LambdaPar = pars[1];
    double& SourceSize = pars[2];
    //set the radius to the fit value
    //the last parameter says its a small change of the radius, which does not require a new computing grid (saves time)
    //however, the last step requires making sure a good initial value of the radius
    KITTY_CATS_WORKFIT_PL->SetAnaSource(0,SourceSize,true);
    //useful tip: this makes CATS to shut up and not flood your screen. Only errors will be displayed. Use nSilent to suppress even those
    KITTY_CATS_WORKFIT_PL->SetNotifications(CATS::nError);
    KITTY_CATS_WORKFIT_PL->KillTheCat();
    return Norm*(LambdaPar*KITTY_CATS_WORKFIT_PL->EvalCorrFun(*x)+1.-LambdaPar);
}


void Worksheet_ProtonLambda(){

    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    CATSparameters POT_PARS_1S0(CATSparameters::tPotential,4,true);
    POT_PARS_1S0.SetParameter(0,0);
    POT_PARS_1S0.SetParameter(1,2137);
    POT_PARS_1S0.SetParameter(2,0.5);
    POT_PARS_1S0.SetParameter(3,0.2);
    CATSparameters POT_PARS_3S1(CATSparameters::tPotential,4,true);
    POT_PARS_3S1.SetParameter(0,1);
    POT_PARS_3S1.SetParameter(1,2137);
    POT_PARS_3S1.SetParameter(2,0.5);
    POT_PARS_3S1.SetParameter(3,0.2);

    CATS Kitty_pL;

    //DEFINE YOUR CATS OBJECT HERE
    //USE the potential function:
    Kitty_pL.SetShortRangePotential(0,0,Basics_Potential_Usmani,POT_PARS_1S0);
    Kitty_pL.SetShortRangePotential(1,0,Basics_Potential_Usmani,POT_PARS_3S1);

    Kitty_pL.KillTheCat();

    KITTY_CATS_WORKFIT_PL = &Kitty_pL;

    TFile* fInput = new TFile("../Files/DummyProtonLambda.root","read");
    TH1F* hInput = (TH1F*)fInput->Get("hDummyProtonLambda");
    TF1* FITTER = new TF1("FITTER",CATS_WORKFIT_PL,kMin,kMax,3);
    FITTER->SetParameter(0,1);FITTER->SetParLimits(0,0.5,2);
    FITTER->SetParameter(1,0.5);FITTER->SetParLimits(1,0,1);
    //! SET THE SOURCE PAR
    //FITTER->SetParameter(2,SOURCE_PARS.GetParameter(0));FITTER->SetParLimits(2,0.5,3);


    //FITTER->SetParameter(3,POT_PARS_1S0.GetParameter(1));FITTER->SetParLimits(3,0.99*POT_PARS_1S0.GetParameter(1),1.01*POT_PARS_1S0.GetParameter(1));
    //FITTER->SetParameter(4,POT_PARS_1S0.GetParameter(2));FITTER->SetParLimits(4,0.99*POT_PARS_1S0.GetParameter(2),1.01*POT_PARS_1S0.GetParameter(2));
    //FITTER->SetParameter(5,POT_PARS_1S0.GetParameter(3));FITTER->SetParLimits(5,0.99*POT_PARS_1S0.GetParameter(3),1.01*POT_PARS_1S0.GetParameter(3));

    //FITTER->FixParameter(1,0.6);
    //FITTER->FixParameter(3,POT_PARS_1S0.GetParameter(1));
    //FITTER->FixParameter(4,POT_PARS_1S0.GetParameter(2));
    //FITTER->FixParameter(5,POT_PARS_1S0.GetParameter(3));


    hInput->Fit(FITTER,"S, N, R, M");

    delete FITTER;
    delete fInput;

}


//
DLM_Ck* Worksheet_SetUp_pL(){

    const unsigned NumMomBins = 80;
    const double kMin = 0;
    const double kMax = 320;

    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    SOURCE_PARS.SetParameter(0,1.3);

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
    double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1(CATSparameters::tPotential,8,true);
    cPotPars3S1.SetParameters(PotPars3S1);

    static CATS Kitty_pL;
    Kitty_pL.SetMomBins(NumMomBins,kMin,kMax);
    Kitty_pL.SetAnaSource(GaussSource,SOURCE_PARS);
    Kitty_pL.SetAnaSource(0,SOURCE_PARS.GetParameter(0));
    Kitty_pL.SetUseAnalyticSource(true);
    Kitty_pL.SetMomentumDependentSource(false);
    Kitty_pL.SetThetaDependentSource(false);
    Kitty_pL.SetExcludeFailedBins(false);
    Kitty_pL.SetQ1Q2(0);
    Kitty_pL.SetPdgId(2212, 3122);
    Kitty_pL.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
    Kitty_pL.SetNumChannels(2);
    Kitty_pL.SetNumPW(0,1);
    Kitty_pL.SetNumPW(1,1);
    Kitty_pL.SetSpin(0,0);
    Kitty_pL.SetSpin(1,1);
    Kitty_pL.SetChannelWeight(0, 1./4.);
    Kitty_pL.SetChannelWeight(1, 3./4.);
    Kitty_pL.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    Kitty_pL.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1);

    Kitty_pL.KillTheCat();

    DLM_Ck* Ck_pL = new DLM_Ck(1,0,Kitty_pL);
    Ck_pL->Update();

    RootFile_DlmCk("Worksheet_DlmCk.root","Ck_pL",Ck_pL);

    return Ck_pL;
}

//potential: NN_AV18
DLM_Ck* Worksheet_SetUp_pp(){

    const unsigned NumMomBins = 80;
    const double kMin = 0;
    const double kMax = 320;

    CATSparameters SOURCE_PARS(CATSparameters::tSource,1,true);
    SOURCE_PARS.SetParameter(0,1.3);

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
    double PotPars1D2[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,2,2};

    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true); cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3P0(CATSparameters::tPotential,8,true); cPotPars3P0.SetParameters(PotPars3P0);
    CATSparameters cPotPars3P1(CATSparameters::tPotential,8,true); cPotPars3P1.SetParameters(PotPars3P1);
    CATSparameters cPotPars3P2(CATSparameters::tPotential,8,true); cPotPars3P2.SetParameters(PotPars3P2);
    CATSparameters cPotPars1D2(CATSparameters::tPotential,8,true); cPotPars1D2.SetParameters(PotPars1D2);


}

//potential: pXim_HALQCD1
DLM_Ck* Worksheet_SetUp_pXi(){
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotParsI0S0[9]={pXim_HALQCD1,12,0,-1,1,0,0,0,0};
    double PotParsI0S1[9]={pXim_HALQCD1,12,0,-1,1,1,0,1,0};
    double PotParsI1S0[9]={pXim_HALQCD1,12,1,1,1,0,0,0,0};
    double PotParsI1S1[9]={pXim_HALQCD1,12,1,1,1,1,0,1,0};




}

//Ck model: Lednicky_gauss_Sigma0 (only 1 source parameter, no interaction pars)
DLM_Ck* Worksheet_SetUp_pSigma(){

}
