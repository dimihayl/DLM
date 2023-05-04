
#include "Basics.h"
#include "ExtendedCk.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"

#include "DLM_Potentials.h"
#include "DLM_Source.h"
#include "DLM_CkDecomposition.h"
#include "DLM_CkModels.h"
#include "DLM_Ck.h"
#include "DLM_Random.h"
#include "DLM_CppTools.h"

#include "TGraph.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TROOT.h"

///An example how to compute the pLambda correlation function using the Usmani potential and the Lednicky model
void Ck_pL_Ledni_Usmani(){
    const unsigned NumMomBins = 60;
    const double kMin = 0;
    const double kMax = 240;

    //object for the parameters to be used by the source function
    CATSparameters SOURCE_PARS(CATSparameters::tSource,2,true);
    //set the source size to 1.2 fm
    SOURCE_PARS.SetParameter(0,1.2);

    //the CATS object to model pLambda (using the Usmani potential)
    CATS Kitty_pL;
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
    //two spin channels (0 and 1)
    Kitty_pL.SetNumChannels(2);
    Kitty_pL.SetNumPW(0,1);
    Kitty_pL.SetNumPW(1,1);
    Kitty_pL.SetSpin(0,0);
    Kitty_pL.SetSpin(1,1);
    //the weights are based on the singlet/triplet configuration
    Kitty_pL.SetChannelWeight(0, 1./4.);
    Kitty_pL.SetChannelWeight(1, 3./4.);

    //the standard input to the predefined potentials using fDlmPot
    //this is the configuration for Usmani, where apart from the flag,
    //the the only relevant parameter is the spin (s)
    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
    double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,0};
    //define a set of parameters for each channel
    CATSparameters cPotPars1S0(CATSparameters::tPotential,8,true);
    cPotPars1S0.SetParameters(PotPars1S0);
    CATSparameters cPotPars3S1(CATSparameters::tPotential,8,true);
    cPotPars3S1.SetParameters(PotPars3S1);

    //give to CATS the relevant information for the interaction in each channel
    Kitty_pL.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    Kitty_pL.SetShortRangePotential(1,0,fDlmPot,cPotPars3S1);
    //Kitty_pL.SetMaxNumThreads(4);
    Kitty_pL.KillTheCat();

    //The DLM_Ck object is in essence a histogram, that can be used to read a CATS object, or some function
    //Here we define the DLM_Ck for the Usmani potential, from the CATS object
    //number of source/potential pars to be "accessed"
    DLM_Ck Ck_Usmani(1,0,Kitty_pL);
    //Ck_Usmani.SetSourcePar(0,2);
    Ck_Usmani.Update();

    //here we define the DLM_Ck object for the Lednicky interaction. The arguments are:
    //#source parameters, #interaction parameters, k min, k max, pointer to the function to be used
    DLM_Ck Ck_Lednicky(1,4,NumMomBins,kMin,kMax,Lednicky_SingletTriplet);
    //set the source and interaction (potential) parameters.
    //for the interaction we have used the scattering lengths (pars 0 and 2) and the effective ranges (1 and 3)
    //for the Usmani potential.
    Ck_Lednicky.SetSourcePar(0,SOURCE_PARS.GetParameter(0));
    Ck_Lednicky.SetPotPar(0,2.88);
    Ck_Lednicky.SetPotPar(1,2.92);
    Ck_Lednicky.SetPotPar(2,1.66);
    Ck_Lednicky.SetPotPar(3,3.78);
    Ck_Lednicky.Update();

    RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_Usmani.root","gCk_pL_Usmani",&Ck_Usmani);
    RootFile_DlmCk("./OutputFiles/Ck_pL_Ledni_Usmani.root","gCk_pL_Lednicky",&Ck_Lednicky);
}
//save a DLM_Ck to a file in the form of a TGraph
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot){
    TFile* RootFile = new TFile(RootFileName,"update");
    if(!RootFile) RootFile = new TFile(RootFileName,"recreate");
    const unsigned NumBins = CkToPlot->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        graph.SetPoint(uBin,CkToPlot->GetBinCenter(0,uBin),CkToPlot->GetBinContent(uBin));
    }
    graph.Write("",TObject::kOverwrite);
    delete RootFile;
}
//save a DLM_CkDecomposition to a file in the form of a TGraph
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_CkDecomposition* CkToPlot, bool PlotIndividualContributions){
    TFile* RootFile = new TFile(RootFileName,"update");
    if(!RootFile) RootFile = new TFile(RootFileName,"recreate");
    const unsigned NumBins = CkToPlot->GetCk()->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        graph.SetPoint(uBin,CkToPlot->GetCk()->GetBinCenter(0,uBin),CkToPlot->EvalCk(CkToPlot->GetCk()->GetBinCenter(0,uBin)));
    }
    graph.Write("",TObject::kOverwrite);

    //this bit is to get the individual contributions (associated to lambda pars) to the total correlation function
    if(PlotIndividualContributions){
        TGraph graphSignal;
        graphSignal.SetName(GraphName+"_signal");
        graphSignal.Set(NumBins);
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0,uBin);
            //the signal is defined as C(k)-1
            graphSignal.SetPoint(uBin,MOMENTUM,CkToPlot->EvalSignalSmeared(MOMENTUM));
        }
        graphSignal.Write("",TObject::kOverwrite);

        TGraph graphMain;
        graphMain.SetName(GraphName+"_main");
        graphMain.Set(NumBins);
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0,uBin);
            //this is the contribution of the main channel only
            graphMain.SetPoint(uBin,MOMENTUM,CkToPlot->EvalSignalSmearedMain(MOMENTUM));
        }
        graphMain.Write("",TObject::kOverwrite);

        for(unsigned uChild=0; uChild<CkToPlot->GetNumChildren(); uChild++){
            TGraph graphChild;
            graphChild.SetName(GraphName+"_child"+uChild);
            graphChild.Set(NumBins);
            for(unsigned uBin=0; uBin<NumBins; uBin++){
                double MOMENTUM = CkToPlot->GetCk()->GetBinCenter(0,uBin);
                //this is the contribution of each secondary channel, e.g. feed-down or misid
                graphChild.SetPoint(uBin,MOMENTUM,CkToPlot->EvalSignalSmearedChild(uChild,MOMENTUM));
            }
            graphChild.Write("",TObject::kOverwrite);
        }
    }
    delete RootFile;
}

//initialize the Gauss source for the CATS object
void CATS_GaussSource(CATS& Kitty, const double& SourceSize){
    //object containing the source or potential parameters. Arguments:
    //(tSource or tPotential, NumberOfParameters, ReadyForMulti-threading)
    CATSparameters cPars(CATSparameters::tSource,1,true);
    //set up the parameters
    cPars.SetParameter(0,SourceSize);
    //set the source to a Gaussian function, using the parameters above. These parameters are copied into CATS,
    //thus it is okay if cPars is deleted afterwards.
    Kitty.SetAnaSource(GaussSource, cPars);
    //this is important, since if its `false` it is assumed that the source will be sampled from a transport model, and the Gaussian function will not be used
    Kitty.SetUseAnalyticSource(true);
    //if true, the source is automatically renormalized in the range 0-64 fm. Nice to dummy proof the source, but problematic for sources with large tails
    //for the Gaussian example above, both options should be completely identical
    Kitty.SetAutoNormSource(false);
}
//initialize the Gauss core + resonances for the CATS object
//note that the core is more generic than that, and is in the form of a Levy-Stable distributions,
//which is a class of distributions to which the Cauchy (Alpha=1) and Gauss (Alpha=2) are special cases.
//fracreso gives the fraction of resonance feeding into your protons, taureso is the avg lifetime of the resonances, massreso is the avg mass
//the cutoff sets the limit for the maximum allowed relative momentum of the daughters (ideally it should be a small number, in reality c.a. 200 MeV works)
//the returned pointer is to the source object, that must not be deleted while the CATS object is used. You should delete it after that, to avoid memory leaks
DLM_CleverMcLevyResoTM* CATS_ResoSource_pp(CATS& Kitty, const double& SourceSize, const double& Alpha, const bool& FixAlpha,
                        const double& fracreso, const double& taureso, const double& massreso, const double& cutoff){

    //In the DLM_Source.h there is a special class created to model this scenario, called DLM_CleverMcLevyResoTM
    //Note, that if you want to create some fancy source class as well, you can do this by following the skeleton
    //provided by the base class 'CatsSource', from which you should inherit (see the implementation of DLM_CleverMcLevyResoTM as an example)

    //create the source object
    //this object must NOT be deleted while the program is running
    DLM_CleverMcLevyResoTM* CleverMcLevyResoTM_pp = new DLM_CleverMcLevyResoTM();;

    //standard initialization of DLM_CleverMcLevyResoTM in the following:
    //sets up a discrete grid size, on which only a finite amount of source evaluations will be performed to save CPU time
    //the source function will be extrapolated to the desired parameters.
    if(FixAlpha) CleverMcLevyResoTM_pp->InitStability(1,Alpha-1e-6,Alpha+1e-6);//only a single grid point for the Alpha
    else CleverMcLevyResoTM_pp->InitStability(21,1,2);
    //this sets the possible source sizes between 0.4 and 3 fm. Typically should provide good cover of all ranges, but take care you do not go outside
    CleverMcLevyResoTM_pp->InitScale(52,0.4,3.0);
    CleverMcLevyResoTM_pp->InitRad(400,0,64);
    //do not change this ever, no questions asked
    CleverMcLevyResoTM_pp->InitType(2);
    //for different particle species, e.g. pL, you will need to setup the 0th and 1st particle differently
    CleverMcLevyResoTM_pp->SetUpReso(0,fracreso);
    CleverMcLevyResoTM_pp->SetUpReso(1,fracreso);

    Float_t k_D,fP1,fP2,fM1,fM2,Tau1,Tau2,AngleRcP1,AngleRcP2,AngleP1P2;
    DLM_Random RanGen(11);
    double RanVal1,RanVal2;

    //this is a TNtuple containing the relevant information from your transport model, which is to be used to model the kinematics.
    //The variables are those used in eq.19 of the analysis note here: https://alice-notes.web.cern.ch/node/891
    //if you do not have access there, write me an email
    //N.B. ALL these observables should be evaluated in the CM frame of the two final particle you study!
    //To set up such an TNtuple: be creative :) In case of difficulties (or if you need EPOS output) write me or Max an email and we can discuss.
    //we actually have two TNtuple, one for the case of primary-secondary protons, and one for secondary-secondary. Note that for non-identical particles
    //you will also have to consider separately secondary-primary.
    TFile* F_EposDisto_p_pReso = new TFile("../Files/EposDisto_p_pReso.root");
    TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
    T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    //iterate over the TNtuple to give the source the needed input
    //this is where the magic happens. The idea is, that we give the source some amount of vectors corresponding to s_res1/2 in fig. 20 of the analysis note
    //the source will than build the source by random sampling the core and putting on top some of these vectors,
    //to evaluate the resulting separation between the daughters (final particles of interest). We can either use the information fully from our TNtuple, i.e.
    //take the information about Tau, Mass, Momentum and angles from our transfer model, or use some averaged values. The latter is used since in EPOS we often do not
    //have all resonance we want, so we used some dummy resonances just to obtain a sample with the correct average mass and kinematic properties, although the lifetimes
    //could differ, and the individual masses wrong.
    for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        //the values were provided from the statistical hadronization model, we rewrite the EPOS output
        Tau2 = taureso;
        fM2 = massreso;
        //reject daughters with too high relative momenta
        if(k_D>cutoff) continue;
        //sample randomly the fly length of the resonance
        RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
        //plug in the result for the primordial-resonance case
        CleverMcLevyResoTM_pp->AddBGT_PR(RanVal1,-cos(AngleRcP2));
        //plug in the result for the primordial-resonance case. For identical particle it is symmetric, apart the sign of the cosine.vid
        CleverMcLevyResoTM_pp->AddBGT_RP(RanVal1,cos(AngleRcP2));
    }
    delete F_EposDisto_p_pReso;

    //do exactly the same thing for the reso-reso case
    TFile* F_EposDisto_pReso_pReso = new TFile("../Files/EposDisto_pReso_pReso.root");
    TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
    for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        Tau1 = taureso;
        Tau2 = taureso;
        fM1 = massreso;
        fM2 = massreso;
        if(k_D>cutoff) continue;
        RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
        RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
        CleverMcLevyResoTM_pp->AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_pReso;

    //the Gauss and Cauchy cases are evaluated much faster, so we can allow ourselves more iterations
    //the number of iter. here is how many entries will the final pdf of the source function contain
    if(FixAlpha&&(Alpha==1||Alpha==2)) CleverMcLevyResoTM_pp->InitNumMcIter(1000000);
    else CleverMcLevyResoTM_pp->InitNumMcIter(100000);

    //The CatsSourceForwarder is a technicality that you should not care about,
    //here it is important to know that the source function is passed into cats by a pointer to the object
    //CleverMcLevyResoTM, and the number of parameters of the source.
    Kitty.SetAnaSource(CatsSourceForwarder, CleverMcLevyResoTM_pp, 2);
    //and this, as aways, is the way to communicate with CATS and set the source size
    Kitty.SetAnaSource(0,SourceSize);
    Kitty.SetAnaSource(1,Alpha);

    //this is important, since if its `false` it is assumed that the source will be sampled from a transport model, and the Gaussian function will not be used
    Kitty.SetUseAnalyticSource(true);
    //if true, the source is automatically renormalized in the range 0-64 fm. Nice to dummy proof the source, but problematic for sources with large tails
    //For the core+reso case we could easily have large tails, so this is quite important to set as false
    Kitty.SetAutoNormSource(false);

    return CleverMcLevyResoTM_pp;
}


//Basic initialization of a CATS object for pp, WITHOUT any setup of the source or interaction
void CATS_pp_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetExcludeFailedBins(false);
    //the charge of the pair given in q1 * q2, where q1 and q2 are the charges of the individual particles
    Kitty.SetQ1Q2(1);
    //pid of the particle. Important only to set properly the quantum statistics (on for identical guys)
    //one could also force the QS by: Kitty.SetQuantumStatistics(true/false);
    Kitty.SetPdgId(2212, 2212);
    //reduced mass of the pair
    Kitty.SetRedMass( 0.5*Mass_p );
}
//initialize the interaction for pp, using the AV18 potential.
//by default we have only s-waves included, but one can include the p and d waves if desired
void CATS_pp_AV18(CATS& Kitty, const bool& pwaves, const bool& dwaves){
    //the 4 channels for pp are:
    //s=0: 1S0 + 3D1
    //s=1: 3P0
    //s=1: 3P1
    //s=1: 3P2
    //note that for s=0 the p-waves are Pauli blocked, for s=1 these are the s and d waves
    if(pwaves){
        Kitty.SetNumChannels(4);
        if(dwaves) Kitty.SetNumPW(0,3);
        else Kitty.SetNumPW(0,1);
        Kitty.SetNumPW(1,2);
        Kitty.SetNumPW(2,2);
        Kitty.SetNumPW(3,2);
        Kitty.SetSpin(0,0);
        Kitty.SetSpin(1,1);
        Kitty.SetSpin(2,1);
        Kitty.SetSpin(3,1);
        Kitty.SetChannelWeight(0, 3./12.);
        Kitty.SetChannelWeight(1, 1./12.);
        Kitty.SetChannelWeight(2, 3./12.);
        Kitty.SetChannelWeight(3, 5./12.);
    }
    else{
        //important: even with the p-waves switched off, physics wise the spin 1 state still exists and
        //the p-waves are there, just in there asymptotic state (free wave). To include this in the computation,
        //CATS still needs a second channel, even if it is `empty`!
        Kitty.SetNumChannels(2);
        if(dwaves) Kitty.SetNumPW(0,3);
        else Kitty.SetNumPW(0,1);
        Kitty.SetNumPW(1,0);
        Kitty.SetSpin(0,0);
        Kitty.SetSpin(1,1);
        Kitty.SetChannelWeight(0, 1./4.);
        Kitty.SetChannelWeight(1, 3./4.);
    }

    //to set up the strong interaction, one can use the predefined functions available in DLM_Potentials.h
    //the main idea is to always pass the function fDlmPot but with different input parameters, based on which the interaction is set up automatically
    //this works not only for pp, but many other systems are included.
    //The fDlmPot is only used as an interface to easily use different potentials, that are hard coded in DLM_Potentials.h
    //Feel free to expand the data base of this file if you think others will benefit from it. For further details contact Dimi
    //The input parameters are by default 9 and are defined as follows:
    //0: potential flag (defines which potential to use, see the enumerators in DLM_Potentials.h for more info)
    //1: a second flag, that can be used if needed (depending on the definition of the potential)
    //2: total isospin
    //3: 2 x isospin of particle 1
    //4: 2 x isospin of particle 2
    //5: total spin
    //6: l quantum number
    //7: j quantum number
    CATSparameters cPars_pp_1S0(CATSparameters::tPotential,8,true);
    cPars_pp_1S0.SetParameter(0,NN_ReidV8);//choose the AV18
    cPars_pp_1S0.SetParameter(1,v18_Coupled3P2);//default option, which takes the 3P2 channel from a coupled-channel computation, but in CATS only the first diagonal potential elements is used
    cPars_pp_1S0.SetParameter(2,1);
    cPars_pp_1S0.SetParameter(3,1);
    cPars_pp_1S0.SetParameter(4,1);
    cPars_pp_1S0.SetParameter(5,0);
    cPars_pp_1S0.SetParameter(6,0);
    cPars_pp_1S0.SetParameter(7,0);

    //copy all settings from cPars_pp_1S0, and just change quantum numbers s,l,j
    CATSparameters cPars_pp_3P0(cPars_pp_1S0);
    cPars_pp_3P0.SetParameter(5,1);
    cPars_pp_3P0.SetParameter(6,1);
    cPars_pp_3P0.SetParameter(7,0);

    CATSparameters cPars_pp_3P1(cPars_pp_1S0);
    cPars_pp_3P1.SetParameter(5,1);
    cPars_pp_3P1.SetParameter(6,1);
    cPars_pp_3P1.SetParameter(7,1);

    CATSparameters cPars_pp_3P2(cPars_pp_1S0);
    cPars_pp_3P2.SetParameter(5,1);
    cPars_pp_3P2.SetParameter(6,1);
    cPars_pp_3P2.SetParameter(7,2);

    CATSparameters cPars_pp_1D2(cPars_pp_1S0);
    cPars_pp_1D2.SetParameter(5,0);
    cPars_pp_1D2.SetParameter(6,2);
    cPars_pp_1D2.SetParameter(7,2);

    //plug in the strong potential for each channel and partial wave
    //the arguments are: #WhichChannel,#WhichPartialWave,#PotentialFunction,#PotentialParameters
    Kitty.SetShortRangePotential(0,0,fDlmPot,cPars_pp_1S0);
    if(pwaves){
        Kitty.SetShortRangePotential(1,1,fDlmPot,cPars_pp_3P0);
        Kitty.SetShortRangePotential(2,1,fDlmPot,cPars_pp_3P1);
        Kitty.SetShortRangePotential(3,1,fDlmPot,cPars_pp_3P2);
    }
    if(dwaves){
        Kitty.SetShortRangePotential(0,2,fDlmPot,cPars_pp_1D2);
    }
    //if later on you would like to switch some contribution off, this can be done with:
    //Kitty.RemoveShortRangePotential(#WhichChannel,#WhichPartialWave);
}
void CATS_pL_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );
}
void CATS_pL_Usmani(CATS& Kitty){
    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetChannelWeight(0, 1./4.);
    Kitty.SetChannelWeight(1, 3./4.);

    CATSparameters cPars_pL_Usmani1S0(CATSparameters::tPotential,8,true);
    cPars_pL_Usmani1S0.SetParameter(0,pL_UsmaniOli);
    cPars_pL_Usmani1S0.SetParameter(5,0);
    cPars_pL_Usmani1S0.SetParameter(7,0);

    CATSparameters cPars_pL_Usmani3S1(CATSparameters::tPotential,8,true);
    cPars_pL_Usmani3S1.SetParameter(0,pL_UsmaniOli);
    cPars_pL_Usmani3S1.SetParameter(5,1);
    cPars_pL_Usmani3S1.SetParameter(7,1);

    Kitty.SetShortRangePotential(0,0,fDlmPot,cPars_pL_Usmani1S0);
    Kitty.SetShortRangePotential(1,0,fDlmPot,cPars_pL_Usmani3S1);
}
void CATS_pXim_Basic(CATS& Kitty, const unsigned& NumMomBins, const double& kMin, const double& kMax){
    Kitty.SetMomBins(NumMomBins,kMin,kMax);
    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );
}
void CATS_pXim_Hal(CATS& Kitty){
    //0: S=0 I=0
    //1: S=1 I=0
    //2: S=0 I=1
    //3: S=1 I=1
    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetNumPW(2,1);
    Kitty.SetNumPW(3,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetSpin(2,0);
    Kitty.SetSpin(3,1);
    Kitty.SetChannelWeight(0, 1./8.);
    Kitty.SetChannelWeight(1, 3./8.);
    Kitty.SetChannelWeight(2, 1./8.);
    Kitty.SetChannelWeight(3, 3./8.);

    CATSparameters cPars_pXim_HalI0S0(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI0S0.SetParameter(0,pXim_HALQCD1);
    cPars_pXim_HalI0S0.SetParameter(1,12);
    cPars_pXim_HalI0S0.SetParameter(2,0);
    cPars_pXim_HalI0S0.SetParameter(3,-1);
    cPars_pXim_HalI0S0.SetParameter(4,1);
    cPars_pXim_HalI0S0.SetParameter(5,0);
    cPars_pXim_HalI0S0.SetParameter(6,0);
    cPars_pXim_HalI0S0.SetParameter(7,0);

    CATSparameters cPars_pXim_HalI0S1(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI0S1.SetParameter(0,pXim_HALQCD1);
    cPars_pXim_HalI0S1.SetParameter(1,12);
    cPars_pXim_HalI0S1.SetParameter(2,0);
    cPars_pXim_HalI0S1.SetParameter(3,-1);
    cPars_pXim_HalI0S1.SetParameter(4,1);
    cPars_pXim_HalI0S1.SetParameter(5,1);
    cPars_pXim_HalI0S1.SetParameter(6,0);
    cPars_pXim_HalI0S1.SetParameter(7,1);

    CATSparameters cPars_pXim_HalI1S0(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S0.SetParameter(0,pXim_HALQCD1);
    cPars_pXim_HalI1S0.SetParameter(1,12);
    cPars_pXim_HalI1S0.SetParameter(2,1);
    cPars_pXim_HalI1S0.SetParameter(3,1);
    cPars_pXim_HalI1S0.SetParameter(4,1);
    cPars_pXim_HalI1S0.SetParameter(5,0);
    cPars_pXim_HalI1S0.SetParameter(6,0);
    cPars_pXim_HalI1S0.SetParameter(7,0);

    CATSparameters cPars_pXim_HalI1S1(CATSparameters::tPotential,8,true);
    cPars_pXim_HalI1S1.SetParameter(0,pXim_HALQCD1);
    cPars_pXim_HalI1S1.SetParameter(1,12);
    cPars_pXim_HalI1S1.SetParameter(2,1);
    cPars_pXim_HalI1S1.SetParameter(3,1);
    cPars_pXim_HalI1S1.SetParameter(4,1);
    cPars_pXim_HalI1S1.SetParameter(5,1);
    cPars_pXim_HalI1S1.SetParameter(6,0);
    cPars_pXim_HalI1S1.SetParameter(7,1);

    Kitty.SetShortRangePotential(0,0,fDlmPot,cPars_pXim_HalI0S0);
    Kitty.SetShortRangePotential(1,0,fDlmPot,cPars_pXim_HalI0S1);
    Kitty.SetShortRangePotential(2,0,fDlmPot,cPars_pXim_HalI1S0);
    Kitty.SetShortRangePotential(3,0,fDlmPot,cPars_pXim_HalI1S1);
}


TH2F* GetSmearMatrix(TString filename, TString histoname){
    TFile* FileROOT = new TFile(filename, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}
TH1F* GetExpCorrelation(TString filename, TString histoname){
    TFile* FileROOT = new TFile(filename, "read");
    if(!FileROOT){printf("\033[1;31mERROR:\033[0m The file '%s' does not exist\n",filename.Data());return NULL;}
    TH1F* histo = (TH1F*)FileROOT->Get(histoname);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",histoname.Data(),filename.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT; FileROOT=NULL;
    histoCopy->SetName(Name);
    return histoCopy;
}

DLM_CkDecomposition* Ck_ExampleFitter_pp;
//for this example fitter we have the source size as a free fit parameter [0],
//as long as a pol1 baseline with parameters [1] and [2]
double ExampleFitter_pp(double* x, double* par){
    double& MOM = *x;
    //N.B. this only changes the radius of the p-p, not the feed-down! If the latter is desired
    //it can be additionally set by calling the relevant contributions, e.g. to set the radius of pL->pp one should use:
    //Ck_ExampleFitter_pp->GetContribution("pLambda")->GetCk()->SetSourcePar(0,par[0]);
    //N.B. (2): the SetSourcePar functions assumes we do not significantly change the radius, and when passing this argument into CATS
    //it forces CATS to use the computing grid that is already set (saves a LOT of time). But will only work for source size variations of up to 20-30% (50% if you push it)
    //before running into possible numerical bias. So make sure your starting fit values are reasonable
    Ck_ExampleFitter_pp->GetCk()->SetSourcePar(0,par[0]);
    Ck_ExampleFitter_pp->Update();

    double Correlation = Ck_ExampleFitter_pp->EvalCk(MOM);
    double Baseline = par[1]*(1.+par[2]*MOM);
    return Baseline*Correlation;
}

///an example how compute the pp theoretical correlation
//for the source type, we have either 'Gauss' or 'CoreReso'
void Ck_pp_Decomposition(const TString& SourceType){
    //timer
    DLM_Timer TIMER;
    printf("\033[1;37mExecuting Ck_pp_Decomposition for SourceType=='%s'\033[0m\n\n",SourceType.Data());
    if(SourceType!="Gauss"&&SourceType!="CoreReso"){printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SourceType.Data()); return;}
    TString OutputFileName = "../OutputFiles/Ck_pp_Decomposition_"+SourceType+".root";
    //set the source size slightly smaller for the Gauss+Reso scenario.
    const double SourceSize = SourceType=="Gauss"?1.2:1.1;
    //the source sizes for the feed downs. There we aways use a Gaussian source (to save CPU time)
    //these values are more or less in tone with our HM analysis
    const double SourceSize_pL = 1.3;
    const double SourceSize_pXim = 1.0;
    //these are for the CATS object
    const double NumMomBins = 100;
    const double kMin = 0;
    const double kMax = 400;
    //these are for the fit. Edge effects could occur when the smearing is applied,
    //thus recommended to always keep the fit range below the range of the cats object.
    //typically I leave a buffer of 40-50 MeV
    const double FitMin = 0;
    const double FitMax = 350;

    CATS Kitty_pp_s;
    //initialize a Gaussian source
    if(SourceType=="Gauss") CATS_GaussSource(Kitty_pp_s,SourceSize);
    //initialize a Gaussian core + reso source
    DLM_CleverMcLevyResoTM* CRS_pp = NULL;
    if(SourceType=="CoreReso") CRS_pp = CATS_ResoSource_pp(Kitty_pp_s,SourceSize,2,true,0.6422,1.65,1362,200);

    //you can change the parameters of the source at any time using the function:
    //Kitty_pp_s.SetAnaSource(#WhichParameter,#Value);
    //set up the CATS object for pp (Coulomb and QS) with 100 bins in the range 0-400 MeV
    CATS_pp_Basic(Kitty_pp_s,NumMomBins,kMin,kMax);
    //set up the cats interaction including only s-waves
    CATS_pp_AV18(Kitty_pp_s,false,false);
    //compute the correlation function
    printf(" Executing KillTheCat for Kitty_pp_s\n");
    Kitty_pp_s.KillTheCat();
    //set up a `histogram` for the pp interaction, based on the above CATS object
    //the arguments are the number of source/potential parameters to be controlled by DLM_Ck
    DLM_Ck Ck_pp_s(Kitty_pp_s.GetNumSourcePars(),0,Kitty_pp_s);
    //btw, now you can change the source parameters also by using:
    //Ck_pp_s.SetSourcePar(#WhichParameter,#Value);
    //this function reads the CATS objects and fills the bins of the DLM_Ck object
    Ck_pp_s.Update();
    //save the correlation (theoretical at the moment) in a file
    RootFile_DlmCk(OutputFileName, "Ck_pp_sWaves", &Ck_pp_s);

    //the same for s and p-waves
    CATS Kitty_pp_sp;
    if(SourceType=="Gauss") CATS_GaussSource(Kitty_pp_sp,SourceSize);
    DLM_CleverMcLevyResoTM* CRS_pp_sp = NULL;
    if(SourceType=="CoreReso") CRS_pp_sp = CATS_ResoSource_pp(Kitty_pp_sp,SourceSize,2,true,0.6422,1.65,1362,200);
    CATS_pp_Basic(Kitty_pp_sp,100,0,400);
    CATS_pp_AV18(Kitty_pp_sp,true,false);
    printf(" Executing KillTheCat for Kitty_pp_sp\n");
    //this function can be used to reduce the CATS output if annoying to you.
    //CATS::nWarning only prints Warnings and Errors
    //CATS::nError only prints Errors
    //CATS::nError silent completely removes all messages
    Kitty_pp_sp.SetNotifications(CATS::nWarning);
    Kitty_pp_sp.KillTheCat();
    DLM_Ck Ck_pp_sp(Kitty_pp_sp.GetNumSourcePars(),0,Kitty_pp_sp);
    Ck_pp_sp.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pp_spWaves", &Ck_pp_sp);

    //and now the same including d-waves
    CATS Kitty_pp_spd;
    if(SourceType=="Gauss") CATS_GaussSource(Kitty_pp_spd,SourceSize);
    DLM_CleverMcLevyResoTM* CRS_pp_spd = NULL;
    if(SourceType=="CoreReso") CRS_pp_spd = CATS_ResoSource_pp(Kitty_pp_spd,SourceSize,2,true,0.6422,1.65,1362,200);
    CATS_pp_Basic(Kitty_pp_spd,100,0,400);
    CATS_pp_AV18(Kitty_pp_spd,true,true);
    printf(" Executing KillTheCat for Kitty_pp_spd\n");
    Kitty_pp_spd.SetNotifications(CATS::nWarning);
    Kitty_pp_spd.KillTheCat();
    DLM_Ck Ck_pp_spd(Kitty_pp_spd.GetNumSourcePars(),0,Kitty_pp_spd);
    Ck_pp_spd.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pp_spdWaves", &Ck_pp_spd);

    //we will include feed-downs from Lambda and Xi, thus we need to compute pL and pXi correlations as well
    CATS Kitty_pL_Usmani;
    CATS_GaussSource(Kitty_pL_Usmani,SourceSize_pL);
    CATS_pL_Basic(Kitty_pL_Usmani,100,0,400);
    CATS_pL_Usmani(Kitty_pL_Usmani);
    printf(" Executing KillTheCat for Kitty_pL_Usmani\n");
    Kitty_pL_Usmani.SetNotifications(CATS::nWarning);
    Kitty_pL_Usmani.KillTheCat();
    DLM_Ck Ck_pL_Usmani(Kitty_pL_Usmani.GetNumSourcePars(),0,Kitty_pL_Usmani);
    Ck_pL_Usmani.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pL_Usmani", &Ck_pL_Usmani);

    CATS Kitty_pXim_Hal;
    CATS_GaussSource(Kitty_pXim_Hal,SourceSize_pXim);
    CATS_pXim_Basic(Kitty_pXim_Hal,100,0,400);
    CATS_pXim_Hal(Kitty_pXim_Hal);
    printf(" Executing KillTheCat for Kitty_pXim_Hal\n");
    Kitty_pXim_Hal.SetNotifications(CATS::nWarning);
    Kitty_pXim_Hal.KillTheCat();
    DLM_Ck Ck_pXim_Hal(Kitty_pXim_Hal.GetNumSourcePars(),0,Kitty_pXim_Hal);
    Ck_pXim_Hal.Update();
    RootFile_DlmCk(OutputFileName, "Ck_pXim_Hal", &Ck_pXim_Hal);

    ///next we only concentrate on the last example (s,p,d waves) to include momentum resolution and feed-down effects
    //this is done using the class DLM_CkDecomposition. To set it up, we need the momentum smearing matrix, as well as any residual transformation matrices
    //In short: this object can smear the DLM_Ck histogram, and build the correlation function with the inclusion of feed-down channels. The feed-downs are
    //transfered recursively, i.e. say you define for pL a feed down from pXi and for pp feed down for pL, than pXi->pL->pp is automatically included
    TString FileName_MomSmear = "../Files/pp13TeV_MB.root";
    TString HistoName_MomSmear = "hSigmaMeV_Proton_Proton";
    TH2F* hResolution_pp = GetSmearMatrix(FileName_MomSmear,HistoName_MomSmear);

    TString FileName_Feed_pL_pp = "../Files/DecayMatrices_Oli.root";
    TString HistoName_Feed_pL_pp = "hRes_pp_pL";
    TH2F* hFeedDown_pL_pp = GetSmearMatrix(FileName_Feed_pL_pp,HistoName_Feed_pL_pp);

    TString FileName_Feed_pXim_pL = "../Files/DecayMatrices_Oli.root";
    TString HistoName_Feed_pXim_pL = "hRes_pL_pXim";
    TH2F* hFeedDown_pXim_pL = GetSmearMatrix(FileName_Feed_pXim_pL,HistoName_Feed_pXim_pL);

    //the object is initialized by: a name of your choice, number of contributions to C(k) APART from the primary, DLM_Ck, resolution matrix
    DLM_CkDecomposition CkDec_pp("pp",3,Ck_pp_spd,hResolution_pp);
    DLM_CkDecomposition CkDec_pL("pLambda",3,Ck_pL_Usmani,NULL);
    DLM_CkDecomposition CkDec_pXim("pXim",2,Ck_pXim_Hal,NULL);
    //example lambda pars from Phys. Lett. B 797 (2019) 134822
    double lambda_pL_pp = 0.151;
    double lambda_pp_flat = 0.081;
    double lambda_pp_misid = 0.02;
    //the primary lambda is assumed to be: 1. - lambda_pL_pp - lambda_pp_flat - lambda_pp_misid
    //to add the different contributions to the correlation, we need to specify:
    //#WhichContribution,#Type(cFeedDown or cFake),OPTIONALLY: DLM_CkDecomposition of the feed-down contribution, transformation matrix
    //if the last two entries are not given, it is assumed that the feed-down correlation is flat
    CkDec_pp.AddContribution(0,lambda_pL_pp,DLM_CkDecomposition::cFeedDown,&CkDec_pL,hFeedDown_pL_pp);
    CkDec_pp.AddContribution(1,lambda_pp_flat,DLM_CkDecomposition::cFeedDown);
    CkDec_pp.AddContribution(2,lambda_pp_misid,DLM_CkDecomposition::cFake);

    double lambda_pXim_pL = 0.083;
    double lambda_pL_flat = 0.372;
    double lambda_pL_misid = 0.042;
    CkDec_pL.AddContribution(0,lambda_pXim_pL,DLM_CkDecomposition::cFeedDown,&CkDec_pXim,hFeedDown_pXim_pL);
    CkDec_pL.AddContribution(1,lambda_pL_flat,DLM_CkDecomposition::cFeedDown);
    CkDec_pL.AddContribution(2,lambda_pL_misid,DLM_CkDecomposition::cFake);

    double lambda_pXim_flat = 0.391;
    double lambda_pXim_misid = 0.054;
    CkDec_pXim.AddContribution(0,lambda_pXim_flat,DLM_CkDecomposition::cFeedDown);
    CkDec_pXim.AddContribution(1,lambda_pXim_misid,DLM_CkDecomposition::cFake);

    //computes everything, including the decomposition of the correlation signal C(k)-1,
    //i.e. how much is each contribution to the total correlation
    CkDec_pp.Update(true,true);
    CkDec_pL.Update(true,true);
    CkDec_pXim.Update(true,true);
    RootFile_DlmCk(OutputFileName, "CkDec_pp_spdWaves", &CkDec_pp, true);

    //Now lets fit some correlation function
    //the example is for our analysis of pp 13 TeV MB p-p correlations
    //get the correlation function from a file (x-axis should be in MeV)
    TString DataFileName = "../Files/CFOutput_pp.root";
    TString DataHistoName = "hCk_ReweightedMeV_0";
    TH1F* hExpCk = GetExpCorrelation(DataFileName,DataHistoName);
    //small but important detail: turn off the standard CATS output, leaving only Warnings and Errors
    //this is needed to avoid a ton of output during the fitting!
    Kitty_pp_spd.SetNotifications(CATS::nWarning);
    //set the pointer to the object that we want to use in ExampleFitter_pp
    Ck_ExampleFitter_pp = &CkDec_pp;
    printf(" Fitting begins now!\n");
    //we set up a TF1 fitter, which will use a custom function made such as to use the
    //pointer to the DLM_CkDecomposition of the p-p to perform the fit, which has 3 fit parameters
    TF1* fit_pp = new TF1("fit_pp",ExampleFitter_pp,FitMin,FitMax,3);
    //source size
    fit_pp->SetParameter(0,SourceSize);
    //normalization of the baseline
    fit_pp->SetParameter(1,1.);
    fit_pp->SetParLimits(1,0.5,2.0);
    //slope (pol1) for the baseline
    //fit_pp->SetParameter(2,0.);
    //fit_pp->SetParLimits(2,-1e-3,1e-3);
    //in this example we fit with norm only
    fit_pp->FixParameter(2,0.);
    hExpCk->Fit(fit_pp, "S, N, R, M");

    printf(" The extracted source radius is: %.3f +/- %.3f fm\n",fit_pp->GetParameter(0),fit_pp->GetParError(0));
    printf(" The extracted BL norm is: %.3f +/- %.3f\n",fit_pp->GetParameter(1),fit_pp->GetParError(1));
    printf(" The extracted BL slope is: %.3e +/- %.3e 1/MeV\n",fit_pp->GetParameter(2),fit_pp->GetParError(2));
    printf(" chi2/ndf = %.2f/%i = %.2f\n",fit_pp->GetChisquare(),fit_pp->GetNDF(),fit_pp->GetChisquare()/double(fit_pp->GetNDF()));

    TFile* fOutput = new TFile(OutputFileName,"update");
    hResolution_pp->Write("",TObject::kOverwrite);
    hFeedDown_pL_pp->Write("",TObject::kOverwrite);
    hFeedDown_pXim_pL->Write("",TObject::kOverwrite);
    hExpCk->Write("",TObject::kOverwrite);
    fit_pp->Write("",TObject::kOverwrite);

    delete hResolution_pp;
    delete hFeedDown_pL_pp;
    delete hFeedDown_pXim_pL;
    delete hExpCk;
    delete fit_pp;
    delete fOutput;

    //don't forget to delete the source ones its not needed
    delete CRS_pp;
    delete CRS_pp_sp;
    delete CRS_pp_spd;

    //print out the timer information
    long long ExeTime = TIMER.Stop()/1000.;
    char* strtime = new char [128];
    ShowTime(ExeTime,strtime,0,true,6);
    printf("\nExecution time of Ck_pp_Decomposition for SourceType='%s': %s\n",SourceType.Data(),strtime);
    delete [] strtime;

}
