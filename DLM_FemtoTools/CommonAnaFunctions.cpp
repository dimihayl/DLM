
#include "CommonAnaFunctions.h"
#include "CATS.h"
#include "CATStools.h"
#include "CATSconstants.h"
#include "DLM_Source.h"
#include "DLM_CkModels.h"
#include "DLM_CkDecomposition.h"
#include "DLM_Ck.h"
#include "DLM_Potentials.h"
#include "DLM_WfModel.h"
#include "DLM_Random.h"

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TVector3.h"

DLM_CommonAnaFunctions::DLM_CommonAnaFunctions():NumCleverLevyObjects(6){
    //Simple_Reso = NULL;
    //Simple_Reso = new MS_GaussExp_mT_Simple [NumCleverLevyObjects];
    CleverLevy = NULL;
    CleverLevy = new DLM_CleverLevy [NumCleverLevyObjects];
    CleverMcLevyReso = NULL;
    CleverMcLevyReso = new DLM_CleverMcLevyReso [NumCleverLevyObjects];
    CleverMcLevyResoTM = NULL;
    CleverMcLevyResoTM = new DLM_CleverMcLevyResoTM [NumCleverLevyObjects];
    CatsFilesFolder = new TString();
}

DLM_CommonAnaFunctions::~DLM_CommonAnaFunctions(){
    //if(Simple_Reso){delete[]Simple_Reso;Simple_Reso=NULL;}
    if(CleverLevy){delete[]CleverLevy;CleverLevy=NULL;}
    if(CleverMcLevyReso){delete[]CleverMcLevyReso;CleverMcLevyReso=NULL;}
    if(CleverMcLevyResoTM){delete[]CleverMcLevyResoTM;CleverMcLevyResoTM=NULL;}
    delete CatsFilesFolder;
}

//POT:
//  "AV18", no pot vars so far
//  no sor var so far
void DLM_CommonAnaFunctions::SetUpCats_pp(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){

    CATSparameters* cPars = NULL;

    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3P0 = NULL;
    CATSparameters* cPotPars3P1 = NULL;
    CATSparameters* cPotPars3P2 = NULL;
    CATSparameters* cPotPars1D2 = NULL;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussTheta"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.6)*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.6*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0,sqrt(1.6)*1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[0].InitStability(20,1,2);
        CleverLevy[0].InitScale(35,0.25,2.0);
        CleverLevy[0].InitRad(256,0,64);
        CleverLevy[0].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.6*1.2);
        Kitty.SetAnaSource(1,1.6);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        //printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Gauss_mT_Reso)\n");
        cPars = new CATSparameters(CATSparameters::tSource,11,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.65);//tau
        cPars->SetParameter(2,1.-0.3578);//prim
        cPars->SetParameter(3,1361.52);//reso mass
        cPars->SetParameter(4,Mass_p);
        cPars->SetParameter(5,Mass_pic);
        cPars->SetParameter(6,1.65);
        cPars->SetParameter(7,1.-0.3578);
        cPars->SetParameter(8,1361.52);
        cPars->SetParameter(9,Mass_p);
        cPars->SetParameter(10,Mass_pic);
        Kitty.SetAnaSource(GaussExpTotSimple_2body, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[0].InitStability(21,1,2);
        CleverMcLevyReso[0].InitScale(38,0.15,2.0);
        CleverMcLevyReso[0].InitRad(257,0,64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0,1);
        CleverMcLevyReso[0].InitReso(1,1);
        if(SourceVar==0){
            CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else{
            CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }

        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            //DLM_Histo<double>* HISTO = ConvertThetaAngleHisto(
            //"/home/dmihaylov/Dudek_Ubuntu/MyApps/SomeTests/RandomCos/fCOS_180.root","hTheta2D",400,700);
            CleverMcLevyReso[0].SetUpResoEmission(0,0,HISTO);
            CleverMcLevyReso[0].SetUpResoEmission(1,0,HISTO);
        }
        CleverMcLevyReso[0].InitNumMcIter(200000);

        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[0].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[0].InitScale(38,0.15,2.0);
        CleverMcLevyReso[0].InitRad(257*2,0,64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0,1);
        CleverMcLevyReso[0].InitReso(1,1);
        if(SourceVar==0){
            CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else{
            CleverMcLevyReso[0].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }

        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto(
            "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,600);
            //DLM_Histo<double>* HISTO = ConvertThetaAngleHisto(
            //"/home/dmihaylov/Dudek_Ubuntu/MyApps/SomeTests/RandomCos/fCOS_Flat.root","hTheta2D",400,700);
            CleverMcLevyReso[0].SetUpResoEmission(0,0,HISTO);
            CleverMcLevyReso[0].SetUpResoEmission(1,0,HISTO);
        }

        CleverMcLevyReso[0].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"||SOURCE=="McLevy_ResoTM"){
        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[0].InitStability(1,2-1e-6,2+1e-6);
        else CleverMcLevyResoTM[0].InitStability(21,1,2);
        CleverMcLevyResoTM[0].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[0].InitRad(257*2,0,64);
        CleverMcLevyResoTM[0].InitType(2);
        CleverMcLevyResoTM[0].SetUpReso(0,0.6422);
        CleverMcLevyResoTM[0].SetUpReso(1,0.6422);
        //pure Gauss
        if(SourceVar%100==0){
//printf("Hello 0\n");
        }
        //back-to-back
        else if(SourceVar%100==1){
            CleverMcLevyResoTM[0].AddBGT_PR(490./1362.*1.65,1.);
            CleverMcLevyResoTM[0].AddBGT_RP(490./1362.*1.65,-1.);
            CleverMcLevyResoTM[0].AddBGT_RR(490./1362.*1.65,-1.,490./1362.*1.65,1.,-1.);
        }
        //EPOS, 2 is with fixed mass, 3 is with EPOS mass, 4 is 3 body with fixed mass, 5 is 3 body with EPOS mass
        else{
//printf("Hello 2\n");
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;

            TFile* F_EposDisto_p_pReso;
            if(SourceVar%100==4||SourceVar%100==5) F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/Only3body_p_pReso_3body.root");
            else F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_p_pReso.root");
//printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/").Data());
            TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                Tau2 = 1.65;
                if(SourceVar%100==2){
                    fM2 = 1362;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[0].AddBGT_PR(RanVal1,-cos(AngleRcP2));
                CleverMcLevyResoTM[0].AddBGT_RP(RanVal1,cos(AngleRcP2));
                //CleverMcLevyResoTM[0].AddBGT_PR(RanVal1,cos(AngleRcP2));
                //CleverMcLevyResoTM[0].AddBGT_RP(RanVal1,-cos(AngleRcP2));
            }
            delete F_EposDisto_p_pReso;

            TFile* F_EposDisto_pReso_pReso;
            if(SourceVar%100==4||SourceVar%100==5) F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/Only3body_pReso_pReso_3body.root");
            else F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_pReso_pReso.root");
//printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_pReso_pReso_3body.root").Data());
            TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 1.65;
                if(SourceVar%100==2){
                    fM1 = 1362;
                    fM2 = 1362;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[0].AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[0].InitNumMcIter(1000000);
        else CleverMcLevyResoTM[0].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[0], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        //printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        Kitty.SetUseAnalyticSource(false);
        //AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
        Kitty.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19");
        //AB_pp.SetInputFileName(THERMALFILE.Data());
        Kitty.SetMaxPairsToRead(64e6);
        Kitty.SetMaxPairsPerBin(16000/2);
        Kitty.SetMixingDepth(8);
        Kitty.SetThetaDependentSource(false);
        Kitty.SetTransportRenorm(1);
        Kitty.SetTauCorrection(false);
        //goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="EPOStheta"){
        //printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        Kitty.SetUseAnalyticSource(false);
        //AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
        Kitty.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19");
        //AB_pp.SetInputFileName(THERMALFILE.Data());
        Kitty.SetMaxPairsToRead(64e6);
        Kitty.SetMaxPairsPerBin(16000/2);
        Kitty.SetMixingDepth(8);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetTransportRenorm(1);
        Kitty.SetGridEpsilon(1./8192);
        Kitty.SetGridMaxDepth(8);
        Kitty.SetTauCorrection(false);
        //goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pp;
    }

    if(POT=="AV18"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
        double PotPars3P0[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
        double PotPars3P1[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
        double PotPars3P2[8]={NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};
        double PotPars1D2[8]={NN_AV18,v18_Coupled3P2,1,1,1,0,2,2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if(POT=="ReidV8"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={NN_ReidV8,v18_Coupled3P2,1,1,1,0,0,0};
        double PotPars3P0[8]={NN_ReidV8,v18_Coupled3P2,1,1,1,1,1,0};
        double PotPars3P1[8]={NN_ReidV8,v18_Coupled3P2,1,1,1,1,1,1};
        double PotPars3P2[8]={NN_ReidV8,v18_Coupled3P2,1,1,1,1,1,2};
        double PotPars1D2[8]={NN_ReidV8,v18_Coupled3P2,1,1,1,0,2,2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if(POT=="ReidSC"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={pp_ReidSC,0,1,1,1,0,0,0};
        double PotPars3P0[8]={pp_ReidSC,0,1,1,1,1,1,0};
        double PotPars3P1[8]={pp_ReidSC,0,1,1,1,1,1,1};
        double PotPars3P2[8]={pp_ReidSC,0,1,1,1,1,1,2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3P2->SetParameters(PotPars3P2);
    }
    else if(POT=="Norfolk"){
//printf("SetUpNorfolk...\n");
        SetUpNorfolk("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/NorfolkPotential/OriginalCode/Fine/");
//printf("SetUpNorfolk!\n");
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,3,true);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential,3,true);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential,3,true);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential,3,true);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential,3,true);

        cPotPars1S0->SetParameter(0,PotVar/10);//lpot
        cPotPars3P0->SetParameter(0,PotVar/10);//lpot
        cPotPars3P1->SetParameter(0,PotVar/10);//lpot
        cPotPars3P2->SetParameter(0,PotVar/10);//lpot
        cPotPars1D2->SetParameter(0,PotVar/10);//lpot

        //lemp
        if(PotVar%10==1){
            cPotPars1S0->SetParameter(1,100);
            cPotPars3P0->SetParameter(1,100);
            cPotPars3P1->SetParameter(1,100);
            cPotPars3P2->SetParameter(1,100);
            cPotPars1D2->SetParameter(1,100);
        }
        else if(PotVar%10==2){
            cPotPars1S0->SetParameter(1,101);
            cPotPars3P0->SetParameter(1,101);
            cPotPars3P1->SetParameter(1,101);
            cPotPars3P2->SetParameter(1,101);
            cPotPars1D2->SetParameter(1,101);
        }
        else{
            cPotPars1S0->SetParameter(1,-1);
            cPotPars3P0->SetParameter(1,-1);
            cPotPars3P1->SetParameter(1,-1);
            cPotPars3P2->SetParameter(1,-1);
            cPotPars1D2->SetParameter(1,-1);
        }

        cPotPars1S0->SetParameter(2,100);
        cPotPars3P0->SetParameter(2,310);
        cPotPars3P1->SetParameter(2,311);
        cPotPars3P2->SetParameter(2,312);
        cPotPars1D2->SetParameter(2,122);
    }

    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pp;
    }
    Kitty.SetMomentumDependentSource(false);
    //Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 2212);
    Kitty.SetRedMass( 0.5*Mass_p );

    Kitty.SetNumChannels(4);
    Kitty.SetNumPW(0,3);
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

    if(POT=="Norfolk"){
        if(cPotPars1S0) Kitty.SetShortRangePotential(0,0,pp_Norfolk,*cPotPars1S0);
        if(cPotPars1D2) Kitty.SetShortRangePotential(0,2,pp_Norfolk,*cPotPars1D2);
        if(cPotPars3P0) Kitty.SetShortRangePotential(1,1,pp_Norfolk,*cPotPars3P0);
        if(cPotPars3P1) Kitty.SetShortRangePotential(2,1,pp_Norfolk,*cPotPars3P1);
        if(cPotPars3P2) Kitty.SetShortRangePotential(3,1,pp_Norfolk,*cPotPars3P2);
    }
    else{
        if(cPotPars1S0) Kitty.SetShortRangePotential(0,0,fDlmPot,*cPotPars1S0);
        if(cPotPars1D2) Kitty.SetShortRangePotential(0,2,fDlmPot,*cPotPars1D2);
        if(cPotPars3P0) Kitty.SetShortRangePotential(1,1,fDlmPot,*cPotPars3P0);
        if(cPotPars3P1) Kitty.SetShortRangePotential(2,1,fDlmPot,*cPotPars3P1);
        if(cPotPars3P2) Kitty.SetShortRangePotential(3,1,fDlmPot,*cPotPars3P2);
    }

    CLEAN_SetUpCats_pp: ;
    if(cPars){delete cPars; cPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars1S0){delete cPotPars1S0; cPotPars1S0=NULL;}
    if(cPotPars1D2){delete cPotPars1D2; cPotPars1D2=NULL;}
    if(cPotPars3P0){delete cPotPars3P0; cPotPars3P0=NULL;}
    if(cPotPars3P1){delete cPotPars3P1; cPotPars3P1=NULL;}
    if(cPotPars3P2){delete cPotPars3P2; cPotPars3P2=NULL;}

}

void DLM_CommonAnaFunctions::SetUpCats_pipi(CATS& Kitty, const TString& SOURCE, const int& SourceVar){

    CATSparameters* cPars = NULL;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussTheta"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.6)*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.6*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0,sqrt(1.6)*1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.6*1.2);
        Kitty.SetAnaSource(1,1.6);
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McLevyNolan_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="McGauss_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"||SOURCE=="McLevy_ResoTM"){
        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[4].InitStability(1,2-1e-6,2+1e-6);
        else CleverMcLevyResoTM[4].InitStability(21,1,2);
        CleverMcLevyResoTM[4].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[4].InitRad(257*2,0,64);
        CleverMcLevyResoTM[4].InitType(2);
        CleverMcLevyResoTM[4].SetUpReso(0,0.6422);
        CleverMcLevyResoTM[4].SetUpReso(1,0.6422);
        //pure Gauss
        if(SourceVar%100==0){
        }
        //back-to-back
        else if(SourceVar%100==1){
            printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_ResoTM back-to-back)\n");
            goto CLEAN_SetUpCats_pipi;
        }
        //EPOS, 2(4) is with fixed mass, 3 is with EPOS mass
        //2 is with omega included, 4 is without omega, 8 is neglecting the angular distr. from EPOS and taking them random
        else{
unsigned NumPart=0;
unsigned NumOmega=0;
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;

            TFile* F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/ForMax_pi_piReso_withOmega.root");
            TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
NumPart++;
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                //treat the omega separately
                if(fM2>782&&fM2<783){
NumOmega++;
                    fM2 = 782.6;
                    Tau2 = 23.24;
                    if(SourceVar%100==4) continue;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau2 = 1.5;
                    if(SourceVar%100==2 || SourceVar%100==4){
                        fM2 = 1124;
                    }
                }
                if(k_D>k_CutOff) continue;
//if(fM2>782&&fM2<783){
//printf("omega\n");
//}
                RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
                if(SourceVar%100==8){
                    TVector3 vRCORE(0,0,1.);
                    TVector3 vSR2;
                    vSR2.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1)),RanGen.Uniform(0,2.*Pi));
                    AngleRcP2 = vRCORE.Angle(vSR2);
                }
                CleverMcLevyResoTM[4].AddBGT_PR(RanVal1,-cos(AngleRcP2));
                CleverMcLevyResoTM[4].AddBGT_RP(RanVal1,cos(AngleRcP2));
            }
printf("%u/%u = %f\n",NumOmega,NumPart,double(NumOmega)/double(NumPart));
            delete F_EposDisto_p_pReso;

            TFile* F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/ForMax_piReso_piReso_withOmega.root");
            TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                //treat the omega separately
                if(fM1>782&&fM1<783){
                    fM1 = 782.6;
                    Tau1 = 23.24;
                    if(SourceVar%100==4) continue;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau1 = 1.5;
                    if(SourceVar%100==2 || SourceVar%100==4){
                        fM1 = 1124;
                    }
                }
                //treat the omega separately
                if(fM2>782&&fM2<783){
                    fM2 = 782.6;
                    Tau2 = 23.24;
                    if(SourceVar%100==4) continue;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau2 = 1.5;
                    if(SourceVar%100==2 || SourceVar%100==4){
                        fM2 = 1124;
                    }
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                if(SourceVar%100==8){
                    TVector3 vRCORE(0,0,1.);
                    TVector3 vSR1;
                    vSR1.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1)),RanGen.Uniform(0,2.*Pi));
                    TVector3 vSR2;
                    vSR2.SetMagThetaPhi(1,asin(RanGen.Uniform(-1,1)),RanGen.Uniform(0,2.*Pi));
                    AngleRcP1 = vRCORE.Angle(vSR1);
                    AngleRcP2 = vRCORE.Angle(vSR2);
                    AngleP1P2 = vSR1.Angle(vSR2);
                }
                CleverMcLevyResoTM[4].AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[4].InitNumMcIter(1000000);
        else CleverMcLevyResoTM[4].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[4], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="EPOStheta"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOStheta)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pipi;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(211, 211);
    Kitty.SetRedMass( 0.5*Mass_pic );

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,0);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1.);

    CLEAN_SetUpCats_pipi: ;
    if(cPars){delete cPars; cPars=NULL;}

}
void DLM_CommonAnaFunctions::SetUpCats_pipi_broken(CATS& Kitty, const TString& SOURCE, const int& SourceVar){

    CATSparameters* cPars = NULL;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussTheta"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.6)*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.6*1.2);
        cPars->SetParameter(1,1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0,sqrt(1.6)*1.2);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[4].InitStability(20,1,2);
        CleverLevy[4].InitScale(35,0.25,2.0);
        CleverLevy[4].InitRad(256,0,64);
        CleverLevy[4].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.6*1.2);
        Kitty.SetAnaSource(1,1.6);
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McLevyNolan_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="McGauss_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"||SOURCE=="McLevy_ResoTM"){
        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[4].InitStability(1,2-1e-6,2+1e-6);
        else CleverMcLevyResoTM[4].InitStability(21,1,2);
        CleverMcLevyResoTM[4].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[4].InitRad(257*2,0,64);
        CleverMcLevyResoTM[4].InitType(2);
        CleverMcLevyResoTM[4].SetUpReso(0,0.6422);
        CleverMcLevyResoTM[4].SetUpReso(1,0.6422);
        //pure Gauss
        if(SourceVar%100==0){
        }
        //back-to-back
        else if(SourceVar%100==1){
            printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_ResoTM back-to-back)\n");
            goto CLEAN_SetUpCats_pipi;
        }
        //EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else{
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;
printf("SetUpCats_pipi_broken\n");
            TFile* F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/WithoutOmega/ForMax_pi_piReso.root");
            TNtuple* T_EposDisto_p_pReso = (TNtuple*)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_p_pReso; uEntry++){
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                //treat the omega separately
                if(fM2>782&&fM2<783){
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau2 = 1.5;
                    if(SourceVar%100==2){
                        fM2 = 1124;
                    }
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[4].AddBGT_PR(RanVal1,-cos(AngleRcP2));
                CleverMcLevyResoTM[4].AddBGT_RP(RanVal1,cos(AngleRcP2));
            }
            delete F_EposDisto_p_pReso;

            TFile* F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/WithoutOmega/ForMax_piReso_piReso.root");
            TNtuple* T_EposDisto_pReso_pReso = (TNtuple*)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_pReso; uEntry++){
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                //treat the omega separately
                if(fM2>782&&fM2<783){
                    fM1 = 782.6;
                    Tau1 = 23.24;
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                //the avg. values below should be computed for all resonances besides omega
                else{
                    Tau1 = 1.5;
                    Tau2 = 1.5;
                    if(SourceVar%100==2){
                        fM1 = 1124;
                        fM2 = 1124;
                    }
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[4].AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[4].InitNumMcIter(1000000);
        else CleverMcLevyResoTM[4].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[4], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="EPOStheta"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOStheta)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pipi;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(211, 211);
    Kitty.SetRedMass( 0.5*Mass_pic );

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0,0);
    Kitty.SetSpin(0,0);
    Kitty.SetChannelWeight(0, 1.);

    CLEAN_SetUpCats_pipi: ;
    if(cPars){delete cPars; cPars=NULL;}

}

//POT:
//  "LO"
//  "LO_Coupled_S"
//  "NLO"
//  "NLO_Coupled_S"
//  "Usmani"
//SourceVar:
//  0 = back-to-back
//  1 = flat theta
//  2 = EPOS theta
//PotVar (for Chiral_Coupled_SPD):
//  4 digits ABXXX
//  A is NLO13 (0), NLO19 (1), LO13 (-1)
//  B is: 0 (s waves), 1 (sd waves), 2 (spd waves)
//  XXX are the cutoff
void DLM_CommonAnaFunctions::SetUpCats_pL(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){
    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotPars1S0 = NULL;
    CATSparameters* cPotPars3S1 = NULL;

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    unsigned NumChannels=0;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussTheta"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.2)*1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.2*1.2);
        cPars->SetParameter(1,1.2);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,sqrt(1.2)*1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[1].InitStability(20,1,2);
        CleverLevy[1].InitScale(35,0.25,2.0);
        CleverLevy[1].InitRad(256,0,64);
        CleverLevy[1].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0,0.5*1.2*1.2);
        Kitty.SetAnaSource(1,1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[1].InitStability(21,1,2);
        CleverMcLevyReso[1].InitScale(38,0.15,2.0);
        CleverMcLevyReso[1].InitRad(257,0,64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0,1);
        CleverMcLevyReso[1].InitReso(1,1);
//bool MassSmear = SourceVar%10;
//double MomSmear = SourceVar/100;
        if(SourceVar==0){
            CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else{
            CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }

        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            CleverMcLevyReso[1].SetUpResoEmission(0,0,HISTO);
            CleverMcLevyReso[1].SetUpResoEmission(1,0,HISTO);
        }
        CleverMcLevyReso[1].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[1].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[1].InitScale(35,0.25,2.0);
        CleverMcLevyReso[1].InitRad(257*2,0,64);
        CleverMcLevyReso[1].InitType(2);
        if(SourceVar==0){
            CleverMcLevyReso[1].InitReso(0,1);
            CleverMcLevyReso[1].InitReso(1,1);
            CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else if(SourceVar==1){
            CleverMcLevyReso[1].InitReso(0,1);
            CleverMcLevyReso[1].InitReso(1,1);
            CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }
        if(SourceVar==2){
            CleverMcLevyReso[1].InitReso(0,1);
            CleverMcLevyReso[1].InitReso(1,1);
            CleverMcLevyReso[1].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1,0,1.-0.3562,1462.93,4.69,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double>* HISTO_p = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",400,600);
            DLM_Histo<double>* HISTO_L = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",270,470);

            CleverMcLevyReso[1].SetUpResoEmission(0,0,HISTO_p);
            CleverMcLevyReso[1].SetUpResoEmission(1,0,HISTO_L);
        }
        //improved reso
        else if(SourceVar==3){
            CleverMcLevyReso[1].InitReso(0,2);
            CleverMcLevyReso[1].InitReso(1,2);
            CleverMcLevyReso[1].SetUpReso(0,0,0.4368,1231.98,1.67,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(0,1,0.2054,1636.98,1.62,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1,0,0.4864,1384.54,5.34,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1,1,0.1573,1705.26,2.70,Mass_L,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double>* HISTO_p_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",200,400);
            DLM_Histo<double>* HISTO_p_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",880,1000);
            DLM_Histo<double>* HISTO_L_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",160,360);
            DLM_Histo<double>* HISTO_L_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root","h_rkAngle_Mom2",620,820);

            CleverMcLevyReso[1].SetUpResoEmission(0,0,HISTO_p_0);
            CleverMcLevyReso[1].SetUpResoEmission(0,1,HISTO_p_1);
            CleverMcLevyReso[1].SetUpResoEmission(1,0,HISTO_L_0);
            CleverMcLevyReso[1].SetUpResoEmission(1,1,HISTO_L_1);
        }

        CleverMcLevyReso[1].InitNumMcIter(SourceVar==3?100000:1000000);

        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"||SOURCE=="McLevy_ResoTM"){
        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[1].InitStability(1,2-1e-6,2+1e-6);
        else CleverMcLevyResoTM[1].InitStability(21,1,2);
        CleverMcLevyResoTM[1].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[1].InitRad(257*2,0,64);
        CleverMcLevyResoTM[1].InitType(2);
        CleverMcLevyResoTM[1].SetUpReso(0,0.6422);
        CleverMcLevyResoTM[1].SetUpReso(1,0.6438);
        //pure Gauss
        if(SourceVar%100==0){

        }
        //back-to-back
        else if(SourceVar%100==1){
            CleverMcLevyResoTM[1].AddBGT_PR(360./1462.*4.69,1.);
            CleverMcLevyResoTM[1].AddBGT_RP(490./1362.*1.65,-1.);
            CleverMcLevyResoTM[1].AddBGT_RR(490./1362.*1.65,-1.,360./1462.*4.69,1.,-1.);
        }
        //EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else{
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;

            TFile* F_EposDisto_p_LamReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_p_LamReso.root");
            TNtuple* T_EposDisto_p_LamReso = (TNtuple*)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
            T_EposDisto_p_LamReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_p_LamReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_p_LamReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_p_LamReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_p_LamReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_p_LamReso; uEntry++){
                T_EposDisto_p_LamReso->GetEntry(uEntry);
                Tau1 = 0;
                Tau2 = 4.69;
                if(SourceVar%100==2){
                    fM2 = 1462;
                }
                if(k_D>k_CutOff) continue;
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[1].AddBGT_PR(RanVal2,cos(AngleRcP2));
            }
            delete F_EposDisto_p_LamReso;

            TFile* F_EposDisto_pReso_Lam = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_pReso_Lam.root");
            TNtuple* T_EposDisto_pReso_Lam = (TNtuple*)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
            T_EposDisto_pReso_Lam->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_Lam->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_Lam->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Lam; uEntry++){
                T_EposDisto_pReso_Lam->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if(SourceVar%100==2){
                    fM1 = 1362;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                CleverMcLevyResoTM[1].AddBGT_RP(RanVal1,cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Lam;

            TFile* F_EposDisto_pReso_LamReso = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_pReso_LamReso.root");
            TNtuple* T_EposDisto_pReso_LamReso = (TNtuple*)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
            T_EposDisto_pReso_LamReso->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_LamReso; uEntry++){
                T_EposDisto_pReso_LamReso->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 4.69;
                if(SourceVar%100==2){
                    fM1 = 1362;
                    fM2 = 1462;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
                CleverMcLevyResoTM[1].AddBGT_RR(RanVal1,cos(AngleRcP1),RanVal2,cos(AngleRcP2),cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_LamReso;
        }

        if(SOURCE=="McGauss_ResoTM") CleverMcLevyResoTM[1].InitNumMcIter(1000000);
        else CleverMcLevyResoTM[1].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[1], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pL;
    }

    if(POT=="LO"){
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaLO_600/",
                                Kitty, 0, 600);
        NumChannels=2;
    }
    else if(POT=="LO_Coupled_S"){
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaLO_Coupling/",
                                Kitty, 1, 600);
        NumChannels=4;
    }
    else if(POT=="NLO"){
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaNLO/",
                                Kitty, 10, 600);
        NumChannels=2;
    }
    //s and p waves
    else if(POT=="NLO_sp"){
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaNLO/",
                                Kitty, 12, 600);
        NumChannels=4;
    }
    else if(POT=="NLO_Coupled_S"){
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambdaNLO_Coupling/",
                                Kitty, 11, 600);
        NumChannels=4;
    }
    else if(POT=="Chiral_Coupled_SPD"){
        int CUTOFF = abs(PotVar%1000);
        int TYPE = PotVar/10000;
        if(PotVar==0){CUTOFF=600;TYPE=0;}
//printf("PotVar=%i\n",PotVar);
//printf("CUTOFF=%i\n",CUTOFF);
//printf("TYPE=%i\n",TYPE);
        //original, NLO13 at 600 MeV
        ExternalWF = Init_pL_Haidenbauer2019(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pLambda_Coupled_SD/",
                                    Kitty, TYPE , CUTOFF);
        NumChannels=16;
    }
    else if(POT=="Usmani"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8]={pL_UsmaniOli,0,0,0,0,0,0,0};
        double PotPars3S1[8]={pL_UsmaniOli,0,0,0,0,1,0,1};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,8,true); cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels=2;
    }
    else if(POT=="UsmaniFit"){
        double PotPars1S0[4]={0,2137,0.5,0.2};
        double PotPars3S1[4]={1,2137,0.5,0.2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential,4,true); cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,4,true); cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels=2;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pL;
    }
//printf("NumChannels=%u\n",NumChannels);
//printf("ExternalWF=%p\n",ExternalWF);

    Kitty.SetMomentumDependentSource(false);
    //Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        if(!ExternalWF){
            Kitty.SetNumPW(uCh,1);
            Kitty.SetSpin(uCh, uCh%2==0?0:1);
            Kitty.SetChannelWeight(uCh, uCh%2==0?0.25:0.75);
        }

        if(POT=="UsmaniFit"&&cPotPars1S0&&uCh==0)Kitty.SetShortRangePotential(uCh,0,UsmaniFit,*cPotPars1S0);
        else if(POT=="UsmaniFit"&&cPotPars3S1&&uCh==1)Kitty.SetShortRangePotential(uCh,0,UsmaniFit,*cPotPars3S1);
        else if(cPotPars1S0&&uCh==0)Kitty.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars1S0);
        else if(cPotPars3S1&&uCh==1) Kitty.SetShortRangePotential(uCh,0,fDlmPot,*cPotPars3S1);
        else if(ExternalWF){
            //for(unsigned uMomBin=0; uMomBin<Kitty.GetNumMomBins(); uMomBin++){
                //Kitty.UseExternalWaveFunction(uMomBin,uCh,0,WaveFunctionU[uMomBin][uCh][0], NumRadBins, RadBins, PhaseShifts[uMomBin][uCh][0]);
//printf("Look at that view (%u)!\n",uCh);
                int SPD = (PotVar/1000)%10;
                //unless p-wave
                if( ExternalWF[0][uCh][0].GetDim() ) Kitty.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);
                if(POT=="Chiral_Coupled_SPD"){

                    if(uCh<=6&&SPD==2){
                        Kitty.SetExternalWaveFunction(uCh,1,ExternalWF[0][uCh][1],ExternalWF[1][uCh][1]);
//printf("1: uCh=%u\n",uCh);
                    }
                    if(uCh>=1&&uCh<=3&&SPD!=0){
                        Kitty.SetExternalWaveFunction(uCh,2,ExternalWF[0][uCh][2],ExternalWF[1][uCh][2]);
//printf("uCh=%u\n",uCh);
                    }
                }
                else if(Kitty.GetNumPW(uCh)>=2){
                    Kitty.SetExternalWaveFunction(uCh,1,ExternalWF[0][uCh][1],ExternalWF[1][uCh][1]);
                }

//printf(" --Look at that view (%u)!\n",uCh);
            //}
        }
        else{
printf("PotVar=%i\n",PotVar);
            printf("\033[1;31mERROR:\033[0m SetUpCats_pL says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pL;
        }

    }
//Kitty.KillTheCat();
//printf("------------------------");
    CLEAN_SetUpCats_pL: ;

    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars1S0){delete cPotPars1S0; cPotPars1S0=NULL;}
    if(cPotPars3S1){delete cPotPars3S1; cPotPars3S1=NULL;}
    CleanUpWfHisto(Kitty,ExternalWF);

}

//for the moment only gauss
//models are NLO and ESC16
void DLM_CommonAnaFunctions::SetUpCats_pS0(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){
    CATSparameters* cPars = NULL;

    DLM_Histo<complex<double>>*** ExternalWF=NULL;
    unsigned NumChannels=0;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pS0;
    }

    if(POT=="Chiral"){
        ExternalWF = Init_pSigma0_Haidenbauer(CatsFilesFolder[0]+"/Interaction/Haidenbauer/pSigma0/",Kitty);
        NumChannels=Kitty.GetNumChannels();
    }
    else if(POT=="ESC16"){
        ExternalWF = Init_pS0_ESC16(CatsFilesFolder[0]+"/Interaction/Tom/pSigma0/DimiValeNorm170519/",Kitty);
        NumChannels=Kitty.GetNumChannels();
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pS0;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        if(!ExternalWF){
            Kitty.SetSpin(uCh, uCh%2==0?0:1);
            Kitty.SetChannelWeight(uCh, uCh%2==0?0.25:0.75);
        }
        else if(ExternalWF){
            Kitty.SetExternalWaveFunction(uCh,0,ExternalWF[0][uCh][0],ExternalWF[1][uCh][0]);
        }
        else{
            printf("\033[1;31mERROR:\033[0m SetUpCats_pL says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pS0;
        }
    }
    CLEAN_SetUpCats_pS0: ;

    if(cPars){delete cPars; cPars=NULL;}
    CleanUpWfHisto(Kitty,ExternalWF);

}


//POT:
//  "pXim_Lattice" (the first version)
//  "pXim_HALQCD1" (the second version)
void DLM_CommonAnaFunctions::SetUpCats_pXim(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){

    double POTFLAG;
    switch(PotVar){
        case 11: POTFLAG=11; break;
        case 12: POTFLAG=12; break;
        case 13: POTFLAG=13; break;
        case -11: POTFLAG=-11; break;
        case -12: POTFLAG=-12; break;
        case -13: POTFLAG=-13; break;
        default: POTFLAG=12; break;
    }

    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotParsI0S0 = NULL;
    CATSparameters* cPotParsI0S1 = NULL;
    CATSparameters* cPotParsI1S0 = NULL;
    CATSparameters* cPotParsI1S1 = NULL;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.8)*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.8*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetAnaSource(0,sqrt(1.8)*1.2);
        Kitty.SetAnaSource(1,1.8);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[2].InitStability(20,1,2);
        CleverLevy[2].InitScale(35,0.25,2.0);
        CleverLevy[2].InitRad(256,0,64);
        CleverLevy[2].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.8*1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[2].InitStability(21,1,2);
        CleverMcLevyReso[2].InitScale(38,0.15,2.0);
        CleverMcLevyReso[2].InitRad(257,0,64);
        CleverMcLevyReso[2].InitType(2);
        CleverMcLevyReso[2].InitReso(0,1);
        CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic);
        if(SourceVar==0){
            CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
        }
        else{
            CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
        }
        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            CleverMcLevyReso[2].SetUpResoEmission(0,0,HISTO);
        }

        CleverMcLevyReso[2].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[2], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[2].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[2].InitScale(38,0.15,2.0);
        CleverMcLevyReso[2].InitRad(257,0,64);
        CleverMcLevyReso[2].InitType(2);
        CleverMcLevyReso[2].InitReso(0,1);
        CleverMcLevyReso[2].InitReso(1,1);
        if(SourceVar==0){
            CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[2].SetUpReso(1,0,0.0001,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
        }
        else{
            CleverMcLevyReso[2].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[2].SetUpReso(1,0,0.0001,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
        }
        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            CleverMcLevyReso[2].SetUpResoEmission(0,0,HISTO);
            CleverMcLevyReso[2].SetUpResoEmission(1,0,HISTO);
        }
        CleverMcLevyReso[2].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[2], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"){
        CleverMcLevyResoTM[2].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyResoTM[2].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[2].InitRad(257*2,0,64);
        CleverMcLevyResoTM[2].InitType(2);
        CleverMcLevyResoTM[2].SetUpReso(0,0.6422);
        //pure Gauss
        if(SourceVar%100==0){

        }
        //back-to-back
        else if(SourceVar%100==1){
            CleverMcLevyResoTM[2].AddBGT_RP(490./1362.*1.65,-1.);
        }
        //EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else{
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;

            TFile* F_EposDisto_pReso_Xim = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_pReso_Xim.root");
            TNtuple* T_EposDisto_pReso_Xim = (TNtuple*)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
            T_EposDisto_pReso_Xim->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_Xim->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_Xim->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_Xim->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_Xim->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_Xim->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_Xim->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Xim; uEntry++){
                T_EposDisto_pReso_Xim->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if(SourceVar%100==2){
                    fM1 = 1362;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                CleverMcLevyResoTM[2].AddBGT_RP(RanVal1,cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Xim;
        }

        CleverMcLevyResoTM[2].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[2], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pXim;
    }

    if(POT=="pXim_Lattice"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9]={pXim_Lattice,12,0,-1,1,0,0,0,0};
        double PotParsI0S1[9]={pXim_Lattice,12,0,-1,1,1,0,1,0};
        double PotParsI1S0[9]={pXim_Lattice,6,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_Lattice,6,1,1,1,1,0,1,0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if(POT=="pXim_HALQCD1"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9]={pXim_HALQCD1,POTFLAG,0,-1,1,0,0,0,0};
        double PotParsI0S1[9]={pXim_HALQCD1,POTFLAG,0,-1,1,1,0,1,0};
        double PotParsI1S0[9]={pXim_HALQCD1,POTFLAG,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_HALQCD1,POTFLAG,1,1,1,1,0,1,0};
//printf("POTFLAG = %f\n",POTFLAG);
//printf("PotVar = %i\n",PotVar);
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if(POT=="pXim_HALQCDPaper2020"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9]={pXim_HALQCDPaper2020,POTFLAG,0,-1,1,0,0,0,0};
        double PotParsI0S1[9]={pXim_HALQCDPaper2020,POTFLAG,0,-1,1,1,0,1,0};
        double PotParsI1S0[9]={pXim_HALQCDPaper2020,POTFLAG,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_HALQCDPaper2020,POTFLAG,1,1,1,1,0,1,0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if(POT=="pXim1530"){

    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pXim;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    if(POT=="pXim1530"){
      Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

      Kitty.SetNumChannels(1);
      Kitty.SetNumPW(0,0);
      Kitty.SetSpin(0,0);
      Kitty.SetChannelWeight(0,1.);
    }
    else{
      Kitty.SetRedMass( (Mass_p*Mass_Xim1530)/(Mass_p+Mass_Xim1530) );

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

      if(cPotParsI0S0) Kitty.SetShortRangePotential(0,0,fDlmPot,*cPotParsI0S0);
      if(cPotParsI0S1) Kitty.SetShortRangePotential(1,0,fDlmPot,*cPotParsI0S1);
      if(cPotParsI1S0) Kitty.SetShortRangePotential(2,0,fDlmPot,*cPotParsI1S0);
      if(cPotParsI1S1) Kitty.SetShortRangePotential(3,0,fDlmPot,*cPotParsI1S1);
    }



    CLEAN_SetUpCats_pXim: ;
    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotParsI0S0){delete cPotParsI0S0; cPotParsI0S0=NULL;}
    if(cPotParsI0S1){delete cPotParsI0S1; cPotParsI0S1=NULL;}
    if(cPotParsI1S0){delete cPotParsI1S0; cPotParsI1S0=NULL;}
    if(cPotParsI1S1){delete cPotParsI1S1; cPotParsI1S1=NULL;}
}

void DLM_CommonAnaFunctions::SetUpCats_pXi0(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){

    double POTFLAG;
    switch(PotVar){
        case 11: POTFLAG=11; break;
        case 12: POTFLAG=12; break;
        case 13: POTFLAG=13; break;
        case -11: POTFLAG=-11; break;
        case -12: POTFLAG=-12; break;
        case -13: POTFLAG=-13; break;
        default: POTFLAG=12; break;
    }

    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotParsI1S0 = NULL;
    CATSparameters* cPotParsI1S1 = NULL;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pXi0;
    }

    if(POT=="pXim_Lattice"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9]={pXim_Lattice,6,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_Lattice,6,1,1,1,1,0,1,0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if(POT=="pXim_HALQCD1"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9]={pXim_HALQCD1,POTFLAG,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_HALQCD1,POTFLAG,1,1,1,1,0,1,0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if(POT=="pXim_HALQCDPaper2020"){
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9]={pXim_HALQCDPaper2020,POTFLAG,1,1,1,0,0,0,0};
        double PotParsI1S1[9]={pXim_HALQCDPaper2020,POTFLAG,1,1,1,1,0,1,0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential,9,true); cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pXi0;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3322);
    Kitty.SetRedMass( (Mass_p*Mass_Xi0)/(Mass_p+Mass_Xi0) );

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(1,1);
    Kitty.SetChannelWeight(0, 1./4.);
    Kitty.SetChannelWeight(1, 3./4.);

    if(cPotParsI1S0) Kitty.SetShortRangePotential(0,0,fDlmPot,*cPotParsI1S0);
    if(cPotParsI1S1) Kitty.SetShortRangePotential(1,0,fDlmPot,*cPotParsI1S1);

    CLEAN_SetUpCats_pXi0: ;
    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotParsI1S0){delete cPotParsI1S0; cPotParsI1S0=NULL;}
    if(cPotParsI1S1){delete cPotParsI1S1; cPotParsI1S1=NULL;}
}


void DLM_CommonAnaFunctions::SetUpCats_pOmegam(CATS& Kitty, const TString& POT, const TString& SOURCE, const int& PotVar, const int& SourceVar){
    CATSparameters* cPars = NULL;
    CATSparameters* pPars = NULL;

    CATSparameters* cPotPars3S1 = NULL;
    CATSparameters* cPotPars5S2 = NULL;

    if(SOURCE=="Gauss"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Cauchy"){
        cPars = new CATSparameters(CATSparameters::tSource,1,true);
        cPars->SetParameter(0,1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Nolan"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Single"){
        cPars = new CATSparameters(CATSparameters::tSource,2,true);
        cPars->SetParameter(0,sqrt(1.8)*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="Levy_Diff"){
        cPars->SetParameter(0,0.5*1.8*1.2);
        cPars->SetParameter(1,1.8);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Nolan"){
        CleverLevy[3].InitStability(20,1,2);
        CleverLevy[3].InitScale(35,0.25,2.0);
        CleverLevy[3].InitRad(256,0,64);
        CleverLevy[3].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="CleverLevy_Single"){
        CleverLevy[3].InitStability(20,1,2);
        CleverLevy[3].InitScale(35,0.25,2.0);
        CleverLevy[3].InitRad(256,0,64);
        CleverLevy[3].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetAnaSource(0,sqrt(1.8)*1.2);
        Kitty.SetAnaSource(1,1.8);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="CleverLevy_Diff"){
        CleverLevy[3].InitStability(20,1,2);
        CleverLevy[3].InitScale(35,0.25,2.0);
        CleverLevy[3].InitRad(256,0,64);
        CleverLevy[3].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0,0.5*1.8*1.2);
        Kitty.SetAnaSource(1,1.8);
    }
    else if(SOURCE=="GaussExpTotSimple_2body"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if(SOURCE=="McLevyNolan_Reso"){
        CleverMcLevyReso[3].InitStability(21,1,2);
        CleverMcLevyReso[3].InitScale(38,0.15,2.0);
        CleverMcLevyReso[3].InitRad(257,0,64);
        CleverMcLevyReso[3].InitType(2);
        CleverMcLevyReso[3].InitReso(0,1);
        if(SourceVar==0){
            CleverMcLevyReso[3].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
        }
        else{
            CleverMcLevyReso[3].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
        }
        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            CleverMcLevyReso[3].SetUpResoEmission(0,0,HISTO);
        }
        CleverMcLevyReso[3].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[3], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="McGauss_Reso"){
        CleverMcLevyReso[3].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyReso[3].InitScale(38,0.15,2.0);
        CleverMcLevyReso[3].InitRad(257*2,0,64);
        CleverMcLevyReso[3].InitType(2);
        CleverMcLevyReso[3].InitReso(0,1);
        CleverMcLevyReso[3].InitReso(1,1);
        if(SourceVar==0){
            CleverMcLevyReso[3].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);//rdtBackwards,rdtRandom
            CleverMcLevyReso[3].SetUpReso(1,0,0.0001,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtBackwards);
        }
        else{
            CleverMcLevyReso[3].SetUpReso(0,0,1.-0.3578,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);//rdtBackwards,rdtRandom
            CleverMcLevyReso[3].SetUpReso(1,0,0.0001,1361.52,1.65,Mass_p,Mass_pic,false,false,DLM_CleverMcLevyReso::rdtRandom);
        }

        if(SourceVar==2){
            DLM_Histo<double>* HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root","h_rkAngle_Mom2",400,700);
            CleverMcLevyReso[3].SetUpResoEmission(0,0,HISTO);
            CleverMcLevyReso[3].SetUpResoEmission(1,0,HISTO);
        }

        CleverMcLevyReso[3].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[3], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    //SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if(SOURCE=="McGauss_ResoTM"){
        CleverMcLevyResoTM[3].InitStability(1,2-1e-6,2+1e-6);
        CleverMcLevyResoTM[3].InitScale(38,0.15,2.0);
        CleverMcLevyResoTM[3].InitRad(257*2,0,64);
        CleverMcLevyResoTM[3].InitType(2);
        CleverMcLevyResoTM[3].SetUpReso(0,0.6422);
        //pure Gauss
        if(SourceVar%100==0){

        }
        //back-to-back
        else if(SourceVar%100==1){
            CleverMcLevyResoTM[3].AddBGT_RP(490./1362.*1.65,-1.);
        }
        //EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else{
            const double k_CutOff = int(int(SourceVar)/10)*10.;
            Float_t k_D;
            Float_t fP1;
            Float_t fP2;
            Float_t fM1;
            Float_t fM2;
            Float_t Tau1;
            Float_t Tau2;
            Float_t AngleRcP1;
            Float_t AngleRcP2;
            Float_t AngleP1P2;
            DLM_Random RanGen(11);
            double RanVal1;
            double RanVal2;
            double RanVal3;

            TFile* F_EposDisto_pReso_Omega = new TFile(CatsFilesFolder[0]+"/Source/EposAngularDist/EposDisto_pReso_Omega.root");
            TNtuple* T_EposDisto_pReso_Omega = (TNtuple*)F_EposDisto_pReso_Omega->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Omega = T_EposDisto_pReso_Omega->GetEntries();
            T_EposDisto_pReso_Omega->SetBranchAddress("k_D",&k_D);
            T_EposDisto_pReso_Omega->SetBranchAddress("P1",&fP1);
            T_EposDisto_pReso_Omega->SetBranchAddress("P2",&fP2);
            T_EposDisto_pReso_Omega->SetBranchAddress("M1",&fM1);
            T_EposDisto_pReso_Omega->SetBranchAddress("M2",&fM2);
            T_EposDisto_pReso_Omega->SetBranchAddress("Tau1",&Tau1);
            T_EposDisto_pReso_Omega->SetBranchAddress("Tau2",&Tau2);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP1",&AngleRcP1);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP2",&AngleRcP2);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleP1P2",&AngleP1P2);
            for(unsigned uEntry=0; uEntry<N_EposDisto_pReso_Omega; uEntry++){
                T_EposDisto_pReso_Omega->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if(SourceVar%100==2){
                    fM1 = 1362;
                }
                if(k_D>k_CutOff) continue;
                RanVal1 = RanGen.Exponential(fM1/(fP1*Tau1));
                CleverMcLevyResoTM[3].AddBGT_RP(RanVal1,cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Omega;
        }

        CleverMcLevyResoTM[3].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[3], 2);
        Kitty.SetAnaSource(0,1.0);
        Kitty.SetAnaSource(1,2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if(SOURCE=="EPOS"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if(SOURCE=="EPOSrescaled"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if(SOURCE=="Levy_mT_Reso"){
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n",SOURCE.Data());
        goto CLEAN_SetUpCats_pOmegam;
    }

    if(POT.Contains("pOmega_Lattice")){
        int FLAG = 12;
        if(POT.EndsWith("_11")) FLAG = 11;
        if(POT.EndsWith("_12")) FLAG = 12;
        if(POT.EndsWith("_13")) FLAG = 13;
        if(POT.EndsWith("_14")) FLAG = 14;
        if(POT.EndsWith("_121")) FLAG = 121;
        if(POT.EndsWith("_122")) FLAG = 122;
        if(POT.EndsWith("_123")) FLAG = 123;
        //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars3S1[3]={100000,0.8,0.001};
        double PotPars5S2[9]={pOmega_Lattice,double(FLAG),0,0,0,0,0,0,0};
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential,3,true); cPotPars3S1->SetParameters(PotPars3S1);
        cPotPars5S2 = new CATSparameters(CATSparameters::tPotential,9,true); cPotPars5S2->SetParameters(PotPars5S2);
    }
    else{
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n",POT.Data());
        goto CLEAN_SetUpCats_pOmegam;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3334);
    Kitty.SetRedMass( (Mass_p*MassOmega)/(Mass_p+MassOmega) );

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);

    Kitty.SetSpin(0,1);
    Kitty.SetSpin(1,2);

    Kitty.SetChannelWeight(0, 3./8.);
    Kitty.SetChannelWeight(1, 5./8.);

    if(cPotPars3S1) Kitty.SetShortRangePotential(0,0,RepulsiveCore,*cPotPars3S1);
    if(cPotPars5S2) Kitty.SetShortRangePotential(1,0,fDlmPot,*cPotPars5S2);

    CLEAN_SetUpCats_pOmegam: ;
    if(cPars){delete cPars; cPars=NULL;}
    if(pPars){delete pPars; pPars=NULL;}
    //if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if(cPotPars3S1){delete cPotPars3S1; cPotPars3S1=NULL;}
    if(cPotPars5S2){delete cPotPars5S2; cPotPars5S2=NULL;}
}

DLM_Ck* DLM_CommonAnaFunctions::SetUpLednicky_pL(const unsigned& NumMomBins, const double* MomBins,  const TString& POT){
    double p0,p1,p2,p3;
    if(POT=="Lednicky_ND"){
        p0=1.77,p2=2.06,p1=3.78,p3=3.18;
    }
    else if(POT=="Lednicky_NF"){
        p0=2.18,p2=1.93,p1=3.19,p3=3.358;
    }
    else if(POT=="Lednicky_NSC89"){
        p0=2.73,p2=1.48,p1=2.87,p3=3.04;
    }
    else if(POT=="Lednicky_NSC97a"){
        p0=0.71,p2=2.18,p1=5.86,p3=2.76;
    }
    else if(POT=="Lednicky_NSC97b"){
        p0=0.9,p2=2.13,p1=4.92,p3=2.84;
    }
    else if(POT=="Lednicky_NSC97c"){
        p0=1.2,p2=2.08,p1=4.11,p3=2.92;
    }
    else if(POT=="Lednicky_NSC97d"){
        p0=1.71,p2=1.95,p1=3.46,p3=3.08;
    }
    else if(POT=="Lednicky_NSC97e"){
        p0=2.1,p2=1.86,p1=3.19,p3=3.19;
    }
    else if(POT=="Lednicky_NSC97f"){
        p0=2.51,p2=1.75,p1=3.03,p3=3.32;
    }
    else if(POT=="Lednicky_ESC08"){
        p0=2.7,p2=1.65,p1=2.97,p3=3.63;
    }
    else if(POT=="Lednicky_XeftLO"){
        p0=1.91,p2=1.23,p1=1.4,p3=2.13;
    }
    else if(POT=="Lednicky_XeftNLO"){
        p0=2.91,p2=1.54,p1=2.78,p3=2.72;
    }
    else if(POT=="Lednicky_JulichA"){
        p0=1.56,p2=1.59,p1=1.43,p3=3.16;
    }
    else if(POT=="Lednicky_JulichJ04"){
        p0=2.56,p2=1.66,p1=2.75,p3=2.93;
    }
    else if(POT=="Lednicky_JulichJ04c"){
        p0=2.66,p2=1.57,p1=2.67,p3=3.08;
    }
    else{
        printf("\033[1;31mERROR (SetUpLednicky_pL):\033[0m The pΛ potential '%s' does not exist\n",POT.Data());
    }
    DLM_Ck* DlmCk = new DLM_Ck(1,4,NumMomBins,MomBins,Lednicky_SingletTriplet);//
	DlmCk->SetPotPar(0,p0);
	DlmCk->SetPotPar(1,p1);
	DlmCk->SetPotPar(2,p2);
	DlmCk->SetPotPar(3,p3);
    return DlmCk;
}

void DLM_CommonAnaFunctions::SetUpBinning_pp(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion, const int& MomBinVar, const int& FitRegVar){
    double kMin;
    double kMax;
    double kStep;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        if(MomBinVar==0){
            kMin=0;
            kStep=4;
            kMax=376;//(i.e. max=376 MeV)
        }
        else if(MomBinVar==1){
            kMin=0;
            kStep=4;
            kMax=352;
        }
        else if(MomBinVar==2){
            kMin=0;
            kStep=4;
            kMax=400;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        if(MomBinVar==0){
            kMin=4;
            kStep=4;
            kMax=376;//(i.e. max=376 MeV)
        }
        else if(MomBinVar==1){
            kMin=4;
            kStep=4;
            kMax=352;
        }
        else if(MomBinVar==2){
            kMin=4;
            kStep=4;
            kMax=400;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        if(MomBinVar==0){
            kMin=4;
            kStep=4;
            kMax=376;//(i.e. max=376 MeV)
        }
        else if(MomBinVar==1){
            kMin=4;
            kStep=4;
            kMax=352;
        }
        else if(MomBinVar==2){
            kMin=4;
            kStep=4;
            kMax=400;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        if(MomBinVar==0){
            kMin=4;
            kStep=4;
            kMax=376;//(i.e. max=376 MeV)
        }
        else if(MomBinVar==1){
            kMin=4;
            kStep=4;
            kMax=352;
        }
        else if(MomBinVar==2){
            kMin=4;
            kStep=4;
            kMax=400;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }
    }
    else{
        printf("\033[1;31mERROR (SetUpBinning_pp):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        NumMomBins=0;
        return;
    }
    NumMomBins = floor((kMax-kMin)/(kStep));

    if(MomBins) delete [] MomBins;
    MomBins = new double [NumMomBins+1];
    MomBins[0] = kMin;
    for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
        MomBins[uBin] = MomBins[uBin-1]+kStep;
    }
    if(FitRegion) delete [] FitRegion;
    FitRegion = new double [4];
    FitRegion[0] = MomBins[0];
    FitRegion[1] = MomBins[NumMomBins];
    FitRegion[2] = MomBins[NumMomBins]+kStep;
    FitRegion[3] = MomBins[NumMomBins]+kStep*31.;//till 500

    //if(FitRegVar==0){
    //    FitRegion[0] = MomBins[0];
    //    FitRegion[1] = MomBins[NumMomBins];
    //    FitRegion[2] = MomBins[NumMomBins];
    //    FitRegion[3] = MomBins[NumMomBins];
    //}
    if(FitRegVar==0){
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 576;
    }
    else if(FitRegVar==1){
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 552;
    }
    else if(FitRegVar==2){
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 600;
    }
    else{
        printf("\033[1;31mERROR:\033[0m The FitRegVar '%i' does not exist\n",FitRegVar);
        return;
    }


}
void DLM_CommonAnaFunctions::SetUpBinning_pL(const TString& DataSample, unsigned& NumMomBins, double*& MomBins, double*& FitRegion,
                                             const int& MomBinVar, const int& FitRegVar){

    double kMin;
    double kFineMin;
    double kFineMax;
    double kMax;
    double kCoarseStep;
    double kFineStep;

    if(DataSample=="pp13TeV_MB_Run2paper"){
        if(MomBinVar==0){
            kMin=0;
            kFineMin=336;//272//216
            kFineMax=336;//304
            kMax=336;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==1){
            kMin=0;
            kFineMin=312;//272//216
            kFineMax=312;//304
            kMax=312;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==2){
            kMin=0;
            kFineMin=348;//272//216
            kFineMax=348;//304
            kMax=348;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;//348
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*20.;//588
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        if(MomBinVar==0){
            kMin=0;
            kFineMin=336;//272//216
            kFineMax=336;//304
            kMax=336;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==1){
            kMin=0;
            kFineMin=312;//272//216
            kFineMax=312;//304
            kMax=312;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==2){
            kMin=0;
            kFineMin=348;//272//216
            kFineMax=348;//304
            kMax=348;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==10){
            kMin=0;
            kFineMin=204;//272//216
            kFineMax=204;//304
            kMax=204;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==11){
            kMin=0;
            kFineMin=180;//272//216
            kFineMax=180;//304
            kMax=180;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else if(MomBinVar==12){
            kMin=0;
            kFineMin=228;//272//216
            kFineMax=228;//304
            kMax=228;//336
            kCoarseStep=12;
            kFineStep=12;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n",MomBinVar);
            return;
        }

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];

        if(FitRegVar==0){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 576;
        }
        else if(FitRegVar==1){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 552;
        }
        else if(FitRegVar==2){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 600;
        }
        else if(FitRegVar==10){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 360;//420
            FitRegion[3] = 576;
        }
        else if(FitRegVar==11){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 336;//432
            FitRegion[3] = 552;
        }
        else if(FitRegVar==12){
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 384;//408
            FitRegion[3] = 600;
        }
        //for the Schleching fits
        else if(FitRegVar==50){//only femto range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = MomBins[NumMomBins];
        }
        //only used for pol3
        else if(FitRegVar==51){//extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 456;
        }
        else if(FitRegVar==52){//extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 432;
        }
        else if(FitRegVar==53){//extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins];//336
            FitRegion[3] = 480;
        }
        else{
            printf("\033[1;31mERROR:\033[0m The FitRegVar '%i' does not exist\n",FitRegVar);
            return;
        }
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        kMin=0;
        kFineMin=336;
        kFineMax=336;
        kMax=336;
        kCoarseStep=12;
        kFineStep=12;

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;//348
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*20.;//588
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        kMin=0;
        kFineMin=336;
        kFineMax=336;
        kMax=336;
        kCoarseStep=12;
        kFineStep=12;

        //the number of coarse bins below kFineMin
        //floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        //and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin-kMin)/kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax-double(NumCoarseBinsBelow)*kCoarseStep)/kFineStep);
        //we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax-double(NumCoarseBinsBelow)*kCoarseStep-double(NumFineBins)*kFineStep)/kCoarseStep);

        NumMomBins=NumCoarseBinsBelow+NumFineBins+NumCoarseBinsAbove;

        if(MomBins) delete [] MomBins;
        MomBins = new double [NumMomBins+1];
        MomBins[0] = kMin;
        for(unsigned uBin=1; uBin<=NumMomBins; uBin++){
            if(uBin<=NumCoarseBinsBelow||uBin>NumCoarseBinsBelow+NumFineBins){
                MomBins[uBin] = MomBins[uBin-1]+kCoarseStep;
            }
            else{
                MomBins[uBin] = MomBins[uBin-1]+kFineStep;
            }
        }
        if(FitRegion) delete [] FitRegion;
        FitRegion = new double [4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins]+kCoarseStep;//348
        FitRegion[3] = MomBins[NumMomBins]+kCoarseStep*18.;//564
    }
    else{
        printf("\033[1;31mERROR (SetUpBinning_pL):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        NumMomBins=0;
        return;
    }

}


void DLM_CommonAnaFunctions::GetPurities_p(const TString& DataSample, const int& Variation, double* Purities){
    double PurityProton;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        PurityProton = 0.989859;
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        PurityProton = 0.9943;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityProton = 0.984265;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        PurityProton = 0.984265;
    }
    else if(DataSample=="pp13TeV_HM_BernieSource"){
        PurityProton = 0.9943;
    }
    else{
        printf("\033[1;31mERROR (GetPurities_p):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityProton = 1.0;
    }

    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    //const unsigned NumChannels_p = 4;
    //if(Purities){delete [] Purities; Purities = new double [NumChannels_p];}
    Purities[0] = PurityProton;
    Purities[1] = PurityProton;
    Purities[2] = PurityProton;
    Purities[3] = 1.-PurityProton;
}


void DLM_CommonAnaFunctions::GetPurities_L(const TString& DataSample, const int& Variation, double* Purities){
    double PurityLambda;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        PurityLambda = 0.96768;
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"){
        //printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 is not available yet!\n");
        if(Variation==0) PurityLambda = 0.9595;//the original value for the preliminaries
        else if(Variation==1) PurityLambda = 0.936;//spline fits 4th June 2020
        else if(Variation==2) PurityLambda = 0.936-0.006;//with uncertainties
        else if(Variation==3) PurityLambda = 0.936+0.006;//with uncertainties
        else if(Variation==-1) PurityLambda = 1.0;//use for SB corrected correlations
        else PurityLambda = 0.9595;
    }
    else if(DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        PurityLambda = 1.0;//use for SB corrected correlations
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityLambda = 0.937761;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        PurityLambda = 0.937761;
    }
    else{
        printf("\033[1;31mERROR (GetPurities_L):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityLambda = 1.0;
    }

    //for the Lambda:
    //0 = primary
    //1 = from Sigma0
    //2 = from Xim
    //3 = from Xi0
    //4 = missidentified
    Purities[0] = PurityLambda;
    Purities[1] = PurityLambda;
    Purities[2] = PurityLambda;
    Purities[3] = PurityLambda;
    Purities[4] = 1.-PurityLambda;
}

void DLM_CommonAnaFunctions::GetPurities_Xim(const TString& DataSample, const int& Variation, double* Purities){
    double PurityXim;
    if(DataSample=="pp13TeV_MB_Run2paper"){
        PurityXim = 0.956;
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        PurityXim = 0.956;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        PurityXim = 0.88;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        PurityXim = 0.88;
    }
    else{
        printf("\033[1;31mERROR (GetPurities_Xim):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        PurityXim = 1.0;
    }

    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Xi0(1530)
    //3 = from Omega
    //4 = missidentified
    Purities[0] = PurityXim;
    Purities[1] = PurityXim;
    Purities[2] = PurityXim;
    Purities[3] = PurityXim;
    Purities[4] = 1.-PurityXim;
}

//the last digit of Variation is generic.
//the second to-last digit referes to BernieSource
void DLM_CommonAnaFunctions::GetFractions_p(const TString& DataSample, const int& Variation, double* Fractions){
    //following my lambda pars with the 3 possible modifications
    //for the proton:
    //0 = primary
    //1 = from Lambda
    //2 = other feeddown (flat)
    //3 = missidentified
    double Modify_pp=1;
    switch(Variation%10){
        case 1 : Modify_pp=0.8; break;
        case 2 : Modify_pp=1.2; break;
        default : Modify_pp=1; break;
    }
    double pp_f0;//primary protons
    double pp_f1;//fraction of Lambda
    if(DataSample=="pp13TeV_MB_Run2paper"){
        pp_f0 = 0.87397;
        pp_f1 = 0.0882211;
    }
    else if(DataSample=="pp13TeV_HM_March19"){
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 CROSS CHECK pp_f0 and pp_f1!\n");
        pp_f0 = 0.873;
        pp_f1 = 0.0898;
    }
    else if(DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        pp_f0 = 0.823;
        pp_f1 = 0.125;
    }
    //Variation reflects the mT bin
    else if("pp13TeV_HM_BernieSource"){
      std::vector<double> pp_primary = { 0.82, 0.81, 0.81, 0.81, 0.81, 0.82, 0.83 };
      pp_f0 = pp_primary.at((Variation/10)%10);
      pp_f1 = 0.7*(1.-pp_f0);
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        pp_f0 = 0.862814;
        pp_f1 = 0.09603;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        pp_f0 = 0.862814;
        pp_f1 = 0.09603;
    }
    else{
        printf("\033[1;31mERROR (GetFractions_p):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        pp_f0 = 1.0;
        pp_f1 = 0.0;
    }
    //ratio between feed-down lambdas and all feed downs, it is by default 0.7, i.e. 70% of the feed-down into protons is from lambdas
    double arrayPercLamProton = pp_f1/(1.-pp_f0)*Modify_pp;
    Fractions[0] = pp_f0;
    Fractions[1] = (1.-pp_f0)*(arrayPercLamProton);
    Fractions[2] = (1.-pp_f0)*(1.-arrayPercLamProton);
    Fractions[3] = 1.;
}
//Variation -> use the two digits
//first digit->Modify_SigL
//second digit->Modify_XiL+3*Modify_PrimFrac
//0 -> default; 1 = -20%; 2 = +20%
void DLM_CommonAnaFunctions::GetFractions_L(const TString& DataSample, const int& Variation, double* Fractions){
    double Modify_SigL=1;
    double Modify_XiL=1;
    //the amount of prim lambdas depends on pT, at low pT (and k*)
    //the fraction could be lower. You can modify it here. This is by how much the PrimLambda fractional yield is reduced,
    //compared to the value at the average pT
    double Modify_PrimFrac=1;
    switch(Variation%10){
//the new values will be 0.6,0.84,1.2,1.44
//to keep consistency maybe just use 0.6,0.8,1.0,1.2,1.4
//this will correspond to s/l ratio of 0.2,0.267,0.333,0.4,0.467
        case 0 : Modify_SigL=1; break;
        case 1 : Modify_SigL=0.8;break;
        case 2 : Modify_SigL=1.2; break;
        case 3 : Modify_SigL=0.6; break;
        case 4 : Modify_SigL=1.4; break;
        default : Modify_SigL=1; break;
    }
    switch((Variation/10)%3){
        case 0 : Modify_XiL=1; break;
        case 1 : Modify_XiL=0.8;break;
        case 2 : Modify_XiL=1.2; break;
        default : Modify_XiL=1; break;
    }
    //this only works assuming zero material
    //do not use for newer interations!
    switch((Variation/10)/3){
        case 0 : Modify_PrimFrac=1; break;
        case 1 : Modify_PrimFrac=0.95;break;
        default : Modify_PrimFrac=1; break;
    }
    double pL_f0;//fraction of primary Lambdas
    double pL_f1;//fraction of Sigma0
    double pL_f2;//fractions of Xi0/m (each)
    double pL_fm;//fraction of material
    if(DataSample=="pp13TeV_MB_Run2paper"){
        pL_f0 = 0.601008;
        pL_f1 = 0.200336;
        pL_f2 = 0.099328;
        pL_fm = 0;
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        pL_f0 = 0.576066;
        pL_f1 = 0.192022;
        pL_f2 = 0.115956;
        pL_fm = 0.0;
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        pL_f0 = 0.521433;
        pL_f1 = 0.173811;
        pL_f2 = 0.152378;
        pL_fm = 0.0;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        pL_f0 = 0.521433;
        pL_f1 = 0.173811;
        pL_f2 = 0.152378;
        pL_fm = 0.0;
    }
    else{
        printf("\033[1;31mERROR (GetFractions_L):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        pL_f0 = 1.0;
        pL_f1 = 0.0;
        pL_f2 = 0.0;
        pL_fm = 0.0;
    }

    //this only works assuming zero material
    //do not use for newer interations!
    if(pL_f2){
        pL_f0 *= Modify_PrimFrac;
        pL_f1 *= Modify_PrimFrac;
        pL_f2 = (1.-pL_f0-pL_f1)*0.5;
    }

    double SigLambdaPrimDir = pL_f0+pL_f1;
    //ratio between sigma0 feed down and primary lambdas. By default this should be 1:3
    double arrayPercSigLambda=pL_f1/pL_f0*Modify_SigL;
    //ration between xim feed down to the flat (xi0) feed. By default we assume it is 0.5
    double arrayPercXiLambda=pL_f2/(1.-pL_f0-pL_f1)*Modify_XiL;
    double FracOfLambda = 1./(1.+arrayPercSigLambda);
    //0 is primary
    //1 is from Sigma0
    //2 is is from Xim
    //3 is is the flat feeddown (it should be just xi0, i.e. == xim)
    //  remark: up to Aug 2020, the material contribution, c.a. 1% was distributed amount ALL other
    //  entries, primary, feed etc. Actually it would make most sence to put it as an additional flat contribution!
    //4 is for the missid
    Fractions[0] = SigLambdaPrimDir*FracOfLambda;
    Fractions[1] = SigLambdaPrimDir*(1.-FracOfLambda);
    Fractions[2] = (1.-SigLambdaPrimDir)*(arrayPercXiLambda);
    Fractions[3] = 1.-Fractions[0]-Fractions[1]-Fractions[2];
    Fractions[4] = 1.;
}
void DLM_CommonAnaFunctions::GetFractions_Xim(const TString& DataSample, const int& Variation, double* Fractions){
    //0 = primary
    //1 = from Xi-(1530)
    //2 = from Omega
    //3 = flat
    //4 = missidentified

    //ratio Xi-(1530) to Xi-
    const double Xim1530_to_Xim = 0.32*(1./3.);
    //ratio Xi0(1530) to Xi0 (n=neutral)
    //const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;
    const double OmegamXim_BR = 0.086;
    //the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
    //hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
    Fractions[0] = 1.-3.*Xim1530_to_Xim-Omegam_to_Xim*OmegamXim_BR;
    Fractions[1] = Xim1530_to_Xim;
    Fractions[2] = Omegam_to_Xim*OmegamXim_BR;
    Fractions[3] = 1.-Fractions[2]-Fractions[1]-Fractions[0];
    Fractions[4] = 1.;
}
//0 is primary
//1 is pL->pp
//2 is flat feed
//3 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pp(const TString& DataSample, const int& Variation_p, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    lambda_pars[0] = Purities_p[0]*Fraction_p[0]*Purities_p[0]*Fraction_p[0];
    lambda_pars[1] = Purities_p[0]*Fraction_p[0]*Purities_p[1]*Fraction_p[1]*2.;
    lambda_pars[3] = (Purities_p[0]+Purities_p[0]+Purities_p[3])*Purities_p[3];
    lambda_pars[2] = 1.-lambda_pars[3]-lambda_pars[1]-lambda_pars[0];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<4; uLam++){
    //    printf("λ(pp)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}

//0 is primary
//1 is pSigma0->pL
//2 is pXim->pL
//3 is the flat feeddown
//4 is missid
//Variation_L 4 digits
//the last two digits are passed to the variations of the fraction (last digit is for Lambda, the second to last for Xi)
//the first two digits are passed to the purities.
void DLM_CommonAnaFunctions::SetUpLambdaPars_pL(const TString& DataSample, const int& Variation_p, const int& Variation_L, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_L[5];
    double Fraction_L[5];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    GetPurities_L(DataSample,Variation_L/100,Purities_L);
    GetFractions_L(DataSample,Variation_L%100,Fraction_L);
    lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_L[0]*Fraction_L[0];
    lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_L[1]*Fraction_L[1];
    lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_L[2]*Fraction_L[2];
    lambda_pars[4] =    Purities_p[0]*Purities_L[4]+Purities_p[3]*Purities_L[0]+Purities_p[3]*Purities_L[4];
    lambda_pars[3] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[4];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<5; uLam++){
    //    printf("λ(pΛ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}
//0 is primary
//1 is from Xim1530
//2 is from Xin1530
//3
//4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pXim(const TString& DataSample, const int& Variation_p, const int& Variation_Xim, double* lambda_pars){
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_Xim[5];
    double Fraction_Xim[5];
    GetPurities_p(DataSample,Variation_p,Purities_p);
    GetFractions_p(DataSample,Variation_p,Fraction_p);
    GetPurities_Xim(DataSample,Variation_Xim,Purities_Xim);
    GetFractions_Xim(DataSample,Variation_Xim,Fraction_Xim);
    lambda_pars[0] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[0]*Fraction_Xim[0];
    lambda_pars[1] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[1]*Fraction_Xim[1];
    lambda_pars[2] =    Purities_p[0]*Fraction_p[0]*Purities_Xim[2]*Fraction_Xim[2];
    lambda_pars[4] =    Purities_p[0]*Purities_Xim[4]+Purities_p[3]*Purities_Xim[0]+Purities_p[3]*Purities_Xim[4];
    lambda_pars[3] =    1.-lambda_pars[0]-lambda_pars[1]-lambda_pars[2]-lambda_pars[4];

    //double SUM=0;
    //for(unsigned uLam=0; uLam<5; uLam++){
    //    printf("λ(pΞ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //    SUM+=lambda_pars[uLam]*100.;
    //}
    //printf("SUM: %.1f\n------------\n",SUM);
}

//        DataPeriod=="pp13TeV"?  :
//                                ;

TH2F* DLM_CommonAnaFunctions::GetResolutionMatrix(const TString& DataSample,const TString&& System){
    TString FileName;
    TString HistoName;


    if(DataSample=="pp13TeV_MB_Run2paper"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/ALICE_pp_13TeV.root";
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"||DataSample=="pp13TeV_HM_RotPhiDec19"
            ||DataSample=="pp13TeV_HM_DimiJun20"||DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        //N.B. special rule for pp and pLambda below
        FileName = CatsFilesFolder[0]+"/MomentumSmear/ALICE_pp_13TeV.root";
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/Sample3_MeV_compact.root";;
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/ResolutionMatrices/Sample3_MeV_compact.root";
    }
    else{
        printf("\033[1;31mERROR (GetResolutionMatrix):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    if(System=="pp"){
        HistoName = "hSigmaMeV_Proton_Proton";
    }
    else if(System=="pLambda"){
        //N.B. special rule for pp and pLambda below
        HistoName = "hSigmaMeV_Proton_Lambda";
    }
    else if(System=="LambdaLambda"){
        HistoName = "hSigmaMeV_Lambda_Lambda";
    }
    else if(System=="pXim"){
        HistoName = "hSigmaMeV_Proton_Xim";
    }
    else{
        printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
    }



    if(DataSample=="pp13TeV_HM_DimiJun20"&&System=="pp"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/ALICE_pp_13TeV_MEpp.root";
        HistoName = "h_RESO_pp_MeV";
    }
    if(DataSample=="pp13TeV_HM_DimiJun20"&&System=="pLambda"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/ALICE_pp_13TeV_MEpL.root";
        HistoName = "h_RESO_pL_MeV";
    }
    if(DataSample=="pp13TeV_HM_BernieSource"&&System=="pp"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/Sample6_MeV_compact.root";
        HistoName = "hSigmaMeV_Proton_Proton";
    }
    if(DataSample=="pp13TeV_HM_BernieSource"&&System=="pLambda"){
        FileName = CatsFilesFolder[0]+"/MomentumSmear/Sample6_MeV_compact.root";
        HistoName = "hSigmaMeV_Proton_Lambda";
    }

    //this is the unfolded data
    if(DataSample=="pp13TeV_HM_DimiJul20"||DataSample=="pp13TeV_HM_DimiMay21"||DataSample=="pp13TeV_HM_DimiJun21"){
        return NULL;
    }

    ///FUCKING ROOT SUCKS!!!!!! SUCK MY COCK!!!! SUCK IT YOU BITCH!!!!!!!!!!
    //so we need to copy our histogram, as else we lose it when we delete the file
    //and we need to change to the "central" root directory, as else histoCopy will also be lost
    //and we need to play with the name a little bit, else we are fucked!
    TFile* FileROOT = new TFile(FileName, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",HistoName.Data(),FileName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}
TH2F* DLM_CommonAnaFunctions::GetResidualMatrix(const TString&& FinalSystem, const TString& InitialSystem){
    TString FileName;
    TString HistoName;

    FileName = CatsFilesFolder[0]+"/DecaySmear/run2_decay_matrices_old.root";

    if(FinalSystem=="pp"&&InitialSystem=="pLambda"){
        HistoName = "hRes_pp_pL";
    }
    else if(FinalSystem=="pLambda"&&InitialSystem=="pSigma0"){
        HistoName = "hRes_pL_pSigma0";
    }
    else if(FinalSystem=="pLambda"&&InitialSystem=="pXim"){
        HistoName = "hRes_pL_pXim";
    }
    else if(FinalSystem=="pLambda"&&InitialSystem=="pXi0"){
        FileName = CatsFilesFolder[0]+"/DecaySmear/pXi0_pL.root";
        HistoName = "pXi0_pL";
    }
    else if(FinalSystem=="pXim"&&InitialSystem=="pXim1530"){
        HistoName = "hRes_pXim_pXim1530";
    }
    else{
        printf("\033[1;31mERROR:\033[0m The decay '%s->%s' does not exist\n",InitialSystem.Data(),FinalSystem.Data());
    }
    TFile* FileROOT = new TFile(FileName, "read");
    TH2F* histo = (TH2F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",HistoName.Data(),FileName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F*)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

//mTbin == -1 means we take the integrated function
//iReb = 0 is 4 MeV, 1 is 8, 2 is 12, 3 is 16, 4 is 20
TH1F* DLM_CommonAnaFunctions::GetAliceExpCorrFun(const TString& DataSample,const TString& System, const TString& CutVar, const int& iReb, const bool& AddSyst, const int mTbin){
    TString FileName;
    TString HistoName;
    TString SystFileName="";
    TString SystHistName="";
    TGraph* gRelSyst=NULL;

    if(DataSample=="pp13TeV_MB_Run2paper"){
        if(System=="pp"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pp%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="LambdaLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_LL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pXim"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pXi%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pp13TeV_HM_March19"||DataSample=="pp13TeV_HM_Dec19"){
        if(System=="pp"){
            if(mTbin==-1){
                //buggy for all but _0
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pp%s.root",CutVar.Data());
                if(CutVar=="_0") FileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample10HM_ver2/CFOutput_pp.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else if(mTbin>=0&&mTbin<7){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutputALL_mT_pp_HM%s.root",CutVar.Data());
                //ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if(CutVar=="_0") FileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample10HM_ver2/CFOutputALL_mT_pp_HM.root";
                HistoName = TString::Format("hCkTotNormWeightMeV_mTBin_%i",mTbin);
                SystFileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pp.root";
                SystHistName = "SystErrRel";
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pL%s.root",CutVar.Data());
                //ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if(CutVar=="_0") FileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample10HM/CFOutput_pL.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
                SystFileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pL.root";
                SystHistName = "SystErrRel";
            }
            else if(mTbin>=0&&mTbin<6){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutputALL_mT_pL_HM%s.root",CutVar.Data());
                //ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if(CutVar=="_0") FileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample10HM/CFOutputALL_mT_pL_HM.root";
                HistoName = TString::Format("hCk_RebinnedMeV_%i_mTBin_%i",0,mTbin);
                SystFileName = CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pL.root";
                SystHistName = "SystErrRel";
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="LambdaLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_LL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pXim"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pXi%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    //the correlations are obtained with Dimi's code, normalization to the total yield
    //used for the pLambda paper proposal in June 2020, sideband corrected (folded)
    else if(DataSample=="pp13TeV_HM_DimiJun20"){
        if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/CkSB_pL_%s.root",CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV",float(iReb+1)*4.);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pp13TeV_HM_DimiJul20"){
        if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/Unfolded/090720/CkSB_pL_%s.root",CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV",float(iReb+1)*4.);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pp13TeV_HM_DimiMay21"){
        if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/DimiMay21/UnfoldedNorm240_340/CkSB_pL_%s.root",CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV",float(iReb+1)*4.);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pp13TeV_HM_DimiJun21"){
        if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/DimiJun21/UnfoldedNorm240_340/CkSB_pL_%s.root",CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV",float(iReb+1)*4.);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pp13TeV_HM_RotPhiDec19"){
        if(System=="pp"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pp_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pL_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="LambdaLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_LL_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pXim"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pXi_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    //the data used by Bernie for the source paper. Only pp and pLambda
    else if(DataSample=="pp13TeV_HM_BernieSource"){
        if(System=="pp"){
            if(mTbin>=0&&mTbin<7){
                FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar%s_HM_%i.root",mTbin+1,CutVar.Data(),mTbin);
                HistoName = TString::Format("hCk_RebinnedppVar%sMeV_0",CutVar.Data());
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin>=0&&mTbin<6){
              FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar%s_HM_%i.root",mTbin+1,CutVar.Data(),mTbin);
              HistoName = TString::Format("hCk_RebinnedpLVar%sMeV_0",CutVar.Data());
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pPb5TeV_Run2paper"){
        if(System=="pp"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pp%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="LambdaLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_LL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pXim"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pXi%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else if(DataSample=="pPb5TeV_CPR_Mar19"){
        if(System=="pp"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pp%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="LambdaLambda"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_LL%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else if(System=="pXim"){
            if(mTbin==-1){
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pXi%s.root",CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i",iReb);
            }
            else{
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n",mTbin,DataSample.Data(),System.Data());
            }
        }
        else{
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n",System.Data());
        }
    }
    else{
        printf("\033[1;31mERROR (GetAliceExpCorrFun):\033[0m The data sample '%s' does not exist\n",DataSample.Data());
        FileName = "";
    }

    TFile* FileROOT = new TFile(FileName, "read");
    if(!FileROOT){printf("\033[1;31mERROR:\033[0m The file '%s' does not exist\n",FileName.Data());return NULL;}
    TH1F* histo = (TH1F*)FileROOT->Get(HistoName);
    if(!histo){printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n",HistoName.Data(),FileName.Data());return NULL;}
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F*)histo->Clone("histoCopy");
    delete FileROOT; FileROOT=NULL;
    histoCopy->SetName(Name);

    if(AddSyst){
        FileROOT = new TFile(SystFileName, "read");
        if(FileROOT){
            TH1F* hRelSyst = (TH1F*)FileROOT->Get(SystHistName);
            if(!hRelSyst){printf("\033[1;31mERROR:\033[0m The hRelSyst '%s' does not exist\n",SystHistName.Data());}
            else{
                double MaxMom = hRelSyst->GetXaxis()->GetBinUpEdge(hRelSyst->GetNbinsX());
                gRelSyst = new TGraph();
                gRelSyst->SetName("gRelSyst");
                for(int iBin=1; iBin<=hRelSyst->GetNbinsX(); iBin++){
                    gRelSyst->SetPoint(iBin-1,hRelSyst->GetBinCenter(iBin),hRelSyst->GetBinContent(iBin));
                }
                for(int iBin=1; iBin<=histoCopy->GetNbinsX(); iBin++){
                    if(histoCopy->GetBinCenter(iBin)>MaxMom) break;
                    double StatErr = histoCopy->GetBinError(iBin);
                    double SystErr = histoCopy->GetBinContent(iBin)*gRelSyst->Eval(histoCopy->GetBinCenter(iBin));
                    double TotErr = sqrt(StatErr*StatErr+SystErr*SystErr);
                    histoCopy->SetBinError(iBin,TotErr);
                }
            }
        }
        if(FileROOT) {delete FileROOT; FileROOT=NULL;}
    }

    return histoCopy;
}

DLM_CleverMcLevyReso* DLM_CommonAnaFunctions::GetCleverMcLevyReso_pp(){
    return &CleverMcLevyReso[0];
}

DLM_CleverMcLevyReso* DLM_CommonAnaFunctions::GetCleverMcLevyReso_pL(){
    return &CleverMcLevyReso[1];
}

DLM_CleverMcLevyReso* DLM_CommonAnaFunctions::GetCleverMcLevyReso_pXim(){
    return &CleverMcLevyReso[2];
}

DLM_CleverMcLevyReso* DLM_CommonAnaFunctions::GetCleverMcLevyReso_pOmegam(){
    return &CleverMcLevyReso[3];
}

DLM_CleverMcLevyResoTM* DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pp(){
    return &CleverMcLevyResoTM[0];
}

DLM_CleverMcLevyResoTM* DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pL(){
    return &CleverMcLevyResoTM[1];
}

DLM_CleverMcLevyResoTM* DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pXim(){
    return &CleverMcLevyResoTM[2];
}
DLM_CleverMcLevyResoTM* DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pOmegam(){
    return &CleverMcLevyResoTM[3];
}
DLM_CleverMcLevyResoTM* DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pipi(){
    return &CleverMcLevyResoTM[4];
}
void DLM_CommonAnaFunctions::SetCatsFilesFolder(const TString& folder){
    CatsFilesFolder[0] = folder;
}

DLM_Histo<double>* ConvertThetaAngleHisto(const TString& FileName, const TString& HistoName, const double kMin, const double kMax){
    TFile* InputFile = new TFile(FileName, "read");
    TH2F* InputHisto = NULL;
    if(InputFile){
        InputHisto = (TH2F*)InputFile->Get(HistoName);
    }
    TH1D* Projection=NULL;
    if(InputHisto){
        Projection = InputHisto->ProjectionX("ConvertThetaAngleHisto",InputHisto->GetYaxis()->FindBin(kMin),InputHisto->GetYaxis()->FindBin(kMax));
    }

    const unsigned NumBins = Projection->GetNbinsX();
    Projection->Scale(1./Projection->Integral(1,NumBins));

    DLM_Histo<double> CummDistr;
    CummDistr.SetUp(1);
//printf("NumBins=%u; %f --> %f\n",NumBins,Projection->GetBinLowEdge(1),Projection->GetXaxis()->GetBinUpEdge(NumBins));
//usleep(4e6);
//printf("Calling SetUp\n");
    CummDistr.SetUp(0,NumBins,Projection->GetBinLowEdge(1),Projection->GetXaxis()->GetBinUpEdge(NumBins));
//printf("Ended SetUp\n");
//usleep(4e6);
//printf("Initialize...\n");
    CummDistr.Initialize();
//printf("Ended...\n");
    CummDistr.SetBinContent(unsigned(0),Projection->GetBinContent(1));
    for(unsigned uBin=1; uBin<NumBins; uBin++){
        CummDistr.SetBinContent(uBin,CummDistr.GetBinContent(uBin-1)+Projection->GetBinContent(uBin+1));
    }

    double* MomBins = new double [NumBins+1];
    MomBins[0] = 0;
    //MomBins[1] = CummDistr.GetBinContent(unsigned(0));
    //MomBins[NumBins] = 1;
    for(unsigned uBin=1; uBin<=NumBins; uBin++){
        MomBins[uBin] = (CummDistr.GetBinContent(uBin)+CummDistr.GetBinContent(uBin-1))*0.5;
        if(MomBins[uBin]<=MomBins[uBin-1]) MomBins[uBin] = MomBins[uBin-1]+1e-6;
        //printf("MomBins[%u] = %e\n",uBin,MomBins[uBin]);
    }
//usleep(2e6);
    DLM_Histo<double>* Result = new DLM_Histo<double>();
    Result->SetUp(1);
    Result->SetUp(0,NumBins,MomBins);
    Result->Initialize();
//printf("Result->Initialize\n");
//usleep(2e6);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        if(!Projection) break;
        Result->SetBinCenter(0,uBin,CummDistr.GetBinContent(uBin));//cum. value
        Result->SetBinContent(uBin,CummDistr.GetBinCenter(0,uBin));//momentum value
        //printf("%u: x=%.4f; y=%.4f;\n",uBin,CummDistr.GetBinContent(uBin),CummDistr.GetBinCenter(0,uBin));
    }

    delete [] MomBins;
    delete InputFile;
    return Result;
}


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
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, CATS* Kitty){
    DLM_Ck dlmck(1,0,*Kitty);
    RootFile_DlmCk(RootFileName,GraphName,&dlmck);
}

void RootFile_DlmSource(const TString& RootFileName, const TString& GraphName, CATS* Kitty, const unsigned& NumBins, const double& rMin, const double& rMax, const double& lambda, const bool& FourPi){
    TFile* RootFile = new TFile(RootFileName,"update");
    if(!RootFile) RootFile = new TFile(RootFileName,"recreate");
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    const double rWidth = (rMax-rMin)/double(NumBins);
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        double RAD = rWidth*0.5 + double(uBin)*rWidth;
        if(FourPi) graph.SetPoint(uBin,RAD,Kitty->EvaluateTheSource(10,RAD,0)*lambda);
        else graph.SetPoint(uBin,RAD,Kitty->EvaluateTheSource(10,RAD,0)*lambda/(4.*Pi*(RAD+1e-6)*(RAD+1e-6)));
    }
    graph.Write("",TObject::kOverwrite);
    delete RootFile;
}

/*
void DLM_CommonAnaFunctions::Clean_CommonAnaFunctions(){
    for(unsigned uLevy=0; uLevy<NumCleverLevyObjects; uLevy++){
        delete CleverLevy[uLevy];
        CleverLevy[uLevy] = NULL;
    }
    delete [] CleverLevy;
    CleverLevy=NULL;

    for(unsigned uLevy=0; uLevy<NumCleverLevyObjects; uLevy++){
        delete CleverMcLevyReso[uLevy];
        CleverMcLevyReso[uLevy] = NULL;
    }
    delete [] CleverMcLevyReso;
    CleverMcLevyReso=NULL;
}
*/
