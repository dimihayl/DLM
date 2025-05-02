
#include <cctype> // for tolower
#include <cstring> // for strlen

//#include "DLM_Histo.h"
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
#include "DLM_HistoAnalysis.h"
#include "CECA.h"
#include "TREPNI.h"
#include "DLM_RootWrapper.h"
#include "DLM_RootFit.h"
#include "DLM_CppTools.h"

#include "TString.h"
#include "TH2F.h"
#include "TFile.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TVector3.h"
#include "TF1.h"
#include "TH1F.h"


DLM_CommonAnaFunctions::DLM_CommonAnaFunctions() : NumCleverLevyObjects(9)
{
    // Simple_Reso = NULL;
    // Simple_Reso = new MS_GaussExp_mT_Simple [NumCleverLevyObjects];
    CleverLevy = NULL;
    CleverLevy = new DLM_CleverLevy[NumCleverLevyObjects];
    CleverMcLevyReso = NULL;
    CleverMcLevyReso = new DLM_CleverMcLevyReso[NumCleverLevyObjects];
    CleverMcLevyResoTM = NULL;
    CleverMcLevyResoTM = new DLM_CleverMcLevyResoTM[NumCleverLevyObjects];
    CatsFilesFolder = new TString();
}

DLM_CommonAnaFunctions::~DLM_CommonAnaFunctions()
{
    // if(Simple_Reso){delete[]Simple_Reso;Simple_Reso=NULL;}
    if (CleverLevy)
    {
        delete[] CleverLevy;
        CleverLevy = NULL;
    }
    if (CleverMcLevyReso)
    {
        delete[] CleverMcLevyReso;
        CleverMcLevyReso = NULL;
    }
    if (CleverMcLevyResoTM)
    {
        delete[] CleverMcLevyResoTM;
        CleverMcLevyResoTM = NULL;
    }
    delete CatsFilesFolder;
}

// POT:
//   "AV18", no pot vars so far
//   no sor var so far
void DLM_CommonAnaFunctions::SetUpCats_pp(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    CATSparameters *cPars = NULL;

    CATSparameters *cPotPars1S0 = NULL;
    CATSparameters *cPotPars1P1 = NULL;

    CATSparameters *cPotPars3S1 = NULL;
    CATSparameters *cPotPars3P0 = NULL;
    CATSparameters *cPotPars3P1 = NULL;
    CATSparameters *cPotPars3P2 = NULL;
    CATSparameters *cPotPars1D2 = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    //unsigned NumChannels = 0;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "NULL" || SOURCE == "")
    {
    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "DoubleGauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 3, true);
        cPars->SetParameter(0, 1.0);
        cPars->SetParameter(1, 2.0);
        cPars->SetParameter(2, 0.5);
        Kitty.SetAnaSource(DoubleGaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.6) * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.6 * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[0].InitStability(20, 1, 2);
        CleverLevy[0].InitScale(35, 0.25, 2.0);
        CleverLevy[0].InitRad(256, 0, 64);
        CleverLevy[0].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[0].InitStability(20, 1, 2);
        CleverLevy[0].InitScale(35, 0.25, 2.0);
        CleverLevy[0].InitRad(256, 0, 64);
        CleverLevy[0].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0, sqrt(1.6) * 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[0].InitStability(20, 1, 2);
        CleverLevy[0].InitScale(35, 0.25, 2.0);
        CleverLevy[0].InitRad(256, 0, 64);
        CleverLevy[0].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.6 * 1.2);
        Kitty.SetAnaSource(1, 1.6);
    }
    else if (SOURCE == "GaussExpTotSimple_2body")
    {
        // printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Gauss_mT_Reso)\n");
        cPars = new CATSparameters(CATSparameters::tSource, 11, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.65);        // tau
        cPars->SetParameter(2, 1. - 0.3578); // prim
        cPars->SetParameter(3, 1361.52);     // reso mass
        cPars->SetParameter(4, Mass_p);
        cPars->SetParameter(5, Mass_pic);
        cPars->SetParameter(6, 1.65);
        cPars->SetParameter(7, 1. - 0.3578);
        cPars->SetParameter(8, 1361.52);
        cPars->SetParameter(9, Mass_p);
        cPars->SetParameter(10, Mass_pic);
        Kitty.SetAnaSource(GaussExpTotSimple_2body, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        CleverMcLevyReso[0].InitStability(21, 1, 2);
        CleverMcLevyReso[0].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[0].InitRad(257, 0, 64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0, 1);
        CleverMcLevyReso[0].InitReso(1, 1);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[0].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else
        {
            CleverMcLevyReso[0].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }

        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            // DLM_Histo<double>* HISTO = ConvertThetaAngleHisto(
            //"/home/dmihaylov/Dudek_Ubuntu/MyApps/SomeTests/RandomCos/fCOS_180.root","hTheta2D",400,700);
            CleverMcLevyReso[0].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[0].SetUpResoEmission(1, 0, HISTO);
        }
        CleverMcLevyReso[0].InitNumMcIter(200000);

        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[0].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[0].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[0].InitRad(257 * 2, 0, 64);
        CleverMcLevyReso[0].InitType(2);
        CleverMcLevyReso[0].InitReso(0, 1);
        CleverMcLevyReso[0].InitReso(1, 1);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[0].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else
        {
            CleverMcLevyReso[0].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[0].SetUpReso(1, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }

        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto(
                "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 600);
            // DLM_Histo<double>* HISTO = ConvertThetaAngleHisto(
            //"/home/dmihaylov/Dudek_Ubuntu/MyApps/SomeTests/RandomCos/fCOS_Flat.root","hTheta2D",400,700);
            CleverMcLevyReso[0].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[0].SetUpResoEmission(1, 0, HISTO);
        }

        CleverMcLevyReso[0].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[0], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[0].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[0].InitStability(21, 1, 2);
        CleverMcLevyResoTM[0].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[0].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[0].InitType(2);
        CleverMcLevyResoTM[0].SetUpReso(0, 0.6422);
        CleverMcLevyResoTM[0].SetUpReso(1, 0.6422);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
            // printf("Hello 0\n");
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            CleverMcLevyResoTM[0].AddBGT_PR(490. / 1362. * 1.65, 1.);
            CleverMcLevyResoTM[0].AddBGT_RP(490. / 1362. * 1.65, -1.);
            CleverMcLevyResoTM[0].AddBGT_RR(490. / 1362. * 1.65, -1., 490. / 1362. * 1.65, 1., -1.);
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass, 4 is 3 body with fixed mass, 5 is 3 body with EPOS mass
        else
        {
            // printf("Hello 2\n");
            const double k_CutOff = fabs(int(int(SourceVar) / 10) * 10.);
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

            TFile *F_EposDisto_p_pReso;
            if (SourceVar % 100 == 4 || SourceVar % 100 == 5)
                F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/Only3body_p_pReso_3body.root");
            else
                F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_p_pReso.root");
            // printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/").Data());
            TNtuple *T_EposDisto_p_pReso = (TNtuple *)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_p_pReso; uEntry++)
            {
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                Tau2 = 1.65;
                if (fabs(SourceVar % 100) == 2)
                {
                    fM2 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                // wrong sign (e.g. SourceVar==202)
                if (SourceVar > 0)
                {
                    CleverMcLevyResoTM[0].AddBGT_PR(RanVal1, -cos(AngleRcP2));
                    CleverMcLevyResoTM[0].AddBGT_RP(RanVal1, cos(AngleRcP2));
                }
                // correct sign (e.g. SourceVar==-202)
                else
                {
                    CleverMcLevyResoTM[0].AddBGT_PR(RanVal1, cos(AngleRcP2));
                    CleverMcLevyResoTM[0].AddBGT_RP(RanVal1, -cos(AngleRcP2));
                }
            }
            delete F_EposDisto_p_pReso;

            TFile *F_EposDisto_pReso_pReso;
            if (SourceVar % 100 == 4 || SourceVar % 100 == 5)
                F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/Only3body_pReso_pReso_3body.root");
            else
                F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_pReso.root");
            // printf("F_EposDisto_p_pReso=%p (%s)\n",F_EposDisto_p_pReso,TString(CatsFilesFolder[0]+"/Source/EposAngularDist/Epos3body_pReso_pReso_3body.root").Data());
            TNtuple *T_EposDisto_pReso_pReso = (TNtuple *)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_pReso; uEntry++)
            {
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 1.65;
                if (fabs(SourceVar % 100) == 2)
                {
                    fM1 = 1362;
                    fM2 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[0].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[0].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[0].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[0], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        // printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        Kitty.SetUseAnalyticSource(false);
        // AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
        Kitty.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19");
        // AB_pp.SetInputFileName(THERMALFILE.Data());
        Kitty.SetMaxPairsToRead(64e6);
        Kitty.SetMaxPairsPerBin(16000 / 2);
        Kitty.SetMixingDepth(8);
        Kitty.SetThetaDependentSource(false);
        Kitty.SetTransportRenorm(1);
        Kitty.SetTauCorrection(false);
        // goto CLEAN_SetUpCats_pp;
    }
    else if (SOURCE == "EPOStheta")
    {
        // printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        Kitty.SetUseAnalyticSource(false);
        // AB_pp.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/Scratch9_OSCAR1997_AllLrzProtonFiles.f19");
        Kitty.SetInputFileName("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/EPOS_OUTPUT_FILES/EPOS_LBF_pp91/pp_nist0_030219.f19");
        // AB_pp.SetInputFileName(THERMALFILE.Data());
        Kitty.SetMaxPairsToRead(64e6);
        Kitty.SetMaxPairsPerBin(16000 / 2);
        Kitty.SetMixingDepth(8);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetTransportRenorm(1);
        Kitty.SetGridEpsilon(1. / 8192);
        Kitty.SetGridMaxDepth(8);
        Kitty.SetTauCorrection(false);
        // goto CLEAN_SetUpCats_pp;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pp;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pp;
    }

    if (POT == "AV18")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        int POT_FLAG = v18_Coupled3P2;
        if(PotVar==1) POT_FLAG = v18_SingleChannelMagic;
        double PotPars1S0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 0, 0};
        double PotPars3P0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 0};
        double PotPars3P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 1};
        double PotPars3P2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 1, 1, 2};
        double PotPars1D2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, 1, 0, 2, 2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if (POT == "AV18_pn")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        int POT_FLAG = v18_Coupled3P2;
        if(PotVar==1) POT_FLAG = v18_SingleChannelMagic;
        double PotPars1S0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 0, 0, 0};
        double PotPars1P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 0, 1, 1};
        //const int& Spin, const int& AngMom, const int& TotMom
        double PotPars3S1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 1, 0, 1};
        double PotPars3P0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 1, 1, 0};
        double PotPars3P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 1, 1, 1};
        double PotPars3P2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 1, 1, 2};
        double PotPars1D2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, 1, -1, 0, 2, 2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars1P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1P1->SetParameters(PotPars1P1);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if (POT == "AV18_nn")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        int POT_FLAG = v18_Coupled3P2;
        if(PotVar==1) POT_FLAG = v18_SingleChannelMagic;
        double PotPars1S0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 0, 0, 0};
        double PotPars1P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 0, 1, 1};
        //const int& Spin, const int& AngMom, const int& TotMom
        double PotPars3S1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 1, 0, 1};
        double PotPars3P0[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 1, 1, 0};
        double PotPars3P1[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 1, 1, 1};
        double PotPars3P2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 1, 1, 2};
        double PotPars1D2[8] = {NN_AV18, static_cast<double>(POT_FLAG), 1, -1, -1, 0, 2, 2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars1P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1P1->SetParameters(PotPars1P1);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if(POT == "pp_AV18_Toy"){
      double PotPars1S0[3] = {0,0,0};
      double PotPars3P0[3] = {1,0,0};
      double PotPars3P1[3] = {2,0,0};
      double PotPars3P2[3] = {3,0,0};
      double PotPars1D2[3] = {4,0,0};
      cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 3, true);
      cPotPars1S0->SetParameters(PotPars1S0);
      cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 3, true);
      cPotPars3P0->SetParameters(PotPars3P0);
      cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 3, true);
      cPotPars3P1->SetParameters(PotPars3P1);
      cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 3, true);
      cPotPars3P2->SetParameters(PotPars3P2);
      cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 3, true);
      cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if (POT == "ReidV8")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8] = {NN_ReidV8, v18_Coupled3P2, 1, 1, 1, 0, 0, 0};
        double PotPars3P0[8] = {NN_ReidV8, v18_Coupled3P2, 1, 1, 1, 1, 1, 0};
        double PotPars3P1[8] = {NN_ReidV8, v18_Coupled3P2, 1, 1, 1, 1, 1, 1};
        double PotPars3P2[8] = {NN_ReidV8, v18_Coupled3P2, 1, 1, 1, 1, 1, 2};
        double PotPars1D2[8] = {NN_ReidV8, v18_Coupled3P2, 1, 1, 1, 0, 2, 2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P2->SetParameters(PotPars3P2);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1D2->SetParameters(PotPars1D2);
    }
    else if (POT == "ReidSC")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8] = {pp_ReidSC, 0, 1, 1, 1, 0, 0, 0};
        double PotPars3P0[8] = {pp_ReidSC, 0, 1, 1, 1, 1, 1, 0};
        double PotPars3P1[8] = {pp_ReidSC, 0, 1, 1, 1, 1, 1, 1};
        double PotPars3P2[8] = {pp_ReidSC, 0, 1, 1, 1, 1, 1, 2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P0->SetParameters(PotPars3P0);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P1->SetParameters(PotPars3P1);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3P2->SetParameters(PotPars3P2);
    }
    else if (POT == "Norfolk")
    {
        // printf("SetUpNorfolk...\n");
        SetUpNorfolk("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/NorfolkPotential/OriginalCode/Fine/");
        // printf("SetUpNorfolk!\n");
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 3, true);
        cPotPars3P0 = new CATSparameters(CATSparameters::tPotential, 3, true);
        cPotPars3P1 = new CATSparameters(CATSparameters::tPotential, 3, true);
        cPotPars3P2 = new CATSparameters(CATSparameters::tPotential, 3, true);
        cPotPars1D2 = new CATSparameters(CATSparameters::tPotential, 3, true);

        cPotPars1S0->SetParameter(0, PotVar / 10); // lpot
        cPotPars3P0->SetParameter(0, PotVar / 10); // lpot
        cPotPars3P1->SetParameter(0, PotVar / 10); // lpot
        cPotPars3P2->SetParameter(0, PotVar / 10); // lpot
        cPotPars1D2->SetParameter(0, PotVar / 10); // lpot

        // lemp
        if (PotVar % 10 == 1)
        {
            cPotPars1S0->SetParameter(1, 100);
            cPotPars3P0->SetParameter(1, 100);
            cPotPars3P1->SetParameter(1, 100);
            cPotPars3P2->SetParameter(1, 100);
            cPotPars1D2->SetParameter(1, 100);
        }
        else if (PotVar % 10 == 2)
        {
            cPotPars1S0->SetParameter(1, 101);
            cPotPars3P0->SetParameter(1, 101);
            cPotPars3P1->SetParameter(1, 101);
            cPotPars3P2->SetParameter(1, 101);
            cPotPars1D2->SetParameter(1, 101);
        }
        else
        {
            cPotPars1S0->SetParameter(1, -1);
            cPotPars3P0->SetParameter(1, -1);
            cPotPars3P1->SetParameter(1, -1);
            cPotPars3P2->SetParameter(1, -1);
            cPotPars1D2->SetParameter(1, -1);
        }

        cPotPars1S0->SetParameter(2, 100);
        cPotPars3P0->SetParameter(2, 310);
        cPotPars3P1->SetParameter(2, 311);
        cPotPars3P2->SetParameter(2, 312);
        cPotPars1D2->SetParameter(2, 122);
    }
    else if(POT == "Johann"){
        ExternalWF = Init_pp_Haidenbauer(CatsFilesFolder[0],Kitty,PotVar);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pp potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pp;
    }
    Kitty.SetMomentumDependentSource(false);
    // Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    if(POT == "AV18_pn")
      Kitty.SetQ1Q2(0);
    else if(POT == "AV18_nn")
      Kitty.SetQ1Q2(0);
    else
      Kitty.SetQ1Q2(1);

    if(POT == "AV18_pn")
      Kitty.SetQuantumStatistics(false);
    else if(POT == "AV18_nn")
      Kitty.SetQuantumStatistics(true);
    else
      Kitty.SetPdgId(2212, 2212);


    
    if(POT == "AV18_pn")
      Kitty.SetRedMass((Mass_p*Mass_n)/(Mass_p+Mass_n));
    else if(POT == "AV18_nn")
        Kitty.SetRedMass(0.5 * Mass_n);
    else
      Kitty.SetRedMass(0.5 * Mass_p);

    if(!Kitty.GetNumChannels()){
        Kitty.SetNumChannels(4);
        Kitty.SetNumPW(0, 3);
        Kitty.SetNumPW(1, 2);
        Kitty.SetNumPW(2, 2);
        Kitty.SetNumPW(3, 2);
        Kitty.SetSpin(0, 0);
        Kitty.SetSpin(1, 1);
        Kitty.SetSpin(2, 1);
        Kitty.SetSpin(3, 1);
        Kitty.SetChannelWeight(0, 3. / 12.);
        Kitty.SetChannelWeight(1, 1. / 12.);
        Kitty.SetChannelWeight(2, 3. / 12.);
        Kitty.SetChannelWeight(3, 5. / 12.);
    }


    if (POT == "Norfolk")
    {
        if (cPotPars1S0)
            Kitty.SetShortRangePotential(0, 0, pp_Norfolk, *cPotPars1S0);
        if (cPotPars1D2)
            Kitty.SetShortRangePotential(0, 2, pp_Norfolk, *cPotPars1D2);
        if (cPotPars3P0)
            Kitty.SetShortRangePotential(1, 1, pp_Norfolk, *cPotPars3P0);
        if (cPotPars3P1)
            Kitty.SetShortRangePotential(2, 1, pp_Norfolk, *cPotPars3P1);
        if (cPotPars3P2)
            Kitty.SetShortRangePotential(3, 1, pp_Norfolk, *cPotPars3P2);
    }
    else if(POT == "pp_AV18_Toy"){
      if (cPotPars1S0)
          Kitty.SetShortRangePotential(0, 0, pp_AV18_Toy, *cPotPars1S0);
      if (cPotPars1D2)
          Kitty.SetShortRangePotential(0, 2, pp_AV18_Toy, *cPotPars1D2);
      if (cPotPars3P0)
          Kitty.SetShortRangePotential(1, 1, pp_AV18_Toy, *cPotPars3P0);
      if (cPotPars3P1)
          Kitty.SetShortRangePotential(2, 1, pp_AV18_Toy, *cPotPars3P1);
      if (cPotPars3P2)
          Kitty.SetShortRangePotential(3, 1, pp_AV18_Toy, *cPotPars3P2);
    }
    else if(POT == "Johann"){
        unsigned uCh = 0;
        Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
        if(PotVar==11){
            uCh=1; Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
            uCh=2; Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
            uCh=3; Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
        }
    }
    else
    {
        if (cPotPars1S0)
            Kitty.SetShortRangePotential(0, 0, fDlmPot, *cPotPars1S0);
        if (cPotPars1P1)
            Kitty.SetShortRangePotential(3, 1, fDlmPot, *cPotPars1P1);
        if (cPotPars1D2)
            Kitty.SetShortRangePotential(0, 2, fDlmPot, *cPotPars1D2);
        if (cPotPars3S1){
          Kitty.SetShortRangePotential(1, 1, fDlmPot, *cPotPars3S1);
          if (cPotPars3P1){
            Kitty.SetShortRangePotential(2, 1, fDlmPot, *cPotPars3S1);
          }
          if (cPotPars3P2){
            Kitty.SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);
          }
        }
        if (cPotPars3P0)
            Kitty.SetShortRangePotential(1, 1, fDlmPot, *cPotPars3P0);
        if (cPotPars3P1)
            Kitty.SetShortRangePotential(2, 1, fDlmPot, *cPotPars3P1);
        if (cPotPars3P2)
            Kitty.SetShortRangePotential(3, 1, fDlmPot, *cPotPars3P2);

    }

CLEAN_SetUpCats_pp:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotPars1S0)
    {
        delete cPotPars1S0;
        cPotPars1S0 = NULL;
    }
    if (cPotPars1D2)
    {
        delete cPotPars1D2;
        cPotPars1D2 = NULL;
    }
    if (cPotPars3P0)
    {
        delete cPotPars3P0;
        cPotPars3P0 = NULL;
    }
    if (cPotPars3P1)
    {
        delete cPotPars3P1;
        cPotPars3P1 = NULL;
    }
    if (cPotPars3P2)
    {
        delete cPotPars3P2;
        cPotPars3P2 = NULL;
    }
    if (cPotPars1P1)
    {
        delete cPotPars1P1;
        cPotPars1P1 = NULL;
    }
    if (cPotPars3S1)
    {
        delete cPotPars3S1;
        cPotPars3S1 = NULL;
    }
}

// POT:
// DG_pip_d0 (double gauss potential, p_pi+ with zero effective range)
// DG_pip_d (double gauss potential, p_pi+ with non-zero effective range)
// DG_pim_d0 (double gauss potential, p_pi- with zero effective range)
// DG_pim_d (double gauss potential, p_pi- with non-zero effective range)
//if the SourceVar last digit is zero, we use EPOS, last digit 1 we use CECA (for the ang disto)
//if the number is negative -> a Lambda instead of proton (using though the same ang. distos)
void DLM_CommonAnaFunctions::SetUpCats_ppic(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    CATSparameters *cPars = NULL;
    CATSparameters *cPotPars = NULL;

    double RED_MASS;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[6].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[6].InitStability(21, 1, 2);
        if(SourceVar>=0){
            CleverMcLevyResoTM[6].InitScale(38, 0.15, 2.0);
            CleverMcLevyResoTM[6].InitRad(257 * 2, 0, 64);
            CleverMcLevyResoTM[6].InitType(2);
            CleverMcLevyResoTM[6].SetUpReso(0, 0.6422);
            CleverMcLevyResoTM[6].SetUpReso(1, 0.682);
        }
        else{
            CleverMcLevyResoTM[6].InitScale(38, 0.15, 2.0);
            CleverMcLevyResoTM[6].InitRad(257 * 2, 0, 64);
            CleverMcLevyResoTM[6].InitType(2);
            CleverMcLevyResoTM[6].SetUpReso(0, 0.6438);
            CleverMcLevyResoTM[6].SetUpReso(1, 0.682);            
        }


        const double k_CutOff = int(int(fabs(SourceVar)) / 10) * 10.;
        const int SVAR = fabs(SourceVar % 10);
        int PPid, PRid, RPid, RRid;
        // EPOS
        if (SVAR == 0)
        {
            PPid = 0;
            PRid = 1;
            RPid = 10;
            RPid = 11;
        }
        // CECA
        else if (SVAR == 1)
        {
            PPid = 100;
            PRid = 101;
            RPid = 110;
            RPid = 111;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m Unknown source variation for p-pi\n");
            goto CLEAN_SetUpCats_ppic;
        }

        Float_t Type;
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

        TFile *F_EposDisto_p_pi = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/p_pi_ALL.root");
        TNtuple *T_EposDisto_p_pi = (TNtuple *)F_EposDisto_p_pi->Get("nt_p_pi");
        T_EposDisto_p_pi->SetBranchAddress("Type", &Type);
        T_EposDisto_p_pi->SetBranchAddress("k_D", &k_D);
        T_EposDisto_p_pi->SetBranchAddress("P1", &fP1);
        T_EposDisto_p_pi->SetBranchAddress("P2", &fP2);
        T_EposDisto_p_pi->SetBranchAddress("M1", &fM1);
        T_EposDisto_p_pi->SetBranchAddress("M2", &fM2);
        T_EposDisto_p_pi->SetBranchAddress("Tau1", &Tau1);
        T_EposDisto_p_pi->SetBranchAddress("Tau2", &Tau2);
        T_EposDisto_p_pi->SetBranchAddress("AngleRcP1", &AngleRcP1);
        T_EposDisto_p_pi->SetBranchAddress("AngleRcP2", &AngleRcP2);
        T_EposDisto_p_pi->SetBranchAddress("AngleP1P2", &AngleP1P2);


        for (unsigned uEntry = 0; uEntry < T_EposDisto_p_pi->GetEntries(); uEntry++)
        {
            T_EposDisto_p_pi->GetEntry(uEntry);
            if (Type == PPid)
            {
                continue;
            }
            else if (Type == PRid)
            {
                Tau1 = 0;
                Tau2 = 1.50;
                fM2 = 1124;
                if (k_D > k_CutOff)
                    continue;
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[6].AddBGT_PR(RanVal2, cos(AngleRcP2));
            }
            else if (Type == PRid)
            {
                Tau1 = SourceVar>=0?1.65:4.69;
                Tau2 = 0;
                fM1 = SourceVar>=0?1362:1462;//this was wrong all along for the ppic, had the pionReso mass instead of the protonReso
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                CleverMcLevyResoTM[6].AddBGT_RP(RanVal1, cos(AngleRcP1));
            }
            else if (Type == RRid)
            {
                Tau1 = SourceVar>=0?1.65:4.69;;
                Tau2 = 1.50;
                if (fabs(SourceVar % 100) == 2)
                {
                    fM1 = SourceVar>=0?1362:1462;
                    fM2 = 1124;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[6].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            else
            {
                continue;
            }
        }
        delete F_EposDisto_p_pi;

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[6].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[6].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[6], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_ppic;
    }

    if(POT == ""){
        cPotPars = NULL;
    }
    else if (POT == "DG_pip_d0")
    {
        // f0 = 0.118 fm
        // d0 = -0.02 fm
        cPotPars = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars->SetParameter(0, 2.041351e+02);
        cPotPars->SetParameter(1, 5.924833e-01);
        cPotPars->SetParameter(2, -1.105715e+03);
        cPotPars->SetParameter(3, 4.042254e-01);
    }
    else if (POT == "DG_pip_d")
    {
        cPotPars = new CATSparameters(CATSparameters::tPotential, 4, true);
        // f0 = 0.118 fm
        // d0 = 1.90 fm
        cPotPars->SetParameter(0, 5.680787e+01);
        cPotPars->SetParameter(1, 7.816589e-01);
        cPotPars->SetParameter(2, -3.371110e+02);
        cPotPars->SetParameter(3, 5.692803e-01);
    }
    else if (POT == "DG_pim_d0")
    {
        // f0 = 0.123 fm
        // d0 = 0.02 fm
        cPotPars = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars->SetParameter(0, 2.061667e+02);
        cPotPars->SetParameter(1, 5.982571e-01);
        cPotPars->SetParameter(2, -1.100512e+03);
        cPotPars->SetParameter(3, 4.099134e-01);
    }
    else if (POT == "DG_pim_d")
    {
        // f0 = 0.123 fm
        // d0 = 11.54 fm
        cPotPars = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars->SetParameter(0, -3.227447e+01);
        cPotPars->SetParameter(1, 1.077665e+00);
        cPotPars->SetParameter(2, -2.228376e+02);
        cPotPars->SetParameter(3, 9.697892e-02);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing ppi potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_ppic;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(POT.Contains("pip") ? 1 : POT.Contains("pim") ? -1
                                                                : 0);
    Kitty.SetQuantumStatistics(false);
    RED_MASS = Mass_p * Mass_pic / (Mass_p + Mass_pic);
    if(SourceVar<0) RED_MASS = Mass_L * Mass_pic / (Mass_L + Mass_pic);
    Kitty.SetRedMass(RED_MASS);

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);

    if (cPotPars)
        Kitty.SetShortRangePotential(0, 0, DoubleGaussSum, *cPotPars);

CLEAN_SetUpCats_ppic:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (cPotPars)
    {
        delete cPotPars;
        cPotPars = NULL;
    }
}

void DLM_CommonAnaFunctions::SetUpCats_pipi(CATS &Kitty, const TString &SOURCE, const int &SourceVar)
{

    CATSparameters *cPars = NULL;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.6) * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.6 * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0, sqrt(1.6) * 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.6 * 1.2);
        Kitty.SetAnaSource(1, 1.6);
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McLevyNolan_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "McGauss_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[4].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[4].InitStability(21, 1, 2);
        CleverMcLevyResoTM[4].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[4].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[4].InitType(2);
        CleverMcLevyResoTM[4].SetUpReso(0, 0.682);
        CleverMcLevyResoTM[4].SetUpReso(1, 0.682);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_ResoTM back-to-back)\n");
            goto CLEAN_SetUpCats_pipi;
        }
        // EPOS, 2(4) is with fixed mass, 3 is with EPOS mass
        // 2 is with omega included, 4 is without omega, 8 is neglecting the angular distr. from EPOS and taking them random
        else
        {
            unsigned NumPart = 0;
            unsigned NumOmega = 0;
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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

            TFile *F_EposDisto_p_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForMax_pi_piReso_withOmega.root");
            TNtuple *T_EposDisto_p_pReso = (TNtuple *)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_p_pReso; uEntry++)
            {
                NumPart++;
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                // treat the omega separately
                if (fM2 > 782 && fM2 < 783)
                {
                    NumOmega++;
                    fM2 = 782.6;
                    Tau2 = 23.24;
                    if (SourceVar % 100 == 4)
                        continue;
                }
                // the avg. values below should be computed for all resonances besides omega
                else
                {
                    Tau2 = 1.5;
                    if (SourceVar % 100 == 2 || SourceVar % 100 == 4)
                    {
                        fM2 = 1124;
                    }
                }
                if (k_D > k_CutOff)
                    continue;
                // if(fM2>782&&fM2<783){
                // printf("omega\n");
                // }
                RanVal1 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                if (SourceVar % 100 == 8)
                {
                    TVector3 vRCORE(0, 0, 1.);
                    TVector3 vSR2;
                    vSR2.SetMagThetaPhi(1, asin(RanGen.Uniform(-1, 1)), RanGen.Uniform(0, 2. * Pi));
                    AngleRcP2 = vRCORE.Angle(vSR2);
                }
                CleverMcLevyResoTM[4].AddBGT_PR(RanVal1, -cos(AngleRcP2));
                CleverMcLevyResoTM[4].AddBGT_RP(RanVal1, cos(AngleRcP2));
            }
            printf("%u/%u = %f\n", NumOmega, NumPart, double(NumOmega) / double(NumPart));
            delete F_EposDisto_p_pReso;

            TFile *F_EposDisto_pReso_pReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForMax_piReso_piReso_withOmega.root");
            TNtuple *T_EposDisto_pReso_pReso = (TNtuple *)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_pReso; uEntry++)
            {
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                // treat the omega separately
                if (fM1 > 782 && fM1 < 783)
                {
                    fM1 = 782.6;
                    Tau1 = 23.24;
                    if (SourceVar % 100 == 4)
                        continue;
                }
                // the avg. values below should be computed for all resonances besides omega
                else
                {
                    Tau1 = 1.5;
                    if (SourceVar % 100 == 2 || SourceVar % 100 == 4)
                    {
                        fM1 = 1124;
                    }
                }
                // treat the omega separately
                if (fM2 > 782 && fM2 < 783)
                {
                    fM2 = 782.6;
                    Tau2 = 23.24;
                    if (SourceVar % 100 == 4)
                        continue;
                }
                // the avg. values below should be computed for all resonances besides omega
                else
                {
                    Tau2 = 1.5;
                    if (SourceVar % 100 == 2 || SourceVar % 100 == 4)
                    {
                        fM2 = 1124;
                    }
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                if (SourceVar % 100 == 8)
                {
                    TVector3 vRCORE(0, 0, 1.);
                    TVector3 vSR1;
                    vSR1.SetMagThetaPhi(1, asin(RanGen.Uniform(-1, 1)), RanGen.Uniform(0, 2. * Pi));
                    TVector3 vSR2;
                    vSR2.SetMagThetaPhi(1, asin(RanGen.Uniform(-1, 1)), RanGen.Uniform(0, 2. * Pi));
                    AngleRcP1 = vRCORE.Angle(vSR1);
                    AngleRcP2 = vRCORE.Angle(vSR2);
                    AngleP1P2 = vSR1.Angle(vSR2);
                }
                CleverMcLevyResoTM[4].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[4].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[4].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[4], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "EPOStheta")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOStheta)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pipi;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(211, 211);
    Kitty.SetRedMass(0.5 * Mass_pic);

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);

CLEAN_SetUpCats_pipi:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
}
void DLM_CommonAnaFunctions::SetUpCats_pipi_broken(CATS &Kitty, const TString &SOURCE, const int &SourceVar)
{

    CATSparameters *cPars = NULL;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.6) * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.6 * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetAnaSource(0, sqrt(1.6) * 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[4].InitStability(20, 1, 2);
        CleverLevy[4].InitScale(35, 0.25, 2.0);
        CleverLevy[4].InitRad(256, 0, 64);
        CleverLevy[4].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[4], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.6 * 1.2);
        Kitty.SetAnaSource(1, 1.6);
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McLevyNolan_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "McGauss_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[4].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[4].InitStability(21, 1, 2);
        CleverMcLevyResoTM[4].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[4].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[4].InitType(2);
        CleverMcLevyResoTM[4].SetUpReso(0, 0.6422);
        CleverMcLevyResoTM[4].SetUpReso(1, 0.6422);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (McGauss_ResoTM back-to-back)\n");
            goto CLEAN_SetUpCats_pipi;
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else
        {
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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
            TFile *F_EposDisto_p_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/WithoutOmega/ForMax_pi_piReso.root");
            TNtuple *T_EposDisto_p_pReso = (TNtuple *)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
            T_EposDisto_p_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_p_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_p_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_p_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_p_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_p_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_p_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_p_pReso; uEntry++)
            {
                T_EposDisto_p_pReso->GetEntry(uEntry);
                Tau1 = 0;
                // treat the omega separately
                if (fM2 > 782 && fM2 < 783)
                {
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                // the avg. values below should be computed for all resonances besides omega
                else
                {
                    Tau2 = 1.5;
                    if (SourceVar % 100 == 2)
                    {
                        fM2 = 1124;
                    }
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[4].AddBGT_PR(RanVal1, -cos(AngleRcP2));
                CleverMcLevyResoTM[4].AddBGT_RP(RanVal1, cos(AngleRcP2));
            }
            delete F_EposDisto_p_pReso;

            TFile *F_EposDisto_pReso_pReso = new TFile("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/Max/WithoutOmega/ForMax_piReso_piReso.root");
            TNtuple *T_EposDisto_pReso_pReso = (TNtuple *)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
            T_EposDisto_pReso_pReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_pReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_pReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_pReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_pReso; uEntry++)
            {
                T_EposDisto_pReso_pReso->GetEntry(uEntry);
                // treat the omega separately
                if (fM2 > 782 && fM2 < 783)
                {
                    fM1 = 782.6;
                    Tau1 = 23.24;
                    fM2 = 782.6;
                    Tau2 = 23.24;
                }
                // the avg. values below should be computed for all resonances besides omega
                else
                {
                    Tau1 = 1.5;
                    Tau2 = 1.5;
                    if (SourceVar % 100 == 2)
                    {
                        fM1 = 1124;
                        fM2 = 1124;
                    }
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[4].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_pReso;
        }

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[4].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[4].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[4], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "EPOStheta")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOStheta)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pipi;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pipi;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(211, 211);
    Kitty.SetRedMass(0.5 * Mass_pic);

    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, 0);
    Kitty.SetSpin(0, 0);
    Kitty.SetChannelWeight(0, 1.);

CLEAN_SetUpCats_pipi:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
}

// POT:
//   "LO"
//   "LO_Coupled_S"
//   "NLO"
//   "NLO_Coupled_S"
//   "Usmani"
//   "UsmaniPLB" as used Phys.Lett.B 850 (2024) 138550. 
//              The PotVar corresponds to the point number of Tab. 1
//              The NLO13(600) is 13, NLO19(600) is 19 and N2LO is 20
//   "Chiral_Coupled_SPD"
//   "Chiral2023" N.B. it has to start with that, following some stuff related to the actual variation
//      "_INFOon1S0_and_INFOon3S1"
//      see Init_pL_Haidenbauer2023 for information on the naming convension
// SourceVar:
//   0 = back-to-back
//   1 = flat theta
//   2 = EPOS theta
// PotVar (for Chiral_Coupled_SPD):
//   4 digits ABXXX
//   A is NLO13 (0), NLO19 (1), LO13 (-1)
//   B is: 0 (s waves), 1 (sd waves), 2 (spd waves)
//   XXX are the cutoff
void DLM_CommonAnaFunctions::SetUpCats_pL(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{
    CATSparameters *cPars = NULL;
    CATSparameters *pPars = NULL;

    CATSparameters *cPotPars1S0 = NULL;
    CATSparameters *cPotPars3S1 = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    unsigned NumChannels = 0;

    std::string singlet;
    std::string triplet;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE == "NULL" || SOURCE == ""){

    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.2) * 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.2 * 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, sqrt(1.2) * 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, 0.5 * 1.2 * 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussExpTotSimple_2body")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        CleverMcLevyReso[1].InitStability(21, 1, 2);
        CleverMcLevyReso[1].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[1].InitRad(257, 0, 64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0, 1);
        CleverMcLevyReso[1].InitReso(1, 1);
        // bool MassSmear = SourceVar%10;
        // double MomSmear = SourceVar/100;
        if (SourceVar == 0)
        {
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else
        {
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }

        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO);
        }
        CleverMcLevyReso[1].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[1].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[1].InitScale(35, 0.25, 2.0);
        CleverMcLevyReso[1].InitRad(257 * 2, 0, 64);
        CleverMcLevyReso[1].InitType(2);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else if (SourceVar == 1)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }
        if (SourceVar == 2)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double> *HISTO_p = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 600);
            DLM_Histo<double> *HISTO_L = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 270, 470);

            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO_p);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO_L);
        }
        // improved reso
        else if (SourceVar == 3)
        {
            CleverMcLevyReso[1].InitReso(0, 2);
            CleverMcLevyReso[1].InitReso(1, 2);
            CleverMcLevyReso[1].SetUpReso(0, 0, 0.4368, 1231.98, 1.67, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(0, 1, 0.2054, 1636.98, 1.62, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1, 0, 0.4864, 1384.54, 5.34, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1, 1, 0.1573, 1705.26, 2.70, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double> *HISTO_p_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 200, 400);
            DLM_Histo<double> *HISTO_p_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 880, 1000);
            DLM_Histo<double> *HISTO_L_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 160, 360);
            DLM_Histo<double> *HISTO_L_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 620, 820);

            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO_p_0);
            CleverMcLevyReso[1].SetUpResoEmission(0, 1, HISTO_p_1);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO_L_0);
            CleverMcLevyReso[1].SetUpResoEmission(1, 1, HISTO_L_1);
        }

        CleverMcLevyReso[1].InitNumMcIter(SourceVar == 3 ? 100000 : 1000000);

        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[1].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[1].InitStability(21, 1, 2);
        CleverMcLevyResoTM[1].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[1].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[1].InitType(2);
        CleverMcLevyResoTM[1].SetUpReso(0, 0.6422);
        CleverMcLevyResoTM[1].SetUpReso(1, 0.6438);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            CleverMcLevyResoTM[1].AddBGT_PR(360. / 1462. * 4.69, 1.);
            CleverMcLevyResoTM[1].AddBGT_RP(490. / 1362. * 1.65, -1.);
            CleverMcLevyResoTM[1].AddBGT_RR(490. / 1362. * 1.65, -1., 360. / 1462. * 4.69, 1., -1.);
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else
        {
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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

            TFile *F_EposDisto_p_LamReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_p_LamReso.root");
            TNtuple *T_EposDisto_p_LamReso = (TNtuple *)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
            T_EposDisto_p_LamReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_p_LamReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_p_LamReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_p_LamReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_p_LamReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_p_LamReso; uEntry++)
            {
                T_EposDisto_p_LamReso->GetEntry(uEntry);
                Tau1 = 0;
                Tau2 = 4.69;
                if (SourceVar % 100 == 2)
                {
                    fM2 = 1462;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[1].AddBGT_PR(RanVal2, cos(AngleRcP2));
            }
            delete F_EposDisto_p_LamReso;

            TFile *F_EposDisto_pReso_Lam = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_Lam.root");
            TNtuple *T_EposDisto_pReso_Lam = (TNtuple *)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
            T_EposDisto_pReso_Lam->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_Lam->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_Lam->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Lam; uEntry++)
            {
                T_EposDisto_pReso_Lam->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                CleverMcLevyResoTM[1].AddBGT_RP(RanVal1, cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Lam;

            TFile *F_EposDisto_pReso_LamReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_LamReso.root");
            TNtuple *T_EposDisto_pReso_LamReso = (TNtuple *)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
            T_EposDisto_pReso_LamReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_LamReso; uEntry++)
            {
                T_EposDisto_pReso_LamReso->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 4.69;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                    fM2 = 1462;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[1].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_LamReso;
        }

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[1].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[1].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pL;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pL;
    }

    if (POT == "LO")
    {
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaLO_600/",
                                         Kitty, 0, 600);
        NumChannels = 2;
    }
    else if (POT == "LO_Coupled_S")
    {
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaLO_Coupling/",
                                         Kitty, 1, 600);
        NumChannels = 4;
    }
    else if (POT == "NLO")
    {
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaNLO/",
                                         Kitty, 10, 600);
        NumChannels = 2;
    }
    // s and p waves
    else if (POT == "NLO_sp")
    {
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaNLO/",
                                         Kitty, 12, 600);
        NumChannels = 4;
    }
    else if (POT == "NLO_Coupled_S")
    {
        ExternalWF = Init_pL_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaNLO_Coupling/",
                                         Kitty, 11, 600);
        NumChannels = 4;
    }
    else if (POT == "Chiral_Coupled_SPD")
    {
        int CUTOFF = abs(PotVar % 1000);
        int TYPE = PotVar / 10000;
        if (PotVar == 0)
        {
            CUTOFF = 600;
            TYPE = 0;
        }
        // printf("PotVar=%i\n",PotVar);
        // printf("CUTOFF=%i\n",CUTOFF);
         //printf("TYPE=%i\n",TYPE);
        // original, NLO13 at 600 MeV
        if(TYPE==0 || TYPE==1 || TYPE==-1){
          ExternalWF = Init_pL_Haidenbauer2019(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambda_Coupled_SD/",
                                               Kitty, TYPE, CUTOFF);
        }
        //
        else if(TYPE==2){
          ExternalWF = Init_pL_Haidenbauer2019(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pLambdaN2LO/",
                                               Kitty, TYPE, CUTOFF);
        }
        else{
          printf("\033[1;31mERROR:\033[0m Non-existing pL potential '%s'\n", POT.Data());
          goto CLEAN_SetUpCats_pL;
        }

        NumChannels = 16;
    }
    else if(POT.BeginsWith("Chiral2023_")){
        std::string inputString = POT.Data();
        // Find the positions of "_" and "and"
        size_t underscorePos = inputString.find("_");
        size_t andPos = inputString.find("_and_");
        singlet = inputString.substr(underscorePos + 1, andPos - underscorePos - 1);
        triplet = inputString.substr(andPos + 5);
        ExternalWF = Init_pL_Haidenbauer2023(CatsFilesFolder[0],Kitty,singlet.c_str(),triplet.c_str());
        NumChannels = 16;        
    }
    else if (POT == "Usmani")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars1S0[8] = {pL_UsmaniOli, 0, 0, 0, 0, 0, 0, 0};
        double PotPars3S1[8] = {pL_UsmaniOli, 0, 0, 0, 0, 1, 0, 1};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    else if (POT == "UsmaniFit")
    {
        double PotPars1S0[4] = {0, 2137, 0.5, 0.2};
        double PotPars3S1[4] = {1, 2137, 0.5, 0.2};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    else if (POT == "UsmaniPLB")
    {
        double PotPars1S0[4] = {0, 2137, 0.5, 0.2};
        double PotPars3S1[4] = {1, 2137, 0.5, 0.2};
        if(PotVar==1){
            PotPars1S0[1] = 2029.21;
            PotPars1S0[2] = 0.516704;
            PotPars1S0[3] = 0.201386;
            PotPars3S1[1] = 1900.49;
            PotPars3S1[2] = 0.447223;
            PotPars3S1[3] = 0.238315;
        }
        else if(PotVar==2){
            PotPars1S0[1] = 2084.27;
            PotPars1S0[2] = 0.501872;
            PotPars1S0[3] = 0.204962;
            PotPars3S1[1] = 2069.51;
            PotPars3S1[2] = 0.442597;
            PotPars3S1[3] = 0.228846;
        }
        else if(PotVar==3){
            PotPars1S0[1] = 2079.81;
            PotPars1S0[2] = 0.519965;
            PotPars1S0[3] = 0.197214;
            PotPars3S1[1] = 1983.95;
            PotPars3S1[2] = 0.456457;
            PotPars3S1[3] = 0.227331;
        }
        else if(PotVar==4){
            PotPars1S0[1] = 2111.71;
            PotPars1S0[2] = 0.512425;
            PotPars1S0[3] = 0.197126;
            PotPars3S1[1] = 2024.79;
            PotPars3S1[2] = 0.428883;
            PotPars3S1[3] = 0.238985;
        }        
        else if(PotVar==5){
            PotPars1S0[1] = 2099.52;
            PotPars1S0[2] = 0.518209;
            PotPars1S0[3] = 0.195196;
            PotPars3S1[1] = 1932.33;
            PotPars3S1[2] = 0.456657;
            PotPars3S1[3] = 0.231787;
        }       
        else if(PotVar==6){
            PotPars1S0[1] = 2169.98;
            PotPars1S0[2] = 0.506783;
            PotPars1S0[3] = 0.196690;
            PotPars3S1[1] = 2129.60;
            PotPars3S1[2] = 0.438970;
            PotPars3S1[3] = 0.227010;
        }    
        else if(PotVar==7){
            PotPars1S0[1] = 1976.88;
            PotPars1S0[2] = 0.520216;
            PotPars1S0[3] = 0.199594;
            PotPars3S1[1] = 2057.49;
            PotPars3S1[2] = 0.443663;
            PotPars3S1[3] = 0.230643;
        }    
        else if(PotVar==8){
            PotPars1S0[1] = 2098.16;
            PotPars1S0[2] = 0.512709;
            PotPars1S0[3] = 0.195045;
            PotPars3S1[1] = 1957.61;
            PotPars3S1[2] = 0.434795;
            PotPars3S1[3] = 0.241739;
        }    
        else if(PotVar==13){
            PotPars1S0[1] = 2155.97;
            PotPars1S0[2] = 0.514819;
            PotPars1S0[3] = 0.192401;
            PotPars3S1[1] = 2176.83;
            PotPars3S1[2] = 0.432072;
            PotPars3S1[3] = 0.227458;
        }   
        else if(PotVar==19){
            PotPars1S0[1] = 1950.83;
            PotPars1S0[2] = 0.532720;
            PotPars1S0[3] = 0.195344;
            PotPars3S1[1] = 2009.28;
            PotPars3S1[2] = 0.449633;
            PotPars3S1[3] = 0.230304;
        }  
        else if(PotVar==20){
            PotPars1S0[1] = 2098.16;
            PotPars1S0[2] = 0.512709;
            PotPars1S0[3] = 0.195045;
            PotPars3S1[1] = 1957.61;
            PotPars3S1[2] = 0.434795;
            PotPars3S1[3] = 0.241739;
        }  

        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    else if (POT == "UsmaniFitAll")
    {
        double PotPars1S0[8] = {0, 2137, 0.5, 0.2, 6.2, 0.25, 0.7, 2.0};
        double PotPars3S1[8] = {1, 2137, 0.5, 0.2, 6.2, 0.25, 0.7, 2.0};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pL potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pL;
    }
    // printf("NumChannels=%u\n",NumChannels);
    // printf("ExternalWF=%p\n",ExternalWF);

    Kitty.SetMomentumDependentSource(false);
    // Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    Kitty.SetRedMass((Mass_p * Mass_L) / (Mass_p + Mass_L));

    Kitty.SetNumChannels(NumChannels);
    for (unsigned uCh = 0; uCh < NumChannels; uCh++)
    {
        if (!ExternalWF)
        {
            Kitty.SetNumPW(uCh, 1);
            Kitty.SetSpin(uCh, uCh % 2 == 0 ? 0 : 1);
            Kitty.SetChannelWeight(uCh, uCh % 2 == 0 ? 0.25 : 0.75);
        }

        if ( (POT == "UsmaniFit"||POT == "UsmaniPLB")  && cPotPars1S0 && uCh == 0)
            Kitty.SetShortRangePotential(uCh, 0, UsmaniFit, *cPotPars1S0);
        else if ( (POT == "UsmaniFit"||POT == "UsmaniPLB") && cPotPars3S1 && uCh == 1)
            Kitty.SetShortRangePotential(uCh, 0, UsmaniFit, *cPotPars3S1);
        else if (POT == "UsmaniFitAll" && cPotPars1S0 && uCh == 0)
            Kitty.SetShortRangePotential(uCh, 0, UsmaniFitAll, *cPotPars1S0);
        else if (POT == "UsmaniFitAll" && cPotPars3S1 && uCh == 1)
            Kitty.SetShortRangePotential(uCh, 0, UsmaniFitAll, *cPotPars3S1);
        else if (cPotPars1S0 && uCh == 0)
            Kitty.SetShortRangePotential(uCh, 0, fDlmPot, *cPotPars1S0);
        else if (cPotPars3S1 && uCh == 1)
            Kitty.SetShortRangePotential(uCh, 0, fDlmPot, *cPotPars3S1);
        else if (ExternalWF)
        {
            if(false){
                    //this produces same Ck for both !!!!
                    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
                    //if (uCh <= 6)//main channels, p-waves
                    //{
                        //this produces same Ck for both !!!!
                        //Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
                    //}
                    //if (uCh >= 1 && uCh <= 3)//the d-waves for the S=1 main channel
                    //{
                        //this produces same Ck for both !!!!
                        //Kitty.SetExternalWaveFunction(uCh, 2, ExternalWF[0][uCh][2], ExternalWF[1][uCh][2]);
                    //}                   

            }
            else if(POT.BeginsWith("Chiral2023_")){

                // Use strtok to split the string
                std::string separator = "_";

                // Create a stringstream to tokenize the string
                std::stringstream ss1S0(singlet);
                std::string token1S0;
                // Tokenize the string using the separator
                std::vector<std::string> tokens1S0;
                while (std::getline(ss1S0, token1S0, separator[0])) {
                    tokens1S0.push_back(token1S0);
                }

                // Create a stringstream to tokenize the string
                std::stringstream ss3S1(triplet);
                std::string token3S1;
                // Tokenize the string using the separator
                std::vector<std::string> tokens3S1;
                while (std::getline(ss3S1, token3S1, separator[0])) {
                    tokens3S1.push_back(token3S1);
                }
                
                if(token1S0<4){
                    printf("\033[1;31mERROR:\033[0m Something wrong with the pL set up (%s)!\n",singlet.c_str());
                }
                if(token3S1<4){
                    printf("\033[1;31mERROR:\033[0m Something wrong with the pL set up (%s)!\n",singlet.c_str());
                }
                //N.B. all CC are listed as s-waves for some technical reason. So it is a bit of a mess,
                //but some are actually s-waves, others, p or d. This is taken care of below

                bool ElasticOnly1S0 = false;
                bool ElasticOnly3S1 = false;
                if(tokens1S0.size()>=5){
                    for(unsigned uLen=4; uLen<tokens1S0.size(); uLen++)
                    if(tokens1S0.at(uLen)=="elastic"){
                        ElasticOnly1S0 = true;
                    }
                }

                if(tokens3S1.size()>=5){
                    for(unsigned uLen=4; uLen<tokens3S1.size(); uLen++)
                    if(tokens3S1.at(uLen)=="elastic"){
                        ElasticOnly3S1 = true;
                    }
                }  

    //main channels:
    //0: 1S0+1P1
    //1: 3S1+3P0+3D1
    //2: 3S1+3P1+3D1
    //3: 3S1+3P2+3D1
    //4: 3S1+3P0
    //5: 3S1+3P1
    //6: 3S1+3P2
    //coupled channels:
    //7: 1S0 SN(s) -> LN(s)
    //8: 3S1 SN(s) -> LN(s)
    //9: 3S1 LN(d) -> LN(s)
    //10: 3S1 SN(d) -> LN(s)
    //11: 3P0 SN(p) -> LN(p)//
    //12: 3P2 SN(p) -> LN(p)//
    //13: 3D1 SN(d) -> LN(d)
    //14: 3D1 LN(s) -> LN(d)
    //15: 3D1 SN(s) -> LN(d)

//Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);

                // for(unsigned uMomBin=0; uMomBin<Kitty.GetNumMomBins(); uMomBin++){
                // Kitty.UseExternalWaveFunction(uMomBin,uCh,0,WaveFunctionU[uMomBin][uCh][0], NumRadBins, RadBins, PhaseShifts[uMomBin][uCh][0]);
                // printf("Look at that view (%u)!\n",uCh);


                //the s-wave channels for the singlet
                if(tokens1S0.at(3).find("s")!=std::string::npos){
                    if(uCh==0)
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}
                    if(uCh==7 && !ElasticOnly1S0)
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}                  
                }
                //the s-wave channels for the triplet
                if(tokens3S1.at(3).find("s")!=std::string::npos){
                    if(uCh>=1 && uCh<=6)//main channels
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}
                    if(uCh>=8 && uCh<=10 && !ElasticOnly3S1)//coupled channels
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}                   
                }

                //the p-wave channels for the singlet
                if(tokens1S0.at(3).find("p")!=std::string::npos){
                    if(uCh==0)
                        {Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);}
                }
                //the p-wave channels for the triplet
                if(tokens3S1.at(3).find("p")!=std::string::npos){
                    if(uCh>=1 && uCh<=3)
                        {Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);}
                    if( (uCh==11 || uCh==12) && !ElasticOnly3S1 )
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}
                }

                //the d-wave channels for the singlet
                if(tokens1S0.at(3).find("d")!=std::string::npos){

                }
                //the d-wave channels for the triplet
                if(tokens3S1.at(3).find("d")!=std::string::npos){
                    if(uCh>=1 && uCh<=3)
                        {Kitty.SetExternalWaveFunction(uCh, 2, ExternalWF[0][uCh][2], ExternalWF[1][uCh][2]);}
                    if(uCh>=13 && uCh<=15 && !ElasticOnly3S1)
                        {Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);}
                }

            }
            else{
                // for(unsigned uMomBin=0; uMomBin<Kitty.GetNumMomBins(); uMomBin++){
                // Kitty.UseExternalWaveFunction(uMomBin,uCh,0,WaveFunctionU[uMomBin][uCh][0], NumRadBins, RadBins, PhaseShifts[uMomBin][uCh][0]);
                //printf("Look at that view (%u)!\n",uCh);
                int SPD = (PotVar / 1000) % 10;
                // unless p-wave
                if (ExternalWF[0][uCh][0].GetDim()){//sets all s-waves (main and coupled channels)
                    Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
                    //printf("%u s\n",uCh);
                }
                if (POT == "Chiral_Coupled_SPD")
                {
                    if (uCh <= 6 && SPD == 2)//main channels, p-waves
                    {
                        Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
                        //printf("%u p\n",uCh);
                        // printf("1: uCh=%u\n",uCh);
                    }
                    if (uCh >= 1 && uCh <= 3 && SPD != 0)//the d-waves for the S=1 main channel
                    {
                        Kitty.SetExternalWaveFunction(uCh, 2, ExternalWF[0][uCh][2], ExternalWF[1][uCh][2]);
                        //printf("%u d\n",uCh);
                        // printf("uCh=%u\n",uCh);
                    }
                }
                else if (Kitty.GetNumPW(uCh) >= 2)
                {
                    Kitty.SetExternalWaveFunction(uCh, 1, ExternalWF[0][uCh][1], ExternalWF[1][uCh][1]);
                    //printf("???\n");
                }

                // printf(" --Look at that view (%u)!\n",uCh);
                // }
            }
        }
        else
        {
            printf("PotVar=%i\n", PotVar);
            printf("\033[1;31mERROR:\033[0m SetUpCats_pL says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pL;
        }
    }
    // Kitty.KillTheCat();
    // printf("------------------------");
CLEAN_SetUpCats_pL:;

    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotPars1S0)
    {
        delete cPotPars1S0;
        cPotPars1S0 = NULL;
    }
    if (cPotPars3S1)
    {
        delete cPotPars3S1;
        cPotPars3S1 = NULL;
    }
    CleanUpWfHisto(Kitty, ExternalWF);
}


void DLM_CommonAnaFunctions::SetUpCats_pSp(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{
    CATSparameters *cPars = NULL;
    CATSparameters *pPars = NULL;

    CATSparameters *cPotPars1S0 = NULL;
    CATSparameters *cPotPars3S1 = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    unsigned NumChannels = 0;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE == "NULL" || SOURCE == ""){

    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.2) * 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.2 * 1.2);
        cPars->SetParameter(1, 1.2);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, sqrt(1.2) * 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[1].InitStability(20, 1, 2);
        CleverLevy[1].InitScale(35, 0.25, 2.0);
        CleverLevy[1].InitRad(256, 0, 64);
        CleverLevy[1].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[1], 2);
        Kitty.SetAnaSource(0, 0.5 * 1.2 * 1.2);
        Kitty.SetAnaSource(1, 1.2);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussExpTotSimple_2body")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pSp;
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        CleverMcLevyReso[1].InitStability(21, 1, 2);
        CleverMcLevyReso[1].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[1].InitRad(257, 0, 64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0, 1);
        CleverMcLevyReso[1].InitReso(1, 1);
        // bool MassSmear = SourceVar%10;
        // double MomSmear = SourceVar/100;
        if (SourceVar == 0)
        {
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else
        {
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }

        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO);
        }
        CleverMcLevyReso[1].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[1].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[1].InitScale(35, 0.25, 2.0);
        CleverMcLevyReso[1].InitRad(257 * 2, 0, 64);
        CleverMcLevyReso[1].InitType(2);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else if (SourceVar == 1)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }
        if (SourceVar == 2)
        {
            CleverMcLevyReso[1].InitReso(0, 1);
            CleverMcLevyReso[1].InitReso(1, 1);
            CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double> *HISTO_p = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 600);
            DLM_Histo<double> *HISTO_L = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 270, 470);

            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO_p);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO_L);
        }
        // improved reso
        else if (SourceVar == 3)
        {
            CleverMcLevyReso[1].InitReso(0, 2);
            CleverMcLevyReso[1].InitReso(1, 2);
            CleverMcLevyReso[1].SetUpReso(0, 0, 0.4368, 1231.98, 1.67, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(0, 1, 0.2054, 1636.98, 1.62, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1, 0, 0.4864, 1384.54, 5.34, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
            CleverMcLevyReso[1].SetUpReso(1, 1, 0.1573, 1705.26, 2.70, Mass_L, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);

            DLM_Histo<double> *HISTO_p_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 200, 400);
            DLM_Histo<double> *HISTO_p_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 880, 1000);
            DLM_Histo<double> *HISTO_L_0 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 160, 360);
            DLM_Histo<double> *HISTO_L_1 = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/EPOS_HM_4PI/pp/DimiPhi_pp.root", "h_rkAngle_Mom2", 620, 820);

            CleverMcLevyReso[1].SetUpResoEmission(0, 0, HISTO_p_0);
            CleverMcLevyReso[1].SetUpResoEmission(0, 1, HISTO_p_1);
            CleverMcLevyReso[1].SetUpResoEmission(1, 0, HISTO_L_0);
            CleverMcLevyReso[1].SetUpResoEmission(1, 1, HISTO_L_1);
        }

        CleverMcLevyReso[1].InitNumMcIter(SourceVar == 3 ? 100000 : 1000000);

        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[1].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[1].InitStability(21, 1, 2);
        CleverMcLevyResoTM[1].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[1].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[1].InitType(2);
        CleverMcLevyResoTM[1].SetUpReso(0, 0.6422);
        CleverMcLevyResoTM[1].SetUpReso(1, 0.6438);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            CleverMcLevyResoTM[1].AddBGT_PR(360. / 1462. * 4.69, 1.);
            CleverMcLevyResoTM[1].AddBGT_RP(490. / 1362. * 1.65, -1.);
            CleverMcLevyResoTM[1].AddBGT_RR(490. / 1362. * 1.65, -1., 360. / 1462. * 4.69, 1., -1.);
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else
        {
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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

            TFile *F_EposDisto_p_LamReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_p_LamReso.root");
            TNtuple *T_EposDisto_p_LamReso = (TNtuple *)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
            T_EposDisto_p_LamReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_p_LamReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_p_LamReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_p_LamReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_p_LamReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_p_LamReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_p_LamReso; uEntry++)
            {
                T_EposDisto_p_LamReso->GetEntry(uEntry);
                Tau1 = 0;
                Tau2 = 4.69;
                if (SourceVar % 100 == 2)
                {
                    fM2 = 1462;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[1].AddBGT_PR(RanVal2, cos(AngleRcP2));
            }
            delete F_EposDisto_p_LamReso;

            TFile *F_EposDisto_pReso_Lam = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_Lam.root");
            TNtuple *T_EposDisto_pReso_Lam = (TNtuple *)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
            T_EposDisto_pReso_Lam->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_Lam->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_Lam->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_Lam->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Lam; uEntry++)
            {
                T_EposDisto_pReso_Lam->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                CleverMcLevyResoTM[1].AddBGT_RP(RanVal1, cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Lam;

            TFile *F_EposDisto_pReso_LamReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_LamReso.root");
            TNtuple *T_EposDisto_pReso_LamReso = (TNtuple *)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
            T_EposDisto_pReso_LamReso->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_LamReso; uEntry++)
            {
                T_EposDisto_pReso_LamReso->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 4.69;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                    fM2 = 1462;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
                CleverMcLevyResoTM[1].AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
            }
            delete F_EposDisto_pReso_LamReso;
        }

        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[1].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[1].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[1], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pSp;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pSp;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pSp;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pSp;
    }

    //a double gaussian, that reproduces the correct scattering parameters
    if (POT == "DG_NLO19")
    {
        double PotPars1S0[4] = {-82.20693060436437, 1.3392323926534693, 2274.508479434418, 0.4182867242441114};
        double PotPars3S1[4] = {6.295180245126864, 1.995145226219902, -58.68689895288922, 0.32017712049353114};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    //a double gaussian, that reproduces the correct scattering parameters
    else if (POT == "DG_N2LO")
    {
        double PotPars1S0[4] = {-79.07298063891983, 1.3716237629795942, 2287.907638142889, 0.43371332814881547};
        double PotPars3S1[4] = {6.6923386649252645, 1.970250601599477, -48.30709539991798, 0.3549784459798298};
        cPotPars1S0 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars1S0->SetParameters(PotPars1S0);
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 4, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        NumChannels = 2;
    }
    //the WF obtained from Johann
    else if (POT == "N2LO")
    {
        ExternalWF = Init_pSigmaPlus_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pSigmaPlus/",
                                         Kitty, 0);
        NumChannels = 2;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pSigmaPlus potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pSp;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass((Mass_p * Mass_Sch) / (Mass_p + Mass_Sch));

    Kitty.SetNumChannels(NumChannels);
    for (unsigned uCh = 0; uCh < NumChannels; uCh++)
    {
        if (!ExternalWF)
        {
            Kitty.SetNumPW(uCh, 1);
            Kitty.SetSpin(uCh, uCh % 2 == 0 ? 0 : 1);
            Kitty.SetChannelWeight(uCh, uCh % 2 == 0 ? 0.25 : 0.75);
        }

        if (POT.Contains("DG_") && cPotPars1S0 && uCh == 0)
            Kitty.SetShortRangePotential(uCh, 0, DoubleGaussSum, *cPotPars1S0);
        else if (POT.Contains("DG_") && cPotPars3S1 && uCh == 1)
            Kitty.SetShortRangePotential(uCh, 0, DoubleGaussSum, *cPotPars3S1);
        else if (ExternalWF)
        {
            Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
        }
        else
        {
            printf("PotVar=%i\n", PotVar);
            printf("\033[1;31mERROR:\033[0m SetUpCats_pSp says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pSp;
        }
    }
CLEAN_SetUpCats_pSp:;

    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    if (cPotPars1S0)
    {
        delete cPotPars1S0;
        cPotPars1S0 = NULL;
    }
    if (cPotPars3S1)
    {
        delete cPotPars3S1;
        cPotPars3S1 = NULL;
    }
    CleanUpWfHisto(Kitty, ExternalWF);
}



// for the moment only gauss
// models are NLO and ESC16
void DLM_CommonAnaFunctions::SetUpCats_pS0(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{
    CATSparameters *cPars = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    unsigned NumChannels = 0;

    Kitty.SetThetaDependentSource(false);

    if(SOURCE == "NULL" || SOURCE == ""){

    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pS0;
    }

    if (POT == "Chiral")
    {
        ExternalWF = Init_pSigma0_Haidenbauer(CatsFilesFolder[0] + "/Interaction/Haidenbauer/pSigma0/", Kitty);
        NumChannels = Kitty.GetNumChannels();
    }
    else if (POT == "ESC16")
    {
        ExternalWF = Init_pS0_ESC16(CatsFilesFolder[0] + "/Interaction/Tom/pSigma0/DimiValeNorm170519/", Kitty);
        NumChannels = Kitty.GetNumChannels();
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pS0 potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pS0;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    for (unsigned uCh = 0; uCh < NumChannels; uCh++)
    {
        if (!ExternalWF)
        {
            Kitty.SetSpin(uCh, uCh % 2 == 0 ? 0 : 1);
            Kitty.SetChannelWeight(uCh, uCh % 2 == 0 ? 0.25 : 0.75);
        }
        else if (ExternalWF)
        {
            Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m SetUpCats_pL says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pS0;
        }
    }
CLEAN_SetUpCats_pS0:;

    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    CleanUpWfHisto(Kitty, ExternalWF);
}

// POT:
//   "pXim_Lattice" (the first version)
//   "pXim_HALQCD1" (the second version)
void DLM_CommonAnaFunctions::SetUpCats_pXim(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    double POTFLAG;
    switch (PotVar)
    {
    case 11:
        POTFLAG = 11;
        break;
    case 12:
        POTFLAG = 12;
        break;
    case 13:
        POTFLAG = 13;
        break;
    case -11:
        POTFLAG = -11;
        break;
    case -12:
        POTFLAG = -12;
        break;
    case -13:
        POTFLAG = -13;
        break;
    default:
        POTFLAG = 12;
        break;
    }

    CATSparameters *cPars = NULL;
    CATSparameters *pPars = NULL;

    CATSparameters *cPotParsI0S0 = NULL;
    CATSparameters *cPotParsI0S1 = NULL;
    CATSparameters *cPotParsI1S0 = NULL;
    CATSparameters *cPotParsI1S1 = NULL;

    if(SOURCE == "NULL" || SOURCE == ""){

    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.8) * 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.8 * 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[2].InitStability(20, 1, 2);
        CleverLevy[2].InitScale(35, 0.25, 2.0);
        CleverLevy[2].InitRad(256, 0, 64);
        CleverLevy[2].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.8);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[2].InitStability(20, 1, 2);
        CleverLevy[2].InitScale(35, 0.25, 2.0);
        CleverLevy[2].InitRad(256, 0, 64);
        CleverLevy[2].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetAnaSource(0, sqrt(1.8) * 1.2);
        Kitty.SetAnaSource(1, 1.8);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[2].InitStability(20, 1, 2);
        CleverLevy[2].InitScale(35, 0.25, 2.0);
        CleverLevy[2].InitRad(256, 0, 64);
        CleverLevy[2].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[2], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.8 * 1.2);
        Kitty.SetAnaSource(1, 1.8);
    }
    else if (SOURCE == "GaussExpTotSimple_2body")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        CleverMcLevyReso[2].InitStability(21, 1, 2);
        CleverMcLevyReso[2].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[2].InitRad(257, 0, 64);
        CleverMcLevyReso[2].InitType(2);
        CleverMcLevyReso[2].InitReso(0, 1);
        CleverMcLevyReso[2].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[2].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
        }
        else
        {
            CleverMcLevyReso[2].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
        }
        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[2].SetUpResoEmission(0, 0, HISTO);
        }

        CleverMcLevyReso[2].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[2], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[2].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[2].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[2].InitRad(257, 0, 64);
        CleverMcLevyReso[2].InitType(2);
        CleverMcLevyReso[2].InitReso(0, 1);
        CleverMcLevyReso[2].InitReso(1, 1);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[2].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[2].SetUpReso(1, 0, 0.0001, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);      // rdtBackwards,rdtRandom
        }
        else
        {
            CleverMcLevyReso[2].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[2].SetUpReso(1, 0, 0.0001, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);      // rdtBackwards,rdtRandom
        }
        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[2].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[2].SetUpResoEmission(1, 0, HISTO);
        }
        CleverMcLevyReso[2].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[2], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM")
    {
        CleverMcLevyResoTM[2].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyResoTM[2].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[2].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[2].InitType(2);
        CleverMcLevyResoTM[2].SetUpReso(0, 0.6422);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            CleverMcLevyResoTM[2].AddBGT_RP(490. / 1362. * 1.65, -1.);
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else
        {
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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

            TFile *F_EposDisto_pReso_Xim = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_Xim.root");
            TNtuple *T_EposDisto_pReso_Xim = (TNtuple *)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
            T_EposDisto_pReso_Xim->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_Xim->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_Xim->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_Xim->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_Xim->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_Xim->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_Xim->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Xim; uEntry++)
            {
                T_EposDisto_pReso_Xim->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                CleverMcLevyResoTM[2].AddBGT_RP(RanVal1, cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Xim;
        }

        CleverMcLevyResoTM[2].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[2], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pXim;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pXim;
    }

    if (POT == "pXim_Lattice")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9] = {pXim_Lattice, 12, 0, -1, 1, 0, 0, 0, 0};
        double PotParsI0S1[9] = {pXim_Lattice, 12, 0, -1, 1, 1, 0, 1, 0};
        double PotParsI1S0[9] = {pXim_Lattice, 6, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_Lattice, 6, 1, 1, 1, 1, 0, 1, 0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if (POT == "pXim_HALQCD1")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9] = {pXim_HALQCD1, POTFLAG, 0, -1, 1, 0, 0, 0, 0};
        double PotParsI0S1[9] = {pXim_HALQCD1, POTFLAG, 0, -1, 1, 1, 0, 1, 0};
        double PotParsI1S0[9] = {pXim_HALQCD1, POTFLAG, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_HALQCD1, POTFLAG, 1, 1, 1, 1, 0, 1, 0};
        // printf("POTFLAG = %f\n",POTFLAG);
        // printf("PotVar = %i\n",PotVar);
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if (POT == "pXim_HALQCDPaper2020")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI0S0[9] = {pXim_HALQCDPaper2020, POTFLAG, 0, -1, 1, 0, 0, 0, 0};
        double PotParsI0S1[9] = {pXim_HALQCDPaper2020, POTFLAG, 0, -1, 1, 1, 0, 1, 0};
        double PotParsI1S0[9] = {pXim_HALQCDPaper2020, POTFLAG, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_HALQCDPaper2020, POTFLAG, 1, 1, 1, 1, 0, 1, 0};
        cPotParsI0S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S0->SetParameters(PotParsI0S0);
        cPotParsI0S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI0S1->SetParameters(PotParsI0S1);
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if (POT == "pXim1530")
    {
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pXi potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pXim;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    if (POT == "pXim1530")
    {
        Kitty.SetRedMass((Mass_p * Mass_Xim1530) / (Mass_p + Mass_Xim1530));

        Kitty.SetNumChannels(1);
        Kitty.SetNumPW(0, 0);
        Kitty.SetSpin(0, 0);
        Kitty.SetChannelWeight(0, 1.);
    }
    else
    {
        Kitty.SetRedMass((Mass_p * Mass_Xim) / (Mass_p + Mass_Xim));

        Kitty.SetNumChannels(4);
        Kitty.SetNumPW(0, 1);
        Kitty.SetNumPW(1, 1);
        Kitty.SetNumPW(2, 1);
        Kitty.SetNumPW(3, 1);
        Kitty.SetSpin(0, 0);
        Kitty.SetSpin(1, 1);
        Kitty.SetSpin(2, 0);
        Kitty.SetSpin(3, 1);
        Kitty.SetChannelWeight(0, 1. / 8.);
        Kitty.SetChannelWeight(1, 3. / 8.);
        Kitty.SetChannelWeight(2, 1. / 8.);
        Kitty.SetChannelWeight(3, 3. / 8.);

        if (cPotParsI0S0)
            Kitty.SetShortRangePotential(0, 0, fDlmPot, *cPotParsI0S0);
        if (cPotParsI0S1)
            Kitty.SetShortRangePotential(1, 0, fDlmPot, *cPotParsI0S1);
        if (cPotParsI1S0)
            Kitty.SetShortRangePotential(2, 0, fDlmPot, *cPotParsI1S0);
        if (cPotParsI1S1)
            Kitty.SetShortRangePotential(3, 0, fDlmPot, *cPotParsI1S1);
    }

CLEAN_SetUpCats_pXim:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotParsI0S0)
    {
        delete cPotParsI0S0;
        cPotParsI0S0 = NULL;
    }
    if (cPotParsI0S1)
    {
        delete cPotParsI0S1;
        cPotParsI0S1 = NULL;
    }
    if (cPotParsI1S0)
    {
        delete cPotParsI1S0;
        cPotParsI1S0 = NULL;
    }
    if (cPotParsI1S1)
    {
        delete cPotParsI1S1;
        cPotParsI1S1 = NULL;
    }
}

void DLM_CommonAnaFunctions::SetUpCats_pXi0(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    double POTFLAG;
    switch (PotVar)
    {
    case 11:
        POTFLAG = 11;
        break;
    case 12:
        POTFLAG = 12;
        break;
    case 13:
        POTFLAG = 13;
        break;
    case -11:
        POTFLAG = -11;
        break;
    case -12:
        POTFLAG = -12;
        break;
    case -13:
        POTFLAG = -13;
        break;
    default:
        POTFLAG = 12;
        break;
    }

    CATSparameters *cPars = NULL;
    CATSparameters *pPars = NULL;

    CATSparameters *cPotParsI1S0 = NULL;
    CATSparameters *cPotParsI1S1 = NULL;

    if(SOURCE == "NULL" || SOURCE == ""){

    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pXi0;
    }

    if (POT == "pXim_Lattice")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9] = {pXim_Lattice, 6, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_Lattice, 6, 1, 1, 1, 1, 0, 1, 0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if (POT == "pXim_HALQCD1")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9] = {pXim_HALQCD1, POTFLAG, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_HALQCD1, POTFLAG, 1, 1, 1, 1, 0, 1, 0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else if (POT == "pXim_HALQCDPaper2020")
    {
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotParsI1S0[9] = {pXim_HALQCDPaper2020, POTFLAG, 1, 1, 1, 0, 0, 0, 0};
        double PotParsI1S1[9] = {pXim_HALQCDPaper2020, POTFLAG, 1, 1, 1, 1, 0, 1, 0};
        cPotParsI1S0 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S0->SetParameters(PotParsI1S0);
        cPotParsI1S1 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotParsI1S1->SetParameters(PotParsI1S1);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pXi potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pXi0;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3322);
    Kitty.SetRedMass((Mass_p * Mass_Xi0) / (Mass_p + Mass_Xi0));

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1. / 4.);
    Kitty.SetChannelWeight(1, 3. / 4.);

    if (cPotParsI1S0)
        Kitty.SetShortRangePotential(0, 0, fDlmPot, *cPotParsI1S0);
    if (cPotParsI1S1)
        Kitty.SetShortRangePotential(1, 0, fDlmPot, *cPotParsI1S1);

CLEAN_SetUpCats_pXi0:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotParsI1S0)
    {
        delete cPotParsI1S0;
        cPotParsI1S0 = NULL;
    }
    if (cPotParsI1S1)
    {
        delete cPotParsI1S1;
        cPotParsI1S1 = NULL;
    }
}

void DLM_CommonAnaFunctions::SetUpCats_pOmegam(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{
    CATSparameters *cPars = NULL;
    CATSparameters *pPars = NULL;

    CATSparameters *cPotPars3S1 = NULL;
    CATSparameters *cPotPars5S2 = NULL;

    if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.8) * 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.8 * 1.2);
        cPars->SetParameter(1, 1.8);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[3].InitStability(20, 1, 2);
        CleverLevy[3].InitScale(35, 0.25, 2.0);
        CleverLevy[3].InitRad(256, 0, 64);
        CleverLevy[3].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.8);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[3].InitStability(20, 1, 2);
        CleverLevy[3].InitScale(35, 0.25, 2.0);
        CleverLevy[3].InitRad(256, 0, 64);
        CleverLevy[3].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetAnaSource(0, sqrt(1.8) * 1.2);
        Kitty.SetAnaSource(1, 1.8);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[3].InitStability(20, 1, 2);
        CleverLevy[3].InitScale(35, 0.25, 2.0);
        CleverLevy[3].InitRad(256, 0, 64);
        CleverLevy[3].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[3], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.8 * 1.2);
        Kitty.SetAnaSource(1, 1.8);
    }
    else if (SOURCE == "GaussExpTotSimple_2body")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (GaussExpTotSimple_2body)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if (SOURCE == "McLevyNolan_Reso")
    {
        CleverMcLevyReso[3].InitStability(21, 1, 2);
        CleverMcLevyReso[3].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[3].InitRad(257, 0, 64);
        CleverMcLevyReso[3].InitType(2);
        CleverMcLevyReso[3].InitReso(0, 1);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[3].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
        }
        else
        {
            CleverMcLevyReso[3].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
        }
        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[3].SetUpResoEmission(0, 0, HISTO);
        }
        CleverMcLevyReso[3].InitNumMcIter(200000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[3], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[3].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[3].InitScale(38, 0.15, 2.0);
        CleverMcLevyReso[3].InitRad(257 * 2, 0, 64);
        CleverMcLevyReso[3].InitType(2);
        CleverMcLevyReso[3].InitReso(0, 1);
        CleverMcLevyReso[3].InitReso(1, 1);
        if (SourceVar == 0)
        {
            CleverMcLevyReso[3].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards); // rdtBackwards,rdtRandom
            CleverMcLevyReso[3].SetUpReso(1, 0, 0.0001, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtBackwards);
        }
        else
        {
            CleverMcLevyReso[3].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom); // rdtBackwards,rdtRandom
            CleverMcLevyReso[3].SetUpReso(1, 0, 0.0001, 1361.52, 1.65, Mass_p, Mass_pic, false, false, DLM_CleverMcLevyReso::rdtRandom);
        }

        if (SourceVar == 2)
        {
            DLM_Histo<double> *HISTO = ConvertThetaAngleHisto("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/MixedEvents/AngleStudy_1/DimiPhi_pp.root", "h_rkAngle_Mom2", 400, 700);
            CleverMcLevyReso[3].SetUpResoEmission(0, 0, HISTO);
            CleverMcLevyReso[3].SetUpResoEmission(1, 0, HISTO);
        }

        CleverMcLevyReso[3].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[3], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    // SourceVar last digit is 0-9 the type
    //(SourceVar/10)*10 is the cutoff value (e.g. 192 is cutoff value of 190 and type 2)
    else if (SOURCE == "McGauss_ResoTM")
    {
        CleverMcLevyResoTM[3].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyResoTM[3].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[3].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[3].InitType(2);
        CleverMcLevyResoTM[3].SetUpReso(0, 0.6422);
        // pure Gauss
        if (SourceVar % 100 == 0)
        {
        }
        // back-to-back
        else if (SourceVar % 100 == 1)
        {
            CleverMcLevyResoTM[3].AddBGT_RP(490. / 1362. * 1.65, -1.);
        }
        // EPOS, 2 is with fixed mass, 3 is with EPOS mass
        else
        {
            const double k_CutOff = int(int(SourceVar) / 10) * 10.;
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

            TFile *F_EposDisto_pReso_Omega = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_Omega.root");
            TNtuple *T_EposDisto_pReso_Omega = (TNtuple *)F_EposDisto_pReso_Omega->Get("InfoTuple_ClosePairs");
            unsigned N_EposDisto_pReso_Omega = T_EposDisto_pReso_Omega->GetEntries();
            T_EposDisto_pReso_Omega->SetBranchAddress("k_D", &k_D);
            T_EposDisto_pReso_Omega->SetBranchAddress("P1", &fP1);
            T_EposDisto_pReso_Omega->SetBranchAddress("P2", &fP2);
            T_EposDisto_pReso_Omega->SetBranchAddress("M1", &fM1);
            T_EposDisto_pReso_Omega->SetBranchAddress("M2", &fM2);
            T_EposDisto_pReso_Omega->SetBranchAddress("Tau1", &Tau1);
            T_EposDisto_pReso_Omega->SetBranchAddress("Tau2", &Tau2);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP1", &AngleRcP1);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP2", &AngleRcP2);
            T_EposDisto_pReso_Omega->SetBranchAddress("AngleP1P2", &AngleP1P2);
            for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Omega; uEntry++)
            {
                T_EposDisto_pReso_Omega->GetEntry(uEntry);
                Tau1 = 1.65;
                Tau2 = 0;
                if (SourceVar % 100 == 2)
                {
                    fM1 = 1362;
                }
                if (k_D > k_CutOff)
                    continue;
                RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
                CleverMcLevyResoTM[3].AddBGT_RP(RanVal1, cos(AngleRcP1));
            }
            delete F_EposDisto_pReso_Omega;
        }

        CleverMcLevyResoTM[3].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[3], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "EPOS")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOS)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if (SOURCE == "EPOSrescaled")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (EPOSrescaled)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else if (SOURCE == "Levy_mT_Reso")
    {
        printf("\033[1;33mWARNING:\033[0m The CommonAnaFunction is still under construction (Levy_mT_Reso)\n");
        goto CLEAN_SetUpCats_pOmegam;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pOmegam;
    }

    if (POT.Contains("pOmega_Lattice"))
    {
        int FLAG = 12;
        if (POT.EndsWith("_11"))
            FLAG = 11;
        if (POT.EndsWith("_12"))
            FLAG = 12;
        if (POT.EndsWith("_13"))
            FLAG = 13;
        if (POT.EndsWith("_14"))
            FLAG = 14;
        if (POT.EndsWith("_121"))
            FLAG = 121;
        if (POT.EndsWith("_122"))
            FLAG = 122;
        if (POT.EndsWith("_123"))
            FLAG = 123;
        // #,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
        double PotPars3S1[3] = {100000, 0.8, 0.001};
        double PotPars5S2[9] = {pOmega_Lattice, double(FLAG), 0, 0, 0, 0, 0, 0, 0};
        cPotPars3S1 = new CATSparameters(CATSparameters::tPotential, 3, true);
        cPotPars3S1->SetParameters(PotPars3S1);
        cPotPars5S2 = new CATSparameters(CATSparameters::tPotential, 9, true);
        cPotPars5S2->SetParameters(PotPars5S2);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pOmega potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pOmegam;
    }
    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3334);
    Kitty.SetRedMass((Mass_p * MassOmega) / (Mass_p + MassOmega));

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);

    Kitty.SetSpin(0, 1);
    Kitty.SetSpin(1, 2);

    Kitty.SetChannelWeight(0, 3. / 8.);
    Kitty.SetChannelWeight(1, 5. / 8.);

    if (cPotPars3S1)
        Kitty.SetShortRangePotential(0, 0, RepulsiveCore, *cPotPars3S1);
    if (cPotPars5S2)
        Kitty.SetShortRangePotential(1, 0, fDlmPot, *cPotPars5S2);

CLEAN_SetUpCats_pOmegam:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotPars3S1)
    {
        delete cPotPars3S1;
        cPotPars3S1 = NULL;
    }
    if (cPotPars5S2)
    {
        delete cPotPars5S2;
        cPotPars5S2 = NULL;
    }
}

// Xi - Kaon ONLY COULOMB
void DLM_CommonAnaFunctions::SetUpCats_XiKCoulomb(CATS &Kitty, const TString &POT, const TString &SOURCE,
                                                  const TString &DataSample)
{

    CATSparameters *cPars = NULL;
    double radius;
    if (DataSample == "pp13TeV_MB_BBar")
        radius = 1.188;
    if (DataSample == "pp13TeV_HM_BBar")
        radius = 1.28;

    if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, radius);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "McGauss_Reso")
    {
        CleverMcLevyReso[1].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        CleverMcLevyReso[1].InitScale(35, 0.25, 2.0);
        CleverMcLevyReso[1].InitRad(256, 0, 64);
        CleverMcLevyReso[1].InitType(2);
        CleverMcLevyReso[1].InitReso(0, 1);
        CleverMcLevyReso[1].InitReso(1, 1);
        CleverMcLevyReso[1].SetUpReso(0, 0, 1. - 0.3578, 1361.52, 1.65, Mass_p, Mass_pic);
        CleverMcLevyReso[1].SetUpReso(1, 0, 1. - 0.3562, 1462.93, 4.69, Mass_L, Mass_pic);
        CleverMcLevyReso[1].InitNumMcIter(1000000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyReso[1], 2);
        Kitty.SetAnaSource(0, 0.87);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
        printf("SetUpCats_XiKCoulomb should NOT use this source set up!!!\n");
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pApCoulomb;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, -2212);
    Kitty.SetRedMass((Mass_p * Mass_p) / (Mass_p + Mass_p));

    Kitty.SetNumChannels(2);
    Kitty.SetNumPW(0, 1);
    Kitty.SetNumPW(1, 1);
    Kitty.SetSpin(0, 0);
    Kitty.SetSpin(1, 1);
    Kitty.SetChannelWeight(0, 1. / 4.);
    Kitty.SetChannelWeight(1, 3. / 4.);

CLEAN_SetUpCats_pApCoulomb:;
}

void DLM_CommonAnaFunctions::SetUpCats_LKVidana(CATS &Kitty, const TString &SOURCE, const TString &DataSample)
{
    CATSparameters *cPars = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    unsigned NumChannels = 6;

    double rad1;
    double rad2;
    double w_source;
    double lam_source;
    if (DataSample == "pp13TeV_HM_LKMAC")
    {
        rad1 = 1.202;
        rad2 = 2.33;
        w_source = 0.7993;
        lam_source = 0.9806;
        ExternalWF = Init_LAntiK_Vidana("/Users/sartozza/cernbox/CATSFiles/Models_WFs_Theoreticians/Vidana/", Kitty, 0);
    }
    else if (DataSample == "pp13TeV_HM_LK")
    {
        rad1 = 1.202;
        rad2 = 2.33;
        w_source = 0.7993;
        lam_source = 0.9806;
        ExternalWF = Init_LAntiK_Vidana("/home/valentina/thor/cernbox/CATSFiles/Models_WFs_Theoreticians/Vidana/", Kitty, 0);
    }
    if (SOURCE == "Gauss")
    {
        cout << "Gaussian source" << endl;
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, rad1);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "DoubleGauss")
    {
        cout << "Double Gaussian source" << endl;
        cPars = new CATSparameters(CATSparameters::tSource, 4, true);
        cPars->SetParameter(0, rad1);
        cPars->SetParameter(1, rad2);
        cPars->SetParameter(2, w_source);
        cPars->SetParameter(3, lam_source);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(NormDoubleGaussSource, *cPars);
        Kitty.SetAutoNormSource(false);  // MUST ALWAYS BE FALSE!!
        Kitty.SetNormalizedSource(true); // do not touch the source, set to true so CATS adds (1-λs), if false it changes the distribution
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pApHaide;
    }

    NumChannels = 6;

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(3122, -321);
    Kitty.SetRedMass((Mass_L * Mass_Kch) / (Mass_L + Mass_Kch));

    Kitty.SetNumChannels(NumChannels);

    for (unsigned uCh = 0; uCh < NumChannels; uCh++)
    {
        if (ExternalWF)
        {
            Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m SetUpCats_LKVidana says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pApHaide;
        }

    } // end of for

CLEAN_SetUpCats_pApHaide:;

    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    CleanUpWfHisto(Kitty, ExternalWF);
}

void DLM_CommonAnaFunctions::SetUpCats_Lcp_LQCDe(CATS &Kitty, const TString &SOURCE, const int &CUTOFF)
{
    CATSparameters *cPars = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    unsigned NumChannels = 8;

    double rad1;
    double rad2;
    double w_source;
    double lam_source;

    rad1 = 1.2;
    ExternalWF = Init_Lcp_Haidenbauer("/Users/sartozza/cernbox/CATSFiles/Models_WFs_Theoreticians/Johann/Lc_Proton/", Kitty, CUTOFF);

    if (SOURCE == "Gauss")
    {
        cout << "Gaussian source" << endl;
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, rad1);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "DoubleGauss")
    {
        cout << "Double Gaussian source" << endl;
        cPars = new CATSparameters(CATSparameters::tSource, 4, true);
        cPars->SetParameter(0, rad1);
        cPars->SetParameter(1, rad2);
        cPars->SetParameter(2, w_source);
        cPars->SetParameter(3, lam_source);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(NormDoubleGaussSource, *cPars);
        Kitty.SetAutoNormSource(false);  // MUST ALWAYS BE FALSE!!
        Kitty.SetNormalizedSource(true); // do not touch the source, set to true so CATS adds (1-λs), if false it changes the distribution
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pApHaide;
    }

    Kitty.SetMomentumDependentSource(false);
    Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(4122, 2212);
    Kitty.SetRedMass((Mass_Lcp * Mass_p) / (Mass_Lcp + Mass_p));

    Kitty.SetNumChannels(NumChannels);

    for (unsigned uCh = 0; uCh < NumChannels; uCh++)
    {
        if (ExternalWF)
        {
            Kitty.SetExternalWaveFunction(uCh, 0, ExternalWF[0][uCh][0], ExternalWF[1][uCh][0]);
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m SetUpCats_Lcp says that you should NEVER see this message! BIG BUG!\n");
            goto CLEAN_SetUpCats_pApHaide;
        }

    } // end of for

CLEAN_SetUpCats_pApHaide:;

    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    CleanUpWfHisto(Kitty, ExternalWF);
}

DLM_Ck *DLM_CommonAnaFunctions::SetUpLednicky_pL(const unsigned &NumMomBins, const double *MomBins, const TString &POT)
{
    double p0, p1, p2, p3;
    if (POT == "Lednicky_ND")
    {
        p0 = 1.77, p2 = 2.06, p1 = 3.78, p3 = 3.18;
    }
    else if (POT == "Lednicky_NF")
    {
        p0 = 2.18, p2 = 1.93, p1 = 3.19, p3 = 3.358;
    }
    else if (POT == "Lednicky_NSC89")
    {
        p0 = 2.73, p2 = 1.48, p1 = 2.87, p3 = 3.04;
    }
    else if (POT == "Lednicky_NSC97a")
    {
        p0 = 0.71, p2 = 2.18, p1 = 5.86, p3 = 2.76;
    }
    else if (POT == "Lednicky_NSC97b")
    {
        p0 = 0.9, p2 = 2.13, p1 = 4.92, p3 = 2.84;
    }
    else if (POT == "Lednicky_NSC97c")
    {
        p0 = 1.2, p2 = 2.08, p1 = 4.11, p3 = 2.92;
    }
    else if (POT == "Lednicky_NSC97d")
    {
        p0 = 1.71, p2 = 1.95, p1 = 3.46, p3 = 3.08;
    }
    else if (POT == "Lednicky_NSC97e")
    {
        p0 = 2.1, p2 = 1.86, p1 = 3.19, p3 = 3.19;
    }
    else if (POT == "Lednicky_NSC97f")
    {
        p0 = 2.51, p2 = 1.75, p1 = 3.03, p3 = 3.32;
    }
    else if (POT == "Lednicky_ESC08")
    {
        p0 = 2.7, p2 = 1.65, p1 = 2.97, p3 = 3.63;
    }
    else if (POT == "Lednicky_XeftLO")
    {
        p0 = 1.91, p2 = 1.23, p1 = 1.4, p3 = 2.13;
    }
    else if (POT == "Lednicky_XeftNLO")
    {
        p0 = 2.91, p2 = 1.54, p1 = 2.78, p3 = 2.72;
    }
    else if (POT == "Lednicky_JulichA")
    {
        p0 = 1.56, p2 = 1.59, p1 = 1.43, p3 = 3.16;
    }
    else if (POT == "Lednicky_JulichJ04")
    {
        p0 = 2.56, p2 = 1.66, p1 = 2.75, p3 = 2.93;
    }
    else if (POT == "Lednicky_JulichJ04c")
    {
        p0 = 2.66, p2 = 1.57, p1 = 2.67, p3 = 3.08;
    }
    else
    {
        printf("\033[1;31mERROR (SetUpLednicky_pL):\033[0m The pΛ potential '%s' does not exist\n", POT.Data());
    }
    DLM_Ck *DlmCk = new DLM_Ck(1, 4, NumMomBins, MomBins, Lednicky_SingletTriplet); //
    DlmCk->SetPotPar(0, p0);
    DlmCk->SetPotPar(1, p1);
    DlmCk->SetPotPar(2, p2);
    DlmCk->SetPotPar(3, p3);
    return DlmCk;
}

void DLM_CommonAnaFunctions::SetUpBinning_pp(const TString &DataSample, unsigned &NumMomBins, double *&MomBins, double *&FitRegion, const int &MomBinVar, const int &FitRegVar)
{
    double kMin;
    double kMax;
    double kStep;
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        if (MomBinVar == 0)
        {
            kMin = 0;
            kStep = 4;
            kMax = 376; //(i.e. max=376 MeV)
        }
        else if (MomBinVar == 1)
        {
            kMin = 0;
            kStep = 4;
            kMax = 352;
        }
        else if (MomBinVar == 2)
        {
            kMin = 0;
            kStep = 4;
            kMax = 400;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        if (MomBinVar == 0)
        {
            kMin = 4;
            kStep = 4;
            kMax = 376; //(i.e. max=376 MeV)
        }
        else if (MomBinVar == 1)
        {
            kMin = 4;
            kStep = 4;
            kMax = 352;
        }
        else if (MomBinVar == 2)
        {
            kMin = 4;
            kStep = 4;
            kMax = 400;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        if (MomBinVar == 0)
        {
            kMin = 4;
            kStep = 4;
            kMax = 376; //(i.e. max=376 MeV)
        }
        else if (MomBinVar == 1)
        {
            kMin = 4;
            kStep = 4;
            kMax = 352;
        }
        else if (MomBinVar == 2)
        {
            kMin = 4;
            kStep = 4;
            kMax = 400;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        if (MomBinVar == 0)
        {
            kMin = 4;
            kStep = 4;
            kMax = 376; //(i.e. max=376 MeV)
        }
        else if (MomBinVar == 1)
        {
            kMin = 4;
            kStep = 4;
            kMax = 352;
        }
        else if (MomBinVar == 2)
        {
            kMin = 4;
            kStep = 4;
            kMax = 400;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }
    }
    else
    {
        printf("\033[1;31mERROR (SetUpBinning_pp):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        NumMomBins = 0;
        return;
    }
    NumMomBins = floor((kMax - kMin) / (kStep));

    if (MomBins)
        delete[] MomBins;
    MomBins = new double[NumMomBins + 1];
    MomBins[0] = kMin;
    for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
    {
        MomBins[uBin] = MomBins[uBin - 1] + kStep;
    }
    if (FitRegion)
        delete[] FitRegion;
    FitRegion = new double[4];
    FitRegion[0] = MomBins[0];
    FitRegion[1] = MomBins[NumMomBins];
    FitRegion[2] = MomBins[NumMomBins] + kStep;
    FitRegion[3] = MomBins[NumMomBins] + kStep * 31.; // till 500

    // if(FitRegVar==0){
    //     FitRegion[0] = MomBins[0];
    //     FitRegion[1] = MomBins[NumMomBins];
    //     FitRegion[2] = MomBins[NumMomBins];
    //     FitRegion[3] = MomBins[NumMomBins];
    // }
    if (FitRegVar == 0)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 576;
    }
    else if (FitRegVar == 1)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 552;
    }
    else if (FitRegVar == 2)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 600;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m The FitRegVar '%i' does not exist\n", FitRegVar);
        return;
    }
}
void DLM_CommonAnaFunctions::SetUpBinning_pL(const TString &DataSample, unsigned &NumMomBins, double *&MomBins, double *&FitRegion,
                                             const int &MomBinVar, const int &FitRegVar)
{

    double kMin;
    double kFineMin;
    double kFineMax;
    double kMax;
    double kCoarseStep;
    double kFineStep;

    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        if (MomBinVar == 0)
        {
            kMin = 0;
            kFineMin = 336; // 272//216
            kFineMax = 336; // 304
            kMax = 336;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 1)
        {
            kMin = 0;
            kFineMin = 312; // 272//216
            kFineMax = 312; // 304
            kMax = 312;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 2)
        {
            kMin = 0;
            kFineMin = 348; // 272//216
            kFineMax = 348; // 304
            kMax = 348;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }

        // the number of coarse bins below kFineMin
        // floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        // and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin - kMin) / kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax - double(NumCoarseBinsBelow) * kCoarseStep) / kFineStep);
        // we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax - double(NumCoarseBinsBelow) * kCoarseStep - double(NumFineBins) * kFineStep) / kCoarseStep);

        NumMomBins = NumCoarseBinsBelow + NumFineBins + NumCoarseBinsAbove;

        if (MomBins)
            delete[] MomBins;
        MomBins = new double[NumMomBins + 1];
        MomBins[0] = kMin;
        for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
        {
            if (uBin <= NumCoarseBinsBelow || uBin > NumCoarseBinsBelow + NumFineBins)
            {
                MomBins[uBin] = MomBins[uBin - 1] + kCoarseStep;
            }
            else
            {
                MomBins[uBin] = MomBins[uBin - 1] + kFineStep;
            }
        }
        if (FitRegion)
            delete[] FitRegion;
        FitRegion = new double[4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins] + kCoarseStep;       // 348
        FitRegion[3] = MomBins[NumMomBins] + kCoarseStep * 20.; // 588
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        if (MomBinVar == 0)
        {
            kMin = 0;
            kFineMin = 336; // 272//216
            kFineMax = 336; // 304
            kMax = 336;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 1)
        {
            kMin = 0;
            kFineMin = 312; // 272//216
            kFineMax = 312; // 304
            kMax = 312;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 2)
        {
            kMin = 0;
            kFineMin = 348; // 272//216
            kFineMax = 348; // 304
            kMax = 348;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 10)
        {
            kMin = 0;
            kFineMin = 204; // 272//216
            kFineMax = 204; // 304
            kMax = 204;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 11)
        {
            kMin = 0;
            kFineMin = 180; // 272//216
            kFineMax = 180; // 304
            kMax = 180;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else if (MomBinVar == 12)
        {
            kMin = 0;
            kFineMin = 228; // 272//216
            kFineMax = 228; // 304
            kMax = 228;     // 336
            kCoarseStep = 12;
            kFineStep = 12;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }

        // the number of coarse bins below kFineMin
        // floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        // and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin - kMin) / kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax - double(NumCoarseBinsBelow) * kCoarseStep) / kFineStep);
        // we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax - double(NumCoarseBinsBelow) * kCoarseStep - double(NumFineBins) * kFineStep) / kCoarseStep);

        NumMomBins = NumCoarseBinsBelow + NumFineBins + NumCoarseBinsAbove;

        if (MomBins)
            delete[] MomBins;
        MomBins = new double[NumMomBins + 1];
        MomBins[0] = kMin;
        for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
        {
            if (uBin <= NumCoarseBinsBelow || uBin > NumCoarseBinsBelow + NumFineBins)
            {
                MomBins[uBin] = MomBins[uBin - 1] + kCoarseStep;
            }
            else
            {
                MomBins[uBin] = MomBins[uBin - 1] + kFineStep;
            }
        }
        if (FitRegion)
            delete[] FitRegion;
        FitRegion = new double[4];

        if (FitRegVar == 0)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 576;
        }
        else if (FitRegVar == 1)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 552;
        }
        else if (FitRegVar == 2)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 600;
        }
        else if (FitRegVar == 10)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 360; // 420
            FitRegion[3] = 576;
        }
        else if (FitRegVar == 11)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 336; // 432
            FitRegion[3] = 552;
        }
        else if (FitRegVar == 12)
        {
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = 384; // 408
            FitRegion[3] = 600;
        }
        // for the Schleching fits
        else if (FitRegVar == 50)
        { // only femto range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = MomBins[NumMomBins];
        }
        // only used for pol3
        else if (FitRegVar == 51)
        { // extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 456;
        }
        else if (FitRegVar == 52)
        { // extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 432;
        }
        else if (FitRegVar == 53)
        { // extended fit range
            FitRegion[0] = MomBins[0];
            FitRegion[1] = MomBins[NumMomBins];
            FitRegion[2] = MomBins[NumMomBins]; // 336
            FitRegion[3] = 480;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The FitRegVar '%i' does not exist\n", FitRegVar);
            return;
        }
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        kMin = 0;
        kFineMin = 336;
        kFineMax = 336;
        kMax = 336;
        kCoarseStep = 12;
        kFineStep = 12;

        // the number of coarse bins below kFineMin
        // floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        // and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin - kMin) / kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax - double(NumCoarseBinsBelow) * kCoarseStep) / kFineStep);
        // we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax - double(NumCoarseBinsBelow) * kCoarseStep - double(NumFineBins) * kFineStep) / kCoarseStep);

        NumMomBins = NumCoarseBinsBelow + NumFineBins + NumCoarseBinsAbove;

        if (MomBins)
            delete[] MomBins;
        MomBins = new double[NumMomBins + 1];
        MomBins[0] = kMin;
        for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
        {
            if (uBin <= NumCoarseBinsBelow || uBin > NumCoarseBinsBelow + NumFineBins)
            {
                MomBins[uBin] = MomBins[uBin - 1] + kCoarseStep;
            }
            else
            {
                MomBins[uBin] = MomBins[uBin - 1] + kFineStep;
            }
        }
        if (FitRegion)
            delete[] FitRegion;
        FitRegion = new double[4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins] + kCoarseStep;       // 348
        FitRegion[3] = MomBins[NumMomBins] + kCoarseStep * 20.; // 588
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        kMin = 0;
        kFineMin = 336;
        kFineMax = 336;
        kMax = 336;
        kCoarseStep = 12;
        kFineStep = 12;

        // the number of coarse bins below kFineMin
        // floor/ceil combination makes sure that we include the WHOLE region we want in the fine binning,
        // and if there is rounding needed it is done so that we make our region larger, not smaller!
        unsigned NumCoarseBinsBelow = floor((kFineMin - kMin) / kCoarseStep);
        unsigned NumFineBins = ceil((kFineMax - double(NumCoarseBinsBelow) * kCoarseStep) / kFineStep);
        // we floor the highest point, to make sure we do not run out of the range provided by experimental data
        unsigned NumCoarseBinsAbove = floor((kMax - double(NumCoarseBinsBelow) * kCoarseStep - double(NumFineBins) * kFineStep) / kCoarseStep);

        NumMomBins = NumCoarseBinsBelow + NumFineBins + NumCoarseBinsAbove;

        if (MomBins)
            delete[] MomBins;
        MomBins = new double[NumMomBins + 1];
        MomBins[0] = kMin;
        for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
        {
            if (uBin <= NumCoarseBinsBelow || uBin > NumCoarseBinsBelow + NumFineBins)
            {
                MomBins[uBin] = MomBins[uBin - 1] + kCoarseStep;
            }
            else
            {
                MomBins[uBin] = MomBins[uBin - 1] + kFineStep;
            }
        }
        if (FitRegion)
            delete[] FitRegion;
        FitRegion = new double[4];
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins] + kCoarseStep;       // 348
        FitRegion[3] = MomBins[NumMomBins] + kCoarseStep * 18.; // 564
    }
    else
    {
        printf("\033[1;31mERROR (SetUpBinning_pL):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        NumMomBins = 0;
        return;
    }
}

void DLM_CommonAnaFunctions::SetUpBinning_LKmin(const TString &DataSample, unsigned &NumMomBins, double *&MomBins, double *&FitRegion, const int &MomBinVar, const int &FitRegVar)
{
    double kMin;
    double kMax;
    double kStep;
    if (DataSample == "pp13TeV_HM_LK")
    {
        if (MomBinVar == 0)
        {
            kMin = 0;
            kStep = 4;
            kMax = 362; //(i.e. max=376 MeV)
        }
        else if (MomBinVar == 1)
        {
            kMin = 0;
            kStep = 4;
            kMax = 602;
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The MomBinVar '%i' does not exist\n", MomBinVar);
            return;
        }
    }
    else
    {
        printf("\033[1;31mERROR (SetUpBinning_pp):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        NumMomBins = 0;
        return;
    }
    NumMomBins = floor((kMax - kMin) / (kStep));

    if (MomBins)
        delete[] MomBins;
    MomBins = new double[NumMomBins + 1];
    MomBins[0] = kMin;
    for (unsigned uBin = 1; uBin <= NumMomBins; uBin++)
    {
        MomBins[uBin] = MomBins[uBin - 1] + kStep;
    }
    if (FitRegion)
        delete[] FitRegion;
    FitRegion = new double[4];
    FitRegion[0] = MomBins[0];
    FitRegion[1] = MomBins[NumMomBins];
    FitRegion[2] = MomBins[NumMomBins] + kStep;
    FitRegion[3] = MomBins[NumMomBins] + kStep * 31.; // till 500

    // if(FitRegVar==0){
    //     FitRegion[0] = MomBins[0];
    //     FitRegion[1] = MomBins[NumMomBins];
    //     FitRegion[2] = MomBins[NumMomBins];
    //     FitRegion[3] = MomBins[NumMomBins];
    // }
    if (FitRegVar == 0)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 576;
    }
    else if (FitRegVar == 1)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 552;
    }
    else if (FitRegVar == 2)
    {
        FitRegion[0] = MomBins[0];
        FitRegion[1] = MomBins[NumMomBins];
        FitRegion[2] = MomBins[NumMomBins];
        FitRegion[3] = 600;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m The FitRegVar '%i' does not exist\n", FitRegVar);
        return;
    }
}

void DLM_CommonAnaFunctions::GetPurities_p(const TString &DataSample, const int &Variation, double *Purities)
{
    double PurityProton;
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        PurityProton = 0.989859;
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        PurityProton = 0.9943;
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        PurityProton = 0.984265;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        PurityProton = 0.984265;
    }
    else if (DataSample == "pp13TeV_HM_BernieSource")
    {
        PurityProton = 0.9943;
    }
    else
    {
        printf("\033[1;31mERROR (GetPurities_p):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        PurityProton = 1.0;
    }

    // following my lambda pars with the 3 possible modifications
    // for the proton:
    // 0 = primary
    // 1 = from Lambda
    // 2 = other feeddown (flat)
    // 3 = missidentified
    // const unsigned NumChannels_p = 4;
    // if(Purities){delete [] Purities; Purities = new double [NumChannels_p];}
    Purities[0] = PurityProton;
    Purities[1] = PurityProton;
    Purities[2] = PurityProton;
    Purities[3] = 1. - PurityProton;
}

void DLM_CommonAnaFunctions::GetPurities_L(const TString &DataSample, const int &Variation, double *Purities)
{
    double PurityLambda;
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        PurityLambda = 0.96768;
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19")
    {
        // printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 is not available yet!\n");
        if (Variation == 0)
            PurityLambda = 0.9595; // the original value for the preliminaries
        else if (Variation == 1)
            PurityLambda = 0.936; // spline fits 4th June 2020
        else if (Variation == 2)
            PurityLambda = 0.936 - 0.006; // with uncertainties
        else if (Variation == 3)
            PurityLambda = 0.936 + 0.006; // with uncertainties
        else if (Variation == 4)
            PurityLambda = 0.953; // the purity as published in PLB (lower limit)
        else if (Variation == 5)
            PurityLambda = 0.963; // the purity as published in PLB (upper limit)
        else if (Variation == -1)
            PurityLambda = 1.0; // use for SB corrected correlations
        else
            PurityLambda = 0.9595;
    }
    else if (DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        PurityLambda = 1.0; // use for SB corrected correlations
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        PurityLambda = 0.937761;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        PurityLambda = 0.937761;
    }
    else
    {
        printf("\033[1;31mERROR (GetPurities_L):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        PurityLambda = 1.0;
    }

    // for the Lambda:
    // 0 = primary
    // 1 = from Sigma0
    // 2 = from Xim
    // 3 = from Xi0
    // 4 = missidentified
    Purities[0] = PurityLambda;
    Purities[1] = PurityLambda;
    Purities[2] = PurityLambda;
    Purities[3] = PurityLambda;
    Purities[4] = 1. - PurityLambda;
}

void DLM_CommonAnaFunctions::GetPurities_L_Vale(const TString &DataSample, const int &Variation, double *Purities, int SamplePurity)
{
    double PurityLambdaVars[45] = {0.9465, 0.965985, 0.955872, 0.966306, 0.956546,
                                   0.965578, 0.965118, 0.967649, 0.951028, 0.958578, 0.958875, 0.958578, 0.956233,

                                   0.958878, 0.953028, 0.969392, 0.956932, 0.96478, 0.954015, 0.967371, 0.953029,
                                   0.954015, 0.965445, 0.959869, 0.965097, 0.958929, 0.952171, 0.960708, 0.96445,
                                   0.956484, 0.969944, 0.965535, 0.965945, 0.964151, 0.953982, 0.96738, 0.961568,
                                   0.965563, 0.964181, 0.962324, 0.960303, 0.966329, 0.970294, 0.951154, 0.968579};
    double PurityLambda;
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        PurityLambda = 0.96768;
    }
    else if (DataSample == "pp13TeV_HM_LK")
    {
        PurityLambda = PurityLambdaVars[SamplePurity];
    }
    else
    {
        printf("\033[1;31mERROR (GetPurities_L):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        PurityLambda = 1.0;
    }

    // for the Lambda:
    // 0 = primary
    // 1 = from Sigma0
    // 2 = from Xim
    // 3 = from Xi0
    // 4 = missidentified
    Purities[0] = PurityLambda;
    Purities[1] = PurityLambda;
    Purities[2] = PurityLambda;
    Purities[3] = PurityLambda;
    Purities[4] = 1. - PurityLambda;
}

void DLM_CommonAnaFunctions::GetPurities_Xim(const TString &DataSample, const int &Variation, double *Purities)
{
    double PurityXim;
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        PurityXim = 0.956;
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        PurityXim = 0.956;
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        PurityXim = 0.88;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        PurityXim = 0.88;
    }
    else
    {
        printf("\033[1;31mERROR (GetPurities_Xim):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        PurityXim = 1.0;
    }

    // 0 = primary
    // 1 = from Xi-(1530)
    // 2 = from Xi0(1530)
    // 3 = from Omega
    // 4 = missidentified
    Purities[0] = PurityXim;
    Purities[1] = PurityXim;
    Purities[2] = PurityXim;
    Purities[3] = PurityXim;
    Purities[4] = 1. - PurityXim;
}

/// PurityKaon[0] os the default pT averaged
void DLM_CommonAnaFunctions::GetPurities_K(const TString &DataSample, const int &Variation, double *Purities, int SamplePurity)
{
    double PurityKaonVars[45] = {0.997177, 0.997259, 0.997129, 0.997177, 0.997093, 0.994909,
                                 0.99508, 0.994934, 0.997216, 0.997183, 0.997093, 0.994909,
                                 0.995998, 0.996075, 0.996066, 0.996075, 0.994954, 0.994928, 0.996075,
                                 0.996031, 0.997259, 0.994892, 0.997093, 0.994903, 0.994985, 0.997158,
                                 0.996066, 0.996079, 0.997259, 0.996009, 0.997098, 0.996042, 0.994996,
                                 0.997139, 0.997104, 0.994965, 0.997062, 0.996085, 0.997221, 0.99616,
                                 0.997152, 0.994909, 0.994997, 0.994934, 0.994934};
    double PurityKaon;

    if (DataSample == "pp13TeV_HM_LK")
    {
        PurityKaon = PurityKaonVars[SamplePurity];
    }
    else
    {
        printf("\033[1;31mERROR (GetPurities_K):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        PurityKaon = 1.0;
    }

    // following my lambda pars with the 3 possible modifications
    // for the kaon:
    // 0 = primary
    // 1 = from Phi
    // 2 = other feeddown (flat)
    // 3 = missidentified
    // const unsigned NumChannels_p = 4;
    // if(Purities){delete [] Purities; Purities = new double [NumChannels_p];}
    Purities[0] = PurityKaon;
    Purities[1] = PurityKaon;
    Purities[2] = PurityKaon;
    Purities[3] = 1. - PurityKaon;
}

// the last digit of Variation is generic.
// the second to-last digit referes to BernieSource
void DLM_CommonAnaFunctions::GetFractions_p(const TString &DataSample, const int &Variation, double *Fractions)
{
    // following my lambda pars with the 3 possible modifications
    // for the proton:
    // 0 = primary
    // 1 = from Lambda
    // 2 = other feeddown (flat)
    // 3 = missidentified
    double Modify_pp = 1;
    switch (Variation % 10)
    {
    case 1:
        Modify_pp = 0.8;
        break;
    case 2:
        Modify_pp = 1.2;
        break;
    default:
        Modify_pp = 1;
        break;
    }
    double pp_f0; // primary protons
    double pp_f1; // fraction of Lambda
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        pp_f0 = 0.87397;
        pp_f1 = 0.0882211;
    }
    else if (DataSample == "pp13TeV_HM_March19")
    {
        printf("\033[1;33mWARNING:\033[0m pp13TeV_HM_March19 CROSS CHECK pp_f0 and pp_f1!\n");
        pp_f0 = 0.873;
        pp_f1 = 0.0898;
    }
    else if (DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        pp_f0 = 0.823;
        pp_f1 = 0.125;
    }
    // Variation reflects the mT bin
    else if (DataSample == "pp13TeV_HM_BernieSource")
    {
        std::vector<double> pp_primary = {0.82, 0.81, 0.81, 0.81, 0.81, 0.82, 0.83};
        pp_f0 = pp_primary.at((Variation / 10) % 10);
        pp_f1 = 0.7 * (1. - pp_f0);
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        pp_f0 = 0.862814;
        pp_f1 = 0.09603;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        pp_f0 = 0.862814;
        pp_f1 = 0.09603;
    }
    else
    {
        printf("\033[1;31mERROR (GetFractions_p):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        pp_f0 = 1.0;
        pp_f1 = 0.0;
    }
    // ratio between feed-down lambdas and all feed downs, it is by default 0.7, i.e. 70% of the feed-down into protons is from lambdas
    double arrayPercLamProton = pp_f1 / (1. - pp_f0) * Modify_pp;
    Fractions[0] = pp_f0;
    Fractions[1] = (1. - pp_f0) * (arrayPercLamProton);
    Fractions[2] = (1. - pp_f0) * (1. - arrayPercLamProton);
    Fractions[3] = 1.;
}
// Variation -> use the two digits
// first digit->Modify_SigL
// second digit->Modify_XiL+3*Modify_PrimFrac
// 0 -> default; 1 = -20%; 2 = +20%
void DLM_CommonAnaFunctions::GetFractions_L(const TString &DataSample, const int &Variation, double *Fractions)
{
    double Modify_SigL = 1;
    double Modify_XiL = 1;
    // the amount of prim lambdas depends on pT, at low pT (and k*)
    // the fraction could be lower. You can modify it here. This is by how much the PrimLambda fractional yield is reduced,
    // compared to the value at the average pT
    double Modify_PrimFrac = 1;
    switch (Variation % 10)
    {
        // the new values will be 0.6,0.84,1.2,1.44
        // to keep consistency maybe just use 0.6,0.8,1.0,1.2,1.4
        // this will correspond to s/l ratio of 0.2,0.267,0.333,0.4,0.467
    case 0:
        Modify_SigL = 1;
        break;
    case 1:
        Modify_SigL = 0.8;
        break;
    case 2:
        Modify_SigL = 1.2;
        break;
    case 3:
        Modify_SigL = 0.6;
        break;
    case 4:
        Modify_SigL = 1.4;
        break;
    default:
        Modify_SigL = 1;
        break;
    }
    // this gives info about ratio of Xim to Xi0, absolute as we now assume it to be 50:50
    switch ((Variation / 10) % 3)
    {
    case 0:
        Modify_XiL = 1;
        break;
    case 1:
        Modify_XiL = 0.8;
        break;
    case 2:
        Modify_XiL = 1.2;
        break;
    default:
        Modify_XiL = 1;
        break;
    }
    // this only works assuming zero material
    // do not use for newer interations!
    // i.e. the 0.95 variation comes from the lower pT, which results in a bit different primary fraction
    switch ((Variation / 10) / 3)
    {
    case 0:
        Modify_PrimFrac = 1;
        break;
    case 1:
        Modify_PrimFrac = 0.95;
        break;
    default:
        Modify_PrimFrac = 1;
        break;
    }
    double pL_f0; // fraction of primary Lambdas
    double pL_f1; // fraction of Sigma0
    double pL_f2; // fractions of Xi0/m (each)
    double pL_fm; // fraction of material
    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        pL_f0 = 0.601008;
        pL_f1 = 0.200336;
        pL_f2 = 0.099328;
        pL_fm = 0;
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        pL_f0 = 0.576066;
        pL_f1 = 0.192022;
        pL_f2 = 0.115956;
        pL_fm = 0.0;
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        pL_f0 = 0.521433;
        pL_f1 = 0.173811;
        pL_f2 = 0.152378;
        pL_fm = 0.0;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        pL_f0 = 0.521433;
        pL_f1 = 0.173811;
        pL_f2 = 0.152378;
        pL_fm = 0.0;
    }
    else
    {
        printf("\033[1;31mERROR (GetFractions_L):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        pL_f0 = 1.0;
        pL_f1 = 0.0;
        pL_f2 = 0.0;
        pL_fm = 0.0;
    }

    // this only works assuming zero material
    // do not use for newer interations!
    if (pL_f2)
    {
        pL_f0 *= Modify_PrimFrac;
        pL_f1 *= Modify_PrimFrac;
        pL_f2 = (1. - pL_f0 - pL_f1) * 0.5;
    }

    double SigLambdaPrimDir = pL_f0 + pL_f1;
    // ratio between sigma0 feed down and primary lambdas. By default this should be 1:3
    double arrayPercSigLambda = pL_f1 / pL_f0 * Modify_SigL;
    // ration between xim feed down to the flat (xi0) feed. By default we assume it is 0.5
    double arrayPercXiLambda = pL_f2 / (1. - pL_f0 - pL_f1) * Modify_XiL;
    double FracOfLambda = 1. / (1. + arrayPercSigLambda);
    // 0 is primary
    // 1 is from Sigma0
    // 2 is is from Xim
    // 3 is is the flat feeddown (it should be just xi0, i.e. == xim)
    //   remark: up to Aug 2020, the material contribution, c.a. 1% was distributed amount ALL other
    //   entries, primary, feed etc. Actually it would make most sence to put it as an additional flat contribution!
    // 4 is for the missid
    Fractions[0] = SigLambdaPrimDir * FracOfLambda;
    Fractions[1] = SigLambdaPrimDir * (1. - FracOfLambda);
    Fractions[2] = (1. - SigLambdaPrimDir) * (arrayPercXiLambda);
    Fractions[3] = 1. - Fractions[0] - Fractions[1] - Fractions[2];
    Fractions[4] = 1.;
}

void DLM_CommonAnaFunctions::GetFractions_L_Vale(const TString &DataSample, const int &Variation, double *Fractions)
{
    double Modify_SigL = 1;
    double Modify_XiL = 1;
    // the amount of prim lambdas depends on pT, at low pT (and k*)
    // the fraction could be lower. You can modify it here. This is by how much the PrimLambda fractional yield is reduced,
    // compared to the value at the average pT
    double Modify_PrimFrac = 1;
    switch (Variation % 10)
    {
    case 0:
        Modify_SigL = 1.;
        break;
    case 1:
        Modify_SigL = 0.9;
        break;
    case 2:
        Modify_SigL = 1.1;
        break;
    default:
        Modify_SigL = 1;
        break;
    }
    switch (Variation / 10)
    {
    case 0:
        Modify_XiL = 1.;
        break;
    case 1:
        Modify_XiL = 0.9;
        break;
    case 2:
        Modify_XiL = 1.1;
        break;
    default:
        Modify_XiL = 1;
        break;
    }

    double pL_f0; // fraction of primary Lambdas
    double pL_f1; // fraction of Sigma0
    double pL_f2; // fractions of Xi0/m (each)
    double pL_fm; // fraction of material

    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        pL_f0 = 0.601008;
        pL_f1 = 0.200336;
        pL_f2 = 0.099328;
        pL_fm = 0;
    }
    else if (DataSample == "pp13TeV_HM_LK")
    {
        pL_f1 = 0.192;
        pL_f2 = 0.232 / 2.;
        pL_f0 = 1. - pL_f1 - pL_f2 * 2;
        pL_fm = 0.0;
    }
    else
    {
        printf("\033[1;31mERROR (GetFractions_L):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        pL_f0 = 1.0;
        pL_f1 = 0.0;
        pL_f2 = 0.0;
        pL_fm = 0.0;
    }

    double SigLambdaPrimDir = pL_f0 + pL_f1;
    // ratio between sigma0 feed down and primary lambdas. By default this should be 1:3
    double arrayPercSigLambda = pL_f1 / pL_f0 * Modify_SigL;
    // ration between xim feed down to the flat (xi0) feed. By default we assume it is 0.5
    double arrayPercXiLambda = pL_f2 / (1. - pL_f0 - pL_f1) * Modify_XiL;
    double FracOfLambda = 1. / (1. + arrayPercSigLambda);
    // 0 is primary
    // 1 is from Sigma0
    // 2 is is from Xim
    // 3 is is the flat feeddown (it should be just xi0, i.e. == xim)
    //   remark: up to Aug 2020, the material contribution, c.a. 1% was distributed amount ALL other
    //   entries, primary, feed etc. Actually it would make most sence to put it as an additional flat contribution!
    // 4 is for the missid
    Fractions[0] = SigLambdaPrimDir * FracOfLambda;
    Fractions[1] = SigLambdaPrimDir * (1. - FracOfLambda);
    Fractions[2] = (1. - SigLambdaPrimDir) * (arrayPercXiLambda);
    Fractions[3] = 1. - Fractions[0] - Fractions[1] - Fractions[2];
    Fractions[4] = 1.;
}

void DLM_CommonAnaFunctions::GetFractions_Xim(const TString &DataSample, const int &Variation, double *Fractions)
{
    // 0 = primary
    // 1 = from Xi-(1530)
    // 2 = from Omega
    // 3 = flat
    // 4 = missidentified

    // ratio Xi-(1530) to Xi-
    const double Xim1530_to_Xim = 0.32 * (1. / 3.);
    // ratio Xi0(1530) to Xi0 (n=neutral)
    // const double Xin1530_to_Xim = 0.32*(2./3.);
    const double Omegam_to_Xim = 0.1;
    const double OmegamXim_BR = 0.086;
    // the ratios that we have for Xis are referred to the total number of Xi particles (which already include all contributions)
    // hence Xi1530_to_Xi indeed is simply the number of Xis that stem from a Xi1530
    Fractions[0] = 1. - 3. * Xim1530_to_Xim - Omegam_to_Xim * OmegamXim_BR;
    Fractions[1] = Xim1530_to_Xim;
    Fractions[2] = Omegam_to_Xim * OmegamXim_BR;
    Fractions[3] = 1. - Fractions[2] - Fractions[1] - Fractions[0];
    Fractions[4] = 1.;
}

void DLM_CommonAnaFunctions::GetFractions_K(const TString &DataSample, const int &Variation, double *Fractions)
{
    // following my lambda pars with the 3 possible modifications
    // for the kaon:
    // 0 = primary
    // 1 = from Phi
    // 2 = other feeddown (flat)
    // 3 = missidentified
    double Modify_pp = 1;
    switch (Variation % 10)
    {
    case 1:
        Modify_pp = 0.9;
        break;
    case 2:
        Modify_pp = 1.1;
        break;
    default:
        Modify_pp = 1;
        break;
    }
    double pp_f0; // primary Kaons
    double pp_f1; // fraction from Phi
    if (DataSample == "pp13TeV_HM_LK")
    {
        pp_f0 = 0.94;
        pp_f1 = 0.06;
    }
    else
    {
        printf("\033[1;31mERROR (GetFractions_p):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        pp_f0 = 1.0;
        pp_f1 = 0.0;
    }
    // ratio between feed-down Phi and all feed downs, it is by default 0.06, i.e. 6% of the feed-down into Kaons is from Phi
    double arrayPercLamProton = pp_f1 / (1. - pp_f0) * Modify_pp;
    Fractions[0] = pp_f0;
    Fractions[1] = (1. - pp_f0) * (arrayPercLamProton);
    Fractions[2] = (1. - pp_f0) * (1. - arrayPercLamProton);
    Fractions[3] = 1.;
}

// 0 is primary
// 1 is pL->pp
// 2 is flat feed
// 3 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pp(const TString &DataSample, const int &Variation_p, double *lambda_pars)
{
    double Purities_p[4];
    double Fraction_p[4];
    GetPurities_p(DataSample, Variation_p, Purities_p);
    GetFractions_p(DataSample, Variation_p, Fraction_p);
    lambda_pars[0] = Purities_p[0] * Fraction_p[0] * Purities_p[0] * Fraction_p[0];
    lambda_pars[1] = Purities_p[0] * Fraction_p[0] * Purities_p[1] * Fraction_p[1] * 2.;
    lambda_pars[3] = (Purities_p[0] + Purities_p[0] + Purities_p[3]) * Purities_p[3];
    lambda_pars[2] = 1. - lambda_pars[3] - lambda_pars[1] - lambda_pars[0];

    // double SUM=0;
    // for(unsigned uLam=0; uLam<4; uLam++){
    //     printf("λ(pp)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //     SUM+=lambda_pars[uLam]*100.;
    // }
    // printf("SUM: %.1f\n------------\n",SUM);
}

// 0 is primary
// 1 is pSigma0->pL
// 2 is pXim->pL
// 3 is the flat feeddown
// 4 is missid
// Variation_L 4 digits
// the last two digits are passed to the variations of the fraction (last digit is for Lambda, the second to last for Xi)
// the first two digits are passed to the purities.
void DLM_CommonAnaFunctions::SetUpLambdaPars_pL(const TString &DataSample, const int &Variation_p, const int &Variation_L, double *lambda_pars)
{
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_L[5];
    double Fraction_L[5];
    GetPurities_p(DataSample, Variation_p, Purities_p);
    GetFractions_p(DataSample, Variation_p, Fraction_p);
    GetPurities_L(DataSample, Variation_L / 100, Purities_L);
    GetFractions_L(DataSample, Variation_L % 100, Fraction_L);
    lambda_pars[0] = Purities_p[0] * Fraction_p[0] * Purities_L[0] * Fraction_L[0];
    lambda_pars[1] = Purities_p[0] * Fraction_p[0] * Purities_L[1] * Fraction_L[1];
    lambda_pars[2] = Purities_p[0] * Fraction_p[0] * Purities_L[2] * Fraction_L[2];
    lambda_pars[4] = Purities_p[0] * Purities_L[4] + Purities_p[3] * Purities_L[0] + Purities_p[3] * Purities_L[4];
    lambda_pars[3] = 1. - lambda_pars[0] - lambda_pars[1] - lambda_pars[2] - lambda_pars[4];

    // double SUM=0;
    // for(unsigned uLam=0; uLam<5; uLam++){
    //     printf("λ(pΛ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //     SUM+=lambda_pars[uLam]*100.;
    // }
    // printf("SUM: %.1f\n------------\n",SUM);
}
// 0 is primary
// 1 is from Xim1530
// 2 is from Xin1530
// 3
// 4 is missid
void DLM_CommonAnaFunctions::SetUpLambdaPars_pXim(const TString &DataSample, const int &Variation_p, const int &Variation_Xim, double *lambda_pars)
{
    double Purities_p[4];
    double Fraction_p[4];
    double Purities_Xim[5];
    double Fraction_Xim[5];
    GetPurities_p(DataSample, Variation_p, Purities_p);
    GetFractions_p(DataSample, Variation_p, Fraction_p);
    GetPurities_Xim(DataSample, Variation_Xim, Purities_Xim);
    GetFractions_Xim(DataSample, Variation_Xim, Fraction_Xim);
    lambda_pars[0] = Purities_p[0] * Fraction_p[0] * Purities_Xim[0] * Fraction_Xim[0];
    lambda_pars[1] = Purities_p[0] * Fraction_p[0] * Purities_Xim[1] * Fraction_Xim[1];
    lambda_pars[2] = Purities_p[0] * Fraction_p[0] * Purities_Xim[2] * Fraction_Xim[2];
    lambda_pars[4] = Purities_p[0] * Purities_Xim[4] + Purities_p[3] * Purities_Xim[0] + Purities_p[3] * Purities_Xim[4];
    lambda_pars[3] = 1. - lambda_pars[0] - lambda_pars[1] - lambda_pars[2] - lambda_pars[4];
    /*
        double SUM=0;
        printf("pp * fp * pxi * fxi = %.3f * %.3f * %.3f * %.3f\n",Purities_p[0],Fraction_p[0],Purities_Xim[0],Fraction_Xim[0]);
        for(unsigned uLam=0; uLam<5; uLam++){
            printf("λ(pΞ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
            SUM+=lambda_pars[uLam]*100.;
        }
        printf("SUM: %.1f\n------------\n",SUM);
        */
}

//        DataPeriod=="pp13TeV"?  :
//                                ;

// 0 is primary
// 1 is XimKmin->LKim
// 2 is from phi
// 3 is the flat feeddown (includes Sigma0KMin->LKmin Xi0Kmin->LKmin + missid)
// Variation_L 4 digits
// the last two digits are passed to the variations of the fraction (last digit is for Lambda, the second to last for Xi)
// the first two digits are passed to the purities.
void DLM_CommonAnaFunctions::SetUpLambdaPars_LKmin(const TString &DataSample, const int &Variation_L, const int &Variation_K, double *lambda_pars, int SamplePurity)
{
    double Purities_K[4]; // Primary, from phi, from Secondaries, mis-id
    double Fraction_K[4]; // Primary, from phi, from Secondaries, mis-id
    double Purities_L[5]; //
    double Fraction_L[5]; // primary, from Sigma0, from Xim, flat, mis-id
    GetPurities_K(DataSample, Variation_K, Purities_K, SamplePurity);
    GetFractions_K(DataSample, Variation_K, Fraction_K);
    GetPurities_L_Vale(DataSample, Variation_L, Purities_L, SamplePurity);
    GetFractions_L_Vale(DataSample, Variation_L, Fraction_L);
    lambda_pars[0] = Purities_K[0] * Fraction_K[0] * Purities_L[0] * Fraction_L[0]; // Primary
    lambda_pars[1] = Purities_K[0] * Fraction_K[0] * Purities_L[2] * Fraction_L[2]; // ΞK
    lambda_pars[2] = Purities_K[1] * Fraction_K[1] * Purities_L[0] * Fraction_L[0]; // ΛΦ
    lambda_pars[3] = 1. - lambda_pars[0] - lambda_pars[1] - lambda_pars[2];         // flat

    // double SUM=0;
    // for(unsigned uLam=0; uLam<5; uLam++){
    //     printf("λ(pΛ)_%u = %.1f\n",uLam,lambda_pars[uLam]*100.);
    //     SUM+=lambda_pars[uLam]*100.;
    // }
    // printf("SUM: %.1f\n------------\n",SUM);
}

TH2F *DLM_CommonAnaFunctions::GetResolutionMatrix(const TString &DataSample, const TString &&System)
{
    TString FileName;
    TString HistoName;

    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/ALICE_pp_13TeV.root";
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19" || DataSample == "pp13TeV_HM_RotPhiDec19" || DataSample == "pp13TeV_HM_DimiJun20" || DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        // N.B. special rule for pp and pLambda below
        FileName = CatsFilesFolder[0] + "/MomentumSmear/ALICE_pp_13TeV.root";
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/Sample3_MeV_compact.root";
        ;
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        FileName = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/ResolutionMatrices/Sample3_MeV_compact.root";
    }
    else if (DataSample == "pp13TeV_HM_LK")
    {
        FileName = "/home/valentina/thor/cernbox/CATSFiles/SystematicsAndCalib/MomentumSmear/ALICE_pp_13TeVHM_MELKmin.root";
    }
    else if (DataSample == "pp13TeV_HM_LKMAC")
    {
        FileName = "/Users/sartozza/cernbox/CATSFiles/SystematicsAndCalib/MomentumSmear/ALICE_pp_13TeVHM_MELKmin.root";
    }
    else
    {
        printf("\033[1;31mERROR (GetResolutionMatrix):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        FileName = "";
    }

    if (System == "pp")
    {
        HistoName = "hSigmaMeV_Proton_Proton";
    }
    else if (System == "pLambda")
    {
        // N.B. special rule for pp and pLambda below
        HistoName = "hSigmaMeV_Proton_Lambda";
    }
    else if (System == "LambdaLambda")
    {
        HistoName = "hSigmaMeV_Lambda_Lambda";
    }
    else if (System == "pXim")
    {
        HistoName = "hSigmaMeV_Proton_Xim";
    }
    else if (System == "LKmin")
    {
        HistoName = "MomentumResolutionME_Particle1_Particle2";
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
    }

    if (DataSample == "pp13TeV_HM_DimiJun20" && System == "pp")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/ALICE_pp_13TeV_MEpp.root";
        HistoName = "h_RESO_pp_MeV";
    }
    if (DataSample == "pp13TeV_HM_DimiJun20" && System == "pLambda")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/ALICE_pp_13TeV_MEpL.root";
        HistoName = "h_RESO_pL_MeV";
    }
    if (DataSample == "pp13TeV_HM_BernieSource" && System == "pp")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/Sample6_MeV_compact.root";
        HistoName = "hSigmaMeV_Proton_Proton";
    }
    if (DataSample == "pp13TeV_HM_BernieSource" && System == "pLambda")
    {
        FileName = CatsFilesFolder[0] + "/MomentumSmear/Sample6_MeV_compact.root";
        HistoName = "hSigmaMeV_Proton_Lambda";
    }

    // this is the unfolded data
    if (DataSample == "pp13TeV_HM_DimiJul20" || DataSample == "pp13TeV_HM_DimiMay21" || DataSample == "pp13TeV_HM_DimiJun21")
    {
        return NULL;
    }

    /// FUCKING ROOT SUCKS!!!!!! SUCK MY COCK!!!! SUCK IT YOU BITCH!!!!!!!!!!
    // so we need to copy our histogram, as else we lose it when we delete the file
    // and we need to change to the "central" root directory, as else histoCopy will also be lost
    // and we need to play with the name a little bit, else we are fucked!
    TFile *FileROOT = new TFile(FileName, "read");
    TH2F *histo = (TH2F *)FileROOT->Get(HistoName);
    if (!histo)
    {
        printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n", HistoName.Data(), FileName.Data());
        return NULL;
    }
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F *)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}
TH2F *DLM_CommonAnaFunctions::GetResidualMatrix(const TString &&FinalSystem, const TString &InitialSystem)
{
    TString FileName;
    TString HistoName;

    FileName = CatsFilesFolder[0] + "/DecaySmear/run2_decay_matrices_old.root";

    if (FinalSystem == "pp" && InitialSystem == "pLambda")
    {
        HistoName = "hRes_pp_pL";
    }
    else if (FinalSystem == "pLambda" && InitialSystem == "pSigma0")
    {
        HistoName = "hRes_pL_pSigma0";
    }
    else if (FinalSystem == "pLambda" && InitialSystem == "pXim")
    {
        HistoName = "hRes_pL_pXim";
    }
    else if (FinalSystem == "pLambda" && InitialSystem == "pXi0")
    {
        FileName = CatsFilesFolder[0] + "/DecaySmear/pXi0_pL.root";
        HistoName = "pXi0_pL";
    }
    else if (FinalSystem == "pXim" && InitialSystem == "pXim1530")
    {
        HistoName = "hRes_pXim_pXim1530";
    }
    else if (FinalSystem == "LKmin" && InitialSystem == "XiminKmin")
    {
        FileName = "/home/valentina/thor/cernbox/CATSFiles/DecaySmear/histDecayKinematics_LK_pp13HM.root";
        HistoName = "KXi_KLambda";
    }
    else if (FinalSystem == "pp" && InitialSystem == "pSigmaPlus")
    {
        FileName = CatsFilesFolder[0]+"/DecaySmear/Decay_matrix_pp_pSp.root";
        HistoName = "hRes_pp_pSp";
    }
    else if (FinalSystem == "LKminMAC" && InitialSystem == "XiminKminMAC")
    {
        FileName = "/Users/sartozza/cernbox/CATSFiles/DecaySmear/histDecayKinematics_LK_pp13HM.root";
        HistoName = "KXi_KLambda";
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m The decay '%s->%s' does not exist\n", InitialSystem.Data(), FinalSystem.Data());
    }
    TFile *FileROOT = new TFile(FileName, "read");
    TH2F *histo = (TH2F *)FileROOT->Get(HistoName);
    if (!histo)
    {
        printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n", HistoName.Data(), FileName.Data());
        return NULL;
    }
    TString Name = histo->GetName();
    gROOT->cd();
    TH2F *histoCopy = (TH2F *)histo->Clone("histoCopy");
    delete FileROOT;
    histoCopy->SetName(Name);
    return histoCopy;
}

// mTbin == -1 means we take the integrated function
// iReb = 0 is 4 MeV, 1 is 8, 2 is 12, 3 is 16, 4 is 20
TH1F *DLM_CommonAnaFunctions::GetAliceExpCorrFun(const TString &DataSample, const TString &System, const TString &CutVar, const int &iReb, const bool &AddSyst, const int mTbin)
{
    TString FileName;
    TString HistoName;
    TString SystFileName = "";
    TString SystHistName = "";
    TGraph *gRelSyst = NULL;

    if (DataSample == "pp13TeV_MB_Run2paper")
    {
        if (System == "pp")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pp%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "LambdaLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_LL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pXim")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pp_13TeV/Sample8/CFOutput_pXi%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pp13TeV_HM_March19" || DataSample == "pp13TeV_HM_Dec19")
    {
        if (System == "pp")
        {
            if (mTbin == -1)
            {
                // buggy for all but _0
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pp%s.root", CutVar.Data());
                if (CutVar == "_0")
                    FileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample10HM_ver2/CFOutput_pp.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else if (mTbin >= 0 && mTbin < 7)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutputALL_mT_pp_HM%s.root", CutVar.Data());
                // ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if (CutVar == "_0")
                    FileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample10HM_ver2/CFOutputALL_mT_pp_HM.root";
                HistoName = TString::Format("hCkTotNormWeightMeV_mTBin_%i", mTbin);
                SystFileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pp.root";
                SystHistName = "SystErrRel";
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pL%s.root", CutVar.Data());
                // ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if (CutVar == "_0")
                    FileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample10HM/CFOutput_pL.root";
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
                SystFileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pL.root";
                SystHistName = "SystErrRel";
            }
            else if (mTbin >= 0 && mTbin < 6)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutputALL_mT_pL_HM%s.root", CutVar.Data());
                // ugly solution to make sure that the default fit combination is the correlation we have showed at the preliminaries preview
                if (CutVar == "_0")
                    FileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample10HM/CFOutputALL_mT_pL_HM.root";
                HistoName = TString::Format("hCk_RebinnedMeV_%i_mTBin_%i", 0, mTbin);
                SystFileName = CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/Systematics_pL.root";
                SystHistName = "SystErrRel";
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "LambdaLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_LL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pXim")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/Sample12HM/CFOutput_pXi%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    // the correlations are obtained with Dimi's code, normalization to the total yield
    // used for the pLambda paper proposal in June 2020, sideband corrected (folded)
    else if (DataSample == "pp13TeV_HM_DimiJun20")
    {
        if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/CkSB_pL_%s.root", CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV", float(iReb + 1) * 4.);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pp13TeV_HM_DimiJul20")
    {
        if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                // FileName = TString::Format(CatsFilesFolder[0]+"/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/Unfolded/090720/CkSB_pL_%s.root",CutVar.Data());
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/CkSB_pL_%s.root", CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV", float(iReb + 1) * 4.);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pp13TeV_HM_DimiMay21")
    {
        if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/DimiMay21/UnfoldedNorm240_340/CkSB_pL_%s.root", CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV", float(iReb + 1) * 4.);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pp13TeV_HM_DimiJun21")
    {
        if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/ALICE_pp_13TeV_HM/DimiJun21/UnfoldedNorm240_340/CkSB_pL_%s.root", CutVar.Data());
                HistoName = TString::Format("hCkS_Norm_%.0fMeV", float(iReb + 1) * 4.);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pp13TeV_HM_RotPhiDec19")
    {
        if (System == "pp")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pp_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pL_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "LambdaLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_LL_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pXim")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/pp13TeV_HM_Baseline/MyResults_Vale/Data/CF/ME_Methods_14122019/CFOutput_pXi_3.root");
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    // the data used by Bernie for the source paper. Only pp and pLambda
    else if (DataSample == "pp13TeV_HM_BernieSource")
    {
        if (System == "pp")
        {
            if (mTbin >= 0 && mTbin < 7)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/Bernie_Source/ppData/mTBin_%i/CFOutput_mT_ppVar%s_HM_%i.root", mTbin + 1, CutVar.Data(), mTbin);
                HistoName = TString::Format("hCk_RebinnedppVar%sMeV_0", CutVar.Data());
                if(AddSyst){
                  SystFileName = (CatsFilesFolder[0] + "/ExpData/Bernie_Source/ppData/DimiSystematics_pp.root");
                  SystHistName = TString::Format("hSyst_mT%i", mTbin + 1);
                }
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin >= 0 && mTbin < 6)
            {
                FileName = TString::Format(CatsFilesFolder[0] + "/ExpData/Bernie_Source/pLData/mTBin_%i/CFOutput_mT_pLVar%s_HM_%i.root", mTbin + 1, CutVar.Data(), mTbin);
                HistoName = TString::Format("hCk_RebinnedpLVar%sMeV_0", CutVar.Data());
                if(AddSyst){
                  SystFileName = (CatsFilesFolder[0] + "/ExpData/Bernie_Source/pLData/DimiSystematics_pL.root");
                  SystHistName = TString::Format("hSyst_mT%i", mTbin + 1);
                }
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pPb5TeV_Run2paper")
    {
        if (System == "pp")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pp%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "LambdaLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_LL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pXim")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample10/CFOutput_pXi%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else if (DataSample == "pPb5TeV_CPR_Mar19")
    {
        if (System == "pp")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pp%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "LambdaLambda")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_LL%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else if (System == "pXim")
        {
            if (mTbin == -1)
            {
                FileName = TString::Format("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/ALICE_pPb_5TeV/Sample11/CFOutput_pXi%s.root", CutVar.Data());
                HistoName = TString::Format("hCk_ReweightedMeV_%i", iReb);
            }
            else
            {
                printf("\033[1;31mERROR:\033[0m The mT bin #%i is not defined for %s (%s)\n", mTbin, DataSample.Data(), System.Data());
            }
        }
        else
        {
            printf("\033[1;31mERROR:\033[0m The system '%s' does not exist\n", System.Data());
        }
    }
    else
    {
        printf("\033[1;31mERROR (GetAliceExpCorrFun):\033[0m The data sample '%s' does not exist\n", DataSample.Data());
        FileName = "";
    }

    TFile *FileROOT = new TFile(FileName, "read");
    if (!FileROOT)
    {
        printf("\033[1;31mERROR:\033[0m The file '%s' does not exist\n", FileName.Data());
        return NULL;
    }
    TH1F *histo = (TH1F *)FileROOT->Get(HistoName);
    if (!histo)
    {
        printf("\033[1;31mERROR:\033[0m The histo '%s' in file '%s' does not exist\n", HistoName.Data(), FileName.Data());
        return NULL;
    }
    TString Name = histo->GetName();
    gROOT->cd();
    TH1F *histoCopy = (TH1F *)histo->Clone("histoCopy");
    delete FileROOT;
    FileROOT = NULL;
    histoCopy->SetName(Name);

    if (AddSyst)
    {
        FileROOT = new TFile(SystFileName, "read");
        if (FileROOT)
        {
            TH1F *hRelSyst = (TH1F *)FileROOT->Get(SystHistName);
            if (!hRelSyst)
            {
                printf("\033[1;31mERROR:\033[0m The hRelSyst '%s' does not exist\n", SystHistName.Data());
            }
            //in this case the systematics are absolute!!!
            else if(SystFileName.Contains("DimiSystematics_")){
              double MaxMom = hRelSyst->GetXaxis()->GetBinUpEdge(hRelSyst->GetNbinsX());
              for (int iBin = 1; iBin <= histoCopy->GetNbinsX(); iBin++){
                if (histoCopy->GetBinCenter(iBin) > MaxMom)
                    break;
                double StatErr = histoCopy->GetBinError(iBin);
                double SystErr = hRelSyst->GetBinError(hRelSyst->FindBin(histoCopy->GetBinCenter(iBin)));
                double TotErr = sqrt(StatErr * StatErr + SystErr * SystErr);
                histoCopy->SetBinError(iBin, TotErr);
//printf("System=%s; iBin=%i; Ck=%.4f; stat=%.4f; syst=%.4f; tot=%.4f\n",System.Data(),iBin,histoCopy->GetBinContent(iBin),StatErr,SystErr,TotErr);
              }
            }
            else
            {
                double MaxMom = hRelSyst->GetXaxis()->GetBinUpEdge(hRelSyst->GetNbinsX());
                gRelSyst = new TGraph();
                gRelSyst->SetName("gRelSyst");
                for (int iBin = 1; iBin <= hRelSyst->GetNbinsX(); iBin++)
                {
                    gRelSyst->SetPoint(iBin - 1, hRelSyst->GetBinCenter(iBin), hRelSyst->GetBinContent(iBin));
                }
                for (int iBin = 1; iBin <= histoCopy->GetNbinsX(); iBin++)
                {
                    if (histoCopy->GetBinCenter(iBin) > MaxMom)
                        break;
                    double StatErr = histoCopy->GetBinError(iBin);
                    double SystErr = histoCopy->GetBinContent(iBin) * gRelSyst->Eval(histoCopy->GetBinCenter(iBin));
                    double TotErr = sqrt(StatErr * StatErr + SystErr * SystErr);
                    histoCopy->SetBinError(iBin, TotErr);
                }
            }
        }
        if (FileROOT)
        {
            delete FileROOT;
            FileROOT = NULL;
        }
    }

    return histoCopy;
}





// POT:
//   double Gaussian based on the scat pars described in https://arxiv.org/abs/2308.16120
// SRC:
//    The RSM source is based on Dmeson -- K epos files
//    The "CECA" based source with SourceVar==1 is the on straight out of the box, using the params from the CECA paper (full fit)
void DLM_CommonAnaFunctions::SetUpCats_Kd(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    CATSparameters *cPars = NULL;

    CATSparameters *cPotPars_AVG = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    //unsigned NumChannels = 0;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "NULL" || SOURCE == "")
    {
    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "DoubleGauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 3, true);
        cPars->SetParameter(0, 1.0);
        cPars->SetParameter(1, 2.0);
        cPars->SetParameter(2, 0.5);
        Kitty.SetAnaSource(DoubleGaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.6) * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.6 * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[7].InitStability(20, 1, 2);
        CleverLevy[7].InitScale(35, 0.25, 2.0);
        CleverLevy[7].InitRad(256, 0, 64);
        CleverLevy[7].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[7].InitStability(20, 1, 2);
        CleverLevy[7].InitScale(35, 0.25, 2.0);
        CleverLevy[7].InitRad(256, 0, 64);
        CleverLevy[7].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[7], 2);
        Kitty.SetAnaSource(0, sqrt(1.6) * 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[7].InitStability(20, 1, 2);
        CleverLevy[7].InitScale(35, 0.25, 2.0);
        CleverLevy[7].InitRad(256, 0, 64);
        CleverLevy[7].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[7], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.6 * 1.2);
        Kitty.SetAnaSource(1, 1.6);
    }
    //last digit of SourceVar is the smooth sampling
    //(SourceVar/10)*10 is the cutoff value
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[7].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[7].InitStability(21, 1, 2);
        CleverMcLevyResoTM[7].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[7].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[7].InitType(2);
        CleverMcLevyResoTM[7].SetUpReso(0, 0.0);//deuterons
        CleverMcLevyResoTM[7].SetUpReso(1, 0.476);//kaons
        // pure Gauss

        // EPOS, 2 is with fixed mass, 3 is with EPOS mass, 4 is 3 body with fixed mass, 5 is 3 body with EPOS mass
        // printf("Hello 2\n");
        const double k_CutOff = fabs(int(int(SourceVar) / 10) * 10.);
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

        double MeanP2 = 0;
        double RanVal2 = 0;


        TFile *F_EposDisto_Kreso_d;
        F_EposDisto_Kreso_d = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/TEMP/SecondTry/ALL_D_KaonReso.root");
        TNtuple *T_EposDisto_Kreso_d = (TNtuple *)F_EposDisto_Kreso_d->Get("InfoTuple_ClosePairs");
        unsigned N_EposDisto_Kreso_d = T_EposDisto_Kreso_d->GetEntries();
        T_EposDisto_Kreso_d->SetBranchAddress("k_D", &k_D);
        T_EposDisto_Kreso_d->SetBranchAddress("P1", &fP1);
        T_EposDisto_Kreso_d->SetBranchAddress("P2", &fP2);
        T_EposDisto_Kreso_d->SetBranchAddress("M1", &fM1);
        T_EposDisto_Kreso_d->SetBranchAddress("M2", &fM2);
        T_EposDisto_Kreso_d->SetBranchAddress("Tau1", &Tau1);
        T_EposDisto_Kreso_d->SetBranchAddress("Tau2", &Tau2);
        T_EposDisto_Kreso_d->SetBranchAddress("AngleRcP1", &AngleRcP1);
        T_EposDisto_Kreso_d->SetBranchAddress("AngleRcP2", &AngleRcP2);
        T_EposDisto_Kreso_d->SetBranchAddress("AngleP1P2", &AngleP1P2);


        gROOT->cd();
        TH1F* hAngle = new TH1F("hAngle","hAngle",32,0,TMath::Pi());
        TH1F* hCos = new TH1F("hCos","hCos",32,-1.,1.);
        TH1F* hFinalAngle = new TH1F("hFinalAngle","hFinalAngle",32,0,TMath::Pi());
        TH1F* hFinalCos = new TH1F("hFinalCos","hFinalCos",32,-1.,1.);
        F_EposDisto_Kreso_d->cd();
        int NumUsefulEntries = 0;
        for(unsigned uEntry=0; uEntry<N_EposDisto_Kreso_d; uEntry++){
          T_EposDisto_Kreso_d->GetEntry(uEntry);
          if(k_D>k_CutOff) continue;
          hAngle->Fill(AngleRcP2);
          hCos->Fill(cos(AngleRcP2));
          MeanP2 += fP2;
          NumUsefulEntries++;
        }
        MeanP2 /= double(NumUsefulEntries);

        hCos->Scale(1./hAngle->Integral(),"width");


      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_Kreso_d; uEntry++){
          //get each entry
          T_EposDisto_Kreso_d->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          //overwrite the value for the lifetime. This is computed from the
          //stat. hadronization model (Vale) or thermal fist (Max)
          //this is the value for the secondary protons
          Tau1 = 0;
          //for primoridials (the Xis) we put 0
          Tau2 = 3.66;
          //put in the average mass of the resonances (again from SHM or TF)
          //this is the value for protons
          fM2 = 1054;
          //generate a random path length for the propagation of the resonances
          //nothing to change!
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          //adds a single entry into the PDF for the angular distribution to be used
          CleverMcLevyResoTM[7].AddBGT_PR(RanVal2,cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }

      delete hAngle;
      delete hCos;
      delete hFinalAngle;
      delete hFinalCos;

        delete F_EposDisto_Kreso_d;


        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[7].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[7].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[7], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_Kd;
    }

//f_best     -0.4691022467008116 //V1
//d_best     1.7127749936682812 //V1
//f_best     -0.4703530995890609//V2
//d_best     1.7532889758457184    //V2
//IMPORTANT!!! The sign of d0 was wrong, also in the PRX, now flipped in the V3:
//f_best     -0.4701902500023152 //V3
//d_best     -1.7527813610447043 //V3
    if (POT == "DG_ER")
    {
        //double PotPars_AVG[4] = {676.1923588303322, 0.7595548889641572, -44.11395844823363, 1.3254079885988337};//V1
        //double PotPars_AVG[4] = {677.6412964551082, 0.7649742727111178, -46.04483165662934, 1.3250780448057782};//V2
        double PotPars_AVG[4] = {121.19992348368771, 0.2945870040375218, 29.909893773423107, 1.3329202449286048};//V3
        
        cPotPars_AVG = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars_AVG->SetParameters(PotPars_AVG);
    }
    //f_best     -0.5399918971499024//V1
    //d_best     0.02053659406423805//V1
    //f_best     -0.5400300923056585//V2
    //d_best     0.003990374220081213//V2
    else if (POT == "DG_FCA")
    {
        //double PotPars_AVG[4] = {369.4094817373952, 0.7847035161927353, -125.73561350111179, 0.7507974357114263};//V1
        double PotPars_AVG[4] = {362.7384262781105, 0.7962158970471339, -130.77196423166328, 0.7685517499557406};//V2
        cPotPars_AVG = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars_AVG->SetParameters(PotPars_AVG);    
    }
    //square well potential,see mail from Johann @ CERN on 27.11.2024
    else if(POT == "SW_ER"){
        double PotPars_AVG[2] = {9.70, 2.14};
        cPotPars_AVG = new CATSparameters(CATSparameters::tPotential, 8, true);
        cPotPars_AVG->SetParameters(PotPars_AVG);   
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing kd potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_Kd;
    }
    Kitty.SetMomentumDependentSource(false);
    // Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(1);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass((Mass_Kch*Mass_d)/(Mass_Kch+Mass_d));

    if(!Kitty.GetNumChannels()){
        Kitty.SetNumChannels(1);
        Kitty.SetNumPW(0, 1);
        Kitty.SetSpin(0, 0);
        Kitty.SetChannelWeight(0, 1.);
    }

    if(POT=="SW_ER"){
        Kitty.SetShortRangePotential(0,0,SquareWell,*cPotPars_AVG);
    }
    else if(POT != "")
        Kitty.SetShortRangePotential(0,0,DoubleGaussSum,*cPotPars_AVG);

CLEAN_SetUpCats_Kd:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotPars_AVG)
    {
        delete cPotPars_AVG;
        cPotPars_AVG = NULL;
    }
}





// SRC:
//    The RSM source is based on Dmeson -- pi epos files
//SourceVar, 4, ABCDE
//CDE = k_cutoff, with E being actually variation on the amount of resonances, 0 = default, 1 = +10%, 2 = -10%
//A is Tau variation, 0 = default, 1 = +10%, 2 = -10%
//B is ResoMass variation, 0 = default, 1 = +10%, 2 = -10%
//if POT=="Gauss":
//  a positive PotVar means pi^+ d
//  a negative means pi^- d
//  -1 is the -0.037 fm of the real part as taken from Nuclear Physics A 872 (2011) 69–116 
void DLM_CommonAnaFunctions::SetUpCats_pi_d(CATS &Kitty, const TString &POT, const TString &SOURCE, const int &PotVar, const int &SourceVar)
{

    CATSparameters *cPars = NULL;

    CATSparameters *cPotPars = NULL;

    DLM_Histo<complex<double>> ***ExternalWF = NULL;
    //unsigned NumChannels = 0;

    Kitty.SetThetaDependentSource(false);

    if (SOURCE == "NULL" || SOURCE == "")
    {
    }
    else if (SOURCE == "Gauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "DoubleGauss")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 3, true);
        cPars->SetParameter(0, 1.0);
        cPars->SetParameter(1, 2.0);
        cPars->SetParameter(2, 0.5);
        Kitty.SetAnaSource(DoubleGaussSource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "GaussTheta")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(GaussSourceTheta, *cPars);
        Kitty.SetThetaDependentSource(true);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Cauchy")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 1, true);
        cPars->SetParameter(0, 1.2);
        Kitty.SetAnaSource(CauchySource, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Nolan")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Single")
    {
        cPars = new CATSparameters(CATSparameters::tSource, 2, true);
        cPars->SetParameter(0, sqrt(1.6) * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_single, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "Levy_Diff")
    {
        cPars->SetParameter(0, 0.5 * 1.6 * 1.2);
        cPars->SetParameter(1, 1.6);
        Kitty.SetAnaSource(LevySource3D_2particle, *cPars);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Nolan")
    {
        CleverLevy[8].InitStability(20, 1, 2);
        CleverLevy[8].InitScale(35, 0.25, 2.0);
        CleverLevy[8].InitRad(256, 0, 64);
        CleverLevy[8].InitType(2);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[0], 2);
        Kitty.SetAnaSource(0, 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Single")
    {
        CleverLevy[8].InitStability(20, 1, 2);
        CleverLevy[8].InitScale(35, 0.25, 2.0);
        CleverLevy[8].InitRad(256, 0, 64);
        CleverLevy[8].InitType(0);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[8], 2);
        Kitty.SetAnaSource(0, sqrt(1.6) * 1.2);
        Kitty.SetAnaSource(1, 1.6);
        Kitty.SetUseAnalyticSource(true);
    }
    else if (SOURCE == "CleverLevy_Diff")
    {
        CleverLevy[8].InitStability(20, 1, 2);
        CleverLevy[8].InitScale(35, 0.25, 2.0);
        CleverLevy[8].InitRad(256, 0, 64);
        CleverLevy[8].InitType(1);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverLevy[8], 2);
        Kitty.SetUseAnalyticSource(true);
        Kitty.SetAnaSource(0, 0.5 * 1.6 * 1.2);
        Kitty.SetAnaSource(1, 1.6);
    }
    //last digit of SourceVar is the smooth sampling
    //(SourceVar/10)*10 is the cutoff value
    else if (SOURCE == "McGauss_ResoTM" || SOURCE == "McLevy_ResoTM")
    {
        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[8].InitStability(1, 2 - 1e-6, 2 + 1e-6);
        else
            CleverMcLevyResoTM[8].InitStability(21, 1, 2);
        CleverMcLevyResoTM[8].InitScale(38, 0.15, 2.0);
        CleverMcLevyResoTM[8].InitRad(257 * 2, 0, 64);
        CleverMcLevyResoTM[8].InitType(2);
        CleverMcLevyResoTM[8].SetUpReso(0, 0.0);//deuterons
        if(SourceVar%10==1) CleverMcLevyResoTM[8].SetUpReso(1, 0.682*1.1);//pions
        else if(SourceVar%10==2) CleverMcLevyResoTM[8].SetUpReso(1, 0.682*0.9);//pions
        else CleverMcLevyResoTM[8].SetUpReso(1, 0.682);//pions

        // pure Gauss

        // EPOS, 2 is with fixed mass, 3 is with EPOS mass, 4 is 3 body with fixed mass, 5 is 3 body with EPOS mass
        // printf("Hello 2\n");
        const double k_CutOff = fabs( (int(int(SourceVar) / 10))%100 * 10.);
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

        double MeanP2 = 0;
        double RanVal2 = 0;


        TFile *F_EposDisto_d_piReso;
        F_EposDisto_d_piReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/TEMP/SecondTry/ALL_D_piReso.root");
        TNtuple *T_EposDisto_d_piReso = (TNtuple *)F_EposDisto_d_piReso->Get("InfoTuple_ClosePairs");
        unsigned N_EposDisto_d_piReso = T_EposDisto_d_piReso->GetEntries();
        T_EposDisto_d_piReso->SetBranchAddress("k_D", &k_D);
        T_EposDisto_d_piReso->SetBranchAddress("P1", &fP1);
        T_EposDisto_d_piReso->SetBranchAddress("P2", &fP2);
        T_EposDisto_d_piReso->SetBranchAddress("M1", &fM1);
        T_EposDisto_d_piReso->SetBranchAddress("M2", &fM2);
        T_EposDisto_d_piReso->SetBranchAddress("Tau1", &Tau1);
        T_EposDisto_d_piReso->SetBranchAddress("Tau2", &Tau2);
        T_EposDisto_d_piReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
        T_EposDisto_d_piReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
        T_EposDisto_d_piReso->SetBranchAddress("AngleP1P2", &AngleP1P2);


        gROOT->cd();
        TH1F* hAngle = new TH1F("hAngle","hAngle",32,0,TMath::Pi());
        TH1F* hCos = new TH1F("hCos","hCos",32,-1.,1.);
        TH1F* hFinalAngle = new TH1F("hFinalAngle","hFinalAngle",32,0,TMath::Pi());
        TH1F* hFinalCos = new TH1F("hFinalCos","hFinalCos",32,-1.,1.);
        F_EposDisto_d_piReso->cd();
        int NumUsefulEntries = 0;
        for(unsigned uEntry=0; uEntry<N_EposDisto_d_piReso; uEntry++){
          T_EposDisto_d_piReso->GetEntry(uEntry);
          if(k_D>k_CutOff) continue;
          hAngle->Fill(AngleRcP2);
          hCos->Fill(cos(AngleRcP2));
          MeanP2 += fP2;
          NumUsefulEntries++;
        }
        MeanP2 /= double(NumUsefulEntries);

        hCos->Scale(1./hAngle->Integral(),"width");


      //iterate over the ntuple
      for(unsigned uEntry=0; uEntry<N_EposDisto_d_piReso; uEntry++){
          //get each entry
          T_EposDisto_d_piReso->GetEntry(uEntry);
          //disregard the entry of you are outside the desired k*
          if(k_D>k_CutOff) continue;
          //overwrite the value for the lifetime. This is computed from the
          //stat. hadronization model (Vale) or thermal fist (Max)
          //this is the value for the secondary protons
          Tau1 = 0;
          //for primoridials (the Xis) we put 0
          Tau2 = 1.50;
          if(SourceVar/10000==1) Tau2*=1.1;
          else if(SourceVar/10000==2) Tau2*=0.9;
          //put in the average mass of the resonances (again from SHM or TF)
          //this is the value for pions
          fM2 = 1124;
          if((SourceVar/1000)%10==1) fM2*=1.1;
          else if((SourceVar/1000)%10==2) fM2*=0.9;
          //generate a random path length for the propagation of the resonances
          //nothing to change!
          RanVal2 = RanGen.Exponential(fM2/(fP2*Tau2));
          //adds a single entry into the PDF for the angular distribution to be used
          CleverMcLevyResoTM[8].AddBGT_PR(RanVal2,cos(AngleRcP2));
          hFinalAngle->Fill(AngleRcP2);
          hFinalCos->Fill(cos(AngleRcP2));
      }

      delete hAngle;
      delete hCos;
      delete hFinalAngle;
      delete hFinalCos;

        delete F_EposDisto_d_piReso;


        if (SOURCE == "McGauss_ResoTM")
            CleverMcLevyResoTM[8].InitNumMcIter(1000000);
        else
            CleverMcLevyResoTM[8].InitNumMcIter(100000);
        Kitty.SetAnaSource(CatsSourceForwarder, &CleverMcLevyResoTM[8], 2);
        Kitty.SetAnaSource(0, 1.0);
        Kitty.SetAnaSource(1, 2.0);
        Kitty.SetUseAnalyticSource(true);
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing source '%s'\n", SOURCE.Data());
        goto CLEAN_SetUpCats_pi_d;
    }

    if(POT==""){

    }
    else if(POT == "Gauss"){
        double PotPars[2] = {273.0806, 0.3778127};
        cPotPars = new CATSparameters(CATSparameters::tPotential, 2, true);
        cPotPars->SetParameters(PotPars);    
    }
    //else if (POT == "DG_ER")
    //{
    //    double PotPars_AVG[4] = {0, 1, 0, 1};
    //    cPotPars_AVG = new CATSparameters(CATSparameters::tPotential, 8, true);
    //    cPotPars_AVG->SetParameters(PotPars_AVG);
    //}
    //else if (POT == "DG_FCA")
    //{
    //    double PotPars_AVG[4] = {0, 1, 0, 1};
    //    cPotPars_AVG = new CATSparameters(CATSparameters::tPotential, 8, true);
    //    cPotPars_AVG->SetParameters(PotPars_AVG);    
    //}
    else
    {
        printf("\033[1;31mERROR:\033[0m Non-existing pi_d potential '%s'\n", POT.Data());
        goto CLEAN_SetUpCats_pi_d;
    }
    Kitty.SetMomentumDependentSource(false);
    // Kitty.SetThetaDependentSource(false);
    Kitty.SetExcludeFailedBins(false);

    Kitty.SetQ1Q2(-1);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass((Mass_pic*Mass_d)/(Mass_pic+Mass_d));

    if(!Kitty.GetNumChannels()){
        Kitty.SetNumChannels(1);
        Kitty.SetNumPW(0, 1);
        Kitty.SetSpin(0, 0);
        Kitty.SetChannelWeight(0, 1.);
    }


    if(cPotPars && POT=="Gauss")
        Kitty.SetShortRangePotential(0,0,SingleGauss,*cPotPars);

CLEAN_SetUpCats_pi_d:;
    if (cPars)
    {
        delete cPars;
        cPars = NULL;
    }
    // if(CleverLevy){delete CleverLevy; CleverLevy=NULL;}
    if (cPotPars)
    {
        delete cPotPars;
        cPotPars = NULL;
    }
}












DLM_CleverMcLevyReso *DLM_CommonAnaFunctions::GetCleverMcLevyReso_pp()
{
    return &CleverMcLevyReso[0];
}

DLM_CleverMcLevyReso *DLM_CommonAnaFunctions::GetCleverMcLevyReso_pL()
{
    return &CleverMcLevyReso[1];
}

DLM_CleverMcLevyReso *DLM_CommonAnaFunctions::GetCleverMcLevyReso_pXim()
{
    return &CleverMcLevyReso[2];
}

DLM_CleverMcLevyReso *DLM_CommonAnaFunctions::GetCleverMcLevyReso_pOmegam()
{
    return &CleverMcLevyReso[3];
}

DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pp()
{
    return &CleverMcLevyResoTM[0];
}

DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pL()
{
    return &CleverMcLevyResoTM[1];
}

DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pXim()
{
    return &CleverMcLevyResoTM[2];
}
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pOmegam()
{
    return &CleverMcLevyResoTM[3];
}
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pipi()
{
    return &CleverMcLevyResoTM[4];
}
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_ppic()
{
    return &CleverMcLevyResoTM[6];
}
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_Kd()
{
    return &CleverMcLevyResoTM[7];
}
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GetCleverMcLevyResoTM_pi_d()
{
    return &CleverMcLevyResoTM[8];
}

DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_LK(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.6438);
    MagicSource->SetUpReso(1, 0.476);

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_LK = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/LK_ALL.root");
    TNtuple *T_EposDisto_LK = (TNtuple *)F_EposDisto_LK->Get("nt_LK");
    T_EposDisto_LK->SetBranchAddress("Type", &Type);
    T_EposDisto_LK->SetBranchAddress("k_D", &k_D);
    T_EposDisto_LK->SetBranchAddress("P1", &fP1);
    T_EposDisto_LK->SetBranchAddress("P2", &fP2);
    T_EposDisto_LK->SetBranchAddress("M1", &fM1);
    T_EposDisto_LK->SetBranchAddress("M2", &fM2);
    T_EposDisto_LK->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_LK->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_LK->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_LK->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_LK->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < T_EposDisto_LK->GetEntries(); uEntry++)
    {
        T_EposDisto_LK->GetEntry(uEntry);
        if (Type == PPid)
        {
            continue;
        }
        else if (Type == PRid)
        {
            Tau1 = 0;
            Tau2 = 3.66;
            fM2 = 1054;
            if (k_D > k_CutOff)
                continue;
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
        }
        else if (Type == PRid)
        {
            Tau1 = 4.69;
            Tau2 = 0;
            fM1 = 1463;
            if (k_D > k_CutOff)
                continue;
            RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
            MagicSource->AddBGT_RP(RanVal1, cos(AngleRcP1));
        }
        else if (Type == RRid)
        {
            Tau1 = 4.69;
            Tau2 = 3.66;
            fM1 = 1463;
            fM2 = 1054;

            if (k_D > k_CutOff)
                continue;
            RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
        }
        else
        {
            continue;
        }
    }
    delete F_EposDisto_LK;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

// Crosscheck on effect of CC in LK correlations
// Assuming same kinematics as ΛΚ but now for Σ we have 0.6265 primordial
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_SigmaK(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.3735); // Res. contrib. for Σ
    MagicSource->SetUpReso(1, 0.476);  // Res. contrib. for Κ.

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_LK = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/LK_ALL.root");
    TNtuple *T_EposDisto_LK = (TNtuple *)F_EposDisto_LK->Get("nt_LK");
    T_EposDisto_LK->SetBranchAddress("Type", &Type);
    T_EposDisto_LK->SetBranchAddress("k_D", &k_D);
    T_EposDisto_LK->SetBranchAddress("P1", &fP1);
    T_EposDisto_LK->SetBranchAddress("P2", &fP2);
    T_EposDisto_LK->SetBranchAddress("M1", &fM1);
    T_EposDisto_LK->SetBranchAddress("M2", &fM2);
    T_EposDisto_LK->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_LK->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_LK->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_LK->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_LK->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < T_EposDisto_LK->GetEntries(); uEntry++)
    {
        T_EposDisto_LK->GetEntry(uEntry);
        if (Type == PPid)
        {
            continue;
        }
        else if (Type == PRid)
        {
            Tau1 = 0;
            Tau2 = 3.66;
            fM2 = 1054;
            if (k_D > k_CutOff)
                continue;
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
        }
        else if (Type == PRid)
        {
            Tau1 = 4.69;
            Tau2 = 0;
            fM1 = 1463;
            if (k_D > k_CutOff)
                continue;
            RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
            MagicSource->AddBGT_RP(RanVal1, cos(AngleRcP1));
        }
        else if (Type == RRid)
        {
            Tau1 = 4.69;
            Tau2 = 3.66;
            fM1 = 1463;
            fM2 = 1054;

            if (k_D > k_CutOff)
                continue;
            RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
            RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
            MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
        }
        else
        {
            continue;
        }
    }
    delete F_EposDisto_LK;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

// Crosscheck on effect of CC in LK correlations
// Assuming kinematics taken from Dπ analysis, which used Ωπ
// ctau and Meff given by pions since other particles is assumed primordian
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_XiPi(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.0);   // Ξ assumed primary
    MagicSource->SetUpReso(1, 0.682); // Res. contrib. for π.
    MagicSource->InitNumMcIter(1000000);

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_XiPi = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ALL_D_piReso.root");
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto_XiPi->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < T_EposDisto->GetEntries(); uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        Tau1 = 0;
        Tau2 = 1.5;
        fM2 = 1124;

        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2)); // here put angle between
    }
    delete F_EposDisto_XiPi;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

// Crosscheck on effect of CC in LK correlations
// Assuming kinematics taken from pΞ analysis
// ctau and Meff given by pions since other particles is assumed primordian
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_XiEta(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.51); // η with resonances
    MagicSource->SetUpReso(1, 0.);   // Ξ assumed primary
    MagicSource->InitNumMcIter(1000000);

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    // TFile *F_EposDisto_XiEta = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/EposDisto_pReso_Xim.root");
    TFile *F_EposDisto_XiEta = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/Vale_etaReso_Xi.root");

    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto_XiEta->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < T_EposDisto->GetEntries(); uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        Tau1 = 2.1;
        Tau2 = 0.;
        fM1 = 1230;

        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource->AddBGT_RP(RanVal1, cos(AngleRcP1)); // here put angle between
    }
    delete F_EposDisto_XiEta;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

// Getting total source for ΞK analysis
// Assuming kinematics taken from DK analysis, which used ΩK
// ctau and Meff given by kaons since other particles is assumed primordian
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_XiK(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.0);   // Ξ assumed primary
    MagicSource->SetUpReso(1, 0.476); // Res. contrib. for K.
    MagicSource->InitNumMcIter(1000000);

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_XiK = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ALL_D_KaonReso.root");
    TNtuple *T_EposDisto = (TNtuple *)F_EposDisto_XiK->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto = T_EposDisto->GetEntries();
    T_EposDisto->SetBranchAddress("k_D", &k_D);
    T_EposDisto->SetBranchAddress("P1", &fP1);
    T_EposDisto->SetBranchAddress("P2", &fP2);
    T_EposDisto->SetBranchAddress("M1", &fM1);
    T_EposDisto->SetBranchAddress("M2", &fM2);
    T_EposDisto->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto->SetBranchAddress("AngleP1P2", &AngleP1P2);

    for (unsigned uEntry = 0; uEntry < T_EposDisto->GetEntries(); uEntry++)
    {
        T_EposDisto->GetEntry(uEntry);
        if (k_D > k_CutOff)
            continue;
        Tau1 = 0;
        Tau2 = 3.66; // Kaons
        fM2 = 1054;  // Kaons

        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2)); // here put angle between
    }
    delete F_EposDisto_XiK;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

// Getting total source for Λπ analysis
// Following the work on Kp coupled channel analysis, resonances from Σπ
DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_LPi(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.6438);//lambda
    MagicSource->SetUpReso(1, 0.682);//Pion

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_pReso_Kaon = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForMaxRamona_pi_SigReso.root");
    TNtuple *T_EposDisto_pReso_Kaon = (TNtuple *)F_EposDisto_pReso_Kaon->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Kaon = T_EposDisto_pReso_Kaon->GetEntries();
    T_EposDisto_pReso_Kaon->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Kaon->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Kaon; uEntry++)
    {
        // get each entry
        T_EposDisto_pReso_Kaon->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 4.69; // Lambda
        Tau2 = 0;
        fM1 = 1463;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource->AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete T_EposDisto_pReso_Kaon;
    delete F_EposDisto_pReso_Kaon;

    TFile *F_EposDisto_p_KaonReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForMaxRamona_piReso_Sig.root");
    TNtuple *T_EposDisto_p_KaonReso = (TNtuple *)F_EposDisto_p_KaonReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_KaonReso = T_EposDisto_p_KaonReso->GetEntries();
    T_EposDisto_p_KaonReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_KaonReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_KaonReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_KaonReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_KaonReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_KaonReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_KaonReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_KaonReso; uEntry++)
    {
        // get each entry
        T_EposDisto_p_KaonReso->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 0;
        Tau2 = 1.5;
        fM2 = 1130;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
    }
    delete T_EposDisto_p_KaonReso;
    delete F_EposDisto_p_KaonReso;

    TFile *F_EposDisto_pReso_KaonReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForMaxRamona_piReso_SigReso.root");
    // TFile* F_EposDisto_pReso_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_dReso.root",GetCernBoxDimi()));

    TNtuple *T_EposDisto_pReso_KaonReso = (TNtuple *)F_EposDisto_pReso_KaonReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_KaonReso = T_EposDisto_pReso_KaonReso->GetEntries();
    T_EposDisto_pReso_KaonReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_KaonReso; uEntry++)
    {
        // get each entry
        T_EposDisto_pReso_KaonReso->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 4.69; // Lambda
        fM1 = 1463;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        Tau2 = 1.5;//pion
        fM2 = 1130;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        // check the signs
        MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
    delete T_EposDisto_pReso_KaonReso;
    delete F_EposDisto_pReso_KaonReso;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}
/*
// the cut off scale in k*, for which the angular distributions from EPOS
// are evaluated. 200 MeV works okay, you can go up to 300 MeV for systematic checks
const double k_CutOff = 200;

// to be used for the NTuple later on
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
// random generator dimi style. The input is incompatible with the ROOT random generator,
// do not mix and match, do not ask me how I know this. Ask Bernie.
// 11 is the seed, you can change that to you favorite number
DLM_Random RanGen(11);
// dummies to save random shit
double RanVal1;
double RanVal2;
double RanVal3;
*/

DLM_CleverMcLevyResoTM *DLM_CommonAnaFunctions::GaussCoreRsm_pK(const int &SourceVar)
{
    DLM_CleverMcLevyResoTM *MagicSource = new DLM_CleverMcLevyResoTM();
    MagicSource->InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource->InitScale(38, 0.15, 2.0);
    MagicSource->InitRad(257 * 2, 0, 64);
    MagicSource->InitType(2);
    MagicSource->SetUpReso(0, 0.6422);
    MagicSource->SetUpReso(1, 0.476);

    const double k_CutOff = int(int(SourceVar) / 10) * 10.;
    const int SVAR = SourceVar % 10;
    int PPid, PRid, RPid, RRid;
    // EPOS
    if (SVAR == 0)
    {
        PPid = 0;
        PRid = 1;
        RPid = 10;
        RPid = 11;
    }
    // CECA
    else if (SVAR == 1)
    {
        PPid = 100;
        PRid = 101;
        RPid = 110;
        RPid = 111;
    }
    else
    {
        printf("\033[1;31mERROR:\033[0m Unknown source variation for LK\n");
        delete MagicSource;
        return NULL;
    }

    Float_t Type;
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

    TFile *F_EposDisto_pReso_Kaon = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForRamona_pReso_Kaon.root");
    TNtuple *T_EposDisto_pReso_Kaon = (TNtuple *)F_EposDisto_pReso_Kaon->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Kaon = T_EposDisto_pReso_Kaon->GetEntries();
    T_EposDisto_pReso_Kaon->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Kaon->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Kaon->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Kaon; uEntry++)
    {
        // get each entry
        T_EposDisto_pReso_Kaon->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 1.65;
        Tau2 = 0;
        fM1 = 1362;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource->AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete T_EposDisto_pReso_Kaon;
    delete F_EposDisto_pReso_Kaon;

    TFile *F_EposDisto_p_KaonReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForRamona_p_KaonReso.root");
    TNtuple *T_EposDisto_p_KaonReso = (TNtuple *)F_EposDisto_p_KaonReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_KaonReso = T_EposDisto_p_KaonReso->GetEntries();
    T_EposDisto_p_KaonReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_KaonReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_KaonReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_KaonReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_KaonReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_KaonReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_KaonReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_KaonReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_KaonReso; uEntry++)
    {
        // get each entry
        T_EposDisto_p_KaonReso->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 0;
        Tau2 = 3.66;
        fM2 = 1054;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource->AddBGT_PR(RanVal2, cos(AngleRcP2));
    }
    delete T_EposDisto_p_KaonReso;
    delete F_EposDisto_p_KaonReso;

    TFile *F_EposDisto_pReso_KaonReso = new TFile(CatsFilesFolder[0] + "/Source/EposAngularDist/ForRamona_pReso_KaonReso.root");
    // TFile* F_EposDisto_pReso_KaonReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_dReso.root",GetCernBoxDimi()));

    TNtuple *T_EposDisto_pReso_KaonReso = (TNtuple *)F_EposDisto_pReso_KaonReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_KaonReso = T_EposDisto_pReso_KaonReso->GetEntries();
    T_EposDisto_pReso_KaonReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_KaonReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    gROOT->cd();
    // iterate over the ntuple
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_KaonReso; uEntry++)
    {
        // get each entry
        T_EposDisto_pReso_KaonReso->GetEntry(uEntry);
        // disregard the entry of you are outside the desired k*
        if (k_D > k_CutOff)
            continue;
        Tau1 = 1.65;
        fM1 = 1362;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        Tau2 = 3.66;
        fM2 = 1054;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        // check the signs
        MagicSource->AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
    delete T_EposDisto_pReso_KaonReso;
    delete F_EposDisto_pReso_KaonReso;

    MagicSource->InitNumMcIter(1000000);

    return MagicSource;
}

void DLM_CommonAnaFunctions::SetCatsFilesFolder(const TString &folder)
{
    CatsFilesFolder[0] = folder;
}

DLM_Histo<double> *ConvertThetaAngleHisto(const TString &FileName, const TString &HistoName, const double kMin, const double kMax)
{
    TFile *InputFile = new TFile(FileName, "read");
    TH2F *InputHisto = NULL;
    if (InputFile)
    {
        InputHisto = (TH2F *)InputFile->Get(HistoName);
    }
    TH1D *Projection = NULL;
    if (InputHisto)
    {
        Projection = InputHisto->ProjectionX("ConvertThetaAngleHisto", InputHisto->GetYaxis()->FindBin(kMin), InputHisto->GetYaxis()->FindBin(kMax));
    }

    const unsigned NumBins = Projection->GetNbinsX();
    Projection->Scale(1. / Projection->Integral(1, NumBins));

    DLM_Histo<double> CummDistr;
    CummDistr.SetUp(1);
    // printf("NumBins=%u; %f --> %f\n",NumBins,Projection->GetBinLowEdge(1),Projection->GetXaxis()->GetBinUpEdge(NumBins));
    // usleep(4e6);
    // printf("Calling SetUp\n");
    CummDistr.SetUp(0, NumBins, Projection->GetBinLowEdge(1), Projection->GetXaxis()->GetBinUpEdge(NumBins));
    // printf("Ended SetUp\n");
    // usleep(4e6);
    // printf("Initialize...\n");
    CummDistr.Initialize();
    // printf("Ended...\n");
    CummDistr.SetBinContent(unsigned(0), Projection->GetBinContent(1));
    for (unsigned uBin = 1; uBin < NumBins; uBin++)
    {
        CummDistr.SetBinContent(uBin, CummDistr.GetBinContent(uBin - 1) + Projection->GetBinContent(uBin + 1));
    }

    double *MomBins = new double[NumBins + 1];
    MomBins[0] = 0;
    // MomBins[1] = CummDistr.GetBinContent(unsigned(0));
    // MomBins[NumBins] = 1;
    for (unsigned uBin = 1; uBin <= NumBins; uBin++)
    {
        MomBins[uBin] = (CummDistr.GetBinContent(uBin) + CummDistr.GetBinContent(uBin - 1)) * 0.5;
        if (MomBins[uBin] <= MomBins[uBin - 1])
            MomBins[uBin] = MomBins[uBin - 1] + 1e-6;
        // printf("MomBins[%u] = %e\n",uBin,MomBins[uBin]);
    }
    // usleep(2e6);
    DLM_Histo<double> *Result = new DLM_Histo<double>();
    Result->SetUp(1);
    Result->SetUp(0, NumBins, MomBins);
    Result->Initialize();
    // printf("Result->Initialize\n");
    // usleep(2e6);
    for (unsigned uBin = 0; uBin < NumBins; uBin++)
    {
        if (!Projection)
            break;
        Result->SetBinCenter(0, uBin, CummDistr.GetBinContent(uBin)); // cum. value
        Result->SetBinContent(uBin, CummDistr.GetBinCenter(0, uBin)); // momentum value
        // printf("%u: x=%.4f; y=%.4f;\n",uBin,CummDistr.GetBinContent(uBin),CummDistr.GetBinCenter(0,uBin));
    }

    delete[] MomBins;
    delete InputFile;
    return Result;
}

void RootFile_DlmCk(const TString &RootFileName, const TString &GraphName, DLM_Ck *CkToPlot)
{
    TFile *RootFile = new TFile(RootFileName, "update");
    if (!RootFile)
        RootFile = new TFile(RootFileName, "recreate");
    const unsigned NumBins = CkToPlot->GetNbins();
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    for (unsigned uBin = 0; uBin < NumBins; uBin++)
    {
        graph.SetPoint(uBin, CkToPlot->GetBinCenter(0, uBin), CkToPlot->GetBinContent(uBin));
    }
    graph.Write("", TObject::kOverwrite);
    delete RootFile;
}
void RootFile_DlmCk(const TString &RootFileName, const TString &GraphName, CATS *Kitty)
{
    DLM_Ck dlmck(1, 0, *Kitty);
    RootFile_DlmCk(RootFileName, GraphName, &dlmck);
}

void RootFile_DlmSource(const TString &RootFileName, const TString &GraphName, CATS *Kitty, const unsigned &NumBins, const double &rMin, const double &rMax, const double &lambda, const bool &FourPi)
{
    TFile *RootFile = new TFile(RootFileName, "update");
    if (!RootFile)
        RootFile = new TFile(RootFileName, "recreate");
    TGraph graph;
    graph.SetName(GraphName);
    graph.Set(NumBins);
    const double rWidth = (rMax - rMin) / double(NumBins);
    for (unsigned uBin = 0; uBin < NumBins; uBin++)
    {
        double RAD = rWidth * 0.5 + double(uBin) * rWidth;
        if (FourPi)
            graph.SetPoint(uBin, RAD, Kitty->EvaluateTheSource(10, RAD, 0) * lambda);
        else
            graph.SetPoint(uBin, RAD, Kitty->EvaluateTheSource(10, RAD, 0) * lambda / (4. * Pi * (RAD + 1e-6) * (RAD + 1e-6)));
    }
    graph.Write("", TObject::kOverwrite);
    delete RootFile;
}

double Get_reff(TH1F *hsource, const float lambda, const float CEI)
{
    TF1 *fit_temp;
    double result = Get_reff_TF1(hsource, fit_temp, lambda, CEI);
    delete fit_temp;
    return result;
}
double Get_reff_TF1(TH1F *hsource, TF1 *&fsource, const float lambda, const float CEI)
{
    TH1F *hfit4325 = (TH1F *)hsource->Clone("hfit4325");
    hfit4325->Scale(1. / hfit4325->Integral(), "width");

    double lowerlimit;
    double upperlimit;
    GetCentralInterval(*hfit4325, CEI, lowerlimit, upperlimit, true);

    fsource = new TF1("fit4325", "[0]*4.*TMath::Pi()*x*x*pow(4.*TMath::Pi()*[1]*[1],-1.5)*exp(-(x*x)/(4.*[1]*[1]))+1.-[0]", lowerlimit, upperlimit);
    fsource->FixParameter(0, lambda);
    fsource->SetParameter(1, hfit4325->GetMean() / 2.3);
    fsource->SetParLimits(1, hfit4325->GetMean() / 10., hfit4325->GetMean() * 2.);

    hfit4325->Fit(fsource, "Q, S, N, R, M");

    // fsource = fit4325;
    delete hfit4325;

    return fsource->GetParameter(1);
}

double GetRcore(DLM_CleverMcLevyResoTM &MagicSource, const double &reff)
{
    const unsigned NumRadBins = 256;
    const double rMin = 0;
    const double rMax = 32;
    double rcore = reff;
    const double Eps = reff * 0.001;
    TH1F *h_rstar = new TH1F("fgsgsdfg", "fgsgsdfg", NumRadBins, rMin, rMax);
    double reff_fit;
    while (true)
    {
        for (unsigned uBin = 0; uBin < NumRadBins; uBin++)
        {
            double rstar = h_rstar->GetBinCenter(uBin + 1);
            double parameters[2];
            parameters[0] = rcore;
            parameters[1] = 2.0;
            double val;
            val = MagicSource.RootEval(&rstar, parameters);
            h_rstar->SetBinContent(uBin + 1, val);
            h_rstar->SetBinError(uBin + 1, 1e-3);
        }
        reff_fit = Get_reff(h_rstar);
        if (fabs(reff_fit - reff) > Eps)
        {
            if (reff_fit > reff)
                rcore *= (reff / reff_fit);
        }
        else
            break;
    }
    delete h_rstar;
    return rcore;
}
double GetReff(DLM_CleverMcLevyResoTM &MagicSource, const double &rcore)
{
    const unsigned NumRadBins = 256;
    const double rMin = 0;
    const double rMax = 32;
    double reff;
    TH1F *h_rstar = new TH1F("fgsgsdfg", "fgsgsdfg", NumRadBins, rMin, rMax);

    for (unsigned uBin = 0; uBin < NumRadBins; uBin++)
    {
        double rstar = h_rstar->GetBinCenter(uBin + 1);
        double parameters[2];
        parameters[0] = rcore;
        parameters[1] = 2.0;
        double val;
        val = MagicSource.RootEval(&rstar, parameters);
        h_rstar->SetBinContent(uBin + 1, val);
        h_rstar->SetBinError(uBin + 1, 1e-3);
    }
    reff = Get_reff(h_rstar);

    delete h_rstar;
    return reff;
}


double ConvertMeanGauss(double mean){
  TF1* fGauss = new TF1("fGauss",GaussSourceTF1,0,64,1);

  unsigned Trials = 1000000;
  double Expected = mean/2.3;
  double MinR = Expected/3;
  double MaxR = Expected*3;
  double Step = (MaxR-MinR)/double(Trials);
  double test_val;

  double best_result = 0;
  double best_nsig = 1e6;
  double current_result;
  double current_nsig;

  bool GettingWorse_up = false;
  bool GettingWorse_low = false;

  for(unsigned uT=0; uT<Trials/2; uT++){
    test_val = Expected - 0.5*Step - double(uT)*Step;
    fGauss->SetParameter(0,test_val);
    current_result = fGauss->Mean(0,64);
    current_nsig = fabs(current_result-mean);
    if(current_nsig<best_nsig){
      best_nsig = current_nsig;
      best_result = current_result;
      GettingWorse_low = false;
    }
    else GettingWorse_low = true;


    test_val = Expected + 0.5*Step + double(uT)*Step;
    fGauss->SetParameter(0,test_val);
    current_result = fGauss->Mean(0,64);
    current_nsig = fabs(current_result-mean);
    if(current_nsig<best_nsig){
      best_nsig = current_nsig;
      best_result = current_result;
      GettingWorse_up = false;
    }
    else GettingWorse_up = true;

    if(GettingWorse_up && GettingWorse_low) break;
  }

  delete fGauss;
  return best_result;
}
double ConvertGaussMean(double reff){
  TF1* fGauss = new TF1("fGauss",GaussSourceTF1,0,reff*10,1);
  fGauss->SetParameter(0,reff);
  double mean = fGauss->Mean(0,reff*10);
  delete fGauss;
  return mean;
}

void SetUp_RSM_Flat(DLM_CleverMcLevyResoTM &MagicSource, const double frac1, const double frac2, const double MassP1, const double MassP2,
                    const double MassR1, const double MassR2, const double TauR1, const double TauR2,
                    const double TIMEOUT, const int flag)
{

    // tested that changing this by 50% makes NO difference!
    const double MomSpread = 850 * 1;
    const double Disp = 0.8;

    BasicSetUp_MS(MagicSource, frac1, frac2);

    TREPNI Database(0);
    Database.SetSeed(11);
    std::vector<TreParticle *> ParticleList;
    ParticleList.push_back(Database.NewParticle("Part1"));
    ParticleList.push_back(Database.NewParticle("Part2"));
    ParticleList.push_back(Database.NewParticle("Pion"));

    if (frac1)
        ParticleList.push_back(Database.NewParticle("Reso1"));
    if (frac2)
        ParticleList.push_back(Database.NewParticle("Reso2"));

    for (TreParticle *prt : ParticleList)
    {
        if (prt->GetName() == "Part1")
        {
            prt->SetMass(MassP1);
            prt->SetAbundance((1. - frac1));
        }
        else if (prt->GetName() == "Part2")
        {
            prt->SetMass(MassP2);
            prt->SetAbundance((1. - frac2));
        }
        else if (prt->GetName() == "Pion")
        {
            prt->SetMass(Mass_pic);
            prt->SetAbundance(0);
        }
        else if (prt->GetName() == "Reso1" && frac1)
        {
            prt->SetMass(MassR1);
            prt->SetAbundance(frac1);
            prt->SetWidth(hbarc / TauR1);

            prt->NewDecay();
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Part1"));
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
            prt->GetDecay(0)->SetBranching(100);
        }
        else if (prt->GetName() == "Reso2" && frac2)
        {
            prt->SetMass(MassR2);
            prt->SetAbundance(frac2);
            prt->SetWidth(hbarc / TauR2);

            prt->NewDecay();
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Part2"));
            prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
            prt->GetDecay(0)->SetBranching(100);
        }
        prt->SetPtPz(prt->GetMass() * MomSpread * 0.001, prt->GetMass() * MomSpread * 0.001);
    }

    std::vector<std::string> ListOfParticles;
    ListOfParticles.push_back("Part1");
    ListOfParticles.push_back("Part2");
    CECA Ivana(Database, ListOfParticles);
    Ivana.SetTargetStatistics(10);
    Ivana.SetEventMult(2);
    Ivana.SetSourceDim(2);
    Ivana.SetThreadTimeout(TIMEOUT);
    Ivana.SetFemtoRegion(100);
    Ivana.EqualizeFsiTime(true);
    Ivana.SetPropagateMother(false);
    Ivana.GHETTO_EVENT = true;

    if (flag % 10 == 0)
    {
        Ivana.SetDisplacement(Disp);
        Ivana.SetHadronization(0);
    }
    else
    {
        Ivana.SetDisplacement(0);
        Ivana.SetHadronization(Disp);
    }

    Ivana.SetUp_RSM = &MagicSource;
    Ivana.GoBabyGo(0);
}

void BasicSetUp_MS(DLM_CleverMcLevyResoTM &MagicSource, double frac1, double frac2)
{
    MagicSource.InitStability(1, 2 - 1e-6, 2 + 1e-6);
    MagicSource.InitScale(38, 0.15, 2.0);
    MagicSource.InitRad(257 * 2, 0, 64);
    MagicSource.InitType(2);
    MagicSource.SetUpReso(0, frac1);
    MagicSource.SetUpReso(1, frac2);
    MagicSource.InitNumMcIter(1000000);
}

void SetUp_RSMflat_pp(DLM_CleverMcLevyResoTM &MagicSource)
{

    const double frac1 = 0.6422;
    const double frac2 = 0.6422;
    BasicSetUp_MS(MagicSource, frac1, frac2);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;
    const double q_CutOff = 200;

    SetUp_RSM_Flat(MagicSource, frac1, frac2, Mass_p, Mass_p,
                   Mass_ProtonReso, Mass_ProtonReso, Tau_Proton, Tau_Proton,
                   16, 0);
}

// const double FracProtonReso = 0.6422*1;
// const double FracLambdaReso = 0.6438*1;
//FLAG = 0 is the CORRECT SIGN !!!
//FLAG = 1 is the WRONG SIGN
void SetUp_RSM_pp(DLM_CleverMcLevyResoTM &MagicSource, const TString InputFolder, const int flag)
{
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

    const double frac1 = 0.6422;
    const double frac2 = 0.6422;
    BasicSetUp_MS(MagicSource, frac1, frac2);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;
    const double q_CutOff = 200;

    TFile *F_EposDisto_p_pReso;
    F_EposDisto_p_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_pReso.root", InputFolder.Data()));
    TNtuple *T_EposDisto_p_pReso = (TNtuple *)F_EposDisto_p_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_pReso = T_EposDisto_p_pReso->GetEntries();
    T_EposDisto_p_pReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_pReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_pReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_pReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_pReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_pReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_pReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_pReso; uEntry++)
    {
        T_EposDisto_p_pReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = Tau_Proton;
        fM2 = Mass_ProtonReso;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        // if(flag/10==0){
        // this is with the correct sign
        //(man, I think its flipped, this should be the wrong sign)
        if (flag % 10 == 1)
        {
            MagicSource.AddBGT_PR(RanVal1, -cos(AngleRcP2));
            MagicSource.AddBGT_RP(RanVal1, cos(AngleRcP2));
        }
        else
        {
            MagicSource.AddBGT_PR(RanVal1, cos(AngleRcP2));
            MagicSource.AddBGT_RP(RanVal1, -cos(AngleRcP2));
        }
    }
    delete F_EposDisto_p_pReso;

    TFile *F_EposDisto_pReso_pReso;
    F_EposDisto_pReso_pReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_pReso.root", InputFolder.Data()));
    TNtuple *T_EposDisto_pReso_pReso = (TNtuple *)F_EposDisto_pReso_pReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_pReso = T_EposDisto_pReso_pReso->GetEntries();
    T_EposDisto_pReso_pReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_pReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_pReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_pReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_pReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_pReso; uEntry++)
    {
        T_EposDisto_pReso_pReso->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = Tau_Proton;
        fM1 = Mass_ProtonReso;
        fM2 = Mass_ProtonReso;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource.AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_pReso;
}

void SetUp_RSMflat_pL(DLM_CleverMcLevyResoTM &MagicSource)
{

    const double frac1 = 0.6422;
    const double frac2 = 0.6438;
    BasicSetUp_MS(MagicSource, frac1, frac2);
    const double Tau_Proton = 1.65;
    const double Tau_Lambda = 4.69; // 4.69
    const double Mass_ProtonReso = 1362;
    const double Mass_LambdaReso = 1462;

    SetUp_RSM_Flat(MagicSource, frac1, frac2, Mass_p, Mass_L,
                   Mass_ProtonReso, Mass_LambdaReso, Tau_Proton, Tau_Lambda,
                   16, 0);
}

void SetUp_RSM_pL(DLM_CleverMcLevyResoTM &MagicSource, const TString InputFolder, const int flag)
{
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

    BasicSetUp_MS(MagicSource, 0.6422, 0.6438);
    const double Tau_Proton = 1.65;
    const double Tau_Lambda = 4.69; // 4.69
    const double Mass_ProtonReso = 1362;
    const double Mass_LambdaReso = 1462;
    const double q_CutOff = 200;

    TFile *F_EposDisto_p_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_p_LamReso.root", InputFolder.Data()));
    TNtuple *T_EposDisto_p_LamReso = (TNtuple *)F_EposDisto_p_LamReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_p_LamReso = T_EposDisto_p_LamReso->GetEntries();
    T_EposDisto_p_LamReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_p_LamReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_p_LamReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_p_LamReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_p_LamReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_p_LamReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_p_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_p_LamReso; uEntry++)
    {
        T_EposDisto_p_LamReso->GetEntry(uEntry);
        Tau1 = 0;
        Tau2 = Tau_Lambda;
        fM2 = Mass_LambdaReso;
        if (k_D > q_CutOff)
            continue;
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource.AddBGT_PR(RanVal2, cos(AngleRcP2));
    }
    delete F_EposDisto_p_LamReso;

    TFile *F_EposDisto_pReso_Lam = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Lam.root", InputFolder.Data()));
    TNtuple *T_EposDisto_pReso_Lam = (TNtuple *)F_EposDisto_pReso_Lam->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Lam = T_EposDisto_pReso_Lam->GetEntries();
    T_EposDisto_pReso_Lam->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Lam->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Lam->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Lam->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Lam->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Lam; uEntry++)
    {
        T_EposDisto_pReso_Lam->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = 0;
        fM1 = Mass_ProtonReso;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource.AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete F_EposDisto_pReso_Lam;

    TFile *F_EposDisto_pReso_LamReso = new TFile(TString::Format("%s/CatsFiles/Source/EposAngularDist/EposDisto_pReso_LamReso.root", InputFolder.Data()));
    TNtuple *T_EposDisto_pReso_LamReso = (TNtuple *)F_EposDisto_pReso_LamReso->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_LamReso = T_EposDisto_pReso_LamReso->GetEntries();
    T_EposDisto_pReso_LamReso->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_LamReso->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_LamReso; uEntry++)
    {
        T_EposDisto_pReso_LamReso->GetEntry(uEntry);
        Tau1 = Tau_Proton;
        Tau2 = Tau_Lambda;
        fM1 = Mass_ProtonReso;
        fM2 = Mass_LambdaReso;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        RanVal2 = RanGen.Exponential(fM2 / (fP2 * Tau2));
        MagicSource.AddBGT_RR(RanVal1, cos(AngleRcP1), RanVal2, cos(AngleRcP2), cos(AngleP1P2));
    }
    delete F_EposDisto_pReso_LamReso;
}

void SetUp_RSMflat_pXi(DLM_CleverMcLevyResoTM &MagicSource)
{

    const double frac1 = 0.6422;
    const double frac2 = 0.0;
    BasicSetUp_MS(MagicSource, frac1, frac2);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;

    SetUp_RSM_Flat(MagicSource, frac1, frac2, Mass_p, Mass_Xim,
                   Mass_ProtonReso, 0, Tau_Proton, 0,
                   16, 0);
}

void SetUp_RSM_pXi(DLM_CleverMcLevyResoTM &MagicSource, const TString InputFolder, const int flag)
{

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

    BasicSetUp_MS(MagicSource, 0.6422, 0.0);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;
    const double q_CutOff = 200;

    TFile *F_EposDisto_pReso_Xim = new TFile(InputFolder + "/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Xim.root");
    TNtuple *T_EposDisto_pReso_Xim = (TNtuple *)F_EposDisto_pReso_Xim->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Xim = T_EposDisto_pReso_Xim->GetEntries();
    T_EposDisto_pReso_Xim->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Xim->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Xim->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Xim->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Xim->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Xim->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Xim->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Xim->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Xim; uEntry++)
    {
        T_EposDisto_pReso_Xim->GetEntry(uEntry);
        Tau1 = 1.65;
        Tau2 = 0;
        fM1 = 1362;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource.AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete F_EposDisto_pReso_Xim;
}

void SetUp_RSMflat_pOmega(DLM_CleverMcLevyResoTM &MagicSource)
{

    const double frac1 = 0.6422;
    const double frac2 = 0.0;
    BasicSetUp_MS(MagicSource, frac1, frac2);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;

    SetUp_RSM_Flat(MagicSource, frac1, frac2, Mass_p, MassOmega,
                   Mass_ProtonReso, 0, Tau_Proton, 0,
                   16, 0);
}

void SetUp_RSM_pOmega(DLM_CleverMcLevyResoTM &MagicSource, const TString InputFolder, const int flag)
{

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

    BasicSetUp_MS(MagicSource, 0.6422, 0.0);
    const double Tau_Proton = 1.65;
    const double Mass_ProtonReso = 1362;
    const double q_CutOff = 200;

    TFile *F_EposDisto_pReso_Omega = new TFile(InputFolder + "/CatsFiles/Source/EposAngularDist/EposDisto_pReso_Omega.root");
    TNtuple *T_EposDisto_pReso_Omega = (TNtuple *)F_EposDisto_pReso_Omega->Get("InfoTuple_ClosePairs");
    unsigned N_EposDisto_pReso_Omega = T_EposDisto_pReso_Omega->GetEntries();
    T_EposDisto_pReso_Omega->SetBranchAddress("k_D", &k_D);
    T_EposDisto_pReso_Omega->SetBranchAddress("P1", &fP1);
    T_EposDisto_pReso_Omega->SetBranchAddress("P2", &fP2);
    T_EposDisto_pReso_Omega->SetBranchAddress("M1", &fM1);
    T_EposDisto_pReso_Omega->SetBranchAddress("M2", &fM2);
    T_EposDisto_pReso_Omega->SetBranchAddress("Tau1", &Tau1);
    T_EposDisto_pReso_Omega->SetBranchAddress("Tau2", &Tau2);
    T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP1", &AngleRcP1);
    T_EposDisto_pReso_Omega->SetBranchAddress("AngleRcP2", &AngleRcP2);
    T_EposDisto_pReso_Omega->SetBranchAddress("AngleP1P2", &AngleP1P2);
    for (unsigned uEntry = 0; uEntry < N_EposDisto_pReso_Omega; uEntry++)
    {
        T_EposDisto_pReso_Omega->GetEntry(uEntry);
        Tau1 = 1.65;
        Tau2 = 0;
        fM1 = 1362;
        if (k_D > q_CutOff)
            continue;
        RanVal1 = RanGen.Exponential(fM1 / (fP1 * Tau1));
        MagicSource.AddBGT_RP(RanVal1, cos(AngleRcP1));
    }
    delete F_EposDisto_pReso_Omega;
}

bool NormalizeSource_rk(DLM_Histo<float> *dlmSource)
{
    if (dlmSource == NULL)
        return false;
    if (dlmSource->GetDim() != 2)
        return false;

    for (unsigned uKstar = 0; uKstar < dlmSource->GetNbins(0); uKstar++)
    {
        double TotInt = 0;
        double BinVal, BinErr;
        double BinWidth;
        for (unsigned uRstar = 0; uRstar < dlmSource->GetNbins(1); uRstar++)
        {
            // scale to the bin width
            BinWidth = dlmSource->GetBinSize(1, uRstar);
            dlmSource->SetBinContent(uKstar, uRstar, dlmSource->GetBinContent(uKstar, uRstar) / BinWidth);
            dlmSource->SetBinError(uKstar, uRstar, dlmSource->GetBinError(uKstar, uRstar) / BinWidth);

            // add to the integral
            BinVal = dlmSource->GetBinContent(uKstar, uRstar);
            TotInt += BinVal;
        }
        for (unsigned uRstar = 0; uRstar < dlmSource->GetNbins(1); uRstar++)
        {
            BinVal = dlmSource->GetBinContent(uKstar, uRstar);
            BinErr = dlmSource->GetBinError(uKstar, uRstar);
            dlmSource->SetBinContent(uKstar, uRstar, dlmSource->GetBinContent(uKstar, uRstar) / TotInt);
            dlmSource->SetBinError(uKstar, uRstar, dlmSource->GetBinError(uKstar, uRstar) / TotInt);
        }
    }

    return true;
}


void SetUpSplPars(TF1*& fitfun){
  fitfun = new TF1("fitfun",DlmTSplineFitPositive,0,8,3+2*10);
  fitfun->FixParameter(0,10);
  fitfun->FixParameter(1,0);
  fitfun->FixParameter(2,0);
  fitfun->FixParameter(3,0.0);
  fitfun->FixParameter(4,0.4);
  fitfun->FixParameter(5,0.8);
  fitfun->FixParameter(6,1.2);
  fitfun->FixParameter(7,1.8);
  fitfun->FixParameter(8,2.4);
  fitfun->FixParameter(9,3.2);
  fitfun->FixParameter(10,4.5);
  fitfun->FixParameter(11,6.0);
  fitfun->FixParameter(12,8.0);
}



//Mode = 0: default
//Mode = 1: reduced: we use only 8 of the Disots, making bigger distance and no allowing sharp peaks (small stdv)
//Mode = 2: reduced+: only 6 Distos
//Mode = 3: reduced+: only 6 Distos, but extra loose limits for long-range distos
void SetUpKdpPars(TF1*& fitfun, int Mode, double KdpFitMax){

  fitfun = new TF1("fKdpPars",PoissonSum,0,KdpFitMax,29);

  if(Mode==1){
    //P0
    fitfun->SetParameter(0,0.4);//mean
    fitfun->SetParLimits(0,0.2,0.6);
    fitfun->SetParameter(1,0.2);//stdv
    fitfun->SetParLimits(1,0.1,0.4);
    fitfun->SetParameter(2,1./10.);//wght
    fitfun->SetParLimits(2,0.0,1.0);

    //P1
    fitfun->SetParameter(3,0.8);//mean
    fitfun->SetParLimits(3,0.6,1.0);
    fitfun->SetParameter(4,0.4);//stdv
    fitfun->SetParLimits(4,0.2,0.8);
    fitfun->SetParameter(5,1./9.);//wght
    fitfun->SetParLimits(5,0.0,1.0);

    //P2
    fitfun->SetParameter(6,1.2);//mean
    fitfun->SetParLimits(6,0.6,1.5);
    fitfun->SetParameter(7,0.6);//stdv
    fitfun->SetParLimits(7,0.3,1.2);
    fitfun->SetParameter(8,1./8.);//wght
    fitfun->SetParLimits(8,0.0,1.0);

    //P3
    fitfun->SetParameter(9,1.8);//mean
    fitfun->SetParLimits(9,1.5,2.1);
    fitfun->SetParameter(10,0.9);//stdv
    fitfun->SetParLimits(10,0.45,1.8);
    fitfun->SetParameter(11,1./7.);//wght
    fitfun->SetParLimits(11,0.0,1.0);

    //P4
    fitfun->SetParameter(12,2.4);//mean
    fitfun->SetParLimits(12,2.1,3.0);
    fitfun->SetParameter(13,1.2);//stdv
    fitfun->SetParLimits(13,0.6,2.4);
    fitfun->SetParameter(14,1./6.);//wght
    fitfun->SetParLimits(14,0.0,1.0);

    //P5
    fitfun->SetParameter(15,3.6);//mean
    fitfun->SetParLimits(15,3.0,4.5);
    fitfun->SetParameter(16,1.8);//stdv
    fitfun->SetParLimits(16,0.9,3.6);
    fitfun->SetParameter(17,1./5.);//wght
    fitfun->SetParLimits(17,0.0,1.0);

    //P6
    fitfun->SetParameter(18,5.4);//mean
    fitfun->SetParLimits(18,4.5,6.2);
    fitfun->SetParameter(19,2.7);//stdv
    fitfun->SetParLimits(19,1.35,5.4);
    fitfun->SetParameter(20,1./4.);//wght
    fitfun->SetParLimits(20,0.0,1.0);

    //P7
    fitfun->SetParameter(21,7.0);//mean
    fitfun->SetParLimits(21,6.2,8.0);
    fitfun->SetParameter(22,3.5);//stdv
    fitfun->SetParLimits(22,1.75,7.0);
    fitfun->SetParameter(23,1./3.);//wght
    fitfun->SetParLimits(23,0.0,1.0);

    //P8
    fitfun->FixParameter(24,0);//mean
    fitfun->FixParameter(25,0);//stdv
    fitfun->FixParameter(26,0);//wght

    //P9
    fitfun->FixParameter(27,0);//mean
    fitfun->FixParameter(28,0);//stdv
  }
  else if(Mode==2){
    //P0
    fitfun->SetParameter(0,0.45);//mean
    fitfun->SetParameter(3,0.9);//mean
    fitfun->SetParameter(6,1.4);//mean
    fitfun->SetParameter(9,2.0);//mean
    fitfun->SetParameter(12,3.2);//mean
    fitfun->SetParameter(15,5.4);//mean

    fitfun->SetParLimits(0,0.225,0.9);
    for(unsigned uP=1; uP<5; uP++){
      double low_lim = 0.55*(fitfun->GetParameter(-3+uP*3)+fitfun->GetParameter(0+uP*3));
      double up_lim = 0.45*(fitfun->GetParameter(3+uP*3)+fitfun->GetParameter(0+uP*3));
      fitfun->SetParLimits(0+uP*3,low_lim,up_lim);
    }
    fitfun->SetParLimits(15,4.3,6.5);

    for(unsigned uP=0; uP<6; uP++){
      double mu = fitfun->GetParameter(0+uP*3);
      fitfun->SetParameter(1+uP*3,mu*0.4);
      fitfun->SetParLimits(1+uP*3,mu*0.3,mu*0.45);
      //fitfun->FixParameter(1+uP*3,mu*0.4);
    }

    for(unsigned uP=0; uP<5; uP++){
      fitfun->SetParameter(2+uP*3,1./(6.-double(uP)));
      fitfun->SetParLimits(2+uP*3,0,1.0);
    }
    fitfun->FixParameter(17,1);

    //P6
    fitfun->FixParameter(18,0);//mean
    fitfun->FixParameter(19,0);//stdv
    fitfun->FixParameter(20,0);//wght

    //P7
    fitfun->FixParameter(21,0);//mean
    fitfun->FixParameter(22,0);//stdv
    fitfun->FixParameter(23,0);//wght

    //P8
    fitfun->FixParameter(24,0);//mean
    fitfun->FixParameter(25,0);//stdv
    fitfun->FixParameter(26,0);//wght

    //P9
    fitfun->FixParameter(27,0);//mean
    fitfun->FixParameter(28,0);//stdv
  }
  else if(Mode==3){
    //P0
    fitfun->SetParameter(0,0.5);//mean
    fitfun->SetParameter(3,1.0);//mean
    fitfun->SetParameter(6,2.0);//mean
    fitfun->SetParameter(9,4.0);//mean
    fitfun->SetParameter(12,8.0);//mean
    fitfun->SetParameter(15,16.0);//mean


    for(unsigned uP=0; uP<6; uP++){
        double low_mean = fitfun->GetParameter(uP*3)*0.25;
        if(low_mean < 0.1) low_mean = 0.1;
        fitfun->SetParLimits(uP*3,low_mean,fitfun->GetParameter(uP*3)*8.);
    }
    
    for(unsigned uP=0; uP<6; uP++){
      double mu = fitfun->GetParameter(0+uP*3);
      fitfun->SetParameter(1+uP*3,mu*0.4);
      fitfun->SetParLimits(1+uP*3,mu*0.1,mu*1.0);
    }

    for(unsigned uP=0; uP<5; uP++){
      fitfun->SetParameter(2+uP*3,1./(6.-double(uP)));
      fitfun->SetParLimits(2+uP*3,0,1.0);
    }
    fitfun->FixParameter(17,1);

    //P6
    fitfun->FixParameter(18,0);//mean
    fitfun->FixParameter(19,0);//stdv
    fitfun->FixParameter(20,0);//wght

    //P7
    fitfun->FixParameter(21,0);//mean
    fitfun->FixParameter(22,0);//stdv
    fitfun->FixParameter(23,0);//wght

    //P8
    fitfun->FixParameter(24,0);//mean
    fitfun->FixParameter(25,0);//stdv
    fitfun->FixParameter(26,0);//wght

    //P9
    fitfun->FixParameter(27,0);//mean
    fitfun->FixParameter(28,0);//stdv
  }

  else{
    //P0
    fitfun->SetParameter(0,0.3);//mean
    fitfun->SetParLimits(0,0.15,0.45);
    fitfun->SetParameter(1,0.15);//stdv
    fitfun->SetParLimits(1,0.0,0.3);
    fitfun->SetParameter(2,1./10.);//wght
    fitfun->SetParLimits(2,0.0,1.0);

    //P1
    fitfun->SetParameter(3,0.6);//mean
    fitfun->SetParLimits(3,0.45,0.75);
    fitfun->SetParameter(4,0.15);//stdv
    fitfun->SetParLimits(4,0.0,0.3);
    fitfun->SetParameter(5,1./9.);//wght
    fitfun->SetParLimits(5,0.0,1.0);

    //P2
    fitfun->SetParameter(6,0.9);//mean
    fitfun->SetParLimits(6,0.75,1.05);
    fitfun->SetParameter(7,0.15);//stdv
    fitfun->SetParLimits(7,0.0,0.3);
    fitfun->SetParameter(8,1./8.);//wght
    fitfun->SetParLimits(8,0.0,1.0);

    //P3
    fitfun->SetParameter(9,1.2);//mean
    fitfun->SetParLimits(9,1.05,1.5);
    fitfun->SetParameter(10,0.3);//stdv
    fitfun->SetParLimits(10,0.0,0.6);
    fitfun->SetParameter(11,1./7.);//wght
    fitfun->SetParLimits(11,0.0,1.0);

    //P4
    fitfun->SetParameter(12,1.8);//mean
    fitfun->SetParLimits(12,1.5,2.1);
    fitfun->SetParameter(13,0.3);//stdv
    fitfun->SetParLimits(13,0.0,0.6);
    fitfun->SetParameter(14,1./6.);//wght
    fitfun->SetParLimits(14,0.0,1.0);

    //P5
    fitfun->SetParameter(15,2.4);//mean
    fitfun->SetParLimits(15,2.1,2.8);
    fitfun->SetParameter(16,0.4);//stdv
    fitfun->SetParLimits(16,0.0,0.8);
    fitfun->SetParameter(17,1./5.);//wght
    fitfun->SetParLimits(17,0.0,1.0);

    //P6
    fitfun->SetParameter(18,3.2);//mean
    fitfun->SetParLimits(18,2.8,3.6);
    fitfun->SetParameter(19,0.4);//stdv
    fitfun->SetParLimits(19,0.0,0.8);
    fitfun->SetParameter(20,1./4.);//wght
    fitfun->SetParLimits(20,0.0,1.0);

    //P7
    fitfun->SetParameter(21,4.0);//mean
    fitfun->SetParLimits(21,3.6,4.6);
    fitfun->SetParameter(22,0.6);//stdv
    fitfun->SetParLimits(22,0.0,1.2);
    fitfun->SetParameter(23,1./3.);//wght
    fitfun->SetParLimits(23,0.0,1.0);

    //P8
    fitfun->SetParameter(24,5.2);//mean
    fitfun->SetParLimits(24,4.6,6.1);
    fitfun->SetParameter(25,0.9);//stdv
    fitfun->SetParLimits(25,0.0,1.8);
    fitfun->SetParameter(26,1./2.);//wght
    fitfun->SetParLimits(26,0.0,1.0);

    //P9
    fitfun->SetParameter(27,7.0);//mean
    fitfun->SetParLimits(27,6.1,8.5);
    fitfun->SetParameter(28,1.5);//stdv
    fitfun->SetParLimits(28,0.0,3.0);

  }

  fitfun->SetNpx(4096);

}

//! N.B. the femto sign convension used (+ scat len = attraction)
//Coulomb case, following Nuclear Physics A253 (1975) 355--364
//the return value is the k*cotg(delta_(c))
//i.e. to get the phase shifts: atan(kstar/EffRangeExp)
//x = kstar (MeV)
//[0] = q1q2 * RedMass (MeV) -> if zero we go back to the default ERE
//[1] = Coulomb corrected inversed scattering length (1/fm)
//[2] = Coulomb corrected effective range (fm)
//[3,4] = the two extra pars in the expansion after effective range (in NU since we dont really care about the values)
double EffRangeExp(double* x, double* par){
    double& kstar = *x;
    double& q1q2rm = par[0];
    double invf_cs = par[1]/FmToNu;
    double d_cs = par[2]*FmToNu;
    double& ere4 = par[3];
    double& ere6 = par[4];

    //this is Eta = RedMass*double(q1q2)*AlphaFS/Momentum, aka Sommerfeld parameter (sometimes called gamma);
    double Eta;
    double hgam;
    //the precision to which we would like to eval the function h(gamma)
    const double eps = 0.0001;
    double precision;
    int nstep = 1;

    if(q1q2rm){
        Eta = AlphaFS*q1q2rm/kstar;
        hgam = -EulerConst - log(Eta);
        precision = 1;
    }
    else{
        Eta = 0;
        hgam = 0;
        precision = 0;
    }

//printf("kstar = %.3e\n",kstar);
//printf("q1q2rm = %.3e\n",q1q2rm);
//printf("Eta = %.3e\n",Eta);
//printf("hgam = %.3e\n",hgam);

    while(precision>eps && nstep<32){
        double dn = nstep;
        double dh = Eta*Eta/(dn*(dn*dn+Eta*Eta));
        hgam += dh;
        //this says we have to take at least the first two terms of the sum
        if(nstep>1) precision = fabs(dh/hgam);
        nstep++;
    }
//printf("HGAM = %.3e\n",hgam);
    double PreFactor;
    double AddFactor;
    if(q1q2rm){
        PreFactor = (exp(2.*Pi*Eta)-1.)/2.*Pi*Eta;
        AddFactor = -2.*kstar*Eta*hgam;
    }
    else{
        PreFactor = 1;
        AddFactor = 0;
    }
//printf("PreFactor = %.3e\n",PreFactor);
//printf("AddFactor = %.3e\n",AddFactor);
    return PreFactor*(AddFactor + invf_cs + 0.5*d_cs*pow(kstar,2.) + ere4*pow(kstar,4.) + ere6*pow(kstar,6.));
}

bool GetScattParameters(CATS& Kitty, double& ScatLen, double& EffRan, TH1F*& hFit, TF1*& fitSP,
  const int Nterms, const bool Fixf0, const bool Fixd0, const unsigned short usCh){
  Kitty.KillTheCat();
  double* MomBins = Kitty.CopyMomBin();
  hFit = new TH1F("hFit52351","hFit52351",Kitty.GetNumMomBins(),MomBins);
  double LAST_POINT;
  double CURRENT_POINT;
  for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
    CURRENT_POINT = Kitty.GetMomentum(uMom)/tan(Kitty.GetPhaseShift(uMom,usCh,0));
    if(uMom){
      if(CURRENT_POINT*LAST_POINT<0&&fabs(CURRENT_POINT-LAST_POINT)>1000&&Kitty.GetMomentum(uMom)<120)
      {fitSP=NULL;delete[]MomBins;return false;}
    }
    hFit->SetBinContent(uMom+1,CURRENT_POINT);
    hFit->SetBinError(uMom+1,1.);
    LAST_POINT = CURRENT_POINT;
  }
  TF1* fitSP2; 
  TF1* fitSP4; 
  TF1* fitSP6;


  fitSP2 = new TF1("fitSP2",EffRangeExp,10,90,5);
  fitSP2->FixParameter(3,0);
  fitSP2->FixParameter(4,0);
  fitSP4 = new TF1("fitSP2",EffRangeExp,10,90,5);
  fitSP4->FixParameter(4,0);
  fitSP6 = new TF1("fitSP2",EffRangeExp,10,90,5);

  fitSP2->FixParameter(0,Kitty.GetQ1Q2()*Kitty.GetRedMass());
  fitSP4->FixParameter(0,Kitty.GetQ1Q2()*Kitty.GetRedMass());
  fitSP6->FixParameter(0,Kitty.GetQ1Q2()*Kitty.GetRedMass());
  

  double inv_f0 = ScatLen==0?0:1./ScatLen;
  if(Fixf0) {fitSP2->FixParameter(1,inv_f0);fitSP4->FixParameter(1,inv_f0);fitSP6->FixParameter(1,inv_f0);}
  else {fitSP2->SetParameter(1,inv_f0);fitSP2->SetParLimits(1,-100,100);
        fitSP4->SetParameter(1,inv_f0);fitSP4->SetParLimits(1,-100,100);
        fitSP6->SetParameter(1,inv_f0);fitSP6->SetParLimits(1,-100,100);}
  if(Fixd0) { fitSP2->FixParameter(2,EffRan);
              fitSP4->FixParameter(2,EffRan);
              fitSP6->FixParameter(2,EffRan);}
  else {fitSP2->SetParameter(2,EffRan);fitSP2->SetParLimits(2,-50,50);
        fitSP4->SetParameter(2,EffRan);fitSP4->SetParLimits(2,-50,50);
        fitSP6->SetParameter(2,EffRan);fitSP6->SetParLimits(2,-50,50);}
  fitSP4->SetParameter(3,0);fitSP6->SetParameter(3,0);
  fitSP6->SetParameter(4,0);

  double Chi2_Old = 1e64;

  //hFit->Fit(fitSP2, "Q, S, N, R, M");
  ROOT::Math::MinimizerOptions MinOpt;
  MinOpt.SetMinimizerType("Minuit2");
  MinOpt.SetPrintLevel(0);
  DLM_FitHisto(hFit, fitSP2, "Q, S, N, R, M", "", &MinOpt);
  //printf("f0 %f\n", 1./fitSP2->GetParameter(0));
  ScatLen = 1./fitSP2->GetParameter(1);
  EffRan = fitSP2->GetParameter(2);
  if(Nterms<=2){delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}

  //hFit->Fit(fitSP4, "Q, S, N, R, M");
  DLM_FitHisto(hFit, fitSP4, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP4->GetChisquare()/fitSP2->GetChisquare()>0.8)
  {delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}
  ScatLen = 1./fitSP4->GetParameter(1);
  EffRan = fitSP4->GetParameter(2);
  if(Nterms<=3){delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}

  //hFit->Fit(fitSP6, "Q, S, N, R, M");
  DLM_FitHisto(hFit, fitSP6, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP6->GetChisquare()/fitSP4->GetChisquare()>0.8)
  {delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}
  ScatLen = 1./fitSP6->GetParameter(1);
  EffRan = fitSP6->GetParameter(2);
  delete fitSP2; delete fitSP4; fitSP=fitSP6; delete[]MomBins; return true;
}

bool GetComplexScattParameters(CATS &KittyComplex, complex<double> &ScatLen, complex<double> &EffRan, TH1F *&hFitR, TH1F *&hFitI, TF1 *&fitSPR, TF1 *&fitSPI, const int Nterms, const bool Fixf0, const bool Fixd0, const unsigned short usCh, const unsigned short usPW)
{
    KittyComplex.KillTheCat();
    double ScatLenR = ScatLen.real();
    double ScatLenI = ScatLen.imag();
    double EffRanR = EffRan.real();
    double EffRanI = EffRan.imag();
    double* MomBins = KittyComplex.CopyMomBin();
    hFitR = new TH1F("hFit52351R", "hFit52351R", KittyComplex.GetNumMomBins(), MomBins);
    hFitI = new TH1F("hFit52351I", "hFit52351I", KittyComplex.GetNumMomBins(), MomBins);

    double LAST_POINT_R, LAST_POINT_I;
    double CURRENT_POINT_R, CURRENT_POINT_I;
    for (unsigned uMom = 0; uMom < KittyComplex.GetNumMomBins(); uMom++)
    {
        // double PhaseShiftR = 0.5 * arg(KittyComplex.EvalScatteringMatrix(KittyComplex.GetMomentum(uMom), usCh, usPW));
        // /// putting a - sign to ensure that we always have a positive imaginary phase shift as notation suggests
        // double PhaseShiftI = - 0.5 * log(abs(KittyComplex.EvalScatteringMatrix(KittyComplex.GetMomentum(uMom), usCh, usPW)));
        double PhaseShiftR = real(KittyComplex.EvalComplexPhaseShifts(KittyComplex.GetMomentum(uMom), usCh, usPW));
        double PhaseShiftI = imag(KittyComplex.EvalComplexPhaseShifts(KittyComplex.GetMomentum(uMom), usCh, usPW));
        CURRENT_POINT_R = KittyComplex.GetMomentum(uMom) / tan(PhaseShiftR);
        CURRENT_POINT_I = KittyComplex.GetMomentum(uMom) / tan(PhaseShiftI);
        if (uMom)
        {
            if (CURRENT_POINT_R * LAST_POINT_R < 0 && fabs(CURRENT_POINT_R - LAST_POINT_R) > 1000 && KittyComplex.GetMomentum(uMom) < 120)
            {
                fitSPR = NULL;
                delete[] MomBins;
                return false;
            }
            if (CURRENT_POINT_I * LAST_POINT_I < 0 && fabs(CURRENT_POINT_I - LAST_POINT_I) > 1000 && KittyComplex.GetMomentum(uMom) < 120)
            {
                fitSPI = NULL;
                delete[] MomBins;
                return false;
            }
        }
        hFitR->SetBinContent(uMom + 1, CURRENT_POINT_R);
        hFitR->SetBinError(uMom + 1, 1.);
        LAST_POINT_R = CURRENT_POINT_R;
        hFitI->SetBinContent(uMom + 1, CURRENT_POINT_I);
        hFitI->SetBinError(uMom + 1, 1.);
        LAST_POINT_I = CURRENT_POINT_I;
    }
    TF1 *fitSP2;
    TF1 *fitSP4;
    TF1 *fitSP6;
    fitSP2 = new TF1("fitSP2", EffRangeExp, 10, 90, 5);
    fitSP2->FixParameter(3, 0);
    fitSP2->FixParameter(4, 0);
    fitSP4 = new TF1("fitSP2", EffRangeExp, 10, 90, 5);
    fitSP4->FixParameter(4, 0);
    fitSP6 = new TF1("fitSP2", EffRangeExp, 10, 90, 5);
    fitSP2->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());
    fitSP4->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());
    fitSP6->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());

    TF1 *fitSP2I;
    TF1 *fitSP4I;
    TF1 *fitSP6I;
    fitSP2I = new TF1("fitSP2I", EffRangeExp, 10, 90, 5);
    fitSP2I->FixParameter(3, 0);
    fitSP2I->FixParameter(4, 0);
    fitSP4I = new TF1("fitSP2I", EffRangeExp, 10, 90, 5);
    fitSP4I->FixParameter(4, 0);
    fitSP6I = new TF1("fitSP2I", EffRangeExp, 10, 90, 5);
    fitSP2I->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());
    fitSP4I->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());
    fitSP6I->FixParameter(0, KittyComplex.GetQ1Q2() * KittyComplex.GetRedMass());

    double inv_f0R = ScatLenR == 0 ? 0 : 1. / ScatLenR;
    double inv_f0I = ScatLenI == 0 ? 0 : 1. / ScatLenI;

    if (Fixf0)
    {
        fitSP2->FixParameter(1, inv_f0R);
        fitSP4->FixParameter(1, inv_f0R);
        fitSP6->FixParameter(1, inv_f0R);

        fitSP2I->FixParameter(1, inv_f0I);
        fitSP4I->FixParameter(1, inv_f0I);
        fitSP6I->FixParameter(1, inv_f0I);
    }
    else
    {
        fitSP2->SetParameter(1, inv_f0R);
        fitSP2->SetParLimits(1, -100, 100);
        fitSP4->SetParameter(1, inv_f0R);
        fitSP4->SetParLimits(1, -100, 100);
        fitSP6->SetParameter(1, inv_f0R);
        fitSP6->SetParLimits(1, -100, 100);

        fitSP2I->SetParameter(1, inv_f0I);
        fitSP2I->SetParLimits(1, 0., 100);
        fitSP4I->SetParameter(1, inv_f0I);
        fitSP4I->SetParLimits(1, 0, 100);
        fitSP6I->SetParameter(1, inv_f0I);
        fitSP6I->SetParLimits(1, 0, 100);
    }
    if (Fixd0)
    {
        fitSP2->FixParameter(2, EffRanR);
        fitSP4->FixParameter(2, EffRanR);
        fitSP6->FixParameter(2, EffRanR);

        fitSP2I->FixParameter(2, EffRanI);
        fitSP4I->FixParameter(2, EffRanI);
        fitSP6I->FixParameter(2, EffRanI);
    }
    else
    {
        fitSP2->SetParameter(2, EffRanR);
        fitSP2->SetParLimits(2, -50, 50);
        fitSP4->SetParameter(2, EffRanR);
        fitSP4->SetParLimits(2, -50, 50);
        fitSP6->SetParameter(2, EffRanR);
        fitSP6->SetParLimits(2, -50, 50);

        fitSP2I->SetParameter(2, EffRanI);
        fitSP2I->SetParLimits(2, 0., 50);
        fitSP4I->SetParameter(2, EffRanI);
        fitSP4I->SetParLimits(2, 0, 50);
        fitSP6I->SetParameter(2, EffRanI);
        fitSP6I->SetParLimits(2, 0, 50);
    }
    fitSP4->SetParameter(3, 0);
    fitSP6->SetParameter(3, 0);
    fitSP6->SetParameter(4, 0);
    fitSP4I->SetParameter(3, 0);
    fitSP6I->SetParameter(3, 0);
    fitSP6I->SetParameter(4, 0);

    double Chi2_Old = 1e64;
    ROOT::Math::MinimizerOptions MinOpt;
    MinOpt.SetMinimizerType("Minuit2");
    MinOpt.SetPrintLevel(0);
    DLM_FitHisto(hFitR, fitSP2, "Q, S, N, R, M", "", &MinOpt);
    ScatLenR = 1. / fitSP2->GetParameter(1);
    EffRanR = fitSP2->GetParameter(2);
    // printf("f0 (R)= %f\n", ScatLenR);
    // printf("d0 (R)= %f\n", EffRanR);

    DLM_FitHisto(hFitI, fitSP2I, "Q, S, N, R, M", "", &MinOpt);
    ScatLenI = 1. / fitSP2I->GetParameter(1);
    EffRanI = fitSP2I->GetParameter(2);
    // printf("f0 (I)= %f\n", ScatLenI);
    // printf("d0 (I)= %f\n", EffRanI);
    // printf("ChiRed = %f\n", fitSP2->GetChisquare() / fitSP2->GetNDF());
    // printf("ChiRed = %f\n", fitSP2I->GetChisquare() / fitSP2I->GetNDF());
    ScatLen = complex<double>(ScatLenR, ScatLenI);
    EffRan = complex<double>(EffRanR, EffRanI);

    if (Nterms <= 2)
    {
        delete fitSP4;
        delete fitSP6;
        delete fitSP4I;
        delete fitSP6I;
        fitSPR = fitSP2;
        fitSPI = fitSP2I;
        delete[] MomBins;
        return true;
    }

    DLM_FitHisto(hFitR, fitSP4, "Q, S, N, R, M", "", &MinOpt);
    if (fitSP4->GetChisquare() / fitSP2->GetChisquare() > 0.8)
    {
        delete fitSP4;
        delete fitSP6;
        fitSPR = fitSP2;
        delete[] MomBins;
        return true;
    }
    ScatLenR = 1. / fitSP4->GetParameter(1);
    EffRanR = fitSP4->GetParameter(2);
    // printf("f0 (R)= %f\n", ScatLenR);
    // printf("d0 (R)= %f\n", EffRanR);

    DLM_FitHisto(hFitI, fitSP2I, "Q, S, N, R, M", "", &MinOpt);
    if (fitSP4I->GetChisquare() / fitSP2I->GetChisquare() > 0.8)
    {
        delete fitSP4I;
        delete fitSP6I;
        fitSPI = fitSP2I;
        delete[] MomBins;
        return true;
    }
    ScatLenI = 1. / fitSP2I->GetParameter(1);
    EffRanI = fitSP2I->GetParameter(2);
    // printf("f0 (I)= %f\n", ScatLenI);
    // printf("d0 (I)= %f\n", EffRanI);

    ScatLen = complex<double>(ScatLenR, ScatLenI);
    EffRan = complex<double>(EffRanR, EffRanI);
    if (Nterms <= 3)
    {
        delete fitSP2;
        delete fitSP6;
        delete fitSP2I;
        delete fitSP6I;
        fitSPI = fitSP4I;
        delete[] MomBins;
        return true;
    }

    // hFit->Fit(fitSP6, "Q, S, N, R, M");
    DLM_FitHisto(hFitR, fitSP6, "Q, S, N, R, M", "", &MinOpt);
    if (fitSP6->GetChisquare() / fitSP4->GetChisquare() > 0.8)
    {
        delete fitSP2;
        delete fitSP6;
        fitSPR = fitSP4;
        delete[] MomBins;
        return true;
    }
    ScatLenR = 1. / fitSP6->GetParameter(1);
    EffRanR = fitSP6->GetParameter(2);
    DLM_FitHisto(hFitI, fitSP6I, "Q, S, N, R, M", "", &MinOpt);
    if (fitSP6I->GetChisquare() / fitSP4I->GetChisquare() > 0.8)
    {
        delete fitSP2I;
        delete fitSP6I;
        fitSPI = fitSP4I;
        delete[] MomBins;
        return true;
    }
    ScatLenI = 1. / fitSP6I->GetParameter(1);
    EffRanI = fitSP6I->GetParameter(2);

    ScatLen = complex<double>(ScatLenR, ScatLenI);
    EffRan = complex<double>(EffRanR, EffRanI);

    delete fitSP2;
    delete fitSP4;
    fitSPR = fitSP6;
    delete[] MomBins;
    return true;
}

bool PotentialDesignerEngine(char* BaseFileName){
  TString RootFileName = TString(BaseFileName)+".root";//output (WILL BE REWRITEN!!!!)
  TString TextFileName = TString(BaseFileName)+".fit";//input

  double Mass1 = 0;
  double Mass2 = 0;
  unsigned short usPW = 0;
  int Q1Q2 = 0;

  double kMin=  0;    //min momentum for the fit to the PS (0)
  double kMax = 100;    //max momentum for the fit to the PS (100)
  unsigned kBin = 50;    //num mommentum bins for the fit to the PS (50)
  unsigned nPar = 3;    //number of parameters for the effective range expansion when fitting the PS (3)
  double eps = 1e-8;     //the numerical precision parameter (1e-8)
  bool out = true;     //do we save the fit in the .root output file, allowed values are 0 and 1 (default is yes = 1)
  char* pot = new char [64];     //Gauss, DoubleGauss (default), Yukawa, YukawaDLM (modified to have some repulsive core to avoid singulartities)
  strcpy(pot,"DoubleGauss");
  char* pw = new char [4];      //in which partial wave should we place the potential. Possible values are s,p,d,f , s is the default
  strcpy(pw,"s");
  unsigned short NumPW;
  int coulomb = 0; //do we include the coulomb in the evaluation (give q1*q2). By default 0 (no).
  double par1 = 0;
  double par2 = 0.5;
  double par3 = 0;
  double par4 = 0.5;

  double rMin = 0;
  double rMax = 8;
  unsigned rBin = 1024;

  //this is for the output
  TH1F* hPhaseShifts = NULL;
  TH1F* hPotential = NULL;
  CATSparameters* pPars = NULL;
  CATSparameters cPars(CATSparameters::tSource,1,true);
  CATS Kitty;

  FILE* InFile;
  InFile = fopen(TextFileName.Data(), "r");
  if(InFile == nullptr){
    printf("\033[1;31mERROR:\033[0m Cannot open the file %s!\n",TextFileName.Data());
    delete [] pw;
    delete [] pot;
    return false;
  }

  char* ch_dummy_1 = new char [128];
  char* ch_dummy_2 = new char [128];

  while(!feof(InFile)){
    if(!fscanf(InFile, "%s %s", ch_dummy_1, ch_dummy_2)){
      printf("\033[1;33mWARNING (1)!\033[0m Possible bad input-file, error when reading from %s!\n",TextFileName.Data());
    }

    // Consume remaining characters on the line
    int chr;
    //std::cout << "XXX\n";
    do{
    chr = fgetc(InFile);
    //std::cout << chr;
    }
    while (chr != EOF && chr != '\n');
    //std::cout << "\n";

    printf("%s %s\n",ch_dummy_1, ch_dummy_2);

    //convert to only lower letters
    for (int ich = 0; ich < strlen(ch_dummy_1); ich++) {
        ch_dummy_1[ich] = tolower(ch_dummy_1[ich]);
    }
    for (int ich = 0; ich < strlen(ch_dummy_2); ich++) {
        ch_dummy_2[ich] = tolower(ch_dummy_2[ich]);
    }

    if(strcmp(ch_dummy_1,"m1")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic m1 value!\n",TextFileName.Data());
      }
      else{
        Mass1 = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"m2")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic m2 value!\n",TextFileName.Data());
      }
      else{
        Mass2 = atof(ch_dummy_2);
      }
    }
    /*
    if(strcmp(ch_dummy_1,"f_goal")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic scattering length!\n",TextFileName.Data());
      }
      else{
        f_goal = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"f_err")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic error of the scattering length!\n",TextFileName.Data());
      }
      else{
        f_err = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"d_goal")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic effective range!\n",TextFileName.Data());
      }
      else{
        d_goal = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"d_err")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic error of the error of the effective range!\n",TextFileName.Data());
      }
      else{
        d_err = atof(ch_dummy_2);
      }
    }
    */
    if(strcmp(ch_dummy_1,"kmin")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for kMin!\n",TextFileName.Data());
      }
      else{
        kMin = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"kmax")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for kMin!\n",TextFileName.Data());
      }
      else{
        kMax = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"kbin")==0){
      if(!isInteger(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for kBin!\n",TextFileName.Data());
      }
      else{
        kBin = atoi(ch_dummy_2);
      }
    }
    //if(strcmp(ch_dummy_1,"npar")==0){
    //  if(!isInteger(ch_dummy_2)){
    //    printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for nPar!\n",TextFileName.Data());
    //  }
    //  else{
    //    nPar = atoi(ch_dummy_2);
    //  }
    //}
    if(strcmp(ch_dummy_1,"eps")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for eps!\n",TextFileName.Data());
      }
      else{
        eps = atof(ch_dummy_2);
      }
    }
    //if(strcmp(ch_dummy_1,"out")==0){
    //  if(!isInteger(ch_dummy_2)){
    //    printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for 'out'!\n",TextFileName.Data());
    //  }
    //  else{
    //    out = bool(atoi(ch_dummy_2));
    //  }
    //}
    if(strcmp(ch_dummy_1,"pot")==0){
      if(strcmp(ch_dummy_2,"gauss")&&strcmp(ch_dummy_2,"doublegauss")&&strcmp(ch_dummy_2,"yukawa")&&strcmp(ch_dummy_2,"yukawadlm")&&strcmp(ch_dummy_2,"usmanicore")){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set an unknown potential %s!\n",TextFileName.Data(),ch_dummy_2);
      }
      else{
        strcpy(pot,ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"pw")==0){
      if(strcmp(ch_dummy_2,"s")&&strcmp(ch_dummy_2,"p")&&strcmp(ch_dummy_2,"d")&&strcmp(ch_dummy_2,"f")){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set an unknown partial wave %s!\n",TextFileName.Data(),ch_dummy_2);
      }
      else{
        strcpy(pw,ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"coulomb")==0){
      if(!isInteger(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for Q1*Q2 (Coulomb)!\n",TextFileName.Data());
      }
      else{
        coulomb = atoi(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"v_par1")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par1!\n",TextFileName.Data());
      }
      else{
        par1 = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"v_par2")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par2!\n",TextFileName.Data());
      }
      else{
        par2 = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"v_par3")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par3!\n",TextFileName.Data());
      }
      else{
        par3 = atof(ch_dummy_2);
      }
    }
    if(strcmp(ch_dummy_1,"v_par4")==0){
      if(!isNumber(ch_dummy_2)){
        printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par4!\n",TextFileName.Data());
      }
      else{
        par4 = atof(ch_dummy_2);
      }
    }
  }
  fclose(InFile);

  if(Mass1<=0 || Mass2<=0){
    printf("\033[1;31mERROR:\033[0m The mass had to be positive and non-zero!\n");
    delete [] ch_dummy_1;
    delete [] ch_dummy_2;
    delete [] pw;
    delete [] pot;
    return false;
  }
  if(strcmp(pw,"s")==0) NumPW = 1;
  if(strcmp(pw,"p")==0) NumPW = 2;
  if(strcmp(pw,"d")==0) NumPW = 3;
  if(strcmp(pw,"f")==0) NumPW = 4;

  //if(f_goal && f_err==0){
  //  f_err = 0.01*fabs(f_goal);
  //}
  //if(d_goal && d_err==0){
  //  d_err = 0.01*fabs(d_goal);
  //}
  //if(f_goal && d_goal==0 && d_err==0){
  //  d_err = 0.05;
  //}

  Kitty.SetMomBins(kBin,kMin,kMax);
  Kitty.SetThetaDependentSource(false);
  cPars.SetParameter(0,1.0);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetAutoNormSource(false);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);

  Kitty.SetExcludeFailedBins(false);
  Kitty.SetQ1Q2(coulomb);
  Kitty.SetQuantumStatistics(false);
  Kitty.SetRedMass( (Mass1*Mass2)/(Mass1+Mass2) );
  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,NumPW);
  Kitty.SetSpin(0,0);
  Kitty.SetChannelWeight(0, 1);

  if(strcmp(pot,"gauss")==0){
    pPars = new CATSparameters(CATSparameters::tPotential,2,true);
    Kitty.SetShortRangePotential(0,NumPW-1,SingleGauss,*pPars);
    Kitty.SetShortRangePotential(0,NumPW-1,0,par1);
    Kitty.SetShortRangePotential(0,NumPW-1,1,par2);
  }
  if(strcmp(pot,"doublegauss")==0){
    pPars = new CATSparameters(CATSparameters::tPotential,4,true);
    Kitty.SetShortRangePotential(0,NumPW-1,DoubleGaussSum,*pPars);
    Kitty.SetShortRangePotential(0,NumPW-1,0,par1);
    Kitty.SetShortRangePotential(0,NumPW-1,1,par2);
    Kitty.SetShortRangePotential(0,NumPW-1,2,par3);
    Kitty.SetShortRangePotential(0,NumPW-1,3,par4);
    //printf("p1 = %f\n",par1);
    //printf("p2 = %f\n",par2);
    //printf("p3 = %f\n",par3);
    //printf("p4 = %f\n",par4);
    //printf("kBin = %u\n",kBin);
  }
  if(strcmp(pot,"yukawa")==0){
    pPars = new CATSparameters(CATSparameters::tPotential,2,true);
    Kitty.SetShortRangePotential(0,NumPW-1,Yukawa,*pPars);
    Kitty.SetShortRangePotential(0,NumPW-1,0,par1);
    Kitty.SetShortRangePotential(0,NumPW-1,1,par2);
  }
  if(strcmp(pot,"yukawadlm")==0){
    pPars = new CATSparameters(CATSparameters::tPotential,2,true);
    Kitty.SetShortRangePotential(0,NumPW-1,YukawaDimiSmooth,*pPars);
    Kitty.SetShortRangePotential(0,NumPW-1,0,par1);
    Kitty.SetShortRangePotential(0,NumPW-1,1,par2);
  }
  if(strcmp(pot,"usmanicore")==0){
    pPars = new CATSparameters(CATSparameters::tPotential,4,true);
    Kitty.SetShortRangePotential(0,NumPW-1,UsmaniFit,*pPars);
    Kitty.SetShortRangePotential(0,NumPW-1,0,par1);
    Kitty.SetShortRangePotential(0,NumPW-1,1,par2);
    Kitty.SetShortRangePotential(0,NumPW-1,2,par3);
    Kitty.SetShortRangePotential(0,NumPW-1,3,par4);
  }

  Kitty.SetEpsilonConv(eps);
  Kitty.SetEpsilonProp(eps);

  Kitty.SetNotifications(CATS::nWarning);
  Kitty.KillTheCat();

  hPhaseShifts = new TH1F("hPhaseShifts","hPhaseShifts",kBin,kMin,kMax);
  for(unsigned uBin=0; uBin<kBin; uBin++){
    hPhaseShifts->SetBinContent(uBin+1,Kitty.GetPhaseShift(uBin,0,NumPW-1));
    if(Kitty.GetPhaseShift(uBin,0,NumPW-1)>TMath::Pi()*0.5){
      hPhaseShifts->SetBinContent(uBin+1,Kitty.GetPhaseShift(uBin,0,NumPW-1)-TMath::Pi());
    }
  }

  hPotential = new TH1F("hPotential","hPotential",rBin,rMin,rMax);
  for(unsigned uRad=0; uRad<rBin; uRad++){
    double RAD = hPotential->GetBinCenter(uRad+1);
    hPotential->SetBinContent(uRad+1,Kitty.EvaluateThePotential(0,NumPW-1,0.5*(kMin+kMax),RAD));
  }


  TFile fOutput(RootFileName,"recreate");
  hPhaseShifts->Write();
  hPotential->Write();

  delete [] ch_dummy_1;
  delete [] ch_dummy_2;
  if(hPhaseShifts){delete hPhaseShifts; hPhaseShifts=NULL;}
  if(hPotential){delete hPotential; hPotential=NULL;}
  if(pPars){delete pPars; pPars=NULL;}
  return true;
}

bool PotentialDesignerImEngine(char *BaseFileName)
{
    TString RootFileName = TString(BaseFileName) + ".root"; // output (WILL BE REWRITEN!!!!)
    TString TextFileName = TString(BaseFileName) + ".fit";  // input

    double Mass1 = 0;
    double Mass2 = 0;
    unsigned short usPW = 0;
    int Q1Q2 = 0;

    double kMin = 0;          // min momentum for the fit to the PS (0)
    double kMax = 100;        // max momentum for the fit to the PS (100)
    unsigned kBin = 50;       // num mommentum bins for the fit to the PS (50)
    unsigned nPar = 3;        // number of parameters for the effective range expansion when fitting the PS (3)
    double eps = 1e-8;        // the numerical precision parameter (1e-8)
    bool out = true;          // do we save the fit in the .root output file, allowed values are 0 and 1 (default is yes = 1)
    char *pot = new char[64]; // Gauss, DoubleGauss (default), Yukawa, YukawaDLM (modified to have some repulsive core to avoid singulartities)
    strcpy(pot, "ComplexGaussian");
    char *pw = new char[4]; // in which partial wave should we place the potential. Possible values are s,p,d,f , s is the default
    strcpy(pw, "s");
    unsigned short NumPW;
    int coulomb = 0; // do we include the coulomb in the evaluation (give q1*q2). By default 0 (no).
    double par1 = 0;
    double par2 = 0.5;
    double par3 = 0;
    double par4 = 0.5;
    double par5 = 0;
    double par6 = 0.5;

    double rMin = 0;
    double rMax = 8;
    unsigned rBin = 1024;

    // this is for the output
    TH1F *hPhaseShifts = NULL;
    TH1F *hScattPars = NULL;
    TH1F *hPotentialReal = NULL;
    TH1F *hPotentialImag = NULL;

    CATSparameters *pPars = NULL;
    CATSparameters cPars(CATSparameters::tSource, 1, true);
    CATS Kitty;

    FILE *InFile;
    InFile = fopen(TextFileName.Data(), "r");
    if (InFile == nullptr)
    {
        printf("\033[1;31mERROR:\033[0m Cannot open the file %s!\n", TextFileName.Data());
        delete[] pw;
        delete[] pot;
        return false;
    }

    char *ch_dummy_1 = new char[128];
    char *ch_dummy_2 = new char[128];

    while (!feof(InFile))
    {
        if (!fscanf(InFile, "%s %s", ch_dummy_1, ch_dummy_2))
        {
            printf("\033[1;33mWARNING (1)!\033[0m Possible bad input-file, error when reading from %s!\n", TextFileName.Data());
        }

        // Consume remaining characters on the line
        int chr;
        // std::cout << "XXX\n";
        do
        {
            chr = fgetc(InFile);
            // std::cout << chr;
        } while (chr != EOF && chr != '\n');
        // std::cout << "\n";

        printf("%s %s\n", ch_dummy_1, ch_dummy_2);

        // convert to only lower letters
        for (int ich = 0; ich < strlen(ch_dummy_1); ich++)
        {
            ch_dummy_1[ich] = tolower(ch_dummy_1[ich]);
        }
        for (int ich = 0; ich < strlen(ch_dummy_2); ich++)
        {
            ch_dummy_2[ich] = tolower(ch_dummy_2[ich]);
        }

        if (strcmp(ch_dummy_1, "m1") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic m1 value!\n", TextFileName.Data());
            }
            else
            {
                Mass1 = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "m2") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic m2 value!\n", TextFileName.Data());
            }
            else
            {
                Mass2 = atof(ch_dummy_2);
            }
        }
        /*
        if(strcmp(ch_dummy_1,"f_goal")==0){
          if(!isNumber(ch_dummy_2)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic scattering length!\n",TextFileName.Data());
          }
          else{
            f_goal = atof(ch_dummy_2);
          }
        }
        if(strcmp(ch_dummy_1,"f_err")==0){
          if(!isNumber(ch_dummy_2)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic error of the scattering length!\n",TextFileName.Data());
          }
          else{
            f_err = atof(ch_dummy_2);
          }
        }
        if(strcmp(ch_dummy_1,"d_goal")==0){
          if(!isNumber(ch_dummy_2)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic effective range!\n",TextFileName.Data());
          }
          else{
            d_goal = atof(ch_dummy_2);
          }
        }
        if(strcmp(ch_dummy_1,"d_err")==0){
          if(!isNumber(ch_dummy_2)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic error of the error of the effective range!\n",TextFileName.Data());
          }
          else{
            d_err = atof(ch_dummy_2);
          }
        }
        */
        if (strcmp(ch_dummy_1, "kmin") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for kMin!\n", TextFileName.Data());
            }
            else
            {
                kMin = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "kmax") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for kMin!\n", TextFileName.Data());
            }
            else
            {
                kMax = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "kbin") == 0)
        {
            if (!isInteger(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for kBin!\n", TextFileName.Data());
            }
            else
            {
                kBin = atoi(ch_dummy_2);
            }
        }
        // if(strcmp(ch_dummy_1,"npar")==0){
        //   if(!isInteger(ch_dummy_2)){
        //     printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for nPar!\n",TextFileName.Data());
        //   }
        //   else{
        //     nPar = atoi(ch_dummy_2);
        //   }
        // }
        if (strcmp(ch_dummy_1, "eps") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numberic value for eps!\n", TextFileName.Data());
            }
            else
            {
                eps = atof(ch_dummy_2);
            }
        }
        // if(strcmp(ch_dummy_1,"out")==0){
        //   if(!isInteger(ch_dummy_2)){
        //     printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for 'out'!\n",TextFileName.Data());
        //   }
        //   else{
        //     out = bool(atoi(ch_dummy_2));
        //   }
        // }
        if (strcmp(ch_dummy_1, "pot") == 0)
        {
            if (strcmp(ch_dummy_2, "complexgaussian"))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set an unknown potential %s!\n", TextFileName.Data(), ch_dummy_2);
            }
            else
            {
                strcpy(pot, ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "pw") == 0)
        {
            if (strcmp(ch_dummy_2, "s") && strcmp(ch_dummy_2, "p") && strcmp(ch_dummy_2, "d") && strcmp(ch_dummy_2, "f"))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set an unknown partial wave %s!\n", TextFileName.Data(), ch_dummy_2);
            }
            else
            {
                strcpy(pw, ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "coulomb") == 0)
        {
            if (!isInteger(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-integer value for Q1*Q2 (Coulomb)!\n", TextFileName.Data());
            }
            else
            {
                coulomb = atoi(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "v_par1") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par1!\n", TextFileName.Data());
            }
            else
            {
                par1 = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "v_par2") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par2!\n", TextFileName.Data());
            }
            else
            {
                par2 = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "v_par3") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par3!\n", TextFileName.Data());
            }
            else
            {
                par3 = atof(ch_dummy_2);
                
            }
        }
        if (strcmp(ch_dummy_1, "v_par4") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par4!\n", TextFileName.Data());
            }
            else
            {
                par4 = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "v_par5") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par5!\n", TextFileName.Data());
            }
            else
            {
                par5 = atof(ch_dummy_2);
            }
        }
        if (strcmp(ch_dummy_1, "v_par6") == 0)
        {
            if (!isNumber(ch_dummy_2))
            {
                printf("\033[1;33mWARNING!\033[0m Possible bad input-file (%s), trying to set a non-numerical value for par6!\n", TextFileName.Data());
            }
            else
            {
                par6 = atof(ch_dummy_2);
            }
        }
    }
    fclose(InFile);

    if (Mass1 <= 0 || Mass2 <= 0)
    {
        printf("\033[1;31mERROR:\033[0m The mass had to be positive and non-zero!\n");
        delete[] ch_dummy_1;
        delete[] ch_dummy_2;
        delete[] pw;
        delete[] pot;
        return false;
    }
    if (strcmp(pw, "s") == 0)
        NumPW = 1;
    if (strcmp(pw, "p") == 0)
        NumPW = 2;
    if (strcmp(pw, "d") == 0)
        NumPW = 3;
    if (strcmp(pw, "f") == 0)
        NumPW = 4;

    // if(f_goal && f_err==0){
    //   f_err = 0.01*fabs(f_goal);
    // }
    // if(d_goal && d_err==0){
    //   d_err = 0.01*fabs(d_goal);
    // }
    // if(f_goal && d_goal==0 && d_err==0){
    //   d_err = 0.05;
    // }

    Kitty.SetMomBins(kBin, kMin, kMax);
    Kitty.SetThetaDependentSource(false);
    cPars.SetParameter(0, 1.0);
    Kitty.SetAnaSource(GaussSource, cPars);
    Kitty.SetAutoNormSource(false);
    Kitty.SetUseAnalyticSource(true);
    Kitty.SetMomentumDependentSource(false);

    Kitty.SetExcludeFailedBins(false);
    Kitty.SetQ1Q2(coulomb);
    Kitty.SetQuantumStatistics(false);
    Kitty.SetRedMass((Mass1 * Mass2) / (Mass1 + Mass2));
    Kitty.SetNumChannels(1);
    Kitty.SetNumPW(0, NumPW);
    Kitty.SetChannelWeight(0, 1);
    Kitty.SetUsingiCATS(true);
    if (strcmp(pot, "complexgaussian") == 0)
    {
        pPars = new CATSparameters(CATSparameters::tPotential, 3, true);
        pPars->SetParameter(0, par1);
        pPars->SetParameter(1, par2);
        pPars->SetParameter(2, par3);/// MUST BE SETTED otherwise BUG!!
        Kitty.SetShortRangePotential(0, NumPW - 1, ComplexGaussian, *pPars);
        Kitty.SetShortRangePotential(0, NumPW - 1, 0, par1);
        Kitty.SetShortRangePotential(0, NumPW - 1, 1, par2);
        Kitty.SetShortRangePotential(0, NumPW - 1, 2, par3);
        cout << "Inside Engine" <<endl;
        printf("p1 = %f\n",par1);
        printf("p2 = %f\n",par2);
        printf("p3 = %f\n",par3);
    }
    Kitty.SetEpsilonConv(eps);
    Kitty.SetEpsilonProp(eps);

    Kitty.SetNotifications(CATS::nWarning);
    Kitty.KillTheCat();

    array<complex<double>, 2> ScattPars = Kitty.EvalComplexScatPars(0, NumPW - 1);
    complex<double> ScattLen = ScattPars[0];
    complex<double> EffRan = ScattPars[1];

    cout << "Inside Engine getting Scatt Pars" << endl;
    cout << "ScattLen" << ScattLen << endl;
    cout << "EffRan" << EffRan << endl;

    hScattPars = new TH1F("hScattPars", "hScattPars", 4, 0.5, 4.5);
    hScattPars->SetBinContent(1, real(ScattPars[0]));
    hScattPars->SetBinContent(2, imag(ScattPars[0]));
    hScattPars->SetBinContent(3, real(ScattPars[1]));
    hScattPars->SetBinContent(4, imag(ScattPars[1]));

    hPotentialReal = new TH1F("hPotentialReal", "hPotentialReal", rBin, rMin, rMax);
    hPotentialImag = new TH1F("hPotentialImag", "hPotentialImag", rBin, rMin, rMax);
    for (unsigned uRad = 0; uRad < rBin; uRad++)
    {
        double RAD = hPotentialReal->GetBinCenter(uRad + 1);
        hPotentialReal->SetBinContent(uRad + 1, real(Kitty.EvaluateTheComplexPotential(0, NumPW - 1, 0.5 * (kMin + kMax), RAD)));
        hPotentialImag->SetBinContent(uRad + 1, imag(Kitty.EvaluateTheComplexPotential(0, NumPW - 1, 0.5 * (kMin + kMax), RAD)));
    }

    TFile fOutput(RootFileName, "recreate");
    hScattPars->Write();
    hPotentialReal->Write();
    hPotentialImag->Write();

    delete[] ch_dummy_1;
    delete[] ch_dummy_2;
    if (hPhaseShifts)
    {
        delete hPhaseShifts;
        hPhaseShifts = NULL;
    }
    if (hScattPars)
    {
        delete hScattPars;
        hScattPars = NULL;
    }
    if (hPotentialReal)
    {
        delete hPotentialReal;
        hPotentialReal = NULL;
    }
    if (hPotentialImag)
    {
        delete hPotentialImag;
        hPotentialImag = NULL;
    }
    if (pPars)
    {
        delete pPars;
        pPars = NULL;
    }
    return true;
}

TF1* fit_source_kdp(TH1F* hSrc, KdpPars& SrcPar, double& Chi2, double Chi2_Limit, double KdpFitMax, int KdpMode){

  double lowerlimit, upperlimit;
  GetCentralInterval(*hSrc, 0.98, lowerlimit, upperlimit, true);
  //if(lowerlimit>5) lowerlimit = 5;
  lowerlimit = 0;
  if(upperlimit>KdpFitMax+2) upperlimit = KdpFitMax+2;


  //const double Chi2_Limit = 3;
  //const double KdpFitMax = 8;

  TF1* fSrc;
  SetUpKdpPars(fSrc,KdpMode,KdpFitMax);

  //for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
  //  fSrc->SetParLimits(1+uP*3,0.1,fSrc->GetParameter(0+uP*3));
  //}
  //fSrc->FixParameter(0,1);
  //fSrc->FixParameter(1,0.5);
  //fSrc->FixParameter(2,0.5);

  //fSrc->FixParameter(3,6);
  //fSrc->FixParameter(4,2);
  //fSrc->FixParameter(5,1);

  hSrc->Fit(fSrc,"Q, S, N, R, M","",lowerlimit,upperlimit);
  //TFile fWTF("fWTF.root","recreate");
  //hSrc->Write();
  //fSrc->Write();


  double Integral = fSrc->Integral(0,KdpFitMax*4);
  if(fabs(Integral-1)>1e-2){
    printf("SUPER BIG BUG WITH THE KDP (%f)!!!\n",Integral);
  }

  Chi2 = 0;
  double NDPts_chi2 = 0;

  //up to 8 fm
  for(unsigned uRad=0; uRad<hSrc->GetNbinsX(); uRad++){
    double Rad = hSrc->GetBinCenter(uRad+1);
    if(Rad>KdpFitMax) break;
    double dst = hSrc->GetBinContent(uRad+1)-fSrc->Eval(Rad);
    double err = hSrc->GetBinError(uRad+1);
    if(hSrc->GetBinContent(uRad+1)){
      Chi2 += (dst*dst)/(err*err);
      NDPts_chi2++;
    }
  }
  Chi2 /= NDPts_chi2;

  if(Chi2>Chi2_Limit){
    printf("WARNING: BadFitWarning (Chi2/ndf = %.2f)\n",Chi2);
  }

  for(unsigned uP=0; uP<KdpPars::NumDistos; uP++){
    SrcPar.mean[uP] = fSrc->GetParameter(0+uP*3);
    SrcPar.stdv[uP] = fSrc->GetParameter(1+uP*3);
    if(uP!=KdpPars::NumDistos-1)
      SrcPar.wght[uP] = fSrc->GetParameter(2+uP*3);
  }

  return fSrc;
}

//takes a 3D histo of Mt Kstar Rstar and returns a 2D kdp histo of Mt Kstar
DLM_Histo<KdpPars>* Convert_3Dsource_Kdp(DLM_Histo<float>& dlmMtKstarRstar, const bool CecaStyle, double Chi2_Limit, double KdpFitMax, int KdpMode){
TString BaseFileName = TString::Format("/home/dimihayl/Software/LocalFemto/Output/HighMtProblem/CECA_kstar_vs_integrated/");
TFile fTest(BaseFileName + "fTEST.root","recreate");
    DLM_Histo<KdpPars>* KdpResult = NULL;

    if(dlmMtKstarRstar.GetDim()!=3){
        printf("\033[1;31mERROR:\033[0m Convert_3Dsource_Kdp has to be a 3D histo\n");
        return KdpResult;
    }

    const unsigned NumMtBins = dlmMtKstarRstar.GetNbins(CecaStyle?2:0);
    double* MtBinRange = dlmMtKstarRstar.GetBinRange(CecaStyle?2:0);
    const unsigned NumKstarBins = dlmMtKstarRstar.GetNbins(CecaStyle?0:1);
    double* KstarBinRange = dlmMtKstarRstar.GetBinRange(CecaStyle?0:1);
    const unsigned NumRadBins = dlmMtKstarRstar.GetNbins(CecaStyle?1:2);
    double* RadBinRange = dlmMtKstarRstar.GetBinRange(CecaStyle?1:2);
    const double rMin = dlmMtKstarRstar.GetLowEdge(CecaStyle?1:2);
    const double rMax = dlmMtKstarRstar.GetLowEdge(CecaStyle?1:2);
//printf("%i %i %i\n",NumKstarBins,NumRadBins,NumMtBins);
    KdpResult = new DLM_Histo<KdpPars> ();
    KdpResult->SetUp(2);
    KdpResult->SetUp(0,NumMtBins,MtBinRange);
    KdpResult->SetUp(1,NumKstarBins,KstarBinRange);
    KdpResult->Initialize();

    for(unsigned uMt=0; uMt<NumMtBins; uMt++){
        for(unsigned uKstar=0; uKstar<NumKstarBins; uKstar++){
            //printf("cc %u %u\n",uMt,uKstar);
            TH1F* hToFit = new TH1F("hToFit_c3dskdp","hToFit_c3dskdp",NumRadBins,RadBinRange);
            for(unsigned uRad=0; uRad<NumRadBins; uRad++){
                double src_val = CecaStyle?dlmMtKstarRstar.GetBinContent(uKstar,uRad,uMt):dlmMtKstarRstar.GetBinContent(uMt,uKstar,uRad);
                double src_err;
                if(src_val==0){src_err = 1;}
                else{
                    src_err = CecaStyle?dlmMtKstarRstar.GetBinError(uKstar,uRad,uMt):dlmMtKstarRstar.GetBinError(uMt,uKstar,uRad);
                }
                //printf("src_val = %.4f +/- %.4f\n",src_val,src_err);
                hToFit->SetBinContent(uRad+1, src_val);
                hToFit->SetBinError(uRad+1, src_err);
            }
            hToFit->Scale(1./hToFit->Integral(),"width");
            double Chi2=0;
            KdpPars my_kdp;
            TF1* fit_ptr = fit_source_kdp(hToFit, my_kdp, Chi2, Chi2_Limit, KdpFitMax, KdpMode);
            fTest.cd();
            //if(Chi2>3){
                static int counter = 0;

                hToFit->SetName(TString::Format("hToFit_%u_%.0f",uMt,dlmMtKstarRstar.GetBinCenter(CecaStyle?0:1,uKstar)));
                fit_ptr->SetName(TString::Format("fit_ptr_%u_%.0f",uMt,dlmMtKstarRstar.GetBinCenter(CecaStyle?0:1,uKstar)));
                
                hToFit->Write();
                fit_ptr->Write();

                counter++;
            //}

            KdpResult->SetBinContent(uMt, uKstar, my_kdp);
            delete hToFit;
            delete fit_ptr;
        }
    }

    delete [] MtBinRange;
    delete [] KstarBinRange;
    delete [] RadBinRange;
    return KdpResult;
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
