

#include "Ana_pi_d_RUN2.h"

#include "CATS.h"
#include "CATStools.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
#include "CommonAnaFunctions.h"
#include "CATSconstants.h"
#include "DLM_Histo.h"
#include "DLM_RootWrapper.h"
#include "DLM_MathFunctions.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomposition.h"
#include "DLM_HistoAnalysis.h"

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGenPhaseSpace.h"

//path to the analysis code input folder
TString INPUT_FOLDER;
void SetAnalysisFolderPath(char* folder_path){
    INPUT_FOLDER = folder_path;
}

TString OUTPUT_FOLDER = "./";
void SetOutputFolderPath(char* folder_path){
    OUTPUT_FOLDER = folder_path;
}

void clean_graph(TGraph* gInput){
  int num_pts = gInput->GetN();
  for(int iPts=num_pts-1; iPts>=0; iPts--){
    gInput->RemovePoint(iPts);
  }
}

void pi_d_source(){

  double reff_min=1000; double reff_max=0; double reff_def;

  unsigned NumRadBins = 1024;
  double rMin = 0;
  double rMax = 12;
  TH1F* hSourceUpper = new TH1F("hSourceUpper","hSourceUpper",NumRadBins,rMin,rMax);
  TH1F* hSourceLower = new TH1F("hSourceLower","hSourceLower",NumRadBins,rMin,rMax);
  TH1F* hSource = new TH1F("hSource","hSource",NumRadBins,rMin,rMax);

  for(unsigned uRad=0; uRad<NumRadBins; uRad++){
    hSourceUpper->SetBinContent(uRad+1,0);
    hSourceLower->SetBinContent(uRad+1,1e6);
  }

  std::vector<float> core_radii_mTint = {1.036,1.075,1.112};
  std::vector<float> core_radii_mT0 = {1.163,1.216,1.271};
  std::vector<float> core_radii_mT1 = {1.068,1.109,1.150};
  std::vector<float> core_radii_mT2 = {0.996,1.033,1.068};
  std::vector<float> core_radii_mT3 = {0.932,0.968,1.003};
  std::vector<float> core_radii_mT4 = {0.782,0.831,0.880};
  std::vector<float> vectors[] = {core_radii_mTint, core_radii_mT0, core_radii_mT1, core_radii_mT2, core_radii_mT3, core_radii_mT4};

  for(int uA=0; uA<1; uA++){for(int uB=0; uB<1; uB++){for(int uE=0; uE<3; uE++){
        //#for the sake of simplicity, make it 1.074 +/- 0.39 and min/avg/max = 1.035/1.074/1.113
        //std::vector<float> core_radii = {1.035,1.074,1.113};//mt integrated
        std::vector<float> core_radii = {0.782,0.831,0.880};

        DLM_CommonAnaFunctions AnalysisObject;
        AnalysisObject.SetCatsFilesFolder(INPUT_FOLDER.Data());

        unsigned NumMomBins = 60;
        double kMin = 0;
        double kMax = 300;
        CATS Kitty;
        Kitty.SetMomBins(NumMomBins, kMin, kMax);
        AnalysisObject.SetUpCats_pi_d(Kitty, "", "McGauss_ResoTM", 0, uA*10000+uB*1000+200+uE);
        DLM_CleverMcLevyResoTM* MagicSource = AnalysisObject.GetCleverMcLevyResoTM_pi_d();
        printf("%i %i %i\n",uA,uB,uE);
        
        for(unsigned uSrc=0; uSrc<core_radii.size(); uSrc++){
          if(uSrc%3==0){
            printf("mT%i-----------------\n",int(uSrc)/3-1);
          }
          double r_eff = GetReff(*MagicSource, core_radii.at(uSrc));
          if(uA==0&&uB==0&&uE==0&&uSrc==1) reff_def = r_eff;
          printf("rcore = %.3f; reff = %.3f\n", core_radii.at(uSrc), r_eff);
          if(reff_min>r_eff) reff_min = r_eff;
          if(reff_max<r_eff) reff_max = r_eff;
        }

        for(unsigned uSrc=0; uSrc<core_radii.size(); uSrc++){
          Kitty.SetAnaSource(0, core_radii.at(uSrc));
          for(unsigned uRad=0; uRad<NumRadBins; uRad++){
            double rstar = hSource->GetBinCenter(uRad+1);
            double sr = Kitty.EvaluateTheSource(10, rstar, 0);
            double sr_lower = hSourceLower->GetBinContent(uRad+1);
            double sr_upper = hSourceLower->GetBinContent(uRad+1);
            if(uA==0&&uB==0&&uE==0&&uSrc==1){
              hSource->SetBinContent(uRad+1, sr);
            }
            if(sr>sr_upper){
              hSourceUpper->SetBinContent(uRad+1,sr);
            }
            if(sr<sr_lower){
              hSourceLower->SetBinContent(uRad+1,sr);
            }
          }
        }
      }
    }
  }

  TF1* fEffSrc = new TF1("fEffSrc",GaussSourceTF1,rMin,rMax,1);
  fEffSrc->SetParameter(0,reff_def);

  TFile fOutput(TString::Format("%s/Source/pi_d_source.root",INPUT_FOLDER.Data()),"recreate");
  hSource->Write();
  hSourceLower->Write();
  hSourceUpper->Write();
  fEffSrc->Write();

  printf("reff_min = %.3f; reff_max = %.3f\n", reff_min, reff_max);
  printf("=> %.3f +/- %.3f\n", (reff_max+reff_min)*0.5, (reff_max-reff_min)*0.5);

  delete hSource;
  delete hSourceLower;
  delete hSourceUpper;

}



//par[0] = source size
//par[1] = alpha par
//par[2] = lambda_genuine
//par[3] = fraction_d_from_Delta = frac_D
//par[4] = delta_amplitude
//par[5] = CkCutOff
//par[6] = CkConv
//par[7] = treat it as a flag for the fit. 0 = we do NOT multiply C_Delta with C_femto, and 1 where we do so
//par[8] = Delta normalization, should be a dummy fixed to -1e6, giving instructions to the code to renorm
//par[9] = mass of daughter 1 (say the pion)
//par[10] = mass of daughter 2 (say the proton)
//par[11] = mass of the delta
//par[12] = width of the delta
//par[13] = avg pT of the daughters
//par[14] = effective temperature
//par[15] = just in case (placeholder)
//par[16] = norm
//par[17] = 0 to remove the linear term of the pol3
//par[18] = position of the max of the pol3
//par[19] = p3 parameter of the pol3
//par[20] = -1e6 to switch off the pol4
//NEW in JULY 2025:
//par[21] = amplitude of extra sill reso
//par[22] = mass of daughter 1 (d)
//par[23] = mass of daughter 1 (pi)
//par[24] = mass of an extra sill reso
//par[25] = width of an extra sill reso
//the total fit function is:
//Baseline*(lambda_femto*(1-frac_d)*Ck_femto + norm_delta*Ck_delta*(either 1 or Ck_femto) + lamda_flat)
DLM_CkDecomposition* ALICE_CkDec=NULL;
TGraph* Kstar_Modifier = NULL;
double ALICE_fit(double* kstar, double* par){
  if(!ALICE_CkDec) return 0;
  double& MOM = kstar[0];
  
  ALICE_CkDec->GetCk()->SetSourcePar(0, par[0]);
  ALICE_CkDec->GetCk()->SetCutOff(par[5],par[6]);
  ALICE_CkDec->Update();

  double Baseline = DLM_Baseline(kstar, &par[16]);
  double modified_kstar;
  modified_kstar = kstar[0]/Kstar_Modifier->Eval(kstar[0]);
  
  double Delta;
  Delta = SillBoltzmann_kstar(&modified_kstar, &par[8]);
  double Femto = ALICE_CkDec->EvalCk(MOM);
  //N.B. the Delta goes to zero, hence does not count as a proper lambda par!
  if(TMath::Nint(par[7])==0){
    return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta + (1.-par[2]));
  }
  else{
    return Baseline*(par[2]*(1.-par[3])*Femto + par[4]*Delta*Femto + (1.-par[2]));
  }
}

//if the |mT_bin|>100 => we deal with pi^- d
//if its >0, we deal with pi^+ d
//-1 = pip mt int
//-101 = pim mt int
void MainAnalysisFunction(char* Description, std::vector<int> fit_types, int mt_bin, int NumIter, bool Bootstrap, bool DataVar, bool FitVar, int SEED){
  if(mt_bin%100<-1 || mt_bin%100>5){
    printf("Bad mT bin\n");
    abort();
  }
  //these are the first file you got from Bhawani. On 29th Oct 2024 this was updated and now completely
  //converd by the if statements below, these paths here are for bookkeeping only and have no affect right now!
  TString InputDataFileName = "";
  TString DefaultHistoName = "histDefault";
  TString VarDirName = "Raw";
  TString VarHistoBaseName = "histVar_";
  TString MeListName1 = "PairDist";
  TString MeListName2 = "PairReweighted";
  TString MeHistName = "hTotalMEpairs _Rebinned_5_Reweighted";
  TString SeHistName = "";


  //variation we will resample
  int NumDataVars = 44;

  const double max_kstar = 700;
  const unsigned max_mom_bins = 35;
  std::vector<float> fit_range = {600, 500};
  std::vector<float> eff_source_radii;// = {1.51, 1.39, 1.63};
  std::vector<float> lambda_gen_relvar = {1.0, 0.95, 1.05};
  //the ratio of the lambda pars (yields) of pairs counting as primary
  //or decay products of the delta. If negative, the value is fitted, with 
  //starting value of fabs(1st elements) and limits fabs(2nd and 3rd element)
  std::vector<float> amplitude_delta = {-1e8, -1, -1e14};
  //these are the scenario B and A (B = 2/3 became default it seems, the other is systematics), how the pT of the p-pi relates to d-pi
  std::vector<float> pT_scale = {2./3., 3./4.};
  std::vector<TString> InputConversionFileName = 
    { TString::Format("%s/kstar_conversion/ResonanceOutput_B.root",INPUT_FOLDER.Data()), 
      TString::Format("%s/kstar_conversion/ResonanceOutput_A.root",INPUT_FOLDER.Data())};
  float CKCUTOFF = 380;
  const double CkConv = 700;

  //these are all per mt bin
  std::vector<std::vector<float>> mT_low;
  std::vector<std::vector<float>> mT_up;
  std::vector<float> mT_low_pip_v1 = {1030,1160,1280,1400,1550};
  std::vector<float> mT_up_pip_v1 = {1160,1280,1400,1550,2400};
  std::vector<float> mT_low_pim_v1 = {1030,1160,1240,1380,1520};
  std::vector<float> mT_up_pim_v1 = {1160,1240,1380,1520,2400};

  mT_low.push_back(mT_low_pip_v1);
  mT_low.push_back(mT_low_pim_v1);

  mT_up.push_back(mT_up_pip_v1);
  mT_up.push_back(mT_up_pim_v1);

  double avg_mt = 1270;
  if(mt_bin>=0 && mt_bin<1000){
    avg_mt = 0.5*(mT_up.at(mt_bin/100).at(mt_bin%100)+mT_low.at(mt_bin/100).at(mt_bin%100));
  }
  double avg_mass_pid = 0.5*(Mass_d+Mass_pic);
  double kay_tee = sqrt(avg_mt*avg_mt-avg_mass_pid*avg_mass_pid);
  
  std::vector<float> lambda_gen_avg = {0.850*0.88, 0.874*0.88, 0.888*0.88, 0.897*0.88, 0.918*0.88};
  double lambda_gen_mTint = 0.88*0.88;//avg_mt of 1.27

  //take 0.75 or 0 weight for the Delta++, based on the system (pip or pim)
  //N.B. for pim we should actually have 1/2 D0 and 1/2 D-, but we dont have access to it
  double Dpp_weight = 0.75;
  if(fabs(mt_bin)>=100 && mt_bin<1000) Dpp_weight = 0;
  std::vector<float> pip_gamma_D0 = {74.11, 74.06, 74.02, 70.19, 90.29};
  std::vector<float> pip_temp_D0 = {18.51, 18.27, 18.04, 16.49, 17.41};
  std::vector<float> pip_gamma_Dpp = {109.63, 106.28, 103.21, 99.56, 86.01};
  std::vector<float> pip_temp_Dpp = {26.49, 26.71, 26.92, 26.30, 23.25};
  //these values are avaraged over Delta++ and 0, with weight of 0.75 for the ++
  std::vector<float> pip_gamma = {100.75, 98.23, 95.92, 92.22, 87.53};
  std::vector<float> pip_temp = {24.49, 24.60, 24.70, 23.84, 21.79};
  float pip_avg_gamma = 95;
  float pip_avg_temp = 24.5;

  std::vector<float> reff_avg_val = {1.673, 1.552, 1.464, 1.390, 1.226};
  std::vector<float> reff_avg_err = {0.132, 0.119, 0.115, 0.118, 0.141};
  double reff_mTint_val = 1.513;//avg_mt of 1.27
  double reff_mTint_err = 0.117;
  if(mt_bin>=0 && mt_bin<1000){
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100));
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100)-reff_avg_err.at(mt_bin%100));
    eff_source_radii.push_back(reff_avg_val.at(mt_bin%100)+reff_avg_err.at(mt_bin%100));
  }
  else{
    eff_source_radii.push_back(reff_mTint_val);
    eff_source_radii.push_back(reff_mTint_val-reff_mTint_err);
    eff_source_radii.push_back(reff_mTint_val+reff_mTint_err);
  }

  TH2F* hResolution_pd = NULL;
  TH2F* hTemp = NULL;

  TFile* fReso = new TFile(TString::Format("%s/ReconstructionFiles/AnalysisResults3996.root",INPUT_FOLDER.Data()), "read");
  TDirectoryFile* fDir_reso = NULL;
  TList* fList1_reso = NULL;
  TList* fList2_reso = NULL;
  TList* fList3_reso = NULL;

  fDir_reso = (TDirectoryFile*)(fReso->FindObjectAny("MBResultsQA0"));
  fDir_reso->GetObject("MBResultsQA0", fList1_reso);
  fList2_reso = (TList*)fList1_reso->FindObject("PairQA");
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle0_Particle2");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle0_Particle2");
  gROOT->cd();
  hResolution_pd = (TH2F*)hTemp->Clone("hResolution_pd");
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle1_Particle3");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle1_Particle3");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle0_Particle3");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle0_Particle3");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
  fList3_reso = (TList*)fList2_reso->FindObject("QA_Particle1_Particle2");
  hTemp = (TH2F*)fList3_reso->FindObject("MomentumResolutionSE_Particle1_Particle2");
  gROOT->cd();
  hResolution_pd->Add(hTemp);
  hResolution_pd->GetXaxis()->SetLimits(hTemp->GetXaxis()->GetXmin()*1000.,hTemp->GetXaxis()->GetXmax()*1000.);
  hResolution_pd->GetYaxis()->SetLimits(hTemp->GetYaxis()->GetXmin()*1000.,hTemp->GetYaxis()->GetXmax()*1000.);
  gROOT->cd();
  
  delete fReso;

  TRandom3 rangen(SEED);

  //GET THE PTs
  TFile fInput_pT(TString::Format("%s/ReconstructionFiles/AnalysisResults_for_pT.root",INPUT_FOLDER.Data()), "read");
  TDirectory* fDir_temp = NULL;
  fDir_temp = (TDirectoryFile*)(fInput_pT.FindObjectAny("HMResults0"));
  TList* fList_temp1;
  fDir_temp->GetObject("HMResults0", fList_temp1);

  TList* fList_P02 = (TList*)fList_temp1->FindObject("Particle0_Particle2");
  TList* fList_P13 = (TList*)fList_temp1->FindObject("Particle1_Particle3");
  TList* fList_P03 = (TList*)fList_temp1->FindObject("Particle0_Particle3");
  TList* fList_P12 = (TList*)fList_temp1->FindObject("Particle1_Particle2");
  int which_mt_bin = abs(mt_bin%100) + 1;
  if(mt_bin<0) {which_mt_bin = 3;}

  printf("which_mt_bin = %i\n",which_mt_bin);
  TH2F* h_P02 = NULL;
  TH2F* h_P13 = NULL;
  TH2F* h_P03 = NULL;
  TH2F* h_P12 = NULL;

  h_P02 = (TH2F*)fList_P02->FindObject(TString::Format("MEmT_%i_pT_PionNucleon0_2_vs_kStar", which_mt_bin));
  h_P13 = (TH2F*)fList_P13->FindObject(TString::Format("MEmT_%i_pT_PionNucleon1_3_vs_kStar", which_mt_bin));
  h_P03 = (TH2F*)fList_P03->FindObject(TString::Format("MEmT_%i_pT_PionNucleon0_3_vs_kStar", which_mt_bin));
  h_P12 = (TH2F*)fList_P12->FindObject(TString::Format("MEmT_%i_pT_PionNucleon1_2_vs_kStar", which_mt_bin));
  

  double mean_pNucl_pT_below300=0;
  TH2F* h_pT_Nucl = NULL;
  if(abs(mt_bin)>=100){
    h_pT_Nucl = (TH2F*)h_P03->Clone("h_pT_Nucl");
    h_pT_Nucl->Add(h_P12);
  }
  else{
    h_pT_Nucl = (TH2F*)h_P02->Clone("h_pT_Nucl");
    h_pT_Nucl->Add(h_P13);
  }

  TH1F* hProj = NULL;

  hProj = (TH1F*)h_pT_Nucl->ProjectionY(TString::Format("h_pT_Nucl_Proj"),h_pT_Nucl->GetXaxis()->FindBin(0.5*0.001),h_pT_Nucl->GetXaxis()->FindBin(299.5*0.001));
  mean_pNucl_pT_below300 = hProj->GetMean()*1000.;
  
  delete hProj;

  TGraphErrors  gData(max_mom_bins);
  gData.SetName("gData");
  gData.SetLineColor(kBlack);
  gData.SetLineWidth(4);

  TGraph gFit(max_mom_bins);
  gFit.SetName("gFit");
  gFit.SetLineColor(kPink-9);
  gFit.SetLineWidth(4);

  TGraph gBaseline(max_mom_bins);
  gBaseline.SetName("gBaseline");
  gBaseline.SetLineColor(kGreen-5);
  gBaseline.SetLineWidth(4);

  TGraph gFemto(max_mom_bins);
  gFemto.SetName("gFemto");
  gFemto.SetLineColor(kBlue+2);
  gFemto.SetLineWidth(4);

  TGraph gDelta(max_mom_bins);
  gDelta.SetName("gDelta");
  gDelta.SetLineColor(kOrange-1);
  gDelta.SetLineWidth(4);

  float MT_LOW;
  float MT_UP;
  float RADIUS;
  int DATA_TYPE;
  int FIT_TYPE;
  int MT_ID = mt_bin;
  float FIT_RANGE;
  float LAM_GEN;
  float FRAC_D;
  float NUM_DELTAS;
  float AMP_DELTA;
  float DELTA_MASS;
  float DELTA_WIDTH;
  float PS_PT;
  float PS_TEMP;
  float CHI2_500;
  float PT_SCALE;


  if(mt_bin==-1){
    MT_LOW = mT_low.at(0).at(0);
    MT_UP = mT_up.at(0).at(mT_up.at(0).size()-1);
  }
  else if(mt_bin==-101){
    MT_LOW = mT_low.at(1).at(0);
    MT_UP = mT_up.at(1).at(mT_up.at(0).size()-1);
  }
  else if(mt_bin%100>=0 && mt_bin%100<=4){
    MT_LOW = mT_low.at(mt_bin/100).at(mt_bin%100);
    MT_UP = mT_up.at(mt_bin/100).at(mt_bin%100);    
  }
  else{
    printf("Silly mT bin\n");
    abort();
  }

  //int SEED, int mt_bin, int NumIter, bool Bootstrap=true, bool DataVar=true, bool FitVar=true
  TFile fOutputFile(TString::Format("%s/%s_mT%i_B%i_DV%i_FV%i_S%i.root",
    OUTPUT_FOLDER.Data(), Description, mt_bin, Bootstrap, DataVar, FitVar, SEED), "recreate");
  TTree* pi_d_Tree = new TTree("pi_d_Tree","pi_d_Tree");
  pi_d_Tree->Branch("gData","TGraphErrors",&gData,32000,0);//
  pi_d_Tree->Branch("gFit","TGraph",&gFit,32000,0);//
  pi_d_Tree->Branch("gBaseline","TGraph",&gBaseline,32000,0);//
  pi_d_Tree->Branch("gFemto","TGraph",&gFemto,32000,0);//
  pi_d_Tree->Branch("gDelta","TGraph",&gDelta,32000,0);//
  pi_d_Tree->Branch("mT_low",&MT_LOW,"mT_low/F");//
  pi_d_Tree->Branch("mT_up",&MT_UP,"mT_up/F");//
  pi_d_Tree->Branch("mT_id",&MT_ID,"mT_id/I");//
  pi_d_Tree->Branch("rad",&RADIUS,"rad/F");//
  pi_d_Tree->Branch("data_type",&DATA_TYPE,"data_type/I");//ABC, A=0 => no boot, A=1 = bootstrap, BC = variation ID. If negative => default
  pi_d_Tree->Branch("fit_type",&FIT_TYPE,"fit_type/I");
  pi_d_Tree->Branch("fit_range",&FIT_RANGE,"fit_range/F");//
  pi_d_Tree->Branch("lam_gen",&LAM_GEN,"lam_gen/F");//
  pi_d_Tree->Branch("frac_D",&FRAC_D,"frac_D/F");//
  pi_d_Tree->Branch("num_deltas",&NUM_DELTAS,"num_deltas/F");//
  pi_d_Tree->Branch("amp_delta",&AMP_DELTA,"amp_delta/F");//
  pi_d_Tree->Branch("CkCutOff",&CKCUTOFF,"CkCutOff/F");//
  pi_d_Tree->Branch("delta_mass",&DELTA_MASS,"delta_mass/F");//
  pi_d_Tree->Branch("delta_width",&DELTA_WIDTH,"delta_width/F");//
  pi_d_Tree->Branch("ps_pT",&PS_PT,"ps_pT/F");//
  pi_d_Tree->Branch("ps_Temp",&PS_TEMP,"ps_Temp/F");//
  pi_d_Tree->Branch("pT_scale",&PT_SCALE,"pT_scale/F");//
  pi_d_Tree->Branch("chi2_500",&CHI2_500,"chi2_500/F");//

  unsigned NumMomBins = 20;
  double kCatMin = 0;
  double kCatMax = 400;
  CATS Kitty;
  Kitty.SetMomBins(NumMomBins, kCatMin, kCatMax);
  DLM_CommonAnaFunctions AnalysisObject;
  AnalysisObject.SetCatsFilesFolder(TString::Format("%s/CatsFiles",INPUT_FOLDER.Data()).Data());

  if(fabs(mt_bin)>=100){
    AnalysisObject.SetUpCats_pi_d(Kitty, "Gauss", "Gauss", -1, 0);
    Kitty.SetQ1Q2(-1);
  }
  else{
    AnalysisObject.SetUpCats_pi_d(Kitty, "", "Gauss", 0, 0);
    Kitty.SetQ1Q2(+1);    
  }
  Kitty.SetAnaSource(0, eff_source_radii.at(0));
  if(Kitty.GetNumSourcePars()>1) Kitty.SetAnaSource(1, 2);
  Kitty.KillTheCat();
  Kitty.SetNotifications(CATS::nWarning);

  DLM_Ck CkKitty_Default(Kitty.GetNumSourcePars(),0,Kitty,max_mom_bins,0,max_kstar);
  CkKitty_Default.SetSourcePar(0,Kitty.GetAnaSourcePar(0));
  if(Kitty.GetNumSourcePars()>1) CkKitty_Default.SetSourcePar(1,Kitty.GetAnaSourcePar(1));
  DLM_Ck& CkKitty = CkKitty_Default;
  
  CkKitty.SetCutOff(CKCUTOFF,CkConv);
  CkKitty.Update();

  DLM_CkDecomposition CkDecKitty("pi_d",0,CkKitty,hResolution_pd);
  CkDecKitty.Update(true, true);
  ALICE_CkDec = &CkDecKitty;

  double Progress = 0;
  bool data_saved = false;
  for(int uIter=0; uIter<NumIter; uIter++){
    printf("uIter = %u\n",uIter);
    TFile* fInputData = NULL;
    TH1F* hData = NULL;
    TDirectoryFile* fDir = NULL;
    TList* fList1 = NULL;
    TList* fList2 = NULL;
    TH1F* hME = NULL;
    TH1F* hSE = NULL;

    //-1 is default, -2 is default without bootstrap
    int iDataVar = -1;
    do{
      if(DataVar==true){
        //if we get a -1, it is the default
        iDataVar = rangen.Integer(NumDataVars+1)-1;
      }
      if(!data_saved){
        iDataVar = -2;
        data_saved = true;
      }

      //mt int pip_d used for aprovals at CF and PF
      if(mt_bin == -1){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/InputFiles/CF_Piondeuteron_Var0.root",INPUT_FOLDER.Data());
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
          SeHistName = "hTotalSE _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/InputFiles/CF_Pid_Var%i.root",INPUT_FOLDER.Data(), iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
          SeHistName = TString::Format("hTotalSE_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }

      }
      //mt int pim_d used for approvals at CF and PF
      else if(mt_bin == -101){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/InputFiles/CF_AntiPiondeuteron_Var0.root",INPUT_FOLDER.Data());
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
          SeHistName = "hTotalSE _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/InputFiles/CF_AntiPid_Var%i.root",INPUT_FOLDER.Data(), iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);    
          SeHistName = TString::Format("hTotalSE_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);      
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      }
      //mt diff pip_d used for approvals at PF
      else if(mt_bin>=0 && mt_bin<100){
        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/InputFiles/mTBin%i/CF_Piondeuteron_Var0.root",INPUT_FOLDER.Data(),abs(mt_bin)%100+1);
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/InputFiles/mTBin%i/CF_Pid_Var%i.root",
            INPUT_FOLDER.Data(),abs(mt_bin)%100+1,iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      } 
      //mt diff pim_d used for approvals at PF
      else if(mt_bin>=100 && mt_bin<200){

        if(iDataVar==-1 || iDataVar==-2){
          InputDataFileName = TString::Format("%s/InputFiles/mTBin%i/CF_AntiPiondeuteron_Var0.root",
            INPUT_FOLDER.Data(),abs(mt_bin)%100+1);
          DefaultHistoName = "hCk_ReweightedPidVar0MeV_0";
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = "hTotalME _Rebinned_5_Reweighted";
        }
        else{
          InputDataFileName = TString::Format("%s/InputFiles/mTBin%i/CF_AntiPid_Var%i.root",
            INPUT_FOLDER.Data(),abs(mt_bin)%100+1,iDataVar+1);
          DefaultHistoName = TString::Format("hCk_ReweightedReweightedPidVar%iMeV_1MeV_0", iDataVar+1);
          MeListName1 = "PairDist";
          MeListName2 = "PairReweighted";
          MeHistName = TString::Format("hTotalME_PidVar%i_Rebinned_5_Reweighted", iDataVar+1);
        }
        std::ifstream file(InputDataFileName.Data());
        if(file.good()) fInputData = new TFile(InputDataFileName, "read");
        if(fInputData) hData = (TH1F*)fInputData->Get(DefaultHistoName);

        if(fInputData) fList1 = (TList*)(fInputData->FindObjectAny(MeListName1));
        if(fList1) fList2 = (TList*)fList1->FindObject(MeListName2);
        if(fList2) hME = (TH1F*)fList2->FindObject(MeHistName);
        if(hME){
          hME->GetXaxis()->SetLimits(hME->GetXaxis()->GetXmin()*1000.,hME->GetXaxis()->GetXmax()*1000.);
        }
        if(fList2) hSE = (TH1F*)fList2->FindObject(SeHistName);
        if(hSE){
          hSE->GetXaxis()->SetLimits(hSE->GetXaxis()->GetXmin()*1000.,hSE->GetXaxis()->GetXmax()*1000.);
        }
      }
    }
    while(hData==NULL);

    if(!hData){
      printf("!hData\n");
      abort();
    }

    gROOT->cd();
    TH1F* hDataToFit = (TH1F*)hData->Clone("hDataToFit");;
    if(Bootstrap==true && iDataVar!=-2){
      double MOM = 0;
      unsigned uMom = 0;
      while(MOM<max_kstar){
        MOM = hData->GetBinCenter(uMom+1);
        double new_value;
        do{
          new_value = rangen.Gaus(hData->GetBinContent(uMom+1), hData->GetBinError(uMom+1));
        }
        while(new_value<0);
        hDataToFit->SetBinContent(uMom+1, new_value);
        uMom++;
      }
    }

    int pt_scale_int = 0;
    FIT_RANGE = fit_range.at(0);
    RADIUS = eff_source_radii.at(0);
    PT_SCALE = pT_scale.at(pt_scale_int);
    FIT_TYPE = fit_types.at(0);
    if(mt_bin<0) LAM_GEN = lambda_gen_mTint;
    else LAM_GEN = lambda_gen_avg.at(mt_bin%100);
    AMP_DELTA = amplitude_delta.at(0);
    if(FitVar){
      int rndint = rangen.Integer(fit_range.size());
      FIT_RANGE = fit_range.at(rndint);
      RADIUS = eff_source_radii.at(rangen.Integer(eff_source_radii.size()));
      if(LAM_GEN>=0){
        if(mt_bin<0) LAM_GEN = lambda_gen_mTint*lambda_gen_relvar.at(rangen.Integer(lambda_gen_relvar.size()));
        else LAM_GEN = lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(rangen.Integer(lambda_gen_relvar.size()));
      }
      if(AMP_DELTA>=0) AMP_DELTA = amplitude_delta.at(rangen.Integer(amplitude_delta.size()));

      rndint = rangen.Integer(pT_scale.size());
      PT_SCALE = pT_scale.at(rndint);


      rndint = rangen.Integer(fit_types.size());
      FIT_TYPE = fit_types.at(rndint);
    }

    TFile fConversion(InputConversionFileName.at(pt_scale_int), "read");
    Kstar_Modifier = (TGraph*)fConversion.Get("gR_ppi_dpi_c");

    ALICE_CkDec->AddPhaseSpace(hME);

    DATA_TYPE = abs(iDataVar);
    if(Bootstrap==true) DATA_TYPE += 100;
    if(iDataVar<0) DATA_TYPE = -DATA_TYPE;
    //overrides the other stuff, this is default data with no boot
    if(iDataVar==-2) DATA_TYPE = -2;

    //par[0] = source size
    //par[1] = alpha par
    //par[2] = lambda_genuine
    //par[3] = lambda_d
    //par[4] = delta_amplitude
    //par[5] = CkCutOff
    //par[6] = CkConv
    //par[7] = just in case
    //par[8] = Delta normalization, should be a dummy fixed to -1e6, giving instructions to the code to renorm
    //par[9] = mass of daughter 1 (say the pion)
    //par[10] = mass of daughter 2 (say the proton)
    //par[11] = mass of the delta
    //par[12] = width of the delta
    //par[13] = avg pT of the daughters
    //par[14] = effective temperature
    //par[15] = just in case
    //par[16] = norm
    //par[17] = 0 to castrate the pol3
    //par[18] = position of the max of the pol3
    //par[19] = p3 parameter of the pol3
    //par[20] = -1e6 to switch off the pol4
    //[21-25] parameters of the extra resonances added during the review at Nature, see ALICE_fit for details 
    const int NumFitPar = 26;
    TF1* fData = new TF1("fData",ALICE_fit,0,FIT_RANGE,NumFitPar);
    TF1* fDataDummy = new TF1("fDataDummy",ALICE_fit,0,FIT_RANGE,NumFitPar);

    //minimum is one, if you whant to refine the FRAC_D recursively, increase it
    int RepeatFit = 3;
    int iRepeat = 0;
    FRAC_D = 1.5e-2;
    do{

      fData->FixParameter(0, RADIUS);
      fData->FixParameter(1, 2);
      if(LAM_GEN>=0) fData->FixParameter(2, LAM_GEN);
      else{
        if(mt_bin<0){
          fData->SetParameter(2, fabs(lambda_gen_mTint));
          fData->SetParLimits(2, fabs(lambda_gen_mTint*lambda_gen_relvar.at(0)), fabs(lambda_gen_mTint*lambda_gen_relvar.at(1)));        
        }
        else{
          fData->SetParameter(2, fabs(lambda_gen_avg.at(mt_bin%100)));
          fData->SetParLimits(2, fabs(lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(0)), fabs(lambda_gen_avg.at(mt_bin%100)*lambda_gen_relvar.at(1)));        
        }
      }

      if(AMP_DELTA>=0) fData->FixParameter(4, AMP_DELTA);
      else{
        fData->SetParameter(4, fabs(amplitude_delta.at(0)));
        fData->SetParLimits(4, fabs(amplitude_delta.at(1)), fabs(amplitude_delta.at(2)));
      }
      fData->FixParameter(5, CKCUTOFF);
      fData->FixParameter(6, CkConv);

      fData->FixParameter(7, FIT_TYPE);

      fData->FixParameter(8, -1e6);
      fData->FixParameter(9, Mass_pic);
      fData->FixParameter(10, Mass_p);
      //fData->SetParameter(11,1232);
      //fData->SetParLimits(11,1180,1260);
      fData->FixParameter(11,1215);
      fData->SetParameter(12,100);
      fData->SetParLimits(12,40,200);

      double kT_Delta;

      fData->SetParameter(13,mean_pNucl_pT_below300);
      fData->SetParLimits(13,mean_pNucl_pT_below300*0.8,mean_pNucl_pT_below300*1.2);

      fData->SetParameter(14,15);
      fData->SetParLimits(14,5,60);

      if(mt_bin<0){
        fData->FixParameter(12,pip_avg_gamma);//gamma
        fData->SetParameter(14,pip_avg_temp);//temp
      }
      else{
        fData->FixParameter(12,pip_gamma_Dpp.at(mt_bin%100)*Dpp_weight + pip_gamma_D0.at(mt_bin%100)*(1.-Dpp_weight));//gamma
        fData->SetParameter(14,pip_temp_Dpp.at(mt_bin%100)*Dpp_weight + pip_temp_D0.at(mt_bin%100)*(1.-Dpp_weight));//temp
      }

      fData->FixParameter(15, 0);

      fData->SetParameter(16,1);
      fData->SetParLimits(16,0.9,1.1);
      fData->FixParameter(17, 0);
      fData->SetParameter(18,400);
      fData->SetParLimits(18,1,100000);
      fData->SetParameter(19,0);
      fData->SetParLimits(19,-1e-8,1e-8);
      fData->FixParameter(20, -1e6);

      if(FIT_TYPE<10){
        fData->FixParameter(21, 0);
        fData->FixParameter(22, 1);
        fData->FixParameter(23, 1);
        fData->FixParameter(24, 1);
        fData->FixParameter(25, 1);
      }
      else{
        fData->SetParameter(21, 1e-4);
        fData->SetParLimits(21, 0,10);
        fData->FixParameter(22, Mass_d);
        fData->FixParameter(23, Mass_pic);
        fData->FixParameter(24, 2114);
        fData->FixParameter(25, 20);
      }
      fData->FixParameter(3,FRAC_D);
      hDataToFit->Fit(fData,"Q, S, N, R, M");   

      RADIUS = fData->GetParameter(0);
      LAM_GEN = fData->GetParameter(2);
      AMP_DELTA = fData->GetParameter(4);
      CKCUTOFF = fData->GetParameter(5);
      DELTA_MASS = fData->GetParameter(11);
      DELTA_WIDTH = fData->GetParameter(12);
      PS_PT = fData->GetParameter(13);
      PS_TEMP = fData->GetParameter(14);
      CHI2_500 = 0;
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>500.001) break;
        if(hDataToFit->GetBinError(uBin+1))
          CHI2_500 += pow((hDataToFit->GetBinContent(uBin+1) - fData->Eval(MOM))/hDataToFit->GetBinError(uBin+1),2.);
        else
          CHI2_500 += pow((hDataToFit->GetBinContent(uBin+1) - fData->Eval(MOM))/1.,2.);
      }

      clean_graph(&gData);
      clean_graph(&gFit);
      clean_graph(&gFemto);
      clean_graph(&gBaseline);
      clean_graph(&gDelta);

      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gData.SetPoint(uBin, MOM, hDataToFit->GetBinContent(uBin+1));
        gData.SetPointError(uBin, hDataToFit->GetBinWidth(uBin+1)*0.5, hDataToFit->GetBinError(uBin+1));
        gFit.SetPoint(uBin, MOM, fData->Eval(MOM));
      }

      fOutputFile.cd();
      for(int iPar=0; iPar<NumFitPar; iPar++){
        fDataDummy->FixParameter(iPar, fData->GetParameter(iPar));
      }
      hDataToFit->SetName(hDataToFit->GetName()+TString::Format("_uIter%i",uIter));
      fData->SetName(fData->GetName()+TString::Format("_uIter%i",uIter));

      double par_bl0 = fData->GetParameter(16);
      double par_bl1 = fData->GetParameter(17);
      double par_bl2 = fData->GetParameter(18);
      double par_bl3 = fData->GetParameter(19);
      double par_bl4 = fData->GetParameter(20);
      double par_fmt = fData->GetParameter(2);//lambda_gen
      double par_dlt = fData->GetParameter(4);
      double par_frc = fData->GetParameter(3);
      double par_er_n = fData->GetParameter(21);
      double par_er_Mm = fData->GetParameter(22);
      double par_er_Md1 = fData->GetParameter(23);
      double par_er_Md2 = fData->GetParameter(24);
      double par_er_w = fData->GetParameter(25);

      fDataDummy->FixParameter(2, 0);
      fDataDummy->FixParameter(4, 0);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gBaseline.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      fDataDummy->FixParameter(16, 1);
      fDataDummy->FixParameter(17, 0);
      fDataDummy->FixParameter(18, 0);
      fDataDummy->FixParameter(19, 0);
      fDataDummy->FixParameter(20, 0);
      fDataDummy->FixParameter(21, 0);
      fDataDummy->FixParameter(2, par_fmt);
      fDataDummy->FixParameter(4, 0);
      fDataDummy->FixParameter(3, 0);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gFemto.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      fDataDummy->FixParameter(2, 0);
      fDataDummy->FixParameter(3, par_frc);
      fDataDummy->FixParameter(4, par_dlt);
      for(unsigned uBin=0; uBin<hDataToFit->GetNbinsX(); uBin++){
        double MOM = hDataToFit->GetBinCenter(uBin+1);
        if(MOM>FIT_RANGE) break;
        gDelta.SetPoint(uBin, MOM, fDataDummy->Eval(MOM));
      }

      fDataDummy->FixParameter(4, 0);
      fDataDummy->FixParameter(21, par_er_n);


      //evaluate the integrals of the SE related to d from Delta and d not from Delta
      //the norm of the integral is not relevant, as we only care about the ratio FRAC_D
      double Integral_Primary = 0;
      double Integral_FromDelta = 0;
      FRAC_D = 0;
      if(hME){
        for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
          double MOM = hME->GetBinCenter(uBin+1);
          if(MOM<FIT_RANGE){
            Integral_Primary += hME->GetBinContent(uBin+1) * (gFemto.Eval(MOM) - (1.-par_fmt));//original
            Integral_FromDelta += hME->GetBinContent(uBin+1) * (gDelta.Eval(MOM)-1);
          }
          else{
            //we assume gFemto is 1, while gDelta is zero
            Integral_Primary += hME->GetBinContent(uBin+1);
          }
        }   
        FRAC_D = Integral_FromDelta/(Integral_FromDelta+Integral_Primary);   
      }
      iRepeat++;
    }
    while(iRepeat<RepeatFit);
    
    double IntME = hME->Integral();
    double IntSE = hSE?hSE->Integral():0;
    NUM_DELTAS = 0;
    for(unsigned uBin=0; uBin<hME->GetNbinsX(); uBin++){
      //checked, this is correct
      double MOM = hME->GetBinCenter(uBin+1);
      NUM_DELTAS += hME->GetBinContent(uBin+1) * (gDelta.Eval(MOM)-1);
    }

    NUM_DELTAS *= IntSE;
    NUM_DELTAS /= IntME;

    fOutputFile.cd();
    pi_d_Tree->Fill();

    delete fData;
    delete fDataDummy;
    delete hDataToFit;
    delete fInputData;
  }//for iter

  fOutputFile.cd();
  pi_d_Tree->Write();

  fOutputFile.Close();
}

