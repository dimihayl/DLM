#include <iostream>
#include <vector>

#include "DLM_FSI_Wrapper.h"
#include "DLM_Histo.h"

#include "TH2F.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TString.h"
#include "TMath.h"

using namespace std;

void Wrap_pp_Epelbaum(const char* InputFileName, const char* OutputFolder){
    const unsigned short num_channels = 12;

    const unsigned num_mom_bins = 80;
    double kShift = 2.5;
    double kWidth  = 5;
    double* mom_bins = new double [num_mom_bins+1];
    double* mom_bin_center = new double [num_mom_bins];
    mom_bins[0] = 0;
    for(int uMom=0; uMom<num_mom_bins; uMom++){
        //bin center = double(uMom)*kWidth+0.5*kWidth
        //shifted bin center = double(uMom)*kWidth+0.5*kWidth + kShift
        //low edge = double(uMom)*kWidth + kShift
        //up edge = double(uMom)*kWidth+1*kWidth + kShift        
        double(uMom)*kWidth+0.5*kWidth;
        mom_bins[uMom+1] = double(uMom+1)*kWidth + kShift;
        mom_bin_center[uMom] = double(uMom)*kWidth+0.5*kWidth + kShift;
        //printf("r %u %.3f %.3f %.3f\n", uMom, mom_bins[uMom], mom_bin_center[uMom], mom_bins[uMom+1]);
    }
    //return;
    const unsigned num_rad_bins = 999;
    double rShift = -0.025;
    double rWidth  = 0.05;
    double* rad_bins = new double [num_rad_bins+1];
    double* rad_bin_center = new double [num_rad_bins];
    rad_bins[0] = rShift;
    for(int uRad=0; uRad<num_rad_bins; uRad++){
        //bin center = double(uMom)*kWidth+0.5*kWidth
        //shifted bin center = double(uMom)*kWidth+0.5*kWidth + kShift
        //low edge = double(uMom)*kWidth + kShift
        //up edge = double(uMom)*kWidth+1*kWidth + kShift        
        double(uRad)*rWidth+0.5*rWidth;
        rad_bins[uRad+1] = double(uRad+1)*rWidth + rShift;
        rad_bin_center[uRad] = double(uRad)*rWidth+0.5*rWidth + rShift;
        //printf("r %u %.3f %.3f %.3f\n", uRad, rad_bins[uRad], rad_bin_center[uRad], rad_bins[uRad+1]);
    }
    //return;
    unsigned short* num_pw = new unsigned short [num_channels];
    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        switch(usCh){
            case  0: num_pw[usCh]=3; break;
            case  3: num_pw[usCh]=2; break;
            case  6: num_pw[usCh]=2; break;
            case  9: num_pw[usCh]=2; break;
            case 10: num_pw[usCh]=1; break;
            case 11: num_pw[usCh]=1; break;
            default: num_pw[usCh]=4; break;
        }
    }

    const unsigned num_dirs = 9;
    TString* dir_names = new TString [num_dirs];

    dir_names[0]="1S0";
    dir_names[1]="1D2";
    dir_names[2]="3P0";
    dir_names[3]="3P1";
    dir_names[4]="3P2";
    dir_names[5]="3F2";
    dir_names[6]="3F3";
    dir_names[7]="3F2-3P2";
    dir_names[8]="3P2-3F2";

    
    DLM_Histo<complex<double>>*** dlm_histo = NULL;
    dlm_histo = new DLM_Histo<complex<double>>** [2];
    DLM_Histo<complex<double>>**& HistoWF = dlm_histo[0];
    DLM_Histo<complex<double>>**& HistoPS = dlm_histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [num_channels];
    HistoPS = new DLM_Histo<complex<double>> *[num_channels];
    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        HistoWF[usCh] = new DLM_Histo<complex<double>> [num_pw[usCh]];
        HistoPS[usCh] = new DLM_Histo<complex<double>> [num_pw[usCh]];
        for(unsigned short usPw=0; usPw<num_pw[usCh]; usPw++){
            if(usCh==0 && usPw%2==1) continue;
            if((usCh>=1&&usCh<=9) && usPw%2==0) continue;
            
            HistoWF[usCh][usPw].SetUp(2);
            HistoWF[usCh][usPw].SetUp(0,num_mom_bins,mom_bins,mom_bin_center);
            HistoWF[usCh][usPw].SetUp(1,num_rad_bins,rad_bins,rad_bin_center);
            HistoWF[usCh][usPw].Initialize();

            HistoPS[usCh][usPw].SetUp(1);
            HistoPS[usCh][usPw].SetUp(0,num_mom_bins,mom_bins,mom_bin_center);
            HistoPS[usCh][usPw].Initialize();
        }
    }


    TFile fInput(InputFileName, "read");
    for(unsigned uDir=0; uDir<num_dirs; uDir++){
        TDirectoryFile *dir1=(TDirectoryFile*)(fInput.FindObjectAny(dir_names[uDir].Data()));
        if(!dir1){
            printf("\033[1;31mERROR:\033[0m Wrap_pp_Epelbaum cannot open the directory %s\n", dir_names[uDir].Data());
            continue;
        }
        TDirectoryFile *dir2=(TDirectoryFile*)(dir1->FindObjectAny("wave_vs_distance"));
        if(!dir2){
            printf("\033[1;31mERROR:\033[0m Wrap_pp_Epelbaum cannot open the directory wave_vs_distance\n");
            continue;
        }
        for(int uMom=0; uMom<num_mom_bins; uMom++){
            double& mom = mom_bin_center[uMom];
            TH1F *hist_wf_real=(TH1F*)(dir2->FindObjectAny(TString::Format("hist_wf_real_vs_r_p%i",TMath::Nint(mom))));
            TH1F *hist_wf_imag=(TH1F*)(dir2->FindObjectAny(TString::Format("hist_wf_imag_vs_r_p%i",TMath::Nint(mom))));
            if(!hist_wf_real){
                printf("\033[1;31mERROR:\033[0m Wrap_pp_Epelbaum cannot open the histogram %s\n", TString::Format("hist_wf_real_vs_r_p%i",TMath::Nint(mom)).Data());
                continue;
            }
            if(!hist_wf_imag){
                printf("\033[1;31mERROR:\033[0m Wrap_pp_Epelbaum cannot open the histogram %s\n", TString::Format("hist_wf_imag_vs_r_p%i",TMath::Nint(mom)).Data());
                continue;
            }
            for(int uRad=0; uRad<num_rad_bins; uRad++){
                double& rad = rad_bin_center[uRad];
                complex<double> wf_val;// = (hist_wf_real->GetBinContent(uRad+1)*rad, hist_wf_imag->GetBinContent(uRad+1)*rad);
                wf_val.real(hist_wf_real->GetBinContent(uRad+1)*rad);
                wf_val.imag(hist_wf_imag->GetBinContent(uRad+1)*rad);
                //printf("%u %f %f\n",uRad,rad,hist_wf_real->GetBinCenter(uRad+1));
                //usleep(100e3);
                if(fabs(rad-hist_wf_real->GetBinCenter(uRad+1))>1e-3 || fabs(rad-hist_wf_imag->GetBinCenter(uRad+1))>1e-3){
                    printf("\033[1;31mERROR:\033[0m Wrap_pp_Epelbaum found a mismatch of the binning in r*\n");
                    break;
                }
                if(dir_names[uDir]=="1S0"){
                    HistoWF[0][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[0][0].SetBinContent(uMom,0);
                }
                else if(dir_names[uDir]=="1D2"){
                    HistoWF[0][2].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[0][2].SetBinContent(uMom,0);                    
                }
                else if(dir_names[uDir]=="3P0"){
                    HistoWF[1][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[1][1].SetBinContent(uMom,0); 
                    HistoWF[2][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[2][1].SetBinContent(uMom,0);
                    HistoWF[3][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[3][1].SetBinContent(uMom,0);                  
                }
                else if(dir_names[uDir]=="3P1"){
                    HistoWF[4][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[4][1].SetBinContent(uMom,0); 
                    HistoWF[5][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[5][1].SetBinContent(uMom,0);
                    HistoWF[6][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[6][1].SetBinContent(uMom,0);                    
                }
                else if(dir_names[uDir]=="3P2"){
                    HistoWF[7][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[7][1].SetBinContent(uMom,0); 
                    HistoWF[8][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[8][1].SetBinContent(uMom,0);
                    HistoWF[9][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[9][1].SetBinContent(uMom,0);                      
                }
                else if(dir_names[uDir]=="3F2"){
                    HistoWF[1][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[1][3].SetBinContent(uMom,0); 
                    HistoWF[4][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[4][3].SetBinContent(uMom,0);
                    HistoWF[7][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[7][3].SetBinContent(uMom,0);                        
                }
                else if(dir_names[uDir]=="3F3"){
                    HistoWF[2][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[2][3].SetBinContent(uMom,0); 
                    HistoWF[5][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[5][3].SetBinContent(uMom,0);
                    HistoWF[8][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[8][3].SetBinContent(uMom,0);                        
                }
                else if(dir_names[uDir]=="3F4"){
                    //
                }
                else if(dir_names[uDir]=="3F2-3P2"){
                    HistoWF[10][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[10][0].SetBinContent(uMom,0);                     
                }
                else if(dir_names[uDir]=="3P2-3F2"){
                    HistoWF[11][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[11][0].SetBinContent(uMom,0);    
                }
                else{
                    printf("\033[1;31mERROR:\033[0m BUG in Wrap_pp_Epelbaum\n");
                    continue;
                }
            }
        }
    }

    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        for(unsigned short usPw=0; usPw<num_pw[usCh]; usPw++){
            if(HistoWF[usCh][usPw].IsInitialized()){
                TString WF_file_name = TString::Format("%s/wf_epel_usCh%u_usPw%u.dlm",OutputFolder,usCh,usPw);
                TString PS_file_name = TString::Format("%s/ps_epel_usCh%u_usPw%u.dlm",OutputFolder,usCh,usPw);
                HistoWF[usCh][usPw].QuickWrite(WF_file_name.Data(),true);
                HistoPS[usCh][usPw].QuickWrite(PS_file_name.Data(),true);
            }
            //delete dlm_histo[0][usCh][usPw];
            //delete dlm_histo[1][usCh][usPw];
        }
        delete [] dlm_histo[0][usCh];
        delete [] dlm_histo[1][usCh];
    }
    delete [] dlm_histo[0];
    delete [] dlm_histo[1];
    delete [] dlm_histo;

    //CleanUp_Wrap_pp_Epelbaum:;
    delete [] mom_bins;
    delete [] mom_bin_center;
    delete [] rad_bins;
    delete [] rad_bin_center;
    delete [] num_pw;

}



void Wrap_pp_Norfolk(const char* InputFileName, const char* OutputFolder, const char* descriptor){

    const unsigned short num_channels = 13;

    const unsigned num_mom_bins = 250;
    double kShift = 1;
    double kWidth  = 2;
    double* mom_bins = new double [num_mom_bins+1];
    double* mom_bin_center = new double [num_mom_bins];
    mom_bins[0] = 0;
    for(int uMom=0; uMom<num_mom_bins; uMom++){
        //bin center = double(uMom)*kWidth+0.5*kWidth
        //shifted bin center = double(uMom)*kWidth+0.5*kWidth + kShift
        //low edge = double(uMom)*kWidth + kShift
        //up edge = double(uMom)*kWidth+1*kWidth + kShift        
        double(uMom)*kWidth+0.5*kWidth;
        mom_bins[uMom+1] = double(uMom+1)*kWidth + kShift;
        mom_bin_center[uMom] = double(uMom)*kWidth+0.5*kWidth + kShift;
        //printf("r %u %.3f %.3f %.3f\n", uMom, mom_bins[uMom], mom_bin_center[uMom], mom_bins[uMom+1]);
    }
    //return;
    const unsigned num_rad_bins = 3017;
    double rShift = -0.01;
    double rWidth  = 0.02;
    double* rad_bins = new double [num_rad_bins+1];
    double* rad_bin_center = new double [num_rad_bins];
    rad_bins[0] = rShift;
    for(int uRad=0; uRad<num_rad_bins; uRad++){
        //bin center = double(uMom)*kWidth+0.5*kWidth
        //shifted bin center = double(uMom)*kWidth+0.5*kWidth + kShift
        //low edge = double(uMom)*kWidth + kShift
        //up edge = double(uMom)*kWidth+1*kWidth + kShift        
        double(uRad)*rWidth+0.5*rWidth;
        rad_bins[uRad+1] = double(uRad+1)*rWidth + rShift;
        rad_bin_center[uRad] = double(uRad)*rWidth+0.5*rWidth + rShift;
        //printf("r %u %.3f %.3f %.3f\n", uRad, rad_bins[uRad], rad_bin_center[uRad], rad_bins[uRad+1]);
    }
    //return;
    unsigned short* num_pw = new unsigned short [num_channels];
    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        switch(usCh){
            case  0: num_pw[usCh]=3; break;
            case 10: num_pw[usCh]=1; break;
            case 11: num_pw[usCh]=1; break;
            case 12: num_pw[usCh]=1; break;
            default: num_pw[usCh]=4; break;
        }
    }

    const unsigned num_dirs = 11;
    TString* dir_names = new TString [num_dirs];

    dir_names[0]="1S0";
    dir_names[1]="1D2";
    dir_names[2]="3P0";
    dir_names[3]="3P1";
    dir_names[4]="3P2";
    dir_names[5]="3F2";
    dir_names[6]="3F3";
    dir_names[7]="3F4";
    dir_names[8]="3F2-3P2";
    dir_names[9]="3P2-3F2";
    dir_names[10]="3H4-3F4";

    
    DLM_Histo<complex<double>>*** dlm_histo = NULL;
    dlm_histo = new DLM_Histo<complex<double>>** [2];
    DLM_Histo<complex<double>>**& HistoWF = dlm_histo[0];
    DLM_Histo<complex<double>>**& HistoPS = dlm_histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [num_channels];
    HistoPS = new DLM_Histo<complex<double>> *[num_channels];
    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        HistoWF[usCh] = new DLM_Histo<complex<double>> [num_pw[usCh]];
        HistoPS[usCh] = new DLM_Histo<complex<double>> [num_pw[usCh]];
        for(unsigned short usPw=0; usPw<num_pw[usCh]; usPw++){
            if(usCh==0 && usPw%2==1) continue;
            if((usCh>=1&&usCh<=9) && usPw%2==0) continue;
            
            HistoWF[usCh][usPw].SetUp(2);
            HistoWF[usCh][usPw].SetUp(0,num_mom_bins,mom_bins,mom_bin_center);
            HistoWF[usCh][usPw].SetUp(1,num_rad_bins,rad_bins,rad_bin_center);
            HistoWF[usCh][usPw].Initialize();

            HistoPS[usCh][usPw].SetUp(1);
            HistoPS[usCh][usPw].SetUp(0,num_mom_bins,mom_bins,mom_bin_center);
            HistoPS[usCh][usPw].Initialize();
        }
    }

    TFile fInput(InputFileName, "read");
    for(unsigned uDir=0; uDir<num_dirs; uDir++){
        TDirectoryFile *dir1=(TDirectoryFile*)(fInput.FindObjectAny(dir_names[uDir].Data()));
        if(!dir1){
            printf("\033[1;31mERROR:\033[0m Wrap_pp_Norfolk cannot open the directory %s\n", dir_names[uDir].Data());
            continue;
        }
        TDirectoryFile *dir2=(TDirectoryFile*)(dir1->FindObjectAny("wave_vs_distance"));
        if(!dir2){
            printf("\033[1;31mERROR:\033[0m Wrap_pp_Norfolk cannot open the directory wave_vs_distance\n");
            continue;
        }
        for(int uMom=0; uMom<num_mom_bins; uMom++){
            double& mom = mom_bin_center[uMom];
            TH1F *hist_wf_real=(TH1F*)(dir2->FindObjectAny(TString::Format("hist_wf_real_vs_r_p%i",TMath::Nint(mom))));
            TH1F *hist_wf_imag=(TH1F*)(dir2->FindObjectAny(TString::Format("hist_wf_imag_vs_r_p%i",TMath::Nint(mom))));

            if(!hist_wf_real){
                printf("\033[1;31mERROR:\033[0m Wrap_pp_Norfolk cannot open the histogram %s\n", TString::Format("hist_wf_real_vs_r_p%i",TMath::Nint(mom)).Data());
                continue;
            }
            if(!hist_wf_imag){
                printf("\033[1;31mERROR:\033[0m Wrap_pp_Norfolk cannot open the histogram %s\n", TString::Format("hist_wf_imag_vs_r_p%i",TMath::Nint(mom)).Data());
                continue;
            }
            for(int uRad=0; uRad<num_rad_bins; uRad++){
                double& rad = rad_bin_center[uRad];
                complex<double> wf_val;// = (hist_wf_real->GetBinContent(uRad+1)*rad, hist_wf_imag->GetBinContent(uRad+1)*rad);
                wf_val.real(hist_wf_real->GetBinContent(uRad+1)*rad);
                wf_val.imag(hist_wf_imag->GetBinContent(uRad+1)*rad);
                if(fabs(rad-hist_wf_real->GetBinCenter(uRad+1))>1e-3 || fabs(rad-hist_wf_imag->GetBinCenter(uRad+1))>1e-3){
                    printf("\033[1;31mERROR:\033[0m Wrap_pp_Norfolk found a mismatch of the binning in r*\n");
                    break;
                }
                if(dir_names[uDir]=="1S0"){
                    HistoWF[0][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[0][0].SetBinContent(uMom,0);
                }
                else if(dir_names[uDir]=="1D2"){
                    HistoWF[0][2].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[0][2].SetBinContent(uMom,0);                    
                }
                else if(dir_names[uDir]=="3P0"){
                    HistoWF[1][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[1][1].SetBinContent(uMom,0); 
                    HistoWF[2][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[2][1].SetBinContent(uMom,0);
                    HistoWF[3][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[3][1].SetBinContent(uMom,0);                  
                }
                else if(dir_names[uDir]=="3P1"){
                    HistoWF[4][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[4][1].SetBinContent(uMom,0); 
                    HistoWF[5][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[5][1].SetBinContent(uMom,0);
                    HistoWF[6][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[6][1].SetBinContent(uMom,0);                    
                }
                else if(dir_names[uDir]=="3P2"){
                    HistoWF[7][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[7][1].SetBinContent(uMom,0); 
                    HistoWF[8][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[8][1].SetBinContent(uMom,0);
                    HistoWF[9][1].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[9][1].SetBinContent(uMom,0);                      
                }
                else if(dir_names[uDir]=="3F2"){
                    HistoWF[1][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[1][3].SetBinContent(uMom,0); 
                    HistoWF[4][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[4][3].SetBinContent(uMom,0);
                    HistoWF[7][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[7][3].SetBinContent(uMom,0);                        
                }
                else if(dir_names[uDir]=="3F3"){
                    HistoWF[2][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[2][3].SetBinContent(uMom,0); 
                    HistoWF[5][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[5][3].SetBinContent(uMom,0);
                    HistoWF[8][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[8][3].SetBinContent(uMom,0);                        
                }
                else if(dir_names[uDir]=="3F4"){
                    HistoWF[3][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[3][3].SetBinContent(uMom,0); 
                    HistoWF[6][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[6][3].SetBinContent(uMom,0);
                    HistoWF[9][3].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[9][3].SetBinContent(uMom,0);                        
                }
                else if(dir_names[uDir]=="3F2-3P2"){
                    HistoWF[10][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[10][0].SetBinContent(uMom,0);                     
                }
                else if(dir_names[uDir]=="3P2-3F2"){
                    HistoWF[11][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[11][0].SetBinContent(uMom,0);    
                }
                else if(dir_names[uDir]=="3H4-3F4"){
                    HistoWF[12][0].SetBinContent(uMom,uRad,wf_val);
                    HistoPS[12][0].SetBinContent(uMom,0);    
                }
                else{
                    printf("\033[1;31mERROR:\033[0m BUG in Wrap_pp_Epelbaum\n");
                    continue;
                }
            }
        }
    }

    for(unsigned short usCh=0; usCh<num_channels; usCh++){
        for(unsigned short usPw=0; usPw<num_pw[usCh]; usPw++){
            if(HistoWF[usCh][usPw].IsInitialized()){
                TString WF_file_name = TString::Format("%s/wf_%s_usCh%u_usPw%u.dlm",OutputFolder,descriptor,usCh,usPw);
                TString PS_file_name = TString::Format("%s/ps_%s_usCh%u_usPw%u.dlm",OutputFolder,descriptor,usCh,usPw);
                HistoWF[usCh][usPw].QuickWrite(WF_file_name.Data(),true);
                HistoPS[usCh][usPw].QuickWrite(PS_file_name.Data(),true);
            }
            //delete dlm_histo[0][usCh][usPw];
            //delete dlm_histo[1][usCh][usPw];
        }
        delete [] dlm_histo[0][usCh];
        delete [] dlm_histo[1][usCh];
    }
    delete [] dlm_histo[0];
    delete [] dlm_histo[1];
    delete [] dlm_histo;

    //CleanUp_Wrap_pp_Epelbaum:;
    delete [] mom_bins;
    delete [] mom_bin_center;
    delete [] rad_bins;
    delete [] rad_bin_center;    
    delete [] num_pw;

}