
#include "DLM_RootWrapper.h"
#include "DLM_Histo.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TROOT.h"
#include <unistd.h>


DLM_Histo<float>* Convert_TH2F_DlmHisto(const TH2F* input){
    if(!input){
        printf("\033[1;31mERROR:\033[0m Bad input into Convert_TH2F_DlmHisto\n");
        return NULL;
    }
    DLM_Histo<float>* dlmHisto;
    dlmHisto = new DLM_Histo<float> ();
    dlmHisto->SetUp(2);
    unsigned NumBinsX = input->GetXaxis()->GetNbins();
    double* Xaxis_pars = new double [NumBinsX+1];
    for(unsigned uBin=0; uBin<NumBinsX; uBin++){
        Xaxis_pars[uBin] = input->GetXaxis()->GetBinLowEdge(uBin+1);
    }
    Xaxis_pars[NumBinsX] = input->GetXaxis()->GetBinUpEdge(NumBinsX);
    unsigned NumBinsY = input->GetYaxis()->GetNbins();
    double* Yaxis_pars = new double [NumBinsY+1];
    for(unsigned uBin=0; uBin<NumBinsY; uBin++){
        Yaxis_pars[uBin] = input->GetYaxis()->GetBinLowEdge(uBin+1);
    }
    Yaxis_pars[NumBinsY] = input->GetYaxis()->GetBinUpEdge(NumBinsY);
    dlmHisto->SetUp(0,NumBinsX,Xaxis_pars);
    dlmHisto->SetUp(1,NumBinsY,Yaxis_pars);
    dlmHisto->Initialize(true);

    unsigned WhichBin[2];
    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
            WhichBin[0] = uBinX;
            WhichBin[1] = uBinY;
            if( input->GetBinContent(uBinX+1,uBinY+1)!=input->GetBinContent(uBinX+1,uBinY+1) ||
                isinf(input->GetBinContent(uBinX+1,uBinY+1))){
              dlmHisto->SetBinContent(WhichBin,0);
              dlmHisto->SetBinError(WhichBin,0);
            }
            else{
              dlmHisto->SetBinContent(WhichBin,input->GetBinContent(uBinX+1,uBinY+1));
              dlmHisto->SetBinError(WhichBin,input->GetBinError(uBinX+1,uBinY+1));
            }
        }
    }
    delete [] Xaxis_pars;
    delete [] Yaxis_pars;
    return dlmHisto;
}

DLM_Histo<float>* Convert_TH1F_DlmHisto(const TH1F* input){
    if(!input){
        printf("\033[1;31mERROR:\033[0m Bad input into Convert_TH1F_DlmHisto\n");
        return NULL;
    }
    DLM_Histo<float>* dlmHisto;
    dlmHisto = new DLM_Histo<float> ();
    dlmHisto->SetUp(1);
    unsigned NumBinsX = input->GetXaxis()->GetNbins();
    double* Xaxis_pars = new double [NumBinsX+1];
    for(unsigned uBin=0; uBin<NumBinsX; uBin++){
        Xaxis_pars[uBin] = input->GetXaxis()->GetBinLowEdge(uBin+1);
    }
    Xaxis_pars[NumBinsX] = input->GetXaxis()->GetBinUpEdge(NumBinsX);
    dlmHisto->SetUp(0,NumBinsX,Xaxis_pars);
    dlmHisto->Initialize(true);

    unsigned WhichBin[1];
    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        WhichBin[0] = uBinX;
        if( input->GetBinContent(uBinX+1)!=input->GetBinContent(uBinX+1) ||
            isinf(input->GetBinContent(uBinX+1))){
          dlmHisto->SetBinContent(WhichBin,0);
          dlmHisto->SetBinError(WhichBin,0);
        }
        else{
          dlmHisto->SetBinContent(WhichBin,input->GetBinContent(uBinX+1));
          dlmHisto->SetBinError(WhichBin,input->GetBinError(uBinX+1));
        }
    }
    delete [] Xaxis_pars;
    return dlmHisto;
}

DLM_Histo<double>* Convert_TH1F_DoubleDlmHisto(const TH1F* input){
    if(!input){
        printf("\033[1;31mERROR:\033[0m Bad input into Convert_TH1F_DlmHisto\n");
        return NULL;
    }
    DLM_Histo<double>* dlmHisto;
    dlmHisto = new DLM_Histo<double> ();
    dlmHisto->SetUp(1);
    unsigned NumBinsX = input->GetXaxis()->GetNbins();
    double* Xaxis_pars = new double [NumBinsX+1];
    for(unsigned uBin=0; uBin<NumBinsX; uBin++){
        Xaxis_pars[uBin] = input->GetXaxis()->GetBinLowEdge(uBin+1);
    }
    Xaxis_pars[NumBinsX] = input->GetXaxis()->GetBinUpEdge(NumBinsX);
    dlmHisto->SetUp(0,NumBinsX,Xaxis_pars);
    dlmHisto->Initialize(true);

    unsigned WhichBin[1];
    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        WhichBin[0] = uBinX;
        if( input->GetBinContent(uBinX+1)!=input->GetBinContent(uBinX+1) ||
            isinf(input->GetBinContent(uBinX+1))){
          dlmHisto->SetBinContent(WhichBin,0);
          dlmHisto->SetBinError(WhichBin,0);
        }
        else{
          dlmHisto->SetBinContent(WhichBin,input->GetBinContent(uBinX+1));
          dlmHisto->SetBinError(WhichBin,input->GetBinError(uBinX+1));
        }
    }
    delete [] Xaxis_pars;
    return dlmHisto;
}


TH2F* Convert_DlmHisto_TH2F(const DLM_Histo<float>* input, const char* name){
    if(!input || !name || input->GetDim()!=2){
        printf("\033[1;31mERROR:\033[0m Bad input into Convert_DlmHisto_TH2F\n");
        return NULL;
    }
    unsigned NumBinsX = input->GetNbins(0);
    unsigned NumBinsY = input->GetNbins(1);
    double* BinsX = input->GetBinRange(0);
    double* BinsY = input->GetBinRange(1);
    gROOT->cd();
    TH2F* hROOT = new TH2F(name,name,NumBinsX,BinsX,NumBinsY,BinsY);
    unsigned WhichBin[2];
    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        for(unsigned uBinY=0; uBinY<NumBinsY; uBinY++){
            WhichBin[0] = uBinX;
            WhichBin[1] = uBinY;
            if( input->GetBinContent(WhichBin)!=input->GetBinContent(WhichBin) ||
                isinf(input->GetBinContent(WhichBin))){
              hROOT->SetBinContent(uBinX+1,uBinY+1,0);
              hROOT->SetBinError(uBinX+1,uBinY+1,0);
            }
            else{
              hROOT->SetBinContent(uBinX+1,uBinY+1,input->GetBinContent(WhichBin));
              hROOT->SetBinError(uBinX+1,uBinY+1,input->GetBinError(WhichBin));
            }
        }
    }
    return hROOT;
}
TH1F* Convert_DlmHisto_TH1F(const DLM_Histo<float>* input, const char* name){
    if(!input || !name || input->GetDim()!=1){
        printf("\033[1;31mERROR:\033[0m Bad input into Convert_DlmHisto_TH1F\n");
        return NULL;
    }
    unsigned NumBinsX = input->GetNbins(0);
    double* BinsX = input->GetBinRange(0);
    gROOT->cd();
    TH1F* hROOT = new TH1F(name,name,NumBinsX,BinsX);
    unsigned WhichBin[1];
    for(unsigned uBinX=0; uBinX<NumBinsX; uBinX++){
        WhichBin[0] = uBinX;
        //printf("%u %e %e\n",uBinX,input->GetBinContent(WhichBin),input->GetBinError(WhichBin));
        if( input->GetBinContent(WhichBin)!=input->GetBinContent(WhichBin) ||
            isinf(input->GetBinContent(WhichBin))){
          hROOT->SetBinContent(uBinX+1,0);
          hROOT->SetBinError(uBinX+1,0);
        }
        else{
          hROOT->SetBinContent(uBinX+1,input->GetBinContent(WhichBin));
          hROOT->SetBinError(uBinX+1,input->GetBinError(WhichBin));
        }
    }
    return hROOT;
}
