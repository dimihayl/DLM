#include "DLM_Fitters.h"

#include "TString.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

//for test only
#include "TFile.h"

#include "DLM_CkDecomposition.h"

DLM_Fitter1::DLM_Fitter1(const unsigned& maxnumsyst):MaxNumSyst(maxnumsyst),NumPar(15),NumRangePar(4){
    HistoOriginal = new const TH1F* [MaxNumSyst];
    HistoToFit = new TH1F* [MaxNumSyst];
    SystemToFit = new DLM_CkDecomposition* [MaxNumSyst];
    SourceSystems = NULL;
    PotentialSystems = NULL;
    ParentSource = NULL;
    ParentParameter = new int [MaxNumSyst*NumPar];
    for(unsigned uPar=0; uPar<MaxNumSyst*NumPar; uPar++) ParentParameter[uPar] = -1;
    ParentPotential = NULL;
    NumSourceSystems = 0;
    NumPotentialSystems = 0;
    FitRange = new double* [MaxNumSyst];
    NumSourceMapEntries = 0;
    SameSourceMap = NULL;
    NumSameSourcePar = NULL;
    NumPotentialMapEntries = 0;
    SamePotentialMap = NULL;
    NumSamePotentialPar = NULL;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        HistoOriginal[uSyst] = NULL;
        HistoToFit[uSyst] = NULL;
        SystemToFit[uSyst] = NULL;
        FitRange[uSyst] = new double [NumRangePar];
    }

    ParValue = new double* [MaxNumSyst];
    ParDownLimit = new double* [MaxNumSyst];
    ParUpLimit = new double* [MaxNumSyst];
    FixPar = new bool* [MaxNumSyst];
    SeparateBaseLineFit = new bool [MaxNumSyst];
    FitBL = new TF1* [MaxNumSyst];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        ParValue[uSyst] = new double [NumPar];
        ParDownLimit[uSyst] = new double [NumPar];
        ParUpLimit[uSyst] = new double [NumPar];
        FixPar[uSyst] = new bool [NumPar];
        //fix the pot pars
        for(unsigned uPar=0; uPar<NumPar; uPar++){
            FixPar[uSyst][uPar] = false;
        }
        SeparateBaseLineFit[uSyst] = true;
        ParValue[uSyst][p_a] = 1; ParDownLimit[uSyst][p_a] = 0.5; ParUpLimit[uSyst][p_a] = 2;
        ParValue[uSyst][p_b] = 1e-3; ParDownLimit[uSyst][p_b] = 1e-5; ParUpLimit[uSyst][p_b] = 1;
        ParValue[uSyst][p_c] = 0; ParDownLimit[uSyst][p_c] = 0; ParUpLimit[uSyst][p_c] = 0;
        FixPar[uSyst][p_c] = true;
        ParValue[uSyst][p_sor0] = 1.5; ParDownLimit[uSyst][p_sor0] = 1; ParUpLimit[uSyst][p_sor0] = 2.5;
        ParValue[uSyst][p_sor1] = 0; ParValue[uSyst][p_sor2] = 0; ParValue[uSyst][p_sor3] = 0; ParValue[uSyst][p_sor4] = 0; ParValue[uSyst][p_sor5] = 0;
        FixPar[uSyst][p_sor0]=true; FixPar[uSyst][p_sor1]=true; FixPar[uSyst][p_sor2]=true; FixPar[uSyst][p_sor3]=true; FixPar[uSyst][p_sor4]=true; FixPar[uSyst][p_sor5]=true;
        ParValue[uSyst][p_Cl] = 1; ParDownLimit[uSyst][p_Cl] = 0.75; ParUpLimit[uSyst][p_Cl] = 1.25;
        ParValue[uSyst][p_kc] = 300; ParDownLimit[uSyst][p_kc] = 150; ParUpLimit[uSyst][p_kc] = 900;
        ParValue[uSyst][p_pot0]=0; ParValue[uSyst][p_pot1]=0; ParValue[uSyst][p_pot2]=0; ParValue[uSyst][p_pot3]=0;
        FixPar[uSyst][p_pot0]=true; FixPar[uSyst][p_pot1]=true; FixPar[uSyst][p_pot2]=true; FixPar[uSyst][p_pot3]=true;
        //FitBL[uSyst] = new TF1(TString::Format("FitBL%u",uSyst),"[0]+[1]*x+[2]*x*x",FitRange[uSyst][kl],FitRange[uSyst][kmax]);
        FitBL[uSyst] = NULL;
    }

    HistoGlobal = NULL;
    FitGlobal = NULL;
    GlobalToMomentum = NULL;
    NumBinsSyst = NULL;
    CumulativeNumBinsSyst = NULL;

    OutputDirName = "./";
    RemoveNegCk = false;

}


DLM_Fitter1::~DLM_Fitter1(){

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(HistoToFit[uSyst]) {delete HistoToFit[uSyst]; HistoToFit[uSyst]=NULL;}
        if(FitRange[uSyst]) {delete [] FitRange[uSyst]; FitRange[uSyst]=NULL;}
        //if(SameSourceMap[uSyst]) {delete [] SameSourceMap[uSyst]; SameSourceMap[uSyst]=NULL;}
    }

    delete [] HistoOriginal; HistoOriginal=NULL;
    delete [] HistoToFit; HistoToFit=NULL;
    delete [] SystemToFit; SystemToFit=NULL;
    delete [] FitRange; FitRange=NULL;

    if(SameSourceMap){
        for(unsigned uEntry=0; uEntry<NumSourceMapEntries; uEntry++){
            delete [] SameSourceMap[uEntry];
        }
        delete [] SameSourceMap; SameSourceMap=NULL;
    }
    if(NumSameSourcePar) {delete [] NumSameSourcePar; NumSameSourcePar=NULL;}

    if(SamePotentialMap){
        for(unsigned uEntry=0; uEntry<NumPotentialMapEntries; uEntry++){
            delete [] SamePotentialMap[uEntry];
        }
        delete [] SamePotentialMap; SamePotentialMap=NULL;
    }
    if(NumSamePotentialPar) {delete [] NumSamePotentialPar; NumSamePotentialPar=NULL;}

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        delete [] ParValue[uSyst];
        delete [] ParDownLimit[uSyst];
        delete [] ParUpLimit[uSyst];
        delete [] FixPar[uSyst];
        if(FitBL[uSyst]) {delete FitBL[uSyst]; FitBL[uSyst]=NULL;}
    }
    delete [] ParValue; ParValue=NULL;
    delete [] ParDownLimit; ParDownLimit=NULL;
    delete [] ParUpLimit; ParUpLimit=NULL;
    delete [] FixPar; FixPar=NULL;
    delete [] SeparateBaseLineFit; SeparateBaseLineFit=NULL;
    delete FitBL; FitBL=NULL;

    if(HistoGlobal){delete HistoGlobal; HistoGlobal=NULL;}
    if(FitGlobal){delete FitGlobal; FitGlobal=NULL;}
    if(GlobalToMomentum){delete[]GlobalToMomentum; GlobalToMomentum=NULL;}
    if(NumBinsSyst){delete[]NumBinsSyst;NumBinsSyst=NULL;}
    if(CumulativeNumBinsSyst){delete[]CumulativeNumBinsSyst;CumulativeNumBinsSyst=NULL;}
//printf("~SourceSystems=%p\n",SourceSystems);
    if(SourceSystems) {delete[]SourceSystems; SourceSystems=NULL;}
    if(ParentSource) {delete[]ParentSource; ParentSource=NULL;}
    if(ParentParameter) {delete[]ParentParameter; ParentParameter=NULL;}
    if(PotentialSystems) {delete[]PotentialSystems; PotentialSystems=NULL;}
    if(ParentPotential) {delete[]ParentPotential; ParentPotential=NULL;}
}

//void DLM_Fitter1::TEST1(const unsigned& WhichSyst, TH1F* histo, const double& FromMeV ,
//                   DLM_CkDecomposition* decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX){
//HistoToFit[0] = new TH1F(TString::Format("HistoToFit%u",0),TString::Format("HistoToFit%u",0),10,0,1);
//return;
//}

void DLM_Fitter1::SetSystem(const unsigned& WhichSyst, const TH1F& histo, const double& FromMeV ,
                   DLM_CkDecomposition& decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX){


//printf("&decomp = %p\n",&decomp);

    if(WhichSyst>=MaxNumSyst){
        printf("WARNING: SetSystem says WhichSyst>=MaxNumSyst\n");
        return;
    }
    if(KMIN>KFEMTO || KMIN>KLINEAR || KMIN>KMAX ||
       KFEMTO>KLINEAR || KFEMTO>KMAX || KLINEAR>KMAX){
        printf("WARNING: SetSystem says that the input ranges make no sense!\n");
printf("KMIN=%.2f; KFEMTO=%.2f; KLINEAR=%.2f; KMAX=%.2f\n",KMIN,KFEMTO,KLINEAR,KMAX);
        return;
    }
    if(KMIN<histo.GetBinLowEdge(1)/FromMeV || KMAX>histo.GetXaxis()->GetBinUpEdge(histo.GetNbinsX())/FromMeV){
//printf("KMIN=%f --> %f\n",KMIN,histo.GetBinLowEdge(1)/FromMeV);
//printf(" KMAX=%f --> %f\n",KMAX,histo.GetXaxis()->GetBinUpEdge(histo.GetNbinsX())/FromMeV);
        printf("WARNING: SetSystem says that the desired fitting region is outside of the histo range!\n");
        return;
    }

    unsigned BinMin = histo.GetXaxis()->FindBin(KMIN/FromMeV);
    //unsigned BinFemto = histo->FindBin(KFEMTO/FromMeV);
    //unsigned BinLinear = histo->FindBin(KLINEAR/FromMeV);
    unsigned BinMax = histo.GetXaxis()->FindBin(KMAX/FromMeV);

    FitRange[WhichSyst][kmin] = KMIN;
    FitRange[WhichSyst][kf] = KFEMTO;
    FitRange[WhichSyst][kl] = KLINEAR;
    FitRange[WhichSyst][kmax] = KMAX;

    unsigned NumBins = BinMax-BinMin+1;
    double* Bins = new double [NumBins+1];
    for(unsigned uBin=BinMin; uBin<=BinMax; uBin++){
        Bins[uBin-BinMin] = histo.GetBinLowEdge(uBin);
    }
    Bins[NumBins] = histo.GetXaxis()->GetBinUpEdge(BinMax);
    //double kFrom = histo->GetBinLowEdge(BinMin);
    //double kTo = histo->GetXaxis()->GetBinUpEdge(BinMax);

    HistoOriginal[WhichSyst] = &histo;
    SystemToFit[WhichSyst] = &decomp;
//printf("HistoToFit[WhichSyst] = %p\n",HistoToFit[WhichSyst]);
    if(HistoToFit[WhichSyst]){
        delete HistoToFit[WhichSyst];
    }

//printf("HistoToFit = %p\n",HistoToFit);
//printf(" HistoToFit[%u] = %p\n",WhichSyst,HistoToFit[WhichSyst]);
//delete [] HistoToFit;
//HistoToFit = new TH1F* [10];
//TH1F* h1 = new TH1F(TString::Format("HistoToFit%u",WhichSyst),TString::Format("HistoToFit%u",WhichSyst),NumBins,Bins);
    HistoToFit[WhichSyst] = new TH1F(TString::Format("HistoToFit%u_%p",WhichSyst,this),TString::Format("HistoToFit%u_%p",WhichSyst,this),NumBins,Bins);
    //HistoToFit[WhichSyst] = new TH1F(TString::Format("HistoToFit%u",WhichSyst),TString::Format("HistoToFit%u",WhichSyst),10,0,1);

//HistoToFit[WhichSyst] = h1;

    for(unsigned uBin=BinMin; uBin<=BinMax; uBin++){
        HistoToFit[WhichSyst]->SetBinContent(uBin-BinMin+1, histo.GetBinContent(uBin));
        HistoToFit[WhichSyst]->SetBinError(uBin-BinMin+1, histo.GetBinError(uBin));
//if(WhichSyst==2){
//printf("uBin=%u; val=%f\n",uBin,histo.GetBinContent(uBin));
//}
    }
//printf("\n");
    delete [] Bins;
}
int DLM_Fitter1::GetSystem(const TString& System){
    char* buffer = new char [128];
    int SYSTID=-1;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        SystemToFit[uSyst]->GetName(buffer);
        if(strcmp(System.Data(),buffer)==0){
            SYSTID=uSyst;
        }
    }
    if(SYSTID==-1){
        printf("\033[1;33mWARNING:\033[0m Non-existing system %s\n", buffer);
    }
    return SYSTID;
    delete [] buffer;
}
TString DLM_Fitter1::GetSystem(const unsigned& WhichSyst){
    if(WhichSyst>=MaxNumSyst || !SystemToFit[WhichSyst]){
        printf("\033[1;33mWARNING:\033[0m Non-existing system %u\n", WhichSyst);
        return "";
    }
    char* buffer = new char [128];
    SystemToFit[WhichSyst]->GetName(buffer);
    TString Result = buffer;
    delete [] buffer;
    return Result;
}

bool DLM_Fitter1::ChangeCkModel(const unsigned& WhichSyst, DLM_CkDecomposition& decomp){
    if(WhichSyst>=MaxNumSyst){
        printf("WARNING: ChangeCkModel says WhichSyst>=MaxNumSyst\n");
        return false;
    }
    char nameNEW[64];
    decomp.GetName(nameNEW);
    char nameOLD[64];
    SystemToFit[WhichSyst]->GetName(nameOLD);
    if(strcmp(nameNEW,nameOLD)!=0){
        printf("WARNING: ChangeCkModel can only do its job for systems with the same name\n");
        return false;
    }
    SystemToFit[WhichSyst] = &decomp;
    return true;
}

void DLM_Fitter1::AddSameSource(const TString& System, const TString& EqualTo, const int& numpars){

    bool ExistingSystem = false;
    bool ExistingEqualTo = false;
    char* buffer = new char [128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            ExistingSystem += bool(SystemToFit[uSyst]->GetContribution(System.Data()));
            //we can only anchor stuff to main contributions! else it makes no sense
            SystemToFit[uSyst]->GetName(buffer);
            ExistingEqualTo += (strcmp(buffer, EqualTo)==0);
            //ExistingEqualTo += bool(SystemToFit[uSyst]->GetContribution(EqualTo.Data()));
        }
    }
    delete [] buffer;
    if(!ExistingSystem){
        printf("WARNING: AddSameSource says that %s does not exist!\n",System.Data());
        return;
    }
    if(!ExistingEqualTo){
        printf("WARNING: AddSameSource says that %s does not exist!\n",EqualTo.Data());
        return;
    }

    //we do not allow that a single entry is both a template for some other entry, but is
    //itself anchored to some other template entry.
    bool Trouble = false;
    for(unsigned uEntry=0; uEntry<NumSourceMapEntries; uEntry++){
        if(SameSourceMap[uEntry][0]==EqualTo) {Trouble=true; break;}
        //We also make sure that one system is not set equal to some other system more than once
        if(SameSourceMap[uEntry][0]==System) {Trouble=true; break;}
        if(SameSourceMap[uEntry][1]==System) {Trouble=true; break;}
    }
    if(Trouble){
        printf("WARNING: AddSameSource says that there is some Trouble!\n");
        return;
    }

    TString** NewMap = NULL;
    int* NewNumSameSourcePar = NULL;
    //if(SameSourceMap){
    NewMap = new TString*[NumSourceMapEntries+1];
    NewNumSameSourcePar = new int [NumSourceMapEntries+1];

    //copy the old info
    for(unsigned uEntry=0; uEntry<NumSourceMapEntries; uEntry++){
        NewMap[uEntry] = new TString [2];
        NewMap[uEntry][0] = SameSourceMap[uEntry][0];
        NewMap[uEntry][1] = SameSourceMap[uEntry][1];
        NewNumSameSourcePar[uEntry] = NumSameSourcePar[uEntry];
        delete [] SameSourceMap[uEntry];
    }
    if(SameSourceMap) {delete [] SameSourceMap;}
    if(NumSameSourcePar) {delete [] NumSameSourcePar;}
    //}
    NewMap[NumSourceMapEntries] = new TString [2];

    SameSourceMap = NewMap;
    NumSameSourcePar = NewNumSameSourcePar;

    SameSourceMap[NumSourceMapEntries][0] = System;
    SameSourceMap[NumSourceMapEntries][1] = EqualTo;
    NumSameSourcePar[NumSourceMapEntries] = numpars;

    NumSourceMapEntries++;
}

void DLM_Fitter1::RemoveSameSource(const TString& System){
    unsigned WhichEntryToRemove = 4294967295;
    for(unsigned uEntry=0; uEntry<NumSourceMapEntries; uEntry++){
        if(SameSourceMap[uEntry][0]==System) WhichEntryToRemove = uEntry;
    }
    if(WhichEntryToRemove!=4294967295){
        if(NumSourceMapEntries==1){
            delete [] SameSourceMap[0];
            delete [] SameSourceMap;
            SameSourceMap = NULL;
            delete [] NumSameSourcePar;
            NumSameSourcePar = NULL;
            NumSourceMapEntries = 0;
        }
        else{
            TString** NewMap = new TString*[NumSourceMapEntries-1];
            int* NewNumSameSourcePar = new int [NumSourceMapEntries-1];
            unsigned uEntryNew=0;
            for(unsigned uEntry=0; uEntry<NumSourceMapEntries; uEntry++){
                if(uEntry==WhichEntryToRemove) continue;
                NewMap[uEntryNew][0] = SameSourceMap[uEntry][0];
                NewMap[uEntryNew][1] = SameSourceMap[uEntry][1];
                NewNumSameSourcePar[uEntryNew] = NumSameSourcePar[uEntry];
                delete [] SameSourceMap[uEntry];
                uEntryNew++;
            }
            delete [] SameSourceMap;
            delete [] NumSameSourcePar;
            SameSourceMap = NewMap;
            NumSameSourcePar = NewNumSameSourcePar;
            NumSourceMapEntries--;
        }
    }
}


void DLM_Fitter1::AddSamePotential(const TString& System, const TString& EqualTo, const int& numpars){

    bool ExistingSystem = false;
    bool ExistingEqualTo = false;
    char* buffer = new char [128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            ExistingSystem += bool(SystemToFit[uSyst]->GetContribution(System.Data()));
            //we can only anchor stuff to main contributions! else it makes no sense
            SystemToFit[uSyst]->GetName(buffer);
            ExistingEqualTo += (strcmp(buffer, EqualTo)==0);
            //ExistingEqualTo += bool(SystemToFit[uSyst]->GetContribution(EqualTo.Data()));
        }
    }
    delete [] buffer;
    if(!ExistingSystem){
        printf("WARNING: AddSamePotential says that %s does not exist!\n",System.Data());
        return;
    }
    if(!ExistingEqualTo){
        printf("WARNING: AddSamePotential says that %s does not exist!\n",EqualTo.Data());
        return;
    }

    //we do not allow that a single entry is both a template for some other entry, but is
    //itself anchored to some other template entry.
    bool Trouble = false;
    for(unsigned uEntry=0; uEntry<NumPotentialMapEntries; uEntry++){
        if(SamePotentialMap[uEntry][0]==EqualTo) {Trouble=true; break;}
        //We also make sure that one system is not set equal to some other system more than once
        if(SamePotentialMap[uEntry][0]==System) {Trouble=true; break;}
        if(SamePotentialMap[uEntry][1]==System) {Trouble=true; break;}
    }
    if(Trouble){
        printf("WARNING: AddSamePotential says that there is some Trouble!\n");
        return;
    }

    TString** NewMap = NULL;
    int* NewNumSamePotentialPar = NULL;
    //if(SamePotentialMap){
    NewMap = new TString*[NumPotentialMapEntries+1];
    NewNumSamePotentialPar = new int [NumPotentialMapEntries+1];

    //copy the old info
    for(unsigned uEntry=0; uEntry<NumPotentialMapEntries; uEntry++){
        NewMap[uEntry] = new TString [2];
        NewMap[uEntry][0] = SamePotentialMap[uEntry][0];
        NewMap[uEntry][1] = SamePotentialMap[uEntry][1];
        NewNumSamePotentialPar[uEntry] = NumSamePotentialPar[uEntry];
        delete [] SamePotentialMap[uEntry];
    }
    if(SamePotentialMap) {delete [] SamePotentialMap;}
    if(NumSamePotentialPar) {delete [] NumSamePotentialPar;}
    //}
    NewMap[NumPotentialMapEntries] = new TString [2];

    SamePotentialMap = NewMap;
    NumSamePotentialPar = NewNumSamePotentialPar;

    SamePotentialMap[NumPotentialMapEntries][0] = System;
    SamePotentialMap[NumPotentialMapEntries][1] = EqualTo;
    NumSamePotentialPar[NumPotentialMapEntries] = numpars;

//printf("System=%s; EqualTo=%s\n",System.Data(),EqualTo.Data());

    NumPotentialMapEntries++;
}

void DLM_Fitter1::RemoveSamePotential(const TString& System){
    unsigned WhichEntryToRemove = 4294967295;
    for(unsigned uEntry=0; uEntry<NumPotentialMapEntries; uEntry++){
        if(SamePotentialMap[uEntry][0]==System) WhichEntryToRemove = uEntry;
    }
    if(WhichEntryToRemove!=4294967295){
        if(NumPotentialMapEntries==1){
            delete [] SamePotentialMap[0];
            delete [] SamePotentialMap;
            SamePotentialMap = NULL;
            delete [] NumSamePotentialPar;
            NumSamePotentialPar = NULL;
            NumPotentialMapEntries = 0;
        }
        else{
            TString** NewMap = new TString*[NumPotentialMapEntries-1];
            int* NewNumSamePotentialPar = new int [NumPotentialMapEntries-1];
            unsigned uEntryNew=0;
            for(unsigned uEntry=0; uEntry<NumPotentialMapEntries; uEntry++){
                if(uEntry==WhichEntryToRemove) continue;
                NewMap[uEntryNew][0] = SamePotentialMap[uEntry][0];
                NewMap[uEntryNew][1] = SamePotentialMap[uEntry][1];
                NewNumSamePotentialPar[uEntryNew] = NumSamePotentialPar[uEntry];
                delete [] SamePotentialMap[uEntry];
                uEntryNew++;
            }
            delete [] SamePotentialMap;
            delete [] NumSamePotentialPar;
            SamePotentialMap = NewMap;
            NumSamePotentialPar = NewNumSamePotentialPar;
            NumPotentialMapEntries--;
        }
    }
}

void DLM_Fitter1::AddSameParameter(const TString& System, const unsigned& WhichPar, const TString& ParentSystem, const unsigned& ParentPar){
    int WhichSyst = GetSystem(System);
    int WhichParent = GetSystem(ParentSystem);
    if(WhichSyst<0 || WhichParent<0){
        return;
    }
    if(WhichPar>=NumPar){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return;
    }
    if(ParentPar>=NumPar){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", ParentPar);
        return;
    }
    //this guy is -1 if the system has no parent, else it is WhichParent*NumPar+ParentPar (i.e. the id of the parent parameter)
    ParentParameter[WhichSyst*NumPar+WhichPar] = WhichParent*NumPar+ParentPar;

    if(GetBaseParameter(System,WhichPar)==-2){
        printf("\033[1;33mWARNING:\033[0m You are making a closed circle without a base parameter for %s par%u\n", System.Data(), WhichPar);
        ParentParameter[WhichSyst*NumPar+WhichPar] = -1;
    }
}
void DLM_Fitter1::RemoveSameParameter(const TString& System, const unsigned& WhichPar){
    int WhichSyst = GetSystem(System);
    if(WhichSyst<0){
        return;
    }
    if(WhichPar>=NumPar){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return;
    }
    ParentParameter[WhichSyst*NumPar+WhichPar] = -1;
}
int DLM_Fitter1::GetBaseParameter(const int& WhichSyst, const int& WhichPar){
    int sDUMMY;
    int uDUMMY;
    return GetBaseParameter(WhichSyst,WhichPar,sDUMMY,uDUMMY);
}
int DLM_Fitter1::GetBaseParameter(const int& WhichSyst, const int& WhichPar, int& ParentSystem, int& ParentPar){
    return GetBaseParameter(WhichSyst,WhichPar,ParentSystem,ParentPar,WhichSyst,WhichPar);
}
//return -2 is error
int DLM_Fitter1::GetBaseParameter(const int& WhichSyst, const int& WhichPar, int& ParentSystem, int& ParentPar, const int& StartSystem, const int& StartPar){
    if(WhichSyst<0){
        return -2;
    }
    if(WhichPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return -2;
    }
    if(ParentParameter[WhichSyst*NumPar+WhichPar]==-1){
        ParentSystem = WhichSyst;
        ParentPar = WhichPar;
        return WhichSyst*NumPar+WhichPar;
    }
    else{
        //unsigned NextSystem = (WhichSyst*NumPar+WhichPar)/WhichSyst;
        //unsigned NextPar = (WhichSyst*NumPar+WhichPar)%(WhichSyst*NumPar);
        //if(NextSystem==StartSystem && NextPar==StartPar) return false;
        //return GetBaseParameter(GetSystem((WhichSyst*NumPar+WhichPar)/WhichSyst),(WhichSyst*NumPar+WhichPar)%(WhichSyst*NumPar),ParentSystem,ParentPar,StartSystem,StartPar);
        int NextPar = ParentParameter[WhichSyst*NumPar+WhichPar];
        int NextSystem = NextPar/NumPar;
        //StartSystem StartPar is the original configuration, if we reach that one, we are in a loop
        if(NextSystem==StartSystem && NextPar==StartPar) return -2;
        //check the next parent
        return GetBaseParameter(NextSystem,NextPar,ParentSystem,ParentPar,StartSystem,StartPar);
    }
}


int DLM_Fitter1::GetBaseParameter(const TString& System, const int& WhichPar){
    int WhichSyst = GetSystem(System);
    if(WhichSyst<0){
        return -2;
    }
    if(WhichPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return -2;
    }
    return GetBaseParameter(WhichSyst,WhichPar);
}
int DLM_Fitter1::GetBaseParameter(const TString& System, const int& WhichPar, TString& ParentSystem, int& ParentPar){
    int WhichSyst = GetSystem(System);
    if(WhichSyst<0){
        return -2;
    }
    if(WhichPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return -2;
    }
    int WhichParent = GetSystem(ParentSystem);
    if(WhichParent<0){
        return -2;
    }
    if(ParentPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", ParentPar);
        return -2;
    }
    return GetBaseParameter(WhichSyst,WhichPar,WhichParent,ParentPar);
}
//figure out which is the parent parameter (recursively, i.e. if the parent has a parent we return that one)
int DLM_Fitter1::GetBaseParameter(const TString& System, const int& WhichPar, TString& ParentSystem, int& ParentPar, const TString& StartSystem, const int& StartPar){
    int WhichSyst = GetSystem(System);
    if(WhichSyst<0){
        return -2;
    }
    if(WhichPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", WhichPar);
        return -2;
    }
    int WhichParent = GetSystem(ParentSystem);
    if(WhichParent<0){
        return -2;
    }
    if(ParentPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", ParentPar);
        return -2;
    }
    int WhichStart = GetSystem(StartSystem);
    if(WhichStart<0){
        return -2;
    }
    if(StartPar>=int(NumPar)){
        printf("\033[1;33mWARNING:\033[0m Non-existing parameter %u\n", StartPar);
        return -2;
    }
    return GetBaseParameter(WhichSyst,WhichPar,WhichParent,ParentPar,WhichStart,StartPar);
}

void DLM_Fitter1::SetParameter(const unsigned& WhichSyst, const unsigned& WhichPar, const double& Value, const double& ValueDown, const double& ValueUp){
    if(WhichSyst>=MaxNumSyst) return;
    if(WhichPar>=NumPar) return;
    ParValue[WhichSyst][WhichPar] = Value;
    ParDownLimit[WhichSyst][WhichPar] = ValueDown;
    ParUpLimit[WhichSyst][WhichPar] = ValueUp;
    FixPar[WhichSyst][WhichPar] = false;
}
void DLM_Fitter1::SetParameter(const TString& WhichSyst, const unsigned& WhichPar, const double& Value, const double& ValueDown, const double& ValueUp){
    char* buffer = new char[128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            SystemToFit[uSyst]->GetName(buffer);
//printf("buffer=%s\n",buffer);
            if(strcmp(WhichSyst,buffer)==0){
                SetParameter(uSyst, WhichPar, Value, ValueDown, ValueUp);
                delete [] buffer;
                return;
            }
        }
    }
    printf("WARNING: SetParameter says that the system '%s' does not exist!\n", WhichSyst.Data());
    delete [] buffer;
}
void DLM_Fitter1::FixParameter(const unsigned& WhichSyst, const unsigned& WhichPar, const double& Value){
//printf("WhichSyst=%u; WhichPar=%u; Value=%f\n",WhichSyst,WhichPar,Value);
    if(WhichSyst>=MaxNumSyst) return;
    if(WhichPar>=NumPar) return;
    ParValue[WhichSyst][WhichPar] = Value;
    ParDownLimit[WhichSyst][WhichPar] = Value;
    ParUpLimit[WhichSyst][WhichPar] = Value;
    FixPar[WhichSyst][WhichPar] = true;
}
void DLM_Fitter1::FixParameter(const TString& WhichSyst, const unsigned& WhichPar, const double& Value){
    if(WhichPar>=NumPar) return;
    char* buffer = new char[128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(WhichSyst,buffer)==0){
                FixParameter(uSyst, WhichPar, Value);
                delete [] buffer;
                return;
            }
        }
    }
    printf("WARNING: FixParameter says that the system '%s' does not exist!\n", WhichSyst.Data());
    delete [] buffer;
}


void DLM_Fitter1::FixParameter(const unsigned& WhichSyst, const unsigned& WhichPar){
    FixParameter(WhichSyst,WhichPar,ParValue[WhichSyst][WhichPar]);
}
void DLM_Fitter1::FixParameter(const TString& WhichSyst, const unsigned& WhichPar){
    FixParameter(WhichSyst,WhichPar,GetParameter(WhichSyst,WhichPar));
}

double DLM_Fitter1::GetParameter(const unsigned& WhichSyst, const unsigned& WhichPar){
    if(WhichSyst>=MaxNumSyst) return 0;
    if(WhichPar>=NumPar) return 0;
    return FitGlobal->GetParameter(WhichSyst*NumPar+WhichPar);
}
double DLM_Fitter1::GetParError(const unsigned& WhichSyst, const unsigned& WhichPar){
    if(WhichSyst>=MaxNumSyst) return 0;
    if(WhichPar>=NumPar) return 0;
    if(SeparateBaseLineFit[WhichSyst] && HistoToFit[WhichSyst] && WhichPar<3){
        return FitBL[WhichSyst]->GetParError(WhichPar);
    }
    else{
        return FitGlobal->GetParError(WhichSyst*NumPar+WhichPar);
    }
}
double DLM_Fitter1::GetChi2(){
    return FitGlobal->GetChisquare();
}
int DLM_Fitter1::GetNdf(){
    return FitGlobal->GetNDF();
}
double DLM_Fitter1::GetChi2Ndf(){
    return GetChi2()/double(GetNdf());
}
double DLM_Fitter1::GetPval(){
    return FitGlobal->GetProb();
}
/*
double DLM_Fitter1::Eval(const unsigned& WhichSyst, const double& Momentum)(){
    if(WhichSyst>=MaxNumSyst) return 0;


    HistoOriginal[WhichSyst]->FindBin(Momentum);

    return FitGlobal->Eval();
}
*/
void DLM_Fitter1::GetFitGraph(const unsigned& WhichSyst, TGraph& OutGraph){
//printf(" HEY -> GetFitGraph!\n");
    if(WhichSyst>=MaxNumSyst) return;
//printf(" HEY2 -> GetFitGraph!\n");
    OutGraph.Set(NumBinsSyst[WhichSyst]);
    double Momentum;
    double xGlobal;
    unsigned GlobalBinShift=0;
    for(unsigned uSyst=0; uSyst<WhichSyst; uSyst++){
        GlobalBinShift += NumBinsSyst[uSyst];
    }
//printf("GlobalBinShift=%u\n",GlobalBinShift);
    for(unsigned uBin=0; uBin<NumBinsSyst[WhichSyst]; uBin++){
        Momentum = GlobalToMomentum[GlobalBinShift+uBin];
        xGlobal = HistoGlobal->GetBinCenter(GlobalBinShift+uBin+1);
//printf(" k=%.2f; xVal=%.2f\n",Momentum,xGlobal);
        OutGraph.SetPoint(uBin,Momentum,FitGlobal->Eval(xGlobal));
    }
}

double DLM_Fitter1::GetParameter(const TString& WhichSyst, const unsigned& WhichPar){
    char* buffer = new char[128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(WhichSyst,buffer)==0){
                delete [] buffer;
                return GetParameter(uSyst, WhichPar);
            }
        }
    }
    printf("WARNING: SetParameter says that the system '%s' does not exist!\n", WhichSyst.Data());
    delete [] buffer;
    return 0;
}
double DLM_Fitter1::GetParError(const TString& WhichSyst, const unsigned& WhichPar){
    char* buffer = new char[128];
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(WhichSyst,buffer)==0){
                delete [] buffer;
                return GetParError(uSyst, WhichPar);
            }
        }
    }
    printf("WARNING: SetParameter says that the system '%s' does not exist!\n", WhichSyst.Data());
    delete [] buffer;
    return 0;
}

void DLM_Fitter1::SetOutputDir(const TString& outdirname){
    OutputDirName = outdirname;
}
void DLM_Fitter1::SetSeparateBL(const unsigned& WhichSyst, const bool& yesno){
    if(WhichSyst>=MaxNumSyst) return;
    SeparateBaseLineFit[WhichSyst] = yesno;
}
void DLM_Fitter1::RemoveNegativeCk(const bool& yesno){
    RemoveNegCk = yesno;
}
bool DLM_Fitter1::CheckNegativeCk(){
    //printf("1/f0=%f\n",FitGlobal->GetParameter(p_pot0));
    //printf("d0=%f\n",FitGlobal->GetParameter(p_pot1));
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(SystemToFit[uSyst]){
            for(unsigned uBin=0; uBin<SystemToFit[uSyst]->GetCk()->GetNbins(); uBin++){
                //printf(" SystemToFit[uSyst]->GetCk()->GetBinContent(%u)=%f\n",uBin,SystemToFit[uSyst]->GetCk()->GetBinContent(uBin));
                if(SystemToFit[uSyst]->GetCk()->GetBinContent(uBin)<0) return true;
            }
        }
    }
    return false;
}

TF1* DLM_Fitter1::GetFit(){
    return FitGlobal;
}
TF1* DLM_Fitter1::GetBaselineFit(const unsigned& WhichSyst) const{
    if(WhichSyst>=MaxNumSyst) return NULL;
    return FitBL[WhichSyst];
}
const TH1F* DLM_Fitter1::GetGlobalHisto() const{
    return HistoGlobal;
}

void DLM_Fitter1::GoBabyGo(){

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]){
            printf("ERROR: At the moment running GoBabyGo is not possible if not all systems are defined!\n");
            return;
        }
    }

    //the separate fit thing
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SeparateBaseLineFit[uSyst] || !HistoToFit[uSyst]) continue;
        if(!HistoToFit[uSyst]) continue;

        if(FitBL[uSyst]) delete FitBL[uSyst];
        FitBL[uSyst] = new TF1(TString::Format("FitBL%u_%p",uSyst,this),"[0]+[1]*x+[2]*x*x",FitRange[uSyst][kl],FitRange[uSyst][kmax]);
        if(FixPar[uSyst][p_a]) {FitBL[uSyst]->FixParameter(0, ParValue[uSyst][p_a]);}
        else{FitBL[uSyst]->SetParameter(0, ParValue[uSyst][p_a]); FitBL[uSyst]->SetParLimits(0, ParDownLimit[uSyst][p_a], ParUpLimit[uSyst][p_a]);}
        if(FixPar[uSyst][p_b]) {FitBL[uSyst]->FixParameter(1, ParValue[uSyst][p_b]);}
        else{FitBL[uSyst]->SetParameter(1, ParValue[uSyst][p_b]); FitBL[uSyst]->SetParLimits(1, ParDownLimit[uSyst][p_b], ParUpLimit[uSyst][p_b]);}
        if(FixPar[uSyst][p_c]) {FitBL[uSyst]->FixParameter(2, ParValue[uSyst][p_c]);}
        else{FitBL[uSyst]->SetParameter(2, ParValue[uSyst][p_c]); FitBL[uSyst]->SetParLimits(2, ParDownLimit[uSyst][p_c], ParUpLimit[uSyst][p_c]);}
        //HistoToFit[uSyst]->Fit(FitBL[uSyst],"S, N, R, M");
        HistoToFit[uSyst]->Fit(FitBL[uSyst],"Q, S, N, R, M");
        FixParameter(uSyst,p_a,FitBL[uSyst]->GetParameter(0));
        FixParameter(uSyst,p_b,FitBL[uSyst]->GetParameter(1));
        FixParameter(uSyst,p_c,FitBL[uSyst]->GetParameter(2));
        FixParameter(uSyst,p_kc,0);

        /*
        TF1* FitBLold = new TF1(TString::Format("FitBLold%u",uSyst),"[0]+[1]*x+[2]*x*x",FitRange[uSyst][kl],FitRange[uSyst][kmax]);
        if(FixPar[uSyst][p_a]) {FitBLold->FixParameter(0, ParValue[uSyst][p_a]);}
        else{FitBLold->SetParameter(0, ParValue[uSyst][p_a]); FitBLold->SetParLimits(0, ParDownLimit[uSyst][p_a], ParUpLimit[uSyst][p_a]);}
        if(FixPar[uSyst][p_b]) {FitBLold->FixParameter(1, ParValue[uSyst][p_b]);}
        else{FitBLold->SetParameter(1, ParValue[uSyst][p_b]); FitBLold->SetParLimits(1, ParDownLimit[uSyst][p_b], ParUpLimit[uSyst][p_b]);}
        if(FixPar[uSyst][p_c]) {FitBLold->FixParameter(2, ParValue[uSyst][p_c]);}
        else{FitBLold->SetParameter(2, ParValue[uSyst][p_c]); FitBLold->SetParLimits(2, ParDownLimit[uSyst][p_c], ParUpLimit[uSyst][p_c]);}
        HistoToFit[uSyst]->Fit(FitBLold,"S, N, R, M");
        //HistoToFit[uSyst]->Fit(FitBLold,"Q, S, N, R, M");
        FixParameter(uSyst,p_a,FitBLold->GetParameter(0));
        FixParameter(uSyst,p_b,FitBLold->GetParameter(1));
        FixParameter(uSyst,p_c,FitBLold->GetParameter(2));
        FixParameter(uSyst,p_kc,0);
        delete FitBLold;
        */
    }

    double Momentum;
    double CkVal;
    double CkValErr;
    unsigned NumOfBins;

    //this might be needed if I allow for zero pointers to some of the main functions
    //for the time being this is not the case
    unsigned ActualNumSyst = 0;

    char* buffer = new char[128];
    char* buffer2 = new char[128];
//----- NEXT WE HANDLE THE POSSIBILITY TO SHARE SOURCES BETWEEN SYSTEMS -----
    //next we count how many unique source we have within our systems
    NumSourceSystems=0;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        NumSourceSystems += bool(SystemToFit[uSyst]);
    }
    for(unsigned uSourceMap=0; uSourceMap<NumSourceMapEntries; uSourceMap++){
        bool IsNotCkToFit=true;
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SameSourceMap[uSourceMap][0].Data())==0) IsNotCkToFit=false;
        }
        //we add to the counter all entries with a parent source, that are not main contributions.
        //the reason is that the main contributions are already included above.
        NumSourceSystems += IsNotCkToFit;
    }

//printf("NumSourceSystems=%u\n",NumSourceSystems);
    //we reset all memory regarding the fit
    if(SourceSystems){delete[]SourceSystems;SourceSystems=NULL;}
    if(NumSourceSystems) SourceSystems = new DLM_CkDecomposition*[NumSourceSystems];
    if(ParentSource){delete[]ParentSource;ParentSource=NULL;}
    if(NumSourceSystems) ParentSource = new int [NumSourceSystems];
//printf("SourceSystems=%p\n",SourceSystems);

    unsigned uActualSource=0;
    //set the addresses of the main Ck functions to the SourceSystems
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;
        SourceSystems[uActualSource] = SystemToFit[uSyst];
        //by default the parent of the system is the system itself
        ParentSource[uActualSource]=uSyst;
        SystemToFit[uSyst]->GetName(buffer);

        /*
        for(unsigned uSyst2=0; uSyst2<MaxNumSyst; uSyst2++){
            if(!SystemToFit[uSyst2]) continue;
            SystemToFit[uSyst2]->GetName(buffer2);
            //if(strcmp(buffer,SameSourceMap[uActualSource][0].Data())==0) ParentSource[uActualSource]=uSyst2;
            //if(strcmp(buffer,SameSourceMap[uActualSource][1].Data())==0) ParentSource[uActualSource]=uSyst2;
            if(strcmp(buffer,buffer2)==0) ParentSource[uActualSource]=uSyst2;
        }
        */
        //find the source-parent if there is any
        for(unsigned uSyst2=0; uSyst2<MaxNumSyst; uSyst2++){
            if(!SystemToFit[uSyst2]) continue;
            SystemToFit[uSyst2]->GetName(buffer2);
            //in case we find in the map a relation between uSyst and uSyst2 -> we set this
            for(unsigned uSourceMap=0; uSourceMap<NumSourceMapEntries; uSourceMap++){
                if(strcmp(buffer,SameSourceMap[uSourceMap][0].Data())==0 &&
                   strcmp(buffer2,SameSourceMap[uSourceMap][1].Data())==0) ParentSource[uActualSource]=uSyst2;
            }
        }

//SourceSystems[uActualSource]->GetName(buffer);
//printf("SourceSystems[%u] name: %s\n",uActualSource,buffer);
//SystemToFit[ParentSource[uActualSource]]->GetName(buffer);
//printf("     Parent name: %s\n",buffer);

        uActualSource++;
    }

    //set the addresses of the secondary Ck functions to the SourceSystems
    for(unsigned uSourceMap=0; uSourceMap<NumSourceMapEntries; uSourceMap++){
        bool IsNotCkToFit=true;
        DLM_CkDecomposition* SystAddress=NULL;
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SameSourceMap[uSourceMap][0].Data())==0) IsNotCkToFit=false;
            SystAddress = SystemToFit[uSyst]->GetContribution(SameSourceMap[uSourceMap][0].Data());
            //make sure this entry is not yet saved
            for(unsigned uAS=0; uAS<uActualSource; uAS++){
                if(SourceSystems[uAS]==SystAddress){
                    SystAddress=NULL;
                    break;
                }
            }
            if(SystAddress) break;
        }
        //this guy lets only contributions from the map which are not main Ck to pass through
        if(!IsNotCkToFit) continue;
        SourceSystems[uActualSource] = SystAddress;

        //find the parent
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SameSourceMap[uSourceMap][1].Data())==0) ParentSource[uActualSource]=uSyst;
        }

//SystAddress->GetName(buffer);
//printf("SystAddress name: %s\n",buffer);
//SystemToFit[ParentSource[uActualSource]]->GetName(buffer);
//printf("     Parent name: %s\n",buffer);

        uActualSource++;
    }

//printf("QA: %u=%u ?\n",NumSourceSystems,uActualSource);


//----- NEXT WE HANDLE THE POSSIBILITY TO SHARE POTENTIALS BETWEEN SYSTEMS -----
    //next we count how many unique source we have within our systems
    NumPotentialSystems=0;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        NumPotentialSystems += bool(SystemToFit[uSyst]);
    }
    for(unsigned uPotentialMap=0; uPotentialMap<NumPotentialMapEntries; uPotentialMap++){
//printf("uPotentialMap=%u; NumPotentialMapEntries=%u\n" ,uPotentialMap, NumPotentialMapEntries);
        bool IsNotCkToFit=true;
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SamePotentialMap[uPotentialMap][0].Data())==0) IsNotCkToFit=false;
//printf("buffer = %s;\n",buffer);
//printf(" SamePotentialMap = %p;\n",SamePotentialMap);
//printf(" SamePotentialMap[%u] = %p;\n",uPotentialMap,SamePotentialMap[uPotentialMap]);
//printf(" SamePotentialMap[%u][0] = %s;\n",uPotentialMap,SamePotentialMap[uPotentialMap][0].Data());


//printf("buffer = %s; SamePotentialMap[%u][0] = %s\n",buffer,SamePotentialMap[uPotentialMap][0].Data());

        }
        //we add to the counter all entries with a parent source, that are not main contributions.
        //the reason is that the main contributions are already included above.
//printf(" uPotentialMap=%u; IsNotCkToFit=%i\n" ,uPotentialMap, int(IsNotCkToFit));
        NumPotentialSystems += IsNotCkToFit;
    }

//printf("NumPotentialSystems=%u\n",NumPotentialSystems);
    //we reset all memory regarding the fit
    if(PotentialSystems){delete[]PotentialSystems;PotentialSystems=NULL;}
    if(NumPotentialSystems) PotentialSystems = new DLM_CkDecomposition*[NumPotentialSystems];
    if(ParentPotential){delete[]ParentPotential;ParentPotential=NULL;}
    if(NumPotentialSystems) ParentPotential = new int [NumPotentialSystems];
//printf("PotentialSystems=%p\n",PotentialSystems);
//printf("NumPotentialSystems=%u\n",NumPotentialSystems);

    unsigned uActualPotential=0;
    //set the addresses of the main Ck functions to the PotentialSystems
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;
        PotentialSystems[uActualPotential] = SystemToFit[uSyst];
        //by default the parent of the system is the system itself
        ParentPotential[uActualPotential]=uSyst;
        SystemToFit[uSyst]->GetName(buffer);

        /*
        for(unsigned uSyst2=0; uSyst2<MaxNumSyst; uSyst2++){
            if(!SystemToFit[uSyst2]) continue;
            SystemToFit[uSyst2]->GetName(buffer2);
            //if(strcmp(buffer,SamePotentialMap[uActualPotential][0].Data())==0) ParentPotential[uActualPotential]=uSyst2;
            //if(strcmp(buffer,SamePotentialMap[uActualPotential][1].Data())==0) ParentPotential[uActualPotential]=uSyst2;
            if(strcmp(buffer,buffer2)==0) ParentPotential[uActualPotential]=uSyst2;
        }
        */
        //find the source-parent if there is any
        for(unsigned uSyst2=0; uSyst2<MaxNumSyst; uSyst2++){
            if(!SystemToFit[uSyst2]) continue;
            SystemToFit[uSyst2]->GetName(buffer2);
            //in case we find in the map a relation between uSyst and uSyst2 -> we set this
            for(unsigned uPotentialMap=0; uPotentialMap<NumPotentialMapEntries; uPotentialMap++){
                if(strcmp(buffer,SamePotentialMap[uPotentialMap][0].Data())==0 &&
                   strcmp(buffer2,SamePotentialMap[uPotentialMap][1].Data())==0) ParentPotential[uActualPotential]=uSyst2;
            }
        }

//PotentialSystems[uActualPotential]->GetName(buffer);
//printf("PotentialSystems[%u] name: %s\n",uActualPotential,buffer);
//SystemToFit[ParentPotential[uActualPotential]]->GetName(buffer);
//printf("     Parent name: %s\n",buffer);

        uActualPotential++;
    }

    //set the addresses of the secondary Ck functions to the PotentialSystems
    for(unsigned uPotentialMap=0; uPotentialMap<NumPotentialMapEntries; uPotentialMap++){
        bool IsNotCkToFit=true;
        DLM_CkDecomposition* SystAddress=NULL;
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SamePotentialMap[uPotentialMap][0].Data())==0) IsNotCkToFit=false;
            SystAddress = SystemToFit[uSyst]->GetContribution(SamePotentialMap[uPotentialMap][0].Data());
            //make sure this entry is not yet saved
            for(unsigned uAS=0; uAS<uActualPotential; uAS++){
                if(PotentialSystems[uAS]==SystAddress){
                    SystAddress=NULL;
                    break;
                }
            }
            if(SystAddress) break;
        }
        //this guy lets only contributions from the map which are not main Ck to pass through
        if(!IsNotCkToFit) continue;
        PotentialSystems[uActualPotential] = SystAddress;

        //find the parent
        for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
            if(!SystemToFit[uSyst]) continue;
            SystemToFit[uSyst]->GetName(buffer);
            if(strcmp(buffer,SamePotentialMap[uPotentialMap][1].Data())==0) ParentPotential[uActualPotential]=uSyst;
        }

//SystAddress->GetName(buffer);
//printf("SystAddress name: %s\n",buffer);
//SystemToFit[ParentPotential[uActualPotential]]->GetName(buffer);
//printf("     Parent name: %s\n",buffer);

        uActualPotential++;
    }

    delete [] buffer;
    delete [] buffer2;

    unsigned NumBinsGlobal = 0;
    if(NumBinsSyst) delete [] NumBinsSyst;
    NumBinsSyst = new unsigned [MaxNumSyst];
    if(CumulativeNumBinsSyst) delete [] CumulativeNumBinsSyst;
    CumulativeNumBinsSyst = new unsigned [MaxNumSyst];

    ActualNumSyst = 0;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        NumBinsSyst[uSyst]=0;
        CumulativeNumBinsSyst[uSyst]=(uSyst==0)?0:CumulativeNumBinsSyst[uSyst-1];
        if(!HistoToFit[uSyst]) continue;
        NumOfBins = HistoToFit[uSyst]->GetNbinsX();
        bool TakeThisBin;
        for(unsigned uBin=0; uBin<NumOfBins; uBin++){
            TakeThisBin = true;
            Momentum = HistoToFit[uSyst]->GetBinCenter(uBin+1);
            if(Momentum<FitRange[uSyst][kmin]) TakeThisBin=false;
            else if(Momentum>FitRange[uSyst][kf] && Momentum<FitRange[uSyst][kl]) TakeThisBin=false;
            else if(Momentum>FitRange[uSyst][kmax]) TakeThisBin=false;
            else if(SeparateBaseLineFit[uSyst] && Momentum>FitRange[uSyst][kf]) TakeThisBin=false;
            if(!TakeThisBin) continue;
            NumBinsGlobal++;
            NumBinsSyst[uSyst]++;
            CumulativeNumBinsSyst[uSyst]++;
        }
        ActualNumSyst++;
    }

    //double* Bin_kVal = new double[NumBinsGlobal];
    if(GlobalToMomentum) delete[]GlobalToMomentum;
    GlobalToMomentum = new double [NumBinsGlobal];


    if(HistoGlobal) delete HistoGlobal;
    HistoGlobal = new TH1F(TString::Format("HistoGlobal%p",this), TString::Format("HistoGlobal%p",this), NumBinsGlobal, 0, 1);
    if(FitGlobal) delete FitGlobal;
    FitGlobal = new TF1(TString::Format("FitGlobal%p",this), this, &DLM_Fitter1::EvalGlobal,0,1,ActualNumSyst*NumPar,"DLM_Fitter1","EvalGlobal");

    unsigned uBinGlob=0;
    //unsigned uActSyst=0;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
//printf("NumBinsSyst[%u]=%u (%u)\n",uSyst,NumBinsSyst[uSyst],CumulativeNumBinsSyst[uSyst]);
        if(!HistoToFit[uSyst]) continue;
        NumOfBins = HistoToFit[uSyst]->GetNbinsX();
        bool TakeThisBin;
        for(unsigned uBin=0; uBin<NumOfBins; uBin++){
            TakeThisBin = true;
            Momentum = HistoToFit[uSyst]->GetBinCenter(uBin+1);
            //exclude bins below the fit range
            if(Momentum<FitRange[uSyst][kmin]) TakeThisBin=false;
            //exclude bins in the dead region
            else if(Momentum>FitRange[uSyst][kf] && Momentum<FitRange[uSyst][kl]) TakeThisBin=false;
            //exclude bins above the fit range
            else if(Momentum>FitRange[uSyst][kmax]) TakeThisBin=false;
            //exclude bins above the femto region, in case we do perform a separate fit in the high-k region
            else if(SeparateBaseLineFit[uSyst] && Momentum>FitRange[uSyst][kf]) TakeThisBin=false;

            if(!TakeThisBin) continue;
            GlobalToMomentum[uBinGlob] = Momentum;
//printf(" GlobalToMomentum[%u]=%f\n",uBinGlob,GlobalToMomentum[uBinGlob]);
            CkVal = HistoToFit[uSyst]->GetBinContent(uBin+1);
            CkValErr = HistoToFit[uSyst]->GetBinError(uBin+1);
            HistoGlobal->SetBinContent(uBinGlob+1,CkVal);
            HistoGlobal->SetBinError(uBinGlob+1,CkValErr);
            uBinGlob++;
        }
        for(unsigned uPar=0; uPar<NumPar; uPar++){

            if(FixPar[uSyst][uPar]){
//printf("Fix PAR%u = %f\n", uSyst*NumPar+uPar, ParValue[uSyst][uPar]);
                FitGlobal->FixParameter(uSyst*NumPar+uPar, ParValue[uSyst][uPar]);
            }
            else{
//printf("Set PAR%u = %f <%f __ %f>\n", uSyst*NumPar+uPar, ParValue[uSyst][uPar], ParDownLimit[uSyst][uPar], ParUpLimit[uSyst][uPar]);
                FitGlobal->SetParameter(uSyst*NumPar+uPar, ParValue[uSyst][uPar]);
                FitGlobal->SetParLimits(uSyst*NumPar+uPar, ParDownLimit[uSyst][uPar], ParUpLimit[uSyst][uPar]);
            }

            //if the radius should be taken from a parent -> fix the source to a dummy value
            if(int(uPar)>=p_sor0 && int(uPar)<=p_sor5 && ParentSource[uSyst]!=int(uSyst)){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_sor1, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_sor2, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_sor3, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_sor4, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_sor5, -1e6);
//printf("FIX SOR %u\n",uSyst*NumPar+uPar);
            }
//printf("uSyst=%u; ParentPotential[uSyst]=%u\n",uSyst,ParentPotential[uSyst]);
            //if the potential should be taken from a parent -> fix the potential pars to dummy values
            if(int(uPar)>=p_pot0 && int(uPar)<=p_pot3 && ParentPotential[uSyst]!=int(uSyst)){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_pot1, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_pot2, -1e6);
                //FitGlobal->FixParameter(uSyst*NumPar+p_pot3, -1e6);
//printf("FIX POT\n");
//printf(" fixing %u\n",uSyst*NumPar+uPar);
//printf(" fixing %u %u %u %u\n",uSyst*NumPar+uPar,uSyst*NumPar+p_pot1,uSyst*NumPar+p_pot2,uSyst*NumPar+p_pot3);
            }

            if(ParentParameter[uSyst*NumPar+uPar]!=-1){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
            }

            //if kf==kl, than there is an active boundary condition for the linear part of C(k)
            //thus Cl will be determined by the fitter, here we set a dummy value
            if(int(uPar)==p_Cl && FitRange[uSyst][kl]==FitRange[uSyst][kf]){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
//printf("Fix PAR%u = %f\n", uSyst*NumPar+uPar, ParValue[uSyst][uPar]);
            }

            //if we want to fit with a flat long-range C(k)
            if(int(uPar)==p_kc && ParValue[uSyst][p_Cl]<=0){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
//printf("Fix PAR%u = %f\n", uSyst*NumPar+uPar, ParValue[uSyst][uPar]);
//printf("HEY from %u\n",uSyst*NumPar+uPar);
            }
        }

        //if there is no long-range fit, fix the variable that are modeling this part
        if(FitRange[uSyst][kl]==FitRange[uSyst][kmax]){
            FitGlobal->FixParameter(uSyst*NumPar+p_Cl, -1e6);
            FitGlobal->FixParameter(uSyst*NumPar+p_kc, -1e6);
//printf("Fix PAR%u = %f\n", uSyst*NumPar+p_Cl, ParValue[uSyst][p_Cl]);
//printf("Fix PAR%u = %f\n", uSyst*NumPar+p_kc, ParValue[uSyst][p_kc]);
        }
        //uActSyst++;
    }

    //HistoGlobal->Fit(FitGlobal,"V, S, N, R, M");
    HistoGlobal->Fit(FitGlobal,"S, N, R, M");
    //HistoGlobal->Fit(FitGlobal,"Q, S, N, R, M");
    //HistoGlobal->Fit(FitGlobal,"Q, N, R, M");

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        for(unsigned uPar=0; uPar<NumPar; uPar++){
            ParValue[uSyst][uPar] = FitGlobal->GetParameter(uSyst*NumPar+uPar);
        }
    }

/*
if(MaxNumSyst==1){
for(unsigned uBin=1; uBin<=5; uBin++){
printf("HG[%u] = %.2f --> Fit = %.2f --> CkDec(%p) = %.2f\n",
       uBin,
       FitGlobal->Eval(HistoOriginal[0]->GetBinCenter(uBin)),
       &SystemToFit[0][0],
       SystemToFit[0]->EvalCk(HistoOriginal[0]->GetBinCenter(uBin)));
printf("XX(%.2f): CkDec(%p) = %.2f\n",
       HistoOriginal[0]->GetBinCenter(uBin),
       &SystemToFit[0][0],
       SystemToFit[0]->EvalCk(HistoOriginal[0]->GetBinCenter(uBin)));
printf("\n");
}
}
if(MaxNumSyst==4){
for(unsigned uBin=1; uBin<=5; uBin++){
printf("pp(%.2f): CkDec(%p) = %.2f\n",
       HistoOriginal[0]->GetBinCenter(uBin),
       &SystemToFit[0][0],
       SystemToFit[0]->EvalCk(HistoOriginal[0]->GetBinCenter(uBin)));
printf("pΛ(%.2f): CkDec(%p) = %.2f\n",
       HistoOriginal[1]->GetBinCenter(uBin),
       &SystemToFit[1][0],
       SystemToFit[1]->EvalCk(HistoOriginal[1]->GetBinCenter(uBin)));
printf("ΛΛ(%.2f): CkDec(%p) = %.2f\n",
       HistoOriginal[2]->GetBinCenter(uBin),
       &SystemToFit[2][0],
       SystemToFit[2]->EvalCk(HistoOriginal[2]->GetBinCenter(uBin)));
printf("pΞ(%.2f): CkDec(%p) = %.2f\n",
       HistoOriginal[3]->GetBinCenter(uBin),
       &SystemToFit[3][0],
       SystemToFit[3]->EvalCk(HistoOriginal[3]->GetBinCenter(uBin)));
printf("\n");
}
}
*/



    //TFile* tmpFile = new TFile(TString::Format("%stmpFile.root",OutputDirName.Data()),"recreate");
    //HistoGlobal->Write();
    //FitGlobal->Write();
    //delete tmpFile;
//if(FitGlobal->GetNDF()==22){

//printf("NGLOB=%u\n",HistoGlobal->GetNbinsX());
//printf("NFIT=%u\n",FitGlobal->GetNDF());
//printf("NFPar=%u\n",FitGlobal->GetNumberFreeParameters());
//printf("NFPts=%u\n\n",FitGlobal->GetNumberFitPoints());

//FitGlobal->Print();
//}

}

double DLM_Fitter1::EvalGlobal(double* xVal, double* Pars){
//printf("CALLING EVALGLOBAL\n");
    unsigned GlobalBin = HistoGlobal->FindBin(*xVal)-1;
    double Momentum = GlobalToMomentum[GlobalBin];
    unsigned WhichSyst=0;
//printf("Pars0[p_sor0/1/2] = %.2f/%.2f/%.2f\n",Pars[p_sor0],Pars[p_sor1],Pars[p_sor2]);
//printf("Pars1[p_sor0/1/2] = %.2f/%.2f/%.2f\n",Pars[NumPar+p_sor0],Pars[NumPar+p_sor1],Pars[NumPar+p_sor2]);

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!HistoToFit[uSyst]) continue;
        //determine in which system the GlobalBin is
        if(GlobalBin<CumulativeNumBinsSyst[uSyst]) {WhichSyst=uSyst; break;}
    }

    if(!SystemToFit[WhichSyst]) return 0;

    //we make sure ALL Pars are set to a meaningful value (i.e. in case there is a parent, we substitute the value)
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        for(unsigned uPar=0; uPar<NumPar; uPar++){
            int uParentParameter = GetBaseParameter(uSyst,uPar);
            Pars[uSyst*NumPar+uPar]=Pars[uParentParameter];
        }
    }
/*
    for(unsigned uPar=0; uPar<NumPar; uPar++){
        //if(ParentParameter[WhichSyst*NumPar+uPar]!=-1){
        int WhichParentParameter = GetBaseParameter(WhichSyst,uPar);
//printf("hi\n");
//if(WhichSyst==1){
//printf("WhichSyst=%u; WhichSyst*NumPar+uPar=%u; WhichParentParameter=%i; Pars[%u]=%f\n",
//       WhichSyst,WhichSyst*NumPar+uPar,WhichParentParameter,WhichSyst*NumPar+uPar,Pars[WhichSyst*NumPar+uPar]);
//}
        if(WhichParentParameter!=int(WhichSyst*NumPar+uPar)){
            Pars[WhichSyst*NumPar+uPar]=Pars[WhichParentParameter];
//printf("Pars[%u]=Pars[%i]=%f\n",WhichSyst*NumPar+uPar,WhichParentParameter,Pars[WhichParentParameter]);
        }
    }
*/


    //determine if we are in the femto region. If false => we are in the linear or flat region
    bool FemtoRegion = (Momentum<=FitRange[WhichSyst][kf]);
    //loop over all systems that may need to have their source size adjusted
    //adjust the source of their corresponding "source-parent"
//printf("WTF man\n");
//usleep(0.1e6);
    for(unsigned uSource=0; uSource<NumSourceSystems; uSource++){
        for(unsigned uPar=0; uPar<SourceSystems[uSource]->GetCk()->GetNumSourcePar(); uPar++){
//char buffer[32];
//SourceSystems[uSource]->GetName(buffer);
//printf("ssName=%s\n",buffer);
            SourceSystems[uSource]->GetCk()->SetSourcePar(uPar, Pars[ParentSource[uSource]*NumPar+p_sor0+uPar]);
        }
        //GetNumSourcePar

        //SourceSystems[uSource]->GetCk()->SetSourcePar(0, Pars[ParentSource[uSource]*NumPar+p_sor0]);
//this is quite badly implemented! The idea is to check if there are more than 1 source parameters, however at the moment
//I only check the main system, but there is no check if all systems are defined equivalently!
        //if(!FixPar[0][p_sor1])
            //SourceSystems[uSource]->GetCk()->SetSourcePar(1, Pars[ParentSource[uSource]*NumPar+p_sor1]);
        //if(!FixPar[0][p_sor2])
            //SourceSystems[uSource]->GetCk()->SetSourcePar(2, Pars[ParentSource[uSource]*NumPar+p_sor2]);
        //if(!FixPar[0][p_sor3])
            //SourceSystems[uSource]->GetCk()->SetSourcePar(3, Pars[ParentSource[uSource]*NumPar+p_sor3]);
    }
/*
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;
        for(unsigned uPot=0; uPot<NumPar; uPot++){
            if(uPot>=SystemToFit[uSyst]->GetCk()->GetNumPotPar()) break;
            SystemToFit[uSyst]->GetCk()->SetPotPar(uPot, Pars[uSyst*NumPar+p_pot0+uPot]);
//printf("uSyst=%u; uPot=%u; %e\n",uSyst,uPot,Pars[uSyst*NumPar+p_pot0+uPot]);
        }
    }
*/
    for(unsigned uPotential=0; uPotential<NumPotentialSystems; uPotential++){
        for(unsigned uPar=0; uPar<PotentialSystems[uPotential]->GetCk()->GetNumPotPar(); uPar++){
            PotentialSystems[uPotential]->GetCk()->SetPotPar(uPar, Pars[ParentPotential[uPotential]*NumPar+p_pot0+uPar]);
//printf(" uPotential=%u; uPar=%u;  --> %e\n",uPotential,uPar,Pars[ParentPotential[uPotential]*NumPar+p_pot0+uPar]);
        }
    }

//usleep(2e6);

    //update all systems. It is better to do this here than in the loop before,
    //in order to make sure that ALL systems have their sources adjusted before the update
    for(unsigned uSource=0; uSource<NumSourceSystems; uSource++){
        SourceSystems[uSource]->Update(false);
    }
    for(unsigned uPotential=0; uPotential<NumPotentialSystems; uPotential++){
        PotentialSystems[uPotential]->Update(false);
    }
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        SystemToFit[uSyst]->Update(false);
    }

    double CkVal;
    double BlVal = Pars[WhichSyst*NumPar+p_a]+Pars[WhichSyst*NumPar+p_b]*Momentum+Pars[WhichSyst*NumPar+p_c]*Momentum*Momentum;

    double Clin;
    //Oli's way
    if(SeparateBaseLineFit[WhichSyst]){
        if(FemtoRegion){
            CkVal = SystemToFit[WhichSyst]->EvalCk(Momentum)/fabs(Pars[WhichSyst*NumPar+p_Cl]);
        }
        else{
            CkVal = 1;
        }
    }
    else{
        if(FemtoRegion){
            CkVal = SystemToFit[WhichSyst]->EvalCk(Momentum);
        }
        else{
            if(Pars[WhichSyst*NumPar+p_Cl]==-1e6){
                Clin = SystemToFit[WhichSyst]->EvalCk(FitRange[WhichSyst][kf]);
            }
            else{
                Clin = fabs(Pars[WhichSyst*NumPar+p_Cl]);
            }
            if(Pars[WhichSyst*NumPar+p_kc]==-1e6){
                CkVal = Clin;
            }
            else if(Pars[WhichSyst*NumPar+p_kc]==FitRange[WhichSyst][kl]){
                CkVal = 1e6;
            }
            else{
                CkVal = ((Momentum-FitRange[WhichSyst][kl])-Clin*(Momentum-Pars[WhichSyst*NumPar+p_kc]))/
                        (Pars[WhichSyst*NumPar+p_kc]-FitRange[WhichSyst][kl]);
            }
            if(Clin<1 && CkVal>1) CkVal=1;
            if(Clin>1 && CkVal<1) CkVal=1;
        }
    }

    //this is here to make sure that we do not get an unphysical correlation function. If we do we put a HUGE value to CkVal to make the chi2 very bad
    //so that the fitter eventually rejects such combinations.
    if(RemoveNegCk){
        //does not really work yet
        double CkMain = SystemToFit[WhichSyst]->EvalMain(Momentum);
        //if(CkMain<0) CkVal-=fabs(CkMain*1000);
        if(CkMain<0) CkVal-=(exp(-CkMain*32)-1);
        //CkVal=-1e16;

        //printf("SystemToFit[WhichSyst]->EvalMain(%.3f)=%.3f\n", Momentum, SystemToFit[WhichSyst]->EvalMain(Momentum));
    }

    return BlVal*CkVal;
}
