#include "DLM_Fitters.h"


DLM_Fitter1::DLM_Fitter1(const unsigned& maxnumsyst):MaxNumSyst(maxnumsyst),NumPar(9),NumRangePar(4){
    HistoOriginal = new TH1F* [MaxNumSyst];
    HistoToFit = new TH1F* [MaxNumSyst];
    SystemToFit = new DLM_CkDecomposition* [MaxNumSyst];
    SourceSystems = NULL;
    ParentSource = NULL;
    NumSourceSystems = 0;
    FitRange = new double* [MaxNumSyst];
    NumSourceMapEntries = 0;
    SameSourceMap = NULL;
    NumSameSourcePar = NULL;
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
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        ParValue[uSyst] = new double [NumPar];
        ParDownLimit[uSyst] = new double [NumPar];
        ParUpLimit[uSyst] = new double [NumPar];
        FixPar[uSyst] = new bool [NumPar];
        //fix the pot pars
        for(unsigned uPar=0; uPar<NumPar; uPar++){
            if(uPar<5) FixPar[uSyst][uPar] = false;
            else FixPar[uSyst][uPar] = true;
        }
        SeparateBaseLineFit[uSyst] = true;
        ParValue[uSyst][p_a] = 1; ParDownLimit[uSyst][p_a] = 0.5; ParUpLimit[uSyst][p_a] = 2;
        ParValue[uSyst][p_b] = 1e-3; ParDownLimit[uSyst][p_b] = 1e-5; ParUpLimit[uSyst][p_b] = 1;
        ParValue[uSyst][p_r] = 1.5; ParDownLimit[uSyst][p_r] = 1; ParUpLimit[uSyst][p_r] = 2.5;
        ParValue[uSyst][p_Cl] = 1; ParDownLimit[uSyst][p_Cl] = 0.75; ParUpLimit[uSyst][p_Cl] = 1.25;
        ParValue[uSyst][p_kc] = 300; ParDownLimit[uSyst][p_kc] = 150; ParUpLimit[uSyst][p_kc] = 900;
        ParValue[uSyst][p_pot0]=0; ParValue[uSyst][p_pot1]=0; ParValue[uSyst][p_pot2]=0; ParValue[uSyst][p_pot3]=0;
    }

    HistoGlobal = NULL;
    FitGlobal = NULL;
    GlobalToMomentum = NULL;
    NumBinsSyst = NULL;
    CumulativeNumBinsSyst = NULL;

    OutputDirName = "./";

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

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        delete [] ParValue[uSyst];
        delete [] ParDownLimit[uSyst];
        delete [] ParUpLimit[uSyst];
        delete [] FixPar[uSyst];
    }
    delete [] ParValue; ParValue=NULL;
    delete [] ParDownLimit; ParDownLimit=NULL;
    delete [] ParUpLimit; ParUpLimit=NULL;
    delete [] FixPar; FixPar=NULL;
    delete [] SeparateBaseLineFit; SeparateBaseLineFit=NULL;

    if(HistoGlobal){delete HistoGlobal; HistoGlobal=NULL;}
    if(FitGlobal){delete FitGlobal; FitGlobal=NULL;}
    if(GlobalToMomentum){delete[]GlobalToMomentum; GlobalToMomentum=NULL;}
    if(NumBinsSyst){delete[]NumBinsSyst;NumBinsSyst=NULL;}
    if(CumulativeNumBinsSyst){delete[]CumulativeNumBinsSyst;CumulativeNumBinsSyst=NULL;}
//printf("~SourceSystems=%p\n",SourceSystems);
    if(SourceSystems) {delete[]SourceSystems; SourceSystems=NULL;}
    if(ParentSource) {delete[]ParentSource; ParentSource=NULL;}
}

//void DLM_Fitter1::TEST1(const unsigned& WhichSyst, TH1F* histo, const double& FromMeV ,
//                   DLM_CkDecomposition* decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX){
//HistoToFit[0] = new TH1F(TString::Format("HistoToFit%u",0),TString::Format("HistoToFit%u",0),10,0,1);
//return;
//}

void DLM_Fitter1::SetSystem(const unsigned& WhichSyst, TH1F& histo, const double& FromMeV ,
                   DLM_CkDecomposition& decomp, const double& KMIN, const double& KFEMTO, const double& KLINEAR, const double& KMAX){

    if(WhichSyst>=MaxNumSyst){
        printf("WARNING: SetSystem says WhichSyst>=MaxNumSyst\n");
        return;
    }
    if(KMIN>KFEMTO || KMIN>KLINEAR || KMIN>KMAX ||
       KFEMTO>KLINEAR || KFEMTO>KMAX || KLINEAR>KMAX){
        printf("WARNING: SetSystem says that the input ranges make no sense!\n");
        return;
    }
    if(KMIN<histo.GetBinLowEdge(1)/FromMeV || KMAX>histo.GetXaxis()->GetBinUpEdge(histo.GetNbinsX())/FromMeV){
//printf("KMIN=%f --> %f",KMIN,GetBinLowEdge(1)/FromMeV);
//printf(" KMAX=%f --> %f",KMIN,GetBinLowEdge(1)/FromMeV);
        printf("WARNING: SetSystem says that the desired fitting region is outside of the histo range!\n");
        return;
    }

    unsigned BinMin = histo.FindBin(KMIN/FromMeV);
    //unsigned BinFemto = histo->FindBin(KFEMTO/FromMeV);
    //unsigned BinLinear = histo->FindBin(KLINEAR/FromMeV);
    unsigned BinMax = histo.FindBin(KMAX/FromMeV);

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

    if(HistoToFit[WhichSyst]){
        delete [] HistoToFit[WhichSyst];
    }

//printf("HistoToFit = %p\n",HistoToFit);
//printf(" HistoToFit[%u] = %p\n",WhichSyst,HistoToFit[WhichSyst]);
//delete [] HistoToFit;
//HistoToFit = new TH1F* [10];
//TH1F* h1 = new TH1F(TString::Format("HistoToFit%u",WhichSyst),TString::Format("HistoToFit%u",WhichSyst),NumBins,Bins);
    HistoToFit[WhichSyst] = new TH1F(TString::Format("HistoToFit%u",WhichSyst),TString::Format("HistoToFit%u",WhichSyst),NumBins,Bins);
    //HistoToFit[WhichSyst] = new TH1F(TString::Format("HistoToFit%u",WhichSyst),TString::Format("HistoToFit%u",WhichSyst),10,0,1);

//HistoToFit[WhichSyst] = h1;

    for(unsigned uBin=BinMin; uBin<=BinMax; uBin++){
        HistoToFit[WhichSyst]->SetBinContent(uBin-BinMin+1, histo.GetBinContent(uBin));
        HistoToFit[WhichSyst]->SetBinError(uBin-BinMin+1, histo.GetBinError(uBin));
    }

    delete [] Bins;
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
    printf("WARNING: SetParameter says that the system '%s' does not exist!\n", WhichSyst.Data());
    delete [] buffer;
}
double DLM_Fitter1::GetParameter(const unsigned& WhichSyst, const unsigned& WhichPar){
    if(WhichSyst>=MaxNumSyst) return 0;
    if(WhichPar>=NumPar) return 0;
    return FitGlobal->GetParameter(WhichSyst*NumPar+WhichPar);
}
double DLM_Fitter1::GetParError(const unsigned& WhichSyst, const unsigned& WhichPar){
    if(WhichSyst>=MaxNumSyst) return 0;
    if(WhichPar>=NumPar) return 0;
    return FitGlobal->GetParError(WhichSyst*NumPar+WhichPar);
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
double DLM_Fitter1::Eval(const unsigned& WhichSyst, const double& Momentum)(){
    if(WhichSyst>=MaxNumSyst) return 0;




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
        TF1* FitBL = new TF1(TString::Format("FitBL%u",uSyst),"[0]+[1]*x",FitRange[uSyst][kl],FitRange[uSyst][kmax]);
        if(FixPar[uSyst][p_a]) {FitBL->FixParameter(0, ParValue[uSyst][p_a]);}
        else{FitBL->SetParameter(0, ParValue[uSyst][p_a]); FitBL->SetParLimits(0, ParDownLimit[uSyst][p_a], ParUpLimit[uSyst][p_a]);}
        if(FixPar[uSyst][p_b]) {FitBL->FixParameter(1, ParValue[uSyst][p_b]);}
        else{FitBL->SetParameter(1, ParValue[uSyst][p_b]); FitBL->SetParLimits(1, ParDownLimit[uSyst][p_b], ParUpLimit[uSyst][p_b]);}
        HistoToFit[uSyst]->Fit(FitBL,"S, N, R, M");
        FixParameter(uSyst,p_a,FitBL->GetParameter(0));
        FixParameter(uSyst,p_b,FitBL->GetParameter(1));
        FixParameter(uSyst,p_kc,0);
        delete FitBL;
    }


    double Momentum;
    double CkVal;
    double CkValErr;
    unsigned NumOfBins;

    //this might be needed if I allow for zero pointers to some of the main functions
    //for the time being this is not the case
    unsigned ActualNumSyst = 0;

    NumSourceSystems=0;
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        NumSourceSystems += bool(SystemToFit[uSyst]);
    }
    char* buffer = new char[128];
    char* buffer2 = new char[128];
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

printf("NumSourceSystems=%u\n",NumSourceSystems);
    if(SourceSystems){delete[]SourceSystems;SourceSystems=NULL;}
    if(NumSourceSystems) SourceSystems = new DLM_CkDecomposition*[NumSourceSystems];
    if(ParentSource){delete[]ParentSource;ParentSource=NULL;}
    if(NumSourceSystems) ParentSource = new int [NumSourceSystems];
printf("SourceSystems=%p\n",SourceSystems);
    unsigned uActualSource=0;
    //set the addresses of the main Ck functions to the SourceSystems
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;
        SourceSystems[uActualSource] = SystemToFit[uSyst];
        //by default the parent of the system is the system itself
        ParentSource[uActualSource]=uSyst;
        SystemToFit[uSyst]->GetName(buffer);

        //find the parent
        /*
        for(unsigned uSyst2=0; uSyst2<MaxNumSyst; uSyst2++){
            if(!SystemToFit[uSyst2]) continue;
            SystemToFit[uSyst2]->GetName(buffer2);
            //if(strcmp(buffer,SameSourceMap[uActualSource][0].Data())==0) ParentSource[uActualSource]=uSyst2;
            //if(strcmp(buffer,SameSourceMap[uActualSource][1].Data())==0) ParentSource[uActualSource]=uSyst2;
            if(strcmp(buffer,buffer2)==0) ParentSource[uActualSource]=uSyst2;
        }
        */
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
    delete [] buffer;
    delete [] buffer2;
printf("QA: %u=%u ?\n",NumSourceSystems,uActualSource);

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
    HistoGlobal = new TH1F("HistoGlobal", "HistoGlobal", NumBinsGlobal, 0, 1);
    if(FitGlobal) delete FitGlobal;
    FitGlobal = new TF1("FitGlobal", this, &DLM_Fitter1::EvalGlobal,0,1,ActualNumSyst*NumPar,"DLM_Fitter1","EvalGlobal");

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
            if(int(uPar)==p_r && ParentSource[uSyst]!=int(uSyst)){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
            }
            //if kf==kl, than there is an active boundary condition for the linear part of C(k)
            //thus Cl will be determined by the fitter, here we set a dummy value
            if(int(uPar)==p_Cl && FitRange[uSyst][kl]==FitRange[uSyst][kf]){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
            }

            //if we want to fit with a flat long-range C(k)
            if(int(uPar)==p_kc && ParValue[uSyst][p_Cl]<=0){
                FitGlobal->FixParameter(uSyst*NumPar+uPar, -1e6);
//printf("HEY from %u\n",uSyst*NumPar+uPar);
            }
        }

        //if there is no long-range fit, fix the variable that are modeling this part
        if(FitRange[uSyst][kl]==FitRange[uSyst][kmax]){
            FitGlobal->FixParameter(uSyst*NumPar+p_Cl, -1e6);
            FitGlobal->FixParameter(uSyst*NumPar+p_kc, -1e6);
        }

        //uActSyst++;
    }




printf("BEFORE FIT!\n");
    HistoGlobal->Fit(FitGlobal,"S, N, R, M");

//!TEMP///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    TFile* TEMPFILE = new TFile(TString::Format("%sTEMPFILE.root",OutputDirName.Data()), "recreate");
    HistoGlobal->Write();
    FitGlobal->Write();

    TCanvas* TEMPCAN = new TCanvas("TEMPCAN","TEMPCAN",1);
    TEMPCAN->cd(0);
    HistoGlobal->Draw();
    FitGlobal->SetNpx(3*NumBinsGlobal);
    FitGlobal->Draw("same");
    TEMPCAN->Write();
    printf("Chi2/NDF = %.2f/%i\n",FitGlobal->GetChisquare(),FitGlobal->GetNDF());
    printf(" pval=%.4f\n",FitGlobal->GetProb());
    for(unsigned uBin=0; uBin<NumBinsGlobal; uBin++){
        //printf("GH = %f <-> FIT = %f\n",HistoGlobal->GetBinContent(uBin+1),FitGlobal->Eval(HistoGlobal->GetBinCenter(uBin+1)));
    }
*/
//!TEMPEND////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//!TEMP///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //delete TEMPCAN;
    //delete TEMPFILE;
//!TEMPEND////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

double DLM_Fitter1::EvalGlobal(double* xVal, double* Pars){
//printf("CALLING EVALGLOBAL\n");
    unsigned GlobalBin = HistoGlobal->FindBin(*xVal)-1;
    double Momentum = GlobalToMomentum[GlobalBin];
    unsigned WhichSyst;

    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!HistoToFit[uSyst]) continue;
        //determine in which system the GlobalBin is
        if(GlobalBin<CumulativeNumBinsSyst[uSyst]) {WhichSyst=uSyst; break;}
    }

    if(!SystemToFit[WhichSyst]) return 0;

    //determine if we are in the femto region. If false => we are in the linear or flat region
    bool FemtoRegion = (Momentum<=FitRange[WhichSyst][kf]);
//printf("xVal=%f (k=%.1f); WhichSyst=%u; --> FEMTO=%i\n",*xVal,Momentum,WhichSyst,int(FemtoRegion));
    //loop over all systems that may need to have their source size adjusted
    //adjust the source of their corresponding "source-parent"
    for(unsigned uSource=0; uSource<NumSourceSystems; uSource++){
//printf(" set source par [%u]\n",uSource);
//printf(" SourceSystems[uSource]=%p\n",SourceSystems[uSource]);
//char buffer[32];
//SourceSystems[uSource]->GetName(buffer);
//printf(" Name=%s\n",buffer);
//printf(" Adjust the source of %p\n",SourceSystems[uSource]->GetCk());
        SourceSystems[uSource]->GetCk()->SetSourcePar(0, Pars[ParentSource[uSource]*NumPar+p_r]);
//printf(" set source par [%u] done\n",uSource);
    }


    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;
        for(unsigned uPot=0; uPot<NumPar; uPot++){
            if(uPot>=SystemToFit[uSyst]->GetCk()->GetNumPotPar()) break;
//printf("uSyst*NumPar+p_pot0+uPot=%u\n",uSyst*NumPar+p_pot0+uPot);
//printf("   VAL=%f\n",Pars[uSyst*NumPar+p_pot0+uPot]);
            SystemToFit[uSyst]->GetCk()->SetPotPar(uPot, Pars[uSyst*NumPar+p_pot0+uPot]);
        }
    }

    //update all systems. It is better to do this here than in the loop before,
    //in order to make sure that ALL systems have their sources adjusted before the update
    for(unsigned uSource=0; uSource<NumSourceSystems; uSource++){
//printf("uSource=%u; r=%f\n",uSource,Pars[ParentSource[uSource]*NumPar+p_r]);
        SourceSystems[uSource]->Update();
    }
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        SystemToFit[uSyst]->Update();
    }

/*
    for(unsigned uSyst=0; uSyst<MaxNumSyst; uSyst++){
        if(!SystemToFit[uSyst]) continue;

        //at the moment I can have only the radius size
        //if( SystemToFit[uSyst]->GetCk()->GetNumSourcePar() ){


            SystemToFit[uSyst]->GetCk()->SetSourcePar(0, Pars[uSyst*NumPar+p_r]);
        //}
    }

    if( SystemToFit[WhichSyst]->GetCk()->GetNumSourcePar() ){


        double& Radius = Pars[WhichSyst*NumPar+p_r];

    }
*/

//SystemToFit[WhichSyst]->GetCk()->SetSourcePar(0, )
//SetSourcePar

    double CkVal;
    double BlVal = Pars[WhichSyst*NumPar+p_a]+Pars[WhichSyst*NumPar+p_b]*Momentum;

    double Clin;
    //double Cconv;
//printf("SeparateBaseLineFit=%i\n",int(SeparateBaseLineFit));
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
//printf("I am at %u\n",WhichSyst*NumPar+p_kc);
            if(Pars[WhichSyst*NumPar+p_kc]==-1e6){
                CkVal = Clin;
//printf("CkVal=%f\n");
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
//printf(" CkVal=%f\n");
        }
    }

    return BlVal*CkVal;
}
