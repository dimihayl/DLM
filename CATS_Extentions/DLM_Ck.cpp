
#include <stdio.h>
#include <string.h>

#include "DLM_Ck.h"
#include "DLM_CkModels.h"

DLM_Ck::DLM_Ck(const unsigned numbins,const double& kmin, const double& kmax):DLM_Histo(),NumSourcePar(0),NumPotPar(0){
    Kitty = NULL;
    CkFunction = NULL;
    SetUp(1);
    SetUp(0,numbins,kmin,kmax);
    Initialize();
    DefaultConstructor();
    SourceUpToDate = true;
    PotUpToDate = true;
}

DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat):DLM_Histo(),
                NumSourcePar(nSourcePar),NumPotPar(nPotPar){
    SetUp(1);
    for(unsigned uBin=1; uBin<=cat.GetNumMomBins(); uBin++){
        if(cat.GetMomBinLowEdge(uBin)<=cat.GetMomBinLowEdge(uBin-1)){
            printf("\033[1;31mERROR:\033[0m DLM_Ck: The bin ranges should be in ascending order and a bin-width of 0 is not allowed!\n");
            return;
        }
    }
    if(!BinRange[0]) BinRange[0] = new double[cat.GetNumMomBins()+1];
    if(!BinCenter[0]) BinCenter[0] = new double[cat.GetNumMomBins()];
    NumBins[0] = cat.GetNumMomBins();
    Initialize();
    for(unsigned uBin=0; uBin<=NumBins[0]; uBin++){
        BinRange[0][uBin] = cat.GetMomBinLowEdge(uBin);
        if(uBin){
            BinCenter[0][uBin-1] = cat.GetMomentum(uBin-1);
            BinValue[uBin-1] = cat.GetCorrFun(uBin-1);
            BinError[uBin-1] = cat.GetCorrFunErr(uBin-1);
        }
    }
    Kitty = &cat;
    CkFunction = NULL;
    DefaultConstructor();
}

DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double* bins):DLM_Histo(),
                NumSourcePar(nSourcePar),NumPotPar(nPotPar){
    SetUp(1);
    for(unsigned uBin=1; uBin<=cat.GetNumMomBins(); uBin++){
        if(cat.GetMomBinLowEdge(uBin)<=cat.GetMomBinLowEdge(uBin-1)){
            printf("\033[1;31mERROR:\033[0m DLM_Ck: The bin ranges should be in ascending order and a bin-width of 0 is not allowed!\n");
            return;
        }
    }

    SetUp(1);
    SetUp(0,numbin,bins);
    Initialize();

    Kitty = &cat;
    CkFunction = NULL;
    DefaultConstructor();
}

DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double& minMom, const double& maxMom):DLM_Histo(),
                NumSourcePar(nSourcePar),NumPotPar(nPotPar){
    SetUp(1);
    for(unsigned uBin=1; uBin<=cat.GetNumMomBins(); uBin++){
        if(cat.GetMomBinLowEdge(uBin)<=cat.GetMomBinLowEdge(uBin-1)){
            printf("\033[1;31mERROR:\033[0m DLM_Ck: The bin ranges should be in ascending order and a bin-width of 0 is not allowed!\n");
            return;
        }
    }

    SetUp(1);
    SetUp(0,numbin,minMom,maxMom);
    Initialize();

    Kitty = &cat;
    CkFunction = NULL;
    DefaultConstructor();
}

DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
        const unsigned& numbin, const double* bins, double (*CorrFun)(const double&, const double*, const double*)):DLM_Histo(),
        NumSourcePar(nSourcePar),NumPotPar(nPotPar),CkFunction(CorrFun){
    Kitty = NULL;
    SetUp(1);
    SetUp(0,numbin,bins);
    Initialize();
    DefaultConstructor();
}
DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
        const unsigned& numbin, const double& minMom, const double& maxMom, double (*CorrFun)(const double&, const double*, const double*)):DLM_Histo(),
        NumSourcePar(nSourcePar),NumPotPar(nPotPar),CkFunction(CorrFun){
    Kitty = NULL;
    SetUp(1);
    SetUp(0,numbin,minMom,maxMom);
    Initialize();
    DefaultConstructor();
}
DLM_Ck::~DLM_Ck(){
    if(SourcePar) {delete [] SourcePar;}
    if(PotPar) {delete [] PotPar;}
}

void DLM_Ck::DefaultConstructor(){
    SourcePar = NULL;
    PotPar = NULL;
    CutOff = BinRange[0][NumBins[0]];
    if(CkFunction){
        if(NumSourcePar) SourcePar = new double [NumSourcePar];
        if(NumPotPar) PotPar = new double [NumPotPar];
    }
    SourceUpToDate = false;
    PotUpToDate = false;
    CutOff = 1e6;
    CutOff_kc = -1;
}

void DLM_Ck::SetSourcePar(const unsigned& WhichPar, const double& Value){
    if(WhichPar>=NumSourcePar) return;
    if(CkFunction){
        if(SourcePar[WhichPar]==Value) return;
        SourcePar[WhichPar]=Value;
        SourceUpToDate = false;
    }
    else if(Kitty && Kitty->GetUseAnalyticSource()){
        if(Kitty->GetAnaSourcePar(WhichPar)==Value) return;
        Kitty->SetAnaSource(WhichPar, Value, true);
        SourceUpToDate = false;
    }
    else if(Kitty && !Kitty->GetUseAnalyticSource()){
        if(Kitty->GetTransportRenorm()==Value) return;
        Kitty->SetPoorManRenorm(Value);
        SourceUpToDate = false;
    }
    else{
        return;
    }
}
double DLM_Ck::GetSourcePar(const unsigned& WhichPar){
    if(WhichPar>=NumSourcePar) return 0;
    if(CkFunction){
        return SourcePar[WhichPar];
    }
    else if(Kitty && Kitty->GetUseAnalyticSource()){
        return Kitty->GetAnaSourcePar(WhichPar);
    }
    else{
        return 0;
    }
}

unsigned DLM_Ck::GetNumSourcePar(){
    return NumSourcePar;
}

double DLM_Ck::GetPotPar(const unsigned& WhichPar){
    if(WhichPar>=NumPotPar) return 0;
    return PotPar[WhichPar];
}
unsigned DLM_Ck::GetNumPotPar(){
    return NumPotPar;
}
void DLM_Ck::SetCutOff(const double& Momentum, const double& kc){
    if(Momentum<GetBinUpEdge(0,NumBins[0]-1)) CutOff = Momentum;
    else CutOff = GetBinUpEdge(0,NumBins[0]-1);
    CutOff_kc = kc;
}
double DLM_Ck::GetCutOff() const{
    return CutOff;
}
double DLM_Ck::GetCutOff_kc() const{
    return CutOff_kc;
}

void DLM_Ck::SetPotPar(const unsigned& WhichPar, const double& Value){
    if(WhichPar>=NumPotPar) return;
    if(CkFunction){
        if(PotPar[WhichPar]==Value) return;
        PotPar[WhichPar]=Value;
        PotUpToDate = false;
    }
    //!N.B. we change the values for ALL potentials.
    //This is done so for ease when fitting. The idea is that the potential function should be the same for all
    //channels and PWs, but simply should have the s,l,j numbers as arguments.
    //However you must take extra care that you do not change s,l,j, since this will be done for all channels!
    else if(Kitty){
        for(unsigned short usCh=0; usCh<Kitty->GetNumChannels(); usCh++){
            for(unsigned short usPw=0; usPw<Kitty->GetNumPW(usCh); usPw++){
                Kitty->SetShortRangePotential(usCh,usPw,WhichPar,Value);
            }
        }
        PotUpToDate = Kitty->PotentialStatus();
    }
    else{
        return;
    }
}

bool DLM_Ck::Status(){
    return (SourceUpToDate && PotUpToDate);
}

//maybe a bug, if I change Transport to Ana -> no update unless I force it
void DLM_Ck::Update(const bool& FORCE){
//printf("SourceUpToDate = %i; PotUpToDate=%i; FORCE=%i\n",SourceUpToDate,PotUpToDate,FORCE);
    if(SourceUpToDate && PotUpToDate && !FORCE) return;
    if(CkFunction){
        for(unsigned uBin=0; uBin<NumBins[0]; uBin++){
            BinValue[uBin] = CkFunction(GetBinCenter(0,uBin), SourcePar, PotPar);
        }
        SourceUpToDate = true;
        PotUpToDate = true;
    }
    else if(Kitty){
        short NOTIF = Kitty->GetNotifications();
        Kitty->SetNotifications(NOTIF==(CATS::nSilent)?NOTIF:CATS::nError);
        Kitty->KillTheCat();
        Kitty->SetNotifications(NOTIF);
        for(unsigned uBin=0; uBin<NumBins[0]; uBin++){
            if(Kitty->GetMomBinUpEdge(Kitty->GetNumMomBins()-1)>GetBinCenter(0,uBin) &&
               Kitty->GetMomBinLowEdge(0)<GetBinCenter(0,uBin) &&
               CutOff>GetBinCenter(0,uBin)){
               BinValue[uBin] = Kitty->EvalCorrFun(GetBinCenter(0,uBin));
//printf("if %f -> %f\n",GetBinCenter(0,uBin),BinValue[uBin]);
            }
            else if(Kitty->GetMomBinLowEdge(0)>=GetBinCenter(0,uBin)){
//printf("Kitty->GetMomBinLowEdge(0)>=GetBinCenter(0,uBin) : %f>=%f",Kitty->GetMomBinLowEdge(0),GetBinCenter(0,uBin));
                BinValue[uBin] = 0;
            }
            //if we go to higher values than the CATS object, we assume a flat correlation equal to one,
            //unless there is the additional condition for slow linear convergence towards unity after the CutOff
            else{
                //Eval will compute automatically the corresponding value in case of a CutOff
                BinValue[uBin] = Eval(GetBinCenter(0,uBin));
//printf("else %f -> %f\n",GetBinCenter(0,uBin),BinValue[uBin]);
            }
            //BinValue[uBin] = 0;
        }
        SourceUpToDate = true;
        PotUpToDate = true;
    }
    else{
        return;
    }
}
double DLM_Ck::Eval(const double& Momentum){
    if(Momentum<BinRange[0][0]) return 0;
    double kf;
    if(Momentum<CutOff&&Momentum<BinRange[0][NumBins[0]]){
        return DLM_Histo::Eval(&Momentum);
    }
    else if(CutOff>BinRange[0][NumBins[0]]){
        kf = BinRange[0][NumBins[0]];
    }
    //Momentum>CutOff
    else{
        kf = CutOff;
    }
///THIS IS STILL BUGGY, QUICK FIXED !!!!!!!!!!!!
    //const double Cf = DLM_Histo::Eval(&kf);
    const double Cf = BinValue[GetBin(0,kf)-1];
    //printf("Cf = %f\n",Cf);
    if(CutOff_kc<0) return fabs(CutOff_kc);
    if(CutOff_kc<=kf) return Cf;
    double ReturnVal = ((Momentum-kf)-Cf*(Momentum-CutOff_kc))/(CutOff_kc-kf);
    if(ReturnVal>1) ReturnVal=1;
    return ReturnVal;
}
