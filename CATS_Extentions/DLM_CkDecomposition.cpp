
#include "DLM_CkDecomposition.h"


#include "DLM_CkModels.h"


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
    if(Momentum<GetBinUpEdge(0,NumBins[0])) CutOff = Momentum;
    else CutOff = GetBinUpEdge(0,NumBins[0]);
    CutOff_kc = kc;
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
            BinValue[uBin] = Kitty->EvalCorrFun(GetBinCenter(0,uBin));
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
    else if(Momentum>BinRange[0][NumBins[0]]){
        kf = BinRange[0][NumBins[0]];
    }
    //Momentum>CutOff
    else{
        kf = CutOff;
    }
    const double Cf = DLM_Histo::Eval(&kf);
//printf("CutOff_kc=%.4f\n",CutOff_kc);
//printf("CutOff=%.4f\n",CutOff);
    if(CutOff_kc<0) return fabs(CutOff_kc);
    if(CutOff_kc<=kf) return Cf;
    double ReturnVal = ((Momentum-kf)-Cf*(Momentum-CutOff_kc))/(CutOff_kc-kf);
    if(ReturnVal>1) ReturnVal=1;
    return ReturnVal;
}

//the main contribution does not count as a child, i.e. 1 child implies one contribution additionally to the main one!
DLM_CkDecomposition::DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const TH2F* hSigmaMatrix, const bool& InvertAxis):
    ERROR_STATE(!name),NumChildren(numchildren),CkMain(&ckfunction){
DEBUGFLAG=0;
    Child = NULL;
    LambdaPar = NULL;
    Type = NULL;
    RM_Child = NULL;
    RM_MomResolution = NULL;
    CurrentStatus = false;
    DecompositionStatus = false;
    CurrentStatusChild = NULL;
    CkMainSmeared = NULL;
    CkMainFeed = NULL;
    CkSmearedMainFeed = NULL;

    if(ERROR_STATE){
        printf("\033[1;31mERROR:\033[0m The DLM_CkDecomposition got some rubbish input, the object will be broken!\n");
        return;
    }

    RM_MomResolution = new DLM_ResponseMatrix(ckfunction, hSigmaMatrix, NULL, InvertAxis);
    CkMainFeed = new DLM_Histo<double>(CkMain[0]);
    CkSmearedMainFeed = new DLM_Histo<double>(CkMain[0]);

    Name = new char [strlen(name)+1];
    strcpy(Name,name);

    LambdaMain = 1;
    MuPar = 1;

    if(NumChildren){
        Child = new DLM_CkDecomposition* [NumChildren];
        LambdaPar = new double [NumChildren];
        Type = new int [NumChildren];
        RM_Child = new DLM_ResponseMatrix* [NumChildren];
        CurrentStatusChild = new bool [NumChildren];
        CkChildMainFeed = new DLM_Histo<double>* [NumChildren];
        CkSmearedChildMainFeed = new DLM_Histo<double>* [NumChildren];
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            Child[uChild] = NULL;
            LambdaPar[uChild] = 0;
            Type[uChild] = -1;
            RM_Child[uChild] = NULL;
            CurrentStatusChild[uChild] = false;
            CkChildMainFeed[uChild] = NULL;
            CkSmearedChildMainFeed[uChild] = NULL;
        }
    }

    Update(true);

}

DLM_CkDecomposition::~DLM_CkDecomposition(){
    if(ERROR_STATE) return;
    if(NumChildren){
        delete [] Child; Child=NULL;
        delete [] LambdaPar; LambdaPar=NULL;
        delete [] Type; Type=NULL;
        delete [] CurrentStatusChild; CurrentStatusChild=NULL;

        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(RM_Child[uChild]) delete RM_Child[uChild];
            if(CkChildMainFeed[uChild]) delete CkChildMainFeed[uChild];
            if(CkSmearedChildMainFeed[uChild]) delete CkSmearedChildMainFeed[uChild];
        }
        delete [] RM_Child; RM_Child=NULL;
        delete [] CkChildMainFeed; CkChildMainFeed=NULL;
        delete [] CkSmearedChildMainFeed; CkSmearedChildMainFeed=NULL;
    }
    if(Name) {delete [] Name; Name=NULL;}
    if(RM_MomResolution) {delete RM_MomResolution; RM_MomResolution=NULL;}
    if(CkMainSmeared) {delete CkMainSmeared; CkMainSmeared=NULL;}
    if(CkMainFeed) {delete CkMainFeed; CkMainFeed=NULL;}
    if(CkSmearedMainFeed) {delete CkSmearedMainFeed; CkSmearedMainFeed=NULL;}
}

void DLM_CkDecomposition::AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomposition* child,
                                          const TH2F* hResidualMatrix, const bool& InvertedAxis){

    if(ERROR_STATE) return;

    if(WhichCk>=NumChildren){
        return;
    }
    if(type!=cFeedDown && type!=cFake){
        return;
    }

    if(fraction<0 || fraction>1){
        printf("\033[1;33mWARNING:\033[0m Trying the set a contribution with a fraction<0 || fraction>1.\n");
        return;
    }
    if(child && !UniqueName(child->Name)){
        printf("\033[1;33mWARNING:\033[0m Trying to duplicate the name '%s' of a contribution in DLM_CkDecomposition::SetContribution.\n", child->Name);
        return;
    }

    if( Child[WhichCk]==child && LambdaPar[WhichCk]==fraction && Type[WhichCk]==type){
        return;
    }

    if(type==cFake && hResidualMatrix){
        printf("\033[1;33mWARNING:\033[0m Impurities are not allowed to have a transformation matrix.\n");
        hResidualMatrix = NULL;
    }

    Child[WhichCk] = child;
    LambdaPar[WhichCk] = fraction;
    Type[WhichCk] = type;
    DecompositionStatus = false;

    if(RM_Child[WhichCk]) {delete RM_Child[WhichCk]; RM_Child[WhichCk]=NULL;}
    if(child && hResidualMatrix){
        RM_Child[WhichCk] = new DLM_ResponseMatrix(*child->CkMain, NULL, hResidualMatrix, InvertedAxis);
    }

    if(CkChildMainFeed[WhichCk]) {delete CkChildMainFeed[WhichCk]; CkChildMainFeed[WhichCk]=NULL;}
    if(child){
        CkChildMainFeed[WhichCk] = new DLM_Histo<double> (*child->CkMain);
    }

    LambdaMain=1;
    MuPar=1;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        //if(Child[uChild])
        LambdaMain-=LambdaPar[uChild];
        if(Type[uChild]==cFake) MuPar-=LambdaPar[uChild];
    }
}

double DLM_CkDecomposition::EvalCk(const double& Momentum){

    if(ERROR_STATE) return 0;
    if(LambdaMain<0){
        printf("\033[1;31mERROR:\033[0m The λ0 parameter is smaller than zero!\n");
        return 0;
    }
    Update(false);

    double VAL_CkSmearedMainFeed=0;
    double VAL_CkSmearedMainFake=0;
    //this is needed in order to get a smooth transition for the range outside of the histogram
    if(Momentum>CkMain->GetUpEdge(0)){
        VAL_CkSmearedMainFeed = CkMain->Eval(Momentum)*LambdaMain/MuPar;
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]==cFeedDown){
                //printf("Momentum=%.4f\n",Momentum);
                if(Child[uChild]){
                    VAL_CkSmearedMainFeed+=GetChild(uChild)->GetCk()->Eval(Momentum)*LambdaPar[uChild]/MuPar;
                    //printf("(1): Ck=%.4f; lam=%.4f; mu=%.4f\n",GetChild(uChild)->GetCk()->Eval(Momentum),LambdaPar[uChild],MuPar);
                }
                else{
                    VAL_CkSmearedMainFeed+=(LambdaPar[uChild]/MuPar);
                    //printf("(2): lam=%.4f; mu=%.4f\n",LambdaPar[uChild],MuPar);
                }
            }
            else{
                if(Child[uChild]){
                    VAL_CkSmearedMainFake+=GetChild(uChild)->GetCk()->Eval(Momentum)*LambdaPar[uChild];
                }
                else{
                    VAL_CkSmearedMainFake+=(LambdaPar[uChild]);
                }
            }
        }
    }
    else{
        VAL_CkSmearedMainFeed = CkSmearedMainFeed->Eval(&Momentum);
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]!=cFake) continue;
            if(Child[uChild]){
                VAL_CkSmearedMainFake += LambdaPar[uChild]*Child[uChild]->CkSmearedMainFeed->Eval(&Momentum);
            }
            else{
                VAL_CkSmearedMainFake += LambdaPar[uChild];
            }
        }
    }
//if(Momentum>400&&Momentum<405)
//printf("MuPar=%.4f; VAL_CkSmearedMainFeed=%.4f; VAL_CkSmearedMainFake=%.4f\n",MuPar,VAL_CkSmearedMainFeed,VAL_CkSmearedMainFake);
    return MuPar*VAL_CkSmearedMainFeed+VAL_CkSmearedMainFake;
}

double DLM_CkDecomposition::EvalMain(const double& Momentum){
    return CkMain->Eval(Momentum);
}
double DLM_CkDecomposition::EvalSmearedMain(const double& Momentum){
    if(!CkMainSmeared){
        CkMainSmeared = new DLM_Histo<double>(*CkMain);
        Smear(CkMain,RM_MomResolution,CkMainSmeared);
    }
    return CkMainSmeared->Eval(&Momentum);
}
double DLM_CkDecomposition::EvalMainFeed(const double& Momentum){
    return CkMainFeed->Eval(&Momentum);
}
double DLM_CkDecomposition::EvalSmearedMainFeed(const double& Momentum){
    return CkSmearedMainFeed->Eval(&Momentum);
}
double DLM_CkDecomposition::EvalSmearedFeed(const unsigned& WhichChild,const double& Momentum){
    if(WhichChild>=NumChildren) return 0;
    if(!CkSmearedChildMainFeed[WhichChild]&&Child[WhichChild]) {
        CkSmearedChildMainFeed[WhichChild] = new DLM_Histo<double> (*Child[WhichChild]->CkMain);
        Smear(CkChildMainFeed[WhichChild], RM_MomResolution, CkSmearedChildMainFeed[WhichChild]);
    }
    if(CkSmearedChildMainFeed[WhichChild]) return CkSmearedChildMainFeed[WhichChild]->Eval(&Momentum);
    else return 0;
}

unsigned DLM_CkDecomposition::GetNumChildren(){
    return NumChildren;
}


DLM_CkDecomposition* DLM_CkDecomposition::GetChild(const unsigned& WhichChild){
    if(WhichChild>=NumChildren) return NULL;
    else return Child[WhichChild];
}
DLM_CkDecomposition* DLM_CkDecomposition::GetContribution(const char* name){
    DLM_CkDecomposition* RESULT = NULL;
    if(strcmp(name,Name)==0) return this;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Child[uChild]) RESULT = Child[uChild]->GetContribution(name);
        if(RESULT) break;
    }
    return RESULT;
}
DLM_Histo<double>* DLM_CkDecomposition::GetChildContribution(const unsigned& WhichChild, const bool& WithLambda){
    if(WhichChild>=NumChildren) return NULL;
    DLM_Histo<double>* Histo = new DLM_Histo<double> (*CkChildMainFeed[WhichChild]);
    if(Type[WhichChild]!=cFake){
        Smear(CkChildMainFeed[WhichChild], RM_MomResolution, Histo);
    }
    if(WithLambda) Histo->Scale(LambdaPar[WhichChild]);
    return Histo;
}
DLM_Histo<double>* DLM_CkDecomposition::GetChildContribution(const char* name, const bool& WithLambda){
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(!Child[uChild]) continue;
        if(strcmp(name,Child[uChild]->Name)==0){
            return GetChildContribution(uChild,WithLambda);
        }
    }
    return NULL;
}
DLM_Ck* DLM_CkDecomposition::GetCk(){
    return CkMain;
}

void DLM_CkDecomposition::GetName(char* name){
    if(!name) return;
    strcpy(name,Name);
}

bool DLM_CkDecomposition::Status(){
    bool STATUS = CkMain->Status();
    STATUS *= DecompositionStatus;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(!STATUS) break;
        if(Child[uChild]){
            STATUS *= Child[uChild]->Status();
        }
    }
    return STATUS;
}


bool SHOWSTUFF=false;

//the strategy is as follows:
//1) check the status and update the main C(k)
//2) look if the children are up to date. If not, update them iteratively.
//3) smear all correlations that needed an update
void DLM_CkDecomposition::Update(const bool& FORCE_FULL_UPDATE){
    double Momentum;
    CurrentStatus = CkMain->Status();
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Child[uChild]){
            //always false if FORCE_FULL_UPDATE==true
            CurrentStatusChild[uChild] = Status()*(!FORCE_FULL_UPDATE);
        }
        //if the child is not defined -> take a flat C(k)
        else if(!CurrentStatusChild[uChild]){
            CurrentStatusChild[uChild] = true;
        }
    }

    if(!CurrentStatus || FORCE_FULL_UPDATE){
        CkMain->Update(FORCE_FULL_UPDATE);
    }

    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if( !CurrentStatusChild[uChild] ){
            Child[uChild]->Update(FORCE_FULL_UPDATE);
        }
    }

    if(!CurrentStatus || FORCE_FULL_UPDATE){

        CkMainFeed->Copy(CkMain[0]);
        CkMainFeed->Scale(LambdaMain/MuPar);

        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]!=cFeedDown) continue;
            if(Child[uChild]){
                Smear(Child[uChild]->CkMainFeed, RM_Child[uChild], CkChildMainFeed[uChild]);

                for(unsigned uBin=0; uBin<CkMainFeed->GetNbins(); uBin++){
                    Momentum = CkMainFeed->GetBinCenter(0,uBin);
                    CkMainFeed->Add(uBin, CkChildMainFeed[uChild]->Eval(&Momentum)*LambdaPar[uChild]/MuPar);
                }
            }
            else{
                *CkMainFeed+=(LambdaPar[uChild]/MuPar);
            }
        }
        Smear(CkMainFeed, RM_MomResolution, CkSmearedMainFeed);
    }
    DecompositionStatus = true;
}

void DLM_CkDecomposition::Smear(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared){
    if(!SmearMatrix){
        CkSmeared[0] = CkToSmear[0];
    }
    else{
        for(unsigned uBinSmear=0; uBinSmear<CkToSmear->GetNbins(); uBinSmear++){
            CkSmeared->SetBinContent(uBinSmear, 0);
            for(unsigned uBinTrue=0; uBinTrue<CkToSmear->GetNbins(); uBinTrue++){
                //as the response matrix is normalized to the size of the bin, during the integration we multiply for it
                CkSmeared->Add(uBinSmear,   SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue)*
                                            CkToSmear->GetBinSize(uBinTrue)*CkSmeared->GetBinSize(uBinSmear));
            }
        }
    }
}

//check if any of the children has the same name
bool DLM_CkDecomposition::UniqueName(const char* name){
    if(ERROR_STATE) return 0;
    if(strcmp(Name,name)==0) return false;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Child[uChild]){
            if(Child[uChild]->UniqueName(name)==false) return false;
        }
    }
    return true;
}
