
#include "DLM_CkDecomposition.h"


#include "DLM_CkModels.h"


/*
DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat):
    NumSourcePar(nSourcePar),NumPotPar(nPotPar),MomBinCopy(cat.CopyMomBin()),Kitty(&cat),
    CATShisto<double>::CATShisto(cat.GetNumMomBins(),MomBinCopy){
printf("Hell yeah!\n");
    CkFunction = NULL;
    delete [] MomBinCopy; MomBinCopy=NULL;
    DefaultConstructor();
}
*/
DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat):
                CATShisto<double>::CATShisto(cat.GetNumMomBins()),NumSourcePar(nSourcePar),NumPotPar(nPotPar){
    for(unsigned uBin=0; uBin<NumBins; uBin++){
        BinRange[uBin] = cat.GetMomBinLowEdge(uBin);
        BinValue[uBin] = cat.GetCorrFun(uBin);
        BinCenter[uBin] = cat.GetMomentum(uBin);
    }
    BinRange[NumBins] = cat.GetMomBinUpEdge(NumBins-1);

    Kitty = &cat;
    CkFunction = NULL;
    DefaultConstructor();
}

DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
        const unsigned& numbin, const double* bins, double (*CorrFun)(const double&, const double*, const double*)):
        CATShisto<double>::CATShisto(numbin,bins),NumSourcePar(nSourcePar),NumPotPar(nPotPar),CkFunction(CorrFun){
    Kitty = NULL;
    DefaultConstructor();
}
DLM_Ck::DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
        const unsigned& numbin, const double& minMom, const double& maxMom, double (*CorrFun)(const double&, const double*, const double*)):
        CATShisto<double>::CATShisto(numbin,minMom,maxMom),NumSourcePar(nSourcePar),NumPotPar(nPotPar),CkFunction(CorrFun){
    Kitty = NULL;
    DefaultConstructor();
}
DLM_Ck::~DLM_Ck(){
    if(SourcePar) {delete [] SourcePar;}
    if(PotPar) {delete [] PotPar;}
}

void DLM_Ck::DefaultConstructor(){
    SourcePar = NULL;
    PotPar = NULL;
    if(CkFunction){
        if(NumSourcePar) SourcePar = new double [NumSourcePar];
        if(NumPotPar) PotPar = new double [NumPotPar];
//printf("I WAS in %p with %u and %u --> PotPar=%p\n",this,NumSourcePar,NumPotPar,PotPar);
    }
    SourceUpToDate = false;
    PotUpToDate = false;
}

void DLM_Ck::SetSourcePar(const unsigned& WhichPar, const double& Value){
//printf(" WhichPar=%u (%u)\n", WhichPar, NumSourcePar);
    if(WhichPar>=NumSourcePar) return;
    if(CkFunction){
        if(SourcePar[WhichPar]==Value) return;
        SourcePar[WhichPar]=Value;
        SourceUpToDate = false;
//printf(" SOURCE CHANGED %p!\n",this);
    }
    else if(Kitty && Kitty->GetUseAnalyticSource()){
//printf(" Kitty=%p; r=%.2f; SUD=%i\n",Kitty,Kitty->GetAnaSourcePar(WhichPar),int(SourceUpToDate));
        if(Kitty->GetAnaSourcePar(WhichPar)==Value) return;
        Kitty->SetAnaSource(WhichPar, Value, true);
        SourceUpToDate = false;
//printf(" Kitty SOURCE CHANGED %p!\n",this);
    }
    else if(Kitty && !Kitty->GetUseAnalyticSource()){
        if(Kitty->GetTransportRenorm()==Value) return;
        //Kitty->SetTransportRenorm(Value);
        Kitty->SetPoorManRenorm(Value);
        SourceUpToDate = false;
    }
    else{
        return;
//printf(" NO CHANGE :( wtf!!!\n");
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
unsigned DLM_Ck::GetNumPotPar(){
    return NumPotPar;
}

void DLM_Ck::SetPotPar(const unsigned& WhichPar, const double& Value){
//printf("this=%p; WhichPar = %u (NumPotPar=%u) -> PotPar=%p\n",this,WhichPar,NumPotPar,PotPar);
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
//PotUpToDate=true;
//printf("%p: SourceUpToDate=%i; PotUpToDate=%i\n",this,int(SourceUpToDate),int(PotUpToDate));
    if(SourceUpToDate && PotUpToDate && !FORCE) return;
    if(CkFunction){
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            BinValue[uBin] = CkFunction(GetBinCenter(uBin), SourcePar, PotPar);
/*
if(BinValue[uBin]<0) {
printf("BinValue[%.2f]=%.2f\n",GetBinCenter(uBin),BinValue[uBin]);
printf("R=%.2f\n",SourcePar[0]);
printf("PotPar=%.2f and %.2f\n",PotPar[0],PotPar[1]);
printf("%.2f\n",Lednicky_Identical_Singlet(GetBinCenter(uBin), SourcePar, PotPar));
}
*/
        }
        SourceUpToDate = true;
        PotUpToDate = true;
    }
    else if(Kitty){
        short NOTIF = Kitty->GetNotifications();
        Kitty->SetNotifications(NOTIF==(CATS::nSilent)?NOTIF:CATS::nError);
        Kitty->KillTheCat();
        Kitty->SetNotifications(NOTIF);
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            BinValue[uBin] = Kitty->EvalCorrFun(GetBinCenter(uBin));
            //BinValue[uBin] = 0;
        }
        SourceUpToDate = true;
        PotUpToDate = true;
    }
    else{
        return;
    }
}

//the main contribution does not count as a child, i.e. 1 child implies one contribution additionally to the main one!
DLM_CkDecomposition::DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, TH2F* hSigmaMatrix, const bool& InvertAxis):
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
    CkMainFeed = new CATShisto<double>(CkMain[0]);
    CkSmearedMainFeed = new CATShisto<double>(CkMain[0]);

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
        CkChildMainFeed = new CATShisto<double>* [NumChildren];
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            Child[uChild] = NULL;
            LambdaPar[uChild] = 0;
            Type[uChild] = -1;
            RM_Child[uChild] = NULL;
            CurrentStatusChild[uChild] = false;
            CkChildMainFeed[uChild] = NULL;
        }
    }

    Update(true);

}

DLM_CkDecomposition::~DLM_CkDecomposition(){
return;
    if(ERROR_STATE) return;
    if(NumChildren){
        delete [] Child; Child=NULL;
        delete [] LambdaPar; LambdaPar=NULL;
        delete [] Type; Type=NULL;
        delete [] CurrentStatusChild; CurrentStatusChild=NULL;

        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(RM_Child[uChild]) delete RM_Child[uChild];
            if(CkChildMainFeed[uChild]) delete CkChildMainFeed[uChild];
        }
        delete [] RM_Child; RM_Child=NULL;
        delete [] CkChildMainFeed; CkChildMainFeed=NULL;
    }
    if(Name) {delete [] Name; Name=NULL;}
    if(RM_MomResolution) {delete RM_MomResolution; RM_MomResolution=NULL;}
    if(CkMainSmeared) {delete CkMainSmeared; CkMainSmeared=NULL;}
    if(CkMainFeed) {delete CkMainFeed; CkMainFeed=NULL;}
    if(CkSmearedMainFeed) {delete CkSmearedMainFeed; CkSmearedMainFeed=NULL;}
}

void DLM_CkDecomposition::AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomposition* child,
                                          TH2F* hResidualMatrix, const bool& InvertedAxis){
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
        //CkChildMainFeed[WhichCk] = new CATShisto<double> (*CkMain);
        CkChildMainFeed[WhichCk] = new CATShisto<double> (*child->CkMain);

//if(strcmp("pLambda",Name)==0){
//printf("I want to have kiddies! --> %u with name %p\n",WhichCk,CkChildMainFeed[WhichCk]);
//}
//why do I have this???
/*
        for(unsigned uBin=0; uBin<CkMain->GetNbins(); uBin++){
if(strcmp("pLambda",Name)==0){
printf(" k=%.2f MeV; C(k)=%.3f\n", CkMain->GetBinCenter(uBin), child->CkMain->Eval(CkMain->GetBinCenter(uBin)));
}
            CkChildMainFeed[WhichCk]->SetBinContent(uBin, child->CkMain->Eval(CkMain->GetBinCenter(uBin)));
        }

if(strcmp("pLambda",Name)==0){
printf("k=%.2f MeV; C(k)=%.3f\n", 130., CkChildMainFeed[WhichCk]->Eval(130));
}
*/
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

    double RESULT = MuPar*CkSmearedMainFeed->Eval(Momentum);
//printf(" MuPar=%f\n",MuPar);
//CATShisto<double>* CkMainSmeared = new CATShisto<double>(*CkMain);
//Smear(CkMain,RM_MomResolution,CkMainSmeared);
//printf(" RESULT=%f\n",RESULT);
//printf(" CkMain=%f <-> %f\n", CkMain->Eval(Momentum), CkMainSmeared->Eval(Momentum));
//printf(" CkFeed=%f <-> %f\n", Child[0]->CkMainFeed->Eval(Momentum), Child[0]->CkSmearedMainFeed->Eval(Momentum));
//printf(" CkMainFeed=%f <-> %f\n",CkMainFeed->Eval(Momentum), CkSmearedMainFeed->Eval(Momentum));
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Type[uChild]!=cFake) continue;
        if(Child[uChild]){
            RESULT += LambdaPar[uChild]*Child[uChild]->CkSmearedMainFeed->Eval(Momentum);
//printf("CHILD LambdaPar[%u]=%f\n",uChild,LambdaPar[uChild]);
        }
        else{
            RESULT += LambdaPar[uChild];
//printf("LambdaPar[%u]=%f\n",uChild,LambdaPar[uChild]);
        }
//printf("  RESULT=%f\n",RESULT);
    }
    return RESULT;
}

double DLM_CkDecomposition::EvalMain(const double& Momentum){
    return CkMain->Eval(Momentum);
}
double DLM_CkDecomposition::EvalSmearedMain(const double& Momentum){
    if(!CkMainSmeared){
        CkMainSmeared = new CATShisto<double>(*CkMain);
        Smear(CkMain,RM_MomResolution,CkMainSmeared);
    }
    return CkMainSmeared->Eval(Momentum);
}
double DLM_CkDecomposition::EvalMainFeed(const double& Momentum){
    return CkMainFeed->Eval(Momentum);
}
double DLM_CkDecomposition::EvalSmearedMainFeed(const double& Momentum){
    return CkSmearedMainFeed->Eval(Momentum);
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
CATShisto<double>* DLM_CkDecomposition::GetChildContribution(const unsigned& WhichChild){
    if(WhichChild>=NumChildren) return NULL;
    CATShisto<double>* Histo = new CATShisto<double> (*CkChildMainFeed[WhichChild]);
    Smear(CkChildMainFeed[WhichChild], RM_MomResolution, Histo);
    return Histo;
}
CATShisto<double>* DLM_CkDecomposition::GetChildContribution(const char* name){
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(!Child[uChild]) continue;
        if(strcmp(name,Child[uChild]->Name)==0){
            return GetChildContribution(uChild);
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

/*
double DLM_CkDecomposition::EvalMainCk(const double& Momentum, const bool& Renorm){
    CkSmearedMainFeed->Eval(Momentum)*(Renorm?1:MuPar);
}
double DLM_CkDecomposition::EvalFeedCk(const double& Momentum, const bool& Renorm){

}
double DLM_CkDecomposition::EvalFakeCk(const double& Momentum, const bool& Renorm){

}
double DLM_CkDecomposition::EvalSpecificCk(const unsigned& WhichCk, const double& Momentum, const bool& Renorm){

}
*/
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
//if(strcmp("pLambda",Name)==0) printf("UPDATE\n");
//printf("NAME: %s; ",Name);
    //1)
    CurrentStatus = CkMain->Status();
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
//printf("So far so bad %u...\n",uChild);
        if(Child[uChild]){
            //always false if FORCE_FULL_UPDATE==true
            CurrentStatusChild[uChild] = Status()*(!FORCE_FULL_UPDATE);
//printf("CSC at %u: %ix%i\n", uChild, int(Status()),int(!FORCE_FULL_UPDATE));
        }
        //if the child is not defined -> take a flat C(k)
        else if(!CurrentStatusChild[uChild]){
//printf("So far so good %u...\n",uChild);
            //for(unsigned uBin=0; uBin<CkSmearedChild[uChild]->GetNbins(); uBin++){
//printf("So far so good %u %u...\n",uChild,uBin);
                //CkSmearedChild[uChild]->SetBinContent(uBin,1);
            //}
            CurrentStatusChild[uChild] = true;
        }
    }

    //2)
    if(!CurrentStatus || FORCE_FULL_UPDATE){
        CkMain->Update(FORCE_FULL_UPDATE);
//if(FORCE_FULL_UPDATE && DEBUGFLAG==1){
    //printf("CkMain(10)=%.2f\n",CkMain->Eval(10));
//}

    }

    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if( !CurrentStatusChild[uChild] ){
            Child[uChild]->Update(FORCE_FULL_UPDATE);
        }
    }

//if(strcmp("pLambda",Name)==0){
//printf(" now Kid 0 is %p\n", 130., CkChildMainFeed[0]);
//if(CkChildMainFeed[0])
//printf(" Eval(%f) = %f\n", 130., CkChildMainFeed[0]->Eval(130));
//}

    //3)
    if(!CurrentStatus || FORCE_FULL_UPDATE){
//if(strcmp("pLambda",Name)==0){
//    printf(" CkMain->GetNbins()=%u; CkMainFeed->GetNbins()=%u\n",CkMain->GetNbins(),CkMainFeed->GetNbins());
//    printf(" MuPar = %f\n",MuPar);
//}


        CkMainFeed->Copy(CkMain[0], LambdaMain/MuPar);
//if(strcmp("pp",Name)==0) printf("At 50 MeV: CkMain->Eval(50)=%f\n", CkMain->Eval(50));
//if(strcmp("pp",Name)==0) printf("          : CkMainFeed->Eval(50)=%f\n", CkMainFeed->Eval(50));
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]!=cFeedDown) continue;
            if(Child[uChild]){
//if(strcmp("pp",Name)==0){
//SHOWSTUFF=true;
//printf(" BEFORE SMEAR: CkChildMainFeed[%u](50)=%.3f\n",uChild,CkChildMainFeed[uChild]->Eval(40));
//printf(" BEFORE SMEAR: Child[%u]->CkMainFeed(50)=%.3f\n",uChild,Child[uChild]->CkMainFeed->Eval(50));
//}
                Smear(Child[uChild]->CkMainFeed, RM_Child[uChild], CkChildMainFeed[uChild]);
//Smear(Child[uChild]->CkMainFeed, RM_Child[uChild], CkChildMainFeed[uChild]);
                //for(){
//if(strcmp("pp",Name)==0){
//printf(" AFTER SMEAR: CkChildMainFeed[%u](50)=%.3f\n",uChild,CkChildMainFeed[uChild]->Eval(50));
//printf("            : CkChildMainFeed[%u]=%p\n",uChild,CkChildMainFeed[uChild]);
//}

                //}

                for(unsigned uBin=0; uBin<CkMainFeed->GetNbins(); uBin++){
                    //if(strcmp("pp",Name)==0){
                    //    printf(" %u before add: %.2f\n", uBin, CkMainFeed->GetBinContent(uBin));
                    //}
                    CkMainFeed->Add(uBin, CkChildMainFeed[uChild]->Eval(CkMainFeed->GetBinCenter(uBin))*LambdaPar[uChild]/MuPar);
                    //CkMainFeed->Add(uBin, 2);
                    //if(strcmp("pp",Name)==0){
                    //    printf(" BIN %u -> k=%.0f; add %.3f\n",
                    //    uBin, CkMainFeed->GetBinCenter(uBin),
                    //    CkChildMainFeed[uChild]->Eval(CkMainFeed->GetBinCenter(uBin))*LambdaPar[uChild]/MuPar);
                    //    printf(" %u after add: %.2f\n", uBin, CkMainFeed->GetBinContent(uBin));
                    //}
                }
                //CkMainFeed->Add(CkChildMainFeed[uChild][0],LambdaPar[uChild]/MuPar);

                //if(!CkMainFeed->Add(CkChildMainFeed[uChild][0],LambdaPar[uChild]/MuPar)){
                //printf("BIG TROUBLE!\n");
                //}
//if(strcmp("pp",Name)==0) printf("--> uChild=%u;  LambdaMain/MuPar=%f\n", uChild, LambdaMain/MuPar);
//if(strcmp("pp",Name)==0) printf("          : CkMain=%f; CkChildMainFeed[uChild]=%f; LambdaPar[uChild]/MuPar=%f\n",
//       CkMain->Eval(50), CkChildMainFeed[uChild]->Eval(50), LambdaPar[uChild]/MuPar);
//if(strcmp("pp",Name)==0) printf("          : CkMainFeed=%f\n", CkMainFeed->Eval(50));
            }
            else{
                CkMainFeed->Add(LambdaPar[uChild]/MuPar);
            }
//SHOWSTUFF=false;
        }

        //Smear(CkMain, RM_MomResolution, CkMainSmeared);
        Smear(CkMainFeed, RM_MomResolution, CkSmearedMainFeed);
//if(strcmp("pp",Name)==0) printf("         A: CkMainFeed=%f\n", CkMainFeed->Eval(50));
//if(strcmp("pp",Name)==0) printf("         A: CkSmearedMainFeed=%f\n", CkSmearedMainFeed->Eval(50));
    }

    DecompositionStatus = true;

}

void DLM_CkDecomposition::Smear(const CATShisto<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, CATShisto<double>* CkSmeared){
//printf("CALL FOR SMEAR!\n");
    //double MomentumTrue;
    if(!SmearMatrix){
        CkSmeared[0] = CkToSmear[0];
//printf(" BAD ENDING!\n");
    }
    else{
        for(unsigned uBinSmear=0; uBinSmear<CkToSmear->GetNbins(); uBinSmear++){
            CkSmeared->SetBinContent(uBinSmear, 0);
//double THISBIN1=0;
//double THISBIN2=0;
            for(unsigned uBinTrue=0; uBinTrue<CkToSmear->GetNbins(); uBinTrue++){
                //MomentumTrue = CkToSmear->GetBinCenter(uBinTrue);
                CkSmeared->Add(uBinSmear, SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue));
                //CkSmeared->Add(uBinSmear, SmearMatrix->ResponseMatrix[uBinTrue][uBinSmear]*CkToSmear->GetBinContent(uBinTrue));
//THISBIN1+=SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue);
//THISBIN2+=SmearMatrix->ResponseMatrix[uBinTrue][uBinSmear]*CkToSmear->GetBinContent(uBinTrue);
//if(SHOWSTUFF) printf("%f vs %f\n",SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue],SmearMatrix->ResponseMatrix[uBinTrue][uBinSmear]);
//if(SHOWSTUFF) printf(" ADD1: ubs=%u; ubt=%u; rm=%f; ck=%f\n",uBinSmear,uBinTrue,SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue],CkToSmear->GetBinContent(uBinTrue));
//if(SHOWSTUFF) printf(" ADD2: ubs=%u; ubt=%u; rm=%f; ck=%f\n",uBinSmear,uBinTrue,SmearMatrix->ResponseMatrix[uBinTrue][uBinSmear],CkToSmear->GetBinContent(uBinTrue));
            }
//if(SHOWSTUFF) printf("k=%.0f (%.0f); THISBIN1=%.3f; THISBIN2=%.3f; Original=%.3f; CkSmeared=%.3f\n",
//                     CkToSmear->GetBinCenter(uBinSmear),CkSmeared->GetBinCenter(uBinSmear),
//                     THISBIN1,THISBIN2,CkToSmear->GetBinContent(uBinSmear),CkSmeared->GetBinContent(uBinSmear));
//if(SHOWSTUFF) printf(" CkSmeared=%p\n",CkSmeared);
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
