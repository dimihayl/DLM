
#include <stdio.h>
#include <string.h>

#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_CkModels.h"
#include "DLM_ResponseMatrix.h"
#include "DLM_Histo.h"



//the main contribution does not count as a child, i.e. 1 child implies one contribution additionally to the main one!
DLM_CkDecomp::DLM_CkDecomp(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const DLM_Histo<float>* hSigmaMatrix, const bool& InvertAxis):
    ERROR_STATE(!name),NumChildren(numchildren),CkMain(&ckfunction){

    StandardSetup(name);
    AddTheMatrix(hSigmaMatrix,InvertAxis);
    Update(true);

}
DLM_CkDecomp::DLM_CkDecomp(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction):
    ERROR_STATE(!name),NumChildren(numchildren),CkMain(&ckfunction){

    StandardSetup(name);

}
void DLM_CkDecomp::StandardSetup(const char* name){
    Child = NULL;
    LambdaPar = NULL;
    Type = NULL;
    RM_Child = NULL;
    SM_Child = NULL;
    RM_MomResolution = NULL;
    CurrentStatus = false;
    DecompositionStatus = false;
    CurrentStatusChild = NULL;
    CkMainSmeared = NULL;
    CkMainFeed = NULL;
    CkSmearedMainFeed = NULL;
    SignalMain = NULL;
    SignalChild = NULL;
    SignalSmearedMain = NULL;
    SignalSmearedChild = NULL;
    PS_Main = NULL;
    PS_Child = NULL;
    SignalsUpdated = false;
    qa_lam = true;

    if(ERROR_STATE){
        printf("\033[1;31mERROR:\033[0m The DLM_CkDecomp got some rubbish input, the object will be broken!\n");
        return;
    }

    CkMainFeed = new DLM_Histo<double>(CkMain[0]);
    CkSmearedMainFeed = new DLM_Histo<double>(CkMain[0]);
    SignalMain = new DLM_Histo<double>(CkMain[0]);
    SignalSmearedMain = new DLM_Histo<double>(CkMain[0]);

    Name = new char [strlen(name)+1];
    strcpy(Name,name);

    LambdaMain = new DLM_Histo<double>(*CkMain);
    LambdaMain->SetBinContentAll(1);
    MuPar = new DLM_Histo<double>(*LambdaMain);
    MuPar->SetBinContentAll(1);

    if(NumChildren){
        Child = new DLM_CkDecomp* [NumChildren];
        LambdaPar = new DLM_Histo<double>* [NumChildren];
        Type = new int [NumChildren];
        RM_Child = new DLM_ResponseMatrix* [NumChildren];
        SM_Child = new DLM_ResponseMatrix* [NumChildren];
        CurrentStatusChild = new bool [NumChildren];
        CkChildMainFeed = new DLM_Histo<double>* [NumChildren];
        CkSmearedChildMainFeed = new DLM_Histo<double>* [NumChildren];
        SignalChild = new DLM_Histo<double>* [NumChildren];
        SignalSmearedChild = new DLM_Histo<double>* [NumChildren];
        PS_Main = NULL;
        PS_Child = new DLM_Histo<float>* [NumChildren];
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            Child[uChild] = NULL;
            LambdaPar[uChild] = NULL;
            Type[uChild] = -1;
            RM_Child[uChild] = NULL;
            SM_Child[uChild] = NULL;
            CurrentStatusChild[uChild] = false;
            CkChildMainFeed[uChild] = NULL;
            CkSmearedChildMainFeed[uChild] = NULL;
            SignalChild[uChild] = new DLM_Histo<double>(CkMain[0]);
            SignalSmearedChild[uChild] = new DLM_Histo<double>(CkMain[0]);
            PS_Child[uChild] = NULL;
        }
    }

}

void DLM_CkDecomp::AddTheMatrix(const DLM_Histo<float>* hSigmaMatrix, const bool& InvertAxis){
    RM_MomResolution = new DLM_ResponseMatrix(*CkMain, hSigmaMatrix, NULL, InvertAxis);
}

DLM_CkDecomp::~DLM_CkDecomp(){

    if(ERROR_STATE) return;
    if(NumChildren){
        delete [] Child; Child=NULL;
        
        delete [] Type; Type=NULL;
        delete [] CurrentStatusChild; CurrentStatusChild=NULL;

        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(RM_Child[uChild]) delete RM_Child[uChild];
            if(SM_Child[uChild]) delete SM_Child[uChild];
            if(CkChildMainFeed[uChild]) delete CkChildMainFeed[uChild];
            if(CkSmearedChildMainFeed[uChild]) delete CkSmearedChildMainFeed[uChild];
            if(SignalChild[uChild]) delete SignalChild[uChild];
            if(SignalSmearedChild[uChild]) delete SignalSmearedChild[uChild];
            if(PS_Child[uChild]) delete PS_Child[uChild];
            if(LambdaPar[uChild]) delete LambdaPar[uChild];
        }
        delete [] LambdaPar; LambdaPar=NULL;
        delete [] RM_Child; RM_Child=NULL;
        delete [] SM_Child; SM_Child=NULL;
        delete [] CkChildMainFeed; CkChildMainFeed=NULL;
        delete [] CkSmearedChildMainFeed; CkSmearedChildMainFeed=NULL;
        delete [] SignalChild; SignalChild=NULL;
        delete [] SignalSmearedChild; SignalSmearedChild=NULL;
        delete [] PS_Child; PS_Child=NULL;
    }
    if(Name) {delete [] Name; Name=NULL;}
    if(RM_MomResolution) {delete RM_MomResolution; RM_MomResolution=NULL;}
    if(CkMainSmeared) {delete CkMainSmeared; CkMainSmeared=NULL;}
    if(CkMainFeed) {delete CkMainFeed; CkMainFeed=NULL;}
    if(CkSmearedMainFeed) {delete CkSmearedMainFeed; CkSmearedMainFeed=NULL;}
    if(SignalMain) {delete SignalMain; SignalMain=NULL;}
    if(SignalSmearedMain) {delete SignalSmearedMain; SignalSmearedMain=NULL;}
    if(PS_Main) {delete PS_Main; PS_Main=NULL;}
    if(LambdaMain) {delete LambdaMain; LambdaMain=NULL;}
    if(MuPar) {delete MuPar; MuPar=NULL;}
}

void DLM_CkDecomp::AddPhaseSpace(const unsigned& WhichCk, const DLM_Histo<float>* hPhaseSpace){
  if(ERROR_STATE) return;
  if(WhichCk>=NumChildren){
      return;
  }
  if(PS_Child[WhichCk]) {delete PS_Child[WhichCk]; PS_Child[WhichCk]=NULL;}
  if(!hPhaseSpace) return;
  PS_Child[WhichCk] = new DLM_Histo<float>(*hPhaseSpace);
  PS_Child[WhichCk]->ScaleToBinSize();
}

void DLM_CkDecomp::AddPhaseSpace(const DLM_Histo<float>* hPhaseSpace){
  if(ERROR_STATE) return;
  if(PS_Main) {delete PS_Main; PS_Main=NULL;}
  if(!hPhaseSpace) return;
  PS_Main = new DLM_Histo<float>(*hPhaseSpace);
  PS_Main->ScaleToBinSize();
}

void DLM_CkDecomp::AddContribution(unsigned WhichCk, DLM_Histo<double>& fraction, int type, DLM_CkDecomp* child,
                                          const DLM_Histo<float>* hResidualMatrix, const bool& InvertedAxis){
    if(ERROR_STATE) return;

    if(WhichCk>=NumChildren){
        return;
    }
    if(type!=cFeedDown && type!=cFake){
        return;
    }

    if(fraction.GetDim()!=1){
        printf("\033[1;33mWARNING:\033[0m The histogram for the lambda pars has to be 1D.\n");
        return;        
    }
    for(unsigned uMom=0; uMom<fraction.GetNbins(); uMom++){
        if(fraction.GetBinContent(uMom)<0 || fraction.GetBinContent(uMom)>1){
            printf("\033[1;33mWARNING:\033[0m Trying the set a contribution with a fraction<0 || fraction>1.\n");
            return;
        }
    }

    if(child && !UniqueName(child->Name)){
        printf("\033[1;33mWARNING:\033[0m Trying to duplicate the name '%s' of a contribution in DLM_CkDecomp::SetContribution.\n", child->Name);
        return;
    }

    if(LambdaPar[WhichCk]){
        if( Child[WhichCk]==child && *LambdaPar[WhichCk]==fraction && Type[WhichCk]==type){
            return;
        }
        delete LambdaPar[WhichCk];
    }

    if(type==cFake && hResidualMatrix){
        printf("\033[1;33mWARNING:\033[0m Impurities are not allowed to have a transformation matrix.\n");
        hResidualMatrix = NULL;
    }

    Child[WhichCk] = child;
    LambdaPar[WhichCk] = new DLM_Histo<double>(*LambdaMain);
    LambdaPar[WhichCk]->SetEqualTo(fraction);
    Type[WhichCk] = type;
    DecompositionStatus = false;

    if(RM_Child[WhichCk]) {delete RM_Child[WhichCk]; RM_Child[WhichCk]=NULL;}
    if(SM_Child[WhichCk]) {delete SM_Child[WhichCk]; SM_Child[WhichCk]=NULL;}
    if(child && hResidualMatrix){
        RM_Child[WhichCk] = new DLM_ResponseMatrix(*child->CkMain, NULL, hResidualMatrix, InvertedAxis);
    }
    if(child && RM_MomResolution){
        SM_Child[WhichCk] = new DLM_ResponseMatrix(*child->CkMain, RM_MomResolution->hSigmaMatrix, NULL, RM_MomResolution->InvertedAxis);
    }

    if(CkChildMainFeed[WhichCk]) {delete CkChildMainFeed[WhichCk]; CkChildMainFeed[WhichCk]=NULL;}
    if(child){
        CkChildMainFeed[WhichCk] = new DLM_Histo<double> (*child->CkMain);
    }

    LambdaMain->SetBinContentAll(1);
    MuPar->SetBinContentAll(1);
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(LambdaPar[uChild]){
            *LambdaMain -= *LambdaPar[uChild];
            if(Type[uChild]==cFake) *MuPar -= *LambdaPar[uChild];
        }
    }
}


void DLM_CkDecomp::AddContribution(unsigned WhichCk, double fraction, int type, DLM_CkDecomp* child,
                                          const DLM_Histo<float>* hResidualMatrix, const bool& InvertedAxis){
    DLM_Histo<double> fraction_histo(*LambdaMain);
    fraction_histo.SetBinContentAll(fraction);
    AddContribution(WhichCk, fraction_histo, type, child, hResidualMatrix, InvertedAxis);
}

//void DLM_CkDecomp::UnfoldData(const unsigned& nboot){
    //CkMain the main correlation function
    //DLM_CkDecomp** Child; list of all children
    //int* Type; //type of the children


//}

double DLM_CkDecomp::EvalCk(const double& Momentum){
    //printf("MOM %f\n",Momentum);
    if(ERROR_STATE) return 0;

    double MinLambdaPar=0;
    LambdaMain->FindMinima(MinLambdaPar);
    if(MinLambdaPar<0){
        printf("\033[1;31mERROR:\033[0m The λ0 parameter is smaller than zero!\n");
        return 0;
    }
    Update(false);
    double VAL_CkSmearedMainFeed=0;
    double VAL_CkSmearedMainFake=0;
    //flat towards unity
    if(Momentum>CkMain->GetCutOff()){
        double kf;
        if(CkMain->GetCutOff()>CkMain->GetUpEdge(0)){
            kf = CkMain->GetUpEdge(0);
        }
        //Momentum>CutOff
        else{
            kf = CkMain->GetCutOff();
        }
        const double Cf = EvalCk(CkMain->GetCutOff());
        double CutOff_kc = CkMain->GetCutOff_kc();
        if(CutOff_kc<0) return fabs(CutOff_kc);
        if(CutOff_kc<=kf) return Cf;
        double ReturnVal = ((Momentum-kf)-Cf*(Momentum-CutOff_kc))/(CutOff_kc-kf);
        if(ReturnVal>1) ReturnVal=1;
        return ReturnVal;
    }
    //this is needed in order to get a smooth transition for the range outside of the histogram
    else if(Momentum>CkMain->GetUpEdge(0)){
        VAL_CkSmearedMainFeed = CkMain->Eval(Momentum)*LambdaMain->Eval(Momentum)/MuPar->Eval(Momentum);
        //printf("CkMain->Eval(%.0f) = %.3f (l/mu = %.3f)\n",CkMain->Eval(Momentum),Momentum,LambdaMain->Eval(Momentum)/MuPar->Eval(Momentum));
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]==cFeedDown){
                if(Child[uChild]){
                    VAL_CkSmearedMainFeed+=GetChild(uChild)->GetCk()->Eval(Momentum)*LambdaPar[uChild]->Eval(Momentum)/MuPar->Eval(Momentum);
                }
                else if(LambdaPar[uChild]){
                    VAL_CkSmearedMainFeed+=(LambdaPar[uChild]->Eval(Momentum)/MuPar->Eval(Momentum));
                }
            }
            else{
                if(Child[uChild]){
                    VAL_CkSmearedMainFake+=GetChild(uChild)->GetCk()->Eval(Momentum)*LambdaPar[uChild]->Eval(Momentum);
                }
                else if(LambdaPar[uChild]){
                    VAL_CkSmearedMainFake+=(LambdaPar[uChild]->Eval(Momentum));
                }
            }
        }
    }
    else{
        VAL_CkSmearedMainFeed = CkSmearedMainFeed->Eval(Momentum);
        //printf("VAL_CkSmearedMainFeed(%.0f) = %.3f\n",Momentum,VAL_CkSmearedMainFeed);
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(Type[uChild]!=cFake) continue;
            if(Child[uChild]){
                VAL_CkSmearedMainFake += LambdaPar[uChild]->Eval(Momentum)*Child[uChild]->CkSmearedMainFeed->Eval(Momentum);
            }
            else if(LambdaPar[uChild]){
                VAL_CkSmearedMainFake += LambdaPar[uChild]->Eval(Momentum);
            }
        }
    }
    //printf("k*=%.0f; Mu=%.3f; CkSmMF=%.3f; CkSMM=%.3f\n",Momentum,MuPar->Eval(Momentum),VAL_CkSmearedMainFeed,VAL_CkSmearedMainFake);
    return MuPar->Eval(Momentum)*VAL_CkSmearedMainFeed+VAL_CkSmearedMainFake;
}

double DLM_CkDecomp::EvalMain(const double& Momentum){
    return CkMain->Eval(Momentum);
}
double DLM_CkDecomp::EvalSmearedMain(const double& Momentum){
    if(!CkMainSmeared){
        CkMainSmeared = new DLM_Histo<double>(*CkMain);
        Smear(CkMain,RM_MomResolution,CkMainSmeared,PS_Main);
    }
    return CkMainSmeared->Eval(Momentum);
}
double DLM_CkDecomp::EvalMainFeed(const double& Momentum){
    return CkMainFeed->Eval(Momentum);
}
double DLM_CkDecomp::EvalSmearedMainFeed(const double& Momentum){
    return CkSmearedMainFeed->Eval(Momentum);
}
double DLM_CkDecomp::EvalSmearedFeed(const unsigned& WhichChild,const double& Momentum){
    if(WhichChild>=NumChildren) return 0;
    if(!CkSmearedChildMainFeed[WhichChild]&&Child[WhichChild]) {
        CkSmearedChildMainFeed[WhichChild] = new DLM_Histo<double> (*Child[WhichChild]->CkMain);
        Smear(CkChildMainFeed[WhichChild], SM_Child[WhichChild], CkSmearedChildMainFeed[WhichChild], PS_Child[WhichChild]);
    }
    if(CkSmearedChildMainFeed[WhichChild]) return CkSmearedChildMainFeed[WhichChild]->Eval(Momentum);
    else return 0;
}

double DLM_CkDecomp::EvalSignal(const double& Momentum){
	double Result = 0;
	Result += EvalSignalMain(Momentum);
	for(unsigned uChild=0; uChild<NumChildren; uChild++) Result += EvalSignalChild(uChild,Momentum);
	return Result;
}
double DLM_CkDecomp::EvalSignalSmeared(const double& Momentum){
	double Result = 0;
	Result += EvalSignalSmearedMain(Momentum);
	for(unsigned uChild=0; uChild<NumChildren; uChild++) Result += EvalSignalSmearedChild(uChild,Momentum);
	return Result;
}
double DLM_CkDecomp::EvalSignalMain(const double& Momentum){
	if(SignalMain) return SignalMain->Eval(Momentum);
	return 0;
}
double DLM_CkDecomp::EvalSignalSmearedMain(const double& Momentum){
	if(SignalSmearedMain) return SignalSmearedMain->Eval(Momentum);
	return 0;
}
double DLM_CkDecomp::EvalSignalChild(const unsigned& WhichChild,const double& Momentum){
	if(WhichChild>=NumChildren) return 0;
	if(SignalChild[WhichChild]) return SignalChild[WhichChild]->Eval(Momentum);
	return 0;
}
double DLM_CkDecomp::EvalSignalSmearedChild(const unsigned& WhichChild,const double& Momentum){
	if(WhichChild>=NumChildren) return 0;
	if(SignalSmearedChild[WhichChild]) return SignalSmearedChild[WhichChild]->Eval(Momentum);
	return 0;
}

void DLM_CkDecomp::SetLambdaMain(const DLM_Histo<double>& lambda_par){
    double MinLambdaPar=0;
    lambda_par.FindMinima(MinLambdaPar);
    if(MinLambdaPar<0 || MinLambdaPar>1){
        printf("\033[1;33mWARNING:\033[0m Trying the set a contribution with a lambda_par<0 || lambda_par>1.\n");
        return;
    }
    if(*LambdaMain == lambda_par){
        return;
    }
    LambdaMain->SetEqualTo(lambda_par);
    DecompositionStatus = false;
}
void DLM_CkDecomp::SetLambdaMain(const double& lambda_par){
    DLM_Histo<double> lambda_histo(*LambdaMain);
    lambda_histo.SetBinContentAll(lambda_par);
    SetLambdaMain(lambda_histo);
}

void DLM_CkDecomp::SetLambdaChild(const unsigned& WhichChild, const DLM_Histo<double>& lambda_par){
    double MinLambdaPar=0;
    lambda_par.FindMinima(MinLambdaPar);
    if(MinLambdaPar<0 || MinLambdaPar>1){
        printf("\033[1;33mWARNING:\033[0m Trying the set a contribution with a lambda_par<0 || lambda_par>1.\n");
        return;
    }
    if(*LambdaPar[WhichChild]==lambda_par){
        return;
    }
  LambdaPar[WhichChild]->SetEqualTo(lambda_par);
  DecompositionStatus = false;
}
void DLM_CkDecomp::SetLambdaChild(const unsigned& WhichChild, const double& lambda_par){
    DLM_Histo<double> lambda_histo(*LambdaMain);
    lambda_histo.SetBinContentAll(lambda_par);
    SetLambdaChild(WhichChild, lambda_histo);
}

void DLM_CkDecomp::QA_LambdaPar(const bool& yesno){
  qa_lam = yesno;
}

double DLM_CkDecomp::GetLambdaMain(const double Momentum){
	return LambdaMain->Eval(Momentum);
}
double DLM_CkDecomp::GetLambdaChild(const unsigned& WhichChild, const double Momentum){
    if(WhichChild>=NumChildren) return 0;
    if(!LambdaPar[WhichChild]) return 0;
    return LambdaPar[WhichChild]->Eval(Momentum);
}
unsigned DLM_CkDecomp::GetNumChildren(){
    return NumChildren;
}

DLM_CkDecomp* DLM_CkDecomp::GetChild(const unsigned& WhichChild){
    if(WhichChild>=NumChildren) return NULL;
    else return Child[WhichChild];
}
DLM_CkDecomp* DLM_CkDecomp::GetContribution(const char* name){
    DLM_CkDecomp* RESULT = NULL;
    if(strcmp(name,Name)==0) return this;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Child[uChild]) RESULT = Child[uChild]->GetContribution(name);
        if(RESULT) break;
    }
    return RESULT;
}
DLM_Histo<double>* DLM_CkDecomp::GetChildContribution(const unsigned& WhichChild, const bool& WithLambda){
    if(WhichChild>=NumChildren) return NULL;
    if(!Child[WhichChild]) return NULL;
    DLM_Histo<double>* Histo = new DLM_Histo<double> (*CkChildMainFeed[WhichChild]);
    if(Type[WhichChild]!=cFake){
        Histo->SetEqualTo(*CkChildMainFeed[WhichChild]);
        if(WithLambda){
            Histo->MultiplyHisto(*LambdaPar[WhichChild]);
        }
        Smear(Histo, SM_Child[WhichChild], Histo, PS_Child[WhichChild]);
    }
    return Histo;
}
DLM_Histo<double>* DLM_CkDecomp::GetChildContribution(const char* name, const bool& WithLambda){
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(!Child[uChild]) continue;
        if(strcmp(name,Child[uChild]->Name)==0){
            return GetChildContribution(uChild,WithLambda);
        }
    }
    return NULL;
}
DLM_Ck* DLM_CkDecomp::GetCk(){
    return CkMain;
}

const DLM_Histo<float>* DLM_CkDecomp::GetResolutionMatrix(){
    return RM_MomResolution->hSigmaMatrix;
}

void DLM_CkDecomp::GetName(char* name){
    if(!name) return;
    strcpy(name,Name);
}

bool DLM_CkDecomp::Status(){
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
void DLM_CkDecomp::Update(const bool& FORCE_FULL_UPDATE, const bool& UpdateDecomp){
    double Momentum;
    //check the lambda pars
    if(!DecompositionStatus && qa_lam){
      DLM_Histo<double> TotLambda(*LambdaMain);
      for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(LambdaPar[uChild]) TotLambda += *LambdaPar[uChild];
      }
      for(unsigned uMom=0; uMom<TotLambda.GetNbins(); uMom++){
        if(fabs(TotLambda.GetBinContent(uMom)-1.)>1e-6){
            printf("\033[1;33mWARNING:\033[0m The lambda parameters sum to %f. Rescaling them such as to add to unity!\n",TotLambda);
            LambdaMain->SetBinContent(uMom, LambdaMain->GetBinContent(uMom)/TotLambda.GetBinContent(uMom));
            for(unsigned uChild=0; uChild<NumChildren; uChild++){
                if(LambdaPar[uChild])
                    LambdaPar[uChild]->SetBinContent(uMom, LambdaPar[uChild]->GetBinContent(uMom)/TotLambda.GetBinContent(uMom));
            }
        }
      }
    }


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
        CkMainFeed->MultiplyHisto(*LambdaMain);
        CkMainFeed->DivideHisto(*MuPar);

        if(UpdateDecomp){
			SignalMain[0] = CkMain[0];
			SignalMain->MultiplyHisto(*LambdaMain);
		 }
		else{
        	SignalMain->SetBinContentAll(0);
        	SignalSmearedMain->SetBinContentAll(0);
		}

        for(unsigned uChild=0; uChild<NumChildren; uChild++){
			SignalChild[uChild]->SetBinContentAll(0);
			SignalSmearedChild[uChild]->SetBinContentAll(0);
            if(Type[uChild]!=cFeedDown){
                if(Child[uChild]){
                    for(unsigned uBin=0; uBin<SignalSmearedChild[uChild]->GetNbins(); uBin++){
                        Momentum = SignalSmearedChild[uChild]->GetBinCenter(0,uBin);
                        SignalChild[uChild]->SetBinContent(uBin,Child[uChild]->CkMainFeed->Eval(&Momentum));
                        if(LambdaPar[uChild])
                            SignalSmearedChild[uChild]->SetBinContent(uBin,LambdaPar[uChild]->Eval(Momentum)*Child[uChild]->CkSmearedMainFeed->Eval(&Momentum));
                    }
                    if(LambdaPar[uChild]){
                        SignalChild[uChild][0] -= *LambdaPar[uChild];
                        SignalSmearedChild[uChild][0] -= *LambdaPar[uChild];
                    }
                }
                continue;
            }
            if(Child[uChild]){
                Smear(Child[uChild]->CkMainFeed, RM_Child[uChild], CkChildMainFeed[uChild], PS_Child[uChild]);
                for(unsigned uBin=0; uBin<CkMainFeed->GetNbins(); uBin++){
                    Momentum = CkMainFeed->GetBinCenter(0,uBin);
                    if(LambdaPar[uChild]){
                        CkMainFeed->Add(uBin, CkChildMainFeed[uChild]->Eval(&Momentum)*LambdaPar[uChild]->Eval(Momentum)/MuPar->Eval(Momentum));
                        SignalChild[uChild]->SetBinContent(uBin,CkChildMainFeed[uChild]->Eval(&Momentum)*LambdaPar[uChild]->Eval(Momentum));
                    }
                }
                if(UpdateDecomp){
                    Smear(SignalChild[uChild], RM_MomResolution, SignalSmearedChild[uChild], PS_Main);
                    if(LambdaPar[uChild]){
                        SignalChild[uChild][0] -= *LambdaPar[uChild];
                        SignalSmearedChild[uChild][0] -= *LambdaPar[uChild];
                    }
                }
            }
            else{
                *CkMainFeed += (*LambdaPar[uChild] / *MuPar);
            }
        }
        Smear(CkMainFeed, RM_MomResolution, CkSmearedMainFeed, PS_Main);
        if(UpdateDecomp){
			Smear(SignalMain, RM_MomResolution, SignalSmearedMain, PS_Main);
			SignalMain[0] -= *LambdaMain;
			SignalSmearedMain[0] -= *LambdaMain;
		}
    }
    DecompositionStatus = true;
    SignalsUpdated = UpdateDecomp;    
}

void DLM_CkDecomp::Smear(const DLM_Histo<double>* CkToSmear, const DLM_ResponseMatrix* SmearMatrix, DLM_Histo<double>* CkSmeared, DLM_Histo<float>* PhaseSpace){
    if(!SmearMatrix){
        CkSmeared[0] = CkToSmear[0];
    }
    else{
        for(unsigned uBinSmear=0; uBinSmear<CkToSmear->GetNbins(); uBinSmear++){
            CkSmeared->SetBinContent(uBinSmear, 0);
            unsigned FirstBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisFirst];
            unsigned LastBin = SmearMatrix->SparseResponse[uBinSmear][DLM_ResponseMatrix::xAxisLast];
            if(LastBin>=CkToSmear->GetNbins()) LastBin = CkToSmear->GetNbins()-1;
            double PSPACE_NORM = PhaseSpace?0:1;
            for(unsigned uBinTrue=FirstBin; uBinTrue<=LastBin; uBinTrue++){
                //as the response matrix is normalized to the size of the bin, during the integration we multiply for it
                double TRANSF = SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue];
                double CKTRUE = CkToSmear->GetBinContent(uBinTrue);
                double DTRUE = CkToSmear->GetBinSize(uBinTrue);
                double DSMEAR = CkSmeared->GetBinSize(uBinSmear);
                //the phase space is assumed to be counts normalized to the bin width (i.e. function)
                double TrueCrd = CkToSmear->GetBinCenter(0,uBinTrue);
                //unsigned PSBIN = PhaseSpace->GetBin(0,TrueCrd);
                //double PSVOL = 1./PhaseSpace->GetBinSize(PSBIN);
                //this would be converting counts to functional value
                double PSPACE = PhaseSpace?PhaseSpace->Eval(TrueCrd):1;
                //CkSmeared->Add(uBinSmear,   SmearMatrix->ResponseMatrix[uBinSmear][uBinTrue]*CkToSmear->GetBinContent(uBinTrue)*
                //                            CkToSmear->GetBinSize(uBinTrue)*CkSmeared->GetBinSize(uBinSmear));
                CkSmeared->Add(uBinSmear,   TRANSF*PSPACE*CKTRUE*DTRUE*DSMEAR);
                PSPACE_NORM += TRANSF*PSPACE*DTRUE*DSMEAR;
            }
            if(PhaseSpace&&PSPACE_NORM){
              CkSmeared->ScaleBin(uBinSmear,1./PSPACE_NORM);
            }
        }
    }
}

//check if any of the children has the same name
bool DLM_CkDecomp::UniqueName(const char* name){
    if(ERROR_STATE) return 0;
    if(strcmp(Name,name)==0) return false;
    for(unsigned uChild=0; uChild<NumChildren; uChild++){
        if(Child[uChild]){
            if(Child[uChild]->UniqueName(name)==false) return false;
        }
    }
    return true;
}
