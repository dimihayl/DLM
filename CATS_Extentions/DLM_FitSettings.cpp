#include "DLM_FitSettings.h"
#include "DLM_CkDecomposition.h"
#include "CATS.h"
#include "DLM_Source.h"
#include "DLM_Potentials.h"
//#include "DLM_CATSresults2.h"
#include "DLM_CkModels.h"
#include "DLM_WfModel.h"

#include <iostream>
//#include <vector>

using namespace std;

FemtoParticle::FemtoParticle(const TString& name, const unsigned& numfeed, const unsigned& nummisid):
    Name(name),NumFeedChannels(numfeed),NumMisidChannels(nummisid){
    FeedParticle = NULL; if(NumFeedChannels) FeedParticle = new FemtoParticle* [NumFeedChannels];
    for(unsigned uFeed=0; uFeed<NumFeedChannels; uFeed++) FeedParticle[uFeed]=NULL;
    //printf("%s: FeedParticle=%p\n",name.Data(),FeedParticle[0]);
    MisidParticle = NULL; if(NumMisidChannels) MisidParticle = new FemtoParticle* [NumMisidChannels];
    for(unsigned uMis=0; uMis<NumMisidChannels; uMis++) MisidParticle[uMis]=NULL;
    //Mass = 0;
}

FemtoParticle::~FemtoParticle(){
    if(FeedParticle){delete[]FeedParticle; FeedParticle=NULL;}
    if(MisidParticle){delete[]MisidParticle; MisidParticle=NULL;}
}

void FemtoParticle::SetFeed(const unsigned& WhichOne, FemtoParticle& Particle){
    if(WhichOne>=NumFeedChannels) return;
    if(!NameIsAvailable(Particle.Name)) {printf("WARNING: The name %s is taken.\n",Particle.Name.Data()); return;}
    //this makes sure that the Particle is not a "parent" of this
    if(!Particle.NameIsAvailable(Name)) {printf("WARNING: The name %s is taken.\n",Particle.Name.Data()); return;}
    FeedParticle[WhichOne] = &Particle;
}

void FemtoParticle::SetMisid(const unsigned& WhichOne, FemtoParticle& Particle){
    if(WhichOne>=NumMisidChannels) return;
    if(!NameIsAvailable(Particle.Name)) {printf("WARNING: The name %s is taken.\n",Particle.Name.Data()); return;}
    if(!Particle.NameIsAvailable(Name)) {printf("WARNING: The name %s is taken.\n",Particle.Name.Data()); return;}
    MisidParticle[WhichOne] = &Particle;
}

bool FemtoParticle::NameIsAvailable(const TString& name) const{
    if(name==Name) return false;
    for(unsigned uFeed=0; uFeed<NumFeedChannels; uFeed++){
        if(!FeedParticle[uFeed]) continue;
        if(FeedParticle[uFeed]->NameIsAvailable(name)==false) return false;
    }
    for(unsigned uMis=0; uMis<NumMisidChannels; uMis++){
        if(!MisidParticle[uMis]) continue;
        if(MisidParticle[uMis]->NameIsAvailable(name)==false) return false;
    }
    return true;
}
TString FemtoParticle::GetName() const{
    return Name;
}
unsigned FemtoParticle::GetNumFeed() const{
    return NumFeedChannels;
}
unsigned FemtoParticle::GetNumContrib() const{
    return GetNumContrib(true)+1;//the plus one is self counting
}
unsigned FemtoParticle::GetNumContrib(const bool& IncludeMisid) const{
    unsigned Result=NumFeedChannels;
    for(unsigned uFeed=0; uFeed<NumFeedChannels; uFeed++){
        if(!FeedParticle[uFeed]) continue;
        Result += FeedParticle[uFeed]->GetNumContrib(false);
    }
    if(IncludeMisid){
        for(unsigned uMis=0; uMis<NumMisidChannels; uMis++){
            if(!MisidParticle[uMis]) continue;
            Result += MisidParticle[uMis]->GetNumContrib(false);
        }
    }
    return Result;
}

unsigned FemtoParticle::GetNumMisid() const{
    return NumMisidChannels;
}

/*
//this function makes sure we assign a unique id to each feed-down (iteratively)
int FemtoParticle::GetFeedID(const TString& name) const{
    if(name==Name) return 0;
    int Result=1;
    int Result2;
    for(unsigned uFeed=0; uFeed<NumFeedChannels; uFeed++){
        if(!FeedParticle[uFeed]) {Result++; continue;}
        Result2 = FeedParticle[uFeed]->GetFeedID(name);
        if(Result2>=0) return Result+Result2;
        Result += FeedParticle[uFeed]->GetNumContrib();
    }
    return -1;
}
*/
//this function makes sure we assign a unique id to each feed-down
int FemtoParticle::GetFeedID(const TString& name) const{
    if(name==Name) return 0;
    for(unsigned uFeed=0; uFeed<NumFeedChannels; uFeed++){
        if(!FeedParticle[uFeed]) {continue;}
        if(FeedParticle[uFeed]->GetName()==name) return uFeed+1;
    }
    return -1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FemtoExperiment::FemtoExperiment(const TString& name, const unsigned& numpart):
    Name(name),NumParticles(numpart?numpart:1){
    if(!numpart){
        printf("ERROR: You should have at least one particle in your FemtoExperiment\n");
        printf("       NumParticles automatically set to 1!\n");
    }
    Particle = new const FemtoParticle* [NumParticles];
    Purity = new double [NumParticles];
    Fraction = new double* [NumParticles];
    for(unsigned uPart=0; uPart<NumParticles; uPart++){
        Fraction[uPart] = NULL;
    }

}
FemtoExperiment::~FemtoExperiment(){
    delete [] Particle; Particle=NULL;
    delete [] Purity; Purity=NULL;
    for(unsigned uPart=0; uPart<NumParticles; uPart++){
        if(Fraction[uPart]) {delete [] Fraction[uPart]; Fraction[uPart]=NULL;}
    }
    delete [] Fraction; Fraction=NULL;
}

void FemtoExperiment::SetParticle(const unsigned& WhichOne, const FemtoParticle& particle){
    if(WhichOne>=NumParticles) return;
    for(unsigned uPart=0; uPart<NumParticles; uPart++){
        if(&particle==Particle[uPart]) {printf("WARNING: Duplicate particles!\n"); return;}
        if(Particle[uPart]){
            if(particle.GetName()==Particle[uPart]->GetName()) {printf("WARNING: Duplicate particles!\n"); return;}
        }
    }
    Particle[WhichOne] = &particle;
    if(Fraction[WhichOne]) {delete[]Fraction[WhichOne]; Fraction[WhichOne]=NULL;}
    unsigned NumFrac = particle.GetNumFeed()+1;
    if(NumFrac){
        Fraction[WhichOne] = new double [NumFrac];
        for(unsigned uFeed=0; uFeed<NumFrac; uFeed++){
            //Fraction[WhichOne][uFeed] = 1./double(NumFeed);
            Fraction[WhichOne][uFeed] = 0;
        }
    }
}

void FemtoExperiment::SetPurity(const unsigned& WhichOne, const double& pur){
    if(WhichOne>=NumParticles) return;
    if(pur<0 || pur>1){
        printf("WARNING: The purity should be between 0 and 1\n");
        return;
    }
    Purity[WhichOne] = pur;
}

void FemtoExperiment::SetPurity(const TString& WhichOne, const double& pur){
    //if(WhichOne>=NumParticles) return;
    unsigned WichParticle = FindParticle(WhichOne);
    if(WichParticle<0) return;
    SetPurity(WichParticle,pur);
}

TString FemtoExperiment::GetName() const{
    return Name;
}

void FemtoExperiment::SetFraction(const unsigned& WhichMother, const unsigned& WhichFeed, const double& feed){
    if(WhichMother>=NumParticles) return;
    if(WhichFeed>=Particle[WhichMother]->GetNumContrib()) return;
    Fraction[WhichMother][WhichFeed] = feed;
}

void FemtoExperiment::SetFraction(const TString& WhichMother, const TString& WhichFeed, const double& feed){
    unsigned MomID = FindParticle(WhichMother);
    if(MomID<0) return;
    int FeedID = Particle[MomID]->GetFeedID(WhichFeed);
    if(FeedID<0) {printf("WARNING: %s is not a direct feed-down channel onto %s\n",WhichFeed.Data(),WhichMother.Data());}
    SetFraction(MomID, FeedID, feed);
}

int FemtoExperiment::FindParticle(const TString& WhichOne)  const{
    for(unsigned uPart=0; uPart<NumParticles; uPart++){
        if(!Particle[uPart]) continue;
        if(Particle[uPart]->GetName()==WhichOne) {return uPart;}
    }
    printf("WARNING: The particle %s was not found! A crash might follow!\n", WhichOne.Data());
    return -1;
}

FemtoPair::FemtoPair(const FemtoParticle& part1, const FemtoParticle& part2):
    Particle1(part1),Particle2(part2),NumFeed1(Particle1.GetNumFeed()+1),NumFeed2(Particle2.GetNumFeed()+1){
    hFeedDown = NULL;
    hFeedDown = new TH2F** [NumFeed1];
    for(unsigned uFeed1=0; uFeed1<NumFeed1; uFeed1++){
        hFeedDown[uFeed1] = new TH2F* [NumFeed2];
        for(unsigned uFeed2=0; uFeed2<NumFeed2; uFeed2++){
            hFeedDown[uFeed1][uFeed2] = NULL;
        }
    }
    //basically later on we treat the particles to be the same type in case &Particle1==&Particle2
    if(Particle1.GetName()==Particle2.GetName() && &Particle1!=&Particle2){
        printf("ERROR: Two different FemtoPair objects share the same name! This is not allowed, the FemtoPair object will be broken!\n");
    }
}
FemtoPair::~FemtoPair(){
    for(unsigned uFeed1=0; uFeed1<NumFeed1; uFeed1++){
        delete [] hFeedDown[uFeed1];
    }
    delete [] hFeedDown; hFeedDown=NULL;
}

bool FemtoPair::IsIt(const TString& Name1, const TString& Name2) const{
    if(Particle1.GetName()==Name1 && Particle2.GetName()==Name2) return true;
    if(Particle2.GetName()==Name1 && Particle1.GetName()==Name2) return true;
    return false;
}

void FemtoPair::SetFeedDownMatrix(const TString& Mother1, const TString& Mother2,
                           const TString& Child1, const TString& Child2,
                           TH2F* hfeed){

    if(Child1!=Particle1.GetName() && Child1!=Particle2.GetName()){
        printf("The particle %s cannot be found!\n",Child1.Data());
        return;
    }
    if(Child2!=Particle1.GetName() && Child2!=Particle2.GetName()){
        printf("Particle1.GetName()=%s\n",Particle1.GetName().Data());
        printf("Particle2.GetName()=%s\n",Particle2.GetName().Data());
        printf("The particle %s cannot be found!\n",Child2.Data());
        return;
    }

    //const TString& CHILD1 = Child1==Particle1.GetName()?Child1:Child2;
    //const TString& CHILD2 = Child2==Particle2.GetName()?Child2:Child1;
    const TString& MOM1 = Child1==Particle1.GetName()?Mother1:Mother2;
    const TString& MOM2 = Child2==Particle2.GetName()?Mother2:Mother1;

    int FeedID1 = Particle1.GetFeedID(MOM1);
    if(FeedID1<0){
        printf("The feed-down %s -> %s cannot be found!\n",MOM1.Data(),Particle1.GetName().Data());
        return;
    }
    int FeedID2 = Particle2.GetFeedID(MOM2);
    if(FeedID2<0){
        printf("The feed-down %s -> %s cannot be found!\n",MOM2.Data(),Particle2.GetName().Data());
        return;
    }
    hFeedDown[FeedID1][FeedID2] = hfeed;
    if(&Particle1==&Particle2) hFeedDown[FeedID2][FeedID1] = hfeed;
}
void FemtoPair::SetFeedDownMatrix(const TString& Mother1, const TString& Mother2,
                           TH2F* hfeed){
    SetFeedDownMatrix(Mother1,Mother2,Particle1.GetName(),Particle2.GetName(),hfeed);
}
const FemtoParticle* FemtoPair::GetParticle(const TString& NAME) const{
    if(Particle1.GetName()==NAME) return &Particle1;
    if(Particle2.GetName()==NAME) return &Particle2;
    else return NULL;
}
const FemtoParticle* FemtoPair::GetParticle(const unsigned& WhichOne) const{
    if(WhichOne==0) return &Particle1;
    if(WhichOne==1) return &Particle2;
    return NULL;
}

FemtoExpPair::FemtoExpPair(const FemtoPair& ppair, const FemtoExperiment& experiment):
    PartPair(ppair),Experiment(experiment){

    //NumFemtoRegion=0;
    //NumBaselineRegion=0;
    //NumCorrelationModel=0;

    FemtoRegionLow=0;
    FemtoRegionUp=100;
    BlRegionLow=300;
    BlRegionUp=400;

    Kitty=NULL;
    CorrelationModel=NULL;
    //CatSourcePars = new double [32];
    CatSourcePars = new CATSparameters(CATSparameters::tSource,32,true);
    CatPotPars = new CATSparameters* [32];
    for(int i=0; i<32; i++){
        CatPotPars[i] = new CATSparameters(CATSparameters::tPotential,32,true);
    }
    WaveFunctionU = NULL;
    PhaseShifts = NULL;
    RadBins = NULL;

    NumRadBins=0;

    hResolution=NULL;
    //hFeedDown=NULL;
    hData=NULL;
    //hSystErr=NULL;

}

FemtoExpPair::~FemtoExpPair(){
    if(NumRadBins) CleanHaidenbauer(*Kitty,&WaveFunctionU,&PhaseShifts,&RadBins);
    if(Kitty){delete Kitty; Kitty=NULL;}
    if(CorrelationModel){delete CorrelationModel; CorrelationModel=NULL;}
    //if(hFeedDown){delete [] hFeedDown; hFeedDown=NULL;}
    delete CatSourcePars;
    for(int i=0; i<32; i++){
        delete CatPotPars[i];
    }
    delete [] CatPotPars;

}

void FemtoExpPair::SetFemtoRegion(const double& MIN, const double& MAX){
    if(MIN>=MAX){
        printf("WARNING: Bad FemtoRegion input!\n");
        return;
    }
    FemtoRegionLow=MIN;
    FemtoRegionUp=MAX;
}

void FemtoExpPair::SetBlRegion(const double& MIN, const double& MAX){
    if(MIN>=MAX){
        printf("WARNING: Bad BlRegion input!\n");
        return;
    }
    BlRegionLow=MIN;
    BlRegionUp=MAX;
}

//void FemtoExpPair::AddCorrelationModel(DLM_Ck& CorrMod){
//    CorrelationModel.push_back(&CorrMod);
//}

void FemtoExpPair::SetData(TH1F* data){
    hData=data;
    //hSystErr=systerr;
}

void FemtoExpPair::SetResolutionMatrix(TH2F* hres){
    hResolution = hres;
}
void FemtoExpPair::SetLargeMomCk(const int& TYPE){
    LargeMomCk = TYPE;
}
void FemtoExpPair::SetBlType(const int& TYPE){
    BlType = TYPE;
}

void FemtoExpPair::SetStandardInteraction(const TString& Inter, TH1F* TemplateData){
    //DLM_Ck* Correlation;
    //CATS* Cat;
    TH1F* histoDATA = TemplateData;
    if(!histoDATA) histoDATA = hData;
    if(!histoDATA) {printf("\033[1;31mERROR:\033[0m SetStandardInteraction needs a data template to figure out the binning\n"); return;}
//printf("hData = %p\n",hData);
    unsigned FirstBin = histoDATA->FindBin(FemtoRegionLow);
    unsigned LastBin = histoDATA->FindBin(FemtoRegionUp);
    unsigned NumMomBins = LastBin-FirstBin+1;
    double* MomBins = new double [NumMomBins+1];
    for(unsigned uBin=FirstBin; uBin<=LastBin; uBin++){
        MomBins[uBin-FirstBin] = histoDATA->GetBinLowEdge(uBin);
    }
    MomBins[NumMomBins] = histoDATA->GetXaxis()->GetBinUpEdge(LastBin);

    CatSourcePars->SetParameter(0,1.5,false);

    //if statements for the standard interactions
    if(PartPair.IsIt("Proton","Proton")){
        if(Inter=="CATS_AV18"){

            StandardInter = Inter;

            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);

            Kitty->SetExcludeFailedBins(false);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(4);
            Kitty->SetNumPW(0,2);
            Kitty->SetNumPW(1,2);
            Kitty->SetNumPW(2,2);
            Kitty->SetNumPW(3,2);
            Kitty->SetSpin(0,0);
            Kitty->SetSpin(1,1);
            Kitty->SetSpin(2,1);
            Kitty->SetSpin(3,1);
            Kitty->SetChannelWeight(0, 3./12.);//1S0
            Kitty->SetChannelWeight(1, 1./12.);//3P0
            Kitty->SetChannelWeight(2, 3./12.);//3P1
            Kitty->SetChannelWeight(3, 5./12.);//3P2

            Kitty->SetQ1Q2(1);
            Kitty->SetPdgId(2212, 2212);
            const double Mass_p=938.272;
            Kitty->SetRedMass( 0.5*Mass_p );

            CatPotPars[0]->SetParameter(0,NN_AV18);
            CatPotPars[0]->SetParameter(1,v18_Coupled3P2);
            CatPotPars[0]->SetParameter(2,1);
            CatPotPars[0]->SetParameter(3,1);
            CatPotPars[0]->SetParameter(4,1);
            CatPotPars[0]->SetParameter(5,0);
            CatPotPars[0]->SetParameter(6,0);
            CatPotPars[0]->SetParameter(7,0);

            CatPotPars[1]->SetParameter(0,NN_AV18);
            CatPotPars[1]->SetParameter(1,v18_Coupled3P2);
            CatPotPars[1]->SetParameter(2,1);
            CatPotPars[1]->SetParameter(3,1);
            CatPotPars[1]->SetParameter(4,1);
            CatPotPars[1]->SetParameter(5,1);
            CatPotPars[1]->SetParameter(6,1);
            CatPotPars[1]->SetParameter(7,0);

            CatPotPars[2]->SetParameter(0,NN_AV18);
            CatPotPars[2]->SetParameter(1,v18_Coupled3P2);
            CatPotPars[2]->SetParameter(2,1);
            CatPotPars[2]->SetParameter(3,1);
            CatPotPars[2]->SetParameter(4,1);
            CatPotPars[2]->SetParameter(5,1);
            CatPotPars[2]->SetParameter(6,1);
            CatPotPars[2]->SetParameter(7,1);

            CatPotPars[3]->SetParameter(0,NN_AV18);
            CatPotPars[3]->SetParameter(1,v18_Coupled3P2);
            CatPotPars[3]->SetParameter(2,1);
            CatPotPars[3]->SetParameter(3,1);
            CatPotPars[3]->SetParameter(4,1);
            CatPotPars[3]->SetParameter(5,1);
            CatPotPars[3]->SetParameter(6,1);
            CatPotPars[3]->SetParameter(7,2);

            Kitty->SetShortRangePotential(0,0,fDlmPot,*CatPotPars[0]);
            Kitty->SetShortRangePotential(1,1,fDlmPot,*CatPotPars[1]);
            Kitty->SetShortRangePotential(2,1,fDlmPot,*CatPotPars[2]);
            Kitty->SetShortRangePotential(3,1,fDlmPot,*CatPotPars[3]);

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);

        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else if(PartPair.IsIt("Proton","Lambda")){
        if(Inter=="CATS_NLO"){

            StandardInter = Inter;

            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(2);
            Kitty->SetNumPW(0,1);
            Kitty->SetNumPW(1,1);
            Kitty->SetSpin(0,0);
            Kitty->SetSpin(1,1);
            Kitty->SetChannelWeight(0, 0.25);
            Kitty->SetChannelWeight(1, 0.75);

            Kitty->SetQ1Q2(0);
            Kitty->SetPdgId(2212, 3122);
            const double Mass_p=938.272; const double Mass_L=1115.683;
            Kitty->SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

            InitHaidenbauerNLO("/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/CorrelationFiles_2018/Haidenbauer/pLambdaNLO/",Kitty[0],&WaveFunctionU,&PhaseShifts,&RadBins,NumRadBins);
            printf("Please keep in mind that InitHaidenbauerNLO works only for specific binning configurations!\n");

            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                Kitty->UseExternalWaveFunction(uBin,0,0,WaveFunctionU[uBin][0][0], NumRadBins, RadBins, PhaseShifts[uBin][0][0]);
                Kitty->UseExternalWaveFunction(uBin,1,0,WaveFunctionU[uBin][1][0], NumRadBins, RadBins, PhaseShifts[uBin][1][0]);
            }

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);

        }
        else if(Inter=="CATS_Usmani"){

            StandardInter = Inter;

            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);

            Kitty->SetExcludeFailedBins(false);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(2);
            Kitty->SetNumPW(0,1);
            Kitty->SetNumPW(1,1);
            Kitty->SetSpin(0,0);
            Kitty->SetSpin(1,1);
            Kitty->SetChannelWeight(0, 0.25);
            Kitty->SetChannelWeight(1, 0.75);

            Kitty->SetQ1Q2(0);
            Kitty->SetPdgId(2212, 3122);
            const double Mass_p=938.272; const double Mass_L=1115.683;
            Kitty->SetRedMass( (Mass_p*Mass_L)/(Mass_p+Mass_L) );

            CatPotPars[0]->SetParameter(0,pL_UsmaniOli);
            CatPotPars[0]->SetParameter(1,0);
            CatPotPars[0]->SetParameter(2,0);
            CatPotPars[0]->SetParameter(3,0);
            CatPotPars[0]->SetParameter(4,0);
            CatPotPars[0]->SetParameter(5,0);
            CatPotPars[0]->SetParameter(6,0);
            CatPotPars[0]->SetParameter(7,0);

            CatPotPars[1]->SetParameter(0,pL_UsmaniOli);
            CatPotPars[1]->SetParameter(1,0);
            CatPotPars[1]->SetParameter(2,0);
            CatPotPars[1]->SetParameter(3,0);
            CatPotPars[1]->SetParameter(4,0);
            CatPotPars[1]->SetParameter(5,1);
            CatPotPars[1]->SetParameter(6,0);
            CatPotPars[1]->SetParameter(7,1);

            Kitty->SetShortRangePotential(0,0,fDlmPot,*CatPotPars[0]);
            Kitty->SetShortRangePotential(1,0,fDlmPot,*CatPotPars[1]);

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);

        }
        else if(Inter=="Ledni"){

            StandardInter = Inter;

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,4,NumMomBins,MomBins,Lednicky_SingletTriplet);
        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else if(PartPair.IsIt("Lambda","Lambda")){
        if(Inter=="Ledni"){

            StandardInter = Inter;

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,2,NumMomBins,MomBins,Lednicky_Identical_Singlet);
        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else if(PartPair.IsIt("Proton","Xim")){
        if(Inter=="HALQCD_Preliminary"){

            StandardInter = Inter;

            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);

            Kitty->SetExcludeFailedBins(false);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(4);
            Kitty->SetNumPW(0,1);
            Kitty->SetNumPW(1,1);
            Kitty->SetNumPW(2,1);
            Kitty->SetNumPW(3,1);
            Kitty->SetSpin(0,0);//I=0; S=0
            Kitty->SetSpin(1,1);//I=0; S=1
            Kitty->SetSpin(2,0);//I=1; S=0
            Kitty->SetSpin(3,1);//I=1; S=1
            Kitty->SetChannelWeight(0, 1./8.);
            Kitty->SetChannelWeight(1, 3./8.);
            Kitty->SetChannelWeight(2, 1./8.);
            Kitty->SetChannelWeight(3, 3./8.);

            Kitty->SetQ1Q2(-1);
            Kitty->SetPdgId(2212, 3312);
            const double Mass_p=938.272; const double Mass_Xim = 1321.7;
            Kitty->SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

            CatPotPars[0]->SetParameter(0,pXim_Lattice);
            CatPotPars[0]->SetParameter(1,12);
            CatPotPars[0]->SetParameter(2,0);
            CatPotPars[0]->SetParameter(3,-1);
            CatPotPars[0]->SetParameter(4,1);
            CatPotPars[0]->SetParameter(5,0);
            CatPotPars[0]->SetParameter(6,0);
            CatPotPars[0]->SetParameter(7,0);

            CatPotPars[1]->SetParameter(0,pXim_Lattice);
            CatPotPars[1]->SetParameter(1,12);
            CatPotPars[1]->SetParameter(2,0);
            CatPotPars[1]->SetParameter(3,-1);
            CatPotPars[1]->SetParameter(4,1);
            CatPotPars[1]->SetParameter(5,1);
            CatPotPars[1]->SetParameter(6,0);
            CatPotPars[1]->SetParameter(7,1);

            CatPotPars[2]->SetParameter(0,pXim_Lattice);
            CatPotPars[2]->SetParameter(1,6);
            CatPotPars[2]->SetParameter(2,1);
            CatPotPars[2]->SetParameter(3,1);
            CatPotPars[2]->SetParameter(4,1);
            CatPotPars[2]->SetParameter(5,0);
            CatPotPars[2]->SetParameter(6,0);
            CatPotPars[2]->SetParameter(7,0);

            CatPotPars[3]->SetParameter(0,pXim_Lattice);
            CatPotPars[3]->SetParameter(1,6);
            CatPotPars[3]->SetParameter(2,1);
            CatPotPars[3]->SetParameter(3,1);
            CatPotPars[3]->SetParameter(4,1);
            CatPotPars[3]->SetParameter(5,1);
            CatPotPars[3]->SetParameter(6,0);
            CatPotPars[3]->SetParameter(7,1);

            Kitty->SetShortRangePotential(0,0,fDlmPot,*CatPotPars[0]);
            Kitty->SetShortRangePotential(1,0,fDlmPot,*CatPotPars[1]);
            Kitty->SetShortRangePotential(2,0,fDlmPot,*CatPotPars[2]);
            Kitty->SetShortRangePotential(3,0,fDlmPot,*CatPotPars[3]);

            Kitty->SetMaxRad(64);
            Kitty->SetMaxRho(32);

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);
        }
        else if(Inter=="CoulombOnly"){

            StandardInter = Inter;

            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);

            Kitty->SetExcludeFailedBins(false);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(2);
            Kitty->SetNumPW(0,1);
            Kitty->SetNumPW(1,1);
            Kitty->SetSpin(0,0);//S=0
            Kitty->SetSpin(1,1);//S=1
            Kitty->SetChannelWeight(0, 1./4.);
            Kitty->SetChannelWeight(1, 3./4.);

            Kitty->SetQ1Q2(-1);
            Kitty->SetPdgId(2212, 3312);
            const double Mass_p=938.272; const double Mass_Xim = 1321.7;
            Kitty->SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);
        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else if(PartPair.IsIt("Proton","Sigma0")){
        if(Inter=="LedniOliver"){
            StandardInter = Inter;
            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,NumMomBins,MomBins,Lednicky_gauss_Sigma0);
        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else if(PartPair.IsIt("Proton","Xim1530")){
        if(Inter=="CoulombOnly"){
            StandardInter = Inter;
            if(Kitty){delete Kitty;}
            Kitty = new CATS();
            Kitty->SetAnaSource(GaussSource, *CatSourcePars);
            Kitty->SetUseAnalyticSource(true);

            Kitty->SetExcludeFailedBins(false);
            Kitty->SetMomBins(NumMomBins,MomBins);

            Kitty->SetNumChannels(2);
            Kitty->SetNumPW(0,1);
            Kitty->SetNumPW(1,1);
            Kitty->SetSpin(0,0);//S=0
            Kitty->SetSpin(1,1);//S=1
            Kitty->SetChannelWeight(0, 1./4.);
            Kitty->SetChannelWeight(1, 3./4.);

            Kitty->SetQ1Q2(-1);
            //!Maybe wrong?
            Kitty->SetPdgId(2212, 3312);
            const double Mass_p=938.272; const double Mass_Xim1530 = 1535;
            Kitty->SetRedMass( (Mass_p*Mass_Xim1530)/(Mass_p+Mass_Xim1530) );

            Kitty->KillTheCat();

            if(CorrelationModel){delete CorrelationModel;}
            CorrelationModel = new DLM_Ck(1,0,*Kitty);
        }
        else{
            printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
        }
    }
    else{
        printf("WARNING: The interaction %s is unknown!\n",Inter.Data());
    }

    delete [] MomBins;
}

TString FemtoExpPair::GetStandardInteraction(){
    return StandardInter;
}

const FemtoExperiment* FemtoExpPair::GetExperiment() const{
    return &Experiment;
}

const FemtoParticle* FemtoExpPair::GetParticle(const unsigned& WhichOne) const{
    return PartPair.GetParticle(WhichOne);
}
