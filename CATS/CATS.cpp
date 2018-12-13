
#include "CATS.h"

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gsl_sf_coulomb.h"
#include "gsl_sf_bessel.h"
#include "gsl_sf_legendre.h"

#include "DLM_CppTools.h"
#include "CATSconstants.h"

//#include <omp.h>
//#include <unistd.h>

using namespace std;

CATS::CATS():
    NumPotPars(2),NumSourcePars(3)
    {
    IdenticalParticles = false;
    Q1Q2 = 0;
    Gamow = false;
    RedMass = 0;
    pdgID[0] = 0;
    pdgID[1] = 0;
    QuantumStatistics = -1;
    NumCh = 0;
    NumMomBins = 0;
    NumIpBins = 0;
    StartRad = 0.005*FmToNu;
    NumGridPts = 0;
    EpsilonProp = 5e-6;
    EpsilonConv = 5e-6;
    MaxRad = 32.*FmToNu;
    MaxRho = 16;
    MaxPw = 256;
    ExcludeFailedConvergence = true;
    GridMinDepth = 5;
    GridMaxDepth = 0;
    GridEpsilon = 0;
    NumPairs = 0;
    WeightIp = NULL;
    WeightIpError = NULL;
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;

    MaxPairsPerBin = 8e3;
    MaxPairsToRead = 4294967295;
//    MaxPairsToLoad = 4294967295;
    MixingDepth = 1;
    TauCorrection = false;
    UseAnalyticSource = false;
    ThetaDependentSource = false;
    TransportRenorm = 1;
    PoorManRenorm = 1;
    SourceMinRad=-1e9;
    SourceMaxRad=1e9;
    MinTotPairMom = -1;
    MaxTotPairMom = 1e100;
    LoadedMinTotPairMom = -1;
    LoadedMaxTotPairMom = 1e100;
    UseTotMomCut = false;
    Notifications = nAll;

    Spin = NULL;
    NumPW = NULL;
    MomBin = NULL;
    MomBinCenter = NULL;
    IpBin = NULL;
    ChannelWeight = NULL;
    SavedWaveFunBins = NULL;
    PhaseShift = NULL;
    PhaseShiftF = NULL;
    WaveFunRad = NULL;
    WaveFunctionU = NULL;
    MomBinConverged = NULL;
    InputFileName = NULL;
    RelativeMomentum = NULL;
    RelativePosition = NULL;
    RelativeCosTheta = NULL;
    TotalPairMomentum = NULL;
    LoadedPairsPerMomBin = NULL;
    LoadedPairsPerBin = NULL;
    PairMomBin = NULL;
    PairIpBin = NULL;
    GridBoxId = NULL;
    BaseSourceGrid = NULL;
    kSourceGrid = NULL;
    WaveFunction2 = NULL;
    kbSourceGrid = NULL;

    ShortRangePotential = NULL;

    AnalyticSource = NULL;
    ForwardedSource = NULL;

    kCorrFun=NULL;
    kCorrFunErr=NULL;
    kbCorrFun=NULL;
    kbCorrFunErr=NULL;
    PotPar = NULL;
    AnaSourcePar = NULL;

    ExternalWF = NULL;
    NumExtWfRadBins = NULL;
    ExtWfRadBins = NULL;

    RefPartWave = NULL;
    SolvedPartWave = NULL;
    LegPol = NULL;

    CPF = NULL;

    SetIpBins(1, -1000, 1000);
}

CATS::~CATS(){
    DelAll();
    if(Spin) {delete[]Spin; Spin=NULL;}
    if(NumPW) {delete[]NumPW; NumPW=NULL;}
    if(MomBin) {delete[]MomBin; MomBin=NULL;}
    if(MomBinCenter) {delete[]MomBinCenter; MomBinCenter=NULL;}
    if(IpBin) {delete[]IpBin; IpBin=NULL;}
    if(ChannelWeight) {delete[]ChannelWeight; ChannelWeight=NULL;}
    if(InputFileName) {delete[]InputFileName; InputFileName=NULL;}
    if(ShortRangePotential){
        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            if(ShortRangePotential[usCh]){
                delete[]ShortRangePotential[usCh];
                ShortRangePotential[usCh] = NULL;

                delete[]PotPar[usCh];
                PotPar[usCh] = NULL;
            }
        }
        delete[]ShortRangePotential; ShortRangePotential=NULL;
        delete[]PotPar; PotPar=NULL;
    }
    if(BaseSourceGrid){delete BaseSourceGrid; BaseSourceGrid=NULL;}
    if(RefPartWave){delete[]RefPartWave; RefPartWave=NULL;}
    if(SolvedPartWave){delete[]SolvedPartWave; SolvedPartWave=NULL;}
    if(LegPol){delete[]LegPol; LegPol=NULL;}
    if(CPF){delete[]CPF; CPF=NULL;}
}

//N.B. While those guy seemingly do not depend in NumIpBins directly,
//actually the length of the array we reserve for all those guys depends on it, thus
//we will need to reinitialize all of them!
void CATS::DelMomIpMp(){
    if(RelativeMomentum){
        delete [] RelativeMomentum; RelativeMomentum=NULL;
        delete [] RelativePosition; RelativePosition=NULL;
        delete [] RelativeCosTheta; RelativeCosTheta=NULL;
        delete [] PairMomBin; PairMomBin=NULL;
        delete [] PairIpBin; PairIpBin=NULL;
        delete [] GridBoxId; GridBoxId=NULL;
    }

    if(TotalPairMomentum){
        delete [] TotalPairMomentum; TotalPairMomentum=NULL;
        UseTotMomCut = false;
    }
}

void CATS::DelIp(){
    if(!WeightIp) return;
    delete [] WeightIp; WeightIp=NULL;
    delete [] WeightIpError; WeightIpError=NULL;
}

void CATS::DelMomChPw(){
    if(SavedWaveFunBins){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned short usCh=0; usCh<NumCh; usCh++){
                for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
                    if(WaveFunRad[uMomBin][usCh][usPW])
                        delete [] WaveFunRad[uMomBin][usCh][usPW];
                    if(WaveFunctionU[uMomBin][usCh][usPW])
                        delete [] WaveFunctionU[uMomBin][usCh][usPW];
                }
                delete [] SavedWaveFunBins[uMomBin][usCh];
                delete [] PhaseShift[uMomBin][usCh];
                delete [] WaveFunRad[uMomBin][usCh];
                delete [] WaveFunctionU[uMomBin][usCh];
            }
            delete [] SavedWaveFunBins[uMomBin];
            delete [] PhaseShift[uMomBin];
            delete [] WaveFunRad[uMomBin];
            delete [] WaveFunctionU[uMomBin];
        }

        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
                delete [] PhaseShiftF[usCh][usPW];
            }
            delete [] PhaseShiftF[usCh];
        }

        delete [] SavedWaveFunBins; SavedWaveFunBins=NULL;
        delete [] PhaseShift; PhaseShift=NULL;
        delete [] PhaseShiftF; PhaseShiftF=NULL;
        delete [] WaveFunRad; WaveFunRad=NULL;
        delete [] WaveFunctionU; WaveFunctionU=NULL;
    }
    if(ExternalWF){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned short usCh=0; usCh<NumCh; usCh++){
                delete [] ExternalWF[uMomBin][usCh];
                delete [] NumExtWfRadBins[uMomBin][usCh];
                delete [] ExtWfRadBins[uMomBin][usCh];
            }
            delete [] ExternalWF[uMomBin];
            delete [] NumExtWfRadBins[uMomBin];
            delete [] ExtWfRadBins[uMomBin];
        }

        delete [] ExternalWF; ExternalWF=NULL;
        delete [] NumExtWfRadBins; NumExtWfRadBins=NULL;
        delete [] ExtWfRadBins; ExtWfRadBins=NULL;
    }
}

void CATS::DelMomCh(){
    if(WaveFunction2){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                delete [] WaveFunction2[uMomBin][uGrid];
            }
            delete [] WaveFunction2[uMomBin];
        }
        delete [] WaveFunction2; WaveFunction2 = NULL;
    }
}

void CATS::DelMom(){
    if(MomBinConverged) {delete [] MomBinConverged; MomBinConverged=NULL;}
    if(LoadedPairsPerMomBin) {delete [] LoadedPairsPerMomBin; LoadedPairsPerMomBin = NULL;}
    if(kCorrFun){
        delete [] kCorrFun; kCorrFun=NULL;
        delete [] kCorrFunErr; kCorrFunErr=NULL;
    }
    if(kSourceGrid){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            if(kSourceGrid[uMomBin]) delete kSourceGrid[uMomBin];
        }
        delete [] kSourceGrid; kSourceGrid=NULL;
    }
}

void CATS::DelMomIp(){
    if(kbCorrFun){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete [] kbCorrFun[uMomBin];
            delete [] kbCorrFunErr[uMomBin];
        }
        delete [] kbCorrFun; kbCorrFun=NULL;
        delete [] kbCorrFunErr; kbCorrFunErr=NULL;
    }
    if(LoadedPairsPerBin){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete [] LoadedPairsPerBin[uMomBin];
        }
        delete [] LoadedPairsPerBin; LoadedPairsPerBin=NULL;
    }
    if(kbSourceGrid){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                if(kbSourceGrid[uMomBin][uIpBin]) delete kbSourceGrid[uMomBin][uIpBin];
            }
            delete [] kbSourceGrid[uMomBin];
        }
        delete [] kbSourceGrid; kbSourceGrid=NULL;
    }
}

void CATS::DelAllMom(){
    DelMomIpMp();
    DelMomChPw();
    DelMomCh();
    DelMom();
    DelMomIp();
}

void CATS::DelAllIp(){
    DelMomIpMp();
    DelIp();
    DelMomIp();
}

void CATS::DelAll(){
    DelAllMom();
    DelAllIp();
}

void CATS::SetRedMass(const double& redMass){
    if(redMass==RedMass) return;
    RedMass = redMass;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;
}
double CATS::GetRedMass() const{
    return RedMass;
}

void CATS::SetPdgId(const int& id1, const int& id2){
    if(pdgID[0]==id1 && pdgID[1]==id2) return;
    pdgID[0] = id1;
    pdgID[1] = id2;
    //if the QS is set manually, then the IdenticalParticles does not depend on the pdgID
    if(QuantumStatistics==0||QuantumStatistics==1) return;
    if(id1==id2) IdenticalParticles=true;
    else IdenticalParticles=false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
void CATS::GetPdgId(int& id1, int& id2) const{
    id1 = pdgID[0];
    id2 = pdgID[1];
}

void CATS::SetQuantumStatistics(short qs){
    if(qs!=0 && qs!=1) qs=-1;
    if(qs==QuantumStatistics) return;
    QuantumStatistics = qs;
    bool Identical;
    switch(QuantumStatistics){
        case 0 : Identical=false; break;
        case 1 : Identical=true; break;
        default : Identical=(pdgID[0]==pdgID[1]); break;
    }
    if(Identical==IdenticalParticles) return;
    IdenticalParticles = Identical;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
short CATS::GetQuantumStatistics() const{
    return QuantumStatistics;
}

//If the number of channels is changed, all previous input about the
//channels themselves is lost (i.e. NumPW, WhichPartialWave and the potentials are reset!)
void CATS::SetNumChannels(const unsigned short& numCh){
    if(NumCh == numCh) return;
    if(!numCh){
        if(Notifications>=nError){
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetNumChannels(unsigned short numCh)\n");
            printf("         NumCh cannot be zero!\n");
        }
        return;
    }

    if(Spin) {delete[]Spin; Spin=NULL;}
    Spin = new unsigned short [numCh];

    if(NumPW) {delete[]NumPW; NumPW=NULL;}
    NumPW = new unsigned short [numCh];

    if(ShortRangePotential){
        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            if(ShortRangePotential[usCh]){
                delete [] ShortRangePotential[usCh];
                delete [] PotPar[usCh];
            }
        }
        delete[]ShortRangePotential; ShortRangePotential=NULL;
        delete[]PotPar; PotPar=NULL;
    }
    ShortRangePotential = new CatsPotential* [numCh];

    PotPar = new double** [numCh];
    for(unsigned short usCh=0; usCh<numCh; usCh++){
        ShortRangePotential[usCh] = NULL;
        PotPar[usCh] = NULL;
    }

    if(ChannelWeight) {delete[]ChannelWeight; ChannelWeight=NULL;}
    ChannelWeight = new double [numCh];
    for(unsigned short usCh=0; usCh<numCh; usCh++){
        Spin[usCh] = 0;
        NumPW[usCh] = 0;
        ChannelWeight[usCh] = 0;
    }

    DelMomChPw();
    NumCh = numCh;

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
unsigned short CATS::GetNumChannels() const{
    return NumCh;
}

void CATS::SetNumPW(const unsigned short& usCh, const unsigned short& numPW){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetNumPW(unsigned short usCh, unsigned short numPW)\n");
        return;
    }
    if(NumPW[usCh]==numPW) return;
    DelMomChPw();
    NumPW[usCh] = numPW;

    if(ShortRangePotential[usCh]) delete[]ShortRangePotential[usCh];
    ShortRangePotential[usCh] = new CatsPotential [numPW];
    for(unsigned short usPW=0; usPW<numPW; usPW++){
        ShortRangePotential[usCh][usPW] = 0;
    }

    if(PotPar[usCh]) delete[]PotPar[usCh];
    PotPar[usCh] = new double* [numPW];

    //Reserve memory for the output
    SavedWaveFunBins = new unsigned** [NumMomBins];
    PhaseShift = new double** [NumMomBins];
    WaveFunRad = new double*** [NumMomBins];
    WaveFunctionU = new complex<double>*** [NumMomBins];
    ExternalWF = new const complex<double>*** [NumMomBins];
    NumExtWfRadBins = new unsigned** [NumMomBins];
    ExtWfRadBins = new const double*** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        SavedWaveFunBins[uMomBin] = new unsigned* [NumCh];
        PhaseShift[uMomBin] = new double* [NumCh];
        WaveFunRad[uMomBin] = new double** [NumCh];
        WaveFunctionU[uMomBin] = new complex<double>** [NumCh];
        ExternalWF[uMomBin] = new const complex<double>** [NumCh];
        NumExtWfRadBins[uMomBin] = new unsigned* [NumCh];
        ExtWfRadBins[uMomBin] = new const double** [NumCh];
        for(unsigned short usCh=0; usCh<NumCh; usCh++){
            SavedWaveFunBins[uMomBin][usCh] = new unsigned [NumPW[usCh]];
            PhaseShift[uMomBin][usCh] = new double [NumPW[usCh]];
            WaveFunRad[uMomBin][usCh] = new double* [NumPW[usCh]];
            WaveFunctionU[uMomBin][usCh] = new complex<double>* [NumPW[usCh]];
            ExternalWF[uMomBin][usCh] = new const complex<double>* [NumPW[usCh]];
            NumExtWfRadBins[uMomBin][usCh] = new unsigned [NumPW[usCh]];
            ExtWfRadBins[uMomBin][usCh] = new const double* [NumPW[usCh]];
            for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
                SavedWaveFunBins[uMomBin][usCh][usPW] = 0;
                PhaseShift[uMomBin][usCh][usPW] = 0;
                WaveFunRad[uMomBin][usCh][usPW] = NULL;
                WaveFunctionU[uMomBin][usCh][usPW] = NULL;
                ExternalWF[uMomBin][usCh][usPW] = NULL;
                NumExtWfRadBins[uMomBin][usCh][usPW] = 0;
                ExtWfRadBins[uMomBin][usCh][usPW] = NULL;
            }
        }
    }
    PhaseShiftF = new float** [NumCh];
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        PhaseShiftF[usCh] = new float* [NumPW[usCh]];
        for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
            PhaseShiftF[usCh][usPW] = new float [NumMomBins];
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                PhaseShiftF[usCh][usPW][uMomBin] = 0;
            }
        }
    }

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
unsigned short CATS::GetNumPW(const unsigned short& usCh) const{
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::GetNumPW(unsigned short usCh)\n");
        return 0;
    }
    return NumPW[usCh];
}

void CATS::SetSpin(const unsigned short& usCh, const unsigned short& spin){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetSpin(...)\n");
        return;
    }
    if(spin==Spin[usCh]) return;
    Spin[usCh] = spin;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
unsigned short CATS::GetSpin(const unsigned short& usCh) const{
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::GetSpin(...)\n");
        return 0;
    }
    return Spin[usCh];
}

void CATS::SetQ1Q2(const int& q1q2){
    if(Q1Q2==q1q2) return;
    Q1Q2 = q1q2;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;
}
int CATS::GetQ1Q2() const{
    return Q1Q2;
}
void CATS::SetGamow(const bool& gamow){
    if(Gamow==gamow) return;
    Gamow = gamow;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;
}
bool CATS::GetGamow() const{
    return Gamow;
}

unsigned CATS::GetNumMomBins() const{
    return NumMomBins;
}

unsigned CATS::GetNumIpBins() const{
    return NumIpBins;
}

unsigned CATS::GetNumPairs() const{
    return NumPairs;
}

void CATS::SetMomBins(const unsigned& nummombins, const double* mombins, const double* bincenter){
    if(!nummombins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
        return;
    }
    if(!mombins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
        return;
    }
    //check if the momentum bins set are the same as before. If yes, change nothing
    if(nummombins==NumMomBins){
        bool SameBinning = true;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            SameBinning *= (mombins[uMomBin]==MomBin[uMomBin]);
        }
        if(SameBinning) return;
    }
    if(nummombins!=NumMomBins || !MomBin){
        if(MomBin) {delete[]MomBin; MomBin=NULL;}
        if(MomBinCenter) {delete[]MomBinCenter; MomBinCenter=NULL;}
        MomBin = new double [nummombins+1];
        MomBinCenter = new double [nummombins];
        DelAllMom();
        NumMomBins = nummombins;
    }

    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;

    for(unsigned uBin=0; uBin<=NumMomBins; uBin++){
        MomBin[uBin] = mombins[uBin];
        if(uBin){
            if(bincenter) MomBinCenter[uBin-1]=bincenter[uBin-1];
            else MomBinCenter[uBin-1] = 0.5*(mombins[uBin-1]+mombins[uBin]);
        }

        if(MomBin[uBin]<0){
            if(Notifications>=nError){
                printf("\033[1;31mERROR:\033[0m CATS::SetMomBins(const unsigned& nummombins, const double* mombins)\n");
                printf("         The momentum should be positive!\n");
            }
            return;
        }
    }
}
void CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom){
    if(!nummombins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    if(MinMom>MaxMom){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    if(MinMom==MaxMom && nummombins!=1){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetMomBins(const unsigned& nummombins, const double& MinMom, const double& MaxMom)\n");
        return;
    }
    //check if the momentum bins set are the same as before. If yes, change nothing
    double BinWidth = (MaxMom-MinMom)/double(nummombins);
    if(nummombins==NumMomBins && MomBin[0]==MinMom && MomBin[nummombins]==MinMom+double(nummombins)*BinWidth) return;
    if(nummombins!=NumMomBins || !MomBin){
        if(MomBin) {delete[]MomBin; MomBin=NULL;}
        if(MomBinCenter) {delete[]MomBinCenter; MomBinCenter=NULL;}
        MomBin = new double [nummombins+1];
        MomBinCenter = new double [nummombins];
        DelAllMom();
        NumMomBins = nummombins;
    }
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
    GamowCorrected = false;

    for(unsigned uBin=0; uBin<=NumMomBins; uBin++){
        MomBin[uBin] = MinMom+double(uBin)*BinWidth;
        if(uBin!=NumMomBins) MomBinCenter[uBin] = MinMom+double(uBin)*BinWidth+0.5*BinWidth;
    }
}

void CATS::SetIpBins(const unsigned& numBbins, const double* imppar){
    if(!numBbins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetIpBins(const unsigned& numBbins, const double* imppar)\n");
        return;
    }
    if(!imppar){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetIpBins(const unsigned& numBbins, const double* imppar)\n");
        return;
    }
    if(numBbins!=NumIpBins || !IpBin){
        if(IpBin) {delete[]IpBin; IpBin=NULL;}
        IpBin = new double [numBbins+1];
        DelAllIp();
        NumIpBins = numBbins;
    }
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    for(unsigned uBin=0; uBin<=NumIpBins; uBin++){
        IpBin[uBin] = imppar[uBin];
    }
}
void CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar){
    if(!numBbins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(MinImpPar>MaxImpPar){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(MinImpPar==MaxImpPar && numBbins!=1){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetIpBins(const unsigned& numBbins, const double& MinImpPar, const double& MaxImpPar)\n");
        return;
    }
    if(numBbins!=NumIpBins || !IpBin){
        if(IpBin) {delete[]IpBin; IpBin=NULL;}
        IpBin = new double [numBbins+1];
        DelAllIp();
        NumIpBins = numBbins;
    }
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    double BinWidth = (MaxImpPar-MinImpPar)/double(NumIpBins);
    for(unsigned uBin=0; uBin<=NumIpBins; uBin++){
        IpBin[uBin] = MinImpPar+double(uBin)*BinWidth;
    }
}

void CATS::SetChannelWeight(const unsigned short& usCh, const double& weight){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetSpinWeight(const unsigned short& usCh, const double& weight)\n");
        return;
    }
    if(ChannelWeight[usCh]==weight) return;
    ChannelWeight[usCh] = weight;
    ComputedCorrFunction = false;
}

double CATS::GetChannelWeight(const unsigned short& usCh) const{
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::GetSpinWeight(const unsigned short& usCh)\n");
        return 0;
    }
    return ChannelWeight[usCh];
}

void CATS::SetStartRad(const double& srad){
    if(StartRad==fabs(srad)) return;
    StartRad = fabs(srad);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetStartRad() const{
    return StartRad;
}

void CATS::SetEpsilonProp(const double& epsp){
    if(EpsilonProp==fabs(epsp)) return;
    //make sure that EpsilonProp is always non-zero and positive
    EpsilonProp = fabs(epsp)+1e-64;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetEpsilonProp() const{
    return EpsilonProp;
}

void CATS::SetEpsilonConv(const double& epsc){
    if(EpsilonConv==fabs(epsc)) return;
    EpsilonConv = fabs(epsc);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
double CATS::GetEpsilonConv() const{
    return EpsilonConv;
}

void CATS::SetMaxRad(const double& maxrad){
    if(MaxRad==fabs(maxrad*FmToNu)) return;
    MaxRad = fabs(maxrad*FmToNu);
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

double CATS::GetMaxRad() const{
    return MaxRad;
}

void CATS::SetMaxRho(const double& maxrho){
    if(MaxRho==fabs(maxrho)) return;
    MaxRho = fabs(maxrho);
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

double CATS::GetMaxRho() const{
    return MaxRho;
}
void CATS::SetMaxPw(const unsigned short& maxpw){
    if(MaxPw==maxpw) return;
    MaxPw = maxpw;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

unsigned short CATS::GetMaxPw() const{
    return MaxRho;
}

void CATS::SetExcludeFailedBins(const bool& efb){
    if(ExcludeFailedConvergence==efb) return;
    ExcludeFailedConvergence = efb;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
bool CATS::GetExcludeFailedBins() const{
    return ExcludeFailedConvergence;
}

void CATS::SetGridMinDepth(const short& val){
    if(val<0 || val>6) GridMinDepth=0;
    else GridMinDepth=val;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
short CATS::GetGridMinDepth() const{
    return GridMinDepth;
}
void CATS::SetGridMaxDepth(const short& val){
    if(val<5 || val>20) GridMaxDepth=0;
    else GridMaxDepth=val;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}
short CATS::GetGridManDepth() const{
    return GridMaxDepth;
}
void CATS::SetGridEpsilon(const double& val){
    if(val<0 || val>0.125) GridEpsilon=0;
    else GridEpsilon = val;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

}
double CATS::GetGridEpsilon() const{
    return GridEpsilon;
}

void CATS::SetMaxPairsPerBin(unsigned mpp){
    if(!mpp){
        if(Notifications>=nWarning)
            printf("\033[1;33mWARNING:\033[0m MaxPairsPerBin cannot be zero, setting MaxPairsPerBin=1\n");
        mpp=1;
    }
    if(MaxPairsPerBin==mpp) return;
    DelMomIpMp();
    MaxPairsPerBin = mpp;
    LoadedData = false;
    if(!UseAnalyticSource) SourceGridReady = false;
    if(!UseAnalyticSource) SourceUpdated = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned CATS::GetMaxPairsPerBin() const{
    return MaxPairsPerBin;
}

void CATS::SetMaxPairsToRead(unsigned mpp){
    if(!mpp){
        if(Notifications>=nWarning)
            printf("\033[1;33mWARNING:\033[0m MaxPairsToRead cannot be zero, setting MaxPairsToRead=1\n");
        mpp=1;
    }
    if(MaxPairsToRead==mpp) return;
    MaxPairsToRead = mpp;
    LoadedData = false;
    if(!UseAnalyticSource) SourceGridReady = false;
    if(!UseAnalyticSource) SourceUpdated = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned CATS::GetMaxPairsToRead() const{
    return MaxPairsToRead;
}

//void CATS::SetMaxPairsToLoad(unsigned mpp){
//    if(!mpp){
//        printf("WARNING: MaxPairsToLoad cannot be zero, setting MaxPairsToLoad=1\n");
//        mpp=1;
//    }
//    if(MaxPairsToLoad==mpp) return;
//    MaxPairsToLoad = mpp;
//    LoadedData = false;
//    if(!UseAnalyticSource) ComputedCorrFunction = false;
//}
//unsigned CATS::GetMaxPairsToLoad(){
//    return MaxPairsToLoad;
//}
void CATS::SetMixingDepth(const unsigned short& mix){
    if(MixingDepth==mix) return;
    MixingDepth = mix?mix:1;
    LoadedData = false;
    if(!UseAnalyticSource) SourceGridReady = false;
    if(!UseAnalyticSource) SourceUpdated = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
unsigned short CATS::GetMixingDepth() const{
    return MixingDepth;
}
void CATS::SetTauCorrection(const bool& tc){
    if(TauCorrection==tc) return;
    TauCorrection = tc;
    LoadedData = false;
    if(!UseAnalyticSource) SourceGridReady = false;
    if(!UseAnalyticSource) SourceUpdated = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}
bool CATS::GetTauCorrection() const{
    return TauCorrection;
}

void CATS::SetUseAnalyticSource(const bool& val){
    if(UseAnalyticSource==val) return;
    UseAnalyticSource = val;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
bool CATS::GetUseAnalyticSource() const{
    return UseAnalyticSource;
}

void CATS::SetThetaDependentSource(const bool& val){
    if(ThetaDependentSource==val) return;
    ThetaDependentSource = val;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
    if(ThetaDependentSource && !RefPartWave){
        RefPartWave = new complex<double> [MaxPw];
        SolvedPartWave = new complex<double> [MaxPw];
        LegPol = new double [MaxPw];
    }
    if(!ThetaDependentSource && RefPartWave){
        delete [] RefPartWave; RefPartWave=NULL;
        delete [] SolvedPartWave; SolvedPartWave=NULL;
        delete [] LegPol; LegPol=NULL;
    }
}
bool CATS::GetThetaDependentSource() const{
    return ThetaDependentSource;
}

void CATS::SetTransportRenorm(const double& val){
    if(TransportRenorm==fabs(val)) return;
    TransportRenorm = fabs(val);
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
double CATS::GetTransportRenorm() const{
    return TransportRenorm;
}

void CATS::SetPoorManRenorm(const double& val){
    if(PoorManRenorm==fabs(val)) return;
    PoorManRenorm = fabs(val);
    ComputedCorrFunction = false;
}
double CATS::GetPoorManRenorm() const{
    return PoorManRenorm;
}

void CATS::SetSourceMinRange(const double& val){
    if(SourceMinRad==val) return;
    if(val>SourceMaxRad) ComputedCorrFunction = false;
    SourceMinRad = val;
}
void CATS::SetSourceMaxRange(const double& val){
    if(SourceMaxRad==val) return;
    if(val<SourceMaxRad) ComputedCorrFunction = false;
    SourceMaxRad = val;
}
double CATS::GetSourceMinRange() const{
    return SourceMinRad;
}
double CATS::GetSourceMaxRange() const{
    return SourceMaxRad;
}

void CATS::SetTotPairMomCut(const double& minval, const double& maxval){
    if(minval<0 || maxval<minval){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in void CATS::SetTotPairMomCut(const double& minval, const double& maxval)");
        return;
    }
    if(minval==MinTotPairMom && maxval==MaxTotPairMom){
        return;
    }

    MinTotPairMom = minval;
    MaxTotPairMom = maxval;
    UseTotMomCut = true;

    if(!TotalPairMomentum) LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
void CATS::GetTotPairMomCut(double& minval, double& maxval) const{
    minval = MinTotPairMom;
    maxval = MaxTotPairMom;
}
void CATS::RemoveTotPairMomCut(){
    UseTotMomCut = false;
}

void CATS::SetNotifications(const short& notify){
    Notifications = notify;
}
short CATS::GetNotifications() const{
    return Notifications;
}

void CATS::SetInputFileName(const char* fname){
    unsigned StrLen = strlen(fname);
    if(!StrLen){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m The input file name is empty!\n");
        return;
    }
    if(InputFileName){
        //if this file was already loaded before
        if(strcmp(InputFileName, fname)==0 && LoadedData){
            return;
        }
        delete [] InputFileName;
        InputFileName = NULL;
    }

    InputFileName = new char [StrLen+1];
    strcpy(InputFileName, fname);

    LoadedData = false;
    if(!UseAnalyticSource) SourceGridReady = false;
    if(!UseAnalyticSource) SourceUpdated = false;
    if(!UseAnalyticSource) ComputedCorrFunction = false;
}

void CATS::GetInputFileName(char* fname) const{
    if(!InputFileName){
        strcpy(fname, "");
        return;
    }
    strcpy(fname, InputFileName);
}

unsigned CATS::GetNumPairsPerBin(const unsigned& uMomBin, const unsigned& uIpBin) const{
    if(uMomBin>=NumMomBins || uIpBin>=NumIpBins || !LoadedData) return 0;
    return LoadedPairsPerBin[uMomBin][uIpBin];
}

unsigned CATS::GetNumPairsPerBin(const unsigned& uMomBin) const{
    if(uMomBin>=NumMomBins || !LoadedData) return 0;
    return LoadedPairsPerMomBin[uMomBin];
}

void CATS::GetPairInfo(const unsigned& uWhichPair,
                     double& RelMom, double& RelPos, double& RelCosTh, double& TotMom) const{

    RelMom=0;
    RelPos=0;
    RelCosTh=0;
    TotMom=0;
    if(!LoadedData) return;
    if(uWhichPair>=NumPairs) return;
    RelMom=RelativeMomentum[uWhichPair];
    RelPos=RelativePosition[uWhichPair]*NuToFm;
    RelCosTh=RelativeCosTheta[uWhichPair];
    TotMom=UseTotMomCut?TotalPairMomentum[uWhichPair]:0;

}
void CATS::GetPairInfo(const unsigned& uWhichPair, double* Output) const{
    GetPairInfo(uWhichPair,
                Output[0],Output[1],Output[2],Output[3]);
}

unsigned CATS::GetLoadedPairs(const unsigned& WhichMomBin, const unsigned& WhichIpBin) const{
    if(!LoadedData) return 0;
    if(WhichMomBin>=NumMomBins) return 0;
    if(WhichIpBin>=NumIpBins) return 0;
    return LoadedPairsPerBin[WhichMomBin][WhichIpBin];
}
unsigned CATS::GetRelativeMomentum(const unsigned& WhichParticle) const{
    if(!LoadedData) return 0;
    if(WhichParticle>=NumPairs) return 0;
    return RelativeMomentum[WhichParticle];
}
unsigned CATS::GetRelativePosition(const unsigned& WhichParticle) const{
    if(!LoadedData) return 0;
    if(WhichParticle>=NumPairs) return 0;
    return RelativePosition[WhichParticle]*NuToFm;
}
unsigned CATS::GetRelativeCosTheta(const unsigned& WhichParticle) const{
    if(!LoadedData) return 0;
    if(WhichParticle>=NumPairs) return 0;
    return RelativeCosTheta[WhichParticle];
}
unsigned CATS::GetTotalPairMomentum(const unsigned& WhichParticle) const{
    if(!LoadedData) return 0;
    if(WhichParticle>=NumPairs) return 0;
    return TotalPairMomentum[WhichParticle];
}

double CATS::GetCorrFun(const unsigned& WhichMomBin) const{
    if(WhichMomBin>=NumMomBins || !kCorrFun) return 0;
    return kCorrFun[WhichMomBin];
}
double CATS::GetCorrFun(const unsigned& WhichMomBin, double& Momentum) const{
    if(WhichMomBin>=NumMomBins || !kCorrFun) return 0;
    Momentum = WhichMomBin<NumMomBins?MomBinCenter[WhichMomBin]:0;
    //Momentum = WhichMomBin<NumMomBins?GetMomentum(WhichMomBin):0;
    return GetCorrFun(WhichMomBin);
}

//in short: here we want to make linear interpolation. For that we need to find two points.
//In perfect case Momentum should be between the two points. However, in case we are working in the first
//half of the 0th bin, or the last half of the last mom. bin, than we need to take either the two point above or
//the two points below the value of Momentum
double CATS::EvalCorrFun(const double& Momentum) const{
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins]) return 0;
    return EvalBinnedFun(Momentum, NumMomBins, MomBin, MomBinCenter, kCorrFun);
    //return EvalBinnedFun(Momentum, NumMomBins, MomBin, NULL, kCorrFun);
}

double CATS::EvalCorrFunErr(const double& Momentum) const{
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins]) return 0;
    return EvalBinnedFun(Momentum, NumMomBins, MomBin, MomBinCenter, kCorrFunErr);
    //return EvalBinnedFun(Momentum, NumMomBins, MomBin, NULL, kCorrFunErr);
}

double CATS::GetCorrFunErr(const unsigned& WhichMomBin) const{
    if(WhichMomBin>=NumMomBins || !kCorrFunErr) return 0;
    return kCorrFunErr[WhichMomBin];
}

double CATS::GetCorrFunErr(const unsigned& WhichMomBin, double& Momentum) const{
    if(WhichMomBin>=NumMomBins || !kCorrFunErr) return 0;
    Momentum = GetMomentum(WhichMomBin);
    return kCorrFunErr[WhichMomBin];
}

double CATS::GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin) const{
    if(WhichMomBin>=NumMomBins || WhichIpBin>=NumIpBins || !kbCorrFun) return 0;
    return kbCorrFun[WhichMomBin][WhichIpBin];
}

double CATS::GetCorrFunIp(const unsigned& WhichMomBin, const unsigned& WhichIpBin, double& Momentum, double& ImpPar) const{
    if(NumMomBins<=WhichMomBin || NumIpBins<=WhichIpBin || !kCorrFun) return 0;
    Momentum = WhichMomBin<NumMomBins?MomBinCenter[WhichMomBin]:0;
    //Momentum = WhichMomBin<NumMomBins?GetMomentum(WhichMomBin):0;
    ImpPar = WhichIpBin<NumIpBins?(IpBin[WhichIpBin]+IpBin[WhichIpBin+1])*0.5:0;
    return GetCorrFunIp(WhichMomBin, WhichIpBin);
}

double CATS::GetPhaseShift(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW) const{
    if(NumMomBins<=WhichMomBin || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    return PhaseShift[WhichMomBin][usCh][usPW];
}

float CATS::EvalPhaseShift(const double& Momentum, const unsigned short& usCh, const unsigned short& usPW) const{
    if(Momentum<MomBin[0] || Momentum>MomBin[NumMomBins] || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    return EvalBinnedFun(Momentum, NumMomBins, MomBin, MomBinCenter, PhaseShiftF[usCh][usPW]);
    //return EvalBinnedFun(Momentum, NumMomBins, MomBin, NULL, PhaseShiftF[usCh][usPW]);
}

unsigned CATS::GetNumRadialWFpts(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW) const{
    return SavedWaveFunBins[WhichMomBin][usCh][usPW];
}

complex<double> CATS::GetRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const unsigned& WhichRadBin) const{
    if(NumMomBins<=WhichMomBin || NumCh<=usCh || NumPW[usCh]<=usPW || SavedWaveFunBins[WhichMomBin][usCh][usPW]<=WhichRadBin) return 0;
    return WaveFunctionU[WhichMomBin][usCh][usPW][WhichRadBin];
}

complex<double> CATS::EvalRadialWaveFunction(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const double& Radius,
                                    const bool& DivideByR) const{
    if(NumMomBins<=WhichMomBin || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    return EvalWaveFunctionU(WhichMomBin, Radius*FmToNu, usCh, usPW, DivideByR);
}

double CATS::EvalWaveFun2(const unsigned& uMomBin, const double& Radius, const double& CosTheta, const unsigned short& usCh){
    if(NumMomBins<=uMomBin || NumCh<=usCh) return 0;
    return EffectiveFunctionTheta(uMomBin,Radius*FmToNu,CosTheta,usCh);
}
double CATS::EvalWaveFun2(const unsigned& uMomBin, const double& Radius, const unsigned short& usCh){
    if(NumMomBins<=uMomBin || NumCh<=usCh) return 0;
    return EffectiveFunction(uMomBin,Radius*FmToNu,usCh);
}

complex<double> CATS::EvalAsymptoticRadialWF(const unsigned& WhichMomBin, const unsigned short& usCh, const unsigned short& usPW, const double& Radius,
                                    const bool& DivideByR){
    if(NumMomBins<=WhichMomBin || NumCh<=usCh || NumPW[usCh]<=usPW) return 0;
    return EvalWaveFunctionU(WhichMomBin, Radius*FmToNu, usCh, usPW, DivideByR, true);
}

complex<double> CATS::EvalReferenceRadialWF(const unsigned& WhichMomBin, const unsigned short& usPW, const double& Radius,
                                    const bool& DivideByR){
    if(NumMomBins<=WhichMomBin) return 0;
    double MultFactor = DivideByR?1./(Radius*FmToNu+1e-64):1;
    return ReferencePartialWave(Radius*FmToNu, GetMomentum(WhichMomBin), usPW, Q1Q2)*MultFactor;
}

double CATS::GetMomentum(const unsigned& WhichMomBin) const{
    if(NumMomBins<=WhichMomBin) return 0;
    return MomBinCenter[WhichMomBin];
    //return 0.5*(MomBin[WhichMomBin]+MomBin[WhichMomBin+1]);
}

double CATS::GetMomBinLowEdge(const unsigned& WhichMomBin) const{
    if(NumMomBins<WhichMomBin) return 0;
    return MomBin[WhichMomBin];
}

double CATS::GetMomBinUpEdge(const unsigned& WhichMomBin) const{
    if(NumMomBins<=WhichMomBin) return 0;
    return MomBin[WhichMomBin+1];
}

double* CATS::CopyMomBin(){
    double* MomBinCopy = new double [NumMomBins+1];
    for(unsigned uMomBin=0; uMomBin<=NumMomBins; uMomBin++){
        MomBinCopy[uMomBin] = MomBin[uMomBin];
    }
    return MomBinCopy;
}

const double& CATS::FmNu() const{
    return FmToNu;
}

const double& CATS::NuFm() const{
    return NuToFm;
}

void CATS::RemoveShortRangePotential(){
    if(!ShortRangePotential) return;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        for(unsigned short usPW=0; usPW<NumPW[usCh]; usPW++){
            ShortRangePotential[usCh][usPW] = 0;
            PotPar[usCh][usPW] = 0;
        }
    }
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::RemoveShortRangePotential(const unsigned& usCh, const unsigned& usPW){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::RemoveShortRangePotential(...)\n");
        return;
    }
    if(usPW>=NumPW[usCh]){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::RemoveShortRangePotential(...)\n");
        return;
    }
    if(!ShortRangePotential) return;
    if(!ShortRangePotential[usCh]) return;
    ShortRangePotential[usCh][usPW] = 0;
    PotPar[usCh][usPW] = 0;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::SetShortRangePotential(const unsigned& usCh, const unsigned& usPW, double (*pot)(double* Pars), double* Pars){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }
    if(NumPW[usCh] && usPW>=NumPW[usCh]){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }
    if(ShortRangePotential[usCh][usPW]==pot && PotPar[usCh][usPW]==Pars){
        return;
    }

    ShortRangePotential[usCh][usPW] = pot;
    PotPar[usCh][usPW] = Pars;

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;
}

void CATS::SetShortRangePotential(const unsigned& usCh, const unsigned& usPW, const unsigned& WhichPar, const double& Value){
    if(usCh>=NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }
    if(usPW>=NumPW[usCh]){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::SetShortRangePotential(...)\n");
        return;
    }

    if(PotPar[usCh][usPW][NumPotPars+WhichPar]==Value) return;
    if(!ShortRangePotential[usCh][usPW]) return;

    PotPar[usCh][usPW][NumPotPars+WhichPar] = Value;
    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

}

void CATS::RemoveAnaSource(){
    if(!AnaSourcePar) return;
    AnalyticSource = NULL;
    ForwardedSource = NULL;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
void CATS::SetAnaSource(double (*AS)(double*), double* Pars){
    if(AnalyticSource==AS && AnaSourcePar==Pars) return;
    if(!Pars){
        if(Notifications>=nWarning) printf("\033[1;33mWARNING:\033[0m NULL pointer to the source parameters!\n");
        return;
    }
    AnalyticSource = AS;
    ForwardedSource = NULL;
    AnaSourcePar = Pars;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
void CATS::SetAnaSource(double (*FS)(void*, double*), void* context){
    if(ForwardedSource==FS) return;
    if(!FS){
        if(Notifications>=nWarning) printf("\033[1;33mWARNING:\033[0m NULL pointer to the source function!\n");
        return;
    }
    if(!context){
        if(Notifications>=nWarning) printf("\033[1;33mWARNING:\033[0m NULL pointer to the source context!\n");
        return;
    }
    AnalyticSource = NULL;
    ForwardedSource = FS;
    AnaSourcePar = ForwardedSourcePar;
    SourceContext = context;
    SourceGridReady = false;
    SourceUpdated = false;
    ComputedCorrFunction = false;
}
void CATS::SetAnaSource(const unsigned& WhichPar, const double& Value, const bool& SmallChange){
    if(!AnaSourcePar) return;
    //if we use a member function, we assume that there are no parameters to be changed here
    if(ForwardedSource){
        if(Notifications>=nWarning) printf("\033[1;33mWARNING:\033[0m Using a source member function does not allow to set any parameters!\n");
        return;
    }
    if(AnaSourcePar[NumSourcePars+WhichPar]==Value) return;
    AnaSourcePar[NumSourcePars+WhichPar] = Value;
    if(UseAnalyticSource){
        SourceGridReady = SourceGridReady?SmallChange:false;
        SourceUpdated = false;
        ComputedCorrFunction = false;
    }
}
double CATS::GetAnaSourcePar(const unsigned& WhichPar) const{
    if(!UseAnalyticSource) return 0;
    if(ForwardedSource){
        if(Notifications>=nWarning) printf("\033[1;33mWARNING:\033[0m Using a source member function does not allow to get any parameters!\n");
        return 0;
    }
    return AnaSourcePar[NumSourcePars+WhichPar];
}
double CATS::GetPotPar(const unsigned& usCh, const unsigned& usPW, const unsigned& WhichPar) const{
    if(NumCh<=usCh || NumPW[usCh]<=usPW){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::GetPotPar(...)\n");
        return 0;
    }
    return PotPar[usCh][usPW][NumPotPars+WhichPar];
}

void CATS::UseExternalWaveFunction(const unsigned& uMomBin, const unsigned& usCh, const unsigned& usPW,
                                 const complex<double>* RadWF, const unsigned& NumRadBins, const double* RadBins, const double& PHASESHIFT){

    if(NumMomBins<=uMomBin || NumCh<=usCh || NumPW[usCh]<=usPW){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m Bad input in CATS::UseExternalWaveFunction(...)\n");
        return;
    }

    ComputedWaveFunction = false;
    ComputedCorrFunction = false;

    if(!RadWF || NumRadBins==0 || !RadBins){
        ExternalWF[uMomBin][usCh][usPW] = NULL;
        NumExtWfRadBins[uMomBin][usCh][usPW] = 0;
        ExtWfRadBins[uMomBin][usCh][usPW] = NULL;
        return;
    }

    ExternalWF[uMomBin][usCh][usPW] = RadWF;
    NumExtWfRadBins[uMomBin][usCh][usPW] = NumRadBins;
    ExtWfRadBins[uMomBin][usCh][usPW] = RadBins;
    PhaseShift[uMomBin][usCh][usPW] = PHASESHIFT;

}

//!Running the analysis
void CATS::KillTheCat(const int& Options){
    //Check if all needed variables are defined
    if(!UseAnalyticSource && !InputFileName){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m The data cannot be loaded! The path is not set!\n\n");
        return;
    }
    //if(UseAnalyticSource && !AnalyticSourceRad && !AnalyticSourceRadCosTh){
    if(UseAnalyticSource && !AnalyticSource && !ForwardedSource){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m The analytic source function is not set!\n\n");
        return;
    }
    if( (!pdgID[0] || !pdgID[1]) && (!UseAnalyticSource || QuantumStatistics==-1) ){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m The PDG IDs of the particles are not set!\n\n");
        return;
    }
    if(!MomBin){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m The momentum bins are not set!\n\n");
        return;
    }
    if(!IpBin){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m The impact parameter bins are not set!\n\n");
        return;
    }
    if(!NumCh){
        if(Notifications>=nError)
            printf("\033[1;31mERROR!\033[0m There is not a single spin-state set!\n\n");
        return;
    }
    if(Notifications>=nAll)
        printf("\033[1;37mKilling the cat!\033[0m\n----------------\n");

    switch(Options){
    case kSourceChanged:
        LoadedData *= UseAnalyticSource;
        SourceGridReady = false;
        SourceUpdated = false;
        ComputedCorrFunction = false;
        break;
    case kPotentialChanged:
        ComputedWaveFunction = false;
        ComputedCorrFunction = false;
        break;
    case kAllChanged:
        LoadedData *= UseAnalyticSource;
        SourceGridReady = false;
        SourceUpdated = false;
        ComputedWaveFunction = false;
        ComputedCorrFunction = false;
        break;
    default: break;
    }

    bool TotWaveFunEvaluated = SourceGridReady*ComputedWaveFunction;
    //if the poor man renormalization is used, the wave function needs to be reevaluated each time
    if(PoorManRenorm!=0 && UseAnalyticSource==false) TotWaveFunEvaluated*=ComputedCorrFunction;
    bool ReallocateTotWaveFun = !SourceGridReady;

    //in case we have a cut on the total momentum, but the array to save it is not present
    //than the data needs to be reloaded
    if(UseTotMomCut && !TotalPairMomentum) {LoadedData=false; SourceGridReady=false;}

    short LoadExitCode;

    if(Notifications>=nAll) printf("\033[1;37m Stage 1:\033[0m Obtaining the source...");
    if( (!LoadedData) && !UseAnalyticSource ){
        if(Notifications>=nAll) printf(" Loading from the data-file...\n");
        LoadExitCode = LoadData();
        if(LoadedData && LoadExitCode){
            if(LoadExitCode==1){
                if(Notifications>=nAll)
                    printf("          \033[1;37mLoading status:\033[1;32m Success!\033[1;37m\n"
                    "          Number of pairs loaded: %u\033[0m\n",NumPairs);
            }
            else{
                if(Notifications>=nWarning)
                    printf("          \033[1;37mLoading status:\033[1;33m Finished with warnings!\033[1;37m\n");
                if(Notifications>=nAll)
                    printf("          Number of pairs loaded: %u\033[0m\n",NumPairs);
            }
        }
        else{
            if(Notifications>=nError)
                printf("          \033[1;37mLoading status:\033[1;31m Failed!\033[0m\n\n");
            return;
        }
    }
    else if(Notifications>=nAll){
        printf("\n");
    }

    if(Notifications>=nAll){
        if(UseAnalyticSource){
            printf("          Using analytic source.\n");
        }
        else{
            printf("          \033[0mUsing data-defined source.\n");
        }
    }

    if(Notifications>=nAll)
        printf("\033[1;37m Stage 2:\033[0m Setup the computing grid...\n");
    if(!SourceGridReady){
        SetUpSourceGrid();
    }
    if(!GamowCorrected){
        UpdateCPF();
    }
    if(Notifications>=nAll)
        printf("          \033[1;32mDone!\033[0m\n");

    if(Notifications>=nAll)
        printf("\033[1;37m Stage 3:\033[0m Solving the Schroedinger equation...\n");
    if(!ComputedWaveFunction) ComputeWaveFunction();
    if(Notifications>=nAll)
        printf("          \033[1;32mDone!\033[0m\n");

    if(Notifications>=nAll)
        printf("\033[1;37m Stage 4:\033[0m Computing the total wave function...\n");
    if(!TotWaveFunEvaluated) ComputeTotWaveFunction(ReallocateTotWaveFun);
    if(Notifications>=nAll)
        printf("          \033[1;32mDone!\033[0m\n");

    if(Notifications>=nAll)
        printf("\033[1;37m Stage 5:\033[0m Computing the correlation function...\n");
    if(!ComputedCorrFunction) FoldSourceAndWF();
    if(Notifications>=nAll)
        printf("          \033[1;32mDone!\033[0m\n");

    if(Notifications>=nAll) printf("\n");
}
bool CATS::CkStatus(){
    return ComputedCorrFunction;
}
bool CATS::SourceStatus(){
    return SourceUpdated;
}
bool CATS::PotentialStatus(){
    return ComputedWaveFunction;
}

void CATS::ComputeTheRadialWaveFunction(){
    if(!NumMomBins) {if(Notifications>=nError)printf("\033[1;31mERROR:\033[0m The momentum bins are not defined!\n"); return;}
    ComputeWaveFunction();
}

void CATS::ComputeWaveFunction(){

    if(!MomBinConverged){
        MomBinConverged = new bool [NumMomBins];
    }
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        MomBinConverged[uMomBin] = true;
    }

    //here one finds the max num PWs, this is needed later on for
    //the correct mapping of all variables
    unsigned short MaxNumPW=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        if(NumPW[usCh]>MaxNumPW){
            MaxNumPW = NumPW[usCh];
        }
    }
    unsigned TotalNumberOfBins = NumMomBins*NumCh*MaxNumPW;

    unsigned uMomBin;
    unsigned short usCh;
    unsigned short usPW;

    //#pragma omp parallel for private(uMomBin,usCh,usPW)
    //the problem with the omp is that the same PotPar are used, we need separate instance for each thread if we want it to work
    //this is however difficult as currently I only pass a single pointer as PotPar, and I do not know how many arguments there are
    //the same problem should occur for the source
    for(unsigned uMPP=0; uMPP<TotalNumberOfBins; uMPP++){
        //compute to which MomBin, Polarization and PW corresponds this MPP,
        //i.e. map uMomBin, usCh and usPW to uMPP
        uMomBin = uMPP/(NumCh*MaxNumPW);
        usCh = (uMPP%(NumCh*MaxNumPW))/(MaxNumPW);
        usPW = (uMPP%(NumCh*MaxNumPW))%(MaxNumPW);
        //since uMPP is build under the assumption all NumPW are the same
        //one has to check if usPW has a meaningful value!
        if(usPW>=NumPW[usCh]) continue;
        //if the potential for this channel and PW is zero => continue
        if(!ShortRangePotential[usCh][usPW] && !ExternalWF[uMomBin][usCh][usPW]) continue;
        //skip momentum bins which have obtained an error code
        if(!MomBinConverged[uMomBin] && ExcludeFailedConvergence) continue;

        SavedWaveFunBins[uMomBin][usCh][usPW]=0;

        //if s+l is odd, than this partial wave will cancel out during the
        //symmetrization for identical particles
        if( IdenticalParticles && (usPW+Spin[usCh])%2 ) continue;

        double Momentum;
        //the momentum is taken from the center of the bin
        Momentum = GetMomentum(uMomBin);
        //Momentum = 0.5*(MomBin[uMomBin]+MomBin[uMomBin+1]);

        //perform the numerical computation
        if(!ExternalWF[uMomBin][usCh][usPW]){
            double* BufferWaveFunction;
            double* BufferRad;

            int q1q2 = Gamow?0:Q1Q2;
            unsigned NumComputedPoints = 2;//counting the initial two starting points
            //index that keeps track of which value in the arrays correspond to the previous, current and next computing step.
            short kOld;
            short kCurrent;
            short kNew;
            double WaveFun[3];
            double PosRad[3];
            double Rho[3];
            //value of the propagating function
            double PropFunVal[3];
            //the value of the prop. function without strong interaction.
            //in case it is equal to the PropFunVal, than the algorithm must have converged
            double PropFunWithoutSI[3];
            //Step size
            double DeltaRad[3];
            //Step size^2
            double DeltaRad2[3];

            double MaxDeltaRad;
            double MinDeltaRad;

            PropagatingFunction(PropFunWithoutSI[0], PropFunVal[0], StartRad, Momentum, usPW, usCh);
            MinDeltaRad = sqrt(fabs(EpsilonProp/(PropFunVal[0]+1e-64)));
            MaxDeltaRad = sqrt(EpsilonProp/(Momentum*Momentum));

            if(MinDeltaRad>MaxDeltaRad) MinDeltaRad=MaxDeltaRad;

            kOld=0;
            kCurrent=1;
            kNew=2;

            //set the initial step in r
            //!This definition of DeltaRad is used in other functions as well.
            //!It is vital that in case there is a need to change the definition, it is changed
            //!in all other functions as well!
            //DeltaRad[kOld] = RhoStep/Momentum;
            //DeltaRad2 = DeltaRad*DeltaRad;

            PosRad[kOld] = 0;
            Rho[kOld] = Momentum*PosRad[kOld];
            DeltaRad[kOld] = MinDeltaRad;
            DeltaRad2[kOld] = DeltaRad[kOld]*DeltaRad[kOld];
            PropagatingFunction(PropFunWithoutSI[kOld], PropFunVal[kOld], PosRad[kOld], Momentum, usPW, usCh);

            PosRad[kCurrent] = PosRad[kOld]+DeltaRad[kOld];
            Rho[kCurrent] = Momentum*PosRad[kCurrent];
            DeltaRad[kCurrent] = MinDeltaRad;
            DeltaRad2[kCurrent] = DeltaRad[kCurrent]*DeltaRad[kCurrent];
            PropagatingFunction(PropFunWithoutSI[kCurrent], PropFunVal[kCurrent], PosRad[kCurrent], Momentum, usPW, usCh);

            //the initial values for the wave function are set based on the solution without the strong potential.
            //this will of course lead to a wrong normalization in the asymptotic region, but this will be corrected for later on,
            //at this stage it is only important that the algorithm gets a meaningful guess so that we do not encounter
            //overflow problems during the calculation.
            WaveFun[kOld] = ReferencePartialWave(PosRad[kOld], Momentum, usPW, q1q2);
            WaveFun[kCurrent] = ReferencePartialWave(PosRad[kCurrent], Momentum, usPW, q1q2);

            bool Convergence = false;
            bool Converged = false;
            //at which point the convergence criteria first occurred
            double ConvergenceRadius=0;
            //how much after the convergence was the WF propagated
            //(this is needed for the normalization)
            const double ConvIntervalLength=3.14;
            const double ConvIntervalRadLength=ConvIntervalLength/Momentum;
            unsigned ConvIntervalSteps=0;

            unsigned StepOfMaxConvergedNumWF=0;
            double MaxConvergedNumWF=0;
            double MaxConvRho=0;
            //the last radius at which the computation was performed
            double MaxConvRad=0;
            double DeltaRadAtMaxConvPrev;
            double DeltaRadAtMaxConv;
            double ConvPointOldWeight=-1;
            double ConvPointWeight=0;

            const unsigned short MinConvSteps = 32;

            //the next two variables are estimates of the max. required bin numbers in case one fails to converge
            //MaxNumRadSteps assumes the worst case scenario -> we compute always with the worst possible step length
            //the second one is the optimistic -> we always compute with the best possible step length
            //bins that converge have usually less than MinNumRadSteps entries, bins that fail to converge a bit above.
            //EstNumRadSteps is a fair small over-estimation of the actual bins needed. It is computed by assuming we have the worst-case
            //scenario up until rho=1.57 or r=4 fm (whichever occurs first at this Momentum) and assumes that the rest as best case scenario.
            //This number looks to be c.a. 10x bigger than the actual num. bins, so use it for memory allocation.

            unsigned EstNumRadSteps = ceil(( (MaxRad*Momentum)>MaxRho?MaxRad-1.57/Momentum:MaxRho/Momentum+ConvIntervalRadLength-1.57/Momentum)/MaxDeltaRad)+MinConvSteps+1 +
                                        ceil(( 4.*FmToNu<1.57/Momentum?4.*FmToNu:1.57/Momentum )/MinDeltaRad);

            BufferWaveFunction = new double [EstNumRadSteps];
            BufferRad = new double [EstNumRadSteps];

            BufferWaveFunction[0] = WaveFun[kOld];
            BufferWaveFunction[1] = WaveFun[kCurrent];

            BufferRad[0] = PosRad[kOld];
            BufferRad[1] = PosRad[kCurrent];

            //!this is the main loop, which computes each next WF point until the result converges or
            //!a maximum value of rho is reached.
            //N.B. if the result is currently converging, the numerical method should not be interrupted even
            //if the MaxRad condition is met!
            while( (!Converged||(PosRad[kCurrent]<ConvergenceRadius+ConvIntervalRadLength))
                  && (PosRad[kOld]<MaxRad || Rho[kOld]<MaxRho || Converged) ){

                PosRad[kNew] = PosRad[kCurrent] + DeltaRad[kCurrent];
                Rho[kNew] = PosRad[kNew]*Momentum;

                WaveFun[kNew] = WaveFun[kCurrent]*(1.+DeltaRad[kCurrent]/DeltaRad[kOld]) -
                                WaveFun[kOld]*DeltaRad[kCurrent]/DeltaRad[kOld] +
                                PropFunVal[kCurrent]*WaveFun[kCurrent]*DeltaRad2[kCurrent];

                //if we run out of memory...
                if(NumComputedPoints>=EstNumRadSteps){
                    EstNumRadSteps*=2;

                    double* BufferTempWF = new double [EstNumRadSteps];
                    for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                        BufferTempWF[uTmp] = BufferWaveFunction[uTmp];
                    }
                    delete [] BufferWaveFunction;
                    BufferWaveFunction = BufferTempWF;

                    double* BufferTempRad = new double [EstNumRadSteps];
                    for(unsigned uTmp=0; uTmp<NumComputedPoints-1; uTmp++){
                        BufferTempRad[uTmp] = BufferRad[uTmp];
                    }
                    delete [] BufferRad;
                    BufferRad = BufferTempRad;
                }

                BufferWaveFunction[NumComputedPoints] = WaveFun[kNew];

                PropagatingFunction(PropFunWithoutSI[kNew], PropFunVal[kNew], PosRad[kNew], Momentum, usPW, usCh);

                DeltaRad2[kNew] = EpsilonProp/(fabs(PropFunVal[kNew])+1e-64);
                DeltaRad[kNew] = sqrt(DeltaRad2[kNew]);
                if(DeltaRad[kNew]<MinDeltaRad){
                    DeltaRad[kNew]=MinDeltaRad;
                    DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
                }
                else if(DeltaRad[kNew]>MaxDeltaRad){
                    DeltaRad[kNew]=MaxDeltaRad;
                    DeltaRad2[kNew]=DeltaRad[kNew]*DeltaRad[kNew];
                }

                BufferRad[NumComputedPoints] = PosRad[kNew];

                Convergence =   fabs( (PropFunWithoutSI[kOld]-PropFunVal[kOld])/(fabs(PropFunWithoutSI[kOld]+PropFunVal[kOld])+1e-64) ) < EpsilonConv &&
                                fabs( (PropFunWithoutSI[kCurrent]-PropFunVal[kCurrent])/(fabs(PropFunWithoutSI[kCurrent]+PropFunVal[kCurrent])+1e-64) ) < EpsilonConv &&
                                fabs( (PropFunWithoutSI[kNew]-PropFunVal[kNew])/(fabs(PropFunWithoutSI[kNew]+PropFunVal[kNew])+1e-64) ) < EpsilonConv;

                //in case we have detected a Convergence
                //1) make sure the convergence is real and not a local artifact of the potential
                //2) reset all variables used to characterize the convergence region
                if( Convergence && !Converged){
                    Converged = true;
                    ConvergenceRadius = PosRad[kNew];

                    //2):
                    MaxConvergedNumWF = 0;
                    MaxConvRho = 0;
                    MaxConvRad = 0;
                    StepOfMaxConvergedNumWF = 0;
                    DeltaRadAtMaxConvPrev = MinDeltaRad;
                    DeltaRadAtMaxConv = MinDeltaRad;
                }
                //this condition makes sure that 1) is fulfilled.
                else{
                    Converged = Convergence;
                }

                //this will reset the ConvIntervalSteps in case the convergence is reset
                ConvIntervalSteps *= Converged;
                //this will count the number of steps computed after convergence
                ConvIntervalSteps += Converged;

                ConvPointWeight = (PosRad[kNew]-ConvergenceRadius)/ConvIntervalRadLength;

                if(fabs(MaxConvergedNumWF)*ConvPointOldWeight<=fabs(WaveFun[kNew])*ConvPointWeight){
                    MaxConvergedNumWF = WaveFun[kNew];
                    MaxConvRho = Rho[kNew];
                    MaxConvRad = PosRad[kNew];
                    StepOfMaxConvergedNumWF = NumComputedPoints;
                    DeltaRadAtMaxConvPrev = DeltaRad[kCurrent];
                    DeltaRadAtMaxConv = DeltaRad[kNew];
                    ConvPointOldWeight = ConvPointWeight;
                }

                NumComputedPoints++;

                kNew++;
                kNew=kNew%3;
                kCurrent++;
                kCurrent=kCurrent%3;
                kOld++;
                kOld=kOld%3;
            }//while(!Converged && PosRad[kOld]<MaxRad)

            //this is not desired for the later computation
            if(StepOfMaxConvergedNumWF==NumComputedPoints-1){
                StepOfMaxConvergedNumWF--;
                MaxConvergedNumWF = BufferWaveFunction[StepOfMaxConvergedNumWF];
                MaxConvRho -= DeltaRadAtMaxConvPrev*Momentum;
                MaxConvRad -= DeltaRadAtMaxConvPrev;
                DeltaRadAtMaxConv = DeltaRadAtMaxConvPrev;
            }
    //! SOMETIMES AT HIGH MOMENTA THE RESULT SIMPLY DOES NOT CONVERGE TO THE REQUIRED VALUE
    //за това дай на юзера опцията да реши какво да прави с неконвергирали резултати
    //1) изхвърли ги
    //2) запаза ги, като нормировката я направи на база максимума в MaxRho-3.14 (btw. make 3.14 the min. value for MaxRho!!!)

            //in case the method failed to converge, the whole momentum bin is marked as unreliable and no
            //further computations are done. This point will be excluded from the final result
            if(!Converged){
                MomBinConverged[uMomBin] = false;
                if(Notifications>=nWarning)
                    printf("          \033[1;33mWARNING:\033[0m A momentum bin at %.2f MeV failed to converge!\n", Momentum);
                if(ExcludeFailedConvergence && Notifications>=nWarning){
                    printf("                   It will be excluded from the final result!\n");
                }
            }

            //if the maximum wave-function is zero the computation will fail!
            //By design this should not really happen.
            if(!MaxConvergedNumWF && Notifications>=nWarning){
                printf("\033[1;33mWARNING:\033[0m MaxConvergedNumWF is zero, which is not allowed and points to a bug in the code!\n");
                printf("         Please contact the developers and do not trust your current results!\n");
            }

            //!Now follows the normalization of the numerical wave function to the asymptotic solution
            //the individual steps are explained in detail in the official CATS documentation
            if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence){
                double ShiftRadStepLen = 0.1/Momentum;
                const unsigned MaximumShiftIter = 121;

                double ShiftRad;

                double DownShift=0;
                double UpShift=ShiftRadStepLen;
                if((!BufferWaveFunction[StepOfMaxConvergedNumWF] || !BufferWaveFunction[StepOfMaxConvergedNumWF+1])&&Notifications>=nWarning){
                    printf("\033[1;33mWARNING:\033[0m BufferWaveFunction is zero, which is not allowed and points to a bug in the code!\n");
                    printf("         Please contact the developers and do not trust your current results!\n");
                }
                double NumRatio = BufferWaveFunction[StepOfMaxConvergedNumWF+1]/BufferWaveFunction[StepOfMaxConvergedNumWF];
                double DownValue;
                double UpValue;
                double SignProduct;
                double SignRefAtLimits;

                //CHECK IF ZERO IS A VIABLE OPTION!
                DownShift = -0.5*ShiftRadStepLen;
                if(fabs(DownShift)>MaxConvRad) DownShift = -MaxConvRad;
                UpShift = 0.5*ShiftRadStepLen;

                //a potential solution should be locked in a region where the ShiftedReferenceWave*NumericalSol have the same sign.
                //furthermore the ShiftedReferenceWave cannot possibly change sign in this region
                for(unsigned uShift=0; uShift<MaximumShiftIter; uShift++){
                    //positive side
                    CurrentRhoStep = DeltaRadAtMaxConv*Momentum;
                    DownValue = AsymptoticRatio(MaxConvRad+DownShift, Momentum, usPW, q1q2) - NumRatio;
                    UpValue = AsymptoticRatio(MaxConvRad+UpShift, Momentum, usPW, q1q2) - NumRatio;

                    if(DownValue*UpValue<0){
                        ShiftRad = 0.5*(DownShift+UpShift);
                        SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)*MaxConvergedNumWF;
                        SignRefAtLimits = ReferencePartialWave(MaxConvRad+DownShift, Momentum, usPW, q1q2)*ReferencePartialWave(MaxConvRad+UpShift, Momentum, usPW, q1q2);
                        if(SignProduct>0 && SignRefAtLimits>0) break;
                    }

                    //negative side
                    if(UpShift<=MaxConvRad){
                        DownValue = AsymptoticRatio(MaxConvRad-DownShift, Momentum, usPW, q1q2) - NumRatio;
                        UpValue = AsymptoticRatio(MaxConvRad-UpShift, Momentum, usPW, q1q2) - NumRatio;
                        if(DownValue*UpValue<0){
                            ShiftRad = -0.5*(DownShift+UpShift);
                            SignProduct = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)*MaxConvergedNumWF;
                            SignRefAtLimits = ReferencePartialWave(MaxConvRad-DownShift, Momentum, usPW, q1q2)*ReferencePartialWave(MaxConvRad-UpShift, Momentum, usPW, q1q2);
                            if(SignProduct>0 && SignRefAtLimits>0){
                                double Temp = DownShift;
                                DownShift = -UpShift;
                                UpShift = -Temp;
                                break;
                            }
                        }
                    }

                    //the factor 0.95 is there just to make sure that we don't by accident
                    //miss a zero just around the limiting values
                    UpShift += 0.95*ShiftRadStepLen;
                    DownShift = UpShift-ShiftRadStepLen;
                }
                ShiftRad = NewtonRapson(&CATS::AsymptoticRatio,
                        DeltaRadAtMaxConv, usPW, Momentum, q1q2, MaxConvRad+DownShift, MaxConvRad+UpShift, NumRatio) - MaxConvRad;
                double ShiftRho = ShiftRad*Momentum;
                double Norm = ReferencePartialWave(MaxConvRad+ShiftRad, Momentum, usPW, q1q2)/MaxConvergedNumWF;

                PhaseShift[uMomBin][usCh][usPW] = ShiftRho;
                PhaseShiftF[usCh][usPW][uMomBin] = ShiftRho;

                //we save all entries up to the point in which the normalization was performed. For higher values
                //we will later on use the asymptotic form. N.B. in principle one could save a bit of space and
                //save the result at the first moment of detected convergence, however than in case the convergence
                //was not "perfect" there might be some inaccuracy in the interval up to the normalization point.
                //btw. the step size is set to be the MaxDeltaRad

                SavedWaveFunBins[uMomBin][usCh][usPW] = StepOfMaxConvergedNumWF+1;
                unsigned& SWFB = SavedWaveFunBins[uMomBin][usCh][usPW];

                if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
                WaveFunRad[uMomBin][usCh][usPW] = new double [SWFB+1];

                if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
                WaveFunctionU[uMomBin][usCh][usPW] = new complex<double> [SWFB];

                for(unsigned uPoint=0; uPoint<SWFB; uPoint++){
                    WaveFunctionU[uMomBin][usCh][usPW][uPoint] = Norm*CPF[uMomBin]*BufferWaveFunction[uPoint];
                }

                //we set up those as bin-range (i.e. the middle of the bin should be buffer rad)
                //we start from r=0, than proceed by adding the limit between different bins as the mid-way to the next point
                WaveFunRad[uMomBin][usCh][usPW][0] = 0;
                for(unsigned uPoint=1; uPoint<SWFB; uPoint++){
                    WaveFunRad[uMomBin][usCh][usPW][uPoint] = (BufferRad[uPoint]+BufferRad[uPoint-1])*0.5;
                }
                //the very last point we simply set as the double distance between the last bin-limit and the last radius value
                WaveFunRad[uMomBin][usCh][usPW][SWFB] = 2*BufferRad[SWFB-1]-WaveFunRad[uMomBin][usCh][usPW][SWFB-1];

            }//if(MomBinConverged[uMomBin] || !ExcludeFailedConvergence)

            delete [] BufferRad;
            delete [] BufferWaveFunction;
        }//end of the numerical computation
        //the case with external wave function
        else{
            SavedWaveFunBins[uMomBin][usCh][usPW] = NumExtWfRadBins[uMomBin][usCh][usPW];
            unsigned& SWFB = SavedWaveFunBins[uMomBin][usCh][usPW];

            if(WaveFunRad[uMomBin][usCh][usPW]) delete [] WaveFunRad[uMomBin][usCh][usPW];
            WaveFunRad[uMomBin][usCh][usPW] = new double [SWFB+1];

            if(WaveFunctionU[uMomBin][usCh][usPW]) delete [] WaveFunctionU[uMomBin][usCh][usPW];
            WaveFunctionU[uMomBin][usCh][usPW] = new complex<double> [SWFB];

            WaveFunRad[uMomBin][usCh][usPW][0] = 0;
            for(unsigned uPoint=1; uPoint<=SWFB; uPoint++){
                //the input from outside is supposed to be in fermi, hence to conversion
                WaveFunRad[uMomBin][usCh][usPW][uPoint] = ExtWfRadBins[uMomBin][usCh][usPW][uPoint]*FmToNu;
            }

            for(unsigned uPoint=0; uPoint<SWFB; uPoint++){
                WaveFunctionU[uMomBin][usCh][usPW][uPoint] = CPF[uMomBin]*ExternalWF[uMomBin][usCh][usPW][uPoint]*FmToNu;
            }
        }
    }//for(unsigned uMPP=0; uMPP<TotalNumberOfBins; uMPP++)
    ComputedWaveFunction = true;

}

void CATS::ComputeTotWaveFunction(const bool& ReallocateTotWaveFun){
    if(!BaseSourceGrid) return;

    if(WaveFunction2 && ReallocateTotWaveFun){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                delete [] WaveFunction2[uMomBin][uGrid];
            }
            delete [] WaveFunction2[uMomBin];
        }
        delete [] WaveFunction2; WaveFunction2 = NULL;
    }

    NumGridPts = BaseSourceGrid->GetNumEndNodes();

    if(!WaveFunction2){
        WaveFunction2 = new double** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            WaveFunction2[uMomBin] = new double* [NumGridPts];
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                WaveFunction2[uMomBin][uGrid] = new double [NumCh];
            }
        }
    }

    //double Momentum;
    double Radius;
    double CosTheta;

    double Time;
    int pTotal;
    int pTotalOld;
    double Progress;
    char* cdummy = new char [512];
    double TotalSteps = double(NumMomBins)*double(NumGridPts);
    double CurrentStep;
    DLM_Timer dlmTimer;

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        //Momentum = GetMomentum(uMomBin);
        for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
            if(PoorManRenorm!=1 && UseAnalyticSource==false) Radius = BaseSourceGrid->GetParValue(uGrid, 0)*FmToNu*PoorManRenorm;
            else Radius = BaseSourceGrid->GetParValue(uGrid, 0)*FmToNu;
            CosTheta = BaseSourceGrid->GetParValue(uGrid, 1);
            for(unsigned usCh=0; usCh<NumCh; usCh++){
                WaveFunction2[uMomBin][uGrid][usCh] =   ThetaDependentSource?
                                                        EffectiveFunctionTheta(uMomBin, Radius, CosTheta, usCh):
                                                        EffectiveFunction(uMomBin, Radius, usCh);
            }
            CurrentStep = double(uMomBin)*double(NumGridPts)+double(uGrid+1);
            Progress = CurrentStep/TotalSteps;
            pTotal = int(Progress*100);
            if(pTotal!=pTotalOld){
                Time = double(dlmTimer.Stop())/1000000.;
                Time = round((1./Progress-1.)*Time);
                ShowTime((long long)(Time), cdummy, 2, true, 5);
                if(Notifications>=nAll) printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
                cout << flush;
                pTotalOld = pTotal;
            }
        }
    }
    if(Notifications>=nAll) printf("\r\033[K");
}

//! N.B. the units in this function (until the result is saved) are fm and GeV!!!
//Exit codes: -1 = WARNING; 0 = ERROR; 1 = OKAY
short CATS::LoadData(const unsigned short& NumBlankHeaderLines){
    DLM_Timer dlmTimer;
    double Time;

    bool ProgressBar=false;
    LoadedData = false;
    SourceGridReady = false;
    SourceUpdated = false;
    char* cdummy = new char [512];

    unsigned MaxTotNumPairs = MaxPairsPerBin*NumMomBins*NumIpBins;

    if(!RelativeMomentum){
        RelativeMomentum = new double [MaxTotNumPairs];
        RelativePosition = new double [MaxTotNumPairs];
        RelativeCosTheta = new double [MaxTotNumPairs];
        PairMomBin = new unsigned [MaxTotNumPairs];
        PairIpBin = new unsigned [MaxTotNumPairs];
        GridBoxId = new int64_t [MaxTotNumPairs];
    }

    if(UseTotMomCut && !TotalPairMomentum){
        TotalPairMomentum = new double [MaxTotNumPairs];
    }

    if(!LoadedPairsPerBin){
        LoadedPairsPerBin = new unsigned* [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            LoadedPairsPerBin[uMomBin] = new unsigned [NumIpBins];
        }
    }
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            LoadedPairsPerBin[uMomBin][uIpBin] = 0;
        }
    }

    if(!LoadedPairsPerMomBin)
        LoadedPairsPerMomBin = new unsigned [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        LoadedPairsPerMomBin[uMomBin] = 0;
    }

    if(!WeightIp) WeightIp = new double [NumIpBins];
    if(!WeightIpError) WeightIpError = new double [NumIpBins];

    FILE *InFile;
    InFile = fopen(InputFileName, "r");
    if(!InFile){
        if(Notifications>=nError)
            printf("          \033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName);
        return 0;
    }

    fseek ( InFile , 0 , SEEK_END );
    long EndPos;
    EndPos = ftell (InFile);
    fseek ( InFile , 0 , SEEK_SET );
    long CurPos;

    //Read the header lines
    for(unsigned short us=0; us<NumBlankHeaderLines; us++){
        if(!fgets(cdummy, 255, InFile)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
            continue;
        }
    }

    if(feof(InFile)){
        if(Notifications>=nError){
            printf("\033[1;31m          ERROR:\033[0m Trying to read past end of file %s\n", InputFileName);
            printf("         No particle pairs were loaded :(\n");

        }
        return 0;
    }

    NumPairs=0;
    unsigned NumTotalPairs=0;

    int EventNumber;
    int NumPartInEvent;
    double ImpPar;
    double fDummy;

    unsigned* NumSePairsIp = new unsigned [NumIpBins];
    unsigned TotalNumSePairs=0;
    unsigned* NumEvents = new unsigned [NumIpBins];
    unsigned TotNumEvents = 0;
    //unsigned* NumEvPart = new unsigned [NumIpBins];

    CatsParticle KittyParticle;
    CatsDataBuffer** KittyBuffer = new CatsDataBuffer* [NumIpBins];
    unsigned* uBuffer = new unsigned[NumIpBins];
    //CatsDataBuffer KittyBuffer(MixingDepth,pdgID[0],pdgID[1]);

    CatsEvent DummyEvent(pdgID[0],pdgID[1]);
    CatsEvent*** KittyEvent;
    KittyEvent = new CatsEvent** [NumIpBins];

    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        //NumEvPart[uIpBin] = 0;
        NumSePairsIp[uIpBin] = 0;
        NumEvents[uIpBin] = 0;
        uBuffer[uIpBin] = 0;
        KittyBuffer[uIpBin] = new CatsDataBuffer(MixingDepth,pdgID[0],pdgID[1]);
        KittyEvent[uIpBin] = new CatsEvent* [MixingDepth];
        for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
            KittyEvent[uIpBin][uDepth] = new CatsEvent(pdgID[0],pdgID[1]);
        }
    }

    unsigned WhichIpBin;
    unsigned MaxTotPairs = MaxPairsPerBin*NumIpBins*NumMomBins;
//    if(MaxPairsToLoad<MaxTotPairs){
//        MaxTotPairs = MaxPairsToLoad;
//    }

    //progress
    //percentage of the file read. The reading speed should be more or less constant,
    //so this value can always give an accurate maximum ETA
    float pFile;
    //percentage of the required number of pairs found in the file. Unless the file has some special
    //internal structure and the events are saved randomly, the should also be a very accurate
    //estimate of the max ETA
    float pMaxPairsToRead;
    //most of the times this is a bad estimate to the ETA, since the MaxPairsPerBin can change drastically
    //the amount of pairs that is accepted per second, i.e. when the bins with high-statistics are full
//    float pMaxPairsToLoad;
    //this should also give a realistic ETA, if we define it based on the bin with fewest entries
    float pMaxPairsPerBin;
    float pTemp;

    short pTotal;
    short pTotalOld;

    bool SkipThisEvent;
    bool NewInterestingEvent;

    unsigned SelectedSePairs;

    unsigned RejectedHighMultEvents = 0;
    const int HighMultLimit = 128;

    //!---Iteration over all events---
    while(!feof(InFile)){
        SkipThisEvent = false;
        if(NumPairs>=MaxTotPairs) break;
        if(NumTotalPairs>=MaxPairsToRead) break;

        if(!fscanf(InFile,"%i %i %lf %lf",&EventNumber,&NumPartInEvent,&ImpPar,&fDummy)){
            printf("\033[1;33mWARNING!\033[0m Possible bad input-file, error when reading from %s!\n",InputFileName);
            continue;
        }

        NewInterestingEvent = false;

        ImpPar = fabs(ImpPar);
        if(ImpPar<IpBin[0] || ImpPar>IpBin[NumIpBins]) {SkipThisEvent=true; WhichIpBin=0;}
        else WhichIpBin = GetIpBin(ImpPar);
        if(WhichIpBin>=NumIpBins) {SkipThisEvent=true; WhichIpBin=0;}

        TotNumEvents++;
        if(NumPartInEvent>HighMultLimit) RejectedHighMultEvents++;

        //!---Iteration over all particles in this event---
        for(int iPart=0; iPart<NumPartInEvent; iPart++){
            if(NumPairs>=MaxTotPairs) break;
            if(NumTotalPairs>=MaxPairsToRead) break;

            KittyParticle.ReadFromOscarFile(InFile);
            if(TransportRenorm!=1){
                KittyParticle.RenormSpacialCoordinates(TransportRenorm);
            }

            if(NumPartInEvent>HighMultLimit) continue;

            //sometimes one might go beyond the limit of the file
            if(ftell(InFile)>=EndPos)continue;

            if(SkipThisEvent) continue;

            if(KittyParticle.GetE()==0){
                if(Notifications>=nWarning)
                    printf("\033[1;33mWARNING!\033[0m Possible bad input-file, there are particles with zero energy!\n");
                continue;
            }

            if(KittyParticle.GetPid()!=pdgID[0] && KittyParticle.GetPid()!=pdgID[1])
                continue; //don't save this particle if it is of the wrong type

            if(!NewInterestingEvent){
                NewInterestingEvent = true;
                NumEvents[WhichIpBin]++;
            }

            KittyEvent[WhichIpBin][uBuffer[WhichIpBin]]->AddParticle(KittyParticle);

        }//for(int iPart=0; iPart<NumPartInEvent; iPart++)

        KittyEvent[WhichIpBin][uBuffer[WhichIpBin]]->ComputeParticlePairs(TauCorrection);

        KittyBuffer[WhichIpBin]->SetEvent(uBuffer[WhichIpBin], *KittyEvent[WhichIpBin][uBuffer[WhichIpBin]]);

        uBuffer[WhichIpBin]++;

        //if the buffer is full -> empty it!
        //note that if it happens the we leave the while loop before emptying the buffer,
        //uBuffer will be != than zero! use this condition to empty the buffer when exiting the loop!
        if(uBuffer[WhichIpBin]==MixingDepth){
            SelectedSePairs=LoadDataBuffer(WhichIpBin, KittyBuffer[WhichIpBin]);

            for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
                KittyEvent[WhichIpBin][uDepth]->Reset();
            }
            NumTotalPairs+=KittyBuffer[WhichIpBin]->GetNumPairs();
            TotalNumSePairs+=SelectedSePairs;
            NumSePairsIp[WhichIpBin]+=SelectedSePairs;
            uBuffer[WhichIpBin]=0;
        }

        CurPos = ftell (InFile);
        pMaxPairsToRead = double(NumTotalPairs)/double(MaxPairsToRead);//
//        pMaxPairsToLoad = double(NumPairs)/double(MaxTotPairs);
        pFile = double(CurPos)/double(EndPos);//what fraction of the file has been read
        ProgressLoad = pMaxPairsToRead>pFile?pMaxPairsToRead:pFile;
//        ProgressLoad = pMaxPairsToLoad>ProgressLoad?pMaxPairsToLoad:ProgressLoad;

        pMaxPairsPerBin = 1;
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                //select the smallest possible pMaxPairsPerBin
                pTemp = float(LoadedPairsPerBin[uMomBin][uIpBin])/float(MaxPairsPerBin);
                if(pTemp<pMaxPairsPerBin) pMaxPairsPerBin=pTemp;
            }
        }
        ProgressLoad = pMaxPairsPerBin>ProgressLoad?pMaxPairsPerBin:ProgressLoad;

        pTotal = int(ProgressLoad*100);
        if(pTotal!=pTotalOld){
            Time = double(dlmTimer.Stop())/1000000.;
            //EtaPerBin = round((1./pMaxPairsPerBin-1.)*Time);
            //EtaToLoad = round((1./pMaxPairsToLoad-1.)*Time);
            Time = round((1./ProgressLoad-1.)*Time);
            ShowTime((long long)(Time), cdummy, 2, true, 5);
            if(Notifications>=nAll)
                printf("\r\033[K          Progress %3d%%, ETA %s",pTotal,cdummy);
            ProgressBar = true;
            cout << flush;
            pTotalOld = pTotal;
        }

    }//while(!feof(InFile))

    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        //empty the buffer for the last time
        if(uBuffer[uIpBin]){
            //if the buffer needs to be emptied, make sure that we do not have any "leftovers" saved in the buffer!
            for(unsigned uLeftOver=uBuffer[uIpBin]; uLeftOver<MixingDepth; uLeftOver++){
                KittyBuffer[uIpBin]->SetEvent(uLeftOver, DummyEvent);
            }
            SelectedSePairs = LoadDataBuffer(uIpBin, KittyBuffer[uIpBin]);
            NumTotalPairs+=KittyBuffer[uIpBin]->GetNumPairs();
            TotalNumSePairs+=SelectedSePairs;
            NumSePairsIp[uIpBin]+=SelectedSePairs;
        }
    }

    if(ProgressBar && Notifications>=nAll){
        printf("\r\033[K");
    }

    if(RejectedHighMultEvents && Notifications>=nWarning){
        printf("\033[1;33m          WARNING:\033[0m CATS cannot handle very high multiplicity events!\n"
               "                   As a result %u events have been rejected!\n", RejectedHighMultEvents);
        if(double(RejectedHighMultEvents)/double(TotNumEvents)>0.01)
            printf("                   We are sorry for the inconvenience :'(\n"
                   "                   If you get in touch with us we will try to find a fix for you!\n");
    }

    if(!NumPairs){
        if(Notifications>=nWarning)
            printf("\033[1;31m          WARNING:\033[0m There were no pairs loaded! The computation cannot proceed!\n");
    }
    else if(!TotalNumSePairs){
        if(Notifications>=nWarning)
            printf("\033[1;31m          WARNING:\033[0m There were no same-events pairs found! The computation cannot proceed!\n");
    }
    else{
        LoadedData = true;
        if(UseTotMomCut){
            LoadedMinTotPairMom = MinTotPairMom;
            LoadedMaxTotPairMom = MaxTotPairMom;
        }
    }
    if(LoadedData && MixingDepth>1){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            WeightIp[uIpBin] = TotalNumSePairs?double(NumSePairsIp[uIpBin])/double(TotalNumSePairs):0;
            if(NumIpBins==1){
                WeightIpError[uIpBin] = 0;
            }
            else if(NumSePairsIp[uIpBin]){
                WeightIpError[uIpBin] = WeightIp[uIpBin]/sqrt(double(NumSePairsIp[uIpBin]));
            }
            else{
                WeightIpError[uIpBin] = 1000;
            }
        }
    }
    else if(MixingDepth==1){
        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            WeightIp[uIpBin] = 1./double(NumIpBins);
            WeightIpError[uIpBin] = 0;
        }
    }
    fclose(InFile);

/*
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        delete [] pIpMomBin[uMomBin];
    }
    delete [] pIpMomBin;
*/
    for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
        for(unsigned uDepth=0; uDepth<MixingDepth; uDepth++){
            delete KittyEvent[uIpBin][uDepth];
        }
        delete KittyBuffer[uIpBin];
        delete [] KittyEvent[uIpBin];
    }

    delete [] cdummy;

    delete [] NumEvents;
    delete [] NumSePairsIp;

    delete [] KittyBuffer;
    delete [] KittyEvent;
    delete [] uBuffer;

    return RejectedHighMultEvents?-1:1;
}

//!CHECK FOR THE MAX NUM PAIRS!
//return the number of same event pairs that pass our basic selection criteria
unsigned CATS::LoadDataBuffer(const unsigned& WhichIpBin, CatsDataBuffer* KittyBuffer){
    KittyBuffer->GoBabyGo(TauCorrection);

    const CatsParticlePair* PairDif;
    const CatsLorentzVector* PairSum;

    double RelPosCom;
    double RelCosThCom;
    double RelMomCom;
    double RedMomComMeV;

    unsigned WhichMomBin;
    unsigned NumSePairs = KittyBuffer->GetNumPairsSameEvent();
    unsigned GoodSePairs = 0;

    double ParticleVector[2];

    for(unsigned uPair=0; uPair<KittyBuffer->GetNumPairs(); uPair++){
        PairDif = KittyBuffer->GetPair(uPair);
        PairSum = &PairDif->GetSum();
        RelPosCom = PairDif->GetR();
        RelCosThCom = (PairDif->GetPx()*PairDif->GetX()+
                       PairDif->GetPy()*PairDif->GetY()+
                       PairDif->GetPz()*PairDif->GetZ())/
                       (PairDif->GetP()*PairDif->GetR());
        RelMomCom = PairDif->GetP();
        RedMomComMeV = 500.*RelMomCom;

        bool Selected = true;

        if(RelPosCom>MaxRad*NuToFm || RelPosCom!=RelPosCom || RelPosCom==0 || RelMomCom==0){
            Selected = false;
        }

        //only save relevant particle pairs
        if(RedMomComMeV<MomBin[0] || RedMomComMeV>MomBin[NumMomBins]){
            Selected = false;
        }
        else if(Selected){
            WhichMomBin = GetMomBin(RedMomComMeV);
            if(WhichMomBin>=NumMomBins) Selected = false;
        }

        //check the total pair momentum condition
        if(UseTotMomCut && (PairSum->GetP()*1000<MinTotPairMom || PairSum->GetP()*1000>MaxTotPairMom)){
            Selected = false;
        }

        //if at this point the pair is not rejected yet, we count it as a "good" pair, i.e. it would be accepted
        //if not for the MaxPairsPerBin limitation. GoodSePairs is than used to compute the weight of each bin.
        //N.B. WE SHOULD COUNT HERE AND NOT AFTER MaxPairsPerBin, since else one would bias the sample!
        if(uPair<NumSePairs && Selected){
            GoodSePairs++;
        }
        if(Selected){
            if(LoadedPairsPerBin[WhichMomBin][WhichIpBin]>=MaxPairsPerBin){
                Selected = false;
            }
        }
        if(!Selected){
            continue;
        }

        RelativeMomentum[NumPairs] = RedMomComMeV;
        RelativePosition[NumPairs] = RelPosCom*FmToNu;
        RelativeCosTheta[NumPairs] = RelCosThCom;
        if(UseTotMomCut) TotalPairMomentum[NumPairs] = PairSum->GetP()*1000;

        ParticleVector[0] = RelativePosition[NumPairs];
        ParticleVector[1] = RelativeCosTheta[NumPairs];

        PairMomBin[NumPairs] = WhichMomBin;
        PairIpBin[NumPairs] = WhichIpBin;
        GridBoxId[NumPairs] = GetBoxId(ParticleVector);

        LoadedPairsPerBin[WhichMomBin][WhichIpBin]++;
        LoadedPairsPerMomBin[WhichMomBin]++;
        NumPairs++;
    }

    return GoodSePairs;
}

void CATS::FoldSourceAndWF(){
    if(!kbCorrFun){
        kbCorrFun = new double* [NumMomBins];
        kbCorrFunErr = new double* [NumMomBins];
        kCorrFun = new double [NumMomBins];
        kCorrFunErr = new double [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            kbCorrFun[uMomBin] = new double [NumIpBins];
            kbCorrFunErr[uMomBin] = new double [NumIpBins];
            kCorrFun[uMomBin] = 0;
            kCorrFunErr[uMomBin] = 0;
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                kbCorrFun[uMomBin][uIpBin] = 0;
                kbCorrFunErr[uMomBin][uIpBin] = 0;
            }
        }
    }
    else if(!ComputedCorrFunction){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            kCorrFun[uMomBin] = 0;
            kCorrFunErr[uMomBin] = 0;
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                kbCorrFun[uMomBin][uIpBin] = 0;
                kbCorrFunErr[uMomBin][uIpBin] = 0;
            }
        }
    }
    else return;

    //this can happen if we change something about the source, but force CATS to use the same grid.
    //in such a case before folding the date one needs to update the values for the source!
    if(!SourceUpdated) UpdateSourceGrid();

    unsigned NumGridPts;
    double SourceVal;
    double WaveFunVal;
    double Integrand;
    double SourceInt;
    double SourceIntCut;
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        kCorrFun[uMomBin] = 0;
        kCorrFunErr[uMomBin] = 0;

        for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
            //perform the k,b analysis only for a data-source which has at least 2 b-bins
            if(UseAnalyticSource || NumIpBins<=1) continue;

            if(!kbSourceGrid[uMomBin][uIpBin]){
                kbCorrFun[uMomBin][uIpBin] = 0;
                kbCorrFunErr[uMomBin][uIpBin] = 1e6;
                continue;
            }
            kbCorrFun[uMomBin][uIpBin] = 0;
            kbCorrFunErr[uMomBin][uIpBin] = 0;

            NumGridPts = kbSourceGrid[uMomBin][uIpBin]->GetNumEndNodes();
            SourceInt = 0;
            SourceIntCut = 0;
            for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                double Radius = kbSourceGrid[uMomBin][uIpBin]->GetParValue(uGrid, 0)*FmToNu;
                //this should return zero in case of !ThetaDependentSource,
                //maybe worth doing a QA to make sure
                //CosTheta = kbSourceGrid[uMomBin][uIpBin]->GetParValue(uGrid, 1);
                SourceVal = kbSourceGrid[uMomBin][uIpBin]->GetGridValue(uGrid);
                if(!SourceVal) continue;
                SourceInt += SourceVal;
                if(Radius<SourceMinRad || Radius>SourceMaxRad) {continue;}
                SourceIntCut += SourceVal;
                WaveFunVal=0;
                for(unsigned usCh=0; usCh<NumCh; usCh++) WaveFunVal+=WaveFunction2[uMomBin][uGrid][usCh]*ChannelWeight[usCh];
                Integrand = SourceVal*WaveFunVal;
                kbCorrFun[uMomBin][uIpBin] += Integrand;
                kbCorrFunErr[uMomBin][uIpBin] += pow(Integrand*kbSourceGrid[uMomBin][uIpBin]->GetGridError(uGrid),2);
            }//uGrid
            kbCorrFun[uMomBin][uIpBin] *= double(SourceInt)/double(SourceIntCut);
            kbCorrFunErr[uMomBin][uIpBin] = sqrt(kbCorrFunErr[uMomBin][uIpBin])*double(SourceInt)/double(SourceIntCut);

            //in case we use event mixing (data-source have more than a single b-bin), this is how we proceed
            if(MixingDepth>1){
                kCorrFun[uMomBin] += WeightIp[uIpBin]*kbCorrFun[uMomBin][uIpBin];
                kCorrFunErr[uMomBin] += pow(WeightIpError[uIpBin]*kbCorrFunErr[uMomBin][uIpBin],2) +
                                                    pow(WeightIpError[uIpBin]*kbCorrFun[uMomBin][uIpBin],2) +
                                                    pow(WeightIp[uIpBin]*kbCorrFunErr[uMomBin][uIpBin],2);
            }
        }

        if( UseAnalyticSource || NumIpBins<=1 || MixingDepth==1 ){
            if(!kSourceGrid[uMomBin]){
                kCorrFunErr[uMomBin] = 1e6;
            }
            else{
                NumGridPts = kSourceGrid[uMomBin]->GetNumEndNodes();
                SourceInt = 0;
                SourceIntCut = 0;
                for(unsigned uGrid=0; uGrid<NumGridPts; uGrid++){
                    double Radius = kSourceGrid[uMomBin]->GetParValue(uGrid, 0)*FmToNu;
                    //this should return zero in case of !ThetaDependentSource,
                    //maybe worth doing a QA to make sure
                    //CosTheta = kSourceGrid[uMomBin]->GetParValue(uGrid, 1);
                    SourceVal = kSourceGrid[uMomBin]->GetGridValue(uGrid);
                    if(!SourceVal) continue;
                    SourceInt += SourceVal;
                    if(Radius<SourceMinRad || Radius>SourceMaxRad) {continue;}
                    SourceIntCut += SourceVal;
                    WaveFunVal=0;
                    for(unsigned usCh=0; usCh<NumCh; usCh++) WaveFunVal+=WaveFunction2[uMomBin][uGrid][usCh]*ChannelWeight[usCh];
                    Integrand = SourceVal*WaveFunVal;
                    kCorrFun[uMomBin] += Integrand;
                    kCorrFunErr[uMomBin] += pow(Integrand*kSourceGrid[uMomBin]->GetGridError(uGrid),2);
                }//uGrid
                kCorrFun[uMomBin] *= double(SourceInt)/double(SourceIntCut);
                kCorrFunErr[uMomBin] = sqrt(kCorrFunErr[uMomBin])*double(SourceInt)/double(SourceIntCut);
            }
        }
    }
    ComputedCorrFunction = true;
}

void CATS::SortAllData(){
    if(NumPairs<=1) return;
    DLM_MergeSort < int64_t, unsigned > SortTool;
    SortTool.SetData(GridBoxId,NumPairs);
    SortTool.MergeSort();

    SortTool.GetSortedData(GridBoxId,GridBoxId);

    ResortData(RelativeMomentum, SortTool);
    ResortData(RelativePosition, SortTool);
    ResortData(RelativeCosTheta, SortTool);
    if(UseTotMomCut) ResortData(TotalPairMomentum, SortTool);
    ResortData(PairMomBin, SortTool);
    ResortData(PairIpBin, SortTool);
}

void CATS::SetUpSourceGrid(){
    double LIMIT = GridEpsilon;
    if(!LIMIT){
        if(ThetaDependentSource) LIMIT=1./16384.;
        else LIMIT=1./1024.;
    }
    short MAXDEPTH = GridMaxDepth;
    if(!MAXDEPTH){
        if(ThetaDependentSource) MAXDEPTH=10;
        else MAXDEPTH=14;
    }
    if(ThetaDependentSource && MAXDEPTH>12){
        MAXDEPTH = 12;
    }
    short DIM = ThetaDependentSource?2:1;
    unsigned MAXGRIDPTS = uipow(2,MAXDEPTH*DIM);

    double* MEAN = new double[DIM];
    double* LENGTH = new double[DIM];
    MEAN[0] = MaxRad*0.5*NuToFm;
    LENGTH[0] = MaxRad*NuToFm;
    if(ThetaDependentSource){
        MEAN[1] = 0;
        LENGTH[1] = 2;
    }

    if(BaseSourceGrid){
        delete BaseSourceGrid; BaseSourceGrid=NULL;
    }
    if(UseAnalyticSource){
        //for setting up the grid we take the "mean" value of k
        AnaSourcePar[0] = (MomBin[0]+MomBin[NumMomBins])*0.5;
        BaseSourceGrid = new CATSelder(DIM, GridMinDepth, MAXDEPTH, LIMIT, MEAN, LENGTH,
                         this, AnaSourcePar, NULL, 0);
    }
    else if(NumPairs){
        //sorts the whole data according to the GridBoxId (all bins!)
        //this will setup the base grid
        SortAllData();
        BaseSourceGrid = new CATSelder(DIM, GridMinDepth, MAXDEPTH, LIMIT, MEAN, LENGTH,
                                        NULL, NULL, GridBoxId, NumPairs);
    }
    else{
        BaseSourceGrid = NULL;
    }

    int64_t ArrayPosition=0;

    //sorts the data in momentum blocks and sorts them according to the GridBoxId
    for(unsigned uPair=0; uPair<NumPairs; uPair++){
        GridBoxId[uPair] += int64_t(MAXGRIDPTS)*int64_t(PairMomBin[uPair]);
    }
    SortAllData();
    for(unsigned uPair=0; uPair<NumPairs; uPair++){
        GridBoxId[uPair] -= (int64_t(GridBoxId[uPair])/int64_t(MAXGRIDPTS))*int64_t(MAXGRIDPTS);
    }

    //resets kSourceGrid
    if(kSourceGrid){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            if(kSourceGrid[uMomBin]) delete kSourceGrid[uMomBin];
        }
        delete [] kSourceGrid;
    }
    kSourceGrid = new CATSelder* [NumMomBins];

    //sets up kSourceGrid
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        if(UseAnalyticSource){
            kSourceGrid[uMomBin] = new CATSelder(BaseSourceGrid, this, AnaSourcePar, NULL, 0);
        }
        else if(LoadedPairsPerMomBin[uMomBin]){
            kSourceGrid[uMomBin] = new CATSelder(BaseSourceGrid, NULL, NULL,
                                                &GridBoxId[ArrayPosition], LoadedPairsPerMomBin[uMomBin]);
            ArrayPosition += LoadedPairsPerMomBin[uMomBin];
        }
        else{
            kSourceGrid[uMomBin] = NULL;
        }
    }

    if(kbSourceGrid){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                if(kbSourceGrid[uMomBin][uIpBin]) delete kbSourceGrid[uMomBin][uIpBin];
            }
            delete [] kbSourceGrid[uMomBin];
        }
        delete [] kbSourceGrid; kbSourceGrid=NULL;
    }

    //sorts the data in k,b blocks and sorts them according to the GridBoxId
    //this is needed only if we use a data source with and more than 1 b-bin
    if(!UseAnalyticSource && NumIpBins>1){
        ArrayPosition = 0;

        //sorts the data in momentum-b blocks and sorts them according to the GridBoxId
        for(unsigned uPair=0; uPair<NumPairs; uPair++){
            GridBoxId[uPair] += int64_t(MAXGRIDPTS)*int64_t(PairIpBin[uPair])+
                                int64_t(MAXGRIDPTS)*int64_t(NumIpBins)*int64_t(PairMomBin[uPair]);
        }
        SortAllData();
        for(unsigned uPair=0; uPair<NumPairs; uPair++){
            GridBoxId[uPair] -= (int64_t(GridBoxId[uPair])/int64_t(MAXGRIDPTS))*int64_t(MAXGRIDPTS);
        }

        kbSourceGrid = new CATSelder** [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            kbSourceGrid[uMomBin] = new CATSelder* [NumIpBins];
            AnaSourcePar[0] = GetMomentum(uMomBin);
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                if(UseAnalyticSource){
                    kbSourceGrid[uMomBin][uIpBin] = new CATSelder(BaseSourceGrid, this, AnaSourcePar, NULL, 0);
                }
                else if(LoadedPairsPerBin[uMomBin][uIpBin]){
                    kbSourceGrid[uMomBin][uIpBin] = new CATSelder(BaseSourceGrid, NULL, NULL,
                                                                &GridBoxId[ArrayPosition], LoadedPairsPerBin[uMomBin][uIpBin]);
                    ArrayPosition += LoadedPairsPerBin[uMomBin][uIpBin];
                }
                else{
                    kbSourceGrid[uMomBin][uIpBin] = NULL;
                }
            }
        }
    }
    SourceGridReady = true;
    SourceUpdated = true;

    delete [] MEAN;
    delete [] LENGTH;
}

void CATS::UpdateSourceGrid(){
    if(BaseSourceGrid) BaseSourceGrid->Update();
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        if(kSourceGrid && kSourceGrid[uMomBin]) kSourceGrid[uMomBin]->Update();
    }
    if(kbSourceGrid){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            for(unsigned uIpBin=0; uIpBin<NumIpBins; uIpBin++){
                if(kbSourceGrid && kbSourceGrid[uMomBin] && kbSourceGrid[uMomBin][uIpBin])
                    kbSourceGrid[uMomBin][uIpBin]->Update();
            }
        }
    }
    SourceUpdated = true;
}

double CATS::CoulombPotential(const double& Radius){
    return Gamow?0:Q1Q2*AlphaFS/(fabs(Radius)+1e-64);
}
//the differential equation for the Schroedinger equation
void CATS::PropagatingFunction(double& Basic, double& Full,
                                 const double& Radius, const double& Momentum,
                                 const unsigned short& usPW, const unsigned short& usCh){
    //make sure that there is no division by zero by adding 1e-64
    //the Basic result is the Prop.Fun. WITHOUT a short range potential
    Basic = 2*RedMass*CoulombPotential(Radius) + double(usPW)*(double(usPW)+1)/(Radius*Radius+1e-64) - Momentum*Momentum;
    //the Full result is the Prop.Fun. WITH a short range potential
    PotPar[usCh][usPW][0] = Radius*NuToFm; PotPar[usCh][usPW][1] = Momentum;
    //! In principle this should be executed only if ShortRangePotential[usCh][usPW] is defined.
    //! Do note that this function is NEVER called in case this is not true! Make sure that this stays so!
    Full = Basic + 2*RedMass*ShortRangePotential[usCh][usPW](PotPar[usCh][usPW]);
}

double CATS::PlanePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW) const{
    double Rho = Radius*Momentum;
    //if Rho is zero, the gsl function will not work
    if(!Rho){
        return 0;
    }
    //N.B. gsl_sf_bessel_jl are defined for Rho>0, but in principle the bessel functions are symmetric for even l
    //and anti-symmetric for odd l, this is implemented here.
    return Rho>0?(Radius)*gsl_sf_bessel_jl(usPW,Rho):pow(-1,usPW)*(Radius)*gsl_sf_bessel_jl(usPW,-Rho);
}

double CATS::CoulombPartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const{
    double Eta = RedMass*double(q1q2)*AlphaFS/Momentum;
    double Rho = Radius*Momentum;
    double Overflow=0;
    double Result;
    if(Rho==0) return 0;
    gsl_sf_coulomb_wave_F_array (usPW, 1, Eta, fabs(Rho), &Result, &Overflow);

    Result /= Momentum;

    //N.B. gsl_sf_coulomb_wave_F_array are defined for Rho>0, but in principle the Coulomb functions are symmetric for odd l
    //and anti-symmetric for even l. However here I assume that Momentum>0, and if rho is negative so is the Momentum. Thus
    //the final result for u_l should be symmetric for even l and antisymmetric for odd l. This is implemented here.
    if(Rho<0 && usPW%2==1){
        Result = -Result;
    }
    return Result;
}

//note that this function is never corrected for the Gamow factor, as in the asymptotic region one could just take the exact solution
//instead of correcting a plane wave by Gamow factor
double CATS::ReferencePartialWave(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const{
    return q1q2 ? CoulombPartialWave(Radius,Momentum,usPW,q1q2) : PlanePartialWave(Radius,Momentum,usPW);
    //return (Q1Q2&&!Gamow) ? CoulombPartialWave(Radius,Momentum,usPW) : PlanePartialWave(Radius,Momentum,usPW);
}

double CATS::AsymptoticRatio(const double& Radius, const double& Momentum, const unsigned short& usPW, const int& q1q2) const{
    return ReferencePartialWave(Radius+CurrentRhoStep/Momentum, Momentum, usPW, q1q2)/(ReferencePartialWave(Radius, Momentum, usPW, q1q2)+1e-64);
}

double CATS::NewtonRapson(double (CATS::*Function)(const double&, const double&, const unsigned short&, const int&) const,
                          const double& EpsilonX, const unsigned short& usPW, const double& Momentum, const int& q1q2,
                          const double&  xMin, const double&  xMax, const double& fValShift) const{

    const unsigned maxIter = 256;

    double DeltaX;
    double xVal=(xMax+xMin)*0.5;
    double fVal;
    double DeltaF;

    for(unsigned iIter=0; iIter<maxIter; iIter++){
        fVal = (this->*Function)(xVal, Momentum, usPW, q1q2)-fValShift;
        DeltaF = (this->*Function)(xVal+EpsilonX, Momentum, usPW, q1q2)-fValShift - fVal;
        if(!DeltaF){
            DeltaX=1;
            if(Notifications>=nWarning)
                printf("\033[1;33mWARNING:\033[0m Something is fishy with the NewtonRapson solver! Might be a bug! Please contact the developers!\n");
        }
        else DeltaX = -fVal*EpsilonX/DeltaF;
        xVal += DeltaX;
        int counter = 0;
        double fValNew = (this->*Function)(xVal, Momentum, usPW, q1q2)-fValShift;
        while( ( (fabs(fValNew)>fabs(fVal) && counter<16)
              || (xVal<xMin || xVal>xMax) ) ){
            counter++;
            xVal -= DeltaX*pow(2., -counter);
            fValNew = (this->*Function)(xVal, Momentum, usPW, q1q2)-fValShift;
        }
        if(counter==16 && fabs(fValNew)>fabs(fVal)){
            if(Notifications>=nWarning)
                printf("\033[1;33mWARNING:\033[0m The backtracking of the NewtonRapson root-finder failed!\n");
        }
        if(fabs(fValNew)<fabs(DeltaF)){
            return xVal;
        }
    }
    if(Notifications>=nWarning)
        printf("\033[1;33mWARNING:\033[0m The NewtonRapson root-finder failed!\n");
    return xVal;
}

template <class Type> void CATS::ResortData(Type* input, DLM_MergeSort <int64_t, unsigned>& Sorter){
    unsigned NumOfEl = Sorter.GetNumOfEl();
    Type* Temp;
    Temp = new Type[NumOfEl];
    for(unsigned uEl=0; uEl<NumOfEl; uEl++){
        Temp[uEl] = input[Sorter.GetKey()[uEl]];
    }
    for(unsigned uEl=0; uEl<NumOfEl; uEl++){
        input[uEl] = Temp[uEl];
    }
    delete [] Temp;
}

//btw, if the range is outside the limits, the return value will be equal
//to the NumberOfBoxes. Used somewhere else this might lead to potential segmentation faults, so
//make sure to take care of that!
unsigned CATS::GetBoxId(double* particle){
    const short Dim = ThetaDependentSource?2:1;
    double ParentMean[Dim];
    double ParentLen[Dim];

    ParentMean[0] = MaxRad*0.5;
    ParentLen[0] = MaxRad;
    if(ThetaDependentSource){
        ParentMean[1] = 0;
        ParentLen[1] = 2;
    }

    short MAXDEPTH = GridMaxDepth;
    if(!MAXDEPTH){
        if(ThetaDependentSource) MAXDEPTH=10;
        else MAXDEPTH=14;
    }
    if(ThetaDependentSource && MAXDEPTH>12){
        MAXDEPTH = 12;
    }

    short ChildDepth=0;
    unsigned ChildFirstID = 0;
    unsigned ChildLastID = uipow(2,Dim*MAXDEPTH);
    const unsigned NumSubNodes = uipow(2,Dim);
    double ChildMean[Dim];
    double ChildLen[Dim];
    unsigned ChildNumBoxes;

    //we want to divide our total interval in two for each parameter on the grid.
    //in order to keep track in which "quadrant" we are, we introduce a very simple counter WhichPart for each
    //of the parameters, that can only take values 0 or 1. Each time WhichPart[x] is increased to 2, than it is set to zero
    //and WhichPart[x+1] is increased, i.e. we continue to iterate over the next parameter.
    char WhichPart[Dim];
    bool ThisBox;

    while( ChildDepth<MAXDEPTH ){
        ChildNumBoxes = (ChildLastID-ChildFirstID+1)/NumSubNodes;
        //ChildFirstID = ParentFirstID;
        //this is the first node
        for(short sDim=0; sDim<Dim; sDim++){
            WhichPart[sDim] = 0;
            ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25;
            ChildLen[sDim] = ParentLen[sDim]*0.5;
        }
        ChildDepth++;
        for(unsigned uSub=0; uSub<NumSubNodes; uSub++){
            ChildLastID = ChildFirstID+ChildNumBoxes-1;
            ThisBox = true;
            //see if the particles in located in one of the following boxes
            for(short sDim=0; sDim<Dim; sDim++){
                ThisBox *= (particle[sDim]>=ChildMean[sDim]-0.5*ChildLen[sDim] &&
                            particle[sDim]<=ChildMean[sDim]+0.5*ChildLen[sDim]);
            }
            //when we find the correct bin, we break out of the loop
            if(ThisBox){
                break;
            }
            ChildFirstID += ChildNumBoxes;
            for(short sDim=0; sDim<Dim; sDim++){
                WhichPart[sDim] = (WhichPart[sDim]+1)%2;
                ChildMean[sDim] = ParentMean[sDim]-ParentLen[sDim]*0.25+0.5*ParentLen[sDim]*WhichPart[sDim];
                if(WhichPart[sDim]) break;
            }
        }
        //update the values for the parent
        //ParentFirstID = ChildFirstID;
        //ParentLastID = ChildLastID;
        for(short sDim=0; sDim<Dim; sDim++){
            ParentMean[sDim] = ChildMean[sDim];
            ParentLen[sDim] = ChildLen[sDim];
        }
    }
    return ChildFirstID;
}

complex<double> CATS::EvalWaveFunctionU(const unsigned& uMomBin, const double& Radius,
                                const unsigned short& usCh, const unsigned short& usPW, const bool& DivideByR, const bool& Asymptotic) const{

    double Momentum = GetMomentum(uMomBin);
    if(uMomBin>=NumMomBins){
        if(Notifications>=nError)
            printf("\033[1;31mERROR:\033[0m There is a bug inside EvalWaveFunctionU! Contact the developer!");
        return 0;
    }

    double MultFactor = DivideByR?1./(Radius+1e-64):1;

    unsigned& SWFB = SavedWaveFunBins[uMomBin][usCh][usPW];

    unsigned RadBin = GetRadBin(Radius, uMomBin, usCh, usPW);

    if(RadBin<SWFB && !Asymptotic){
        const complex<double>* WFU = WaveFunctionU[uMomBin][usCh][usPW];
        const double* WFR = WaveFunRad[uMomBin][usCh][usPW];
        //the external WF is assumed to be given in fm
        complex<double> Result = EvalBinnedFun(Radius, SWFB, WFR, NULL, WFU);
        if(Result==1e6 && Notifications>=nWarning)
            printf("\033[1;33mWARNING:\033[0m DeltaRad==0, which might point to a bug! Please contact the developers!\n");
        return Result*MultFactor;
    }
    //below the 1st bin ==> assume zero
    else if(RadBin==SWFB+1){
        return 0;
    }
    else{
        return ReferencePartialWave(Radius+PhaseShift[uMomBin][usCh][usPW]/Momentum, Momentum, usPW, Q1Q2)*MultFactor;
    }
}

double CATS::EffectiveFunction(const unsigned& uMomBin, const double& Radius, const unsigned short& usCh){
    complex<double> Result;
    complex<double> OldResult=100;
    double TotalResult=0;
    double Momentum = GetMomentum(uMomBin);
    //here we make the assumption that the individual partial waves are orthogonal to one another!
    //this should be the case in the absence of angular dependence due to the Legendre polynomials
    for(unsigned short usPW=0; usPW<MaxPw; usPW++){
        //wave function symmetrization
        if( IdenticalParticles && (usPW+Spin[usCh])%2 ) continue;
        //numerical solution, no computation result for zero potential
        if(usPW<NumPW[usCh] && (ShortRangePotential[usCh][usPW] || ExternalWF[uMomBin][usCh][usPW])){
        //if(false){
            Result = EvalWaveFunctionU(uMomBin, Radius, usCh, usPW, true);
            //Check this!!! Should it be squared?
            //the integration of Pl itself results in 1/(2l+1), so this should be fine as it is
            TotalResult += double(2*usPW+1)*pow(abs(Result),2);
        }
        else{
            Result = ReferencePartialWave(Radius, Momentum, usPW, Q1Q2)/(Radius+1e-64);
            //Check this!!! Should it be squared?
            //the integration of Pl itself results in 1/(2l+1), so this should be fine as it is
            Result = double(2*usPW+1)*pow(abs(Result),2);
            TotalResult += abs(Result);
            //convergence criteria
            if(usPW>=NumPW[usCh] && abs(OldResult)<1e-7 && abs(Result)<1e-8) break;
            OldResult = Result;
        }
    }

    return TotalResult*(1+IdenticalParticles);
}


double CATS::EffectiveFunction(const unsigned& uMomBin, const double& Radius){
    double TotWF=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        TotWF += EffectiveFunction(uMomBin, Radius, usCh);
    }
    return TotWF;
}

double CATS::EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta, const unsigned short& usCh){
    complex<double> Result;
    complex<double> OldResult=100;
    complex<double> TotalResult=0;
    double Momentum = GetMomentum(uMomBin);

    if(!RefPartWave){
        RefPartWave = new complex<double> [MaxPw];
        SolvedPartWave = new complex<double> [MaxPw];
        LegPol = new double [MaxPw];
    }

    //we set an unrealistic default value to monitor which PW is computed and which not
    for(unsigned short usPW=0; usPW<MaxPw; usPW++){
        RefPartWave[usPW] = 1e6;
        SolvedPartWave[usPW] = 1e6;
        LegPol[usPW] = 1e6;
    }

    for(unsigned short usPW=0; usPW<MaxPw; usPW++){
        //wave function symmetrization
        if( IdenticalParticles && (usPW+Spin[usCh])%2 ) continue;
        if(usPW<NumPW[usCh] && (ShortRangePotential[usCh][usPW] || ExternalWF[uMomBin][usCh][usPW])){
            if(LegPol[usPW]==1e6) LegPol[usPW]=gsl_sf_legendre_Pl(usPW,CosTheta);
            if(SolvedPartWave[usPW]==1e6) SolvedPartWave[usPW]=EvalWaveFunctionU(uMomBin, Radius, usCh, usPW, true);
            Result = double(2*usPW+1)*SolvedPartWave[usPW]*LegPol[usPW];
        }
        else{
            if(LegPol[usPW]==1e6) LegPol[usPW]=gsl_sf_legendre_Pl(usPW,CosTheta);
            if(RefPartWave[usPW]==1e6) RefPartWave[usPW] = ReferencePartialWave(Radius, Momentum, usPW, Q1Q2);
            //please check if a sqrt is needed for 2*l+1. I think not, because above, when we integrate Pl, the integration of Pl itself results in 1/(2l+1)
            Result = pow(i,usPW)*double(2*usPW+1)*RefPartWave[usPW]/(Radius+1e-64)*LegPol[usPW];
            if(usPW>=NumPW[usCh] && abs(OldResult)<3.16e-4 && abs(Result)<1e-4) break;
            OldResult = Result;
        }
        TotalResult += Result;
    }
    return pow(abs(TotalResult),2);
}

double CATS::EffectiveFunctionTheta(const unsigned& uMomBin, const double& Radius, const double& CosTheta){
    double TotWF=0;
    for(unsigned short usCh=0; usCh<NumCh; usCh++){
        TotWF += EffectiveFunctionTheta(uMomBin, Radius, CosTheta, usCh)*ChannelWeight[usCh];
    }
    return TotWF;
}

unsigned CATS::GetBin(const double& Value, const double* Range, const unsigned& NumBins) const{
    if(NumBins<=1) return 0;
    unsigned WhichBin=NumBins/2;
    unsigned BinMod=4;
    unsigned BinStep;
    //makes sure that the Value is in Range. If not, the returned value is either
    //NumBins or NumBins+1, depending on if we have an underflow or overflow
    if(Value<Range[0]) return NumBins;
    if(Value>Range[NumBins-1]) return NumBins+1;
    while(true){
        if(Range[WhichBin]<=Value && Range[WhichBin+1]>=Value){
            return WhichBin;
        }
//!check the logic of this <=, I made it so, since else we crashed on limit values, before that it was <
        else if(Value<=Range[WhichBin]){
            BinStep = NumBins/BinMod;
            WhichBin -= BinStep?BinStep:1;
            BinMod *= 2;
        }
        else{
            BinStep = NumBins/BinMod;
            WhichBin += BinStep?BinStep:1;
            BinMod *= 2;
        }
    }
}

unsigned CATS::GetMomBin(const double& Momentum) const{
    return GetBin(Momentum, MomBin, NumMomBins+1);
}
unsigned CATS::GetIpBin(const double& bVal) const{
    return GetBin(bVal, IpBin, NumIpBins+1);
}
unsigned CATS::GetRadBin(const double& Radius, const unsigned& uMomBin,
                         const unsigned short& usCh, const unsigned short& usPW) const{
    return GetBin(Radius, WaveFunRad[uMomBin][usCh][usPW], SavedWaveFunBins[uMomBin][usCh][usPW]+1);
}

template <class Type> Type CATS::GetBinCenter(const Type* Bins, const unsigned& WhichBin) const{
    return 0.5*(Bins[WhichBin]+Bins[WhichBin+1]);
}

template <class Type> Type CATS::EvalBinnedFun(const double& xVal, const unsigned& NumBins, const double* Bins, const double* BinCent, const Type* Function) const{
    if(xVal<Bins[0] || xVal>Bins[NumBins]) return 0;
    if(NumBins==1) return Function[0];
    unsigned WhichBin = GetBin(xVal,Bins,NumBins+1);

    double Value[3];
    if(BinCent){
        Value[0] = WhichBin?BinCent[WhichBin-1]:-1;
        Value[1] = BinCent[WhichBin];
        Value[2] = WhichBin<(NumBins-1)?BinCent[WhichBin+1]:-1;
    }
    else{
        Value[0] = WhichBin?GetBinCenter(Bins,WhichBin-1):-1;
        Value[1] = GetBinCenter(Bins,WhichBin);
        Value[2] = WhichBin<(NumBins-1)?GetBinCenter(Bins,WhichBin+1):-1;
    }

    double* InterpolRange;
    const Type* FunRange;

    if(Value[0]==-1){
        InterpolRange = &Value[1];
        FunRange = &Function[WhichBin];
    }
    else if(Value[2]==-1){
        InterpolRange = &Value[0];
        FunRange = &Function[WhichBin-1];
    }
    else if(xVal<Value[1]){
        InterpolRange = &Value[0];
        FunRange = &Function[WhichBin-1];
    }
    else if(Value[1]<xVal){
        InterpolRange = &Value[1];
        FunRange = &Function[WhichBin];
    }
    else{//Value[1]==xVal
        return Function[WhichBin];
    }

    if(InterpolRange[1]-InterpolRange[0]){
        return (FunRange[1]*(xVal-InterpolRange[0])-
                FunRange[0]*(xVal-InterpolRange[1]))/
                (InterpolRange[1]-InterpolRange[0]);
    }
    else{
        printf("\033[1;33mWARNING:\033[0m EvalBinnedFun could not properly extrapolate! The output of the Eval functions might be wrong!\n");
        printf("         Please use Get functions to cross-check your output and contact the developers!\n");
        return 1e6;
    }
}

double CATS::EvaluateTheSource(double* Pars) const{
    if(!AnalyticSource && !ForwardedSource){
        if(Notifications>=nError){
            printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported a crash! The source is not defined!\n");
        }
        return 0;
    }
    if(AnalyticSource && ForwardedSource){
        if(Notifications>=nError){
            printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported multiple definitions of the source!\n");
            printf("         This is a bug, please contact the developers!\n");
        }
        return 0;
    }
    if(!Pars){
        if(Notifications>=nError){
            printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported a crash! This is most likely a bug!\n");
            printf("         Please contact the developers, reporting that Pars=NULL\n");
        }
        return 0;
    }
    if(AnalyticSource){
        return AnalyticSource(Pars);
    }
    //this is the case of ForwardedSource
    if(!SourceContext){
        if(Notifications>=nError){
            printf("\033[1;31mERROR:\033[0m EvaluateTheSource reported a crash! This is most likely a bug!\n");
            printf("         Please contact the developers, reporting that SourceContext=NULL\n");
        }
        return 0;
    }
    return ForwardedSource(SourceContext,Pars);
}

void CATS::UpdateCPF(){
    if(CPF){delete[]CPF; CPF=NULL;}
    CPF = new complex<double> [NumMomBins];
    if(Gamow){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            CPF[uMomBin] = GamowCorrection(GetMomentum(uMomBin),RedMass,Q1Q2);
        }
    }
    else{
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            CPF[uMomBin] = 1;
        }
    }
    GamowCorrected = true;
}
