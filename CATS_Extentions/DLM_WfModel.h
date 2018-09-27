#ifndef DLM_WFMODEL_H
#define DLM_WFMODEL_H

#include <complex>

using namespace std;

class CATS;

class DLM_WfModels{
public:

    enum WfModel { ProtonLambdaHaideNLO, KminProtonHaide, KplusProtonHaide, KminProtonTetsuo };


    DLM_WfModels(WfModel Model);
    ~DLM_WfModels();

    void Init(const CATS& Kitty, const char* inFileName);

    complex<double>**** WaveFunctionU;
    double*** PhaseShifts;
    double* RadBins;

    unsigned GetNumRadBins();

private:

    enum ProtonLambdaFiles {pl1S0, pl1P1, pl3S1, pl3P0, pl3P1, pl3P2};
    enum KminProtonFiles {kmpS01, kmpS11, kmpP01, kmpP03, kmpP11, kmpP13};
    enum KplusProtonFiles {kppS11, kppP11, kppP13};

    const WfModel MODEL;
    const unsigned short NumChannels;
    //const unsigned short NumPwPerCh;
    //const unsigned short NumFiles;

    char** InputFileName;

    unsigned NumRadBins;
};




void InitHaidenbauerNLO(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=2, const int& CUTOFF=600);
void InitHaidenbauerKaonMinus(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=2);
void InitHaidenbauerKaonMinus_ver2(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=2);
void InitHaidenbauerKaonPlus(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=2);
void InitTetsuoKaonMinus(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=2);
//TYPE 0 = p k-; //TYPE 1 = n k0
//TYPE 2 = La pi0; //TYPE 3 = Si+ pi-
//TYPE 4 = Si0 pi0; //TYPE 5 = Si- pi+
void InitCatForHaidenbauerKaonProton1(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE=0);
void InitESC08(const char* InputFolder, CATS& Kitty, complex<double>***** WaveFunctionU, double** RadBins, unsigned& NumRadBins,
               const double& MaxMomentum, const int& TYPE);
void CleanHaidenbauer(CATS& Kitty, complex<double>***** WaveFunctionU, double**** PhaseShifts, double** RadBins);



#endif
