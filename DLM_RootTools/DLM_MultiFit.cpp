
#include <iostream>

#include "TH1.h"
#include "TList.h"
#include "HFitInterface.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TH1F.h"
#include "TF1.h"

#include "Fit/Fitter.h"

#include "DLM_MultiFit.h"


using namespace std;

DLM_MultiFit::DLM_MultiFit(){
    num_spectra = 0;
    max_spectra = 0;
    HistoToFit = NULL;
    FitFunction = NULL;
    EqualPar = NULL;
    totalpar = 0;
    ParMap = NULL;
    WMTF1 = NULL;
    DR = NULL;
    BD = NULL;
    Chi2Fun = NULL;
    spec_par = NULL;
    opt = new ROOT::Fit::DataOptions ();
}

DLM_MultiFit::~DLM_MultiFit(){
    ClearMemFitter();
    if(HistoToFit){
        delete [] HistoToFit;
        HistoToFit = NULL; //The spectra are never created here, thus NO need to delete the objects themselves
    }
    if(FitFunction){
        delete [] FitFunction;
        FitFunction = NULL; //The spectra are never created here, thus NO need to delete the objects themselves
    }
    if(EqualPar){
        for(unsigned short i=0; i<max_spectra; i++){
            delete [] EqualPar[i];
            EqualPar[i] = NULL;
        }
        delete [] EqualPar;
        EqualPar = NULL;
    }
    if(ParMap){
        for(unsigned short i=0; i<num_spectra; i++){
            delete [] ParMap[i];
            ParMap[i] = NULL;
        }
        delete [] ParMap;
        ParMap = NULL;
    }
    delete opt;
}

void DLM_MultiFit::AllocMemPar(){
    if(totalpar){
        for(unsigned int i=0; i<totalpar; i++){
            delete [] EqualPar[i];
        }
        delete [] EqualPar;
        totalpar = 0;
    }

    for(unsigned short i=0; i<num_spectra; i++){
        totalpar += FitFunction[i]->GetNpar();
    }

    EqualPar = new bool * [totalpar];
    //initially there are no relations between the parameters, i.e.
    //they are only equall to themselves!
    for(unsigned int i=0; i<totalpar; i++){
        EqualPar[i] = new bool [totalpar];
    }
    ResetParRelations();
}

void DLM_MultiFit::AllocMemFitter(){
    ClearMemFitter();
    WMTF1 = new ROOT::Math::WrappedMultiTF1 * [num_spectra];
    DR = new ROOT::Fit::DataRange[num_spectra];
    BD = new ROOT::Fit::BinData * [num_spectra];
    Chi2Fun = new ROOT::Fit::Chi2Function * [num_spectra];
}

void DLM_MultiFit::ClearMemFitter(){
    for(unsigned short i=0; i<num_spectra; i++){
        if(WMTF1 && WMTF1[i])
            {delete WMTF1[i]; WMTF1[i] = NULL;}
        if(BD && BD[i])
            {delete BD[i]; BD[i] = NULL;}
        if(Chi2Fun && Chi2Fun[i])
            {delete Chi2Fun[i]; Chi2Fun[i] = NULL;}
    }
    if(WMTF1) {delete [] WMTF1; WMTF1 = NULL;}
    if(BD) {delete [] BD; BD = NULL;}
    if(Chi2Fun) {delete [] Chi2Fun; Chi2Fun = NULL;}
    if(DR) {delete [] DR; DR = NULL;}
}

void DLM_MultiFit::CreateParMap(){

    if(ParMap){
        for(unsigned short i=0; i<NumUniquePar; i++){
            delete [] ParMap;
        }
        delete [] ParMap;
    }

    NumUniquePar = 0;

    //I loop over GNS, but need the "coordinates" in UNS as well
    int iiUNS=0;
    int ijUNS=0;

    ParMap = new unsigned int * [num_spectra];
    for(unsigned short iGNS=0; iGNS<totalpar; iGNS++){
        if(!ijUNS){
            ParMap[iiUNS] = new unsigned int [FitFunction[iiUNS]->GetNpar()];
        }
        int jiUNS=0;
        int jjUNS=0;
        for(unsigned short jGNS=0; jGNS<iGNS; jGNS++){
            //if the current par. iGNS is equal to one of the preveous parameters,
            //than there is NO new unique par.
            if(EqualPar[iGNS][jGNS]){
                ParMap[iiUNS][ijUNS] = ParMap[jiUNS][jjUNS];
                goto EndOfLoop;
            }
            //calc. the jGNS in UNS
            jjUNS++;
            if(jjUNS==FitFunction[jiUNS]->GetNpar()){
                jiUNS++;
                jjUNS=0;
            }
        }
        ParMap[iiUNS][ijUNS] = NumUniquePar;
        NumUniquePar++;

        EndOfLoop:
        //calc. the iGNS in UNS
        ijUNS++;
        if(ijUNS==FitFunction[iiUNS]->GetNpar()){
            iiUNS++;
            ijUNS=0;
        }
    }

}

void DLM_MultiFit::ReallocMemSpectra(unsigned short m){

    if(num_spectra > m){
        max_spectra = m;
        num_spectra = max_spectra;
        Printf("WARNING: ReallocMem had to remove some spectra from the list! Potential loss of data "
               "and(or) memory leak!!!");
    }
    if(!m) {ClearSpectra(); return;}
    max_spectra = m;

    TH1F ** histo_temp = NULL;
    TF1 ** fit_temp = NULL;

    if(num_spectra){
        histo_temp = new TH1F * [num_spectra];
        for(unsigned short i=0; i<num_spectra; i++){
            histo_temp[i] = HistoToFit[i];
        }
        delete [] HistoToFit;

        fit_temp = new TF1 * [num_spectra];
        for(unsigned short i=0; i<num_spectra; i++){
            fit_temp[i] = FitFunction[i];
        }
        delete [] FitFunction;
    }
    HistoToFit = new TH1F * [max_spectra];
    FitFunction = new TF1 * [max_spectra];
    spec_par = new double* [max_spectra];
    for(unsigned short i=0; i<max_spectra; i++){
      spec_par[i] = NULL;
    }
    if(num_spectra){
        for(unsigned short i=0; i<num_spectra; i++){
            HistoToFit[i] = histo_temp[i];
            FitFunction[i] = fit_temp[i];
            //for(int ipar=0; ipar<FitFunction[i]->GetNpar(); ipar++){
            //  spec_par[i][ipar]
            //}
        }
        delete [] histo_temp;
        delete [] fit_temp;
    }
}

bool DLM_MultiFit::AddSpectrum(TH1F* histo, TF1* fit){
    if(!histo || !fit){
        return false;
    }
    if(num_spectra==max_spectra){
        ReallocMemSpectra(num_spectra+1);
    }
    HistoToFit[num_spectra] = histo;
    FitFunction[num_spectra] = fit;
    num_spectra++;
    ResetParRelations();
    AllocMemPar();
    return true;
}


void DLM_MultiFit::ClearSpectra(){
    if(HistoToFit){
        delete [] HistoToFit;
        delete [] FitFunction;
        for(int iSpec=0; iSpec<max_spectra; iSpec++) delete [] spec_par[iSpec];
        delete [] spec_par;
        HistoToFit = NULL;
        FitFunction = NULL;
        spec_par = NULL;
        num_spectra = 0;
        max_spectra = 0;
    }
}

void DLM_MultiFit::ResetParRelations(){
    for(unsigned int i=0; i<totalpar; i++){
        for(unsigned int j=0; j<totalpar; j++){
            EqualPar[i][j] = i==j;
        }
    }
}

bool DLM_MultiFit::SetEqualPar(unsigned short spec1, int par1, unsigned short spec2, int par2){
    if(spec1>=num_spectra || spec2>=num_spectra){
        return false;
    }
    if(par1>=FitFunction[spec1]->GetNpar() || par2>=FitFunction[spec2]->GetNpar()){
        return false;
    }

    //the input is in form : the par1-th parameter from the spec1-th spectrum
    //this needs to be translated into the internal numeration of the number of parameters
    unsigned int whichpar1 = 0;
    for(int i=0; i<spec1; i++) whichpar1 += FitFunction[i]->GetNpar();
    whichpar1 += par1;

    unsigned int whichpar2 = 0;
    for(int i=0; i<spec2; i++) whichpar2 += FitFunction[i]->GetNpar();
    whichpar2 += par2;

    SetEqualPar(whichpar1, whichpar2);

    return true;
}

bool DLM_MultiFit::SetEqualPar(unsigned int whichpar1, unsigned int whichpar2){
    //Printf("totalpar = %i", totalpar);
    if(whichpar1>=totalpar || whichpar2>=totalpar) return false;
    EqualPar[whichpar1][whichpar2] = true;
    EqualPar[whichpar2][whichpar1] = true;
    for(unsigned int i=0; i<totalpar; i++){
        //if whichpar1 is equal to a 3-rd parameter, and there is yet no relation between whichpar2 and
        //the 3-rd par in question, than whichpar2 is set equal to that parameter
        if(EqualPar[whichpar1][i] && i!=whichpar1 && i!=whichpar2 && !EqualPar[whichpar2][i]){
            SetEqualPar(whichpar2, i);
        }
        //the same is repeated for whichpar2
        if(EqualPar[whichpar2][i] && i!=whichpar1 && i!=whichpar2 && !EqualPar[whichpar1][i]){
            SetEqualPar(whichpar1, i);
        }
    }
    return true;
}




ROOT::Fit::FitResult DLM_MultiFit::PerformGlobalFit(bool printinfo){
    AllocMemFitter();
    CreateParMap();

    int DataSize = 0;

    for(unsigned short i=0; i<num_spectra; i++){
        WMTF1[i] = new ROOT::Math::WrappedMultiTF1(*FitFunction[i],1);
        double Min, Max; FitFunction[i]->GetRange(Min, Max);
        DR[i].SetRange(Min, Max);
        BD[i] = new ROOT::Fit::BinData(*opt,DR[i]);
        ROOT::Fit::FillData(*BD[i], HistoToFit[i]);
        Chi2Fun[i] = new ROOT::Fit::Chi2Function(*BD[i], *WMTF1[i]);
        DataSize += BD[i]->Size();
    }
    ROOT::Fit::Fitter fitter;

    double * parameter;
    parameter = new double [NumUniquePar];

    fitter.Config().SetParamsSettings(NumUniquePar,parameter);
    for(unsigned int i=0; i<NumUniquePar; i++){
        double parmin, parmax;
        for(unsigned int ispec=0; ispec<num_spectra; ispec++){
            for(int ipar=0; ipar<FitFunction[ispec]->GetNpar(); ipar++){
                if(ParMap[ispec][ipar]==i){
                    fitter.Config().ParSettings(i).SetValue(FitFunction[ispec]->GetParameter(ipar));
                    FitFunction[ispec]->GetParLimits(ipar, parmin, parmax);
                    if(parmin==parmax && parmin){//TF1 condition for fixed parameter
                        fitter.Config().ParSettings(i).Fix();
                    }
                    else if(parmin<parmax){//if the limits are reasonable
                        fitter.Config().ParSettings(i).SetLimits(parmin, parmax);
                    }
                    goto EndOfLoop;
                }
            }
        }
        EndOfLoop:;
    }

    fitter.Config().MinimizerOptions().SetPrintLevel(-1);
    fitter.Config().SetMinimizer("Minuit","Migrad");

    fitter.FitFCN(NumUniquePar,*this,0,DataSize,true);

    ROOT::Fit::FitResult myfitresult = fitter.Result();
    if(printinfo)
        myfitresult.Print(std::cout);

    bool CorrectCovMatrix = myfitresult.CovMatrixStatus()>2;
    if(!CorrectCovMatrix){
        //PrintError("Please note, that the covariance matrix is wrong, hense the fit might be incorrect!", "DLM_MultiFit::PerformGlobalFit");
    }

    //the parameters are set within the Spectrum1D
    for(unsigned int ispec=0; ispec<num_spectra; ispec++){
        for(int ipar=0; ipar<FitFunction[ispec]->GetNpar(); ipar++){
            FitFunction[ispec]->SetParameter(ipar, myfitresult.Parameter(ParMap[ispec][ipar]));
            FitFunction[ispec]->SetParError(ipar, myfitresult.ParError(ParMap[ispec][ipar]));
        }
        //one could here set the cov matrix as well
    }

    delete [] parameter;
    return myfitresult;
}

double DLM_MultiFit::operator() (const double *par){
    double result = 0;
    for(unsigned int ispec=0; ispec<num_spectra; ispec++){
        if(!spec_par[ispec]) spec_par[ispec] = new double [FitFunction[ispec]->GetNpar()];
        for(int ipar=0; ipar<FitFunction[ispec]->GetNpar(); ipar++){
            spec_par[ispec][ipar] = par[ParMap[ispec][ipar]];
        }
        result += (*Chi2Fun[ispec])(spec_par[ispec]);
    }
    return result;
}
