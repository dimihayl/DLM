#ifndef DLM_MILTIFIT_H
#define DLM_MILTIFIT_H

#ifndef ROOT_Math_WrappedMultiTF1
#include "Math/WrappedMultiTF1.h"
#endif

#ifndef ROOT_Fit_DataOptions
#include "Fit/DataOptions.h"
#endif

#ifndef ROOT_Fit_DataRange
#include "Fit/DataRange.h"
#endif

#ifndef ROOT_Fit_BinData
#include "Fit/BinData.h"
#endif

#ifndef ROOT_Fit_Chi2FCN
#include "Fit/Chi2FCN.h"
#endif

class TH1;
class TF1;

class DLM_MultiFit{

protected:
//public:
    //!NUMBERING OF THE PARAMETERS!
    //there are several numbering system, and it is important to know when which is used:
    //!Spectrum numbering system (SNS):
    //[i][j] -> the j-th par of the i-th spectrum
    //!Global numbering system (GNS):
    //[i] -> the i-the global par., where all parameters are listed one ofter the other,
    //the firts few parameters are from spectrum1, than from spectrum2 etc.
    //!Unique numbering system (UNS):
    //[i] -> since some of the par. should be equal to one another, the UNS transforms the GNS in
    //a way, that there are only independent parameters. ParMap gives the conversion from SNS to UNS

    unsigned short num_spectra;
    unsigned short max_spectra;
    unsigned int totalpar;
    //DLM_Spectrum1D ** Spectra;
    TH1** HistoToFit;
    TF1** FitFunction;
    double** spec_par;

    //[i][j] : tells you if the i-th par is equal to
    //j-th par. The numbering follows the GNS
    bool ** EqualPar;
    //the number of independent parameters (i.e. in UNS)
    unsigned int NumUniquePar;
    //[i][j] = n : the j-th par of the i-th spectrum corresponds to the n-th independent parameter
    //i.e. transformation from SNS to UNS
    unsigned int ** ParMap;
    //finds NumUniquePar and creates the ParMap (based on EqualPar)
    void CreateParMap();

    ROOT::Math::WrappedMultiTF1 ** WMTF1;
    ROOT::Fit::DataOptions* opt;
    ROOT::Fit::DataRange * DR;
    ROOT::Fit::BinData ** BD;
    ROOT::Fit::Chi2Function ** Chi2Fun;

    void ReallocMemSpectra(unsigned short m);
    void AllocMemPar();

    void AllocMemFitter();
    void ClearMemFitter();


public:

    DLM_MultiFit(); //done
    ~DLM_MultiFit(); //done

    bool AddSpectrum(TH1* histo, TF1* fit); //done

    void ClearSpectra();
    void ResetParRelations();

    //!IMPORTANT! The code is optimized for setting the paramters AFTER all spectra are added.
    //This means, that if any spectra are added later on, the whole information about relations between
    //parameters will be lost!
    bool SetEqualPar(unsigned short spec1, int par1, unsigned short spec2, int par2); //in UNS
    bool SetEqualPar(unsigned int whichpar1, unsigned int whichpar2); //in GNS

    unsigned int GetNumTotalPar() {return totalpar;}
    unsigned int GetNumUniquePar() {return NumUniquePar;}

    //one potential problem with this whole procedure are the initial values, limits, fixations etc. of the
    //parameters. At the moment the whole procedure
    //!applies the limits and values of each parameter according to the first spectrum that contains
    //!the parameter in question (i.e. the first one that has been added using AddSpectrum)!
    ROOT::Fit::FitResult PerformGlobalFit(bool printinfo=true);

    double operator() (const double *par);

};

#endif
