
#ifndef DLM_UNFOLD_H
#define DLM_UNFOLD_H

#include "DLM_Histo.h"
#include "TString.h"

class TH1F;
class TH2F;
class TF1;

class DLM_Unfold{
public:
    DLM_Unfold();
    ~DLM_Unfold();
    //supposed to be in counts
    void SetData(TH1F* data);
    //supposed to be a function, i.e. normalized to the bin width
    void SetResponse(TH2F* response);
    void SetUnfoldPrecision(const double& precision, const double& worst);
    void SetUnfoldBootstrap(const unsigned& numiter, const int& seedmin);
    void SetUnfoldSeconds(const double& seconds);
    void SetUnfoldMinutes(const double& minutes);
    void SetUnfoldRange(const double& min, const double& max);
    void SetSilentMode(const bool& silentmode);
    void SetPcOutput(const bool& pcoutput) {PcOutput=pcoutput;}
    //void SetOutput(const int& level, const TString& ouputfilename);
    DLM_Histo<float>* Unfold();
    //DLM_Histo<float>* UnfoldNew();
    DLM_Histo<float>* Fold();
    double GetReachedPrecision() const {return ReachedPrecision;}
    TF1* GetFitFoldOriginal() const {return FitFoldOriginal;}
    TF1* GetFitFoldFinal() const {return FitFoldFinal;}
    TF1* GetFitUnfoldFinal() const {return FitUnfoldFinal;}
    std::vector<float>* GetUnfoldedResult(){return BinByBinUnfolded;}
    std::vector<float>* GetFoldedResult(){return BinByBinFolded;}
private:
    void Fold(const DLM_Histo<float>& DATA, DLM_Histo<float>& RESULT);
    TH1F* hData;
    TH2F* hResponse;
    DLM_Histo<float>* dlmData;
    DLM_Histo<float>* dlmResponse;
    unsigned NumIter;
    double Precision;
    double PrecisionWorst;
    double ReachedPrecision;
    double MaxTime;
    double RangeMin;
    double RangeMax;
    int SEEDmin;
    bool Silent;
    bool PcOutput;
    double SplineFit(double* xVal, double* pars);
    double SplineFitFolded(double* xVal, double* pars);
    double SplineFitUnfolded(double* xVal, double* pars);
    double DimiFit(double* xVal, double* pars);
    double DimiFitFolded(double* xVal, double* pars);
    DLM_Histo<float>* dlmData_reb2;
    double* DLM_FITTER_ARRAY_SPLINE3_X;
    double* DLM_FITTER_ARRAY_SPLINE3_Y;
    unsigned* SparseFirst;
    unsigned* SparseLast;
    int OutputLevel;
    TString OutputFileName;
    TF1* FitFoldOriginal;
    TF1* FitFoldFinal;
    TF1* FitUnfoldFinal;
    //TH1F* hFit_folded;
    //TH1F* hData_Unfolded;
    DLM_Histo<float>* dlmFoldedWorkHorse;
    DLM_Histo<float>* dlmUnfoldedWorkHorse;

    DLM_Histo<float>* dlmFoldedDimiWorkHorse;
    DLM_Histo<float>* dlmUnfoldedDimiWorkHorse;
    DLM_Histo<float>* dlmDimiWorkHorse;

    std::vector<float>* BinByBinFolded;
    std::vector<float>* BinByBinUnfolded;
};

#endif
