#ifndef DLM_SmearedCatsH
#define DLM_SmearedCatsH

class CATS;
class TH2F;
class DLM_ResponseMatrix;



//! NOTE THAT THIS THING SEEMS TO PERFORM VERY BADLY IN CASE THE RESPONSE MATRICES ARE
//NOT WITH THE SAME BINNING AS THE DATA, SO TRY TO ALWAYS USE THE SAME BINNING
class DLM_SmearedCats{

public:
    DLM_SmearedCats(CATS** InCat, const unsigned& numCk);
//    DLM_SmearedCats(const char* InputFile);
    ~DLM_SmearedCats();

    void SetResolutionMatrix(TH2F* resolution);
    void SetResidualMatrix(const unsigned& WhichNr, TH2F* residual);
    void SetLambda(const unsigned& WhichNr, const double& lam);
    void SetUseLednicky(const unsigned& WhichNr, const int& val,
                                     const double& grad,
                                     const double& ScattLen1, const double& EffRan1,
                                     const double& ScattLen3, const double& EffRan3,
                                     const double& ares, const double& arad, const double& lambda0);
    void UseCustomModel(const unsigned& WhichNr, double (*fModel)(const double&, const double*), double* Pars);
    void UseCats(const unsigned& WhichNr);
    double GetCorrectedCk(const unsigned& WhichBin);
    double GetCorrectedCkErr(const unsigned& WhichBin);
    double EvalCorrectedCk(const double& Momentum);
    double EvalCorrectedCkErr(const double& Momentum);
    CATS* GetTheKitty(const unsigned& WhichOne);

    void Correct(const bool& NewBinning);

    double CkLednicky(const double& Momentum, const bool& SinOnly, const bool& QS, const bool& WithLambda=false);

    //without extension, automatically the extension is set to .dlmsc
//    void SaveAs(const char* FileName);

    bool GetStatus();

private:
    double* LambdaCoeff;
    double* CorrectedCk;
    double* CorrectedCkErr;
    TH2F* hResolution;
    TH2F** hResidual;
    DLM_ResponseMatrix** RespMatrix;
    CATS** cat;
    const unsigned NumCk;

    //[WhichCk], the values:
    //0=use cats
    //1=Lednicky(sin+tri) without QS
    //2=Lednicky(sin+tri) with QS -> may be wrong, check
    //3=Lednicky(sin only) without QS -> may be wrong, check
    //4=Lednicky(sin only) with QS
    //-1=Custom Model
    int* EvalUsingEquation;

    double ScattLenSin;
    double EffRangeSin;
    double ScattLenTri;
    double EffRangeTri;
    double GaussR;
    double ResidualR;
    double aResidual;
    double LambdaLedni;

    bool* LedniStatus;

    //Momentum, Parameters
    double (*ModelFunction)(const double&, const double*);
    double* ModelParameter;
};



#endif
