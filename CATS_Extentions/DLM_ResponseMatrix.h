#ifndef DLM_ResponseMatrixH
#define DLM_ResponseMatrixH

#include "CATS.h"
#include "DLM_Histo.h"

#include "TH2F.h"

class DLM_ResponseMatrix{

public:
    DLM_ResponseMatrix(DLM_Histo<double>& ch, const TH2F* hs, const TH2F* hr, const bool& ia=false);
    DLM_ResponseMatrix(CATS& ab, const TH2F* hs, const TH2F* hr, const bool& ia=false);
    ~DLM_ResponseMatrix();

    DLM_Histo<double>* CatHisto;
    bool CatHistoIsMyOwn;
    const TH2F* hSigmaMatrix;
    const TH2F* hResidualMatrix;

    //by default the input matrices are supposed to be [Y][X]
    //with the X axis being the the original (unsmeared) momentum
    const bool& InvertedAxis;
    const int NumMomBins;

//! The residual matrix is actually the transverse. The reason is that this makes a better
//memory ordering for matrix multiplication
    //[Y][X]
    double** SigmaMatrix;
    //the sparse info says which is the first and the last non-zero
    //entry along the X/Y axis
    int** SparseSigma;
    //[X][Y]
    double** ResidualMatrix;
    int** SparseResidual;

    //[Y][X]
    double** ResponseMatrix;
    int** SparseResponse;

private:

    void DefaultConstructor();

    void AllocateMatrix(const int& WhichMatr);
    void DeleteMatrix(const int& WhichMatr);

    double** GetMatrix(const int& WhichMatr);
    double*** GetMatrixAddress(const int& WhichMatr);
    int** GetSparse(const int& WhichMatr);
    int*** GetSparseAddress(const int& WhichMatr);
    void MakeUnitMatrix(const int& WhichMatr);

//converts the TH2F matrix to double** format
    void ConvertMatrix(const int& WhichMatr, const TH2F* input, const bool& InvAxis);

    double BilinearInterpolation(const double& x0, const double& y0,
                                 const double& x1, const double& y1,
                                 const double& f00, const double& f01, const double& f10, const double& f11,
                                 const double& xVal, const double& yVal);

    void NormalizeMatrix(const int& WhichMatr);

    enum enumSparse { xAxisFirst, xAxisLast, yAxisFirst, yAxisLast };
    enum WhichMatrix { mSigma, mResidual, mResponse };
};


void RespMatrTest1();

#endif

