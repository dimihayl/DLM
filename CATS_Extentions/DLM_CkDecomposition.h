#ifndef DLM_CKDECOMPOSITION_H
#define DLM_CKDECOMPOSITION_H

#include "CATS.h"
#include "CATStools.h"
#include "DLM_ResponseMatrix.h"
#include "DLM_Ck.h"
#include "DLM_CkDecomp.h"
#include "DLM_Histo.h"

class TH2F;
class TH1F;
//class DLM_Histo<float>;

class DLM_CkDecomposition : public DLM_CkDecomp{

public:
    enum CkDecType { cFeedDown, cFake };

    DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const TH2F* hSigmaMatrix, const bool& InvertAxis=false);
    ~DLM_CkDecomposition();

    void AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomposition* child=NULL,
                         const TH2F* hResidualMatrix=NULL, const bool& InvertedAxis=false);
    void AddPhaseSpace(const unsigned& WhichCk, const TH1F* hPhaseSpace);
    //const TH2F* GetResolutionMatrix();

private:
    DLM_Histo<float>* dlmSigmaMatrix;
    DLM_Histo<float>** dlmFeedMatrix;
    DLM_Histo<float>** dlmPhaseSpace;
};

#endif
