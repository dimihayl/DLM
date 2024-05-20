#ifndef DLM_CKDECOMPOSITION_H
#define DLM_CKDECOMPOSITION_H

#include "DLM_CkDecomp.h"

class TH2F;
class TH1F;
class DLM_Ck;
class DLM_CkDecomp;
template <class Type> class DLM_Histo;

class DLM_CkDecomposition : public DLM_CkDecomp{

public:
    enum CkDecType { cFeedDown, cFake };

    DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const TH2F* hSigmaMatrix, const bool& InvertAxis=false);
    ~DLM_CkDecomposition();

    void AddContribution(unsigned WhichCk, DLM_Histo<double>& fraction, int type, DLM_CkDecomposition* child=NULL,
                         const TH2F* hResidualMatrix=NULL, const bool& InvertedAxis=false);
    void AddContribution(unsigned WhichCk, double fraction, int type, DLM_CkDecomposition* child=NULL,
                         const TH2F* hResidualMatrix=NULL, const bool& InvertedAxis=false);
    void AddPhaseSpace(const unsigned& WhichCk, const TH1F* hPhaseSpace);
    void AddPhaseSpace(const TH1F* hPhaseSpace);
    //const TH2F* GetResolutionMatrix();

private:
    DLM_Histo<float>* dlmSigmaMatrix;
    DLM_Histo<float>** dlmFeedMatrix;
    DLM_Histo<float>** dlmPhaseSpace;
    DLM_Histo<float>* dlmPhaseSpaceMain;
};

#endif
