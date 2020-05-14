
#include <stdio.h>
#include <string.h>

#include "DLM_CkDecomposition.h"
#include "DLM_RootWrapper.h"

#include "TH2F.h"

DLM_CkDecomposition::DLM_CkDecomposition(const char* name, const unsigned& numchildren, DLM_Ck& ckfunction, const TH2F* hSigmaMatrix, const bool& InvertAxis):
    DLM_CkDecomp(name, numchildren, ckfunction){

    dlmSigmaMatrix = NULL;
    dlmFeedMatrix = NULL;
    if(hSigmaMatrix){
        dlmSigmaMatrix = Convert_TH2F_DlmHisto(hSigmaMatrix);
        AddTheMatrix(dlmSigmaMatrix,InvertAxis);
    }
    Update(true);
}

DLM_CkDecomposition::~DLM_CkDecomposition(){
    if(dlmSigmaMatrix) {delete dlmSigmaMatrix; dlmSigmaMatrix=NULL;}
    if(dlmFeedMatrix){
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            if(dlmFeedMatrix[uChild]){
                delete dlmFeedMatrix[uChild];
                dlmFeedMatrix[uChild] = NULL;
            }
        }
        delete [] dlmFeedMatrix;
        dlmFeedMatrix = NULL;
    }
}

void DLM_CkDecomposition::AddContribution(const unsigned& WhichCk, const double& fraction, const int& type, DLM_CkDecomposition* child,
                         const TH2F* hResidualMatrix, const bool& InvertedAxis){
    if(WhichCk>=NumChildren){
        return;
    }
    if(!dlmFeedMatrix){
        dlmFeedMatrix = new DLM_Histo<float>* [NumChildren];
        for(unsigned uChild=0; uChild<NumChildren; uChild++){
            dlmFeedMatrix[uChild] = NULL;
        }
    }
    if(dlmFeedMatrix[WhichCk]){
        delete dlmFeedMatrix[WhichCk];
        dlmFeedMatrix[WhichCk] = NULL;
    }
    if(hResidualMatrix) dlmFeedMatrix[WhichCk] = Convert_TH2F_DlmHisto(hResidualMatrix);
    DLM_CkDecomp::AddContribution(WhichCk,fraction,type,child,dlmFeedMatrix[WhichCk],InvertedAxis);
}

