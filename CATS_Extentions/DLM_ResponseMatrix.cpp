#include "DLM_ResponseMatrix.h"


DLM_ResponseMatrix::DLM_ResponseMatrix(CATS& ab, TH2F* hs, TH2F* hr, const bool& ia):
    Kitty(ab),hSigmaMatrix(hs),hResidualMatrix(hr),InvertedAxis(ia),NumMomBins(Kitty.GetNumMomBins()){

    //in case we have a resolution matrix, we convert it to
    //the double** format and normalize it
    if(hSigmaMatrix){
        AllocateMatrix(mSigma);
        ConvertMatrix(mSigma, hSigmaMatrix, InvertedAxis);
    }
    //the same for the residual matrix
    if(hResidualMatrix){
        AllocateMatrix(mResidual);
        ConvertMatrix(mResidual, hResidualMatrix, InvertedAxis);
    }

    //in case we have no input matrices, we return the unit matrix as the response
    //The ResponseMatrix has its own memory allocation!
    if(!hSigmaMatrix && !hResidualMatrix){
        AllocateMatrix(mResponse);
        MakeUnitMatrix(mResponse);
        return;
    }
    //in case we have only the residual matrix we do not initialize the response matrix,
    //but simply point to the ResidualMatrix
    else if(!hSigmaMatrix){
        ResponseMatrix = ResidualMatrix;
        SparseResponse = SparseResidual;
    }
    //similar in the case of sigma matrix only
    else if(!ResidualMatrix){
        ResponseMatrix = SigmaMatrix;
        SparseResponse = SparseSigma;
    }
    //in case we have both matrices at hand, we multiply them and save the result to the ResponseMatrix
    //The ResponseMatrix has its own memory allocation!
    else{
        AllocateMatrix(mResponse);
        for(int iBin=0; iBin<NumMomBins; iBin++){
            SparseResponse[iBin][xAxisFirst] = -1;
            SparseResponse[iBin][xAxisLast] = -1;
            SparseResponse[iBin][yAxisFirst] = -1;
            SparseResponse[iBin][yAxisLast] = -1;
        }

        for(int iBinY=0; iBinY<NumMomBins; iBinY++){
            for(int iBinX=0; iBinX<NumMomBins; iBinX++){
                ResponseMatrix[iBinY][iBinX] = 0;
                for(int iBinK=0; iBinK<NumMomBins; iBinK++){
                    ResponseMatrix[iBinY][iBinX] += SigmaMatrix[iBinY][iBinK]*ResidualMatrix[iBinK][iBinX];
                }

                if(ResponseMatrix[iBinY][iBinX]!=0 && SparseResponse[iBinY][xAxisFirst]==-1){
                    SparseResponse[iBinY][xAxisFirst] = iBinX;
                }
                if(ResponseMatrix[iBinY][iBinX]!=0 && SparseResponse[iBinX][yAxisFirst]==-1){
                    SparseResponse[iBinX][yAxisFirst] = iBinY;
                }
                if(ResponseMatrix[iBinY][iBinX]!=0){
                    SparseResponse[iBinY][xAxisLast] = iBinX;
                }
                if(ResponseMatrix[iBinY][iBinX]!=0){
                    SparseResponse[iBinX][yAxisLast] = iBinY;
                }

            }
        }
    }

}

DLM_ResponseMatrix::~DLM_ResponseMatrix(){
    if(hSigmaMatrix){
        DeleteMatrix(mSigma);
    }
    if(hResidualMatrix){
        DeleteMatrix(mResidual);
    }
    //if both are defined or both are not defined
    if( !(bool(hSigmaMatrix)^bool(hResidualMatrix)) ){
        DeleteMatrix(mResponse);
    }

}

//allocates both the matrix and the sparse
void DLM_ResponseMatrix::AllocateMatrix(const int& WhichMatr){
    double*** Matrix = GetMatrixAddress(WhichMatr);
    int*** Sparse = GetSparseAddress(WhichMatr);
    Matrix[0] = new double* [NumMomBins];
    Sparse[0] = new int* [NumMomBins];
    for(int iBin=0; iBin<NumMomBins; iBin++){
        Matrix[0][iBin] = new double [NumMomBins];
        Sparse[0][iBin] = new int [4];
    }
}

//deletes both the matrix and the sparse
void DLM_ResponseMatrix::DeleteMatrix(const int& WhichMatr){
    double** Matrix = GetMatrix(WhichMatr);
    int** Sparse = GetSparse(WhichMatr);
    for(int iBin=0; iBin<NumMomBins; iBin++){
        delete [] Matrix[iBin];
        delete [] Sparse[iBin];
    }
    delete [] Matrix;
    delete [] Sparse;
}

double** DLM_ResponseMatrix::GetMatrix(const int& WhichMatr){
    switch(WhichMatr){
    case mSigma :
        return SigmaMatrix;
    case mResidual :
        return ResidualMatrix;
    case mResponse :
        return ResponseMatrix;
    default : printf("ERROR: DLM_ResponseMatrix::GetSparse says bad luck!\n"); return NULL;
    }
}

double*** DLM_ResponseMatrix::GetMatrixAddress(const int& WhichMatr){
    switch(WhichMatr){
    case mSigma :
        return &SigmaMatrix;
    case mResidual :
        return &ResidualMatrix;
    case mResponse :
        return &ResponseMatrix;
    default : printf("ERROR: DLM_ResponseMatrix::GetSparse says bad luck!\n"); return NULL;
    }
}

int** DLM_ResponseMatrix::GetSparse(const int& WhichMatr){
    switch(WhichMatr){
    case mSigma :
        return SparseSigma;
    case mResidual :
        return SparseResidual;
    case mResponse :
        return SparseResponse;
    default : printf("ERROR: DLM_ResponseMatrix::GetSparse says bad luck!\n"); return NULL;
    }
}

int*** DLM_ResponseMatrix::GetSparseAddress(const int& WhichMatr){
    switch(WhichMatr){
    case mSigma :
        return &SparseSigma;
    case mResidual :
        return &SparseResidual;
    case mResponse :
        return &SparseResponse;
    default : printf("ERROR: DLM_ResponseMatrix::GetSparse says bad luck!\n"); return NULL;
    }
}

void DLM_ResponseMatrix::MakeUnitMatrix(const int& WhichMatr){
    int** Sparse = GetSparse(WhichMatr);
    double** Matrix = GetMatrix(WhichMatr);
    for(int iBinX=0; iBinX<NumMomBins; iBinX++){
        for(int iBinY=0; iBinY<NumMomBins; iBinY++){
            Matrix[iBinY][iBinX] = (iBinX==iBinY);
        }
        Sparse[iBinX][xAxisFirst] = iBinX;
        Sparse[iBinX][xAxisLast] = iBinX;
        Sparse[iBinX][yAxisFirst] = iBinX;
        Sparse[iBinX][yAxisLast] = iBinX;
    }
}

void DLM_ResponseMatrix::ConvertMatrix(const int& WhichMatr, TH2F* input, const bool& InvAxis){

    int** Sparse = GetSparse(WhichMatr);
    double** Matrix = GetMatrix(WhichMatr);


    double MomLowEdgeX;
    double MomUpEdgeX;
    double MomLowEdgeY;
    double MomUpEdgeY;
    int WhichBinAtLowEdgeX;
    int WhichBinAtUpEdgeX;
    int WhichBinAtLowEdgeY;
    int WhichBinAtUpEdgeY;
    int NumOldBinsX;
    int NumOldBinsY;
    const int MaxBufferSize = 256;
    double* WeightX = new double [MaxBufferSize];
    double* WeightY = new double [MaxBufferSize];
    double* WeightXY = new double [MaxBufferSize];

    //if the axis are inverted
    TAxis* Xaxis = InvAxis?input->GetYaxis():input->GetXaxis();
    TAxis* Yaxis = InvAxis?input->GetXaxis():input->GetYaxis();
    for(int iBin=0; iBin<NumMomBins; iBin++){
        Sparse[iBin][xAxisFirst] = -1;
        Sparse[iBin][xAxisLast] = -2;
        Sparse[iBin][yAxisFirst] = -1;
        Sparse[iBin][yAxisLast] = -2;
    }

    for(int iBinX=0; iBinX<NumMomBins; iBinX++){
        MomLowEdgeX = Kitty.GetMomBinLowEdge(iBinX);
        MomUpEdgeX = Kitty.GetMomBinUpEdge(iBinX);
        WhichBinAtLowEdgeX = Xaxis->FindBin(MomLowEdgeX);
        WhichBinAtUpEdgeX = Xaxis->FindBin(MomUpEdgeX);
        NumOldBinsX = WhichBinAtUpEdgeX - WhichBinAtLowEdgeX + 1;
//Printf("MomLowEdgeX=%f", MomLowEdgeX);
//Printf("MomUpEdgeX=%f", MomUpEdgeX);
//Printf("WhichBinAtLowEdgeX=%i", WhichBinAtLowEdgeX);
//Printf("WhichBinAtUpEdgeX=%i", WhichBinAtUpEdgeX);
//Printf("Xaxis->GetBinLowEdge(%i)=%f", iBinX, Xaxis->GetBinLowEdge(iBinX));
//usleep(2000e3);
        if(NumOldBinsX>MaxBufferSize){
            printf("ERROR: DLM_ResponseMatrix::ConvertMatrix hates you since NumOldBinsX>MaxBufferSize\nA CRASH SHOULD FOLLOW!\n");
        }
        for(int iBinY=0; iBinY<NumMomBins; iBinY++){
            MomLowEdgeY = Kitty.GetMomBinLowEdge(iBinY);
            MomUpEdgeY = Kitty.GetMomBinUpEdge(iBinY);
            WhichBinAtLowEdgeY = Xaxis->FindBin(MomLowEdgeY);
            WhichBinAtUpEdgeY = Xaxis->FindBin(MomUpEdgeY);
            NumOldBinsY = WhichBinAtUpEdgeY - WhichBinAtLowEdgeY + 1;
            if(NumOldBinsY>MaxBufferSize){
                printf("ERROR: DLM_ResponseMatrix::ConvertMatrix hates you since NumOldBinsX>MaxBufferSize\nA CRASH SHOULD FOLLOW!\n");
            }

            for(int iOldBinX=1; iOldBinX<NumOldBinsX-1; iOldBinX++){
                WeightX[iOldBinX] = 1;
            }
            WeightX[0] = (Xaxis->GetBinUpEdge(WhichBinAtLowEdgeX)-MomLowEdgeX)/Xaxis->GetBinWidth(WhichBinAtLowEdgeX);
            if(NumOldBinsX>1){
                WeightX[NumOldBinsX-1] = (MomUpEdgeX-Xaxis->GetBinLowEdge(WhichBinAtUpEdgeX))/Xaxis->GetBinWidth(WhichBinAtUpEdgeX);
            }

            for(int iOldBinY=1; iOldBinY<NumOldBinsY-1; iOldBinY++){
                WeightY[iOldBinY] = 1;
            }
            WeightY[0] = (Yaxis->GetBinUpEdge(WhichBinAtLowEdgeY)-MomLowEdgeY)/Yaxis->GetBinWidth(WhichBinAtLowEdgeY);
            if(NumOldBinsY>1){
                WeightY[NumOldBinsY-1] = (MomUpEdgeY-Yaxis->GetBinLowEdge(WhichBinAtUpEdgeY))/Yaxis->GetBinWidth(WhichBinAtUpEdgeY);
            }

            Matrix[iBinY][iBinX] = 0;
            for(int iOldBinX=0; iOldBinX<NumOldBinsX; iOldBinX++){
                for(int iOldBinY=0; iOldBinY<NumOldBinsY; iOldBinY++){
                    Matrix[iBinY][iBinX] += WeightX[iOldBinX]*WeightY[iOldBinY]*
                        input->GetBinContent(WhichBinAtLowEdgeX+iOldBinX, WhichBinAtLowEdgeY+iOldBinY);
                }
            }

            if(Matrix[iBinY][iBinX]!=0 && Sparse[iBinY][xAxisFirst]==-1){
                Sparse[iBinY][xAxisFirst] = iBinX;
            }
            if(Matrix[iBinY][iBinX]!=0 && Sparse[iBinX][yAxisFirst]==-1){
                Sparse[iBinX][yAxisFirst] = iBinY;
            }
            if(Matrix[iBinY][iBinX]!=0){
                Sparse[iBinY][xAxisLast] = iBinX;
            }
            if(Matrix[iBinY][iBinX]!=0){
                Sparse[iBinX][yAxisLast] = iBinY;
            }

        }
    }
    NormalizeMatrix(WhichMatr);

    delete [] WeightX;
    delete [] WeightY;
    delete [] WeightXY;
}

double DLM_ResponseMatrix::BilinearInterpolation(const double& x0, const double& y0,
                                 const double& x1, const double& y1,
                                 const double& f00, const double& f01, const double& f10, const double& f11,
                                 const double& xVal, const double& yVal){
    return (f00*(x1-xVal)*(y1-yVal)+f10*(xVal-x0)*(y1-yVal)+f01*(x1-xVal)*(yVal-y0)+f11*(xVal-x0)*(yVal-y0))/(x1-x0)/(y1-y0);
}


void DLM_ResponseMatrix::NormalizeMatrix(const int& WhichMatr){

    int** Sparse = GetSparse(WhichMatr);
    double** Matrix = GetMatrix(WhichMatr);

    double* Norm = new double [NumMomBins];
    for(int iBinX=0; iBinX<NumMomBins; iBinX++) Norm[iBinX] = 0;

    for(int iBinY=0; iBinY<NumMomBins; iBinY++){
        for(int iBinX=Sparse[iBinY][xAxisFirst]; iBinX<=Sparse[iBinY][xAxisLast]; iBinX++){
            Norm[iBinY] += Matrix[iBinY][iBinX];
        }
    }

    for(int iBinY=0; iBinY<NumMomBins; iBinY++){
        for(int iBinX=Sparse[iBinY][xAxisFirst]; iBinX<=Sparse[iBinY][xAxisLast]; iBinX++){
            if(Norm[iBinY]) Matrix[iBinY][iBinX] /= Norm[iBinY];
            else Matrix[iBinY][iBinX]=0;
        }
    }

    delete [] Norm;

}




void RespMatrTest1(){

    const int NumBins = 3;

    TH2F* hTest1 = new TH2F("hTest1", "hTest1", NumBins, 0, 60, NumBins, 0, 60);
    for(int iBinX=1; iBinX<=NumBins; iBinX++){
        for(int iBinY=1; iBinY<=NumBins; iBinY++){
            hTest1->SetBinContent(iBinX, iBinY, 0);
        }
    }

    CATS Afterburner;
    Afterburner.SetMomBins(3, 0, 60);

    const double Norm = 1;

    hTest1->SetBinContent(1,1,0.8*Norm);
    hTest1->SetBinContent(1,2,0.2*Norm);

    hTest1->SetBinContent(2,1,0.1*Norm);
    hTest1->SetBinContent(2,2,0.8*Norm);
    hTest1->SetBinContent(2,3,0.1*Norm);

    hTest1->SetBinContent(3,2,0.125*Norm);
    hTest1->SetBinContent(3,3,0.75*Norm);


    TH2F* hTest2 = new TH2F("hTest2", "hTest2", NumBins, 0, 60, NumBins, 0, 60);
    hTest2->SetBinContent(1,1,0.9*Norm);
    hTest2->SetBinContent(1,2,0.1*Norm);

    hTest2->SetBinContent(2,1,0.05*Norm);
    hTest2->SetBinContent(2,2,0.9*Norm);
    hTest2->SetBinContent(2,3,0.05*Norm);

    hTest2->SetBinContent(3,2,0.05*Norm);
    hTest2->SetBinContent(3,3,0.95*Norm);

    DLM_ResponseMatrix dlmRespMatr(Afterburner, hTest1, hTest2, false);

    delete hTest1;

}


