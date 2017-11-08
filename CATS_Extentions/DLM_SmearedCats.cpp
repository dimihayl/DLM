#include "DLM_SmearedCats.h"

using namespace std;

DLM_SmearedCats::DLM_SmearedCats(CATS** InCat, const unsigned& numCk):
        cat(InCat),NumCk(numCk){
    LambdaCoeff = new double [NumCk];
    hResolution = NULL;
    hResidual = new TH2F* [NumCk];
    RespMatrix = new DLM_ResponseMatrix* [NumCk];
    EvalUsingEquation = new int [NumCk];
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        hResidual[uCk] = NULL;
        RespMatrix[uCk] = NULL;
        EvalUsingEquation[uCk] = 0;
    }
    CorrectedCk = NULL;
    CorrectedCkErr = NULL;
}

DLM_SmearedCats::~DLM_SmearedCats(){
    delete [] LambdaCoeff;
    delete [] hResidual;
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        if(RespMatrix[uCk]){delete RespMatrix[uCk]; RespMatrix[uCk]=NULL;}
    }
    delete [] RespMatrix;
    if(CorrectedCk) {delete [] CorrectedCk; CorrectedCk=NULL;}
    if(CorrectedCkErr) {delete [] CorrectedCkErr; CorrectedCkErr=NULL;}
}

void DLM_SmearedCats::SetResolutionMatrix(TH2F* resolution){
    hResolution = resolution;
}

void DLM_SmearedCats::SetResidualMatrix(const unsigned& WhichNr, TH2F* residual){
    hResidual[WhichNr] = residual;
}

void DLM_SmearedCats::SetLambda(const unsigned& WhichNr, const double& lam){
    LambdaCoeff[WhichNr] = lam;
}

void DLM_SmearedCats::SetUseLednicky(const unsigned& WhichNr, const int& val, const double& grad,
                                     const double& ScattLen1, const double& EffRan1,
                                     const double& ScattLen3, const double& EffRan3,
                                     const double& ares, const double& arad, const double& lambda0){
    EvalUsingEquation[WhichNr] = val;
    GaussR = grad;
    ScattLenSin = ScattLen1;
    EffRangeSin = EffRan1;
    ScattLenTri = ScattLen3;
    EffRangeTri = EffRan3;
    aResidual = ares;
    ResidualR = arad;
    LambdaLedni = lambda0;
}

void DLM_SmearedCats::Correct(const bool& NewBinning){
    if(!cat[0]){
        printf("ERROR: DLM_SmearedCats::Correct has bad input!\n");
        return;
    }
    for(unsigned uCk=0; uCk<NumCk; uCk++){
        if(NewBinning && RespMatrix[uCk]){
            delete RespMatrix[uCk];
            RespMatrix[uCk] = NULL;
        }
        if(!RespMatrix[uCk] && (hResolution||hResidual) && cat[uCk]){
            RespMatrix[uCk] = new DLM_ResponseMatrix(cat[0][0], hResolution, hResidual[uCk]);
        }
    }
    double MomentumTrue;
    //unsigned WhichMomBin;

    if(CorrectedCk && NewBinning){
        delete [] CorrectedCk; CorrectedCk=NULL;
    }
    if(CorrectedCkErr && NewBinning){
        delete [] CorrectedCkErr; CorrectedCkErr=NULL;
    }
    if(!CorrectedCk){
        CorrectedCk = new double [cat[0]->GetNumMomBins()];
    }
    if(!CorrectedCkErr){
        CorrectedCkErr = new double [cat[0]->GetNumMomBins()];
    }

    double CkVal;
    double CkValErr;

    for(unsigned uBinSmear=0; uBinSmear<cat[0]->GetNumMomBins(); uBinSmear++){
        CorrectedCk[uBinSmear] = 0;
        CorrectedCkErr[uBinSmear] = 0;

        for(unsigned uBinTrue=0; uBinTrue<cat[0]->GetNumMomBins(); uBinTrue++){
            MomentumTrue = cat[0]->GetMomentum(uBinTrue);
            for(unsigned uCk=0; uCk<NumCk; uCk++){
                if(RespMatrix[uCk] && cat[uCk]){
//printf("EvalUsingEquation[uCk]=%i\n", EvalUsingEquation[uCk]);
                    switch(EvalUsingEquation[uCk]){
                        case 0 :    CkVal = cat[uCk]->EvalCorrFun(MomentumTrue);
                                    CkValErr = cat[uCk]->EvalCorrFunErr(MomentumTrue);
                                    break;
                        case 1 :    CkVal = CkLednicky(MomentumTrue, false, false);
                                    CkValErr = 0;
                                    break;
                        case 2 :    CkVal = CkLednicky(MomentumTrue, false, true);
                                    CkValErr = 0;
                                    break;
                        case 3 :    CkVal = CkLednicky(MomentumTrue, true, false);
                                    CkValErr = 0;
                                    break;
                        case 4 :    CkVal = CkLednicky(MomentumTrue, true, true);
                                    CkValErr = 0;
                                    break;
                        default:    CkVal = 0;
                                    CkValErr = 0;
                                    break;
                    }

                    CorrectedCk[uBinSmear] +=   LambdaCoeff[uCk]*
                                    RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]*CkVal;
                    CorrectedCkErr[uBinSmear] +=   LambdaCoeff[uCk]*
                                    RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]*CkValErr;
if(CorrectedCkErr[uBinSmear]!=CorrectedCkErr[uBinSmear]){
    if(CkValErr!=CkValErr) printf("TROUBLE CkValErr -> k=%f; cat[uCk]->GetCorrFunErr(0)=%f!!!\n",MomentumTrue,cat[uCk]->GetCorrFunErr(0));
    if(LambdaCoeff[uCk]!=LambdaCoeff[uCk]) printf("TROUBLE lam!!!\n");
    if(RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]!=RespMatrix[uCk]->ResponseMatrix[uBinSmear][uBinTrue]) printf("TROUBLE RM!!!\n");
    //printf();
}

                }
            }
        }

        for(unsigned uCk=0; uCk<NumCk; uCk++){
            if(!RespMatrix[uCk] || !cat[uCk]){
                CorrectedCk[uBinSmear] += LambdaCoeff[uCk];
            }
        }
    }
}


double DLM_SmearedCats::GetCorrectedCk(const unsigned& WhichBin){
    return CorrectedCk[WhichBin];
}

double DLM_SmearedCats::GetCorrectedCkErr(const unsigned& WhichBin){
    return CorrectedCkErr[WhichBin];
}

double DLM_SmearedCats::EvalCorrectedCk(const double& Momentum){
    unsigned NumMomBins = cat[0]->GetNumMomBins();
    if(Momentum<cat[0]->GetMomBinLowEdge(0) || Momentum>cat[0]->GetMomBinLowEdge(NumMomBins)) return 0;
    if(NumMomBins==1) return CorrectedCk[0];
    unsigned WhichMomBin = cat[0]->GetMomBin(Momentum);

    double RelMom[3];
    RelMom[0] = WhichMomBin?cat[0]->GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = cat[0]->GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin<(NumMomBins-1)?cat[0]->GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double* CkRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCk[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCk[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCk[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCk[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrectedCk[WhichMomBin];
    }

    return (CkRange[1]*(Momentum-InterpolRange[0])-
            CkRange[0]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

double DLM_SmearedCats::EvalCorrectedCkErr(const double& Momentum){
    unsigned NumMomBins = cat[0]->GetNumMomBins();
    if(Momentum<cat[0]->GetMomBinLowEdge(0) || Momentum>cat[0]->GetMomBinLowEdge(NumMomBins)) return 0;
    if(NumMomBins==1) return CorrectedCkErr[0];
    unsigned WhichMomBin = cat[0]->GetMomBin(Momentum);

    double RelMom[3];
    RelMom[0] = WhichMomBin?cat[0]->GetMomentum(WhichMomBin-1):-1;
    RelMom[1] = cat[0]->GetMomentum(WhichMomBin);
    RelMom[2] = WhichMomBin<(NumMomBins-1)?cat[0]->GetMomentum(WhichMomBin+1):-1;

    double* InterpolRange;
    double* CkRange;

    if(RelMom[0]==-1){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCkErr[WhichMomBin];
    }
    else if(RelMom[2]==-1){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCkErr[WhichMomBin-1];
    }
    else if(Momentum<RelMom[1]){
        InterpolRange = &RelMom[0];
        CkRange = &CorrectedCkErr[WhichMomBin-1];
    }
    else if(RelMom[1]<Momentum){
        InterpolRange = &RelMom[1];
        CkRange = &CorrectedCkErr[WhichMomBin];
    }
    else{//RelMom[1]==Momentum
        return CorrectedCkErr[WhichMomBin];
    }

    return (CkRange[1]*(Momentum-InterpolRange[0])-
            CkRange[0]*(Momentum-InterpolRange[1]))/
            (InterpolRange[1]-InterpolRange[0]);
}

CATS* DLM_SmearedCats::GetTheKitty(const unsigned& WhichOne){
    if(WhichOne>=NumCk) return NULL;
    return cat[WhichOne];
}

double DLM_SmearedCats::CkLednicky(const double& Momentum, const bool& SinOnly, const bool& QS, const bool& WithLambda){
    const std::complex<double> i(0,1);
    const double Pi(3.141592653589793);

    double F1 = gsl_sf_dawson(2.*Momentum*GaussR)/(2.*Momentum*GaussR);
    double F2 = (1.-exp(-4.*Momentum*Momentum*GaussR*GaussR))/(2.*Momentum*GaussR);

    complex<double> ScattAmplSin = pow(1./ScattLenSin+0.5*EffRangeSin*Momentum*Momentum-i*Momentum,-1.);

    double CkValue = 0.;
    CkValue +=  0.5*pow(abs(ScattAmplSin)/GaussR,2)*(1-(EffRangeSin)/(2*sqrt(Pi)*GaussR))+
                2*real(ScattAmplSin)*F1/(sqrt(Pi)*GaussR)-imag(ScattAmplSin)*F2/GaussR;
    //so far this is the eq. for singlet only, w/o QS

    //if we need to include the triplet, we add the term with a weight factor of 3 more than the singlet.
    //since however the correct norm. coeff. are 0.25 and 0.75 we need to divide by 4 to get the final result
    if(!SinOnly){
        complex<double> ScattAmplTri = pow(1./ScattLenTri+0.5*EffRangeTri*Momentum*Momentum-i*Momentum,-1.);
        CkValue +=  3*(
                    0.5*pow(abs(ScattAmplTri)/GaussR,2)*(1-(EffRangeTri)/(2*sqrt(Pi)*GaussR))+
                    2*real(ScattAmplTri)*F1/(sqrt(Pi)*GaussR)-imag(ScattAmplTri)*F2/GaussR);
        CkValue *= 0.25;
    }
    //if we have to include QS we need to add a correction factor and normalize by factor of 1/2
    if(QS){
        CkValue -= exp(-GaussR*GaussR*4.*Momentum*Momentum);
        CkValue *= 0.5;
    }

    CkValue += 1;
    CkValue *= LambdaLedni;
    CkValue += (1.-LambdaLedni)*(1+aResidual*exp(-ResidualR*ResidualR*4.*Momentum*Momentum));
    //CkValue += aResidual*exp(-ResidualR*ResidualR*4.*Momentum*Momentum);
//printf("ScattLenSin=%f; EffRangeSin=%f; GaussR=%f; LambdaLedni=%f; aResidual=%f; ResidualR=%f\n",
//       ScattLenSin*197.327, EffRangeSin*197.327, GaussR*197.327, LambdaLedni, aResidual, ResidualR*197.327);
    return CkValue;
}


