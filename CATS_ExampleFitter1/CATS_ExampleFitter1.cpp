
#include "CATS_ExampleFitter1.h"

#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "CATS.h"

CATS* CkCats;

double GaussSrc(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];
    double& Size = Pars[3];
    return 4.*TMath::Pi()*Radius*Radius*pow(4.*TMath::Pi()*Size*Size,-1.5)*exp(-(Radius*Radius)/(4.*Size*Size));
}

double DoubleGaussSrc(double* Pars){
    //double& Momentum = Pars[0];
    double& Radius = Pars[1];
    //double& CosTheta = Pars[2];

    double& Size1 = Pars[3];
    double& Size2 = Pars[4];
    double& Weight1 = Pars[5];

    return      Weight1 *4.*TMath::Pi()*Radius*Radius*pow(4.*TMath::Pi()*Size1*Size1,-1.5)*exp(-(Radius*Radius)/(4.*Size1*Size1))+
            (1.-Weight1)*4.*TMath::Pi()*Radius*Radius*pow(4.*TMath::Pi()*Size2*Size2,-1.5)*exp(-(Radius*Radius)/(4.*Size2*Size2));

}

//(a+b*k)*(lam*C(k)+1-lam)
//the source is either Gauss(r1) or f*Gauss(r1)+(1-f)*Gauss(r2)
double CatsFitFunction(double* mom, double* par){
    double& k=*mom;
    double& a=par[0];
    double& b=par[1];
    double& lam=par[2];
    double& r1=par[3];
    double& r2=par[4];
    double& f=par[5];

    if(CkCats->GetUseAnalyticSource()){
        //the 3rd parameter makes sure that the computing grid is reused and not initialized each time the function is called
        //as the grid is initialized dynamically based on the source size, this only works well for "fine-tuning",
        //i.e. when the parameters do not change much, hence a good choice of starting parameters is important
        CkCats->SetAnaSource(0, r1, true);
        CkCats->SetAnaSource(1, r2, true);
        CkCats->SetAnaSource(2, f, true);
    }
    //one could use the source from a transport model, not implemented here though
    else{
        return (a+b*k);
    }
    CkCats->KillTheCat();

    return (a+b*k)*(lam*CkCats->EvalCorrFun(k)+1.-lam);
}

void ChargedPionPionFitter(TH1F* hData){

    const unsigned NumMomBins=25;
    const double kMin = 0;
    const double kMax = 500;

    const double Mass_Pi=139.57;
    enum SOURCE { sGauss, sDoubleGauss };
    const int WhichSource = sDoubleGauss;

    double SourceParameters[6] = {0,0,0,0.9,3.0,0.4};

    CATS Cat;

    Cat.SetExcludeFailedBins(false);
    Cat.SetMomBins(NumMomBins,kMin,kMax);

    //standard settings for a CATS object which has no strong interaction potential included
    Cat.SetNumChannels(1);
    Cat.SetNumPW(0,0);
    Cat.SetSpin(0,0);
    Cat.SetChannelWeight(0, 1);

    //include the coulomb interaction. Q1Q2 is the multiplied charge numbers of the two particles
    Cat.SetQ1Q2(1);
    //the PdgId is needed when using a source from a transport model
    Cat.SetPdgId(211, 211);

    Cat.SetRedMass( 0.5*Mass_Pi );

    Cat.SetMaxRad(100);
    Cat.SetMaxRho(40);

    //standard settings for CATS depending on the source
    if(WhichSource==sDoubleGauss){
        Cat.SetAnaSource(DoubleGaussSrc, SourceParameters);
        Cat.SetUseAnalyticSource(true);
    }
    else if(WhichSource==sGauss){
        Cat.SetAnaSource(GaussSrc, SourceParameters);
        Cat.SetUseAnalyticSource(true);
    }
    else{
        printf("ERROR - Unknown source.\n");
        return;
    }

    //sets up the whole object
    Cat.KillTheCat();
    //makes sure further reiterations only display errors and warning messages
    Cat.SetNotifications(CATS::nWarning);
    //a pointer that will be used by the fitter
    CkCats = &Cat;

    //(a+b*k)*(lam*C(k)+1-lam)
    //the source is either Gauss(r1) or f*Gauss(r1)+(1-f)*Gauss(r2)
    //Fit parameters:
    //[0] = a
    //[1] = b
    //[2] = lam (the correlation strength)
    //[3] = r1 for Gauss and double Gauss
    //[4] = r2 for double Gauss
    //[5] = f for double Gauss
    TF1* Fitter = new TF1 ("Fitter", CatsFitFunction, kMin, kMax, 6);

    //set up all of the parameters below
    Fitter->SetParameter(0, 1);
    Fitter->SetParLimits(0, -2, 2);

    Fitter->FixParameter(1, 0);

    Fitter->SetParameter(2, 0.5);
    Fitter->SetParLimits(2, 0, 1);

    switch(WhichSource){
    case sGauss:
        Fitter->SetParameter(3, 1);
        Fitter->SetParLimits(3, 0.5, 4);
        Fitter->FixParameter(4, 0);
        Fitter->FixParameter(5, 0);
        break;
    case sDoubleGauss:
        Fitter->SetParameter(3, 0.9);
        Fitter->SetParLimits(3, 0.4, 1.5);
        Fitter->SetParameter(4, 3);
        Fitter->SetParLimits(4, 1.5, 5);
        Fitter->SetParameter(5, 0.4);
        Fitter->SetParLimits(5, 0, 1);
        break;
    default:
        break;
    }

    printf("Fitting the data...\n");
    hData->Fit(Fitter, "S, N, R, M");
    printf("Fitting done!\n");
    printf(" chi2ndf = %.2f\n",Fitter->GetChisquare()/double(Fitter->GetNDF()));
    printf(" r1 = %.2f +/- %.2f\n",Fitter->GetParameter(3),Fitter->GetParError(3));
    if(WhichSource==sDoubleGauss){
        printf(" r2 = %.2f +/- %.2f\n",Fitter->GetParameter(4),Fitter->GetParError(4));
        printf(" f = %.2f +/- %.2f\n",1.-Fitter->GetParameter(5),Fitter->GetParError(5));
    }

    //TEST
    TFile* f1 = new TFile("./f1.root","recreate");
    hData->Write();
    Fitter->Write();
    delete f1;

    delete Fitter;
}


void PionPionExample1(){

    TH1F* hDummyData = new TH1F("hDummyData", "hDummyData", 25, 0, 500);

    hDummyData->SetBinContent(1, 1.4);
    hDummyData->SetBinContent(2, 1.3);
    hDummyData->SetBinContent(3, 1.2);
    hDummyData->SetBinContent(4, 1.15);
    hDummyData->SetBinContent(5, 1.1);
    hDummyData->SetBinContent(6, 1.05);
    hDummyData->SetBinContent(7, 1.03);
    hDummyData->SetBinContent(8, 1.01);
    for(unsigned uBin=9; uBin<=25; uBin++) hDummyData->SetBinContent(uBin, 1.0);

    for(unsigned uBin=1; uBin<=25; uBin++) hDummyData->SetBinError(uBin, 0.005);

    ChargedPionPionFitter(hDummyData);

    delete hDummyData;
}


