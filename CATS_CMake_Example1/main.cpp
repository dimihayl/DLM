#include <iostream>

#include "CATStools.h"
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TFile.h"

#include "DimiFit.h"

using namespace std;

void Cat1(){

	const unsigned NumMomBins=40;
	const double kMin=0;
	const double kMax=160;

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
    double PotPars1S0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,0,0,0};
    double PotPars3P0[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,0};
    double PotPars3P1[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,1};
    double PotPars3P2[10]={0,0,NN_AV18,v18_Coupled3P2,1,1,1,1,1,2};

    //CatsLorentzVector lv1;
    CATS MyCat;
    double Pars[4] = {0,0,0,1.2};
    MyCat.SetAnaSource(GaussSource, Pars);
    MyCat.SetUseAnalyticSource(true);

    MyCat.SetExcludeFailedBins(false);
    MyCat.SetMomBins(NumMomBins,kMin,kMax);

    MyCat.SetNumChannels(4);
    MyCat.SetNumPW(0,2);
    MyCat.SetNumPW(1,2);
    MyCat.SetNumPW(2,2);
    MyCat.SetNumPW(3,2);
    MyCat.SetSpin(0,0);
    MyCat.SetSpin(1,1);
    MyCat.SetSpin(2,1);
    MyCat.SetSpin(3,1);
    MyCat.SetChannelWeight(0, 3./12.);
    MyCat.SetChannelWeight(1, 1./12.);
    MyCat.SetChannelWeight(2, 3./12.);
    MyCat.SetChannelWeight(3, 5./12.);

    MyCat.SetQ1Q2(1);
    MyCat.SetPdgId(2212,2212);
    double Mass1=938.272; double Mass2=938.272;
    MyCat.SetRedMass( (Mass1*Mass2)/(Mass1+Mass2) );

    MyCat.SetShortRangePotential(0,0,fDlmPot,PotPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,PotPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,PotPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,PotPars3P2);

    MyCat.KillTheCat();

	TGraph gr1;
	gr1.SetName("gr1");
	gr1.Set(40);

    for(unsigned uBin=0; uBin<NumMomBins; uBin++){
        printf("C(%.2f)=%.2f\n",MyCat.GetMomentum(uBin),MyCat.GetCorrFun(uBin));
        gr1.SetPoint(uBin,MyCat.GetMomentum(uBin),MyCat.GetCorrFun(uBin));
    }

    TFile* f1 = new TFile("f1.root", "recreate");
    f1->cd();
    gr1.Write();
	delete f1;
}

int main()
{
    cout << "Hello world!" << endl;

    //Cat1();

	CALL_BERNIE_AND_VALE();

    return 0;
}
