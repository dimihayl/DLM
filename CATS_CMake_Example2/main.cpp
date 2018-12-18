#include <iostream>

#include "CATStools.h"
#include "CATS.h"
#include "DLM_Potentials.h"
#include "DLM_Source.h"

#include "TGraph.h"
#include "TFile.h"

using namespace std;

void Cat1(){

	const unsigned NumMomBins=100;
	const double kMin=0;
	const double kMax=400;


	
    //keep in mind the units are always MeV or fm
    CATS MyCat;
    //specify that we want to use an analytic source
    //the alternative would be transport model (EPOS)
    MyCat.SetUseAnalyticSource(true);
	//define the parameters which will be needed by the source function
	//type,num pars,threadsafe
    CATSparameters cPars(CATSparameters::tSource,1,false);
    cPars.SetParameter(0,1.2);
    //pass to cats the source function and the parameters
    MyCat.SetAnaSource(GaussSource, cPars);

    MyCat.SetExcludeFailedBins(false);
    //defines your k-bins
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

    //#,#,POT_ID,POT_FLAG,t_tot,t1,t2,s,l,j
	CATSparameters cPotPars1S0(CATSparameters::tPotential,8,false);
	cPotPars1S0.SetParameter(0,NN_AV18);
	cPotPars1S0.SetParameter(1,v18_Coupled3P2);
	cPotPars1S0.SetParameter(2,1);
	cPotPars1S0.SetParameter(3,1);
	cPotPars1S0.SetParameter(4,1);
	cPotPars1S0.SetParameter(5,0);
	cPotPars1S0.SetParameter(6,0);
	cPotPars1S0.SetParameter(7,0);

	CATSparameters cPotPars3P0(CATSparameters::tPotential,8,false);
	cPotPars3P0.SetParameter(0,NN_AV18);
	cPotPars3P0.SetParameter(1,v18_Coupled3P2);
	cPotPars3P0.SetParameter(2,1);
	cPotPars3P0.SetParameter(3,1);
	cPotPars3P0.SetParameter(4,1);
	cPotPars3P0.SetParameter(5,1);
	cPotPars3P0.SetParameter(6,1);
	cPotPars3P0.SetParameter(7,0);
	
	CATSparameters cPotPars3P1(CATSparameters::tPotential,8,false);
	cPotPars3P1.SetParameter(0,NN_AV18);
	cPotPars3P1.SetParameter(1,v18_Coupled3P2);
	cPotPars3P1.SetParameter(2,1);
	cPotPars3P1.SetParameter(3,1);
	cPotPars3P1.SetParameter(4,1);
	cPotPars3P1.SetParameter(5,1);
	cPotPars3P1.SetParameter(6,1);
	cPotPars3P1.SetParameter(7,1);
	
	CATSparameters cPotPars3P2(CATSparameters::tPotential,8,false);
	cPotPars3P2.SetParameter(0,NN_AV18);
	cPotPars3P2.SetParameter(1,v18_Coupled3P2);
	cPotPars3P2.SetParameter(2,1);
	cPotPars3P2.SetParameter(3,1);
	cPotPars3P2.SetParameter(4,1);
	cPotPars3P2.SetParameter(5,1);
	cPotPars3P2.SetParameter(6,1);
	cPotPars3P2.SetParameter(7,2);
	
	//which channel, which partial wave, which potential function, which parameters
    MyCat.SetShortRangePotential(0,0,fDlmPot,cPotPars1S0);
    MyCat.SetShortRangePotential(1,1,fDlmPot,cPotPars3P0);
    MyCat.SetShortRangePotential(2,1,fDlmPot,cPotPars3P1);
    MyCat.SetShortRangePotential(3,1,fDlmPot,cPotPars3P2);

    MyCat.KillTheCat();

	TGraph gr1;
	gr1.SetName("gr1");
	gr1.Set(NumMomBins);

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

    Cat1();

    return 0;
}
