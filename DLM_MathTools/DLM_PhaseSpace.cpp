
#include "DLM_PhaseSpace.h"
#include "DLM_Random.h"

DLM_DecayGen::DLM_DecayGen(const unsigned& seed):
  SEED(seed),
  MaxNumDaughters(4),
  MinNumDaughters(2){

  RanGen = new DLM_Random(seed);
  DaughterMass = new double [MaxNumDaughters];

  SillyThreshold = 0;

}

~DLM_DecayGen::DLM_DecayGen(){
  if(RanGen){delete RanGen; RanGen=NULL;}
  if(DaughterMass){delete[]DaughterMass; DaughterMass=NULL;}
}

void DLM_DecayGen::SetEpsilon(const double& epsilon, const unsigned& steps){
  EpsilonStart = epsilon;
  //on each step Epsilon is halfed
  EpsSteps = steps;
}

void DLM_DecayGen::AddDaughter(const double& mass){
  DaughterMass[NumDaughters]=mass;
  NumDaughters++;
}

void DLM_DecayGen::SetSillyThreshold(const double& nthr){
  SillyThreshold = nthr;
}

void DLM_DecayGen::DirtyDist(const unsigned& numdecays){
  double Epsilon = EpsilonStart;
  const double EpsilonFinal = EpsilonStart*pow(2.,-double(EpsSteps));





  TLorentzVector Particle[3];
  TLorentzVector NotBad[3];
  TLorentzVector Sum;
  TLorentzVector Mother;

//for mu decay :     Phys.Rev.D 90 (2014) 9, 093002
//also EPJ Web of Conferences 118, 01021 (2016)
  const double MotherMass = Mmu;
  double TotMom;
  double Mass[3];
  //Mass[0] = 450;
  //Mass[1] = 250;
  //Mass[2] = 50;
  Mass[0] = Mel;
  Mass[1] = Mneutrino;
  Mass[2] = Mneutrino;

  const double TotMass = Mass[0]+Mass[1]+Mass[2];
  //const double ExcessE = 550.;//550
  //const double MotherMass = TotMass+ExcessE;
  const double ExcessE = MotherMass-Mass[0]-Mass[1]-Mass[2];
  double mag,costheta,sintheta,theta,phi,pT,eta,xcrd,ycrd,zcrd;

  double MaxMass = Mass[0];
  for(unsigned uPart=1; uPart<3; uPart++){
    if(Mass[uPart]>MaxMass) MaxMass=Mass[uPart];
  }
  const double MaxMom = sqrt(ExcessE*ExcessE+2.*MaxMass*ExcessE);
  printf("MM = %f\n",MotherMass);
  printf("M0 = %f\n",Mass[0]);
  printf("M1 = %f\n",Mass[1]);
  printf("M1 = %f\n",Mass[2]);
  printf("EE = %f\n",ExcessE);
  //printf("MaxMom = %f vs %f\n",MaxMom,sqrt(ExcessE*ExcessE/2./2.+2.*MaxMass*ExcessE/2.));

  Mother.SetXYZM(0,0,0,MotherMass);

  TH1F** hMomentum = new TH1F* [3];
  TH1F** hCosTheta = new TH1F* [3];
  TH1F** hPhi = new TH1F* [3];

  TH1F** hX = new TH1F* [3];
  TH1F** hY = new TH1F* [3];
  TH1F** hZ = new TH1F* [3];
  for(unsigned uPart=0; uPart<3; uPart++){
    TString Name;
    Name = TString::Format("hMomentum_%u_of_%u",uPart,3);
    hMomentum[uPart] = new TH1F(Name,Name,NumBins,0,MaxMom);
    Name = TString::Format("hCosTheta_%u_of_%u",uPart,3);
    hCosTheta[uPart] = new TH1F(Name,Name,NumBins,-1,1);
    Name = TString::Format("hPhi_%u_of_%u",uPart,3);
    hPhi[uPart] = new TH1F(Name,Name,NumBins,-2.*TMath::Pi(),2.*TMath::Pi());

    Name = TString::Format("hX_%u_of_%u",uPart,3);
    hX[uPart] = new TH1F(Name,Name,NumBins,-1.25*ExcessE,1.25*ExcessE);

    Name = TString::Format("hY_%u_of_%u",uPart,3);
    hY[uPart] = new TH1F(Name,Name,NumBins,-1.25*ExcessE,1.25*ExcessE);

    Name = TString::Format("hZ_%u_of_%u",uPart,3);
    hZ[uPart] = new TH1F(Name,Name,NumBins,-1.25*ExcessE,1.25*ExcessE);
  }
  unsigned TrueNumIter=NumIter;
  printf("Starting\n");
  for(unsigned uIter=0; uIter<TrueNumIter; uIter++){
    if(uIter%1000==0)
      printf("%u: ",uIter);
    Epsilon = EpsilonStart;
    //bool FirstSample = true;
bool RED_DEBUG = false;
    //int Counter = 0;
    for(unsigned uPart=0; uPart<3; uPart++){
      NotBad[uPart].SetXYZM(0,0,0,0);
    }
    int StuckCounter = 0;
    const int StuckLimit = 1000;

    do{
        //printf("NB %f %f %f\n", NotBad[0].E(),NotBad[1].E(),NotBad[2].E());
      if(NotBad[0].E()!=0||NotBad[1].E()!=0||NotBad[2].E()!=0){
        StuckCounter++;
        if(StuckCounter>=StuckLimit){
          for(unsigned uPart=0; uPart<3; uPart++) NotBad[uPart].SetXYZM(0,0,0,0);
          StuckCounter = 0;
          Epsilon = EpsilonStart;
          //printf("--reset--\n");
        }
      }
      for(unsigned uPart=0; uPart<3; uPart++){
        //if we dont have a good guess
        if(NotBad[uPart].E()==0){
if(NotBad[0].E()!=0||NotBad[1].E()!=0||NotBad[2].E()!=0){
  printf("VERY BAD\n");
  NotBad[0].Print();
  NotBad[1].Print();
  NotBad[2].Print();
  printf("%e %e %e\n", NotBad[0].E(),NotBad[1].E(),NotBad[2].E());
  usleep(1000e3);
}
          //printf("Working HARD\n");
//if(RED_DEBUG){
//printf("bad\n");
//}
          double VAL = sqrt(ExcessE*ExcessE+2.*Mass[uPart]*ExcessE);
          xcrd = -VAL+2.*VAL*rangen.Uniform();
          ycrd = -VAL+2.*VAL*rangen.Uniform();
          zcrd = -VAL+2.*VAL*rangen.Uniform();

        }
        //if we have a good guess
        else{
          //printf("NOTBAD\n");
          xcrd = rangen.Gaus(NotBad[uPart].Px(),ExcessE*Epsilon*1);
          ycrd = rangen.Gaus(NotBad[uPart].Py(),ExcessE*Epsilon*1);
          zcrd = rangen.Gaus(NotBad[uPart].Pz(),ExcessE*Epsilon*1);
          //}
          //else{
          //  costheta = NotBad[uPart].CosTheta();
          //  sintheta = sqrt(1.-costheta*costheta);
          //  eta = -log(sintheta/(1.+costheta+1e-99));
          //  phi = NotBad[uPart].Phi();
          //}
        }

        Particle[uPart].SetXYZM(xcrd,ycrd,zcrd,Mass[uPart]);
        //printf(" mag %f P %f\n",mag,Particle[uPart].P());
        //usleep(250e3);
      }
      Sum = Particle[0]+Particle[1]+Particle[2];
      //Counter++;
      //for(unsigned uPart=0; uPart<3; uPart++){
      //  hX[uPart]->Fill(Particle[uPart].Px());
      //  hY[uPart]->Fill(Particle[uPart].Py());
      //  hZ[uPart]->Fill(Particle[uPart].Pz());
      //}
      //if(Counter>10000) break;
      if(RED_DEBUG){
        printf("  DB %f vs %f && %f vs %f\n",Sum.M(),Mother.M(),Sum.P(),Mother.P());
        usleep(200e3);
      }

      if( (fabs(Sum.M()-Mother.M())<Epsilon*ExcessE && fabs(Sum.P()-Mother.P())<Epsilon*ExcessE) && Epsilon>EpsilonFinal ){
        Epsilon *= 0.5;
        for(unsigned uPart=0; uPart<3; uPart++){
          NotBad[uPart] = Particle[uPart];
        }
        //printf(" RED %f vs %f && %f vs %f\n",Sum.M(),Mother.M(),Sum.P(),Mother.P());
        //RED_DEBUG = true;
      }
    }
    while( fabs(Sum.M()-Mother.M())>EpsilonFinal*ExcessE || fabs(Sum.P()-Mother.P())>EpsilonFinal*ExcessE );

    if(uIter%1000==0) printf("  %f vs %f && %f vs %f\n",Sum.M(),Mother.M(),Sum.P(),Mother.P());
    //usleep(200e3);

    //make it fit
    unsigned WhichOne = rangen.Integer(3);
    Particle[WhichOne] = (Mother-Particle[(WhichOne+1)%3]-Particle[(WhichOne+2)%3]);
    if(Particle[WhichOne].Mag2()<0 && !IgnoreSilly){
      TrueNumIter++;
      continue;
    }

    for(unsigned uPart=0; uPart<3; uPart++){
      hMomentum[uPart]->Fill(Particle[uPart].P());
      hCosTheta[uPart]->Fill(Particle[uPart].CosTheta());
      hPhi[uPart]->Fill(Particle[uPart].Phi());
        hX[uPart]->Fill(Particle[uPart].Px());
        hY[uPart]->Fill(Particle[uPart].Py());
        hZ[uPart]->Fill(Particle[uPart].Pz());
    }
  }

  printf("Finished with %u silly events\n",TrueNumIter-NumIter);
  TFile fOutput("XyzClever_3body.root","recreate");
  for(unsigned uPart=0; uPart<3; uPart++){

    hMomentum[uPart]->Scale(1./hMomentum[uPart]->Integral(),"width");
    hCosTheta[uPart]->Scale(1./hCosTheta[uPart]->Integral(),"width");
    hPhi[uPart]->Scale(1./hPhi[uPart]->Integral(),"width");
    hX[uPart]->Scale(1./hX[uPart]->Integral(),"width");
    hY[uPart]->Scale(1./hY[uPart]->Integral(),"width");
    hZ[uPart]->Scale(1./hZ[uPart]->Integral(),"width");

    hMomentum[uPart]->Write();
    hCosTheta[uPart]->Write();
    hPhi[uPart]->Write();
    hX[uPart]->Write();
    hY[uPart]->Write();
    hZ[uPart]->Write();
  }


}
