#include "DLM_CRAB_PM.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TString.h"

//this phase maker is based on the code of Scott Pratt: https://web.pa.msu.edu/people/pratts/freecodes/crab/home.html
//and written my Dimitar Mihaylov

DLM_CRAB_PM::DLM_CRAB_PM():Pi(3.14159265358979323844){
    //rx = NULL;
    //ry = NULL;
    //rz = NULL;
    //tau = NULL;
    SpaceCoord = NULL;
    //ShapeRx = NULL;
    //ShapeRy = NULL;
    //ShapeRz = NULL;
    //ShapeTau = NULL;
    Shape = NULL;
    NumEvents = 10000;
    NumPerEvent = 20;
    NumSpecies = 0;
    Temperature = 0;
    ParticleSpecies = NULL;
    ParticleMass = NULL;
    CDF = NULL;
    IntegratedSource = NULL;
    OutputFileName = "";
}

DLM_CRAB_PM::~DLM_CRAB_PM(){
    MemCleanup();
}

void DLM_CRAB_PM::MemInit(const unsigned& num){
    MemCleanup();
    ParticleSpecies = new int [num];
    ParticleMass = new double [num];
    //rx = new double [NumSpecies];
    //ry = new double [NumSpecies];
    //rz = new double [NumSpecies];
    //tau = new double [NumSpecies];
    SpaceCoord = new double* [num];

    CDF = new TH1F** [num];
    //ShapeRx = new char [NumSpecies];
    //ShapeRy = new char [NumSpecies];
    //ShapeRz = new char [NumSpecies];
    //ShapeTau = new char [NumSpecies];
    Shape = new char* [num];
    for(unsigned ui=0; ui<num; ui++){
        SpaceCoord[ui] = new double [4];
        Shape[ui] = new char [4];
        CDF[ui] = new TH1F* [4];
        Shape[ui][0] = 0;
        Shape[ui][1] = 0;
        Shape[ui][2] = 0;
        Shape[ui][3] = 0;
        CDF[ui][0] = NULL;
        CDF[ui][1] = NULL;
        CDF[ui][2] = NULL;
        CDF[ui][3] = NULL;
    }
    NumSpecies = num;
}

void DLM_CRAB_PM::MemCleanup(){
    if(ParticleSpecies){
        delete [] ParticleSpecies; ParticleSpecies = NULL;
    }
    if(ParticleMass){
        delete [] ParticleMass; ParticleMass = NULL;
    }
    if(CDF){
        for(unsigned ui=0; ui<NumSpecies; ui++){
            CDF[ui][0] = NULL;
            CDF[ui][1] = NULL;
            CDF[ui][2] = NULL;
            CDF[ui][3] = NULL;
            delete [] CDF[ui];
        }
        delete [] CDF; CDF = NULL;
    }
    if(IntegratedSource){
        for(unsigned ui=0; ui<NumSpecies; ui++){
            delete IntegratedSource[ui]; IntegratedSource[ui] = NULL;
        }
        delete [] IntegratedSource; IntegratedSource = NULL;
    }
    /*
    if(rx){
        delete [] rx; rx = NULL;
    }
    if(ry){
        delete [] ry; ry = NULL;
    }
    if(rz){
        delete [] rz; rz = NULL;
    }
    if(tau){
        delete [] tau; tau = NULL;
    }
    if(ShapeRx){
        delete [] ShapeRx; ShapeRx = NULL;
    }
    if(ShapeRy){
        delete [] ShapeRy; ShapeRy = NULL;
    }
    if(ShapeRz){
        delete [] ShapeRz; ShapeRz = NULL;
    }
    if(ShapeTau){
        delete [] ShapeTau; ShapeTau = NULL;
    }
    */
    if(SpaceCoord){
        for(unsigned ui=0; ui<NumSpecies; ui++){
            delete [] SpaceCoord[ui];
        }
        delete [] SpaceCoord; SpaceCoord = NULL;
    }
    if(Shape){
        for(unsigned ui=0; ui<NumSpecies; ui++){
            delete [] Shape[ui];
        }
        delete [] Shape; Shape = NULL;
    }
}

void DLM_CRAB_PM::SetNumEvents(unsigned num){
    NumEvents = num;
}

void DLM_CRAB_PM::SetNumPerEvent(unsigned num){
    NumPerEvent = num;
}

void DLM_CRAB_PM::SetOutputFile(string OutputName){
    OutputFileName = OutputName;
}

void DLM_CRAB_PM::SetNumSpecies(unsigned num){
    if(num==NumSpecies) return;
    MemInit(num);
}

void DLM_CRAB_PM::SetParticleMass(unsigned np, double mass){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    ParticleMass[np] = mass;
}

void DLM_CRAB_PM::SetParticleMass(const double* mass){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        ParticleMass[ui] = mass[ui];
    }
}

void DLM_CRAB_PM::SetParticleSpecies(unsigned np, int species){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    ParticleSpecies[np] = species;
}

void DLM_CRAB_PM::SetParticleSpecies(const int* species){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        ParticleSpecies[ui] = species[ui];
    }
}

void DLM_CRAB_PM::SetTemperature(double temp){
    Temperature = temp;
}

void DLM_CRAB_PM::SetRx(unsigned np, double value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    SpaceCoord[np][rx] = value;
}

void DLM_CRAB_PM::SetRx(const double* value){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        SpaceCoord[ui][rx] = value[ui];
    }
}

void DLM_CRAB_PM::SetRy(unsigned np, double value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    SpaceCoord[np][ry] = value;
}

void DLM_CRAB_PM::SetRy(const double* value){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        SpaceCoord[ui][ry] = value[ui];
    }
}

void DLM_CRAB_PM::SetRz(unsigned np, double value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    SpaceCoord[np][rz] = value;
}

void DLM_CRAB_PM::SetRz(const double* value){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        SpaceCoord[ui][rz] = value[ui];
    }
}

void DLM_CRAB_PM::SetRxyz(unsigned np, double value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    SpaceCoord[np][rx] = value;
    SpaceCoord[np][ry] = value;
    SpaceCoord[np][rz] = value;
}

void DLM_CRAB_PM::SetRxyz(const double* value){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        SpaceCoord[ui][rx] = value[ui];
        SpaceCoord[ui][ry] = value[ui];
        SpaceCoord[ui][rz] = value[ui];
    }
}

void DLM_CRAB_PM::SetTau(unsigned np, double value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    SpaceCoord[np][tau] = value;
}

void DLM_CRAB_PM::SetTau(const double* value){
    for(unsigned ui=0; ui<NumSpecies; ui++){
        SpaceCoord[ui][tau] = value[ui];
    }
}

void DLM_CRAB_PM::SetShapeRx(unsigned np, char value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
}

void DLM_CRAB_PM::SetShapeRx(const char* value){
}

void DLM_CRAB_PM::SetShapeRx(unsigned np, TH1F* cdf){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    if(!cdf) {printf("ERROR! Bad CDF (NULL)!\n"); return;}
    CDF[np][rx] = cdf;
}

void DLM_CRAB_PM::SetShapeRy(unsigned np, char value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
}

void DLM_CRAB_PM::SetShapeRy(const char* value){
}

void DLM_CRAB_PM::SetShapeRy(unsigned np, TH1F* cdf){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    if(!cdf) {printf("ERROR! Bad CDF (NULL)!\n"); return;}
    CDF[np][ry] = cdf;
}


void DLM_CRAB_PM::SetShapeRz(unsigned np, char value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
}

void DLM_CRAB_PM::SetShapeRz(const char* value){
}

void DLM_CRAB_PM::SetShapeRz(unsigned np, TH1F* cdf){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    if(!cdf) {printf("ERROR! Bad CDF (NULL)!\n"); return;}
    CDF[np][rz] = cdf;
}


void DLM_CRAB_PM::SetShapeRxyz(unsigned np, char value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
}

void DLM_CRAB_PM::SetShapeRxyz(const char* value){
}

void DLM_CRAB_PM::SetShapeRxyz(unsigned np, TH1F* cdf){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    if(!cdf) {printf("ERROR! Bad CDF (NULL)!\n"); return;}
    CDF[np][rx] = cdf;
    CDF[np][ry] = cdf;
    CDF[np][rz] = cdf;
}


void DLM_CRAB_PM::SetShapeTau(unsigned np, char value){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
}

void DLM_CRAB_PM::SetShapeTau(const char* value){
}

void DLM_CRAB_PM::SetShapeTau(unsigned np, TH1F* cdf){
    if(np>=NumSpecies) {printf("ERROR! There are only %u particle species!\n", NumSpecies); return;}
    if(!cdf) {printf("ERROR! Bad CDF (NULL)!\n"); return;}
    CDF[np][tau] = cdf;
}

void DLM_CRAB_PM::RunPhasemaker(unsigned Seed){

    if(IntegratedSource){
        for(unsigned ui=0; ui<NumSpecies; ui++){
            delete IntegratedSource[ui]; IntegratedSource[ui] = NULL;
        }
        delete [] IntegratedSource; IntegratedSource = NULL;
    }
    IntegratedSource = new TH1F* [NumSpecies];
    for(unsigned ui=0; ui<NumSpecies; ui++){
        IntegratedSource[ui] = new TH1F(TString::Format("IntegratedSource%u",ui),TString::Format("IntegratedSource%u",ui),
                                        int(double(NumPerEvent)*double(NumEvents)*double(NumSpecies)/128.),0,64);
    }

    double p[4],r[4];
    double weight,pmag2,pmag,cthet,sthet,phi;
    unsigned i,id,alpha;
    TRandom3 rangen(Seed);
    FILE *fptr;
    if(OutputFileName!="") fptr=fopen(OutputFileName.c_str(),"w");
    else fptr=NULL;
    if(fptr){
        fprintf(fptr,"OSC1997A\n");
        fprintf(fptr,"final_id_p_x\n");
        fprintf(fptr,"   DLM_CRAB_PM --- Phasemaker ---\n");
    }

    for(unsigned uEvents=0; uEvents<NumEvents; uEvents++){
        if(fptr) fprintf(fptr,"%10u%12u%10.3lf%10.3lf\n",uEvents+1,NumPerEvent,0.0,0.0);
        for(i=0;i<NumPerEvent;i++){
            id=(unsigned)floor(NumSpecies*rangen.Uniform());
            if(ParticleMass[id]/Temperature>8.0){
            pmag2=0.0;
            p[kx] = rangen.Gaus();
            p[ky] = rangen.Gaus();
            p[kz] = rangen.Gaus();
            for(alpha=1;alpha<4;alpha++) {
                p[alpha]=sqrt(Temperature*ParticleMass[id])*p[alpha];
                pmag2=pmag2+p[alpha]*p[alpha];
            }
            p[energy]=sqrt(pmag2+ParticleMass[id]*ParticleMass[id]);
            }
            else{
                do{
                    pmag=-Temperature*log(rangen.Uniform()*rangen.Uniform()*rangen.Uniform());
                    p[energy]=sqrt(ParticleMass[id]*ParticleMass[id]+pmag*pmag);
                    weight=exp(-(p[energy]-pmag)/Temperature);
                }
                while(rangen.Uniform()>weight);
                cthet=1.0-2.0*rangen.Uniform();
                sthet=sqrt(1.0-cthet*cthet);
                phi=2.0*Pi*rangen.Uniform();
                p[kz]=pmag*cthet;
                p[kx]=pmag*sthet*cos(phi);
                p[ky]=pmag*sthet*sin(phi);
            }
    /*
            r[tau] = rangen.Gaus();
            r[rx] = rangen.Gaus();
            r[ry] = rangen.Gaus();
            r[rz] = rangen.Gaus();

            r[tau]=SpaceCoord[id][tau]*r[tau];
            r[rx]=SpaceCoord[id][rx]*r[rx];
            r[ry]=SpaceCoord[id][ry]*r[ry];
            r[rz]=SpaceCoord[id][rz]*r[rz];
    */

            //for(unsigned ui=0; ui<4; ui++){
            //    r[ui] = GetRandSpace(rangen, ui);
            //}
            r[tau] = GetRandSpace(rangen, id, tau);
            r[rx] = GetRandSpace(rangen, id, rx);
            r[ry] = GetRandSpace(rangen, id, ry);
            r[rz] = GetRandSpace(rangen, id, rz);

            if(fptr) fprintf(fptr,"%10u%12i%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e\n",
                i+1,ParticleSpecies[id],p[kx]/1000.0,p[ky]/1000.0,p[kz]/1000.0,p[energy]/1000.0,
                ParticleMass[id]/1000.0,r[rx],r[ry],r[rz],r[tau]);

            IntegratedSource[id]->Fill(sqrt(r[rx]*r[rx]+r[ry]*r[ry]+r[rz]*r[rz]));
        }
    }

    for(unsigned us=0; us<NumSpecies; us++){
        IntegratedSource[us]->Scale(1./IntegratedSource[us]->Integral(1,IntegratedSource[us]->GetNbinsX()),"width");
    }

    if(fptr) fclose(fptr);
}

TH1F* DLM_CRAB_PM::GetIntegratedSource(const unsigned& Species){
    if(Species>=NumSpecies) return NULL;
    return IntegratedSource[Species];
}

double DLM_CRAB_PM::GetRandSpace(TRandom3& rangen, unsigned Species, int Coord){
    while(true){
    if(!CDF[Species][Coord]){
        printf("ERROR! Some pdf is not set!\n");
        return 0;
    }
    int NumBins = CDF[Species][Coord]->GetNbinsX();
    if(NumBins<1){
        printf("ERROR! Some pdf is wrongly defined!\n");
        return 0;
    }

    double ExtrapolatedValueX;
    double MinRange = CDF[Species][Coord]->GetBinLowEdge(1);
    double MaxRange = CDF[Species][Coord]->GetXaxis()->GetBinUpEdge(NumBins);
    double RandPositionY = rangen.Uniform();

    int SecondBin=NumBins+1;

    for(int iBin=1; iBin<=NumBins; iBin++){
        if(RandPositionY<CDF[Species][Coord]->GetBinContent(iBin)){
            SecondBin = iBin;
            break;
        }
    }

    double X1,X2,Y1,Y2;
    if(SecondBin==1){
        X1=MinRange;
        Y1=0;
        X2=CDF[Species][Coord]->GetBinCenter(SecondBin);
        Y2=CDF[Species][Coord]->GetBinContent(SecondBin);
    }
    else if(SecondBin>NumBins){
        X1=CDF[Species][Coord]->GetBinCenter(NumBins);
        Y1=CDF[Species][Coord]->GetBinContent(NumBins);
        X2=MaxRange;
        Y2=1;
    }
    else{
        X1=CDF[Species][Coord]->GetBinCenter(SecondBin-1);
        Y1=CDF[Species][Coord]->GetBinContent(SecondBin-1);
        X2=CDF[Species][Coord]->GetBinCenter(SecondBin);
        Y2=CDF[Species][Coord]->GetBinContent(SecondBin);
    }
    if(RandPositionY<Y1 || RandPositionY>Y2){
        printf("There is something fishy! -- if(RandPositionY<Y1 || RandPositionY>Y2) --\n");
    }

    ExtrapolatedValueX = (X2-X1)*RandPositionY/(Y2-Y1) + (Y2*X1-Y1*X2)/(Y2-Y1);

    if(ExtrapolatedValueX<X1 || ExtrapolatedValueX>X2){
        //printf("There is something fishy! -- if(ExtrapolatedValueX<X1 || ExtrapolatedValueX>X2) --\n");
    }
    else{
       return ExtrapolatedValueX;
    }

    }

}

/*
double DLM_CRAB_PM::GetRandSpace(TRandom3& rangen, unsigned Species, int Coord){
    if(!CDF[Species][Coord]){
        printf("ERROR! Some pdf is not set!\n");
        return 0;
    }
    int NumBins = CDF[Species][Coord]->GetNbinsX();
    if(NumBins<1){
        printf("ERROR! Some pdf is wrongly defined!\n");
        return 0;
    }
    //the pdf is extrapolated around the bin
    //for that one chooses two points. In case we deal with the first half of the first bin
    //or the second half of the last bin then
    //the extrapolation is done between 0 and the value of the bin.
    //for all other cases the extrapolation is done between:
    //1) the current and the previous bin in case the current value sits below the middle of the bin
    //2) the current and the next bin in case the current value sits after the middle of the bin
    double FirstPointX,FirstPointY;
    double SecondPointX,SecondPointY;
    //double Distance;
    double ExtrapolatedValue;

    double MinRange = CDF[Species][Coord]->GetBinLowEdge(1);
    double MaxRange = CDF[Species][Coord]->GetXaxis()->GetBinUpEdge(NumBins);

    double RandPositionX = rangen.Uniform(MinRange, MaxRange);
    int InWhichBin = CDF[Species][Coord]->FindBin(RandPositionX);
    double BinCenter = CDF[Species][Coord]->GetBinCenter(InWhichBin);
    //double PositionY = CDF[Species][Coord]->GetBinContent(InWhichBin);
    if( InWhichBin==1 && RandPositionX<BinCenter){
        FirstPointX = MinRange;
        FirstPointY = 0;
        SecondPointX = BinCenter;
        SecondPointY = CDF[Species][Coord]->GetBinContent(InWhichBin);
    }
    else if(InWhichBin==NumBins && RandPositionX>BinCenter){
        FirstPointX = BinCenter;
        FirstPointY = CDF[Species][Coord]->GetBinContent(InWhichBin);
        SecondPointX = MaxRange;
        SecondPointY = 0;
    }
    else if(RandPositionX<BinCenter){
        FirstPointX = CDF[Species][Coord]->GetBinCenter(InWhichBin-1);
        FirstPointY = CDF[Species][Coord]->GetBinContent(InWhichBin-1);
        SecondPointX = BinCenter;
        SecondPointY = CDF[Species][Coord]->GetBinContent(InWhichBin);
    }
    else{
        FirstPointX = BinCenter;
        FirstPointY = CDF[Species][Coord]->GetBinContent(InWhichBin);
        SecondPointX = CDF[Species][Coord]->GetBinCenter(InWhichBin+1);
        SecondPointY = CDF[Species][Coord]->GetBinContent(InWhichBin+1);
    }

    ExtrapolatedValue = ((SecondPointY-FirstPointY)*RandPositionX+SecondPointX*FirstPointY-FirstPointX*SecondPointY)/(SecondPointX-FirstPointX);

    return ExtrapolatedValue;
}
*/

