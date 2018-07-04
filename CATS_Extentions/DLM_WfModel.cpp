#include <iostream>
#include "stdio.h"
#include <string.h>
#include <math.h>

#include "DLM_WfModel.h"
#include "CATS.h"

using namespace std;

DLM_WfModels::DLM_WfModels(WfModel Model):MODEL(Model),
    NumChannels(Model==ProtonLambdaHaideNLO?4:Model==KminProtonHaide?4:Model==KplusProtonHaide?2:1){
    //NumChannels(Model==ProtonLambdaHaideNLO?4:Model==KminProtonHaide?4:Model==KplusProtonFiles?2:1){

    WaveFunctionU = NULL;
    PhaseShifts = NULL;
    RadBins = NULL;
    NumRadBins = 0;

}

DLM_WfModels::~DLM_WfModels(){
    if(WaveFunctionU){
        delete [] WaveFunctionU;
        WaveFunctionU = NULL;
    }
    if(PhaseShifts){
        delete [] PhaseShifts;
        PhaseShifts = NULL;
    }
    if(RadBins){
        delete [] RadBins;
        RadBins = NULL;
    }

}

unsigned DLM_WfModels::GetNumRadBins(){
    return NumRadBins;
}

void DLM_WfModels::Init(const CATS& Kitty, const char* inFileName){

}

void InitHaidenbauerNLO(const char* InputFolder, CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins){

    unsigned NumMomBins=Kitty.GetNumMomBins();
    //double kMin = Kitty.GetMomBinLowEdge(0);
    //double kMax = Kitty.GetMomBinUpEdge(NumMomBins-1);

    const unsigned short NumChannels = 4;
    const unsigned short NumPwPerCh = 2;
    const unsigned short NumFiles = 6;
    enum HaideFiles {f1S0, f1P1, f3S1, f3P0, f3P1, f3P2};

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
    }

    strcpy(InputFileName[f1S0],InputFolder);
    strcpy(InputFileName[f1P1],InputFolder);
    strcpy(InputFileName[f3S1],InputFolder);
    strcpy(InputFileName[f3P0],InputFolder);
    strcpy(InputFileName[f3P1],InputFolder);
    strcpy(InputFileName[f3P2],InputFolder);

    strcat(InputFileName[f1S0], "N1s0.data");
    strcat(InputFileName[f1P1], "N1p1.data");
    strcat(InputFileName[f3S1], "N3s1.data");
    strcat(InputFileName[f3P0], "N3p0.data");
    strcat(InputFileName[f3P1], "N3p1.data");
    strcat(InputFileName[f3P2], "N3p2.data");

    const unsigned NumBlankTransitionLines = 1;
    const double MinRad = 0.2;
    const double RadStep = 0.1;
    const double MaxRad = 16.1;
    NumRadBins = round((MaxRad-MinRad)/RadStep) + 1;
    bool* RadBinLoaded = new bool [NumRadBins+1];

    WaveFunctionU[0] = new double*** [NumMomBins];
    PhaseShifts[0] = new double** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        WaveFunctionU[0][uMomBin] = new double** [NumChannels];
        PhaseShifts[0][uMomBin] = new double* [NumChannels];
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            WaveFunctionU[0][uMomBin][usCh] = new double* [NumPwPerCh];
            PhaseShifts[0][uMomBin][usCh] = new double [NumPwPerCh];
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                WaveFunctionU[0][uMomBin][usCh][usPw] = new double [NumRadBins];
                PhaseShifts[0][uMomBin][usCh][usPw] = 0;
            }
        }
    }

    RadBins[0] = new double [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        RadBins[0][uRad] = MinRad-RadStep*0.5 + RadStep*double(uRad);
        RadBinLoaded[uRad] = false;
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return;
        }

        fseek ( InFile , 0 , SEEK_END );
        //long EndPos;
        //EndPos = ftell (InFile);
        fseek ( InFile , 0 , SEEK_SET );
        //long CurPos;

        char* cdummy = new char [256];

        //Read the header lines
        for(unsigned short us=0; us<NumBlankTransitionLines; us++){
            fgets(cdummy, 255, InFile);
        }

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return;
        }

        float fRadius;
        float fMomentum;
        float fReWf;
        float fImWf;
        float fReAsWf;
        float fImAsWf;
        float fCatsWf;
        float fDummy;
        float fPhaseShift;

        unsigned RadBin;
        unsigned LastRadBin;
        unsigned MomBin;
        //!---Iteration over all events---
        while(!feof(InFile)){
            fgets(cdummy, 255, InFile);
//printf("I am here - %i!\n", strlen(cdummy));
            if(strlen(cdummy)!=105) continue;
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                   &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsWf,&fDummy,&fPhaseShift);
            LastRadBin = RadBin;
            RadBin = Kitty.GetBin(fRadius,RadBins[0],NumRadBins+1);
            //we filled up all rad bins
            if(RadBin<=LastRadBin){
                for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
                    RadBinLoaded[uRad]=false;
                }

            }
            if(RadBin>=NumRadBins) continue;
            if(RadBinLoaded[RadBin]){
                printf("\033[1;33mWARNING:\033[0m Trying to load RadBin #%u twice!\n", RadBin);
                continue;
            }
            MomBin = Kitty.GetMomBin(fMomentum);
            if(MomBin>=NumMomBins) continue;

            //skip momentum bins that we do not have in the CATS object
            if(fabs(fMomentum-Kitty.GetMomentum(MomBin))>0.1){
//printf("fMomentum=%f <> Kitty.GetMomentum(MomBin)=%f\n",fMomentum,Kitty.GetMomentum(MomBin));

                continue;
            }
//printf("!!!!\n");
            switch(uFile){
            case f1S0:
                WaveFunctionU[0][MomBin][0][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][0] = fPhaseShift*3.14159/180.;
//printf("MomBin=%u; RadBin=%u; WaveFunctionU[0][MomBin][0][0][RadBin]=%f\n",MomBin,RadBin,WaveFunctionU[0][MomBin][0][0][RadBin]);
                break;
            case f1P1:
                WaveFunctionU[0][MomBin][0][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][1] = fPhaseShift*3.14159/180.;
                break;
            case f3S1:
                WaveFunctionU[0][MomBin][1][0][RadBin] = fCatsWf*fRadius;
                WaveFunctionU[0][MomBin][2][0][RadBin] = fCatsWf*fRadius;
                WaveFunctionU[0][MomBin][3][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][0] = fPhaseShift*3.14159/180.;
                PhaseShifts[0][MomBin][2][0] = fPhaseShift*3.14159/180.;
                PhaseShifts[0][MomBin][3][0] = fPhaseShift*3.14159/180.;
                break;
            case f3P0:
                WaveFunctionU[0][MomBin][1][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][1] = fPhaseShift*3.14159/180.;
                break;
            case f3P1:
                WaveFunctionU[0][MomBin][2][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][2][1] = fPhaseShift*3.14159/180.;
                break;
            case f3P2:
                WaveFunctionU[0][MomBin][3][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][3][1] = fPhaseShift*3.14159/180.;
                break;
            default:
                printf("\033[1;31mERROR:\033[0m SWITCH CHANNELS AND PWS\n");
                break;
            }
            RadBinLoaded[RadBin] = true;
        }
        delete [] cdummy;
        fclose(InFile);
    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
    }
    delete [] InputFileName;

    delete [] RadBinLoaded;
//printf("hhhhhhhhhhhh\n");
}

void InitHaidenbauerKaonMinus(const char* InputFolder, CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins){

    unsigned NumMomBins=Kitty.GetNumMomBins();
    //double kMin = Kitty.GetMomBinLowEdge(0);
    //double kMax = Kitty.GetMomBinUpEdge(NumMomBins-1);

    const unsigned short NumChannels = 4;
    const unsigned short NumPwPerCh = 2;
    const unsigned short NumFiles = 6;
    enum HaideFiles {fS01, fS11, fP01, fP03, fP11, fP13};

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
    }
//KaonMinus_SIonly
//KaonMinus10MeV

    strcpy(InputFileName[fS01],InputFolder);
    strcpy(InputFileName[fS11],InputFolder);
    strcpy(InputFileName[fP01],InputFolder);
    strcpy(InputFileName[fP03],InputFolder);
    strcpy(InputFileName[fP11],InputFolder);
    strcpy(InputFileName[fP13],InputFolder);

    strcat(InputFileName[fS01], "ws01.dat");
    strcat(InputFileName[fS11], "ws11.dat");
    strcat(InputFileName[fP01], "wp01.dat");
    strcat(InputFileName[fP03], "wp03.daa");
    strcat(InputFileName[fP11], "wp11.dat");
    strcat(InputFileName[fP13], "wp13.dat");

    const unsigned NumBlankTransitionLines = 1;
    const double MinRad = 0.1;
    const double RadStep = 0.1;
    const double MaxRad = 10.0;
    NumRadBins = round((MaxRad-MinRad)/RadStep) + 1;
    bool* RadBinLoaded = new bool [NumRadBins+1];

    WaveFunctionU[0] = new double*** [NumMomBins];
    PhaseShifts[0] = new double** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        WaveFunctionU[0][uMomBin] = new double** [NumChannels];
        PhaseShifts[0][uMomBin] = new double* [NumChannels];
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            WaveFunctionU[0][uMomBin][usCh] = new double* [NumPwPerCh];
            PhaseShifts[0][uMomBin][usCh] = new double [NumPwPerCh];
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                WaveFunctionU[0][uMomBin][usCh][usPw] = new double [NumRadBins];
                PhaseShifts[0][uMomBin][usCh][usPw] = 0;
            }
        }
    }

    RadBins[0] = new double [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        RadBins[0][uRad] = MinRad-RadStep*0.5 + RadStep*double(uRad);
        RadBinLoaded[uRad] = false;
    }
//printf("Hey!\n");
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return;
        }
//printf("uFile=%u\n",uFile);
        fseek ( InFile , 0 , SEEK_END );
        //long EndPos;
        //EndPos = ftell (InFile);
        fseek ( InFile , 0 , SEEK_SET );
        //long CurPos;

        char* cdummy = new char [256];

        //Read the header lines
        for(unsigned short us=0; us<NumBlankTransitionLines; us++){
            fgets(cdummy, 255, InFile);
        }

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return;
        }

        float fRadius;
        float fMomentum;
        float fReWf;
        float fImWf;
        float fReAsWf;
        float fImAsWf;
        float fCatsReWf;
        float fCatsImWf;
        float fCatsWf;
        float fPhaseShift;

        unsigned RadBin;
        unsigned LastRadBin;
        unsigned MomBin;
        //!---Iteration over all events---
        while(!feof(InFile)){
            fgets(cdummy, 255, InFile);
//printf("I am here - %i!\n", strlen(cdummy));
            if(strlen(cdummy)<105 || strlen(cdummy)>106) continue;
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                   &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsReWf,&fCatsImWf,&fPhaseShift);
            //sscanf(cdummy, " %f %f %f %f",
            //       &fRadius,&fMomentum,&fCatsReWf,&fCatsImWf);
            fCatsWf = sqrt(fCatsReWf*fCatsReWf+fCatsImWf*fCatsImWf);
            //if(fCatsReWf<0) fCatsWf=-fCatsWf;
            //fCatsWf = fReWf;
            LastRadBin = RadBin;
            RadBin = Kitty.GetBin(fRadius,RadBins[0],NumRadBins+1);
            //we filled up all rad bins
            if(RadBin<=LastRadBin){
                for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
                    RadBinLoaded[uRad]=false;
                }

            }
            if(RadBin>=NumRadBins) continue;
            if(RadBinLoaded[RadBin]){
                printf("\033[1;33mWARNING:\033[0m Trying to load RadBin #%u twice!\n", RadBin);
                continue;
            }
            MomBin = Kitty.GetMomBin(fMomentum);
            if(MomBin>=NumMomBins) continue;

            //skip momentum bins that we do not have in the CATS object
            if(fabs(fMomentum-Kitty.GetMomentum(MomBin))>0.1){
printf("fMomentum=%f <> Kitty.GetMomentum(MomBin)=%f\n",fMomentum,Kitty.GetMomentum(MomBin));

                continue;
            }
//printf("!!!!\n");
            switch(uFile){
            case fS01://Spin 0, Isospin 0, s-wave
                WaveFunctionU[0][MomBin][0][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][0] = fPhaseShift*3.14159/180.;
                WaveFunctionU[0][MomBin][1][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][0] = fPhaseShift*3.14159/180.;
                break;
            case fS11:
                WaveFunctionU[0][MomBin][2][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][2][0] = fPhaseShift*3.14159/180.;
                WaveFunctionU[0][MomBin][3][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][3][0] = fPhaseShift*3.14159/180.;
                break;
            case fP01:
                WaveFunctionU[0][MomBin][0][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][1] = fPhaseShift*3.14159/180.;
                break;
            case fP03:
                WaveFunctionU[0][MomBin][1][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][1] = fPhaseShift*3.14159/180.;
                break;
            case fP11:
                WaveFunctionU[0][MomBin][2][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][2][1] = fPhaseShift*3.14159/180.;
                break;
            case fP13:
                WaveFunctionU[0][MomBin][3][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][3][1] = fPhaseShift*3.14159/180.;
                break;
            default:
                printf("\033[1;31mERROR:\033[0m SWITCH CHANNELS AND PWS\n");
                break;
            }
            RadBinLoaded[RadBin] = true;
        }
        delete [] cdummy;
        fclose(InFile);
    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
    }
    delete [] InputFileName;

    delete [] RadBinLoaded;
//printf("hhhhhhhhhhhh\n");
}

//see email in may 2018, this is isospin averaged
//basically we have two channels: s1+p1 and s1+p3
void InitHaidenbauerKaonMinus_ver2(const char* InputFolder, CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins){

    unsigned NumMomBins=Kitty.GetNumMomBins();
    //double kMin = Kitty.GetMomBinLowEdge(0);
    //double kMax = Kitty.GetMomBinUpEdge(NumMomBins-1);

    const unsigned short NumChannels = 4;
    const unsigned short NumPwPerCh = 2;
    const unsigned short NumFiles = 3;
    enum HaideFiles {fS1, fP1, fP3};

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
    }

    strcpy(InputFileName[fS1],InputFolder);
    strcpy(InputFileName[fP1],InputFolder);
    strcpy(InputFileName[fP3],InputFolder);

    strcat(InputFileName[fS1], "ws1.dat");
    strcat(InputFileName[fP1], "wp1.dat");
    strcat(InputFileName[fP3], "wp3.dat");

    const unsigned NumBlankTransitionLines = 1;
    const double MinRad = 0.1;
    const double RadStep = 0.1;
    const double MaxRad = 10.0;
    NumRadBins = round((MaxRad-MinRad)/RadStep) + 1;
    bool* RadBinLoaded = new bool [NumRadBins+1];

    WaveFunctionU[0] = new double*** [NumMomBins];
    PhaseShifts[0] = new double** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        WaveFunctionU[0][uMomBin] = new double** [NumChannels];
        PhaseShifts[0][uMomBin] = new double* [NumChannels];
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            WaveFunctionU[0][uMomBin][usCh] = new double* [NumPwPerCh];
            PhaseShifts[0][uMomBin][usCh] = new double [NumPwPerCh];
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                WaveFunctionU[0][uMomBin][usCh][usPw] = new double [NumRadBins];
                PhaseShifts[0][uMomBin][usCh][usPw] = 0;
            }
        }
    }

    RadBins[0] = new double [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        RadBins[0][uRad] = MinRad-RadStep*0.5 + RadStep*double(uRad);
        RadBinLoaded[uRad] = false;
    }
//printf("Hey!\n");
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return;
        }
//printf("uFile=%u\n",uFile);
        fseek ( InFile , 0 , SEEK_END );
        //long EndPos;
        //EndPos = ftell (InFile);
        fseek ( InFile , 0 , SEEK_SET );
        //long CurPos;

        char* cdummy = new char [256];

        //Read the header lines
        for(unsigned short us=0; us<NumBlankTransitionLines; us++){
            fgets(cdummy, 255, InFile);
        }

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return;
        }

        float fRadius;
        float fMomentum;
        float fReWf;
        float fImWf;
        float fReAsWf;
        float fImAsWf;
        float fCatsReWf;
        float fCatsImWf;
        float fCatsWf;
        float fPhaseShift;

        unsigned RadBin;
        unsigned LastRadBin;
        unsigned MomBin;
        //!---Iteration over all events---
        while(!feof(InFile)){
            fgets(cdummy, 255, InFile);
//printf("I am here - %i!\n", strlen(cdummy));
            if(strlen(cdummy)<105 || strlen(cdummy)>106) continue;
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                   &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsReWf,&fCatsImWf,&fPhaseShift);
            //sscanf(cdummy, " %f %f %f %f",
            //       &fRadius,&fMomentum,&fCatsReWf,&fCatsImWf);
            fCatsWf = sqrt(fCatsReWf*fCatsReWf+fCatsImWf*fCatsImWf);
            //if(fCatsReWf<0) fCatsWf=-fCatsWf;
            //fCatsWf = fReWf;
            LastRadBin = RadBin;
            RadBin = Kitty.GetBin(fRadius,RadBins[0],NumRadBins+1);
            //we filled up all rad bins
            if(RadBin<=LastRadBin){
                for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
                    RadBinLoaded[uRad]=false;
                }

            }
            if(RadBin>=NumRadBins) continue;
            if(RadBinLoaded[RadBin]){
                printf("\033[1;33mWARNING:\033[0m Trying to load RadBin #%u twice!\n", RadBin);
                continue;
            }
            MomBin = Kitty.GetMomBin(fMomentum);
            if(MomBin>=NumMomBins) continue;

            //skip momentum bins that we do not have in the CATS object
            if(fabs(fMomentum-Kitty.GetMomentum(MomBin))>0.1){
printf("fMomentum=%f <> Kitty.GetMomentum(MomBin)=%f\n",fMomentum,Kitty.GetMomentum(MomBin));

                continue;
            }
//printf("!!!!\n");
            switch(uFile){
            case fS1://Spin 0, Isospin 0, s-wave
                WaveFunctionU[0][MomBin][0][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][0] = fPhaseShift*3.14159/180.;
                WaveFunctionU[0][MomBin][1][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][0] = fPhaseShift*3.14159/180.;
                break;
            case fP1:
                WaveFunctionU[0][MomBin][0][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][1] = fPhaseShift*3.14159/180.;
                break;
            case fP3:
                WaveFunctionU[0][MomBin][1][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][1] = fPhaseShift*3.14159/180.;
                break;
            default:
                printf("\033[1;31mERROR:\033[0m SWITCH CHANNELS AND PWS\n");
                break;
            }
            RadBinLoaded[RadBin] = true;
        }
        delete [] cdummy;
        fclose(InFile);
    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
    }
    delete [] InputFileName;

    delete [] RadBinLoaded;
//printf("hhhhhhhhhhhh\n");
}


void InitHaidenbauerKaonPlus(const char* InputFolder, CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins){

    unsigned NumMomBins=Kitty.GetNumMomBins();
    //double kMin = Kitty.GetMomBinLowEdge(0);
    //double kMax = Kitty.GetMomBinUpEdge(NumMomBins-1);

    const unsigned short NumChannels = 4;//they are two, but I simply keep this for mem management
    const unsigned short NumPwPerCh = 2;
    const unsigned short NumFiles = 3;
    enum HaideFiles {fS11, fP11, fP13};

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
    }

    strcpy(InputFileName[fS11],InputFolder);
    strcpy(InputFileName[fP11],InputFolder);
    strcpy(InputFileName[fP13],InputFolder);

    strcat(InputFileName[fS11], "ws11p.dat");
    strcat(InputFileName[fP11], "wp11p.dat");
    strcat(InputFileName[fP13], "wp13p.dat");

    const unsigned NumBlankTransitionLines = 1;
    const double MinRad = 0.1;
    const double RadStep = 0.1;
    const double MaxRad = 10.0;
    NumRadBins = round((MaxRad-MinRad)/RadStep) + 1;
    bool* RadBinLoaded = new bool [NumRadBins+1];

    WaveFunctionU[0] = new double*** [NumMomBins];
    PhaseShifts[0] = new double** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        WaveFunctionU[0][uMomBin] = new double** [NumChannels];
        PhaseShifts[0][uMomBin] = new double* [NumChannels];
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            WaveFunctionU[0][uMomBin][usCh] = new double* [NumPwPerCh];
            PhaseShifts[0][uMomBin][usCh] = new double [NumPwPerCh];
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                WaveFunctionU[0][uMomBin][usCh][usPw] = new double [NumRadBins];
                PhaseShifts[0][uMomBin][usCh][usPw] = 0;
            }
        }
    }

    RadBins[0] = new double [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        RadBins[0][uRad] = MinRad-RadStep*0.5 + RadStep*double(uRad);
        RadBinLoaded[uRad] = false;
    }
printf("Hey!\n");
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return;
        }
printf("uFile=%u\n",uFile);
        fseek ( InFile , 0 , SEEK_END );
        //long EndPos;
        //EndPos = ftell (InFile);
        fseek ( InFile , 0 , SEEK_SET );
        //long CurPos;

        char* cdummy = new char [256];

        //Read the header lines
        for(unsigned short us=0; us<NumBlankTransitionLines; us++){
            fgets(cdummy, 255, InFile);
        }

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return;
        }

        float fRadius;
        float fMomentum;
        float fReWf;
        float fImWf;
        float fReAsWf;
        float fImAsWf;
        float fCatsReWf;
        float fCatsImWf;
        float fCatsWf;
        float fPhaseShift;

        unsigned RadBin;
        unsigned LastRadBin;
        unsigned MomBin;
        //!---Iteration over all events---
        while(!feof(InFile)){
            fgets(cdummy, 255, InFile);
//printf("I am here - %i!\n", strlen(cdummy));
            if(strlen(cdummy)<105 || strlen(cdummy)>106) continue;
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                   &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsReWf,&fCatsImWf,&fPhaseShift);
            fCatsWf = sqrt(fCatsReWf*fCatsReWf+fCatsImWf*fCatsImWf);
            LastRadBin = RadBin;
            RadBin = Kitty.GetBin(fRadius,RadBins[0],NumRadBins+1);
            //we filled up all rad bins
            if(RadBin<=LastRadBin){
                for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
                    RadBinLoaded[uRad]=false;
                }

            }
            if(RadBin>=NumRadBins) continue;
            if(RadBinLoaded[RadBin]){
                printf("\033[1;33mWARNING:\033[0m Trying to load RadBin #%u twice!\n", RadBin);
                continue;
            }
            MomBin = Kitty.GetMomBin(fMomentum);
            if(MomBin>=NumMomBins) continue;

            //skip momentum bins that we do not have in the CATS object
            if(fabs(fMomentum-Kitty.GetMomentum(MomBin))>0.1){
printf("fMomentum=%f <> Kitty.GetMomentum(MomBin)=%f\n",fMomentum,Kitty.GetMomentum(MomBin));

                continue;
            }
//printf("!!!!\n");
            switch(uFile){
            case fS11:
                WaveFunctionU[0][MomBin][0][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][0] = fPhaseShift*3.14159/180.;
                WaveFunctionU[0][MomBin][1][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][0] = fPhaseShift*3.14159/180.;
                break;
            case fP11:
                WaveFunctionU[0][MomBin][0][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][1] = fPhaseShift*3.14159/180.;
                break;
            case fP13:
                WaveFunctionU[0][MomBin][1][1][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][1][1] = fPhaseShift*3.14159/180.;
                break;
            default:
                printf("\033[1;31mERROR:\033[0m SWITCH CHANNELS AND PWS\n");
                break;
            }
            RadBinLoaded[RadBin] = true;
        }
        delete [] cdummy;
        fclose(InFile);
    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
    }
    delete [] InputFileName;

    delete [] RadBinLoaded;
//printf("hhhhhhhhhhhh\n");
}


//[0][MomBin][uCh][uPw]
//TYPE = 0 => take the psi^2
//TYPE = 1 => the dirty trick with flipping the sign based on the sign of the real part
void InitTetsuoKaonMinus(const char* InputFolder, CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins, unsigned& NumRadBins, const int& TYPE){

    unsigned NumMomBins=Kitty.GetNumMomBins();
    //double kMin = Kitty.GetMomBinLowEdge(0);
    //double kMax = Kitty.GetMomBinUpEdge(NumMomBins-1);

    const unsigned short NumChannels = 4;
    const unsigned short NumPwPerCh = 2;
    const unsigned short NumFiles = 1;

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
    }
//wf_full
//wf_strong

    strcpy(InputFileName[0],InputFolder);

    strcat(InputFileName[0], "wf_full.dat");

    const unsigned NumBlankTransitionLines = 1;
    const double MinRad = 0.1;
    const double RadStep = 0.1;
    const double MaxRad = 10.0;
    NumRadBins = round((MaxRad-MinRad)/RadStep) + 1;
    bool* RadBinLoaded = new bool [NumRadBins+1];

    WaveFunctionU[0] = new double*** [NumMomBins];
    PhaseShifts[0] = new double** [NumMomBins];
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        WaveFunctionU[0][uMomBin] = new double** [NumChannels];
        PhaseShifts[0][uMomBin] = new double* [NumChannels];
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            WaveFunctionU[0][uMomBin][usCh] = new double* [NumPwPerCh];
            PhaseShifts[0][uMomBin][usCh] = new double [NumPwPerCh];
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                WaveFunctionU[0][uMomBin][usCh][usPw] = new double [NumRadBins];
                PhaseShifts[0][uMomBin][usCh][usPw] = 0;
            }
        }
    }

    RadBins[0] = new double [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        RadBins[0][uRad] = MinRad-RadStep*0.5 + RadStep*double(uRad);
        RadBinLoaded[uRad] = false;
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return;
        }
        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );
        char* cdummy = new char [256];
        //Read the header lines
        for(unsigned short us=0; us<NumBlankTransitionLines; us++){
            fgets(cdummy, 255, InFile);
        }
        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return;
        }
        float fRadius;
        float fMomentum;
        //float fReWf;
        //float fImWf;
        //float fReAsWf;
        //float fImAsWf;
        float fCatsReWf;
        float fCatsImWf;
        float fCatsWf;
        float fPhaseShift=0;

        unsigned RadBin;
        unsigned LastRadBin;
        unsigned MomBin;
        //!---Iteration over all events---
        while(!feof(InFile)){
            fgets(cdummy, 255, InFile);
            if(strlen(cdummy)<105 || strlen(cdummy)>106) continue;
            sscanf(cdummy, " %f %f %f %f",
                   &fRadius,&fMomentum,&fCatsReWf,&fCatsImWf);
            fCatsWf = sqrt(fCatsReWf*fCatsReWf+fCatsImWf*fCatsImWf);
            if(TYPE==1){
                if(fCatsReWf<0) fCatsWf=-fCatsWf;
                //fCatsWf = fReWf;
            }

            LastRadBin = RadBin;
            RadBin = Kitty.GetBin(fRadius,RadBins[0],NumRadBins+1);
            //we filled up all rad bins
            if(RadBin<=LastRadBin){
                for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
                    RadBinLoaded[uRad]=false;
                }
            }
            if(RadBin>=NumRadBins) continue;
            if(RadBinLoaded[RadBin]){
                printf("\033[1;33mWARNING:\033[0m Trying to load RadBin #%u twice!\n", RadBin);
                continue;
            }
            MomBin = Kitty.GetMomBin(fMomentum);
            if(MomBin>=NumMomBins) continue;

            //skip momentum bins that we do not have in the CATS object
            if(fabs(fMomentum-Kitty.GetMomentum(MomBin))>0.1){
                continue;
            }
            switch(uFile){
            case 0://Spin 0, Isospin 0, s-wave
                WaveFunctionU[0][MomBin][0][0][RadBin] = fCatsWf*fRadius;
                PhaseShifts[0][MomBin][0][0] = fPhaseShift*3.14159/180.;
                break;
            default:
                printf("\033[1;31mERROR:\033[0m SWITCH CHANNELS AND PWS\n");
                break;
            }
            RadBinLoaded[RadBin] = true;
        }
        delete [] cdummy;
        fclose(InFile);
    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
    }
    delete [] InputFileName;

    delete [] RadBinLoaded;
//printf("hhhhhhhhhhhh\n");
}

void CleanHaidenbauer(CATS& Kitty, double***** WaveFunctionU, double**** PhaseShifts, double** RadBins){
    unsigned NumMomBins = Kitty.GetNumMomBins();
    const unsigned short NumChannels = 4;
    const unsigned short NumPwPerCh = 2;
    //const unsigned short NumFiles = 6;

    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        for(unsigned usCh=0; usCh<NumChannels; usCh++){
            for(unsigned usPw=0; usPw<NumPwPerCh; usPw++){
                delete [] WaveFunctionU[0][uMomBin][usCh][usPw];
            }
            delete [] WaveFunctionU[0][uMomBin][usCh];
            delete [] PhaseShifts[0][uMomBin][usCh];
        }
        delete [] WaveFunctionU[0][uMomBin];
        delete [] PhaseShifts[0][uMomBin];
    }

    delete [] WaveFunctionU[0];
    delete [] PhaseShifts[0];
    delete [] RadBins[0];
}












