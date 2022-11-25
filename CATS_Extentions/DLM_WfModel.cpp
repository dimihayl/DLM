#include <stdio.h>
#include "string.h"

#include "DLM_WfModel.h"
#include "CATS.h"
#include "DLM_Histo.h"

#include <unistd.h>

const complex<float>fi(0,1);
const double Pi(3.141592653589793);

//TYPE = 00 is LO
//TYPE = 01 is LO+Coupling
//TYPE = 03 is LO+Coupling, s+d waves
//TYPE = 10 is NLO (s-wave)
//TYPE = 11 is NLO+Coupling
//TYPE = 12 is NLO (s and p-waves)
//TYPE = 13 is NLO+Coupling, s+d waves
//for TYPE==01/11 we make it so, that Kitty has 4 channels
//ch0 and ch1 are the LN->LN 1S0 and 3S1 (as per usual)
//ch2 and ch3 are the SN->LN in the S-wave (1S0 and 3S1), however we ONLY add to the total correlation function the s-wave
//to obtain the total correlation we need to add 1/4*ch0+3/4*ch1+1/4*ch2+3/4*ch3. As the ch2 and ch3 are just "corrections" to
//the primary channels (0 and 1) it is completely fine that our weights do not sum up to 1
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF){
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    if(CUTOFF!=500&&CUTOFF!=600){
        printf("Problem with the CUTOFF in Init_pL_Haidenbauer\n");
        return NULL;
    }

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    //unsigned short NumFiles;
    unsigned NumMomBins = 0;
    unsigned short NumFiles = 10; // (to cover all s,p,d, waves)
    bool* TakeThisFile = new bool [NumFiles];

    if(TYPE==0){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 40;
    }
    else if(TYPE==1){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 43;
    }
    //the channels here are:
    //[0] is 1S0 LN->LN (S wave)
    //[1] is 1S0 SN->LN (S wave)
    //[2] is LN->LN (S/D wave)
    //[3] is SN->LN (S/D wave)
    //[4] is LN->LN (D/S->S/D wave)
    //[5] is SN->LN (D/S->S/D wave)
    //the X/Y for the 2-5 means that in the s-wave we have the X-state, in d-wave the Y-state
    //(the channels 2-5 have both s and d waves, ch 0,1 have only s wave)
    else if(TYPE==3){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 6;
        NumPwPerCh = 3;
        NumFiles = 3;
        //NumMomBins = 43;
    }
//! p-waves not yet implemented!
    else if(TYPE==10&&CUTOFF==600){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.2;
        RadiusMaximum = 16.1;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 40;
    }
    else if(TYPE==10&&CUTOFF==500){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 40;
    }
    else if(TYPE==11){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 43;
    }
    else if(TYPE==12&&CUTOFF==600){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.2;
        RadiusMaximum = 16.1;
        NumChannels = 4;
        NumPwPerCh = 2;
        NumFiles = 6;
        //NumMomBins = 40;
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pL_Haidenbauer\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    if(TYPE==12){
        Kitty.SetNumPW(0,NumPwPerCh);
        Kitty.SetChannelWeight(0,3./12.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,NumPwPerCh);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,NumPwPerCh);
        Kitty.SetChannelWeight(2,3./12.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,NumPwPerCh);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);
    }
    else{
        for(unsigned uCh=0; uCh<NumChannels; uCh++){
            Kitty.SetNumPW(uCh,NumPwPerCh);
            Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
            Kitty.SetSpin(uCh,uCh%2==0?0:1);
        }
    }

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    const double Mass_p = 938.272;
    const double Mass_L = 1115.683;
    Kitty.SetRedMass((Mass_p*Mass_L)/(Mass_p+Mass_L));
    if(TYPE%10==1){
        Kitty.SetOnlyNumericalPw(2,true);
        Kitty.SetOnlyNumericalPw(3,true);
    }

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum HaideFiles {f1S0, f3S1, f1P1, f3P0, f3P1, f3P2, f1D2, f3D1, f3D2, f3D3};
    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    if(TYPE==10 && CUTOFF==600){
        strcat(InputFileName[f1S0], "N1s0.data");
        //strcat(InputFileName[f1P1], "N1p1.data");
        strcat(InputFileName[f3S1], "N3s1.data");
        //strcat(InputFileName[f3P0], "N3p0.data");
        //strcat(InputFileName[f3P1], "N3p1.data");
        //strcat(InputFileName[f3P2], "N3p2.data");
    }
    else if(TYPE==10 && CUTOFF==500){
        strcat(InputFileName[f1S0], "Y1s05.data");
        strcat(InputFileName[f3S1], "Y3s15.data");
    }
    else if(TYPE==0 && CUTOFF==600){
        strcat(InputFileName[f1S0], "W1s06.data");
        strcat(InputFileName[f3S1], "W3s16.data");
    }
    else if(TYPE==11 && CUTOFF==600){
        strcat(InputFileName[f1S0], "fort.13");
        strcat(InputFileName[f3S1], "fort.14");
    }
    else if(TYPE==1 && CUTOFF==600){
        strcat(InputFileName[f1S0], "fort.13");
        strcat(InputFileName[f3S1], "fort.14");
        //strcat(InputFileName[f3S1], "fort.14");
    }
    else if(TYPE==12 && CUTOFF==600){
        strcat(InputFileName[f1S0], "N1s0.data");
        strcat(InputFileName[f1P1], "N1p1.data");
        strcat(InputFileName[f3S1], "N3s1.data");
        strcat(InputFileName[f3P0], "N3p0.data");
        strcat(InputFileName[f3P1], "N3p1.data");
        strcat(InputFileName[f3P2], "N3p2.data");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pL_Haidenbauer\n");
    }

    FILE *InFile;
    InFile = fopen(InputFileName[0], "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[0]);
        return NULL;
    }
    int CurrentLine=0;
    //unsigned WhichMomBin;

    const unsigned MaxNumMomBins = 256;
    double* MomentumBins = new double [MaxNumMomBins+1];
    double* Momentum = new double [MaxNumMomBins];

    char* cdummy = new char [512];
    float fMomentum;

    MomentumBins[0] = 0;
    while(!feof(InFile)){
        if(!fgets(cdummy, 511, InFile)) continue;
        if((CurrentLine)%int(NumRadBins)==0){
            sscanf(cdummy, "%f",&fMomentum);
            Momentum[NumMomBins] = fMomentum;
            if(NumMomBins){
                //set the bin range in between the last two bin centers
                MomentumBins[NumMomBins] = 0.5*(Momentum[NumMomBins]+Momentum[NumMomBins-1]);
            }
            NumMomBins++;
        }
        CurrentLine++;
    }
    fclose(InFile);
    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    MomentumBins[NumMomBins] = 2.*Momentum[NumMomBins-1]-MomentumBins[NumMomBins-1];

    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins,MomentumBins);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins,MomentumBins);
            HistoPS[uHist][uPw].Initialize();
        }
    }
//printf("NumDLM_Histos=%u\n",NumDLM_Histos);
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
//if(uFile>2) break;
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }
        //we have the p-waves only for NLO with cutoff of 600 MeV
        if( (uFile!=f1S0 && uFile!=f3S1 && (TYPE!=10 || CUTOFF!=600)) &&
            (uFile>=6||TYPE!=12||CUTOFF!=600)
           ) {printf("Possible bug in Init_pL_Haidenbauer\n"); continue;}

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //unsigned LastRadBin;

        float fRadius;
        float fReWf;
        float fImWf;
        float fReAsWf;
        float fImAsWf;
        float fCatsWf;
        float fDummy;
        float fPhaseShift;

        float fReWf_LNLN_SS;
        float fImWf_LNLN_SS;
        float fReWf_SNLN_SS;
        float fImWf_SNLN_SS;

        //for the time being we do not use those
        float fReWf_LNLN_DS;
        float fImWf_LNLN_DS;
        float fReWf_SNLN_DS;
        float fImWf_SNLN_DS;

        //MomentumBins[0] = 0;

        int RadBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);

                RadBin++;
                WhichMomBin++;
                continue;
            }

            if(TYPE%10==0||TYPE==12){
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                       &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsWf,&fDummy,&fPhaseShift);
            }
            else if(TYPE%10==1){
                if(uFile==0){
                    sscanf(cdummy, " %f %f %f %f %f %f",
                       &fRadius,&fDummy,&fReWf_LNLN_SS,&fImWf_LNLN_SS,&fReWf_SNLN_SS,&fImWf_SNLN_SS);
                        fReWf_LNLN_DS=0;
                        fImWf_LNLN_DS=0;
                        fReWf_SNLN_DS=0;
                        fImWf_SNLN_DS=0;
                }
                else if(uFile==1){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",
                       &fRadius,&fDummy,&fReWf_LNLN_SS,&fImWf_LNLN_SS,&fReWf_SNLN_SS,&fImWf_SNLN_SS,
                                        &fReWf_LNLN_DS,&fImWf_LNLN_DS,&fReWf_SNLN_DS,&fImWf_SNLN_DS);
                }
                else{
                    printf("Oh man... big bug in Init_pL_Haidenbauer.\n");
                }
            }
            else{
                printf("WTF from Init_pL_Haidenbauer\n");
            }

            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            //we fill up the radius bins with and offset of 1, due to the special zeroth bin
            WhichBin[1] = RadBin+1;

            if(TYPE==0){
                HistoWF[uFile][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                HistoPS[uFile][0].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==1 && CUTOFF==600){
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_LNLN_SS+fi*fImWf_LNLN_SS)*fRadius);
                HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_SNLN_SS+fi*fImWf_SNLN_SS)*fRadius);

                HistoPS[uFile][0].SetBinContent(WhichBin,0);
                HistoPS[uFile+2][0].SetBinContent(WhichBin,0);
            }
            else if(TYPE==10 && CUTOFF==600){
                HistoWF[uFile][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                //HistoWF[uFile].SetBinContent(WhichBin,(fReWf+i*fImWf)*fRadius);
                HistoPS[uFile][0].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==10 && CUTOFF==500){
                HistoWF[uFile][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                HistoPS[uFile][0].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==11 && CUTOFF==600){
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_LNLN_SS+fi*fImWf_LNLN_SS)*fRadius);
                HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_SNLN_SS+fi*fImWf_SNLN_SS)*fRadius);

                //HistoWF[uFile].SetBinContent(WhichBin,(fReWf_LNLN_DS+i*fImWf_LNLN_DS)*fRadius);
                //HistoWF[uFile+2].SetBinContent(WhichBin,(fReWf_SNLN_DS+i*fImWf_SNLN_DS)*fRadius);
//&fReWf_LNLN_DS,&fImWf_LNLN_DS,&fReWf_SNLN_DS,&fImWf_SNLN_DS

//if(fMomentum>280 && fMomentum<290){
//printf("uFile=%u; fMomentum=%.1f; fReWf=%.4f(%.4f); fImWf=%.4f(%.4f)\n",uFile,fMomentum,
//       fReWf_LNLN_SS,fReWf_SNLN_SS,fImWf_LNLN_SS,fImWf_SNLN_SS);
//}
                HistoPS[uFile][0].SetBinContent(WhichBin,0);
                HistoPS[uFile+2][0].SetBinContent(WhichBin,0);
            }
            //enum HaideFiles {f1S0, f3S1, f1P1, f3P0, f3P1, f3P2};
            else if(TYPE==12 && CUTOFF==600){
                if(uFile==f1S0){
                    HistoWF[0][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                    HistoPS[0][0].SetBinContent(WhichBin,fPhaseShift);
                }
                else if(uFile==f3S1){
                    HistoWF[1][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                    HistoPS[1][0].SetBinContent(WhichBin,fPhaseShift);

                    HistoWF[2][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                    HistoPS[2][0].SetBinContent(WhichBin,fPhaseShift);

                    HistoWF[3][0].SetBinContent(WhichBin,fCatsWf*fRadius);
                    HistoPS[3][0].SetBinContent(WhichBin,fPhaseShift);
                }
                else{
                    HistoWF[uFile-2][1].SetBinContent(WhichBin,fCatsWf*fRadius);
                    HistoPS[uFile-2][1].SetBinContent(WhichBin,fPhaseShift);
                }
            }
            else{

            }
            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }
        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TakeThisFile;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF){
    return Init_pL_Haidenbauer(InputFolder,*Kitty,TYPE,CUTOFF);
}





//if the CutOff is 601, than we switch to FineTuneMode (for CutOff 600), where TYPE is 4 digits, ABCD, where:
//  AB is a variation of the 1S0. So far we have:
//    00 default (f_s = 2.91 fm)
//    01 default+CSB
//    02 f_s = 2.2 fm
//    03 f_s = 2.55 fm
//  CD is a variation of the 3S1. So far we have:
//    00 default (f_t = 1.41 fm)
//    01 default+CSB
//    02 f_t = 1.3 fm
//    03 f_t = 1.25 fm
//    04 f_t = 1.2 fm
//    05 f_t = 1.15 fm
//    06 f_t = 1.1 fm
//    (10,11,12): reseved potentiall for 1.45, 1.50, 1.55
//    13 f_t = 1.6 fm
//I also did CutOff 602, since the previous scheme was not getting all covered. Here:
//  S=0
//    00 default (2.91)
//    01 CSB
//    02-09 reserve for other special cases
//    the logic from here on: >=10 we go up in scat. len.
//                            <=99 we go down in scat. len.
//    so far:
//    99 2.80
//    98 2.65
//    97 2.55
//    96 2.45
//    95 2.30
//    94 2.20
//  S=1 (same logic as above)
//    00 default (1.41)
//    01 CSB
//    99 1.35
//    98 1.30
//    97 1.25
//    96 1.20
//    95 1.15
//    94 1.10
//    93 1.05
//    92 1.00

//TYPE==0 is NLO13, TYPE==1 is NLO19, TYPE==-1 is LO 13 TYPE==-2 is LO19
//if TYPE is two digti: XY, X is the setting for 1S0, Y for 3S1,
//where both can be {1,2,3} = {ORIGINAL, CSB, FineTune}
//CSB = charge symmetry breaking included, FineTine is where the 3S1 scattreing ampitude is reduced to 1.3
//To make the tests for the files provided by JH in April 2022, we need:
//131600: test the effect of reducing the f0 in 3S1
//221600: test the effect of CSB
//231600: test the effect of CSB only in 1S0, where 3S1 is the new fine tune
//211600: test the effect of CSB only in 1S0
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2019(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF){
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    const unsigned short NumFiles = 6; // (to cover all s,p,d, waves)
    unsigned* NumMomBins = new unsigned [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        NumMomBins[uFile] = 0;
    }

    //fP1 is both 1P1 and 3P1
    enum HaideFiles {f1S0, f3S1, fP1, f3P0, f3P2, f3D1};
    bool FileAvailable[NumFiles];
    bool InludePwaves=true;
    for(unsigned uFile=0; uFile<NumFiles; uFile++) FileAvailable[uFile]=true;
    //we have the p-waves only for NLO at a cutoff of 600 (and 601,602)
    //if( (CUTOFF/10!=60||TYPE%10!=1) && CUTOFF!=601 ){
    if( CUTOFF/10!=60||TYPE%10!=1){
        FileAvailable[fP1] = false;
        FileAvailable[f3P0] = false;
        FileAvailable[f3P2] = false;
        InludePwaves = false;
    }

    //LIST OF CHANNELS:
    //main channels:
    //0: 1S0+1P1
    //1: 3S1+3P0+3D1
    //2: 3S1+3P1+3D1
    //3: 3S1+3P2+3D1
    //4: 3S1+3P0
    //5: 3S1+3P1
    //6: 3S1+3P2
    //coupled channels:
    //7: 1S0 SN(s) -> LN(s)
    //8: 3S1 SN(s) -> LN(s)
    //9: 3S1 LN(d) -> LN(s)
    //10: 3S1 SN(d) -> LN(s)
    //11: 3P0 SN(p) -> LN(p)
    //12: 3P2 SN(p) -> LN(p)
    //13: 3D1 SN(d) -> LN(d)
    //14: 3D1 LN(s) -> LN(d)
    //15: 3D1 SN(s) -> LN(d)
    const unsigned short NumChannels = 16;

    RadiusStepSize = 0.02;
    RadiusMinimum = 0.02;
    RadiusMaximum = 10.;

    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,2);
    Kitty.SetNumPW(1,3);
    Kitty.SetNumPW(2,3);
    Kitty.SetNumPW(3,3);
    Kitty.SetNumPW(4,2);
    Kitty.SetNumPW(5,2);
    Kitty.SetNumPW(6,2);

    //for the coupled channels we only take the integral of the wave function
    //hence it does not matter if we set the waves to s,p,d or whatever
    //thus for simplification, everything is saved in the s-wave
    Kitty.SetNumPW(7,1);
    Kitty.SetNumPW(8,1);
    Kitty.SetNumPW(9,1);
    Kitty.SetNumPW(10,1);
    Kitty.SetNumPW(11,1);
    Kitty.SetNumPW(12,1);
    Kitty.SetNumPW(13,1);
    Kitty.SetNumPW(14,1);
    Kitty.SetNumPW(15,1);

    Kitty.SetChannelWeight(0,1./4.);
    Kitty.SetChannelWeight(1,1./60.);
    Kitty.SetChannelWeight(2,1./20.);
    Kitty.SetChannelWeight(3,1./12.);
    Kitty.SetChannelWeight(4,1./15.);
    Kitty.SetChannelWeight(5,1./5.);
    Kitty.SetChannelWeight(6,1./3.);

    Kitty.SetChannelWeight(7,1./4.);
    Kitty.SetChannelWeight(8,3./4.);
    Kitty.SetChannelWeight(9,3./4.);
    Kitty.SetChannelWeight(10,3./4.);
    Kitty.SetChannelWeight(11,1./12.*InludePwaves);
    Kitty.SetChannelWeight(12,5./12.*InludePwaves);
    Kitty.SetChannelWeight(13,3./20.);
    Kitty.SetChannelWeight(14,3./20.);
    Kitty.SetChannelWeight(15,3./20.);

    for(unsigned short usCh=0; usCh<NumChannels; usCh++){Kitty.SetSpin(usCh,1);}
    Kitty.SetSpin(0,0);
    Kitty.SetSpin(7,0);

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    const double Mass_p = 938.272;
    const double Mass_L = 1115.683;
    Kitty.SetRedMass((Mass_p*Mass_L)/(Mass_p+Mass_L));
//printf("1\n");
    for(unsigned short usCh=7; usCh<NumChannels; usCh++){Kitty.SetOnlyNumericalPw(usCh,true);}




    //LIST OF CHANNELS:
    //main channels:
    //0: 1S0+1P1
    //1: 3S1+3P0+3D1
    //2: 3S1+3P1+3D1
    //3: 3S1+3P2+3D1
    //4: 3S1+3P0
    //5: 3S1+3P1
    //6: 3S1+3P2
    //coupled channels:
    //7: 1S0 SN(s) -> LN(s)
    //8: 3S1 SN(s) -> LN(s)
    //9: 3S1 LN(d) -> LN(s)
    //10: 3S1 SN(d) -> LN(s)
    //11: 3P0 SN(p) -> LN(p)//
    //12: 3P2 SN(p) -> LN(p)//
    //13: 3D1 SN(d) -> LN(d)
    //14: 3D1 LN(s) -> LN(d)
    //15: 3D1 SN(s) -> LN(d)
    unsigned short** WhichFile = new unsigned short* [NumChannels];
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){WhichFile[usCh] = new unsigned short [Kitty.GetNumPW(usCh)];}
    WhichFile[0][0] = f1S0;
    WhichFile[0][1] = fP1;
    WhichFile[1][0] = f3S1;
    WhichFile[1][1] = f3P0;
    WhichFile[1][2] = f3D1;
    WhichFile[2][0] = f3S1;
    WhichFile[2][1] = fP1;
    WhichFile[2][2] = f3D1;
    WhichFile[3][0] = f3S1;
    WhichFile[3][1] = f3P2;
    WhichFile[3][2] = f3D1;
    WhichFile[4][0] = f3S1;
    WhichFile[4][1] = f3P0;
    WhichFile[5][0] = f3S1;
    WhichFile[5][1] = fP1;
    WhichFile[6][0] = f3S1;
    WhichFile[6][1] = f3P2;
    WhichFile[7][0] = f1S0;
    WhichFile[8][0] = f3S1;
    WhichFile[9][0] = f3S1;
    WhichFile[10][0] = f3S1;
    WhichFile[11][0] = f3P0;
    WhichFile[12][0] = f3P2;
    WhichFile[13][0] = f3D1;
    WhichFile[14][0] = f3D1;
    WhichFile[15][0] = f3D1;

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;
//printf("2\n");

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
        strcat(InputFileName[uFile],"Chiral");
    }

    if(TYPE==0){
        char strCutOff[16];
        sprintf(strCutOff,"13_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++) strcat(InputFileName[uFile], strCutOff);

        strcat(InputFileName[f1S0], "NLO1s0.data");
        strcat(InputFileName[f3S1], "NLO3s1.data");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    //the special case for fine tunes of the 600 cutoff (for NLO19)
    else if(CUTOFF==601){
      //  AB is a variation of the 1S0. So far we have:
      //    00 default (f_s = 2.91 fm)
      //    01 default+CSB
      //    02 f_s = 2.2 fm
      //    03 f_s = 2.55 fm
      //  CD is a variation of the 1S0. So far we have:
      //    00 default (f_t = 1.41 fm)
      //    01 default+CSB
      //    02 f_t = 1.3 fm
      //    03 f_t = 1.25 fm
      //    04 f_t = 1.2 fm
      //    05 f_t = 1.15 fm
      //    06 f_t = 1.10 fm
      //    (10,11,12): reseved potentiall for 1.45, 1.50, 1.55
      //    13 f_t = 1.6 fm
      //printf("TYPE = %i\n",TYPE);
      if(TYPE%10!=1){
        printf("Possible BUG in Init_pL_Haidenbauer2019, or a wrong set up of the pL potential\n");
        printf("The CutOff 601 is only allowed to be used with the NLO19 potential\n");
      }
      int ABCD = TYPE/10;
      for(unsigned uFile=0; uFile<NumFiles; uFile++){
        //1S0 ---------------------------------------------------------
        if(uFile==f1S0&&ABCD/100==1){//this is set differently in this case
          strcat(InputFileName[uFile], "19_CSB/");
          strcat(InputFileName[uFile], "csbN191s0.600");
        }
        else if(uFile==f1S0&&ABCD/100==2){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p20.600");
        }
        else if(uFile==f1S0&&ABCD/100==3){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p55.600");
        }
        // ------------------------------------------------------------
        //3S1 ---------------------------------------------------------
        else if(uFile==f3S1&&ABCD%100==1){
          strcat(InputFileName[uFile], "19_CSB/");
          strcat(InputFileName[uFile], "csbN193s1.600");
        }
        else if(uFile==f3S1&&ABCD%100==2){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p30.600");
        }
        else if(uFile==f3S1&&ABCD%100==3){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p25.600");
        }
        else if(uFile==f3S1&&ABCD%100==4){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p20.600");
        }
        else if(uFile==f3S1&&ABCD%100==5){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p15.600");
        }
        else if(uFile==f3S1&&ABCD%100==6){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p10.600");
        }
        else if(uFile==f3S1&&ABCD%100==13){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p60.600");
        }
        // ------------------------------------------------------------
        else{
          char strCutOff[16];
          sprintf(strCutOff,"19_%i/",600);
          strcat(InputFileName[uFile], strCutOff);
          switch (uFile) {
            case f1S0: strcat(InputFileName[uFile], "NLO1s0.data"); break;
            case f3S1: strcat(InputFileName[uFile], "NLO3s1.data"); break;
            case fP1: strcat(InputFileName[uFile], "NLOPU.data"); break;
            case f3P0: strcat(InputFileName[uFile], "NLO3P0.data"); break;
            case f3P2: strcat(InputFileName[uFile], "NLO3P2.data"); break;
            case f3D1: strcat(InputFileName[uFile], "NLO3d1.data"); break;
            default: printf("Something bad has happend in Init_pL_Haidenbauer2019\n"); break;
          }
        }
      }
    }//CUTOFF 601

    //the special case for fine tunes of the 600 cutoff (for NLO19)
    else if(CUTOFF==602){
      //this is correct, TYPE 1 is 1160*
      if(TYPE%10!=1){
        printf("Possible BUG in Init_pL_Haidenbauer2019, or a wrong set up of the pL potential\n");
        printf("The CutOff 602 is only allowed to be used with the NLO19 potential\n");
        //printf("TYPE = %i\n",TYPE);
      }
      int ABCD = TYPE/10;
      for(unsigned uFile=0; uFile<NumFiles; uFile++){
        //1S0 ---------------------------------------------------------
        if(uFile==f1S0&&ABCD/100==1){//this is set differently in this case
          strcat(InputFileName[uFile], "19_CSB/");
          strcat(InputFileName[uFile], "csbN191s0.600");
        }
        else if(uFile==f1S0&&ABCD/100==99){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p80.600");
        }
        else if(uFile==f1S0&&ABCD/100==98){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p65.600");
        }
        else if(uFile==f1S0&&ABCD/100==97){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p55.600");
        }
        else if(uFile==f1S0&&ABCD/100==96){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p45.600");
        }
        else if(uFile==f1S0&&ABCD/100==95){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p30.600");
        }
        else if(uFile==f1S0&&ABCD/100==94){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_1s0_2p20.600");
        }
        // ------------------------------------------------------------
        //3S1 ---------------------------------------------------------
        else if(uFile==f3S1&&ABCD%100==1){
          strcat(InputFileName[uFile], "19_CSB/");
          strcat(InputFileName[uFile], "csbN193s1.600");
        }
        else if(uFile==f3S1&&ABCD%100==99){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p35.600");
        }
        else if(uFile==f3S1&&ABCD%100==98){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p30.600");
        }
        else if(uFile==f3S1&&ABCD%100==97){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p25.600");
        }
        else if(uFile==f3S1&&ABCD%100==96){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p20.600");
        }
        else if(uFile==f3S1&&ABCD%100==95){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p15.600");
        }
        else if(uFile==f3S1&&ABCD%100==94){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p10.600");
        }
        else if(uFile==f3S1&&ABCD%100==93){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p05.600");
        }
        else if(uFile==f3S1&&ABCD%100==92){
          strcat(InputFileName[uFile], "19_FineTunes/");
          strcat(InputFileName[uFile], "NLO19_3s1_1p00.600");
        }
        // ------------------------------------------------------------
        else{
          char strCutOff[16];
          sprintf(strCutOff,"19_%i/",600);
          strcat(InputFileName[uFile], strCutOff);
          switch (uFile) {
            case f1S0: strcat(InputFileName[uFile], "NLO1s0.data"); break;
            case f3S1: strcat(InputFileName[uFile], "NLO3s1.data"); break;
            case fP1: strcat(InputFileName[uFile], "NLOPU.data"); break;
            case f3P0: strcat(InputFileName[uFile], "NLO3P0.data"); break;
            case f3P2: strcat(InputFileName[uFile], "NLO3P2.data"); break;
            case f3D1: strcat(InputFileName[uFile], "NLO3d1.data"); break;
            default: printf("Something bad has happend in Init_pL_Haidenbauer2019\n"); break;
          }
        }
      }
    }//CUTOFF 602


    //the 1S0 for the fine tune is not changed
    else if(TYPE==1||TYPE==11||TYPE==31){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++) strcat(InputFileName[uFile], strCutOff);

        strcat(InputFileName[f1S0], "NLO1s0.data");
        strcat(InputFileName[f3S1], "NLO3s1.data");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    else if(TYPE==-1){
        char strCutOff[16];
        sprintf(strCutOff,"13_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++) strcat(InputFileName[uFile], strCutOff);

        strcat(InputFileName[f1S0], "LO1s0.data");
        strcat(InputFileName[f3S1], "LO3s1.data");
        strcat(InputFileName[fP1], "LOPU.data");
        strcat(InputFileName[f3P0], "LO3P0.data");
        strcat(InputFileName[f3P2], "LO3P2.data");
        strcat(InputFileName[f3D1], "LO3d1.data");
    }
    else if(TYPE==-2){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++) strcat(InputFileName[uFile], strCutOff);

        strcat(InputFileName[f1S0], "LO1s0.data");
        strcat(InputFileName[f3S1], "LO3s1.data");
        strcat(InputFileName[fP1], "LOPU.data");
        strcat(InputFileName[f3P0], "LO3P0.data");
        strcat(InputFileName[f3P2], "LO3P2.data");
        strcat(InputFileName[f3D1], "LO3d1.data");
    }
    //the 1S0 for the fine tune is not changed
    else if(TYPE==12||TYPE==32){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
          if(uFile==f3S1){//this is set differently in this case
            strcat(InputFileName[uFile], "19_CSB/");
          }
          else strcat(InputFileName[uFile], strCutOff);
        }

        strcat(InputFileName[f1S0], "NLO1s0.data");
        strcat(InputFileName[f3S1], "csbN193s1.600");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    else if(TYPE==21){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
          if(uFile==f1S0){//this is set differently in this case
            strcat(InputFileName[uFile], "19_CSB/");
          }
          else strcat(InputFileName[uFile], strCutOff);
        }

        strcat(InputFileName[f1S0], "csbN191s0.600");
        strcat(InputFileName[f3S1], "NLO3s1.data");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    else if(TYPE==22){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
          if(uFile==f1S0||uFile==f3S1){//this is set differently in this case
            strcat(InputFileName[uFile], "19_CSB/");
          }
          else strcat(InputFileName[uFile], strCutOff);
        }

        strcat(InputFileName[f1S0], "csbN191s0.600");
        strcat(InputFileName[f3S1], "csbN193s1.600");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    //the 1S0 for the fine tune is not changed
    else if(TYPE==13||TYPE==33){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
          if(uFile==f3S1){//this is set differently in this case
            strcat(InputFileName[f3S1], "19_FineTunes/");
          }
          else strcat(InputFileName[uFile], strCutOff);
        }

        strcat(InputFileName[f1S0], "NLO1s0.data");
        strcat(InputFileName[f3S1], "WN193s1.600");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    else if(TYPE==23){
        char strCutOff[16];
        sprintf(strCutOff,"19_%i/",CUTOFF);
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
          if(uFile==f3S1){//this is set differently in this case
            strcat(InputFileName[f3S1], "19_FineTunes/");
          }
          else if(uFile==f1S0){
            strcat(InputFileName[uFile], "19_CSB/");
          }
          else strcat(InputFileName[uFile], strCutOff);
        }

        strcat(InputFileName[f1S0], "csbN191s0.600");
        strcat(InputFileName[f3S1], "WN193s1.600");
        strcat(InputFileName[fP1], "NLOPU.data");
        strcat(InputFileName[f3P0], "NLO3P0.data");
        strcat(InputFileName[f3P2], "NLO3P2.data");
        strcat(InputFileName[f3D1], "NLO3d1.data");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pL_Haidenbauer2019\n");
        return NULL;
    }
//printf("3\n");
    FILE *InFile;
    int CurrentLine=0;
    //unsigned WhichMomBin;

    const unsigned MaxNumMomBins = 256;
    double** MomentumBins = new double* [NumFiles];
    double** Momentum = new double* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        MomentumBins[uFile] = new double [MaxNumMomBins+1];
        Momentum[uFile] = new double [MaxNumMomBins];
    }

    char* cdummy = new char [512];
    float fMomentum;
//printf("4\n");
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        if(!FileAvailable[uFile]) continue;
//printf("uf=%u --> %s\n",uFile,InputFileName[uFile]);
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
printf("1-- FileAvailable[%u]=%u\n",uFile,FileAvailable[uFile]);
            return NULL;
        }

        CurrentLine=0;
        MomentumBins[uFile][0] = 0;
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if((CurrentLine)%int(NumRadBins)==0){
                sscanf(cdummy, "%f",&fMomentum);
                Momentum[uFile][NumMomBins[uFile]] = fMomentum;
                if(NumMomBins[uFile]){
                    //set the bin range in between the last two bin centers
                    MomentumBins[uFile][NumMomBins[uFile]] = 0.5*(Momentum[uFile][NumMomBins[uFile]]+Momentum[uFile][NumMomBins[uFile]-1]);
                }
                NumMomBins[uFile]++;
            }
            CurrentLine++;
        }
        fclose(InFile);
        //set the upper edge of the last bin, where we just add the bin width of the last bin
        //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
        MomentumBins[uFile][NumMomBins[uFile]] = 2.*Momentum[uFile][NumMomBins[uFile]-1]-MomentumBins[uFile][NumMomBins[uFile]-1];
    }
//printf("5\n");
    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uHist)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uHist); uPw++){
            if(!FileAvailable[WhichFile[uHist][uPw]]) continue;
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins[WhichFile[uHist][uPw]],MomentumBins[WhichFile[uHist][uPw]]);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }
//printf("6\n");
    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uHist)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uHist); uPw++){
            if(!FileAvailable[WhichFile[uHist][uPw]]) continue;
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins[WhichFile[uHist][uPw]],MomentumBins[WhichFile[uHist][uPw]]);
            HistoPS[uHist][uPw].Initialize();
        }
    }
//printf("7\n");

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        if(!FileAvailable[uFile]) continue;
//if(uFile>2) break;
//printf("uf=%u\n",uFile);
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
printf("2-- FileAvailable[%u]=%u\n",uFile,FileAvailable[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //unsigned LastRadBin;

        float fDummy;
        float fRadius;
        float fReWf1;
        float fImWf1;
        float fReWf2;
        float fImWf2;
        float fReWf3;
        float fImWf3;
        float fReWf4;
        float fImWf4;

        //MomentumBins[0] = 0;

        int RadBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(WhichMomBin>=int(NumMomBins[uFile])){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins[uFile]);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);

                RadBin++;
                WhichMomBin++;
                continue;
            }


            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            //we fill up the radius bins with and offset of 1, due to the special zeroth bin
            WhichBin[1] = RadBin+1;


            if(uFile==f1S0){
                //Wf1 is LN->LN
                //Wf2 is SN->LN
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2);
                //0: 1S0+1P1
                HistoWF[0][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[0][0].SetBinContent(WhichBin,0);
                //7: 1S0 SN(s) -> LN(s)
                HistoWF[7][0].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[7][0].SetBinContent(WhichBin,0);
            }
            else if(uFile==fP1){
                //Wf1 is 1P1
                //Wf2 is 3P1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2);
                //0: 1S0+1P1
                HistoWF[0][1].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[0][1].SetBinContent(WhichBin,0);
                //2: 3S1+3P1+3D1
                HistoWF[2][1].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[2][1].SetBinContent(WhichBin,0);
                ////5: 3S1+3P1
                HistoWF[5][1].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[5][1].SetBinContent(WhichBin,0);
            }
            else if(uFile==f3S1){
                //Wf1 is LN->LN (s wave)
                //Wf2 is SN->LN (s wave)
                //Wf3 is LN->LN (d->s wave)
                //Wf4 is SN->LN (d->s wave)
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2,&fReWf3,&fImWf3,&fReWf4,&fImWf4);
                //1: 3S1+3P0+3D1
                HistoWF[1][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[1][0].SetBinContent(WhichBin,0);
                //2: 3S1+3P1+3D1
                HistoWF[2][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[2][0].SetBinContent(WhichBin,0);
                //3: 3S1+3P2+3D1
                HistoWF[3][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[3][0].SetBinContent(WhichBin,0);
                //4: 3S1+3P0
                HistoWF[4][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[4][0].SetBinContent(WhichBin,0);
                //5: 3S1+3P1
                HistoWF[5][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[5][0].SetBinContent(WhichBin,0);
                //6: 3S1+3P2
                HistoWF[6][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[6][0].SetBinContent(WhichBin,0);
                //8: 3S1 SN(s) -> LN(s)
                HistoWF[8][0].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[8][0].SetBinContent(WhichBin,0);
                //9: 3S1 LN(d) -> LN(s)
                HistoWF[9][0].SetBinContent(WhichBin,(fReWf3+fi*fImWf3)*fRadius);
                HistoPS[9][0].SetBinContent(WhichBin,0);
                //10: 3S1 SN(d) -> LN(s)
                HistoWF[10][0].SetBinContent(WhichBin,(fReWf4+fi*fImWf4)*fRadius);
                HistoPS[10][0].SetBinContent(WhichBin,0);
            }
            else if(uFile==f3P0){
                //Wf1 is LN->LN
                //Wf2 is SN->LN
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2);
                //1: 3S1+3P0+3D1
                HistoWF[1][1].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[1][1].SetBinContent(WhichBin,0);
                //4: 3S1+3P0
                HistoWF[4][1].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[4][1].SetBinContent(WhichBin,0);
                //11: 3P0 SN(p) -> LN(p)
                HistoWF[11][0].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[11][0].SetBinContent(WhichBin,0);
            }
            else if(uFile==f3P2){
                //Wf1 is LN->LN (p wave)
                //Wf2 is SN->LN (p wave)
                //Wf3 is (??? but small, so we ignore it)
                //Wf4 is (??? but small, so we ignore it)
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2,&fReWf3,&fImWf3,&fReWf4,&fImWf4);
                //3: 3S1+3P2+3D1
                HistoWF[3][1].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[3][1].SetBinContent(WhichBin,0);
                //6: 3S1+3P2
                HistoWF[6][1].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[6][1].SetBinContent(WhichBin,0);
                //12: 3P2 SN(p) -> LN(p)
                HistoWF[12][0].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[12][0].SetBinContent(WhichBin,0);
            }
            else if(uFile==f3D1){
                //Wf1 is LN->LN (d wave) WRONG
                //Wf2 is SN->LN (d wave) WRONG
                //Wf3 is LN->LN (s->d wave) WRONG
                //Wf4 is SN->LN (s->d wave) WRONG
                //in reality it is:
                //Wf1 is LN->LN (s->d wave)
                //Wf2 is SN->LN (s->d wave)
                //Wf3 is LN->LN (d wave)
                //Wf4 is SN->LN (d wave)
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf1,&fImWf1,&fReWf2,&fImWf2,&fReWf3,&fImWf3,&fReWf4,&fImWf4);
                //1: 3S1+3P0+3D1
                HistoWF[1][2].SetBinContent(WhichBin,(fReWf3+fi*fImWf3)*fRadius);
                HistoPS[1][2].SetBinContent(WhichBin,0);
                //2: 3S1+3P1+3D1
                HistoWF[2][2].SetBinContent(WhichBin,(fReWf3+fi*fImWf3)*fRadius);
                HistoPS[2][2].SetBinContent(WhichBin,0);
                //3: 3S1+3P2+3D1
                HistoWF[3][2].SetBinContent(WhichBin,(fReWf3+fi*fImWf3)*fRadius);
                HistoPS[3][2].SetBinContent(WhichBin,0);
                //13: 3D1 SN(d) -> LN(d)
                HistoWF[13][0].SetBinContent(WhichBin,(fReWf4+fi*fImWf4)*fRadius);
                HistoPS[13][0].SetBinContent(WhichBin,0);
                //14: 3D1 LN(s) -> LN(d)
                HistoWF[14][0].SetBinContent(WhichBin,(fReWf1+fi*fImWf1)*fRadius);
                HistoPS[14][0].SetBinContent(WhichBin,0);
                //15: 3D1 SN(s) -> LN(d)
                HistoWF[14][0].SetBinContent(WhichBin,(fReWf2+fi*fImWf2)*fRadius);
                HistoPS[14][0].SetBinContent(WhichBin,0);
            }
            else{
                printf("WTF from Init_pL_Haidenbauer2019\n");
            }

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }

        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins[uFile])){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins[uFile]);
        }

    }//uFile

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] MomentumBins[uFile];
        delete [] Momentum[uFile];
    }
    delete [] NumMomBins;
    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        delete [] WhichFile[usCh];
    }
    delete [] WhichFile;
    return Histo;

}
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer2019(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF){
    return Init_pL_Haidenbauer2019(InputFolder,*Kitty,TYPE,CUTOFF);
}



//TYPE = 0 (no options)
//we have 6 channels in total:
//ch0/ch1 are the S0p -> S0p, 1S0/3S1
//ch2/ch3 are the Lp -> S0p, 1S0/3S1
//ch4/ch5 are the S+n -> S0p, 1S0/3S1
//we add the even (1S0) channels with weight 1/4 and the odd (3S1) with 3/4
//the fact that we go beyond a total weight of 1 is not an issue here, as the channels 2-5
//are just "corrections" to the primary channels (0 and 1)
DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE){
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    unsigned short NumFiles;
    unsigned NumMomBins = 0;

    if(TYPE==0){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 6;
        NumPwPerCh = 1;
        NumFiles = 2;
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pSigma0_Haidenbauer\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3212);
    const double Mass_p = 938.272;
    const double Mass_S0 = 1192.642;
    Kitty.SetRedMass((Mass_p*Mass_S0)/(Mass_p+Mass_S0));

    Kitty.SetOnlyNumericalPw(2,true);
    Kitty.SetOnlyNumericalPw(3,true);
    Kitty.SetOnlyNumericalPw(4,true);
    Kitty.SetOnlyNumericalPw(5,true);

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum HaideFiles {f1S0, f3S1};
    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    if(TYPE==0){
        strcat(InputFileName[f1S0], "S0p1s0.data");
        strcat(InputFileName[f3S1], "S0p3s1.data");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pSigma0_Haidenbauer\n");
    }

    FILE *InFile;
    InFile = fopen(InputFileName[0], "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[0]);
        return NULL;
    }
    int CurrentLine=0;

    const unsigned MaxNumMomBins = 256;
    double* MomentumBins = new double [MaxNumMomBins+1];
    double* Momentum = new double [MaxNumMomBins];

    char* cdummy = new char [512];
    float fMomentum;

    MomentumBins[0] = 0;
    while(!feof(InFile)){
        if(!fgets(cdummy, 511, InFile)) continue;
        if((CurrentLine)%int(NumRadBins)==0){
            sscanf(cdummy, "%f",&fMomentum);
            Momentum[NumMomBins] = fMomentum;
            if(NumMomBins){
                //set the bin range in between the last two bin centers
                MomentumBins[NumMomBins] = 0.5*(Momentum[NumMomBins]+Momentum[NumMomBins-1]);
            }
            NumMomBins++;
        }
        CurrentLine++;
    }
    fclose(InFile);
    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    MomentumBins[NumMomBins] = 2.*Momentum[NumMomBins-1]-MomentumBins[NumMomBins-1];

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
        HistoWF[uCh][0].SetUp(2);
        HistoWF[uCh][0].SetUp(0,NumMomBins,MomentumBins);
        HistoWF[uCh][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh][0].Initialize();
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
        HistoPS[uCh][0].SetUp(1);
        HistoPS[uCh][0].SetUp(0,NumMomBins,MomentumBins);
        HistoPS[uCh][0].Initialize();
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        float fRadius;
        float fDummy;

        float fReWf_Lp_S0p;
        float fImWf_Lp_S0p;
        float fReWf_S0p_S0p;
        float fImWf_S0p_S0p;
        float fReWf_Sn_S0p;
        float fImWf_Sn_S0p;

        int RadBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);

                RadBin++;
                WhichMomBin++;
                continue;
            }

            if(uFile==f1S0){
                sscanf(cdummy, " %f %f %f %f %f %f %f %f",
                    &fRadius,&fDummy,&fReWf_Lp_S0p,&fImWf_Lp_S0p,&fReWf_S0p_S0p,&fImWf_S0p_S0p,&fReWf_Sn_S0p,&fImWf_Sn_S0p);
            }
            else if(uFile==f3S1){
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                    &fRadius,&fDummy,&fReWf_Lp_S0p,&fImWf_Lp_S0p,&fReWf_S0p_S0p,&fImWf_S0p_S0p,&fReWf_Sn_S0p,&fImWf_Sn_S0p,
                                    &fDummy,&fDummy,&fDummy,&fDummy,&fDummy,&fDummy);
            }
            else{
                printf("Oh man... big bug in Init_pSigma0_Haidenbauer.\n");
            }

            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            //we fill up the radius bins with and offset of 1, due to the special zeroth bin
            WhichBin[1] = RadBin+1;

            HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_S0p_S0p+fi*fImWf_S0p_S0p)*fRadius);
            HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_Lp_S0p+fi*fImWf_Lp_S0p)*fRadius);
            HistoWF[uFile+4][0].SetBinContent(WhichBin,(fReWf_Sn_S0p+fi*fImWf_Sn_S0p)*fRadius);

            HistoPS[uFile][0].SetBinContent(WhichBin,0);
            HistoPS[uFile+2][0].SetBinContent(WhichBin,0);
            HistoPS[uFile+4][0].SetBinContent(WhichBin,0);

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }
        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile
    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pSigma0_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE){
    return Init_pSigma0_Haidenbauer(InputFolder,*Kitty,TYPE);
}

//TYPE = 0 (no options)
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF){
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    if(CUTOFF!=500&&CUTOFF!=550&&CUTOFF!=600&&CUTOFF!=650){
        printf("Problem with the CUTOFF in Init_pXi_Haidenbauer\n");
        return NULL;
    }

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    unsigned short NumFiles;
    unsigned NumMomBins = 0;

    if(TYPE==0){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 4;
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pXi_Haidenbauer\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.125:0.375);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);
    const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;
    Kitty.SetRedMass((Mass_p*Mass_Xim)/(Mass_p+Mass_Xim));

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum HaideFiles {fI0S0, fI0S1, fI1S0, fI1S1};
    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    char* cdummy = new char [512];

    if(TYPE==0){
        sprintf(cdummy,"%i.data",CUTOFF);
        strcat(InputFileName[fI0S0], "I0S0");
        strcat(InputFileName[fI0S0], cdummy);
        strcat(InputFileName[fI0S1], "I0S1");
        strcat(InputFileName[fI0S1], cdummy);
        strcat(InputFileName[fI1S0], "I1S0");
        strcat(InputFileName[fI1S0], cdummy);
        strcat(InputFileName[fI1S1], "I1S1");
        strcat(InputFileName[fI1S1], cdummy);
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pL_Haidenbauer\n");
    }

    FILE *InFile;
    InFile = fopen(InputFileName[0], "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[0]);
        return NULL;
    }
    int CurrentLine=0;
    //unsigned WhichMomBin;

    const unsigned MaxNumMomBins = 256;
    double* MomentumBins = new double [MaxNumMomBins+1];
    double* Momentum = new double [MaxNumMomBins];

    float fMomentum;

    MomentumBins[0] = 0;
    while(!feof(InFile)){
        if(!fgets(cdummy, 511, InFile)) continue;
        if((CurrentLine)%int(NumRadBins)==0){
            sscanf(cdummy, "%f",&fMomentum);
            Momentum[NumMomBins] = fMomentum;
            if(NumMomBins){
                //set the bin range in between the last two bin centers
                MomentumBins[NumMomBins] = 0.5*(Momentum[NumMomBins]+Momentum[NumMomBins-1]);
            }
            NumMomBins++;
        }
        CurrentLine++;
    }
    fclose(InFile);
    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    MomentumBins[NumMomBins] = 2.*Momentum[NumMomBins-1]-MomentumBins[NumMomBins-1];

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
        HistoWF[uCh][0].SetUp(2);
        HistoWF[uCh][0].SetUp(0,NumMomBins,MomentumBins);
        HistoWF[uCh][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh][0].Initialize();
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
        HistoPS[uCh][0].SetUp(1);
        HistoPS[uCh][0].SetUp(0,NumMomBins,MomentumBins);
        HistoPS[uCh][0].Initialize();
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        float fRadius;
        float fReWf;
        float fImWf;
        float fDummy;

        int RadBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);

                RadBin++;
                WhichMomBin++;
                continue;
            }

            if(TYPE==0){
                if(uFile==fI0S0){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",&fRadius,&fMomentum,&fDummy,&fDummy,&fReWf,&fImWf,&fDummy,&fDummy,&fDummy,&fDummy);
                }
                else if(uFile==fI0S1){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f",&fRadius,&fMomentum,&fDummy,&fDummy,&fReWf,&fImWf,&fDummy,&fDummy);
                }
                else if(uFile==fI1S0){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",&fRadius,&fMomentum,&fDummy,&fDummy,&fDummy,&fDummy,&fReWf,&fImWf,&fDummy,&fDummy);
                }
                else{
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f",&fRadius,&fMomentum,&fDummy,&fDummy,&fDummy,&fDummy,&fReWf,&fImWf);
                }
            }
            else{
                printf("WTF from Init_pXi_Haidenbauer\n");
            }

            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            //we fill up the radius bins with and offset of 1, due to the special zeroth bin
            WhichBin[1] = RadBin+1;

            if(TYPE==0){
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf+fi*fImWf)*fRadius);
                HistoPS[uFile][0].SetBinContent(WhichBin,0);
            }
            else{
                printf("WTF from Init_pXi_Haidenbauer\n");
            }

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }
        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile
    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pXi_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF){
    return Init_pXi_Haidenbauer(InputFolder,*Kitty,TYPE,CUTOFF);
}


DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS& Kitty, const int& TYPE){

    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    unsigned short NumFiles;
    const unsigned NumSpinChannels = 2;
    //unsigned NumMomBins = 0;

    unsigned NumDummyLines=0;
    unsigned NumBufferLines=0;
    unsigned NumCol=0;
    unsigned NumRow=0;
    unsigned NumFuckedUpEntriesAtTheEnd=0;

    unsigned NumTomMomBins=0;
    double* TomMomentum=NULL;
    double* TomMomentumBins=NULL;
    double* pLab=NULL;

    if(TYPE==0){
        RadiusStepSize = 0.05;
        RadiusMinimum = 0.05;
        RadiusMaximum = 14.95;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 2;
        NumDummyLines=32;
        NumBufferLines=2;
        NumCol=10;
        NumRow=30;
        NumFuckedUpEntriesAtTheEnd=1;
        NumTomMomBins = 21;
        pLab = new double [NumTomMomBins];
        pLab[0] = 5;
        pLab[1] = 50;
        for(unsigned uBin=2; uBin<NumTomMomBins; uBin++){
            pLab[uBin] = pLab[uBin-1]+50;
        }
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pXi_ESC16_IS\n");
        return NULL;
    }

    const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;

    TomMomentum = new double [NumTomMomBins];
    TomMomentumBins = new double [NumTomMomBins+1];
    for(unsigned uBin=0; uBin<NumTomMomBins; uBin++){
        TomMomentum[uBin] = sqrt(pow(Mass_p*pLab[uBin],2)/(pow(Mass_p,2)+pow(Mass_Xim,2)+2*Mass_p*sqrt(pow(Mass_Xim,2)+pow(pLab[uBin],2))));
    }
    TomMomentumBins[0]=0;
    for(unsigned short uBin=0; uBin<NumTomMomBins-1; uBin++){
        TomMomentumBins[uBin+1] = (TomMomentum[uBin]+TomMomentum[uBin+1])*0.5;
    }
    TomMomentumBins[NumTomMomBins] = TomMomentum[NumTomMomBins-1]+(TomMomentum[NumTomMomBins-1]-TomMomentumBins[NumTomMomBins-1]);

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);

    Kitty.SetRedMass((Mass_p*Mass_Xim)/(Mass_p+Mass_Xim));

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum TomFiles {fI0, fI1};
    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[512];
        strcpy(InputFileName[uFile],InputFolder);
    }

    if(TYPE==0){
        strcat(InputFileName[fI0], "xnwave.ich=0.ic=0.l=0.s=0-1");
        strcat(InputFileName[fI1], "xnwave.ich=1.ic=0.l=0.s=0-1");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pXi_ESC16_IS\n");
    }

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
        HistoWF[uCh][0].SetUp(2);
        HistoWF[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoWF[uCh][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoWF[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
        HistoPS[uCh][0].SetUp(1);
        HistoPS[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoPS[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoPS[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    char* cdummy = new char [512];

    unsigned WhichSpin=0;
    unsigned WhichIsospin=0;

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        WhichIsospin = uFile;
        int WhichMomBin=-1;
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //Read the header lines
        for(unsigned short us=0; us<NumDummyLines; us++){
            if(!fgets(cdummy, 511, InFile)) continue;
        }

        float Values[NumCol];
        float fRadius;

        //int RadBin=-1;
        unsigned WhichBin[2];
        WhichBin[0]=0;
        WhichBin[1]=0;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(WhichBin[0]>=NumTomMomBins) {WhichBin[0]=0; WhichSpin=1;}
            for(unsigned uBuffer=0; uBuffer<NumBufferLines; uBuffer++){
                if(!fgets(cdummy, 511, InFile)) continue;
            }

            for(unsigned uRow=0; uRow<NumRow; uRow++){
                if(!fgets(cdummy, 511, InFile)) continue;
                if(uRow==NumRow-1&&NumFuckedUpEntriesAtTheEnd==1){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8]);
                }
                else if(NumFuckedUpEntriesAtTheEnd){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8],&Values[9]);
                }
                else{
                    printf("Issues in Init_pXi_ESC16_IS\n");
                }

                for(unsigned uCol=0; uCol<NumCol; uCol++){
                    WhichBin[1] = uRow*NumCol+uCol;
                    fRadius = HistoWF[WhichIsospin*NumSpinChannels+WhichSpin]->GetBinCenter(1,WhichBin[1]);
//printf("fRadius=%f\n",fRadius);
                    HistoWF[WhichIsospin*NumSpinChannels+WhichSpin]->SetBinContent(WhichBin,Values[uCol]*fRadius);
                }
            }
            HistoPS[WhichIsospin*NumSpinChannels+WhichSpin]->SetBinContent(WhichBin,0);
            WhichBin[0]++;

        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumTomMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumTomMomBins (%u vs %u)\n",WhichMomBin,NumTomMomBins);
        }
    }//uFile
    delete [] RadBinLoaded;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TomMomentum;
    delete [] TomMomentumBins;
    delete [] pLab;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS* Kitty, const int& TYPE){
    return Init_pXi_ESC16_IS(InputFolder,*Kitty,TYPE);
}

//TYPE==0, some old WF,
//TYPE==1, function from 2nd April 2019
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS& Kitty, const int& TYPE){

    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    unsigned short NumFiles;
    //const unsigned NumSpinChannels = 2;
    //unsigned NumMomBins = 0;

    unsigned NumDummyLines=0;
    unsigned NumBufferLines=0;
    unsigned NumCol=0;
    unsigned NumRow=0;
    unsigned NumFuckedUpEntriesAtTheEnd=0;
    unsigned NumFuckedUpMomentumEntries=0;

    unsigned NumTomMomBins=0;
    double* TomMomentum=NULL;
    double* TomMomentumBins=NULL;
    double* pLab=NULL;

    if(TYPE==0){
        RadiusStepSize = 0.05;
        RadiusMinimum = 0.05;
        RadiusMaximum = 14.95;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        NumDummyLines=32;
        NumBufferLines=2;
        NumCol=10;
        NumRow=30;
        NumFuckedUpEntriesAtTheEnd=1;
        NumFuckedUpMomentumEntries=0;
        NumTomMomBins = 36;
        pLab = new double [NumTomMomBins];
        pLab[0] = 5;
        pLab[1] = 10;
        for(unsigned uBin=2; uBin<NumTomMomBins; uBin++){
            pLab[uBin] = pLab[uBin-1]+10;
        }
    }
    else if(TYPE==1){
        RadiusStepSize = 0.05;
        RadiusMinimum = 0.05;
        RadiusMaximum = 14.95;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        NumDummyLines=32;
        NumBufferLines=2;
        NumCol=10;
        NumRow=30;
        NumFuckedUpEntriesAtTheEnd=1;
        NumFuckedUpMomentumEntries=1;
        NumTomMomBins = 200;
        pLab = new double [NumTomMomBins];
        pLab[0] = 5.0;
        pLab[1] = 10.0;
        for(unsigned uBin=2; uBin<NumTomMomBins; uBin++){
            pLab[uBin] = pLab[uBin-1]+5;
        }
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pXi_ESC16_IS\n");
        return NULL;
    }

    const double Mass_p = 938.272;
    const double Mass_Xim = 1321.7;

    TomMomentum = new double [NumTomMomBins];
    TomMomentumBins = new double [NumTomMomBins+1];
    for(unsigned uBin=0; uBin<NumTomMomBins; uBin++){
        TomMomentum[uBin] = sqrt(pow(Mass_p*pLab[uBin],2)/(pow(Mass_p,2)+pow(Mass_Xim,2)+2*Mass_p*sqrt(pow(Mass_Xim,2)+pow(pLab[uBin],2))));
    }
    TomMomentumBins[0]=0;
    for(unsigned short uBin=0; uBin<NumTomMomBins-1; uBin++){
        TomMomentumBins[uBin+1] = (TomMomentum[uBin]+TomMomentum[uBin+1])*0.5;
    }
    TomMomentumBins[NumTomMomBins] = TomMomentum[NumTomMomBins-1]+(TomMomentum[NumTomMomBins-1]-TomMomentumBins[NumTomMomBins-1]);
    printf("N.B. The maximum momentum is %.2f MeV\n",TomMomentumBins[NumTomMomBins]);

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, 3312);

    Kitty.SetRedMass((Mass_p*Mass_Xim)/(Mass_p+Mass_Xim));

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum TomFiles {fS0, fS1};
    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[512];
        strcpy(InputFileName[uFile],InputFolder);
    }

    if(TYPE==0){
        strcat(InputFileName[fS0], "xnwave.ich=2.ic=1.xim-prot.lsj=000");
        strcat(InputFileName[fS1], "xnwave.ich=2.ic=1.xim-prot.lsj=011");
    }
    else if(TYPE==1){
        strcat(InputFileName[fS0], "xnwave.ich=2.ic=1.xim-prot.lsj=000.200");
        strcat(InputFileName[fS1], "xnwave.ich=2.ic=1.xim-prot.lsj=011.200");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pXi_ESC16_IS\n");
    }

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
        HistoWF[uCh][0].SetUp(2);
        HistoWF[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoWF[uCh][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoWF[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
        HistoPS[uCh][0].SetUp(1);
        HistoPS[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoPS[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoPS[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    char* cdummy = new char [512];

    unsigned WhichSpin=0;
    //unsigned WhichIsospin=0;

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        WhichSpin = uFile;
        //int WhichMomBin=-1;
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //Read the header lines
        for(unsigned short us=0; us<NumDummyLines; us++){
            if(!fgets(cdummy, 511, InFile)) continue;
        }

        float Values[NumCol];

        //int RadBin=-1;
        unsigned WhichBin[2];
        float fRadius;
        WhichBin[0]=0;
        WhichBin[1]=0;
        unsigned TotNumMomBins=0;
        //!---Iteration over all events---
        while(!feof(InFile)){
            //if(WhichBin[0]>=NumTomMomBins) {WhichBin[0]=0; WhichSpin=1;}
            for(unsigned uBuffer=0; uBuffer<NumBufferLines; uBuffer++){
                if(!fgets(cdummy, 511, InFile)) continue;
            }

            for(unsigned uRow=0; uRow<NumRow; uRow++){
                if(!fgets(cdummy, 511, InFile)) continue;
                if(uRow==NumRow-1&&NumFuckedUpEntriesAtTheEnd==1){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8]);
                }
                else if(NumFuckedUpEntriesAtTheEnd){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8],&Values[9]);
                }
                else{
                    printf("Issues in Init_pXi_ESC16_IS\n");
                }

                if(TotNumMomBins<NumFuckedUpMomentumEntries) continue;
                else if(uFile==1&&(TotNumMomBins%3!=0)) continue;

                for(unsigned uCol=0; uCol<NumCol; uCol++){
                    WhichBin[1] = uRow*NumCol+uCol;
                    fRadius = HistoWF[WhichSpin]->GetBinCenter(1,WhichBin[1]);
//printf("fRadius=%f\n",fRadius);
                    HistoWF[WhichSpin]->SetBinContent(WhichBin,Values[uCol]*fRadius);
                }
            }

            if(TotNumMomBins<NumFuckedUpMomentumEntries){TotNumMomBins++;continue;}
            else if(uFile==1&&(TotNumMomBins%3!=0)) {TotNumMomBins++;continue;}
            else TotNumMomBins++;

            HistoPS[WhichSpin]->SetBinContent(WhichBin,0);
            WhichBin[0]++;
        }

        fclose(InFile);
        if(WhichBin[0]-1!=NumTomMomBins){
            printf("\033[1;31mERROR:\033[0m WhichBin[0]!=NumTomMomBins (%u vs %u)\n",WhichBin[0],NumTomMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TomMomentum;
    delete [] TomMomentumBins;
    delete [] pLab;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS* Kitty, const int& TYPE){
    return Init_pXi_ESC16_Iavg_Coulomb(InputFolder,*Kitty,TYPE);
}



DLM_Histo<complex<double>>*** Init_pS0_ESC08(const char* InputFolder, CATS& Kitty, const int& TYPE){

    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    unsigned short NumFiles;
    //const unsigned NumSpinChannels = 2;
    //unsigned NumMomBins = 0;

    unsigned NumDummyLines=0;
    unsigned NumBufferLines=0;
    unsigned NumCol=0;
    unsigned NumRow=0;
    unsigned NumFuckedUpEntriesAtTheEnd=0;
    unsigned NumFuckedUpMomentumEntries=0;

    unsigned NumTomMomBins=0;
    double* TomMomentum=NULL;
    double* TomMomentumBins=NULL;
    double* pLab=NULL;

    if(TYPE==0){
        RadiusStepSize = 0.05;
        RadiusMinimum = 0.05;
        RadiusMaximum = 14.95;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        NumDummyLines=32;
        NumBufferLines=2;
        NumCol=10;
        NumRow=30;
        NumFuckedUpEntriesAtTheEnd=1;
        NumFuckedUpMomentumEntries=0;
        NumTomMomBins = 101;
        pLab = new double [NumTomMomBins];
        pLab[0] = 5;
        pLab[1] = 10;
        for(unsigned uBin=2; uBin<NumTomMomBins; uBin++){
            pLab[uBin] = pLab[uBin-1]+10;
        }
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pS0_ESC08\n");
        return NULL;
    }

    const double Mass_p = 938.272;
    const double Mass_S0 = 1192.642;

    TomMomentum = new double [NumTomMomBins];
    TomMomentumBins = new double [NumTomMomBins+1];
    for(unsigned uBin=0; uBin<NumTomMomBins; uBin++){
        TomMomentum[uBin] = sqrt(pow(Mass_p*pLab[uBin],2)/(pow(Mass_p,2)+pow(Mass_S0,2)+2*Mass_p*sqrt(pow(Mass_S0,2)+pow(pLab[uBin],2))));
    }
    TomMomentumBins[0]=0;
    for(unsigned short uBin=0; uBin<NumTomMomBins-1; uBin++){
        TomMomentumBins[uBin+1] = (TomMomentum[uBin]+TomMomentum[uBin+1])*0.5;
    }
    TomMomentumBins[NumTomMomBins] = TomMomentum[NumTomMomBins-1]+(TomMomentum[NumTomMomBins-1]-TomMomentumBins[NumTomMomBins-1]);
    printf("N.B. The maximum momentum is %.2f MeV\n",TomMomentumBins[NumTomMomBins]);

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3212);

    Kitty.SetRedMass((Mass_p*Mass_S0)/(Mass_p+Mass_S0));

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    enum TomFiles {fS0, fS1};
    char** InputFileName = new char* [NumFiles];
    char** InputFileNamePhases = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[512];
        strcpy(InputFileName[uFile],InputFolder);
        InputFileNamePhases[uFile] = new char[512];
        strcpy(InputFileNamePhases[uFile],InputFolder);
    }

    if(TYPE==0){
        strcat(InputFileName[fS0], "ynwave.si0-prot.lsj=000");
        strcat(InputFileName[fS1], "ynwave.si0-prot.lsj=011");
        strcat(InputFileNamePhases[fS0], "phases.si0-prot.jsj=000");
        strcat(InputFileNamePhases[fS1], "phases.si0-prot.jsj=011");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pS0_ESC08\n");
    }

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
        HistoWF[uCh][0].SetUp(2);
        HistoWF[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoWF[uCh][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoWF[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
        HistoPS[uCh][0].SetUp(1);
        HistoPS[uCh][0].SetUp(0,NumTomMomBins,TomMomentumBins);
        HistoPS[uCh][0].Initialize();
        for(unsigned uMomBin=0; uMomBin<NumTomMomBins; uMomBin++){
            HistoPS[uCh][0].SetBinCenter(0,uMomBin,TomMomentum[uMomBin]);
        }
    }

    char* cdummy = new char [512];
    char* buffer = new char [512];

    unsigned WhichSpin=0;
    //unsigned WhichIsospin=0;

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        WhichSpin = uFile;

        //int RadBin=-1;
        unsigned WhichBin[2];
        WhichBin[0]=0;
        WhichBin[1]=0;
        unsigned TotNumMomBins=0;

        //int WhichMomBin=-1;
        FILE *InFile;


        InFile = fopen(InputFileNamePhases[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileNamePhases[uFile]);
            return Histo;
        }
        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileNamePhases[uFile]);
            return Histo;
        }
        float PhaseShiftRe;
        float PhaseShiftIm;
        float pLab;
        complex<float> PhaseShift;

        //!---Iteration over all events---
        WhichBin[0]=0;
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(uFile==0) sscanf(cdummy,"%s %f %s %s %f %s %f %s",buffer,&pLab,buffer,buffer,&PhaseShiftRe,buffer,&PhaseShiftIm,buffer);
            else{
                //imagine: the last line has different format due to a space missing
                if(WhichBin[0]==NumTomMomBins-1){
                    sscanf(cdummy,"%s %s %s %s %s %f %s %f %s %f",buffer,buffer,buffer,buffer,buffer,&PhaseShiftRe,buffer,&PhaseShiftIm,buffer,&PhaseShiftIm);
                }
                else{
                    sscanf(cdummy,"%s %s %s %f %s %s %f %s %f %s %f",buffer,buffer,buffer,&pLab,buffer,buffer,&PhaseShiftRe,buffer,&PhaseShiftIm,buffer,&PhaseShiftIm);
                }

                PhaseShiftIm=0;
            }
//if(uFile==1) {printf("%s; %f\n",cdummy,PhaseShiftRe);}
            PhaseShiftRe *= (Pi/180.);
            PhaseShiftIm *= (Pi/180.);
            PhaseShift.real(PhaseShiftRe);
            PhaseShift.imag(PhaseShiftIm);
            HistoPS[WhichSpin]->SetBinContent(WhichBin,PhaseShift);
            WhichBin[0]++;
        }

        fclose(InFile);


        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //Read the header lines
        for(unsigned short us=0; us<NumDummyLines; us++){
            if(!fgets(cdummy, 511, InFile)) continue;
        }

        float Values[NumCol];
        complex<float> PsVal;
        complex<float> Smatrix;
        complex<float> WfVal;
        float RealPrefactor = 2;
        float fRadius;
        WhichBin[0]=0;
        //!---Iteration over all events---
        while(!feof(InFile)){
            //if(WhichBin[0]>=NumTomMomBins) {WhichBin[0]=0; WhichSpin=1;}
            for(unsigned uBuffer=0; uBuffer<NumBufferLines; uBuffer++){
                if(!fgets(cdummy, 511, InFile)) continue;
            }

            for(unsigned uRow=0; uRow<NumRow; uRow++){
                if(!fgets(cdummy, 511, InFile)) continue;
                if(uRow==NumRow-1&&NumFuckedUpEntriesAtTheEnd==1){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8]);
                }
                else if(NumFuckedUpEntriesAtTheEnd){
                    sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",
                       &Values[0],&Values[1],&Values[2],&Values[3],&Values[4],
                       &Values[5],&Values[6],&Values[7],&Values[8],&Values[9]);
                }
                else{
                    printf("Issues in Init_pXi_ESC16_IS\n");
                }

                if(TotNumMomBins<NumFuckedUpMomentumEntries) continue;
                else if(uFile==1&&(TotNumMomBins%3!=0)) continue;

                for(unsigned uCol=0; uCol<NumCol; uCol++){
                    WhichBin[1] = uRow*NumCol+uCol;
                    PsVal = HistoPS[WhichSpin]->GetBinContent(WhichBin);
                    Smatrix = exp(RealPrefactor*fi*PsVal);
                    Smatrix += 1.;
                    fRadius = HistoWF[WhichSpin]->GetBinCenter(1,WhichBin[1]);
//printf("fRadius=%f\n",fRadius);
                    WfVal = Smatrix*Values[uCol]*fRadius;
                    //WfVal = Values[uCol];
//if(uFile==0)
//printf("k = %.3f; Wf = %.3f * (%.3f + i*%.3f)\n",HistoWF[WhichSpin]->GetBinCenter(0,WhichBin[0]),Values[uCol],Smatrix.real(),Smatrix.imag());
                    HistoWF[WhichSpin]->SetBinContent(WhichBin,WfVal);
                }
            }

            if(TotNumMomBins<NumFuckedUpMomentumEntries){TotNumMomBins++;continue;}
            else if(uFile==1&&(TotNumMomBins%3!=0)) {TotNumMomBins++;continue;}
            else TotNumMomBins++;

            //HistoPS[WhichSpin]->SetBinContent(WhichBin,0);
            WhichBin[0]++;
        }

        fclose(InFile);



        if(WhichBin[0]-1!=NumTomMomBins){
            printf("\033[1;31mERROR:\033[0m WhichBin[0]!=NumTomMomBins (%u vs %u)\n",WhichBin[0],NumTomMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] cdummy;
    delete [] buffer;
    delete [] RadBins;
    delete [] TomMomentum;
    delete [] TomMomentumBins;
    delete [] pLab;

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        delete [] InputFileName[uFile];
        delete [] InputFileNamePhases[uFile];
    }
    delete [] InputFileName;
    delete [] InputFileNamePhases;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pS0_ESC08(const char* InputFolder, CATS* Kitty, const int& TYPE){
    return Init_pS0_ESC08(InputFolder,*Kitty,TYPE);
}



DLM_Histo<complex<double>>*** Init_pS0_ESC16(const char* InputFolder, CATS& Kitty, const int& TYPE){

    const unsigned NumFiles=2;
    const unsigned NumChannels=2;
    char** FileNames = new char* [NumFiles];
    if(TYPE==0){
        for(unsigned uFile=0; uFile<NumFiles; uFile++){
            FileNames[uFile] = new char [512];
        }
        strcpy(FileNames[0],InputFolder);
        strcat(FileNames[0],"lsj000.dat");

        strcpy(FileNames[1],InputFolder);
        strcat(FileNames[1],"lsj011.dat");
    }
    else{
        printf("\033[1;31mERROR:\033[0m Init_pS0_ESC16 says you used a wrong TYPE\n");
        return NULL;
    }

    const double Mass_p = 938.272;
    const double Mass_S0 = 1192.642;

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,1);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3212);

    Kitty.SetRedMass((Mass_p*Mass_S0)/(Mass_p+Mass_S0));

    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];
    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];

    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];

    char* cdummy = new char [512];

    float fMomentum;
    float fMomentumOld=-1000;
    float fRadius;
    float fRadiusOld=-1000;
    float fNormRe;
    float fNormIm;
    float fWfRe;
    float fWfIm;

    FILE *InFile;
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        unsigned MaxNumMomBins = 1000;
        unsigned NumMomBins = 0;
        double* MomBins = new double [MaxNumMomBins+1];
        double* MomBinsCenter = new double [MaxNumMomBins];

        unsigned MaxNumRadBins = 1000;
        unsigned NumRadBins = 0;
        double* RadBins = new double [MaxNumRadBins+1];
        double* RadBinsCenter = new double [MaxNumRadBins];

        InFile = fopen(FileNames[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", FileNames[uFile]);
            abort();
        }
        MomBins[0] = 0;
        //we iterate only to get the binning
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            sscanf(cdummy, "%f %f %f %f %f %f",&fMomentum,&fRadius,&fNormRe,&fNormIm,&fWfRe,&fWfIm);
            if(fMomentum!=fMomentumOld){
                if(NumMomBins==0){
                    fMomentumOld = fMomentum;
                    MomBinsCenter[NumMomBins] = fMomentum;
                    NumMomBins++;
                }
                else{
                    MomBins[NumMomBins] = (fMomentumOld+fMomentum)*0.5;
                    MomBinsCenter[NumMomBins] = fMomentum;
                    fMomentumOld = fMomentum;
                    NumMomBins++;
                }
            }

            if(NumRadBins==0){
                fRadiusOld = fRadius;
                RadBinsCenter[NumRadBins] = fRadius;
                NumRadBins++;
            }
            else if(NumMomBins<=1){
                RadBins[NumRadBins] = (fRadiusOld+fRadius)*0.5;
                RadBinsCenter[NumRadBins] = fRadius;
                fRadiusOld = fRadius;
                NumRadBins++;
            }

        }
        fclose(InFile);

        MomBins[NumMomBins] = MomBins[NumMomBins-1] + 2.*(MomBinsCenter[NumMomBins-1]-MomBins[NumMomBins-1]);
        RadBins[NumRadBins] = RadBins[NumRadBins-1] + 2.*(RadBinsCenter[NumRadBins-1]-RadBins[NumRadBins-1]);

        HistoWF[uFile] = new DLM_Histo<complex<double>> [1];
        HistoWF[uFile][0].SetUp(2);
        HistoWF[uFile][0].SetUp(0,NumMomBins,MomBins);
        HistoWF[uFile][0].SetUp(1,NumRadBins,RadBins);
        HistoWF[uFile][0].Initialize();

        HistoPS[uFile] = new DLM_Histo<complex<double>> [1];
        HistoPS[uFile][0].SetUp(1);
        HistoPS[uFile][0].SetUp(0,NumMomBins,MomBins);
        HistoPS[uFile][0].Initialize();
/*
        if(uFile==0){
            printf("NumMomBins=%u\n",NumMomBins);
            for(unsigned uBin=0; uBin<NumMomBins; uBin++){
                printf(" %.3f -> %.3f -> %.3f\n",MomBins[uBin],MomBinsCenter[uBin],MomBins[uBin+1]);
            }
            printf("NumRadBins=%u\n",NumRadBins);
            for(unsigned uRad=0; uRad<NumRadBins; uRad++){
                printf(" %.3f -> %.3f -> %.3f\n",RadBins[uRad],RadBinsCenter[uRad],RadBins[uRad+1]);
            }
        }
*/


        unsigned WhichBin[2];
        complex<float> WfVal;
        InFile = fopen(FileNames[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", FileNames[uFile]);
            abort();
        }
        //we iterate a second time, in order to actually read the wave function
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            sscanf(cdummy, "%f %f %f %f %f %f",&fMomentum,&fRadius,&fNormRe,&fNormIm,&fWfRe,&fWfIm);
            WfVal.real(fWfRe);
            WfVal.imag(fWfIm);
            WfVal *= fRadius;
            WhichBin[0] = HistoWF[uFile][0].GetBin(0,fMomentum);
            WhichBin[1] = HistoWF[uFile][0].GetBin(1,fRadius);
            HistoWF[uFile][0].SetBinContent(WhichBin,WfVal);
        }
        fclose(InFile);

        if(MomBins) delete [] MomBins;
        if(MomBinsCenter) delete [] MomBinsCenter;
        if(RadBins) delete [] RadBins;
        if(RadBinsCenter) delete [] RadBinsCenter;
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        if(FileNames&&FileNames[uFile]) delete [] FileNames[uFile];
    }
    if(FileNames) delete [] FileNames;
    delete [] cdummy;

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pS0_ESC16(const char* InputFolder, CATS* Kitty, const int& TYPE){
    return Init_pS0_ESC16(InputFolder,*Kitty,TYPE);
}

//TYPE 0 = ppbar->ppbar
//TYPE 1 = + nnbar->ppbar
DLM_Histo<complex<double>>*** Init_pantip_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE){
    //SET UP Q1Q2... Kitty
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    unsigned short NumChannels;
    //unsigned short NumPwPerCh;
    //unsigned short NumFiles;
    unsigned NumMomBins = 0;
    unsigned short NumFiles = 7; // (to cover all s,p,d, waves)
    bool* TakeThisFile = new bool [NumFiles];

    //The TYPE selects the different configurations of describing the p-antip CF
    //TYPE=0 -> only the p-antip->p-antip wf
    //TYPE=1 -> p-natip->p-antip + n-antin->p-antip
    //TYPE=2 -> includes the coupled-channel to 2 mesons in an approximate way in 1S0 state:
    //          p-antip->p-antip + n-antin->p-antip + 2M->p-antip
    //TYPE=3 -> includes the coupled-channel to 2 mesons in an approximate way in 3S1 state:
    //          p-antip->p-antip + n-antin->p-antip + 2M->p-antip
    //TYPE=4 -> includes the coupled-channel to 2 mesons in an approximate way in 1S0 + 3S1 states
    //TYPE=5 -> includes the coupled-channel to 2 mesons in an approximate way in 1S0 + 3S1 + 3P0 states
    //TYPE=6 -> includes the coupled-channel to 2 mesons in an approximate way in 1S0 + 3S1 + 3P0 +1P1 states
    //TYPE=7 -> includes the coupled-channel to 2 mesons in EXACT way with wfs in 1S0 states
//Update 07.10.2020 based on Haidenbauer suggestions on Annih. cross-sections in PWs (see Annxs.data file)
    //TYPE=8 -> includes the coupled-channel to 2 mesons in an approximate way in 1S0 + 3S1 + 3P2 +3P0 states
    //TYPE=9 -> includes the coupled-channel to 2 mesons in an approximate way in ALL PARTIAL WAVES


   // VERY IMPORTANT: whenever you add a new configuration, increase of +1 this number
    int NumTYPE = 9;
    if(TYPE==0){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 4;
        NumFiles = 7;
    }
    else if(TYPE==1){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 8;
        NumFiles = 7;
    } else if(TYPE==2){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 9;//additional "2pi" channel only for 1S0
        NumFiles = 7;
    } else if(TYPE==3){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 9;//additional "2pi" channel only for 3S1
        NumFiles = 7;
    }  else if(TYPE==4){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 10;//additional "2pi" channel for 1S0 and 3S1
        NumFiles = 7;
    }  else if(TYPE==5){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 11;//additional "2pi" channel for 1S0, 3S1,3P0
        NumFiles = 7;
    }  else if(TYPE==6){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 12;//additional "2pi" channel for 1S0,1P1,3S1,3P0
        NumFiles = 7;
    } else if(TYPE==7){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 9;//real 2pi channel for 1S0
        NumFiles = 7;//additional file with WFS for 2mesons
    } else if(TYPE==8){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 12;//additional "2pi" channel for 1S0, 3S1,3P2, 3P0
        NumFiles = 7;
    } else if(TYPE==9){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 14;//additional "2pi" channel for 1S0, 3S1, 1P1, 3P0, 3P1, 3P2
        // NumChannels = 9;//additional "2pi" channel for 1S0, 3S1, 1P1, 3P0, 3P1, 3P2
        NumFiles = 7;
    }
    else{
        printf("HELLO MY DEAR VALENTINA: YOU BROKE SOMETHING in Init_pantip_Haidenbauer at the very beginning\n");
        return NULL;
    }
// channels for CATS object
// ch 0: 1S0 -> 1/4
// ch 1: 3S1, 3P0 -> 3/4 x 1/9 = 1/12
// ch 2: 3S1 -> 3/4 x 8/9 = 2/3
// Update 13.07.2020
// ch 0: 1S0 + 1P1 -> 1/4
// ch 1: 3S1 + 3P0 -> 3/4 x 1/9 = 1/12
// ch 2: 3S1 + 3P1 -> 3/4 x 1/3 = 1/4
// ch 3: 3S1 + 3P2 -> 3/4 x 5/9 = 5/12
// The integer TYPE selects different
//Different selections of potentials setup:
//HAIDE_1 -> only p-antip->p-antip (TYPE 0)
//HAIDE_2 -> p-antip->p-antip + n-antin->p-antip (TYPE 1)
//HAIDE_3 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 1S0 (TYPE 2)
//HAIDE_4 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 3S1 (TYPE 3)
//HAIDE_5 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 1S0 + 3S1 (TYPE 4)
//HAIDE_6 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 1S0 + 3S1 + 3P0 (TYPE 5)
//HAIDE_7 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 1S0 + 3S1 + 1P1 + 3P0 (TYPE 6)
//Update 04.10.2020
//HAIDE_8 -> p-antip->p-antip + n-antin->p-antip + REAL 2M -> p-antip in 1S0 (TYPE 7)
//Update 07.10.2020 inclusion of 3P1 and 3P2 for the Migdal-Watson approximation
//HAIDE_9 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in 1S0 + 3S1 + 3P2 + 3P0 (TYPE 8)
//Update 16.10.2020 inclusion of ALL PWs for the Migdal-Watson approximation
//HAIDE_10 -> p-antip->p-antip + n-antin->p-antip + 2M -> p-antip in ALL PWS (TYPE 9)


    Kitty.SetNumChannels(NumChannels);
    if(TYPE==0){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);
    }
    else if(TYPE==1){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);
    }
    else if (TYPE==2){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);
        //"2pi"->ppbar
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);
    }
    else if (TYPE==3){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);
        //"2pi"->ppbar
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,3./4.);
        Kitty.SetSpin(8,1);
    }
    else if (TYPE==4){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);

        //"2pi"->ppbar in 1S0 and 3S1
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);

        Kitty.SetNumPW(9,1);
        Kitty.SetChannelWeight(9,3./4.);
        Kitty.SetSpin(9,1);
    }
        else if (TYPE==5){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);

        //"2pi"->ppbar in 1S0
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);
        //"2pi"->ppbar in 3S1
        Kitty.SetNumPW(9,1);
        Kitty.SetChannelWeight(9,3./4.);
        Kitty.SetSpin(9,1);
        //"2pi"->ppbar in 3P0
        Kitty.SetNumPW(10,1);
        Kitty.SetChannelWeight(10,1./12.);
        Kitty.SetSpin(10,1);
    }
        else if (TYPE==6){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);

        //"2pi"->ppbar in 1S0
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);
        //"2pi"->ppbar in 3S1
        Kitty.SetNumPW(9,1);
        Kitty.SetChannelWeight(9,3./4.);
        Kitty.SetSpin(9,1);
        //"2pi"->ppbar in 1P1
        Kitty.SetNumPW(10,1);
        Kitty.SetChannelWeight(10,1./4.);
        Kitty.SetSpin(10,0);
        //"2pi"->ppbar in 3P0
        Kitty.SetNumPW(11,1);
        Kitty.SetChannelWeight(11,1./12.);//3/4 x 1/9
        Kitty.SetSpin(11,1);
    }    else if (TYPE==7){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);
        //Real 2pi->ppbar
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);
    }
    else if (TYPE==8){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);

        //"2pi"->ppbar in 1S0
        Kitty.SetNumPW(8,1);
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);
        //"2pi"->ppbar in 3S1
        Kitty.SetNumPW(9,1);
        Kitty.SetChannelWeight(9,3./4.);
        Kitty.SetSpin(9,1);
        //"2pi"->ppbar in 3P0
        Kitty.SetNumPW(10,1);
        Kitty.SetChannelWeight(10,1./12.);//3/4 x 1/9
        Kitty.SetSpin(10,1);
        //"2pi"->ppbar in 3P2
        Kitty.SetNumPW(11,1);
        Kitty.SetChannelWeight(11,5./12.);
        Kitty.SetSpin(11,1);
    }  else if(TYPE==9){
        //ppbar->ppbar
        Kitty.SetNumPW(0,2);
        Kitty.SetChannelWeight(0,1./4.);
        Kitty.SetSpin(0,0);

        Kitty.SetNumPW(1,2);
        Kitty.SetChannelWeight(1,1./12.);
        Kitty.SetSpin(1,1);

        Kitty.SetNumPW(2,2);
        Kitty.SetChannelWeight(2,1./4.);
        Kitty.SetSpin(2,1);

        Kitty.SetNumPW(3,2);
        Kitty.SetChannelWeight(3,5./12.);
        Kitty.SetSpin(3,1);

        //nnbar->ppbar
        Kitty.SetNumPW(4,2);
        Kitty.SetChannelWeight(4,1./4.);
        Kitty.SetSpin(4,0);

        Kitty.SetNumPW(5,2);
        Kitty.SetChannelWeight(5,1./12.);
        Kitty.SetSpin(5,1);

        Kitty.SetNumPW(6,2);
        Kitty.SetChannelWeight(6,1./4.);
        Kitty.SetSpin(6,1);

        Kitty.SetNumPW(7,2);
        Kitty.SetChannelWeight(7,5./12.);
        Kitty.SetSpin(7,1);

        //X->ppbar
        Kitty.SetNumPW(8,1);//1S0
        Kitty.SetChannelWeight(8,1./4.);
        Kitty.SetSpin(8,0);

        Kitty.SetNumPW(9,1);//3S1
        Kitty.SetChannelWeight(9,3./4.);
        Kitty.SetSpin(9,1);

        Kitty.SetNumPW(10,1);//1P1
        Kitty.SetChannelWeight(10,1./4.);
        Kitty.SetSpin(10,1);

        Kitty.SetNumPW(11,1);//3P0
        Kitty.SetChannelWeight(11,1./12.);
        Kitty.SetSpin(11,1);

        Kitty.SetNumPW(12,1);//3P1
        Kitty.SetChannelWeight(12,1./4.);
        Kitty.SetSpin(12,1);

        Kitty.SetNumPW(13,1);//3P2
        Kitty.SetChannelWeight(13,5./12.);
        Kitty.SetSpin(13,1);
    }
    else{
        printf("Да му еба майката!\n");
    }


    Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, -2212);
    const double Mass_p = 938.272;
    Kitty.SetRedMass((Mass_p*Mass_p)/(Mass_p+Mass_p));
    //only take the s-waves into account, i.e. C(k)->0 as for coupled channels
    // CATS does not do anything for the asymptotic solution
    if(TYPE%10==1){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
    } else if (TYPE%10==2 || TYPE%10==7){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
    } else if (TYPE%10==3){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
    } else if (TYPE%10==4){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
        Kitty.SetOnlyNumericalPw(9,true);
    } else if (TYPE%10==5){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
        Kitty.SetOnlyNumericalPw(9,true);
        Kitty.SetOnlyNumericalPw(10,true);
    } else if (TYPE%10==6 || TYPE%10==8){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
        Kitty.SetOnlyNumericalPw(9,true);
        Kitty.SetOnlyNumericalPw(10,true);
        Kitty.SetOnlyNumericalPw(11,true);
    } else if (TYPE%10==9){
        Kitty.SetOnlyNumericalPw(4,true);
        Kitty.SetOnlyNumericalPw(5,true);
        Kitty.SetOnlyNumericalPw(6,true);
        Kitty.SetOnlyNumericalPw(7,true);
        Kitty.SetOnlyNumericalPw(8,true);
        Kitty.SetOnlyNumericalPw(9,true);
        Kitty.SetOnlyNumericalPw(10,true);
        Kitty.SetOnlyNumericalPw(11,true);
        Kitty.SetOnlyNumericalPw(12,true);
        Kitty.SetOnlyNumericalPw(13,true);
    }

    //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
    //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=1; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }


    //char* buffer = new char [512];
    if(TYPE<=NumTYPE){
        //1-10 is singlet
        //11-20 is triplet
        //for(unsigned uFile=0; uFile<NumFiles; uFile++){
        //    sprintf(buffer, "w%s%u.%s",uFile==9?"":"0",uFile<10?uFile+1:uFile+1-10,uFile<10?"31":"32");
        //    strcat(InputFileName[uFile], buffer);
        //    printf("%s\n",InputFileName[uFile]);
        //}
        strcat(InputFileName[0], "w.31");
        strcat(InputFileName[1], "w.32");
        strcat(InputFileName[2], "w_29102020.30");//3P0 new 29.10.2020
        strcat(InputFileName[3], "w.11");//1P1
        strcat(InputFileName[4], "w.12");//3P1
        strcat(InputFileName[5], "w_2007.13");//3P2 new 20.07.2020
        strcat(InputFileName[6], "wX1s0.data");//real 2M 1S0 new 04.10.2020
      //  for(unsigned uFile=0; uFile<NumFiles; uFile++){
      //    printf("%s\n",InputFileName[uFile]);
      //  }
    }
    else{
        printf("HELLO MY DEAR VALENTINA: YOU BROKE SOMETHING in Init_pantip_Haidenbauer\n");
    }

    //unsigned WhichMomBin;

    const unsigned MaxNumMomBins = 256;
    double* MomentumBins = new double [MaxNumMomBins+1];
    double* Momentum = new double [MaxNumMomBins];

    char* cdummy = new char [512];
    float fMomentum;

    FILE *InFile;
    InFile = fopen(InputFileName[0], "r");//since it has #mombins=180 and the others have 200!!
    // InFile = fopen(InputFileName[6], "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[0]);
        return NULL;
    }
    int CurrentLine=0;

    MomentumBins[0] = 0;
    while(!feof(InFile)){
        if(!fgets(cdummy, 511, InFile)) continue;
        if((CurrentLine)%int(NumRadBins)==0){
            sscanf(cdummy, "%f",&fMomentum);
            Momentum[NumMomBins] = fMomentum;
            // printf("cdummy = %s\n",cdummy);
            // printf("Momentum[%u]=%.2f\n",NumMomBins,fMomentum);
            if(NumMomBins){
                //set the bin range in between the last two bin centers
                MomentumBins[NumMomBins] = 0.5*(Momentum[NumMomBins]+Momentum[NumMomBins-1]);
//printf("MomentumBins[%u]=%f\n",NumMomBins,MomentumBins[NumMomBins]);
            }
            NumMomBins++;
        }
        CurrentLine++;
      // printf("TOTAL NUMMOMBINS = %i\n",NumMomBins);
      // printf("MAX MOM = %.2f\n",MomentumBins[NumMomBins]);
    }
    fclose(InFile);
      // printf("TOTAL NUMMOMBINS = %i\n",NumMomBins);
      // printf("MAX MOM = %.2f\n",MomentumBins[NumMomBins]);

    // return 0;
    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    MomentumBins[NumMomBins] = 2.*Momentum[NumMomBins-1]-MomentumBins[NumMomBins-1];


    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uHist)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uHist); uPw++){
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins,MomentumBins);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uHist)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uHist); uPw++){
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins,MomentumBins);
            HistoPS[uHist][uPw].Initialize();
        }
    }

//printf("NumDLM_Histos=%u\n",NumDLM_Histos);
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
//if(uFile>2) break;
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        // if(uFile<=5)InFile = fopen(InputFileName[uFile], "r");
        // else if(uFile==6) InFile = fopen(InputFileName[6], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        //unsigned LastRadBin;

        float fRadius;
        float fReWf_ppbar;
        float fImWf_ppbar;
        float fReWf_nnbar;
        float fImWf_nnbar;
        float fReWf_2M;
        float fImWf_2M;
        float fDummy;

        int RadBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            if(!fgets(cdummy, 511, InFile)) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Reading file:: %s\n",InputFileName[uFile]);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);

                RadBin++;
                WhichMomBin++;
                continue;
            }

            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            // printf("NumMomBins Haide = ");
            //we fill up the radius bins with and offset of 1, due to the special zeroth bin
            WhichBin[1] = RadBin+1;

            if(TYPE==0){
              if(uFile==0){//1S0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[0][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){//3S1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[1][0].SetBinContent(WhichBin,0);
                HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[2][0].SetBinContent(WhichBin,0);
                HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[3][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[1][1].SetBinContent(WhichBin,0);
              } else if(uFile==3){//1P1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[0][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[2][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
                HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
                HistoPS[3][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==1){
              if(uFile==0){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);
              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==2){
              if(uFile==0){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              //only 1S0 for channel 8 "2pi"->p-antip
              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);
              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==3){
              if(uFile==0){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              //only 3S1 for channel 8 "2pi"->p-antip
              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);
              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==4){
              if(uFile==0){//1S0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){//3S1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              //"2M"->p-antip
              HistoWF[9][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[9][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);

              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==5){
              if(uFile==0){//1S0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){//3S1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              //"2M"->p-antip
              HistoWF[9][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[9][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);

              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);

              HistoWF[10][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[10][0].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            } else if(TYPE==6){
              if(uFile==0){//1S0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){//3S1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              //"2M"->p-antip
              HistoWF[9][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[9][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);

              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);

              HistoWF[11][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[11][0].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);

              HistoWF[10][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[10][0].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              }
            }  else if(TYPE==7){
              if(uFile==0){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);
              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);
              } else if(uFile==6){//2M in 1S0
              sscanf(cdummy, " %f %f %f %f",&fRadius,&fDummy,&fReWf_2M,&fImWf_2M);
              //only 1S0 for channel 8 real 2M->p-antip
              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_2M+fi*fImWf_2M)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              }
            }  else if(TYPE==8){
              if(uFile==0){//1S0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){//3S1
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              //"2M"->p-antip
              HistoWF[9][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[9][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);

              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);

              //"2M"->p-antip
              HistoWF[10][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fReWf_ppbar)*fRadius);
              HistoPS[10][0].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);

              HistoWF[11][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[11][0].SetBinContent(WhichBin,0);
              }
            }  else if(TYPE==9){
              if(uFile==0){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][0].SetBinContent(WhichBin,0);

              HistoWF[4][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][0].SetBinContent(WhichBin,0);

              //only 1S0 for channel 8 "2pi"->p-antip
              HistoWF[8][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[8][0].SetBinContent(WhichBin,0);
              } else if(uFile==1){
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);
              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);

              HistoWF[5][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][0].SetBinContent(WhichBin,0);
              HistoWF[6][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][0].SetBinContent(WhichBin,0);
              HistoWF[7][0].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][0].SetBinContent(WhichBin,0);

              HistoWF[9][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[9][0].SetBinContent(WhichBin,0);
              } else if(uFile==2){//3P0
                sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[1][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[1][1].SetBinContent(WhichBin,0);
              HistoWF[5][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[5][1].SetBinContent(WhichBin,0);

              HistoWF[11][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[11][0].SetBinContent(WhichBin,0);
              }else if(uFile==3){//1P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[0][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[0][1].SetBinContent(WhichBin,0);

              HistoWF[4][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[4][1].SetBinContent(WhichBin,0);

              HistoWF[10][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[10][0].SetBinContent(WhichBin,0);
              } else if(uFile==4){//3P1
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[2][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[2][1].SetBinContent(WhichBin,0);

              HistoWF[6][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[6][1].SetBinContent(WhichBin,0);

              HistoWF[12][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[12][0].SetBinContent(WhichBin,0);
              } else if(uFile==5){//3P2
              sscanf(cdummy, " %f %f %f %f %f %f",&fRadius,&fDummy,&fReWf_ppbar,&fImWf_ppbar,&fReWf_nnbar,&fImWf_nnbar);
              HistoWF[3][1].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[3][1].SetBinContent(WhichBin,0);

              HistoWF[7][1].SetBinContent(WhichBin,(fReWf_nnbar+fi*fImWf_nnbar)*fRadius);
              HistoPS[7][1].SetBinContent(WhichBin,0);

              HistoWF[13][0].SetBinContent(WhichBin,(fReWf_ppbar+fi*fImWf_ppbar)*fRadius);
              HistoPS[13][0].SetBinContent(WhichBin,0);
              }
            }


            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }
        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TakeThisFile;

    return Histo;
}


//TYPE 0 is strong only
//TYPE 1 is with coulomb
DLM_Histo<complex<double>>*** Init_pKminus_Kyoto2019(const char* InputFolder, CATS& Kitty, const int& TYPE){
    //SET UP Q1Q2... Kitty
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    //unsigned short NumFiles;
    unsigned NumMomBins = 0;
    unsigned short NumFiles = 0;
    bool* TakeThisFile = new bool [NumFiles];

    if(TYPE==0){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.0;
        RadiusMaximum = 20.;
        NumChannels = 6;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 78;
    }
    else if(TYPE==1){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.0;
        RadiusMaximum = 64.;
        NumChannels = 6;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 78;
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pKminus_Kyoto2019\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,NumPwPerCh);
    Kitty.SetChannelWeight(0,1.);
    Kitty.SetSpin(0,0);

    for(unsigned uCh=1; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,1.);
        Kitty.SetSpin(uCh,0);
        Kitty.SetOnlyNumericalPw(uCh,true);
    }

    if(TYPE==0) Kitty.SetQ1Q2(0);
    else Kitty.SetQ1Q2(-1);
    Kitty.SetPdgId(2212, -321);
    const double Mass_p = 938.272;
    const double Mass_KaonCh = 493.677;
    Kitty.SetRedMass((Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh));

    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1;
//printf("NumRadBins=%u\n",NumRadBins);
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad)*RadiusStepSize-0.5*RadiusStepSize;
        //printf("RadBins[%u] = %.4f\n",uRad,RadBins[uRad]);
        RadBinLoaded[uRad] = false;
    }

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    //char* buffer = new char [512];
    if(TYPE==0){
        strcat(InputFileName[0], "wf_Kyoto_cc_woC.dat");
    }
    else if(TYPE==1){
        strcat(InputFileName[0], "wf64_Kyoto_cc_C.dat");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pKminus_Kyoto2019\n");
    }


    double* MomentumBins = new double [NumMomBins+1];
    double* Momentum = new double [NumMomBins];
    MomentumBins[0] = 0;
    Momentum[0] = 2;
    MomentumBins[1] = 3;
    for(unsigned uBin=1; uBin<60; uBin++){
        MomentumBins[uBin+1] = MomentumBins[uBin]+2;
        Momentum[uBin] = Momentum[uBin-1]+2;
    }
    Momentum[60] = 130;
    MomentumBins[61] = 135;
    for(unsigned uBin=61; uBin<NumMomBins; uBin++){
        MomentumBins[uBin+1] = MomentumBins[uBin]+10;
        Momentum[uBin] = Momentum[uBin-1]+10;
    }
    //for(unsigned uBin=0; uBin<NumMomBins; uBin++){
//printf("%.2f - %.2f - %.2f\n",MomentumBins[uBin],Momentum[uBin],MomentumBins[uBin+1]);
    //}

    char* cdummy = new char [512];
    //float fMomentum;




    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoPS[uHist][uPw].Initialize();
        }
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        int WhichMomBin=0;
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }


        //unsigned LastRadBin;

        float fMomentum;
        float fRadius;
        float fReWf_Kp;
        float fImWf_Kp;
        float fReWf_Kn;
        float fImWf_Kn;
        float fReWf_pipSm;
        float fImWf_pipSm;
        float fReWf_pi0S0;
        float fImWf_pi0S0;
        float fReWf_pimSp;
        float fImWf_pimSp;
        float fReWf_pi0L;
        float fImWf_pi0L;

        int RadBin=0;
        WhichMomBin=0;
        //!---Iteration over all events---
        while(!feof(InFile)){
            //ignore empty lines between the momentum bins
            if(!fgets(cdummy, 511, InFile)||strlen(cdummy)<10) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }
//printf("READING THE NEXT LINE ---> \n");
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                    &fMomentum, &fRadius,
                    &fReWf_Kp,&fImWf_Kp,
                    &fReWf_Kn,&fImWf_Kn,
                    &fReWf_pipSm,&fImWf_pipSm,
                    &fReWf_pi0S0,&fImWf_pi0S0,
                    &fReWf_pimSp,&fImWf_pimSp,
                    &fReWf_pi0L,&fImWf_pi0L);

           // printf("-------------(%i %i)\n%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
           //         RadBin,WhichMomBin,fMomentum, fRadius,
           //         fReWf_Kp,fImWf_Kp,
           //         fReWf_Kn,fImWf_Kn,
           //         fReWf_pipSm,fImWf_pipSm,
           //         fReWf_pi0S0,fImWf_pi0S0,
           //         fReWf_pimSp,fImWf_pimSp,
           //         fReWf_pi0L,fImWf_pi0L);
//usleep(62.5e3);
            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            WhichBin[1] = RadBin;

            HistoWF[0][0].SetBinContent(WhichBin,(fReWf_Kp+fi*fImWf_Kp)*fRadius);
            HistoPS[0][0].SetBinContent(WhichBin,0);

            HistoWF[1][0].SetBinContent(WhichBin,(fReWf_Kn+fi*fImWf_Kn)*fRadius);
            HistoPS[1][0].SetBinContent(WhichBin,0);

            HistoWF[2][0].SetBinContent(WhichBin,(fReWf_pipSm+fi*fImWf_pipSm)*fRadius);
            HistoPS[2][0].SetBinContent(WhichBin,0);

            HistoWF[3][0].SetBinContent(WhichBin,(fReWf_pi0S0+fi*fImWf_pi0S0)*fRadius);
            HistoPS[3][0].SetBinContent(WhichBin,0);

            HistoWF[4][0].SetBinContent(WhichBin,(fReWf_pimSp+fi*fImWf_pimSp)*fRadius);
            HistoPS[4][0].SetBinContent(WhichBin,0);

            HistoWF[5][0].SetBinContent(WhichBin,(fReWf_pi0L+fi*fImWf_pi0L)*fRadius);
            HistoPS[5][0].SetBinContent(WhichBin,0);

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)){
                RadBin=0;
                WhichMomBin++;
//printf("SWITCH\n");
//usleep(500e3);
            }
        }

        fclose(InFile);
        if(WhichMomBin!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TakeThisFile;

    return Histo;
}


//Type == 0 K0p model 1
//Type == 1 K0p model 2
//Type == 10 Kbar0p
DLM_Histo<complex<double>>*** Init_pK0_Kyoto2019(const char* InputFolder, CATS& Kitty, const int& TYPE){
    //SET UP Q1Q2... Kitty
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    //unsigned short NumFiles;
    unsigned NumMomBins = 0;
    unsigned short NumFiles = 0;
    bool* TakeThisFile = new bool [NumFiles];

    if(TYPE==0){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.0;
        RadiusMaximum = 20.;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 78;
    }
    else if(TYPE==1){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.0;
        RadiusMaximum = 20.;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 78;
    }
    else if(TYPE==10){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.0;
        RadiusMaximum = 20.;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 78;
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pK0_Kyoto2019\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,NumPwPerCh);
    Kitty.SetChannelWeight(0,1.);
    Kitty.SetSpin(0,0);

    for(unsigned uCh=1; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,1.);
        Kitty.SetSpin(uCh,0);
        Kitty.SetOnlyNumericalPw(uCh,true);
    }

    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212,TYPE<10?311:-311);
    const double Mass_p = 938.272;
    const double Mass_K0 = 497.648;
    Kitty.SetRedMass((Mass_p*Mass_K0)/(Mass_p+Mass_K0));

    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1;
//printf("NumRadBins=%u\n",NumRadBins);
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad)*RadiusStepSize-0.5*RadiusStepSize;
        //printf("RadBins[%u] = %.4f\n",uRad,RadBins[uRad]);
        RadBinLoaded[uRad] = false;
    }

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    //char* buffer = new char [512];
    if(TYPE==0){
        strcat(InputFileName[0], "wf_K0p_model1.dat");
    }
    else if(TYPE==1){
        strcat(InputFileName[0], "wf_K0p_model2.dat");
    }
    else if(TYPE==10){
        strcat(InputFileName[0], "wf_Kbar0p.dat");
    }
    else{
        printf("YOU BROKE SOMETHING in Init_pK0_Kyoto2019\n");
    }


    double* MomentumBins = new double [NumMomBins+1];
    double* Momentum = new double [NumMomBins];
    MomentumBins[0] = 0;
    Momentum[0] = 2;
    MomentumBins[1] = 3;
    for(unsigned uBin=1; uBin<60; uBin++){
        MomentumBins[uBin+1] = MomentumBins[uBin]+2;
        Momentum[uBin] = Momentum[uBin-1]+2;
    }
    Momentum[60] = 130;
    MomentumBins[61] = 135;
    for(unsigned uBin=61; uBin<NumMomBins; uBin++){
        MomentumBins[uBin+1] = MomentumBins[uBin]+10;
        Momentum[uBin] = Momentum[uBin-1]+10;
    }
    //for(unsigned uBin=0; uBin<NumMomBins; uBin++){
//printf("%.2f - %.2f - %.2f\n",MomentumBins[uBin],Momentum[uBin],MomentumBins[uBin+1]);
    //}

    char* cdummy = new char [512];
    //float fMomentum;

    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoPS[uHist][uPw].Initialize();
        }
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        int WhichMomBin=0;
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }

        float fMomentum;
        float fRadius;
        float fReWf_Kp;
        float fImWf_Kp;
        float fReWf_Kn;
        float fImWf_Kn;
        float fReWf_pipS0;
        float fImWf_pipS0;
        float fReWf_pi0Sp;
        float fImWf_pi0Sp;
        float fReWf_pipL;
        float fImWf_pipL;

        int RadBin=0;
        WhichMomBin=0;
        //!---Iteration over all events---
        while(!feof(InFile)){
            //ignore empty lines between the momentum bins
            if(!fgets(cdummy, 511, InFile)||strlen(cdummy)<10) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }
//printf("READING THE NEXT LINE ---> \n");
            if(TYPE<10) sscanf(cdummy, " %f %f %f %f %f %f",
                                &fMomentum, &fRadius,
                                &fReWf_Kp,&fImWf_Kp,
                                &fReWf_Kn,&fImWf_Kn);
            else        sscanf(cdummy, " %f %f %f %f %f %f %f %f %f %f",
                                &fMomentum, &fRadius,
                                &fReWf_Kp,&fImWf_Kp,
                                &fReWf_pipS0,&fImWf_pipS0,
                                &fReWf_pi0Sp,&fImWf_pi0Sp,
                                &fReWf_pipL,&fImWf_pipL);


           // printf("-------------(%i %i)\n%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
           //         RadBin,WhichMomBin,fMomentum, fRadius,
           //         fReWf_Kp,fImWf_Kp,
           //         fReWf_Kn,fImWf_Kn,
           //         fReWf_pipSm,fImWf_pipSm,
           //         fReWf_pi0S0,fImWf_pi0S0,
           //         fReWf_pimSp,fImWf_pimSp,
           //         fReWf_pi0L,fImWf_pi0L);
//usleep(62.5e3);
            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }
//printf(" fRadius = %f\n",fRadius);
            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            WhichBin[1] = RadBin;

            HistoWF[0][0].SetBinContent(WhichBin,(fReWf_Kp+fi*fImWf_Kp)*fRadius);
            HistoPS[0][0].SetBinContent(WhichBin,0);
            if(TYPE<10){
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_Kn+fi*fImWf_Kn)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);
            }
            else{
              HistoWF[1][0].SetBinContent(WhichBin,(fReWf_pipS0+fi*fImWf_pipS0)*fRadius);
              HistoPS[1][0].SetBinContent(WhichBin,0);

              HistoWF[2][0].SetBinContent(WhichBin,(fReWf_pi0Sp+fi*fImWf_pi0Sp)*fRadius);
              HistoPS[2][0].SetBinContent(WhichBin,0);

              HistoWF[3][0].SetBinContent(WhichBin,(fReWf_pipL+fi*fImWf_pipL)*fRadius);
              HistoPS[3][0].SetBinContent(WhichBin,0);
            }

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)){
                RadBin=0;
                WhichMomBin++;
//printf("SWITCH\n");
//usleep(500e3);
            }
        }

        fclose(InFile);
        if(WhichMomBin!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TakeThisFile;

    return Histo;
}




//TYPE 0 without Coulomb
//TYPE 1 with Coulomb (Gamow?)
DLM_Histo<complex<double>>*** InitHaidenbauerKaonPlus(const char* InputFolder, CATS& Kitty, const int& TYPE){

    //SET UP Q1Q2... Kitty
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    unsigned short NumChannels;
    unsigned short NumPwPerCh;
    //unsigned short NumFiles;
    unsigned NumMomBins = 0;
    unsigned short NumFiles = 0;
    bool* TakeThisFile = new bool [NumFiles];

    if(TYPE==0){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 1;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 30;
    }
    else if(TYPE==1){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 1;
        NumPwPerCh = 1;
        NumFiles = 1;
        NumMomBins = 30;
    }
    else{
        printf("YOU BROKE SOMETHING in InitHaidenbauerKaonPlus\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);

    Kitty.SetNumPW(0,NumPwPerCh);
    Kitty.SetChannelWeight(0,1.);
    Kitty.SetSpin(0,0);

    if(TYPE==0) Kitty.SetQ1Q2(0);
    else Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 321);
    const double Mass_p = 938.272;
    const double Mass_KaonCh = 493.677;
    Kitty.SetRedMass((Mass_p*Mass_KaonCh)/(Mass_p+Mass_KaonCh));

    const unsigned NumRadBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
//printf("NumRadBins=%u\n",NumRadBins);
    double* RadBins = new double [NumRadBins+1];
    bool* RadBinLoaded = new bool [NumRadBins+1];
    for(unsigned uRad=0; uRad<=NumRadBins; uRad++){
        //uRad-1 as we have special treatment of the zeroth bin
        //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
        //to define the bin edges
        RadBins[uRad] = RadiusMinimum+double(uRad-1)*RadiusStepSize-0.5*RadiusStepSize;
        RadBinLoaded[uRad] = false;
    }
    //we want to have the very first bin (center) exactly at zero!
    RadBins[0] = -RadBins[1];
    RadBinLoaded[0] = false;

    char** InputFileName = new char* [NumFiles];
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[256];
        strcpy(InputFileName[uFile],InputFolder);
    }

    //char* buffer = new char [512];
    if(TYPE==0){
        strcat(InputFileName[0], "ws11.dat");
    }
    else if(TYPE==1){
        strcat(InputFileName[0], "ws11p.dat");
    }
    else{
        printf("YOU BROKE SOMETHING in InitHaidenbauerKaonPlus\n");
    }


    double* MomentumBins = new double [NumMomBins+1];
    double* Momentum = new double [NumMomBins];
    MomentumBins[0] = 0;
    Momentum[0] = 10;
    MomentumBins[1] = 15;
    for(unsigned uBin=1; uBin<NumMomBins; uBin++){
        MomentumBins[uBin+1] = MomentumBins[uBin]+10;
        Momentum[uBin] = Momentum[uBin-1]+10;
    }

    //for(unsigned uBin=0; uBin<NumMomBins; uBin++){
//printf("%.2f - %.2f - %.2f\n",MomentumBins[uBin],Momentum[uBin],MomentumBins[uBin+1]);
    //}

    char* cdummy = new char [512];
    //float fMomentum;

    //const unsigned NumDLM_Histos = NumFiles>NumChannels?NumFiles:NumChannels;
    const unsigned NumDLM_Histos = NumChannels;
    //the first one is WF, second is PS
    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];

    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];
    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>>* [NumDLM_Histos];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoWF[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoWF[uHist][uPw].SetUp(2);
            HistoWF[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoWF[uHist][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uHist][uPw].Initialize();
        }
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>>* [NumDLM_Histos];
    for(unsigned uHist=0; uHist<NumDLM_Histos; uHist++){
        HistoPS[uHist] = new DLM_Histo<complex<double>> [NumPwPerCh];
        for(unsigned uPw=0; uPw<NumPwPerCh; uPw++){
            HistoPS[uHist][uPw].SetUp(1);
            HistoPS[uHist][uPw].SetUp(0,NumMomBins,MomentumBins,Momentum);
            HistoPS[uHist][uPw].Initialize();
        }
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        FILE *InFile;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }

        fseek ( InFile , 0 , SEEK_END );
        fseek ( InFile , 0 , SEEK_SET );

        if(feof(InFile)){
            printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n", InputFileName[uFile]);
            return Histo;
        }


        //unsigned LastRadBin;

        float fMomentum;
        float fRadius;
        float fReWf;
        float fImWf;
        float fReWfAs;
        float fImWfAs;
        float fReWfRot;
        float fImWfRot;
        float fPhaseShift;

        int RadBin=-1;
        int WhichMomBin=-1;
        //!---Iteration over all events---
        while(!feof(InFile)){
            //ignore empty lines between the momentum bins
            if(!fgets(cdummy, 511, InFile)||strlen(cdummy)<10) continue;
            if(WhichMomBin>=int(NumMomBins)){
                printf("\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",NumMomBins);
                printf(" Buffer reads: %s\n",cdummy);
                break;
            }

            //the first line contains info about the momentum
            if(RadBin==-1) {
                sscanf(cdummy, "%f",&fMomentum);
                RadBin++;
                WhichMomBin++;
//printf("WhichMomBin=%i\n",WhichMomBin);
                continue;
            }

//printf("READING THE NEXT LINE ---> \n");
            sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                    &fRadius, &fMomentum,
                    &fReWf,&fImWf,
                    &fReWfAs,&fImWfAs,
                    &fReWfRot,&fImWfRot,
                    &fPhaseShift);

           // printf("-------------(%i %i)\n%.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n",
           //         RadBin,WhichMomBin,fMomentum, fRadius,
           //         fReWf_Kp,fImWf_Kp,
           //         fReWf_Kn,fImWf_Kn,
           //         fReWf_pipSm,fImWf_pipSm,
           //         fReWf_pi0S0,fImWf_pi0S0,
           //         fReWf_pimSp,fImWf_pimSp,
           //         fReWf_pi0L,fImWf_pi0L);
//usleep(62.5e3);
            if(WhichMomBin<0){
                printf("\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
                continue;
            }

            unsigned WhichBin[2];
            WhichBin[0] = unsigned(WhichMomBin);
            WhichBin[1] = RadBin;

            HistoWF[0][0].SetBinContent(WhichBin,(fReWf+fi*fImWf)*fRadius);
            HistoPS[0][0].SetBinContent(WhichBin,fPhaseShift);

            RadBin++;
            //if we have are passed the last radius bin in this momentum bin => we start over.
            //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
            if(RadBin==int(NumRadBins)-1){
                RadBin=-1;
            }
        }

        fclose(InFile);
        if(WhichMomBin+1!=int(NumMomBins)){
            printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",WhichMomBin,NumMomBins);
        }
    }//uFile

    delete [] RadBinLoaded;
    delete [] MomentumBins;
    delete [] Momentum;
    delete [] cdummy;
    delete [] RadBins;
    delete [] TakeThisFile;

    return Histo;
}












DLM_Histo<complex<double>>*** Init_pd_Sebastian(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF){

    const unsigned short NumMomBins = 26;
    const double kMin = 5;
    const double kSize = 5;
    double* MomBinCenter = new double [NumMomBins];
    double* MomBins = new double [NumMomBins+1];
    MomBins[0] = 0;
    MomBins[1] = kMin+0.5*kSize;
    for(unsigned uMomBin=2; uMomBin<=NumMomBins; uMomBin++){
        MomBins[uMomBin] = MomBins[uMomBin-1]+kSize;
    }
    for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
        MomBinCenter[uMomBin] = kMin+kSize*double(uMomBin);
    }

    const unsigned short NumChannels = 2;
    char* cdummy = new char [512];
    char*** FileName;
    FileName = new char** [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        FileName[uCh] = new char* [NumMomBins];
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            FileName[uCh][uMomBin] = new char [512];
            strcpy(FileName[uCh][uMomBin],InputFolder);
            if(uCh==0) strcat(FileName[uCh][uMomBin],"Wf-pd-D/");
            else strcat(FileName[uCh][uMomBin],"Wf-pd-Q/");
            strcat(FileName[uCh][uMomBin],"Wf-pd-");
            if(uCh==0) strcat(FileName[uCh][uMomBin],"D-k");
            else strcat(FileName[uCh][uMomBin],"Q-k");
            sprintf(cdummy,"%.0f-",MomBinCenter[uMomBin]);
            strcat(FileName[uCh][uMomBin],cdummy);
            sprintf(cdummy,"L%i-",CUTOFF);
            strcat(FileName[uCh][uMomBin],cdummy);
            if(TYPE==0){strcat(FileName[uCh][uMomBin],"l0p5.dat");}
            else {strcat(FileName[uCh][uMomBin],"l1p0.dat");}
        }
    }

    //how many lines to skip at the beginning of each file
    const unsigned NumDummyLines = 17;
    unsigned CurrentLine;
    float fRadius;
    float fWfRe;
    float fWfIm;

    Kitty.SetNumChannels(NumChannels);
    Kitty.SetNumPW(0,1);
    Kitty.SetNumPW(1,1);

    Kitty.SetChannelWeight(0,1./3.);
    Kitty.SetChannelWeight(1,2./3.);

    //dummy values
    Kitty.SetSpin(0,1);
    Kitty.SetSpin(1,3);

    //btw, none of this matters
    Kitty.SetQ1Q2(1);
    Kitty.SetPdgId(2212, 100000);
    const double Mass_p = 938.272;
    const double Mass_d = 1875.613;
    Kitty.SetRedMass((Mass_p*Mass_d)/(Mass_p+Mass_d));

    const unsigned MaxNumRadBins = 256;
    unsigned NumRadBins=0;
    double* RadBins = new double [MaxNumRadBins+1];
    double* RadBinCenter = new double [MaxNumRadBins];

    FILE *InFile;
    InFile = fopen(FileName[0][0], "r");
    if(!InFile){
        printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", FileName[0][0]);
        return NULL;
    }

    CurrentLine=0;
    RadBins[0] = 0;
    while(!feof(InFile)){//line by line
        if(!fgets(cdummy, 511, InFile)) continue;//get a single line into cdummy
        if(CurrentLine<NumDummyLines){CurrentLine++; continue;}

        sscanf(cdummy, "%f %f %f",&fRadius,&fWfRe,&fWfIm);
        RadBinCenter[NumRadBins] = fRadius;
        if(NumRadBins){
            //set the bin range in between the last two bin centers
            RadBins[NumRadBins] = 0.5*(RadBinCenter[NumRadBins]+RadBinCenter[NumRadBins-1]);
        }
        NumRadBins++;

        CurrentLine++;
    }
    fclose(InFile);
    //set the upper edge of the last bin, where we just add the bin width of the last bin
    //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
    RadBins[NumRadBins] = 2.*RadBinCenter[NumRadBins-1]-RadBins[NumRadBins-1];

    DLM_Histo<complex<double>>*** Histo = new DLM_Histo<complex<double>>** [2];
    DLM_Histo<complex<double>>**& HistoWF = Histo[0];
    DLM_Histo<complex<double>>**& HistoPS = Histo[1];

    HistoWF = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uCh)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uCh); uPw++){
            //2dim, k* and radius
            HistoWF[uCh][uPw].SetUp(2);
            HistoWF[uCh][uPw].SetUp(0,NumMomBins,MomBins);
            HistoWF[uCh][uPw].SetUp(1,NumRadBins,RadBins);
            HistoWF[uCh][uPw].Initialize();
        }
    }

    HistoPS = new DLM_Histo<complex<double>>* [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh] = new DLM_Histo<complex<double>> [Kitty.GetNumPW(uCh)];
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uCh); uPw++){
            //2dim, k* and radius
            HistoPS[uCh][uPw].SetUp(2);
            HistoPS[uCh][uPw].SetUp(0,NumMomBins,MomBins);
            HistoPS[uCh][uPw].SetUp(1,NumRadBins,RadBins);
            HistoPS[uCh][uPw].Initialize();
        }
    }

    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        //for the future if we have more PWs, for now it will be a single iteration
        for(unsigned uPw=0; uPw<Kitty.GetNumPW(uCh); uPw++){
            for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
                //open the corresponding file
                InFile = fopen(FileName[uCh][uMomBin], "r");
                if(!InFile){
                    printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", FileName[uCh][uMomBin]);
                    return NULL;
                }

                CurrentLine=0;
                unsigned WhichBin[2];
                WhichBin[0] = uMomBin;
                WhichBin[1] = 0;
                while(!feof(InFile)){//line by line
                    if(!fgets(cdummy, 511, InFile)) continue;//get a single line into cdummy
                    if(CurrentLine<NumDummyLines){CurrentLine++; continue;}

                    sscanf(cdummy, "%f %f %f",&fRadius,&fWfRe,&fWfIm);
                    //figure out which k* and radius bins we are filling

                    HistoWF[uCh][uPw].SetBinContent(WhichBin,(fWfRe+fi*fWfIm)*fRadius);
                    HistoPS[uCh][uPw].SetBinContent(WhichBin,0);

                    WhichBin[1]++;
                    CurrentLine++;
                }
                fclose(InFile);
            }
        }
    }

    delete [] MomBinCenter;
    delete [] MomBins;
    delete [] cdummy;
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        for(unsigned uMomBin=0; uMomBin<NumMomBins; uMomBin++){
            delete [] FileName[uCh][uMomBin];
        }
        delete [] FileName[uCh];
    }
    delete [] FileName;
    delete [] RadBinCenter;
    delete [] RadBins;

    return Histo;
}

DLM_Histo<complex<double>>*** Init_pDminus_Haidenbauer(const char *InputFolder,
                                                       CATS &Kitty,
                                                       const int& TYPE) {

  double RadiusStepSize;
  double RadiusMinimum;
  double RadiusMaximum;

  unsigned short NumChannels;
  unsigned short NumPwPerCh;
  //unsigned short NumFiles;
  unsigned NumMomBins = 0;

  if (TYPE == 0) {
    RadiusStepSize = 0.02;
    RadiusMinimum = 0.02;
    RadiusMaximum = 10.;
    NumChannels = 2;
    NumPwPerCh = 1;
    //NumFiles = 1;
  } else {
    printf("YOU BROKE SOMETHING in Init_pDminus_Haidenbauer\n");
    return NULL;
  }

  Kitty.SetNumChannels(NumChannels);
  for (unsigned uCh = 0; uCh < NumChannels; uCh++) {
    Kitty.SetNumPW(uCh, NumPwPerCh);
    Kitty.SetChannelWeight(uCh, 1);
    Kitty.SetSpin(uCh, 0);
  }
  Kitty.SetQ1Q2(0); // Coulomb is already included in the wave function!
  Kitty.SetPdgId(2212, -411);
  const double Mass_p = 938.272;
  const double Mass_Dminus = 1869.65;
  Kitty.SetRedMass((Mass_p * Mass_Dminus) / (Mass_p + Mass_Dminus));

  Kitty.SetOnlyNumericalPw(1, true);

  //we always add 1 bin at zero, so if we have e.g. 0.1 to 0.3, these are 3 bins + 1 extra
  //NumRadiusBins = round((RadiusMaximum-RadiusMinimum)/RadiusStepSize)+1+1;
  const unsigned NumRadBins = round(
      (RadiusMaximum - RadiusMinimum) / RadiusStepSize) + 1 + 1;
  double *RadBins = new double[NumRadBins + 1];
  bool *RadBinLoaded = new bool[NumRadBins + 1];
  for (unsigned uRad = 1; uRad <= NumRadBins; uRad++) {
    //uRad-1 as we have special treatment of the zeroth bin
    //the -0.5*RadiusStepSize comes from the fact, that else we compute the bin center, while we would like
    //to define the bin edges
    RadBins[uRad] = RadiusMinimum + double(uRad - 1) * RadiusStepSize
        - 0.5 * RadiusStepSize;
    RadBinLoaded[uRad] = false;
  }
  //we want to have the very first bin (center) exactly at zero!
  RadBins[0] = -RadBins[1];
  RadBinLoaded[0] = false;

  char *InputFileName = new char[256];
  strcpy(InputFileName, InputFolder);

  if (TYPE == 0) {
    strcat(InputFileName, "kmpc.data");
  } else {
    printf("YOU BROKE SOMETHING in Init_pDminus_Haidenbauer\n");
  }

  FILE *InFile;
  InFile = fopen(InputFileName, "r");
  if (!InFile) {
    printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n",
           InputFileName);
    return NULL;
  }
  int CurrentLine = 0;

  const unsigned MaxNumMomBins = 256;
  double *MomentumBins = new double[MaxNumMomBins + 1];
  double *Momentum = new double[MaxNumMomBins];

  char *cdummy = new char[512];
  float fMomentum;

  MomentumBins[0] = 0;
  while (!feof(InFile)) {
    if (!fgets(cdummy, 511, InFile))
      continue;
    if ((CurrentLine) % int(NumRadBins) == 0) {
      sscanf(cdummy, "%f", &fMomentum);
      Momentum[NumMomBins] = fMomentum;
      if (NumMomBins) {
        //set the bin range in between the last two bin centers
        MomentumBins[NumMomBins] = 0.5
            * (Momentum[NumMomBins] + Momentum[NumMomBins - 1]);
      }
      NumMomBins++;
    }
    CurrentLine++;
  }
  fclose(InFile);
  //set the upper edge of the last bin, where we just add the bin width of the last bin
  //i.e. if we have l(low) c(center) u(up), we have that u=c+(c-l)=2c-l
  MomentumBins[NumMomBins] = 2. * Momentum[NumMomBins - 1]
      - MomentumBins[NumMomBins - 1];

  //the first one is WF, second is PS
  DLM_Histo<complex<double>> ***Histo = new DLM_Histo<complex<double>>**[2];

  DLM_Histo<complex<double>> **&HistoWF = Histo[0];
  DLM_Histo<complex<double>> **&HistoPS = Histo[1];
  HistoWF = new DLM_Histo<complex<double>>*[NumChannels];
  for (unsigned uCh = 0; uCh < NumChannels; uCh++) {
    HistoWF[uCh] = new DLM_Histo<complex<double>> [1];
    HistoWF[uCh][0].SetUp(2);
    HistoWF[uCh][0].SetUp(0, NumMomBins, MomentumBins);
    HistoWF[uCh][0].SetUp(1, NumRadBins, RadBins);
    HistoWF[uCh][0].Initialize();
  }

  HistoPS = new DLM_Histo<complex<double>>*[NumChannels];
  for (unsigned uCh = 0; uCh < NumChannels; uCh++) {
    HistoPS[uCh] = new DLM_Histo<complex<double>> [1];
    HistoPS[uCh][0].SetUp(1);
    HistoPS[uCh][0].SetUp(0, NumMomBins, MomentumBins);
    HistoPS[uCh][0].Initialize();
  }

  int WhichMomBin = -1;
  InFile = fopen(InputFileName, "r");
  if (!InFile) {
    printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n",
           InputFileName);
    return Histo;
  }

  fseek(InFile, 0, SEEK_END);
  fseek(InFile, 0, SEEK_SET);

  if (feof(InFile)) {
    printf("\033[1;31mERROR:\033[0m Trying to read past end of file %s\n",
           InputFileName);
    return Histo;
  }

  float fRadius;
  float fMom;

  float fReWf_pDm_pDm;
  float fImWf_pDm_pDm;
  float fReWf_nDo_pDm;
  float fImWf_nDo_pDm;

  int RadBin = -1;
  //!---Iteration over all events---
  while (!feof(InFile)) {
    if (!fgets(cdummy, 511, InFile))
      continue;
    if (WhichMomBin >= int(NumMomBins)) {
      printf(
          "\033[1;31mERROR:\033[0m Trying to read more momentum bins than set up (%u)!\n",
          NumMomBins);
      printf(" Buffer reads: %s\n", cdummy);
      break;
    }

    //the first line contains info about the momentum
    if (RadBin == -1) {
      sscanf(cdummy, "%f", &fMomentum);

      RadBin++;
      WhichMomBin++;
      continue;
    }

    sscanf(cdummy, " %f %f %f %f %f %f", &fRadius, &fMom, &fReWf_pDm_pDm,
           &fImWf_pDm_pDm, &fReWf_nDo_pDm, &fImWf_nDo_pDm);

    if (WhichMomBin < 0) {
      printf(
          "\033[1;33mWARNING:\033[0m WhichMomBin==-1, possible bug, please contact the developers!\n");
      continue;
    }

    unsigned WhichBin[2];
    WhichBin[0] = unsigned(WhichMomBin);
    //we fill up the radius bins with and offset of 1, due to the special zeroth bin
    WhichBin[1] = RadBin + 1;

    HistoWF[0][0].SetBinContent(
        WhichBin, (fReWf_pDm_pDm + fi * fImWf_pDm_pDm) * fRadius);
    HistoWF[1][0].SetBinContent(
        WhichBin, (fReWf_nDo_pDm + fi * fImWf_nDo_pDm) * fRadius);

    HistoPS[0][0].SetBinContent(WhichBin, 0);
    HistoPS[1][0].SetBinContent(WhichBin, 0);

    RadBin++;
    //if we have are passed the last radius bin in this momentum bin => we start over.
    //the -1 is due to the special case of the zeroth bin (which we do NOT read from the file)
    if (RadBin == int(NumRadBins) - 1) {
      RadBin = -1;
    }
  }

  fclose(InFile);
  if (WhichMomBin + 1 != int(NumMomBins)) {
    printf("\033[1;31mERROR:\033[0m WhichMomBin!=NumMomBins (%u vs %u)\n",
           WhichMomBin, NumMomBins);
  }
  delete[] RadBinLoaded;
  delete[] MomentumBins;
  delete[] Momentum;
  delete[] cdummy;
  delete[] RadBins;

  return Histo;
}

DLM_Histo<complex<double>>*** Init_pDminus_Haidenbauer(const char *InputFolder,
                                                       CATS *Kitty,
                                                       const int& TYPE) {
  return Init_pDminus_Haidenbauer(InputFolder, *Kitty, TYPE);
}


void CleanUpWfHisto(const unsigned short& NumChannels, DLM_Histo<complex<double>>***& Histo){
    for(unsigned short usCh=0; usCh<NumChannels; usCh++){
        if(Histo&&Histo[0]&&Histo[0][usCh]){delete[]Histo[0][usCh];Histo[0][usCh]=NULL;}
        if(Histo&&Histo[1]&&Histo[1][usCh]){delete[]Histo[1][usCh];Histo[1][usCh]=NULL;}
    }
    if(Histo&&Histo[0]){delete[]Histo[0];Histo[0]=NULL;}
    if(Histo&&Histo[1]){delete[]Histo[1];Histo[1]=NULL;}
    if(Histo){delete[]Histo;Histo=NULL;}
}
void CleanUpWfHisto(const CATS& Kitty, DLM_Histo<complex<double>>***& Histo){
    CleanUpWfHisto(Kitty.GetNumChannels(),Histo);
}
