#include <stdio.h>
#include "string.h"

//#include "DLM_WfModel.h"
#include "CATS.h"

//#include <unistd.h>

const complex<float>i(0,1);

//TYPE = 00 is LO
//TYPE = 01 is LO+Coupling
//TYPE = 10 is NLO (s-wave)
//TYPE = 11 is NLO+Coupling
//TYPE = 12 is NLO (s and p-waves)
//TYPE
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
    unsigned short NumFiles;
    unsigned NumMomBins = 0;

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

    enum HaideFiles {f1S0, f3S1, f1P1, f3P0, f3P1, f3P2};
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
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_LNLN_SS+i*fImWf_LNLN_SS)*fRadius);
                HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_SNLN_SS+i*fImWf_SNLN_SS)*fRadius);

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
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_LNLN_SS+i*fImWf_LNLN_SS)*fRadius);
                HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_SNLN_SS+i*fImWf_SNLN_SS)*fRadius);

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

    return Histo;
}
DLM_Histo<complex<double>>*** Init_pL_Haidenbauer(const char* InputFolder, CATS* Kitty, const int& TYPE, const int& CUTOFF){
    return Init_pL_Haidenbauer(InputFolder,*Kitty,TYPE,CUTOFF);
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

            HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf_S0p_S0p+i*fImWf_S0p_S0p)*fRadius);
            HistoWF[uFile+2][0].SetBinContent(WhichBin,(fReWf_Lp_S0p+i*fImWf_Lp_S0p)*fRadius);
            HistoWF[uFile+4][0].SetBinContent(WhichBin,(fReWf_Sn_S0p+i*fImWf_Sn_S0p)*fRadius);

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
                HistoWF[uFile][0].SetBinContent(WhichBin,(fReWf+i*fImWf)*fRadius);
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


DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS& Kitty, const int& TYPE=0){

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
                    HistoWF[WhichIsospin*NumSpinChannels+WhichSpin]->SetBinContent(WhichBin,Values[uCol]);
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
DLM_Histo<complex<double>>*** Init_pXi_ESC16_IS(const char* InputFolder, CATS* Kitty, const int& TYPE=0){
    return Init_pXi_ESC16_IS(InputFolder,*Kitty,TYPE);
}

//TYPE==0, some old WF,
//TYPE==1, function from 2nd April 2019
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS& Kitty, const int& TYPE=0){

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
                    HistoWF[WhichSpin]->SetBinContent(WhichBin,Values[uCol]);
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
DLM_Histo<complex<double>>*** Init_pXi_ESC16_Iavg_Coulomb(const char* InputFolder, CATS* Kitty, const int& TYPE=0){
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
    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        InputFileName[uFile] = new char[512];
        strcpy(InputFileName[uFile],InputFolder);
    }

    if(TYPE==0){
        strcat(InputFileName[fS0], "ynwave.si0-prot.lsj=000");
        strcat(InputFileName[fS1], "ynwave.si0-prot.lsj=011");
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
                    HistoWF[WhichSpin]->SetBinContent(WhichBin,Values[uCol]);
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
