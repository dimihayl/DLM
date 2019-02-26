#include <stdio.h>
#include "string.h"

//#include "DLM_WfModel.h"
#include "CATS.h"

const complex<float>i(0,1);

//TYPE = 0 is LO
//TYPE = 1 is NLO
//TYPE = 2 is NLO+Coupling
//for TYPE==2 we make it so, that Kitty has 4 channels
//ch0 and ch1 are the LN->LN 1S0 and 3S1 (as per usual)
//ch2 and ch2 are the SN->LN in the S-wave (1S0 and 3S1), however we ONLY add to the total correlation function the s-wave
//to obtain the total correlation we need to add 1/4*ch0+3/4*ch1+1/4*ch2+3/4*ch3. As the ch2 and ch3 are just "corrections" to
//the primary channels (0 and 1) it is completely fine that our weights do not sum up to 1
DLM_Histo<complex<double>>** Init_pL_Haidenbauer(const char* InputFolder, CATS& Kitty, const int& TYPE, const int& CUTOFF){
    double RadiusStepSize;
    double RadiusMinimum;
    double RadiusMaximum;
    //unsigned NumRadiusBins;

    if(CUTOFF!=500&&CUTOFF!=600){
        printf("Problem with the CUTOFF in InitInterpolHaide_pL\n");
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
//! p-waves not yet implemented!
    else if(TYPE==1&&CUTOFF==600){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.2;
        RadiusMaximum = 16.1;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 40;
    }
    else if(TYPE==1&&CUTOFF==500){
        RadiusStepSize = 0.1;
        RadiusMinimum = 0.1;
        RadiusMaximum = 10.;
        NumChannels = 2;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 40;
    }
    else if(TYPE==2){
        RadiusStepSize = 0.02;
        RadiusMinimum = 0.02;
        RadiusMaximum = 10.;
        NumChannels = 4;
        NumPwPerCh = 1;
        NumFiles = 2;
        //NumMomBins = 43;
    }
    else{
        printf("YOU BROKE SOMETHING in InitInterpolHaide_pL\n");
        return NULL;
    }

    Kitty.SetNumChannels(NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        Kitty.SetNumPW(uCh,NumPwPerCh);
        Kitty.SetChannelWeight(uCh,uCh%2==0?0.25:0.75);
        Kitty.SetSpin(uCh,uCh%2==0?0:1);
    }
    Kitty.SetQ1Q2(0);
    Kitty.SetPdgId(2212, 3122);
    const double Mass_p = 938.272;
    const double Mass_L = 1115.683;
    Kitty.SetRedMass((Mass_p*Mass_L)/(Mass_p+Mass_L));
    if(TYPE==2){
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

    if(TYPE==1 && CUTOFF==600){
        strcat(InputFileName[f1S0], "N1s0.data");
        //strcat(InputFileName[f1P1], "N1p1.data");
        strcat(InputFileName[f3S1], "N3s1.data");
        //strcat(InputFileName[f3P0], "N3p0.data");
        //strcat(InputFileName[f3P1], "N3p1.data");
        //strcat(InputFileName[f3P2], "N3p2.data");
    }
    else if(TYPE==1 && CUTOFF==500){
        strcat(InputFileName[f1S0], "Y1s05.data");
        strcat(InputFileName[f3S1], "Y3s15.data");
    }
    else if(TYPE==0 && CUTOFF==600){
        strcat(InputFileName[f1S0], "W1s06.data");
        strcat(InputFileName[f3S1], "W3s16.data");
    }
    else if(TYPE==2 && CUTOFF==600){
        strcat(InputFileName[f1S0], "fort.10");
        strcat(InputFileName[f3S1], "fort.11");
    }
    else{
        printf("YOU BROKE SOMETHING in InitInterpolHaide_pL\n");
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

    //the first one is WF, second is PS
    DLM_Histo<complex<double>>** Histo = new DLM_Histo<complex<double>>* [2];
    DLM_Histo<complex<double>>*& HistoWF = Histo[0];
    DLM_Histo<complex<double>>*& HistoPS = Histo[1];

    //DLM_Histo<complex<double>> HistoWF[NumChannels];
    HistoWF = new DLM_Histo<complex<double>> [NumChannels];
//printf("NumChannels = %u\n",NumChannels);
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoWF[uCh].SetUp(2);
        HistoWF[uCh].SetUp(0,NumMomBins,MomentumBins);
        HistoWF[uCh].SetUp(1,NumRadBins,RadBins);
        HistoWF[uCh].Initialize();
    }

    //DLM_Histo<complex<double>> HistoPS[NumChannels];
    HistoPS = new DLM_Histo<complex<double>> [NumChannels];
    for(unsigned uCh=0; uCh<NumChannels; uCh++){
        HistoPS[uCh].SetUp(1);
        HistoPS[uCh].SetUp(0,NumMomBins,MomentumBins);
        HistoPS[uCh].Initialize();
    }

    for(unsigned uFile=0; uFile<NumFiles; uFile++){
        int WhichMomBin=-1;
        InFile = fopen(InputFileName[uFile], "r");
        if(!InFile){
            printf("\033[1;31mERROR:\033[0m The file\033[0m %s cannot be opened!\n", InputFileName[uFile]);
            return Histo;
        }
        //we have the p-waves only for NLO with cutoff of 600 MeV
        if(uFile!=f1S0 && uFile!=f3S1 && (TYPE!=1 || CUTOFF!=600)) {printf("Possible bug in InitInterpolHaide_pL\n"); continue;}

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

            if(TYPE==0||TYPE==1){
                sscanf(cdummy, " %f %f %f %f %f %f %f %f %f",
                       &fRadius,&fMomentum,&fReWf,&fImWf,&fReAsWf,&fImAsWf,&fCatsWf,&fDummy,&fPhaseShift);
            }
            else if(TYPE==2){
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
                    printf("Oh man... big bug in InitInterpolHaide_pL.\n");
                }
            }
            else{
                printf("WTF from InitInterpolHaide_pL\n");
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
                HistoWF[uFile].SetBinContent(WhichBin,fCatsWf*fRadius);
                HistoPS[uFile].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==1 && CUTOFF==600){
                HistoWF[uFile].SetBinContent(WhichBin,fCatsWf*fRadius);
                //HistoWF[uFile].SetBinContent(WhichBin,(fReWf+i*fImWf)*fRadius);
                HistoPS[uFile].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==1 && CUTOFF==500){
                HistoWF[uFile].SetBinContent(WhichBin,fCatsWf*fRadius);
                HistoPS[uFile].SetBinContent(WhichBin,fPhaseShift);
            }
            else if(TYPE==2 && CUTOFF==600){
                HistoWF[uFile].SetBinContent(WhichBin,(fReWf_LNLN_SS+i*fImWf_LNLN_SS)*fRadius);
                HistoWF[uFile+2].SetBinContent(WhichBin,(fReWf_SNLN_SS+i*fImWf_SNLN_SS)*fRadius);

                //HistoWF[uFile].SetBinContent(WhichBin,(fReWf_LNLN_DS+i*fImWf_LNLN_DS)*fRadius);
                //HistoWF[uFile+2].SetBinContent(WhichBin,(fReWf_SNLN_DS+i*fImWf_SNLN_DS)*fRadius);
//&fReWf_LNLN_DS,&fImWf_LNLN_DS,&fReWf_SNLN_DS,&fImWf_SNLN_DS

//if(fMomentum>280 && fMomentum<290){
//printf("uFile=%u; fMomentum=%.1f; fReWf=%.4f(%.4f); fImWf=%.4f(%.4f)\n",uFile,fMomentum,
//       fReWf_LNLN_SS,fReWf_SNLN_SS,fImWf_LNLN_SS,fImWf_SNLN_SS);
//}
                HistoPS[uFile].SetBinContent(WhichBin,0);
                HistoPS[uFile+2].SetBinContent(WhichBin,0);
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
