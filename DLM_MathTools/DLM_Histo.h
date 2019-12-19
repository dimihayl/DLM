
#ifndef DLM_HISTO_H
#define DLM_HISTO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "DLM_CppTools.h"
//#include <stdint.h>
//#include <complex>

template <class Type> class DLM_Histo1D{
public:
    DLM_Histo1D(const unsigned& numbin, const Type* bins):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        BinCenter = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        BinCenter = new Type [NumBins];
        for(unsigned uBin=0; uBin<=NumBins; uBin++){
            BinRange[uBin] = bins[uBin];
        }
        //for(unsigned uBin=0; uBin<NumBins; uBin++){
            //BinValue[uBin] = GetBinCenter(uBin);
        //}
    }
    DLM_Histo1D(const unsigned& numbin, const Type& xmin, const Type& xmax):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        BinCenter = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        BinCenter = new Type [NumBins];
        if(NumBins==1){
            BinRange[0] = xmin; BinRange[1] = xmax;
            //BinValue[0] = (xmin+xmax)*0.5;
            BinCenter[0] = 1.331e128;
        }
        else{
            Type BinWidth = (xmax-xmin)/Type(NumBins);
//printf("NumBins=%u; BinWidth=%f;\n",NumBins,BinWidth);
            for(unsigned uBin=0; uBin<=NumBins; uBin++){
                BinRange[uBin] = xmin + Type(uBin)*BinWidth;
                if(uBin!=NumBins) BinCenter[uBin] = 1.331e128;
//printf(" uBin=%u -> BR=%f\n",uBin, BinRange[uBin]);
            }
            //for(unsigned uBin=0; uBin<NumBins; uBin++){
            //    BinValue[uBin] = GetBinCenter(uBin);
            //}
        }

    }
    DLM_Histo1D(const DLM_Histo1D& other):NumBins(other.NumBins){
        BinRange = NULL;
        BinValue = NULL;
        BinCenter = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        BinCenter = new Type [NumBins];
        operator=(other);
    }
    ~DLM_Histo1D(){
        if(BinRange) {delete [] BinRange; BinRange=NULL;}
        if(BinValue) {delete [] BinValue; BinValue=NULL;}
        if(BinCenter) {delete [] BinCenter; BinCenter=NULL;}
    }

    void SetBinContent(const unsigned& WhichBin, const Type& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]=Val;
    }
    void SetBinCenter(const unsigned& WhichBin, const Type& Val){
        if(WhichBin>=NumBins) return;
        BinCenter[WhichBin]=Val;
    }
    void Add(const unsigned& WhichBin, const Type& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]+=Val;
    }
    void SetBinAt(const Type& xVal, const Type& Val){
        unsigned WhichBin = GetBin(xVal);
        SetBinContent(WhichBin, Val);
    }
    void AddAt(const Type& xVal, const Type& Val){
        unsigned WhichBin = GetBin(xVal);
        Add(WhichBin, Val);
    }

    unsigned GetBin(const Type& xVal) const{
        if(NumBins<=1) return 0;
        unsigned WhichBin=(NumBins+1)/2;
        unsigned BinMod=4;
        unsigned BinStep;
        //makes sure that the xVal is in BinRange. If not, the returned value is either
        //NumBins or NumBins+1, depending on if we have an underflow or overflow
        if(xVal<BinRange[0]) return NumBins;
        if(xVal>BinRange[NumBins]) return NumBins+1;
        while(true){
            if(BinRange[WhichBin]<=xVal && BinRange[WhichBin+1]>=xVal){
                return WhichBin;
            }
    //!check the logic of this <=, I made it so, since else we crashed on limit values, before that it was <
            else if(xVal<=BinRange[WhichBin]){
                BinStep = (NumBins+1)/BinMod;
                WhichBin -= BinStep?BinStep:1;
                BinMod *= 2;
            }
            else{
                BinStep = (NumBins+1)/BinMod;
                WhichBin += BinStep?BinStep:1;
                BinMod *= 2;
            }
        }
    }
    Type GetBinCenter(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        if(BinCenter[WhichBin]!=1.331e128) return BinCenter[WhichBin];
        return 0.5*(BinRange[WhichBin]+BinRange[WhichBin+1]);
    }
    Type GetBinContent(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return BinValue[WhichBin];
    }
    Type GetBinLowEdge(const unsigned& WhichBin) const{
        //i.e. GetBinLowEdge(NumBins) will return the most upper limit
        if(WhichBin>=NumBins+1) return 0;
        return BinRange[WhichBin];
    }
    Type GetBinUpEdge(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return BinRange[WhichBin+1];
    }
    Type GetBinWidth(const unsigned& WhichBin) const{
        if(WhichBin>=NumBins) return 0;
        return (BinRange[WhichBin+1]-BinRange[WhichBin]);
    }
    Type Eval(const Type& xVal) const{
        if(xVal<BinRange[0] || xVal>BinRange[NumBins]) return 0;
        if(NumBins==1) return BinValue[0];
        unsigned WhichBin = GetBin(xVal);

        Type Value[3];
        Value[0] = WhichBin?GetBinCenter(WhichBin-1):-1;
        Value[1] = GetBinCenter(WhichBin);
        Value[2] = WhichBin<(NumBins-1)?GetBinCenter(WhichBin+1):-1;

        Type* InterpolRange;
        const Type* FunRange;

        if(Value[0]==-1){
            InterpolRange = &Value[1];
            FunRange = &BinValue[WhichBin];
        }
        else if(Value[2]==-1){
            InterpolRange = &Value[0];
            FunRange = &BinValue[WhichBin-1];
        }
        else if(xVal<Value[1]){
            InterpolRange = &Value[0];
            FunRange = &BinValue[WhichBin-1];
        }
        else if(Value[1]<xVal){
            InterpolRange = &Value[1];
            FunRange = &BinValue[WhichBin];
        }
        else{//Value[1]==xVal
            return BinValue[WhichBin];
        }

        if(InterpolRange[1]-InterpolRange[0]){
            return (FunRange[1]*(xVal-InterpolRange[0])-
                    FunRange[0]*(xVal-InterpolRange[1]))/
                    (InterpolRange[1]-InterpolRange[0]);
        }
        else{
            printf("WTF!!! at bin %u\n",WhichBin);
            return 1e6;
        }
    }

    unsigned GetNbins() const{
        return NumBins;
    }

    bool Copy(const DLM_Histo1D& other, const Type& Scale){
        if(NumBins!=other.NumBins) return false;
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            BinRange[uBin] = other.BinRange[uBin];
            BinValue[uBin] = other.BinValue[uBin]*Scale;
            BinCenter[uBin] = other.BinCenter[uBin];
        }
        BinRange[NumBins] = other.BinRange[NumBins];
/*
        for(unsigned uBin=0; uBin<NumBins; uBin++){
            if(other.BinValue[uBin]>0.95452 && other.BinValue[uBin]<0.95453){
                printf("I found it: the copy is %f\n",BinValue[uBin]);
            }
            printf("%u: OBV=%f; BV=%f; Scale=%f\n", uBin, other.BinValue[uBin], BinValue[uBin], Scale);
        }
*/
        return true;
    }
    bool Add(const DLM_Histo1D& other, const Type& Scale){
        if(NumBins!=other.NumBins) return false;
        for(unsigned uBin=0; uBin<=NumBins; uBin++) if(BinRange[uBin]!=other.BinRange[uBin]) return false;
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] += other.BinValue[uBin]*Scale;
        return true;
    }
    void Add(const Type& Value){
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] += Value;
    }
    void Scale(const Type& Value){
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] *= Value;
    }
    void ScaleToBinWidth(){
        for(unsigned uBin=0; uBin<NumBins; uBin++) BinValue[uBin] /= GetBinWidth(uBin);
    }

    bool operator=(const DLM_Histo1D& other){
        return Copy(other, 1);
    }
    bool operator+=(const DLM_Histo1D& other){
        return Add(other, 1);
    }

protected:

    DLM_Histo1D(const unsigned& numbin):NumBins(numbin){
        BinRange = NULL;
        BinValue = NULL;
        BinCenter = NULL;
        if(!NumBins){
            return;
        }
        BinRange = new Type [NumBins+1];
        BinValue = new Type [NumBins];
        BinCenter = new Type [NumBins];
    }

    const unsigned NumBins;
    Type* BinRange;
    Type* BinValue;
    Type* BinCenter;
};



















template <class Type> class DLM_Histo{
public:
    DLM_Histo(){
//printf("DLM_Histo() %p\n",this);
        BinRange = NULL;
        BinValue = NULL;
        BinError = NULL;
        BinCenter = NULL;
        NumBins = NULL;
        TotNumBins=0;
        PER = NULL;
        Initialized=false;
        Dim = 0;
//printf("BinValue=%p\n",BinValue);
        //xValue = NULL;
        //xBin = NULL;
        //FunctionValue = NULL;
        //Normalization = NULL;
        //TotalNorm=0;
    }
    DLM_Histo(const DLM_Histo& other):DLM_Histo(){
//printf("DLM_Histo(const DLM_Histo& other)\n");
//printf("BinValue=%p\n",BinValue);
        operator=(other);
    }
    ~DLM_Histo(){
//printf("~DLM_Histo()\n");
        CleanUp();
    }
    //Warning, using those functions will delete all previous contents of the bins
    void SetUp(const unsigned short& dim){
        if(!dim){
            printf("\033[1;31mERROR:\033[0m DLM_Histo cannot have 0 dimensions!\n");
            return;
        }
        if(Dim==dim) return;

        CleanUp();
        Dim = dim;
        BinRange = new double* [Dim];
        BinValue = NULL;
        BinError = NULL;
        BinCenter = new double* [Dim];
        NumBins = new unsigned [Dim];
        //xValue = new Type* [Dim];
        //xBin = new unsigned* [Dim];
        //FunctionValue = new Type* [Dim];
        //Normalization = new Type* [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            BinRange[sDim]=NULL;
            //BinValue[sDim]=0;
            //BinError[sDim]=0;
            BinCenter[sDim]=NULL;
            NumBins[sDim]=0;
            //xValue[sDim] = new Type[2];
            //xBin[sDim] = new unsigned[2];
            //FunctionValue[sDim] = new Type[2];
            //Normalization[sDim] = new Type[2];
        }
        InitPER();
    }
    void SetUp(const unsigned short& sDim, const unsigned& numbins, const double* bins, const double* bincenter=NULL){
        if(sDim>=Dim) return;
        if(!numbins){
            printf("\033[1;31mERROR:\033[0m DLM_Histo cannot have 0 bins (holds true for each axis)!\n");
            return;
        }
        for(unsigned uBin=0; uBin<numbins; uBin++){
            if(bins[uBin+1]<=bins[uBin]){
//printf("\033[1;31msDim=%u; uBin=%u; %f --> %f\033[0m\n",sDim,uBin,bins[uBin],bins[uBin+1]);
                printf("\033[1;31mERROR:\033[0m DLM_Histo: The bin ranges should be in ascending order and a bin-width of 0 is not allowed!\n");
                return;
            }
//else
//printf("sDim=%u; %f --> %f\n",sDim,bins[uBin],bins[uBin+1]);
        }
        if(NumBins[sDim]!=numbins){
            CleanUp(sDim);
        }
        if(!BinRange[sDim]) BinRange[sDim] = new double[numbins+1];
        if(!BinCenter[sDim]) BinCenter[sDim] = new double[numbins];
        NumBins[sDim] = numbins;
        for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
            BinRange[sDim][uBin] = bins[uBin];
            if(uBin){
                if(bincenter){
                    BinCenter[sDim][uBin-1] = bincenter[uBin-1];
                }
                else{
                    BinCenter[sDim][uBin-1] = (bins[uBin-1]+bins[uBin])*0.5;
                }
            }

        }
    }
    void SetUp(const unsigned short& sDim, const unsigned& numbins, const double& xmin, const double& xmax){
        if(sDim>=Dim) return;
        if(!numbins){
            printf("\033[1;31mERROR:\033[0m DLM_Histo cannot have 0 bins (holds true for each axis)!\n");
            return;
        }
        if(xmin>=xmax){
            printf("\033[1;31mERROR:\033[0m DLM_Histo cannot accept xmin>=xmax!\n");
            return;
        }
        if(NumBins[sDim]!=numbins){
            CleanUp(sDim);
        }
        if(!BinRange[sDim]) BinRange[sDim] = new double[numbins+1];
        if(!BinCenter[sDim]) BinCenter[sDim] = new double[numbins];
        NumBins[sDim] = numbins;
        if(NumBins[sDim]==1){
            BinRange[sDim][0] = xmin;
            BinRange[sDim][1] = xmax;
            BinCenter[sDim][0] = (xmin+xmax)*0.5;
        }
        else{
            double BinWidth = (xmax-xmin)/double(NumBins[sDim]);
            for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
                BinRange[sDim][uBin] = xmin + double(uBin)*BinWidth;
                if(uBin) BinCenter[sDim][uBin-1] = (BinRange[sDim][uBin-1]+BinRange[sDim][uBin])*0.5;
            }
        }
    }
    bool Initialize(const bool& ZeroElements=true){
//printf("Initialize 1\n");
        if(!Dim||!NumBins||!BinRange||!BinCenter) return false;
        TotNumBins=1;
//printf("Initialize 2\n");
        for(unsigned short sDim=0; sDim<Dim; sDim++){
//printf("sDim=%u; NumBins[sDim]=%u; BinRange[sDim]=%p; BinCenter[sDim]=%p\n",sDim,NumBins[sDim],BinRange[sDim],BinCenter[sDim]);
            if(!NumBins[sDim]||!BinRange[sDim]||!BinCenter[sDim]) {TotNumBins=0;return false;}
            TotNumBins*=NumBins[sDim];
        }
        if(!BinValue){delete [] BinValue; BinValue=NULL;}
        if(!BinError){delete [] BinError; BinError=NULL;}
        BinValue = new Type[TotNumBins+2];
        BinError = new Type[TotNumBins];
        if(ZeroElements){
            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                BinValue[uBin]=0;
                BinError[uBin]=0;
            }
            BinValue[TotNumBins]=0;
            BinValue[TotNumBins+1]=0;
        }
        Initialized = true;
//printf("Initialize 3\n");
        return true;
    }
    unsigned GetTotBin(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        unsigned TotBin=0;
        unsigned Base=1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            if(WhichBin[sDim]>=NumBins[sDim]) return TotNumBins;
            TotBin+=Base*WhichBin[sDim];
            Base*=NumBins[sDim];
        }
        return TotBin;
    }
    //BinIdPerAxis[sDim] -> saves for each dim the bin number
    void GetBinCoordinates(const unsigned& WhichTotBin, unsigned* BinIdPerAxis) const{
        if(!Initialized) {InitWarning(); return;}
        //unsigned Division;
        //unsigned Reminder=WhichTotBin;
        ldiv_t divresult;
        divresult.rem = WhichTotBin;
        unsigned short isDim;
        unsigned Base=TotNumBins;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            isDim = Dim-sDim-1;
            Base/=NumBins[isDim];
            divresult = ldiv(divresult.rem,Base);
            BinIdPerAxis[isDim] = divresult.quot;
        }
    }
    bool Copy(const DLM_Histo& other){
        //if(!Initialized) {InitWarning(); return false;}
        return operator=(other);
    }
    bool AddHisto(const DLM_Histo& other, const bool witherror=true, const Type& scale=1){
        if(!Initialized) {InitWarning(); return false;}
        if(!SameStructure(other)||!Initialized) return false;
        for(unsigned uBin=0; uBin<TotNumBins; uBin++) BinValue[uBin] += scale*other.BinValue[uBin];
        if(witherror)
            for(unsigned uBin=0; uBin<TotNumBins; uBin++)
                BinError[uBin] = sqrt(BinError[uBin]*BinError[uBin]+scale*scale*other.BinError[uBin]*other.BinError[uBin]);
        return true;
    }
    bool AddWeightedHisto(const DLM_Histo& other, const Type& scale=1){
        if(!Initialized) {InitWarning(); return false;}
        if(!Initialized) return false;
        unsigned* BinId = new unsigned [Dim];
        double* xVal = new double [Dim];
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
            }
            Type VAR = 1./(1./BinError[uBin]/BinError[uBin] + 1./other.Eval(xVal,1)/other.Eval(xVal,1));
            BinValue[uBin] = BinValue[uBin]/BinError[uBin]/BinError[uBin] + other.Eval(xVal)/other.Eval(xVal,1)/other.Eval(xVal,1);
            BinValue[uBin] *= VAR;
            BinError[uBin] = sqrt(VAR);
        }
        delete [] BinId;
        delete [] xVal;
        return true;
    }
    bool MultiplyHisto(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!SameStructure(other)||!Initialized) return false;
        for(unsigned uBin=0; uBin<TotNumBins; uBin++) BinValue[uBin] *= other.BinValue[uBin];
        if(witherror)
            for(unsigned uBin=0; uBin<TotNumBins; uBin++)
                BinError[uBin] = sqrt(  pow(BinError[uBin]*other.BinError[uBin],2.)+
                                        pow(BinError[uBin]*other.BinValue[uBin],2.)+
                                        pow(BinValue[uBin]*other.BinError[uBin],2.));
        return true;
    }
    /*
    bool DivideHisto(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!SameStructure(other)||!other.Initialized) return false;
        if(witherror)
            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                if(!other.BinValue[uBin]) BinError[uBin] = 1e128;
                else BinError[uBin] = pow(other.BinValue[uBin],-2.)*
                                            sqrt(   pow(BinError[uBin]*other.BinError[uBin],2.)+
                                                    pow(BinError[uBin]*other.BinValue[uBin],2.)+
                                                    pow(BinValue[uBin]*other.BinError[uBin],2.));
            }
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            if(!other.BinValue[uBin]) BinValue[uBin] = 1e128;
            else BinValue[uBin] /= other.BinValue[uBin];
        }
        return true;
    }
    */
    bool DivideHisto(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        unsigned* BinId = new unsigned [Dim];
        double* xVal = new double [Dim];
        if(witherror)
            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                GetBinCoordinates(uBin,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                }
                if(!other.Eval(xVal)) BinError[uBin] = 1e128;
                else BinError[uBin] = pow(other.Eval(xVal),-2.)*
                                            sqrt(   pow(BinError[uBin]*other.Eval(xVal,true),2.)+
                                                    pow(BinError[uBin]*other.Eval(xVal),2.)+
                                                    pow(BinValue[uBin]*other.Eval(xVal,true),2.));
            }
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
            }
            if(!other.Eval(xVal)) BinValue[uBin] = 1e128;
            else BinValue[uBin] /= other.Eval(xVal);
        }
        delete [] BinId;
        delete [] xVal;
        return true;
    }
    void AddToAll(const Type& value, const Type& error=0){
        if(!Initialized) {InitWarning(); return;}
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin]+=value;
            BinError[uBin]=Type(sqrt(double(BinError[uBin]*BinError[uBin])+double(error*error)));
        }
    }
    void Scale(const Type& scale){
        if(!Initialized) {InitWarning(); return;}
        if(!scale) return;
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin] *= scale;
            BinError[uBin] *= scale;
        }
        BinValue[TotNumBins] *= scale;
        BinValue[TotNumBins+1] *= scale;
    }
    void ScaleToBinSize(){
        if(!Initialized) {InitWarning(); return;}
        Type BinSize;
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinSize = GetBinSize(uBin);
            BinValue[uBin] /= BinSize;
            BinError[uBin] /= BinSize;
        }
    }
    void RescaleAxis(const unsigned short& sDim, const double& scale, const bool& binwidth){
        if(!Initialized) {InitWarning(); return;}
        if(sDim>=Dim) return;
        if(!scale) return;
        for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
            BinRange[sDim][uBin] *= scale;
            BinCenter[sDim][uBin] *= scale;
        }
        BinRange[sDim][NumBins[sDim]] *= scale;
        if(binwidth){
            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                BinValue[uBin] /= scale;
                BinError[uBin] /= scale;
            }
        }
    }
    //xMin and xMax have size Dim, and each entry represents the min/max value to which the
    //corresponding dimension should be considered
    Type Integral(const double* xMin, const double* xMax){
        Type RESULT = 0;
        if(!Initialized) {InitWarning(); return RESULT;}
        unsigned* BinId = new unsigned [Dim];
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            double FractionInside=1;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                //in case both xMin and xMax are outside of the bin
                if(BinRange[sDim][BinId[sDim]]>=xMin[sDim] && BinRange[sDim][BinId[sDim]+1]<xMax[sDim])
                {continue;}
                //in case both xMin and xMax are within the bin
                else if(BinRange[sDim][BinId[sDim]]<xMin[sDim] && BinRange[sDim][BinId[sDim]+1]>=xMax[sDim]){
                    FractionInside *= ( (xMax[sDim]-xMin[sDim])/GetBinSize(sDim,BinId[sDim]) );
                }
                //in case xMax is within the bin
                else if(BinRange[sDim][BinId[sDim]]>=xMin[sDim] && BinRange[sDim][BinId[sDim]+1]>=xMax[sDim]){
                    FractionInside *= ( (xMax[sDim]-BinRange[sDim][BinId[sDim]])/GetBinSize(sDim,BinId[sDim]) );
                }
                //in case xMin is within the bin
                else if(BinRange[sDim][BinId[sDim]]<xMin[sDim] && BinRange[sDim][BinId[sDim]+1]<xMin[sDim]){
                    FractionInside *= ( (BinRange[sDim][BinId[sDim]+1]-xMin[sDim])/GetBinSize(sDim,BinId[sDim]) );
                }
                else{
                    FractionInside = 0;
                    break;
                }
            }
            RESULT+=FractionInside*BinValue[uBin];
        }
        delete [] BinId;
        return RESULT;
    }
    unsigned short GetDim() const{
        return Dim;
    }
    unsigned GetNbins() const{
        if(!Initialized) {InitWarning(); return 0;}
        return TotNumBins;
    }
    unsigned GetNbins(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        return NumBins?NumBins[sDim]:0;
    }

    void SetBinContent(const unsigned& WhichTotBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;
        BinValue[WhichTotBin]=Val;
    }
     void SetBinError(const unsigned& WhichTotBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;
        BinError[WhichTotBin]=Val;
    }
    void SetBinContentAll(const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin]=Val;
        }
    }
     void SetBinErrorAll(const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinError[uBin]=Val;
        }
    }

    void Add(const unsigned& WhichTotBin, const Type& Val, const Type& Err=0){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;
        BinValue[WhichTotBin]+=Val;
        BinError[WhichTotBin]=Type(sqrt(double(BinError[WhichTotBin]*BinError[WhichTotBin])+double(Err*Err)));
    }
    void AddAt(const double* AxisValue, const Type& Val=1){
        if(!Initialized) {InitWarning(); return;}
        unsigned* WhichBin = new unsigned [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            WhichBin[sDim] = GetBin(sDim,AxisValue[sDim]);
        }
        BinValue[GetTotBin(WhichBin)] += Val;
    }

    void SetBinContent(const unsigned* WhichBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        SetBinContent(GetTotBin(WhichBin),Val);
    }
    void SetBinError(const unsigned* WhichBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        SetBinError(GetTotBin(WhichBin),Val);
    }
    void ComputeError(){
        if(!Initialized) {InitWarning(); return;}
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinError[uBin] = sqrt(BinValue[uBin]);
//printf("BinValue[%u]=%f+/-%f\n",uBin,BinValue[uBin],BinError[uBin]);
        }
    }

    void SetBinCenter(const unsigned short& sDim, const unsigned& WhichBin, const double& Val){
        if(!Initialized) {InitWarning(); return;}
        if(sDim>=Dim)return;
        if(WhichBin>=NumBins[sDim])return;
        if(!BinCenter)return;
        if(!BinCenter[sDim])return;
        BinCenter[sDim][WhichBin] = Val;
    }

    double GetBinCenter(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim)return 0;
        if(WhichBin>=NumBins[sDim])return 0;
        if(!BinCenter)return 0;
        if(!BinCenter[sDim])return 0;
        return BinCenter[sDim][WhichBin];
    }

    Type GetBinContent(const unsigned& WhichTotBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(WhichTotBin>=TotNumBins+2) return 0;
        return BinValue[WhichTotBin];
    }
    Type GetBinContent(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        return GetBinContent(GetTotBin(WhichBin));
    }

    Type GetBinError(const unsigned& WhichTotBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(WhichTotBin>=TotNumBins) return 0;
        return BinError[WhichTotBin];
    }
    Type GetBinError(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        return GetBinError(GetTotBin(WhichBin));
    }
    //Type GetBinError(const unsigned short& sDim, const unsigned& WhichBin) const{
    //    return GetBinError(sDim,GetBin(sDim,WhichBin));
    //}

    double GetLowEdge(const unsigned short& sDim) const{
        return GetBinLowEdge(sDim,0);
    }
    double GetUpEdge(const unsigned short& sDim) const{
        return GetBinUpEdge(sDim,NumBins[sDim]-1);
    }
    double GetBinLowEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        //i.e. GetBinLowEdge(NumBins) will return the most upper limit
        if(WhichBin>=NumBins[sDim]+1) return 0;
        return BinRange[sDim][WhichBin];
    }
    double GetBinUpEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
        return BinRange[sDim][WhichBin+1];
    }
    double GetBinSize(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
//printf("sDim=%u; WhichBin=%u; BR1=%f; BR0=%f\n",sDim,WhichBin,BinRange[sDim][WhichBin+1],BinRange[sDim][WhichBin]);
        return (BinRange[sDim][WhichBin+1]-BinRange[sDim][WhichBin]);
    }
    double GetBinSize(const unsigned& WhichTotBin) const{
        double BinSize = 1;
        unsigned* WhichBin = new unsigned [Dim];
        GetBinCoordinates(WhichTotBin,WhichBin);
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            BinSize *= GetBinSize(sDim,WhichBin[sDim]);
        }
        delete [] WhichBin;
        return BinSize;
    }

    unsigned GetBin(const unsigned short& sDim, const double& xVal) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        if(NumBins[sDim]<=1) return 0;
        unsigned WhichBin=(NumBins[sDim]+1)/2;
        unsigned BinMod=4;
        unsigned BinStep;
        //makes sure that the xVal is in BinRange. If not, the returned value is either
        //NumBins or NumBins+1, depending on if we have an underflow or overflow
        if(xVal<BinRange[sDim][0]) return NumBins[sDim];
        if(xVal>BinRange[sDim][NumBins[sDim]]) return NumBins[sDim]+1;
        while(true){
            if(BinRange[sDim][WhichBin]<=xVal && BinRange[sDim][WhichBin+1]>=xVal){
                return WhichBin;
            }
    //!check the logic of this <=, I made it so, since else we crashed on limit values, before that it was <
            else if(xVal<=BinRange[sDim][WhichBin]){
                BinStep = (NumBins[sDim]+1)/BinMod;
                WhichBin -= BinStep?BinStep:1;
                BinMod *= 2;
            }
            else{
                BinStep = (NumBins[sDim]+1)/BinMod;
                WhichBin += BinStep?BinStep:1;
                BinMod *= 2;
            }
        }
    }
    double* GetBinRange(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        double* BINRANGE = new double [NumBins[sDim]+1];
        for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
            BINRANGE[uBin] = BinRange[sDim][uBin];
        }
        return BINRANGE;
    }
    double* GetBinCenters(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        double* BINCENTER = new double [NumBins[sDim]+1];
        for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
            BINCENTER[uBin] = BinCenter[sDim][uBin];
        }
        return BINCENTER;
    }

    Type Eval(const double* xVal, const bool& EvalTheError=false) const{
        if(!Initialized) {InitWarning(); return 0;}
        //this is here to make it thread-safe, but maybe hinders performance???
        double* xValue1 = new double [Dim];
        double* xValue2 = new double [Dim];
        double* DeltaX1 = new double [Dim];
        double* DeltaX2 = new double [Dim];
        unsigned* xBin1 = new unsigned [Dim];
        unsigned* xBin2 = new unsigned [Dim];
        //Type* FunctionValue1 = new Type [Dim];
        //Type* FunctionValue2 = new Type [Dim];
        //Type Normalization = new Type* [Dim];
        //Type TotalNorm=0;

        for(unsigned short sDim=0; sDim<Dim; sDim++){
            //xBin[sDim] = new Type [2];
            //xValue[sDim] = unsigned Type [2];
            //FunctionValue[sDim] = new Type [2];
            //Normalization[sDim] = new Type [2];

            xBin1[sDim] = GetBin(sDim,xVal[sDim]);
            xValue1[sDim] = GetBinCenter(sDim,xBin1[sDim]);
            if(NumBins[sDim]==1){
                //set as dummy only one of the x-values with weight zero,
                //this way the value of x2 will actually play no role
                //despite of that just in case we set x1=x2
                xBin1[sDim]=0;
                xBin2[sDim]=0;
                xValue1[sDim]=GetBinCenter(sDim,0);
                xValue2[sDim]=GetBinCenter(sDim,0);
                DeltaX1[sDim]=0;
                DeltaX2[sDim]=1;
            }
            else if(xVal[sDim]==xValue1[sDim]){
                xBin2[sDim]=xBin1[sDim];
                xValue2[sDim]=xValue1[sDim];
                DeltaX1[sDim]=0;
                DeltaX2[sDim]=1;
            }
            //we take this and the previous bin
            //this is the case when:
            //a) xVal is smaller than the bin center and there is a lower bin available
            //b) in case we are in the last bin
            else if( (xVal[sDim]<xValue1[sDim]&&xBin1[sDim]>0&&xBin1[sDim]<NumBins[sDim]) || xBin1[sDim]==NumBins[sDim]-1 ){
                xBin2[sDim]=xBin1[sDim];
                xBin1[sDim]=xBin1[sDim]-1;
                xValue1[sDim] = GetBinCenter(sDim,xBin1[sDim]);
                xValue2[sDim] = GetBinCenter(sDim,xBin2[sDim]);
                DeltaX1[sDim]=xVal[sDim]-xValue1[sDim];
                DeltaX2[sDim]=xValue2[sDim]-xVal[sDim];
            }
            //we take this and the next bin
            //this is the case when:
            //a) xVal is larger than the bin center and there is an upper bin available
            //b) in case we are in the first bin or below
            else if( (xVal[sDim]>xValue1[sDim]&&xBin1[sDim]<NumBins[sDim]-1) || xBin1[sDim]==0 ){
                xBin2[sDim]=xBin1[sDim]+1;
                xValue1[sDim] = GetBinCenter(sDim,xBin1[sDim]);
                xValue2[sDim] = GetBinCenter(sDim,xBin2[sDim]);
                DeltaX1[sDim]=xVal[sDim]-xValue1[sDim];
                DeltaX2[sDim]=xValue2[sDim]-xVal[sDim];
            }
            //if we are in the overflow bin
            else if(xBin1[sDim]==NumBins[sDim]+1){
                xBin1[sDim]=NumBins[sDim]-2;
                xBin2[sDim]=NumBins[sDim]-1;
                xValue1[sDim] = GetBinCenter(sDim,xBin1[sDim]);
                xValue2[sDim] = GetBinCenter(sDim,xBin2[sDim]);
                DeltaX1[sDim]=xVal[sDim]-xValue1[sDim];
                DeltaX2[sDim]=xValue2[sDim]-xVal[sDim];
            }
            //if we are in the underflow bin
            else if(xBin1[sDim]==NumBins[sDim]){
                xBin1[sDim]=0;
                xBin2[sDim]=1;
                xValue1[sDim] = GetBinCenter(sDim,xBin1[sDim]);
                xValue2[sDim] = GetBinCenter(sDim,xBin2[sDim]);
                DeltaX1[sDim]=xVal[sDim]-xValue1[sDim];
                DeltaX2[sDim]=xValue2[sDim]-xVal[sDim];
            }
            else{
                printf("This should not happen, unless there is a bug in DLM_Histo::Eval()\n");
            }
//if(TEMP){
//printf("sDim=%u\n",sDim);
//printf(" xVal=%f\n",xVal[sDim]);
//printf(" xBin1=%u\n",xBin1[sDim]);
//printf(" xBin2=%u\n",xBin2[sDim]);
//printf(" xValue1=%f\n",xValue1[sDim]);
//printf(" xValue2=%f\n",xValue2[sDim]);
//printf(" DeltaX1=%f\n",DeltaX1[sDim]);
//printf(" DeltaX2=%f\n",DeltaX2[sDim]);
//}
        }
        Type Result=0;
        unsigned* BinArray = new unsigned [Dim];
        Type Weight=1;
        Type Norm=0;
        for(unsigned uPer=0; uPer<NumPermutations; uPer++){
            Weight=1;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                if(PER[uPer][sDim]==0){
                    BinArray[sDim]=xBin1[sDim];
                    Weight*=DeltaX2[sDim];
                }
                else{
                    BinArray[sDim]=xBin2[sDim];
                    Weight*=DeltaX1[sDim];
                }
            }
            if(EvalTheError){Result += GetBinError(BinArray)*Weight;}
            else{Result += GetBinContent(BinArray)*Weight;}
            Norm += Weight;
        }
        Result /= Norm;
//if(TEMP){
//unsigned BINARRAY[Dim];
//BINARRAY[0] = GetBin(0,xVal[0]);
//BINARRAY[1] = GetBin(1,xVal[1]);
//printf(" --> RESULT=%f (%f); NORM=%f\n",abs(Result),abs(GetBinContent(BINARRAY)),abs(Norm));
//}
        delete [] xValue1;
        delete [] xValue2;
        delete [] xBin1;
        delete [] xBin2;
        delete [] DeltaX1;
        delete [] DeltaX2;
        delete [] BinArray;

        return Result;
    }

    Type EvalError(const double* xVal) const{
        return Eval(xVal,true);
    }

    bool operator=(const DLM_Histo& other){
//printf("operator=\n");
        //if(!Initialized) {InitWarning(); return false;}
//printf("BinValue=%p\n",BinValue);
        SetUp(other.Dim);
//printf("other.Dim=%u; Dim=%u\n",other.Dim,Dim);
//printf("Dim done\n");
        for(unsigned short sDim=0; sDim<Dim; sDim++) SetUp(sDim,other.NumBins[sDim],other.BinRange[sDim]);
//printf("SetUp done\n");
        if(!Initialize(false)) return false;
//printf("Initialize done\n");
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
                BinRange[sDim][uBin] = other.BinRange[sDim][uBin];
                BinCenter[sDim][uBin] = other.BinCenter[sDim][uBin];
            }
            BinRange[sDim][NumBins[sDim]] = other.BinRange[sDim][NumBins[sDim]];
        }
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin] = other.BinValue[uBin];
            BinError[uBin] = other.BinError[uBin];
        }
        BinValue[TotNumBins] = other.BinValue[TotNumBins];
        BinValue[TotNumBins+1] = other.BinValue[TotNumBins+1];
        return true;
    }
//! FOR DIVISION/MULT:
//http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
//https://apps.dtic.mil/dtic/tr/fulltext/u2/785623.pdf
    bool operator+=(const DLM_Histo& other){
        return AddHisto(other,true,1);
    }
    bool operator-=(const DLM_Histo& other){
        return AddHisto(other,true,-1);
    }
    bool operator*=(const DLM_Histo& other){
        return MultiplyHisto(other,true);
    }
    bool operator/=(const DLM_Histo& other){
        return DivideHisto(other,true);
    }
    //note that with this type of construct, DLM_Histo<Type> histo=histo1+histo2 invokes the constructor only ones!
    DLM_Histo<Type> operator+(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result+=other;
//printf("Result=%p\n",&Result);
        return Result;
    }
    DLM_Histo<Type> operator-(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result-=other;
        return Result;
    }
    DLM_Histo<Type> operator*(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result*=other;
        return Result;
    }
    DLM_Histo<Type> operator/(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result/=other;
        return Result;
    }

    bool operator+=(const Type& Value){
        AddToAll(Value);
        return true;
    }
    bool operator-=(const Type& Value){
        AddToAll(-Value);
        return true;
    }
    bool operator*=(const Type& Value){
        Scale(Value);
        return true;
    }
    bool operator/=(const Type& Value){
        Scale(1./Value);
        return true;
    }

    DLM_Histo<Type> operator+(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result+=Value;
        return Result;
    }
    DLM_Histo<Type> operator-(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result-=Value;
        return Result;
    }
    DLM_Histo<Type> operator*(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result*=Value;
        return Result;
    }
    DLM_Histo<Type> operator/(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result/=Value;
        return Result;
    }

protected:
    void CleanUp(){
//printf("CU: BinRange\n");
        if(BinRange){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                delete [] BinRange[sDim]; BinRange[sDim]=NULL;
            }
            delete [] BinRange; BinRange=NULL;
        }
//printf("CU: BinValue\n");
        if(BinValue){
            delete [] BinValue; BinValue=NULL;
        }
//printf("CU: BinError\n");
        if(BinError){
            delete [] BinError; BinError=NULL;
        }
//printf("CU: BinCenter\n");
        if(BinCenter){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                delete [] BinCenter[sDim]; BinCenter[sDim]=NULL;
            }
            delete [] BinCenter; BinCenter=NULL;
        }
//printf("CU: NumBins\n");
        delete [] NumBins; NumBins=NULL;
//printf("CU: PER\n");
        if(PER){
            for(unsigned uPer=0; uPer<NumPermutations; uPer++){
                delete [] PER[uPer]; PER[uPer]=NULL;
            }
            delete [] PER; PER=NULL;
        }
        Initialized=false;
    }
    void CleanUp(const unsigned short sDim){
        if(BinRange){
            if(BinRange[sDim]){
                delete [] BinRange[sDim];
                BinRange[sDim]=NULL;
            }
        }
        if(BinCenter){
            if(BinCenter[sDim]){
                delete [] BinCenter[sDim];
                BinCenter[sDim]=NULL;
            }
        }
    }
    bool SameStructure(const DLM_Histo& other) const{
        if(Dim!=other.Dim) return false;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            if(NumBins[sDim]!=other.NumBins[sDim]) return false;
            for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
                if(BinRange[sDim][uBin]!=other.BinRange[sDim][uBin]) return false;
                if(BinCenter[sDim][uBin]!=other.BinCenter[sDim][uBin]) return false;
            }
            if(BinRange[sDim][NumBins[sDim]]!=other.BinRange[sDim][NumBins[sDim]]) return false;
        }
        if(TotNumBins!=other.TotNumBins) return false;
        return true;
    }

    void InitPER(){
        NumPermutations = 1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            NumPermutations *= 2;
        }

        PER = new char* [NumPermutations];
        //char UpdatedPer = 0;
        for(unsigned uPer=0; uPer<NumPermutations; uPer++){
            PER[uPer] = new char [Dim];
            bool Increase = true;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                if(!uPer) PER[uPer][sDim]=0;
                else{
                    PER[uPer][sDim]=PER[uPer-1][sDim];
                    if(Increase && PER[uPer][sDim]==0){
                        PER[uPer][sDim]=1;
                        Increase=false;
                    }
                    else if(Increase){
                        PER[uPer][sDim]=0;
                        Increase=true;
                    }
                }
            }
        }
    }

    void InitWarning() const{
        printf("\033[1;33mWARNING:\033[0m DLM_Histo cannot be used until fully SetUp and Initialized!\n");
    }

    unsigned short Dim;
    unsigned TotNumBins;
    unsigned* NumBins;
    double** BinRange;
    //the last two bins are under/overflow
    Type* BinValue;
    Type* BinError;
    double** BinCenter;

    unsigned NumPermutations;
    char** PER;

    bool Initialized;
};
#endif
