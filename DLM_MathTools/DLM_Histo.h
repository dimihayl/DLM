
#ifndef DLM_HISTO_H
#define DLM_HISTO_H

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <math.h>
#include "DLM_CppTools.h"
#include "DLM_Random.h"
#include <fstream>
#include <iostream>
#include <stdint.h>

//to use usleep in tests
#include <unistd.h>
#include <limits>
#include <type_traits>

using namespace std;

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
            BinCenter[0] = 1.331e121;
        }
        else{
            Type BinWidth = (xmax-xmin)/Type(NumBins);
//printf("NumBins=%u; BinWidth=%f;\n",NumBins,BinWidth);
            for(unsigned uBin=0; uBin<=NumBins; uBin++){
                BinRange[uBin] = xmin + Type(uBin)*BinWidth;
                if(uBin!=NumBins) BinCenter[uBin] = 1.331e121;
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
        if(BinCenter[WhichBin]!=1.331e121) return BinCenter[WhichBin];
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

//UPDATE-OKAY
    DLM_Histo():DefaultMaxNumBins(4200000000){
        ConstructorState();
    }
//UPDATE-OKAY
    DLM_Histo(const DLM_Histo& other):DLM_Histo(){
        ConstructorState();
        Copy(other,true);
    }
    ~DLM_Histo(){
        CleanUp();
    }

//UPDATE-OKAY
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
        CumulativeValue = NULL;
        rangen = NULL;
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
        Initialized = false;
        CumUpdated = false;
    }
//UPDATE-OKAY
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
        Initialized = false;
        CumUpdated = false;
    }
//UPDATE-OKAY
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
        Initialized = false;
        CumUpdated = false;
    }
//UPDATE-OKAY
    //sets up the sDim in the same way as the sDim_other of other
    void SetUp(const unsigned short& sDim, const DLM_Histo<Type>& other, const unsigned short& sDim_other){
        SetUp(sDim, other.GetNbins(sDim_other), other.BinRange[sDim_other]);
    }

    //sparse_limit == 0 means we have a normal histogram, with ALL bins initialized
    //sparse_limit == 1 means we have a sparse histogram (saves memory, but more CPU usage)
    //sparse_limit > 1 tells the code the the exact number of initialized bins at the start. 
    //The input is in MILLIONS (10^6)
    //Default value: 23, Minimum value: 2, Maximum value: 4200, but for >~ 1000 there is almost 100% certainty of running out of memory given hardware in 2026.
    //If the total number of bins is smaller than sparse_limit, normal matrix will be used. Use sparse_limit > 1 only if you know what you are doing. It should be a number above 2M
    //
    //elements_per_bucket controls roughly speaking how much memory we save at what extra CPU cost
    //The input is in THOUSANDS (10^3)
    //the lower the number the more memory we save, but CPU overhead gets larger
    //the higher the number the less memory we save, but CPU overhead gets lower.
    //recommened: do not change the default, if you have to, vary it:
    //Default value: 230, Minimum value: 10, Maximum value: 1000.
//UPDATE-OKAY
    void SetUpSparse(unsigned sparse_limit = 1, unsigned elements_per_bucket = 230){

//printf("hello from sparse\n");
        if(sparse_limit>4200 || elements_per_bucket<10 || elements_per_bucket>1000){
            printf("\033[1;33mWARNING:\033[0m DLM_Histo SetUpSparse(...) function failed, check your arguments!\n");
            return;            
        }
        if(sparse_limit>1 && sparse_limit*1000 < 2*elements_per_bucket){
            printf("\033[1;33mWARNING:\033[0m DLM_Histo SetUpSparse(...) function failed, check your arguments!\n");
            return;
        }

        if(sparse_limit==0){
            MaxNumBins = DefaultMaxNumBins;
            Initialized = false;
            CumUpdated = false;
            return;             
        }

        unsigned nsb;
        unsigned mnb;
        if(sparse_limit==1){
            mnb = 23000000;
        }
        else{
            mnb = sparse_limit * 1000000;
        }
        elements_per_bucket *= 1000;

        //if we cannot get the exact number of elements, we reduce the buckets by 1 (safety first, this uses less memory)
        nsb = mnb/(elements_per_bucket) - bool(mnb%(elements_per_bucket));
        //re-evalute the actual MaxNumBins
        mnb = (elements_per_bucket) * nsb;

        //no change
        //if(mnb == MaxNumBins && NumDesiredSparseBuckets == nsb){
        if(mnb == MaxNumBins && NumElPerSparseBucket == elements_per_bucket){
            return;
        }

        MaxNumBins = mnb;
        //NumDesiredSparseBuckets = nsb;
        NumElPerSparseBucket = (elements_per_bucket);

//printf("MaxNumBins = %u\n",MaxNumBins);
//printf("NumElPerSparseBucket = %u\n",NumElPerSparseBucket);
//usleep(1000e3);

        Initialized = false;
        CumUpdated = false;     
//printf("goodbye from sparse (nsb %u)\n",NumDesiredSparseBuckets);   
    }



//UPDATE-OKAY
    bool Initialize(const bool& ZeroElements=true){
//printf("Initialize 1: %i\n",Initialized);
        if(!Dim||!NumBins||!BinRange||!BinCenter) return false;
        if(!Initialized){
            TotNumBins=1;
//printf("Initialize 2\n");
//printf(" BinValue = %p\n",BinValue);
//printf(" BinError = %p\n",BinError);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                if(!NumBins[sDim]||!BinRange[sDim]||!BinCenter[sDim]) {TotNumBins=0;return false;}
                TotNumBins*=NumBins[sDim];
//printf("sDim=%u; NumBins[sDim]=%u; BinRange[sDim]=%p; BinCenter[sDim]=%p | TNB = %lu\n",sDim,NumBins[sDim],BinRange[sDim],BinCenter[sDim],TotNumBins);
            }
            if(BinValue){delete [] BinValue; BinValue=NULL;}
            BucketId.clear();
            BucketOrder.clear();
            if(CumulativeValue){delete [] CumulativeValue; CumulativeValue=NULL;}
            if(BinError){delete [] BinError; BinError=NULL;}

            if(TotNumBins<=MaxNumBins){
                MaxNumBins = TotNumBins;
            }


//printf("MaxNumBins = %u\n",MaxNumBins); 
//usleep(1000e3);
            BinValue = new Type[MaxNumBins+2];
            BinError = new Type[MaxNumBins];

        }
        if(ZeroElements){
//printf("Initialize 2\n");
            for(unsigned ej=0; ej<MaxNumBins; ej++){
                BinValue[ej]=NULL_ELEM;
                BinError[ej]=NULL_ELEM;
            }
            BinValue[MaxNumBins]=NULL_ELEM;
            BinValue[MaxNumBins+1]=NULL_ELEM;            
        }
        Initialized = true;
        CumUpdated = false;
//printf("Initialize 3\n");
        return true;
    }
    bool IsInitialized(){
        return Initialized;
    }
    uint64_t GetTotBin(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        uint64_t TotBin=0;
        uint64_t Base=1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            if(WhichBin[sDim]>=NumBins[sDim]) return TotNumBins;
            TotBin+=Base*WhichBin[sDim];
            Base*=NumBins[sDim];
        }
        return TotBin;
    }
//UPDATE-OKAY
    uint64_t GetTotBin(const unsigned& WhichX, const unsigned& WhichY) const{
        unsigned which_bin[2] = {WhichX, WhichY};
        return GetTotBin(which_bin);
    }
//UPDATE-OKAY
    uint64_t GetTotBin(const unsigned& WhichX, const unsigned& WhichY, const unsigned& WhichZ) const{
        unsigned which_bin[3] = {WhichX, WhichY, WhichZ};
        return GetTotBin(which_bin);
    }   
//UPDATE-OKAY
    uint64_t FindTotBin(const double* AxisValue) const{
        if(!Initialized) {InitWarning(); return 0;}
        unsigned* WhichBin = new unsigned [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            WhichBin[sDim] = GetBin(sDim,AxisValue[sDim]);
        }
        return GetTotBin(WhichBin);
    }
//UPDATE-OKAY
    Type& GetBinElement(const uint64_t& WhichTotBin){
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(WhichTotBin>=TotNumBins) {return NULL_ELEM;}
        unsigned ej = smid(WhichTotBin);
        return ej==-1?NULL_ELEM:BinValue[ej];
    }
//UPDATE-OKAY
    //this gives the GLOBAL axis coordinates (idependent of using a sparse matrix)
    //BinIdPerAxis[sDim] -> saves for each dim the bin number
    void GetBinCoordinates(const uint64_t& WhichTotBin, unsigned* BinIdPerAxis) const{
        if(!Initialized) {InitWarning(); return;}
        //underflow
        if(WhichTotBin==TotNumBins){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                BinIdPerAxis[sDim] = NumBins[sDim];
            }
            return;
        }
        //overflow
        else if(WhichTotBin==TotNumBins+1){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                BinIdPerAxis[sDim] = NumBins[sDim]+1;
            }
            return;
        }
        //n/a
        else if(WhichTotBin>=TotNumBins+2){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                BinIdPerAxis[sDim] = -1;
            }
            return;
        }
        uint64_t rem = WhichTotBin;
        unsigned short isDim;
        uint64_t Base=TotNumBins;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            isDim = Dim-sDim-1;
            Base/=NumBins[isDim];
            if(Base > 0){
                BinIdPerAxis[isDim] = (unsigned)(rem / Base);
                rem %= Base;
            } else { BinIdPerAxis[isDim] = 0; }
        }
    }
//UPDATE-OKAY
    void GetBinAxisCenters(const uint64_t& WhichTotBin, double* BinAxisCenter){
        if(!Initialized) {InitWarning(); return;}
        unsigned* BinId = new unsigned [Dim];
        GetBinCoordinates(WhichTotBin, BinId);
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            //underflow
            if(BinId[sDim]==NumBins[sDim]){
                BinAxisCenter[sDim] = BinRange[sDim][0];
            }
            //overflow
            else if(BinId[sDim]==NumBins[sDim]+1){
                BinAxisCenter[sDim] = BinRange[sDim][NumBins[sDim]];
            }
            //n/a
            else if(BinId[sDim]==-1){
                BinAxisCenter[sDim] = std::numeric_limits<double>::quiet_NaN();
            }
            //normal situation
            else{
                BinAxisCenter[sDim] = BinCenter[sDim][BinId[sDim]];
            }
        }
        delete [] BinId;
    }
//UPDATE-OKAY
    //the full flag means if we should copy the full structure including MaxNumBins etc.
    bool Copy(const DLM_Histo& other, const bool full=false){
        //if(!Initialized) {InitWarning(); return false;}
        //return operator=(other);

        bool sstr = SameStructure(other);
        bool alles_gut = sstr && MaxNumBins==other.MaxNumBins;
//printf("sstr%u alles_gut%u\n",sstr,alles_gut);

        //we want to initialize the histogram if we have different structure 
        //or in case we have a full copy option, but differences in either structure or bin setup
        if(!sstr || (!alles_gut&&full) ){
            SetUp(other.Dim);
            for(unsigned short sDim=0; sDim<Dim; sDim++) SetUp(sDim,other.NumBins[sDim],other.BinRange[sDim]);
            MaxNumBins = other.MaxNumBins;
            if(!Initialize(false)) return false;
            //at that point we ensure everything is identical
            alles_gut = true;
            sstr = true;
        }
        //here we copy the bin ranges of the axes
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
                BinRange[sDim][uBin] = other.BinRange[sDim][uBin];
                BinCenter[sDim][uBin] = other.BinCenter[sDim][uBin];
            }
            BinRange[sDim][NumBins[sDim]] = other.BinRange[sDim][NumBins[sDim]];
        }

        //here we make a flat copy, in case we have same structure and binning
        if(alles_gut){
            for(unsigned ej=0; ej<MaxNumBins; ej++){
                BinValue[ej] = other.BinValue[ej];
                BinError[ej] = other.BinError[ej];
            }
            BinValue[MaxNumBins] = other.BinValue[other.MaxNumBins];
            BinValue[MaxNumBins+1] = other.BinValue[other.MaxNumBins+1];
            
            NumElPerSparseBucket = other.NumElPerSparseBucket;
            BucketId = other.BucketId;
            BucketOrder = other.BucketOrder;            
        }
        //if there are some differences, we have to be careful
        else{
//printf("hello there\n");
            //a caveat: the MaxNumBins CAN be sometimes different.
            //In addition, one histo can be sparse, the other not.
            //for that reason, we have to equalize the elements using smid.
            //but first thing firs,t we set all elements to zero
            for(unsigned ej=0; ej<MaxNumBins; ej++){
                BinValue[ej] = NULL_ELEM;
                BinError[ej] = NULL_ELEM;
            }
            
            for(unsigned oej=0; oej<other.MaxNumBins; oej++){
                uint64_t oei = other.ismid(oej);
//printf("oej %u --> %lu !!!\n", oej, oei);
                if(oei==-1) continue;
                if(other.BinValue[oej]==NULL_ELEM && other.BinError[oej]==NULL_ELEM) continue;
                //at this point we have ensured the same structure
                uint64_t ei = oei;
                unsigned ej = smid(ei,true);
//printf("oej %u eoi %lu ei %lu ej %u | %e %e\n", oej,oei,ei,ej,other.BinValue[oej],other.BinError[oej]);
                if(ej==-1) break;
                BinValue[ej] = other.BinValue[oej];
                BinError[ej] = other.BinError[oej];            
            }
    
            BinValue[MaxNumBins] = other.BinValue[other.MaxNumBins];
            BinValue[MaxNumBins+1] = other.BinValue[other.MaxNumBins+1];
        }
        
        CumUpdated = false;
        //TotNumBins = other.TotNumBins;
        //NumDesiredSparseBuckets = other.NumDesiredSparseBuckets;

        return true;
    }
//UPDATE-OKAY
    bool AddHisto(const DLM_Histo& other, const bool witherror=true, const double& scale=1, const bool& weighted_avg=false){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(other.Dim!=Dim) {return false;}
        bool sstr = SameStructure(other);

        //b = bucket
        //e = element
        //i = global
        //j = in array (local)
        uint64_t ei,oei;
        unsigned ej;
        unsigned oej;

        //overflow bins
        if(sstr){
            BinValue[MaxNumBins] += scale*other.BinValue[MaxNumBins];
            BinValue[MaxNumBins+1] += scale*other.BinValue[MaxNumBins+1];
        }

        unsigned* BinId = sstr?NULL:new unsigned [Dim];
        double* xVal = sstr?NULL:new double [Dim];

        Type VAR = NULL_ELEM;
        Type bin_err = NULL_ELEM;
        Type obin_err = NULL_ELEM;

        //iterate over all bins of THIS
        for(ej=0; ej<MaxNumBins; ej++){
            ei = ismid(ej);
            //zero bin, skip for now
            if(ei==-1) continue;
            if(sstr){
                //due to the same structure
                oei = ei;
                oej = smid(oei);
                //zero bin, nothing to do as both bin value and error are kept the same
                if(oej==-1) continue;

                //note, that depending on the construction of the matrix we can still end up having zero values and zero erros
                //thus they should be handled (important for the weighted case).
                //in the weighted case, if both errors are zero -> just take the average value
                //if one error is zero -> take this value as a final result

                if(weighted_avg){
                    if(BinError[ej]==NULL_ELEM && other.BinError[oej]==NULL_ELEM){
                        BinValue[ej] = (BinValue[ej] + other.BinValue[oej])*0.5;
                        BinError[ej] = 0;
                    }
                    else if(BinError[ej]==NULL_ELEM){
                        //do nothing
                    }
                    else if(other.BinError[oej]==NULL_ELEM){
                        BinValue[ej] = other.BinValue[oej];
                        BinError[ej] = other.BinError[oej];;
                    }
                    else{
                        VAR = 1./(1./BinError[ej]/BinError[ej] + 1./other.BinError[oej]/other.BinError[oej]);
                        BinValue[ej] = BinValue[ej]/BinError[ej]/BinError[ej] + other.BinValue[oej]/other.BinError[oej]/other.BinError[oej];
                        BinValue[ej] *= VAR;
                        BinError[ej] = sqrt(VAR);
                    }
                }                
                else{
                    BinValue[ej] += scale*other.BinValue[oej];
                    if(witherror && ej<MaxNumBins){
                        BinError[ej] = sqrt(BinError[ej]*BinError[ej]+scale*scale*other.BinError[oej]*other.BinError[oej]);
                    }
                }

            }
            //no same structure
            else{
                GetBinCoordinates(ei,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                }

                if(weighted_avg){
                    bin_err = BinError[ej];
                    obin_err = other.Eval(xVal,1);
                    if(bin_err==NULL_ELEM && obin_err==NULL_ELEM){
                        BinValue[ej] = (BinValue[ej] + other.Eval(xVal))*0.5;
                        BinError[ej] = 0;
                    }
                    else if(bin_err==NULL_ELEM){
                        //do nothing
                    }                    
                    else if(obin_err==NULL_ELEM){
                        BinValue[ej] = (BinValue[ej] + other.Eval(xVal))*0.5;
                        BinError[ej] = 0;
                    }
                    else{
                        VAR = 1./(1./bin_err/bin_err + 1./obin_err/obin_err);
                        BinValue[ej] = BinValue[ej]/bin_err/bin_err + other.Eval(xVal)/obin_err/obin_err;
                        BinValue[ej] *= VAR;
                        BinError[ej] = sqrt(VAR);
                    }
                }
                else{
                    BinValue[ej] += scale*other.Eval(xVal);
                    if(witherror){
                        BinError[ej] = sqrt(BinError[ej]*BinError[ej]+scale*scale*other.Eval(xVal,true)*other.Eval(xVal,true));
                    }     
                }
           
            }
        }

        //second iteration over the bins of OTHER
        //needed in order to add non-zero OTHER elements to zero THIS elements
        //this issue can occur only for sparse THIS
        if(TotNumBins>MaxNumBins){
            for(oej=0; oej<other.MaxNumBins; oej++){
                oei = ismid(oej);
                //zero bin, all good
                if(oei==-1) continue;
                if(sstr){
                    ei = oei;
                    ej = smid(ei);
                    //we found what we were looking for, zero THIS and non-zero OTHER
                    //we set it equal to the OTHER even when using avg. weight
                    if(ej==-1){
                        ej = smid(ei,true);
                        if(ej==-1) break;
                        BinValue[ej] = scale*other.BinValue[oej];
                        BinError[ej] = scale*other.BinError[oej];
                    }
                }
                else{
                    GetBinCoordinates(oei,BinId);
                    for(unsigned short sDim=0; sDim<Dim; sDim++){
                        xVal[sDim] = other.BinCenter[sDim][BinId[sDim]];
                    }
                    ei = FindTotBin(xVal);
                    ej = smid(ei);
                    //we found what we were looking for, zero THIS and non-zero OTHER
                    //we set it equal to the OTHER even when using avg. weight
                    if(ej==-1){
                        GetBinCoordinates(ei,BinId);
                        for(unsigned short sDim=0; sDim<Dim; sDim++){
                            xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                        }       
                        ej = smid(ei,true);
                        if(ej==-1) break;
                        BinValue[ej] = scale*other.Eval(xVal);
                        BinError[ej] = scale*other.Eval(xVal,true);
                    }                                              
                }
            }
        }

        if(BinId) delete [] BinId;
        if(xVal) delete [] xVal;





        /*
        //we iterate over the possibly non-zero bins of the other histo, and add them to our histo
        for(oej=0; oej<other.MaxNumBins+2; oej++){

            //sparse
            if(other.TotNumBins>other.MaxNumBins){
                oei = other.ismid(oej);
                //zero bin (both value and error are zero)
                if(oei==-1) continue;
            }
            else{
                oei = oej;
            }

            if(sstr){
                //due to the same structure
                ei = oei;
                if(TotNumBins>MaxNumBins){
                    ej = smid(ei, true);
                }
                else{
                    //if we have the same structure, and this is not sparse, regardless of if other is sparse or not ei = oei = ej, so this hold
                    ej = oei;
                }
if(ei==6)
printf("%u %u %e + %e x %e\n",ej,oej,BinValue[ej],scale,other.BinValue[oej]);
                BinValue[ej] += scale*other.BinValue[oej];
if(ei==6)
printf(" -> %e\n",BinValue[ej]);
                if(witherror && ej<MaxNumBins){
                    BinError[ej] = sqrt(BinError[ej]*BinError[ej]+scale*scale*other.BinError[oej]*other.BinError[oej]);
                }
            }
            else{
printf("NOSSTR\n");
                unsigned* BinId = new unsigned [Dim];
                double* xVal = new double [Dim];


                //this will partially work sparse or not, HOWEVER, if we have THIS as sparse
                //it may well be that we have missed some non-zero bins in other ...
                for(ej=0; ej<MaxNumBins; ej++){
                    ei = other.ismid(ej);
                    if(ei==-1) continue;//zero bin
                    GetBinCoordinates(ei,BinId);
                    for(unsigned short sDim=0; sDim<Dim; sDim++){
                        xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                    }
                    BinValue[ej] += scale*other.Eval(xVal);
                    if(witherror){
                        BinError[ej] = sqrt(BinError[ej]*BinError[ej]+scale*scale*other.Eval(xVal,true)*other.Eval(xVal,true));
                    }
                }


                //... for that reason here we have to try and add to THIS the non-zero OTHER bins that we missed
                if(TotNumBins>MaxNumBins){
printf("sparse correction\n");
                    for(oej=0; oej<other.MaxNumBins; oej++){
                        oei = other.ismid(oej);
                        if(oei==-1) continue;//zero bin
                        GetBinCoordinates(oei,BinId);
                        for(unsigned short sDim=0; sDim<Dim; sDim++){
                            xVal[sDim] = other.BinCenter[sDim][BinId[sDim]];
                        }
                        ei = FindTotBin(xVal);
                        ej = smid(ei);
                        //i.e. we have bin center in the other bin that corresponds to a zero in this
                        if(ej==-1){
                            //get the bin centre coordinates of this
                            GetBinCoordinates(ei,BinId);
                            for(unsigned short sDim=0; sDim<Dim; sDim++){
                                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                            }
                            SetBinContent(ei, scale*other.Eval(xVal));
                            SetBinError(ei, scale*other.Eval(xVal,true));
                        }
                    }

                }

                delete [] BinId;
                delete [] xVal; 
             
            }


        }//oej
        */

    
        CumUpdated = false;
        return true;
    }

//UPDATE-LEGACY-DELETE
    bool AddHistoLEGACY(const DLM_Histo& other, const bool witherror=true, const Type& scale=1){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(SameStructure(other)){
          for(unsigned uBin=0; uBin<TotNumBins+2; uBin++){
            BinValue[uBin] += scale*other.BinValue[uBin];
          }
          if(witherror)
              for(unsigned uBin=0; uBin<TotNumBins; uBin++)
                  BinError[uBin] = sqrt(BinError[uBin]*BinError[uBin]+scale*scale*other.BinError[uBin]*other.BinError[uBin]);
        }
        else{
          if(other.Dim!=Dim) {return false;}
          unsigned* BinId = new unsigned [Dim];
          double* xVal = new double [Dim];
          for(unsigned uBin=0; uBin<TotNumBins; uBin++){
              GetBinCoordinates(uBin,BinId);
              for(unsigned short sDim=0; sDim<Dim; sDim++){
                  xVal[sDim] = BinCenter[sDim][BinId[sDim]];
              }
              BinValue[uBin] += scale*other.Eval(xVal);
              if(witherror){
                BinError[uBin] = sqrt(BinError[uBin]*BinError[uBin]+scale*scale*other.Eval(xVal,true)*other.Eval(xVal,true));
              }
          }
          delete [] BinId;
          delete [] xVal;
        }
        CumUpdated = false;
        return true;
    }

//UPDATE-LEGACY-DELETE
    bool AddWeightedHistoLEGACY(const DLM_Histo& other, const Type& scale=1){
        if(!Initialized) {InitWarning(); return false;}
        if(!Initialized) return false;
        unsigned* BinId = new unsigned [Dim];
        double* xVal = new double [Dim];

        for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
            }
            Type VAR = 1./(1./BinError[uBin]/BinError[uBin] + 1./other.Eval(xVal,1)/other.Eval(xVal,1));
            BinValue[uBin] = BinValue[uBin]/BinError[uBin]/BinError[uBin] + other.Eval(xVal)/other.Eval(xVal,1)/other.Eval(xVal,1);
            BinValue[uBin] *= VAR;
            BinError[uBin] = sqrt(VAR);
        }
        CumUpdated = false;
        delete [] BinId;
        delete [] xVal;
        return true;
    }

//UPDATE-LEGACY-KEEP
    //a weighted average (based on errors) of the histograms
    bool AddWeightedHisto(const DLM_Histo& other, const double& scale=1){
        return AddHisto(other,true,scale,true);
        /*
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(other.Dim!=Dim) {return false;}
        bool sstr = SameStructure(other);

        //b = bucket
        //e = element
        //i = global
        //j = in array (local)
        uint64_t ei,oei;
        unsigned ej;
        unsigned oej;
        Type VAR;

        Type other_bin_val;
        Type other_bin_err;        

        //we iterate over the possibly non-zero bins of the other histo, and add them to our histo
        for(oej=0; oej<other.MaxNumBins+2; oej++){
            //sparse
            if(other.TotNumBins>other.MaxNumBins){
                oei = other.ismid(oej);
                //zero bin (both value end error)
                //we continue, as we assume this is a non-set bin that plays no role.
                if(oei==-1) continue;
            }
            else{
                oei = oej;
            }
            //we continue, as we assume this is a non-set bin that plays no role
            //this makes sence for the weighted sum, as zero-error should not really occur anyway
            if(other.BinValue[oej]==NULL_ELEM && other.BinError[oej]==NULL_ELEM){
                continue;
            }

            if(sstr){
                //due to the same structure
                ei = oei;
                if(TotNumBins>MaxNumBins){
                    ej = smid(ei, true);
                }
                else{
                    //if we have the same structure, and this is not sparse, regardless of if other is sparse or not ei = oei = ej, so this hold
                    ej = oei;
                }

                if(BinError[ej]==NULL_ELEM && other.BinError[oej]==NULL_ELEM){
                    BinValue[ej] = 0.5*(BinValue[ej]+other.BinValue[oej]);
                }
                else if(BinError[ej]==NULL_ELEM){
                    continue;
                }
                else if(other.BinError[oej]==NULL_ELEM){
                    BinValue[ej] = other.BinValue[oej];
                    BinError[ej] = other.BinError[oej];
                }
                else{
                    VAR = 1./(1./BinError[ej]/BinError[ej] + 1./BinError[oej]/BinError[oej]);
                    BinValue[ej] = BinValue[ej]/BinError[ej]/BinError[ej] + other.BinValue[oej]/other.BinError[oej]/other.BinError[oej];
                    BinValue[ej] *= VAR;
                    BinError[ej] = sqrt(VAR);
                }

            }
            else{
                unsigned* BinId = new unsigned [Dim];
                double* xVal = new double [Dim];

                //this will partially work sparse or not, HOWEVER, if we have THIS as sparse
                //it may well be that we have missed some non-zero bins in other ...
                for(ej=0; ej<MaxNumBins; ej++){
                    ei = other.ismid(ej);
                    if(ei==-1) continue;//zero bin
                    GetBinCoordinates(ei,BinId);

                    for(unsigned short sDim=0; sDim<Dim; sDim++){
                        xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                    }

                    Type& bin_val = BinValue[ej];
                    Type& bin_err = BinError[ej];
                    
                    other_bin_val = other.Eval(xVal);
                    other_bin_err = other.Eval(xVal,true);

                    if(bin_err==NULL_ELEM && other_bin_err==NULL_ELEM){
                        bin_val = 0.5*(bin_val+other_bin_val);
                    }
                    else if(bin_err==NULL_ELEM){
                        continue;
                    }
                    else if(other_bin_err==NULL_ELEM){
                        bin_val = other_bin_err;
                        bin_err = other_bin_err;
                    }
                    else{
                        VAR = 1./(1./bin_err/bin_err + 1./other_bin_err/other_bin_err);
                        bin_val = bin_val/bin_err/bin_err + other_bin_val/other_bin_err/other_bin_err;
                        bin_val *= VAR;
                        bin_err = sqrt(VAR);
                    }

                }

                //... for that reason here we have to try and add to THIS the non-zero OTHER bins that we missed
                if(TotNumBins>MaxNumBins){
printf("sparse correction\n");
                    for(unsigned uBin=0; uBin<other.MaxNumBins; uBin++){
                        other.GetBinCoordinates(uBin,BinId);
                        for(unsigned short sDim=0; sDim<Dim; sDim++){
                            xVal[sDim] = other.BinCenter[sDim][BinId[sDim]];
                        }
                        ei = FindTotBin(xVal);
                        ej = smid(ei);
                        //i.e. we have bin center in the other bin that corresponds to a zero in this
                        if(ej==-1){
                            //get the bin centre coordinates of this
                            GetBinCoordinates(ei,BinId);
                            for(unsigned short sDim=0; sDim<Dim; sDim++){
                                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                            }
                            SetBinContent(ei, scale*other.Eval(xVal));
                            SetBinError(ei, scale*other.Eval(xVal,true));
                        }
                    }

                }

                delete [] BinId;
                delete [] xVal; 
             
            }

        }//oej

        CumUpdated = false;
        return true;    
        */
    }
//UPDATE-LEGACY-DELETE
    bool MultiplyHistoLEGACY(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        //do the operation bin by bin
        if(SameStructure(other)){

          for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin] *= other.BinValue[uBin];
          }
          if(witherror)
              for(uint64_t uBin=0; uBin<TotNumBins; uBin++)
                  BinError[uBin] = sqrt(  pow(BinError[uBin]*other.BinError[uBin],2.)+
                                          pow(BinError[uBin]*other.BinValue[uBin],2.)+
                                          pow(BinValue[uBin]*other.BinError[uBin],2.));
        }
        //not tested
        else{
          if(other.Dim!=Dim) {return false;}
          unsigned* BinId = new unsigned [Dim];
          double* xVal = new double [Dim];

          for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
              GetBinCoordinates(uBin,BinId);
              for(unsigned short sDim=0; sDim<Dim; sDim++){
                  xVal[sDim] = BinCenter[sDim][BinId[sDim]];
              }
              BinValue[uBin] *= other.Eval(xVal);
              if(witherror){
                BinError[uBin] = sqrt(  pow(BinError[uBin]*other.Eval(xVal,true),2.)+
                                        pow(BinError[uBin]*other.Eval(xVal),2.)+
                                        pow(BinValue[uBin]*other.Eval(xVal,true),2.));
              }
          }
          delete [] BinId;
          delete [] xVal;
        }
        CumUpdated = false;
        return true;
    }
//UPDATE-OKAY
    //an important caveat: a zero element in ANY of the two histograms should result in a zero
    //i.e. the zeros should not be skipped here
    bool MultiplyHisto(const DLM_Histo& other, const bool witherror=true){

        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(other.Dim!=Dim) {return false;}
        bool sstr = SameStructure(other);

        //b = bucket
        //e = element
        //i = global
        //j = in array (local)
        uint64_t ei,oei;
        unsigned ej;
        unsigned oej;

        unsigned* BinId = sstr?NULL:new unsigned [Dim];
        double* xVal = sstr?NULL:new double [Dim];

        //iterate over all bins of THIS
        for(ej=0; ej<MaxNumBins; ej++){
            ei = ismid(ej);
            //zero bin (with zero error), we know it will give zero anyway, no need to do anything
            if(ei==-1){
                continue;
            }
            if(sstr){
                //due to the same structure
                oei = ei;
                oej = smid(oei);
                //zero bin, set the corrsponding bin in THIS to zero
                if(oej==-1){
                    //printf("ZERO %u\n",ej);
                    BinValue[ej] = NULL_ELEM;
                    BinError[ej] = NULL_ELEM;
                    continue;  
                }
                
                BinValue[ej] *= other.BinValue[oej];
                if(witherror && ej<MaxNumBins){
                    BinError[ej] = sqrt(  pow(BinError[ej]*other.BinError[oej],2.)+
                                            pow(BinError[ej]*other.BinValue[oej],2.)+
                                            pow(BinValue[ej]*other.BinError[oej],2.));
                }
            }
            //no same structure
            else{
                GetBinCoordinates(ei,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                }

                BinValue[ej] *= other.Eval(xVal);
                if(witherror){
                    BinError[ej] = sqrt(  pow(BinError[ej]*other.Eval(xVal,true),2.)+
                                            pow(BinError[ej]*other.Eval(xVal),2.)+
                                            pow(BinValue[ej]*other.Eval(xVal,true),2.));
                }     
           
            }
        }

        if(BinId) delete [] BinId;
        if(xVal) delete [] xVal;

        CumUpdated = false;
        return true;

    }
//UPDATE-LEGACY-DELETE
    bool DivideHistoLEGACY(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        unsigned* BinId = new unsigned [Dim];
        double* xVal = new double [Dim];
        if(witherror)

            for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
                GetBinCoordinates(uBin,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                }
                if(!other.Eval(xVal)) BinError[uBin] = 1e121;
                else BinError[uBin] = pow(other.Eval(xVal),-2.)*
                                            sqrt(   pow(BinError[uBin]*other.Eval(xVal,true),2.)+
                                                    pow(BinError[uBin]*other.Eval(xVal),2.)+
                                                    pow(BinValue[uBin]*other.Eval(xVal,true),2.));
            }

        for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][BinId[sDim]];
            }
            if(!other.Eval(xVal)) BinValue[uBin] = 1e121;
            else BinValue[uBin] /= other.Eval(xVal);
        }
        delete [] BinId;
        delete [] xVal;
        CumUpdated = false;
        return true;
    }
//UPDATE-OKAY
    bool DivideHisto(const DLM_Histo& other, const bool witherror=true){
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(other.Dim!=Dim) {return false;}
        bool sstr = SameStructure(other);

        //b = bucket
        //e = element
        //i = global
        //j = in array (local)
        uint64_t ei,oei;
        unsigned ej;
        unsigned oej;

        unsigned* BinId = sstr?NULL:new unsigned [Dim];
        double* xVal = sstr?NULL:new double [Dim];

        //iterate over all bins of THIS
        for(ej=0; ej<MaxNumBins; ej++){
            ei = ismid(ej);
            //zero bin (with zero error), we know it will give zero anyway, no need to do anything
            //technically speaking, if OTHER is also zero perhaps we have a NAN, but at the moment this is not considered
            if(ei==-1){
                continue;
            }
            if(sstr){
                //due to the same structure
                oei = ei;
                oej = smid(oei);
                //zero bin, set the corrsponding bin in THIS to NAN
                if(oej==-1){
                    BinValue[ej] = NULL_ELEM;
                    BinError[ej] = NULL_ELEM;
                    continue;  
                }
                
                if(other.BinValue[oej]==NULL_ELEM){
                    if(BinValue[ej]==NULL_ELEM){
                        BinValue[ej] = NAN_ELEM;
                        BinError[ej] = NAN_ELEM;
                    }
                    else{
                        BinValue[ej] = INF_ELEM;
                        BinError[ej] = INF_ELEM;
                    }
                }
                else{
                    if(witherror && ej<MaxNumBins){
                        BinError[ej] = pow(other.BinValue[oej],-2.)*
                                                sqrt(   pow(BinError[ej]*other.BinError[oej],2.)+
                                                        pow(BinError[ej]*other.BinValue[oej],2.)+
                                                        pow(BinValue[ej]*other.BinError[oej],2.));
                    }
                    BinValue[ej] /= other.BinValue[oej];
                }
            }

            //no same structure
            else{
                GetBinCoordinates(ei,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    xVal[sDim] = BinCenter[sDim][BinId[sDim]];
                }

                if(other.BinValue[oej]==NULL_ELEM){
                    if(BinValue[ej]==NULL_ELEM){
                        BinValue[ej] = NAN_ELEM;
                        BinError[ej] = NAN_ELEM;
                    }
                    else{
                        BinValue[ej] = INF_ELEM;
                        BinError[ej] = INF_ELEM;
                    }
                }
                else{
                    BinValue[ej] /= other.Eval(xVal);
                    if(witherror && ej<MaxNumBins){
                        BinError[ej] = pow(other.Eval(xVal),-2.)*
                                                sqrt(   pow(BinError[ej]*other.Eval(xVal,true),2.)+
                                                        pow(BinError[ej]*other.Eval(xVal),2.)+
                                                        pow(BinValue[ej]*other.Eval(xVal,true),2.));
                    }
                }    
        
            }
        }

        if(BinId) delete [] BinId;
        if(xVal) delete [] xVal;

        CumUpdated = false;
        return true;

    }

//UPDATE-LEGACY-KEEP    
    //identical to the one before, but performs the operation faster
    //only valid if the two histograms are of the same structure (else use the function above)
    //LEGACY function, its all now integrated into DivideHisto, which automatically decides what to do
    bool DivideHistoBinByBin(const DLM_Histo& other, const bool witherror=true){
        return DivideHisto(other,witherror);
        //
        if(!Initialized) {InitWarning(); return false;}
        if(!other.Initialized) {InitWarning(); return false;}
        if(!SameStructure(other)) {return false;}
        unsigned* BinId = new unsigned [Dim];
        if(witherror)
            for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
                GetBinCoordinates(uBin,BinId);
                if(!other.BinValue[uBin]) BinError[uBin] = 1e121;
                else BinError[uBin] = pow(other.BinValue[uBin],-2.)*
                                            sqrt(   pow(BinError[uBin]*other.BinError[uBin],2.)+
                                                    pow(BinError[uBin]*other.BinValue[uBin],2.)+
                                                    pow(BinValue[uBin]*other.BinError[uBin],2.));
            }
        for(uint64_t uBin=0; uBin<TotNumBins; uBin++){
            GetBinCoordinates(uBin,BinId);
            if(!other.BinValue[uBin]) BinValue[uBin] = 1e121;
            else BinValue[uBin] /= other.BinValue[uBin];
        }
        delete [] BinId;
        CumUpdated = false;
        return true;
    }
//UPDATE-OKAY
    void AddToAll(const Type& value, const Type& error=0){
        if(!Initialized) {InitWarning(); return;}
        if(TotNumBins>DefaultMaxNumBins){
            printf("\033[1;31mERROR:\033[0m DLM_Histo::AddToAll cannot be executed, as the maximum number of total bins (%.3f+E9) is surpassed!\n",float(DefaultMaxNumBins)/1000000000.);
            return;
        }
        unsigned ej;
        for(uint64_t ei=0; ei<TotNumBins; ei++){
            ej = TotNumBins>MaxNumBins?smid(ei,true):ei;
            static bool err_sparse_mem = false;
            if(ej>=MaxNumBins){
                if(!err_sparse_mem) printf("\033[1;31mERROR:\033[0m The sparse DLM_Histo ran out of memory. No bins will be set from now on. Try to change the SetUpSparse() settings.\n"); 
                err_sparse_mem = true;
                continue;
            }
            //we continue, as while we will not set any new bins, there is a chance we modified perviously initialized bins
            if(ej==-1) continue;
            BinValue[ej]+=value;
            BinError[ej]=Type(sqrt(double(BinError[ej]*BinError[ej])+double(error*error)));
        }
        CumUpdated = false;
    }
//UPDATE-OKAY
    void Scale(const Type& scale){
        if(!Initialized) {InitWarning(); return;}
        if(!scale) return;
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            BinValue[ej] *= scale;
            BinError[ej] *= scale;
        }
        BinValue[MaxNumBins] *= scale;
        BinValue[MaxNumBins+1] *= scale;
        CumUpdated = false;
    }
//UPDATE-OKAY
    void ScaleToIntegral(const bool& UnderOverFlow=false){
      if(!Initialized) {InitWarning(); return;}
      Type TotInt = 0;
      for(unsigned ej=0; ej<MaxNumBins+2*UnderOverFlow; ej++){
        TotInt += BinValue[ej];
      }
      for(unsigned ej=0; ej<MaxNumBins+2*UnderOverFlow; ej++){
        BinValue[ej] /= TotInt;
        BinError[ej] /= TotInt;
      }
      CumUpdated = false;
    }
//UPDATE-OKAY
    void ScaleBin(const uint64_t& WhichTotBin, const Type& Value){
        if(!Initialized) {InitWarning(); return;}
        unsigned ej = smid(WhichTotBin);
        if(ej==-1) return;
        BinValue[ej] *= Value;
        BinError[ej] *= Value;
        CumUpdated = false;
    }
//UPDATE-OKAY
    void ScaleToBinSize(){
        if(!Initialized) {InitWarning(); return;}
        Type BinSize;
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            uint64_t ei = ismid(ej);
            if(ei==-1) continue;
            BinSize = GetBinSize(ei);
            BinValue[ej] /= BinSize;
            BinError[ej] /= BinSize;
        }
        CumUpdated = false;
    }
//UPDATE-OKAY
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
            for(unsigned ej=0; ej<MaxNumBins; ej++){
                BinValue[ej] /= scale;
                BinError[ej] /= scale;
            }
            CumUpdated = false;
        }
    }
//UPDATE-OKAY
    void SetSeed(const unsigned& seed){
      if(rangen) delete rangen;
      rangen = new DLM_Random(seed);
    }
    
//UPDATE-NEED-TO-UPDATE Sample(...)
    //samples randomly assuming that the histogram contains yield (non-normalized to bin width)
    bool SampleYield(double* axisValues, const bool& UnderOverFlow=false, DLM_Random* RanGen=NULL){
      return Sample(axisValues,UnderOverFlow,true,RanGen);
    }
    double SampleYield(const bool& UnderOverFlow=false) const{
      //return Sample(UnderOverFlow,true);
    }
    //samples randomly assuming that the histogram is a pdf, i.e. normalized to the bin width
    bool SamplePdf(double* axisValues, const bool& UnderOverFlow=false, DLM_Random* RanGen=NULL){
      return Sample(axisValues,UnderOverFlow,false,RanGen);
    }
    double SamplePdf(const bool& UnderOverFlow=false) const{
      //return Sample(UnderOverFlow,false);
    }
    
//UPDATE-OKAY
    //the total integral
    Type Integral(const bool& OverUnderFlow=true){
      Type RESULT = 0;
      INT_ERROR = 0;
      if(!Initialized) {InitWarning(); return RESULT;}
      //for(uint64_t uBin=0; uBin<TotNumBins+2*OverUnderFlow; uBin++){
      //  RESULT += BinValue[uBin];
      //  INT_ERROR += BinError[uBin]*BinError[uBin];
      //}
      for(unsigned ej=0; ej<MaxNumBins+2*OverUnderFlow; ej++){
        RESULT += BinValue[ej];
        INT_ERROR += BinError[ej]*BinError[ej];
      }
      INT_ERROR = sqrt(INT_ERROR);
      return RESULT;
    }

//UPDATE-TO-TEST
    //xMin and xMax have size Dim, and each entry represents the min/max value to which the
    //corresponding dimension should be considered
    Type Integral(const double* xMin, const double* xMax, const bool& Normalized=false){
        Type RESULT = 0;
        INT_ERROR = 0;
        if(!Initialized) {InitWarning(); return RESULT;}
        unsigned* BinId = new unsigned [Dim];
        double PhaseSpaceSize = 1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
          PhaseSpaceSize *= (xMax[sDim]-xMin[sDim]);
        }
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            uint64_t ei = ismid(ej);
            if(ei==-1) continue;            
            GetBinCoordinates(ei,BinId);
            double FractionInside=1;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                const double& x_min = xMin[sDim];
                const double& x_max = xMax[sDim];
                double& b_low = BinRange[sDim][BinId[sDim]];
                double& b_up = BinRange[sDim][BinId[sDim]+1];
                //in case both xMin and xMax are outside of the bin
                if(x_min<b_low && x_max>b_up){
                  continue;
                }
                //in case both xMin and xMax are within the bin
                else if(x_min>=b_low && x_max<=b_up){
                    FractionInside *= ( (xMax[sDim]-xMin[sDim])/GetBinSize(sDim,BinId[sDim]) );
                }
                //in case xMax is within the bin, xMin is outside
                else if(x_min<b_low && x_max>=b_low && x_max<=b_up){
                    FractionInside *= ( (xMax[sDim]-BinRange[sDim][BinId[sDim]])/GetBinSize(sDim,BinId[sDim]) );
                }
                //in case xMin is within the bin, xMax is outside
                else if(x_min>=b_low && x_min<=b_up && x_max>b_up){
                    FractionInside *= ( (BinRange[sDim][BinId[sDim]+1]-xMin[sDim])/GetBinSize(sDim,BinId[sDim]) );
                }
                else{
                    FractionInside = 0;
                    break;
                }
            }
            INT_ERROR+=FractionInside*BinError[ej]*(Normalized?GetBinSize(ei):1.);
            RESULT+=FractionInside*BinValue[ej]*(Normalized?GetBinSize(ei):1.);
        }//ej
        delete [] BinId;
        INT_ERROR /= (Normalized?PhaseSpaceSize:1.);
        return RESULT/(Normalized?PhaseSpaceSize:1.);
    }
//UPDATE-OKAY
    Type IntegralError(){
      return INT_ERROR;
    }

//UPDATE-TO-TEST
    //WORK IN PROGRESS
    //does not extrapolate on the edges, takes full bins
    Type FastIntegral(const double* xMin, const double* xMax, const bool& Normalized=false){
        Type RESULT = 0;
        if(!Initialized) {InitWarning(); return RESULT;}
        unsigned* BinId = new unsigned [Dim];
        Type VALUE;
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            uint64_t ei = ismid(ej);
            if(ei==-1) continue;
            GetBinCoordinates(ei,BinId);
            bool OutsideRange = true;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
              OutsideRange *= (BinCenter[sDim][BinId[sDim]]<xMin[sDim] || BinCenter[sDim][BinId[sDim]]>xMax[sDim]);
            }
            if(OutsideRange) continue;
            VALUE = BinValue[ej];
            RESULT+=(BinValue[ej]*(Normalized?GetBinSize(ei):1));
        }//ej
        delete [] BinId;
        return RESULT;
    }


    //the following two rebin function are based on projecting the old histogram onto the new bins.
    //it is rather slow operation, thus not advisable for very large or multidimensional histograms
    //honestly, this function makes no sence, better drop it for now
    //for the future, what may be useful, is to copy other into THIS keeping the binning of THIS (SetEqualTo seems to be doing that...)
//UPDATE-NEXT    
    bool Rebin(DLM_Histo<Type>& other, const bool& Normalized=false){
      if(other.Dim!=Dim){
        printf("\033[1;31mERROR:\033[0m DLM_Histo cannot perform the rebin (wrong dimension of the output)!\n");
        return false;
      }
      double* xMin = new double [Dim];
      double* xMax = new double [Dim];
      unsigned* BinId = new unsigned [Dim];
      for(unsigned uNewBin=0; uNewBin<other.TotNumBins; uNewBin++){
        other.GetBinCoordinates(uNewBin,BinId);
        for(unsigned short sDim=0; sDim<Dim; sDim++){
          xMin[sDim] = other.BinRange[sDim][BinId[sDim]];
          xMax[sDim] = other.BinRange[sDim][BinId[sDim]+1];
        }
        other.BinValue[uNewBin] = Integral(xMin,xMax,Normalized);
//printf("%u in %f : %f = %f\n",uNewBin,xMin[0],xMax[0],other.BinValue[uNewBin]);

        other.BinError[uNewBin] = INT_ERROR;
      }
      other.CumUpdated = false;
      delete [] xMin;
      delete [] xMax;
      delete [] BinId;
      return true;
    }

//UPDATE-NEXT  
    //in each dim
    //not really tested
    void Rebin(const unsigned* RebFactor){
      DLM_Histo<Type> RebinnedHisto;
      RebinnedHisto.SetUp(Dim);
      for(unsigned short sDim=0; sDim<Dim; sDim++){
        double* BinRangeNew = new double[NumBins[sDim]/RebFactor[sDim]+2];
        unsigned BinCounter = 0;
        for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
          if(uBin%RebFactor[sDim]==0){
            BinRangeNew[BinCounter] = BinRange[sDim][uBin];
            BinCounter++;
          }
        }
        BinRangeNew[BinCounter] = GetUpEdge(sDim);
        RebinnedHisto.SetUp(sDim,BinCounter,BinRangeNew);
        delete [] BinRangeNew;
      }
      RebinnedHisto.Initialize();
      Rebin(RebinnedHisto, false);
      Copy(RebinnedHisto);
    }
    unsigned short GetDim() const{
        return Dim;
    }
//UPDATE-NEXT  
    void Rebin(const unsigned dimX, const unsigned dimY=0, const unsigned dimZ=0){
        unsigned* reb_fact = new unsigned [Dim];

        if(dimX==0) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}
        if(dimY==0 && dimZ!=0) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}

        if(dimY==0 && Dim>1) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}
        if(dimY!=0 && Dim==1) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}

        if(dimZ==0 && Dim>2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}
        if(dimZ!=0 && Dim==2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Rebin(...) function failed, check your arguments!\n"); return;}

        reb_fact[0] = dimX;
        if(Dim>1) reb_fact[1] = dimY;
        if(Dim>2) reb_fact[2] = dimZ;

        Rebin(reb_fact);

        delete [] reb_fact;
    }
//UPDATE-NEXT  
    void RebinAll(const unsigned RebFactor){
      unsigned* reb_fact = new unsigned [Dim];
      for(unsigned short sDim=0; sDim<Dim; sDim++){
        reb_fact[sDim] = RebFactor;
      }
      Rebin(reb_fact);
      delete [] reb_fact;
    }





//A FRESH ATTEMPT OT MAKE REBIN A USEFUL FUNCTION (2026)
/*
    //rebin a single dimension by a factor
    //we have to decide what to do if NumBins[sDim]%RebFactor != 0
    bool Rebin(const unsigned sDim, const unsigned RebFactor){
        if(!Initialized) {InitWarning(); return false;}
        if(sDim>=Dim) return false;
        if(!RebFactor) return false;
        //nothing to do
        if(NumBins[sDim]/RebFactor <= 1) return true;



        return true;
    }

    //rebin all dims by the same factor
    void Rebin(const unsigned RebFactor){

    }
    void RebinX(const unsigned RebFactor){

    }  
    void RebinY(const unsigned RebFactor){

    } 
    void RebinZ(const unsigned RebFactor){

    }
    void RebinXY(const unsigned RebFactorX, const unsigned RebFactorY){

    }
    void RebinXYZ(const unsigned RebFactorX, const unsigned RebFactorY, const unsigned RebFactorZ){

    }
    //rebins this based on other.
    //we need to have the same dimensions
    void Rebin(DLM_Histo<Type>& other){

    }

    //a function that copies other into THIS, by keeping the binning of this the same (interpolation based) 
    void GentleCopy(DLM_Histo<Type>& other){

    }
*/



//UPDATE-OKAY
    uint64_t GetNbins() const{
        if(!Initialized) {InitWarning(); return 0;}
        return TotNumBins;
    }
//UPDATE-OKAY
    unsigned GetMaxNbins() const{
        if(!Initialized) {InitWarning(); return 0;}
        return MaxNumBins;
    }
//UPDATE-OKAY
    unsigned GetNbins(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        return NumBins?NumBins[sDim]:0;
    }

    //grep -rnF --include=\*.{cpp,h,hpp,c} "unsigned(0)"
    //grep -rlF --include=\*.{cpp,h,hpp,c} "unsigned(0)" . | xargs -r sed -i 's/unsigned(0)/+0U/g'
    //N.B. even if the Val == 0 and the matrix sparse, the value WILL be saved an use memory.
    //the reason is a possible non-zero error.
//UPDATE-OKAY
    void SetBinContent(const uint64_t& WhichTotBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;
//if(WhichTotBin>=22770000) printf(" setter: %lu\n",WhichTotBin);
        unsigned ej = smid(WhichTotBin,true);

        if(ej==-1) return;
        static bool err_sparse_mem = false;
//printf("WhichTotBin %lu; ej %u; MaxNumBins %u\n",WhichTotBin,ej,MaxNumBins);
        //if(ej==-1){
        //    printf("\033[1;31mERROR:\033[0m DLM_Histo::SetBinContent(...) reports a bug, please contact the developer.\n");
        //    return;
        //}
        if(ej>=MaxNumBins){
            //printf("WhichTotBin %lu; ej %u; MaxNumBins %u\n",WhichTotBin,ej,MaxNumBins);
            if(!err_sparse_mem) printf("\033[1;31mERROR:\033[0m The sparse DLM_Histo ran out of memory. No bins will be set from now on. Try to change the SetUpSparse() settings.\n"); 
            err_sparse_mem = true;
            return;
        }
//if(WhichTotBin>=22770000) printf("  ej: %lu\n",ej);
        //if(smid(WhichTotBin) == -1){
            //BucketId.emplace(WhichTotBin,BucketId.size());
            //bi,bj (key,value)
        //    BucketId.try_emplace(WhichTotBin/NumElPerSparseBucket,BucketId.size());
        //}
        //else{

            //BinValue[WhichTotBin]=Val;
        BinValue[ej]=Val;
        //}
        CumUpdated = false;
    }
//UPDATE-OKAY
    void SetBinContent(const unsigned& WhichX, const unsigned& WhichY, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(Dim!=2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo SetBinContent function failed, this set up works only for Dim=2!\n"); return;}
        if(WhichX>=NumBins[0]) return;
        if(WhichY>=NumBins[1]) return;
        unsigned WhichBin[2];
        WhichBin[0]=WhichX;
        WhichBin[1]=WhichY;
        SetBinContent(GetTotBin(WhichBin),Val);
    }
//UPDATE-OKAY
    void SetBinContent(const unsigned& WhichX, const unsigned& WhichY, const unsigned& WhichZ, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(Dim!=3) {printf("\033[1;33mWARNING:\033[0m DLM_Histo SetBinContent function failed, this set up works only for Dim=3!\n"); return;}
        if(WhichX>=NumBins[0]) return;
        if(WhichY>=NumBins[1]) return;
        if(WhichZ>=NumBins[2]) return;
        unsigned WhichBin[3];
        WhichBin[0]=WhichX;
        WhichBin[1]=WhichY;
        WhichBin[2]=WhichZ;
        SetBinContent(GetTotBin(WhichBin),Val);
    }
//UPDATE-OKAY
     void SetBinError(const uint64_t& WhichTotBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;

        unsigned ej = smid(WhichTotBin,true);
        if(ej==-1) return;
        //static bool err_sparse_mem = false;
        //if(ej==-1){
        //    printf("\033[1;31mERROR:\033[0m DLM_Histo::SetBinError(...) reports a bug, please contact the developer.\n");
        //    return;
        //}
        //if(ej>=MaxNumBins){
        //    if(!err_sparse_mem) printf("\033[1;31mERROR:\033[0m The sparse DLM_Histo ran out of memory. No bins will be set from now on. Try to change the SetUpSparse() settings.\n"); 
        //    err_sparse_mem = true;
        //    return;
        //}

        BinError[ej]=fabs(Val);
    }
//UPDATE-OKAY
    void SetBinError(const unsigned& WhichX, const unsigned& WhichY, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(Dim!=2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo SetBinError function failed, this set up works only for Dim=2!\n"); return;}
        if(WhichX>=NumBins[0]) return;
        if(WhichY>=NumBins[1]) return;
        unsigned WhichBin[2];
        WhichBin[0]=WhichX;
        WhichBin[1]=WhichY;
        SetBinError(GetTotBin(WhichBin),Val);
    }
//UPDATE-OKAY
    void SetBinError(const unsigned& WhichX, const unsigned& WhichY, const unsigned& WhichZ, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        if(Dim!=3) {printf("\033[1;33mWARNING:\033[0m DLM_Histo SetBinError function failed, this set up works only for Dim=2!\n"); return;}
        if(WhichX>=NumBins[0]) return;
        if(WhichY>=NumBins[1]) return;
        if(WhichZ>=NumBins[2]) return;
        unsigned WhichBin[3];
        WhichBin[0]=WhichX;
        WhichBin[1]=WhichY;
        WhichBin[2]=WhichZ;
        SetBinError(GetTotBin(WhichBin),Val);
    }
//UPDATE-OKAY
    void SetBinContentAll(const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        for(uint64_t ei=0; ei<TotNumBins; ei++){
            SetBinContent(ei,Val);
        }
        CumUpdated = false;
    }
//UPDATE-OKAY
    void SetBinErrorAll(const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        for(uint64_t ei=0; ei<TotNumBins; ei++){
            SetBinError(ei,Val);
        }
    }

//UPDATE-TO-TEST
    void Add(const uint64_t& WhichTotBin, const Type& Val, const Type& Err=0){
        if(!Initialized) {InitWarning(); return;}
        if(WhichTotBin>=TotNumBins) return;
        unsigned ej = smid(WhichTotBin,true);
        if(ej==-1) return;
        //if(ej==-1){
        //    printf("\033[1;31mERROR:\033[0m DLM_Histo::Add(...) reports a bug, please contact the developer.\n");
        //    return;
        //}
        BinValue[ej]+=Val;
        BinError[ej]=Type(sqrt(double(BinError[ej]*BinError[ej])+double(Err*Err)));
        CumUpdated = false;
    }
//UPDATE-TO-TEST
    void AddAt(const double* AxisValue, const Type& Val=1){
        if(!Initialized) {InitWarning(); return;}
        unsigned* WhichBin = new unsigned [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            WhichBin[sDim] = GetBin(sDim,AxisValue[sDim]);
        }
        uint64_t ei = GetTotBin(WhichBin);
        if(ei>=TotNumBins || ei==-1) return;
        unsigned ej = smid(ei,true);
        if(ej==-1) {delete [] WhichBin; return;}
        BinValue[ej] += Val;
        delete [] WhichBin;
        CumUpdated = false;
    }
//UPDATE-OKAY
    void AddAt(const double& xValue, const Type& Val){
      if(!Initialized) {InitWarning(); return;}
      if(Dim!=1) {printf("\033[1;33mWARNING:\033[0m DLM_Histo AddAt/Fill function failed, this set up works only for Dim=1!\n"); return;}
      AddAt(&xValue,Val);
    }
//UPDATE-OKAY
    void AddAt(const double& xValue, const double& yValue, const Type& Val){
      if(!Initialized) {InitWarning(); return;}
      if(Dim!=2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo AddAt/Fill function failed, this set up works only for Dim=2!\n"); return;}
      double axis[2];
      axis[0] = xValue;
      axis[1] = yValue;
      AddAt(axis,Val);
    }
//UPDATE-OKAY
    void AddAt(const double& xValue, const double& yValue, const double& zValue, const Type& Val){
      if(!Initialized) {InitWarning(); return;}
      if(Dim!=3) {printf("\033[1;33mWARNING:\033[0m DLM_Histo AddAt/Fill function failed, this set up works only for Dim=3!\n"); return;}
      double axis[3];
      axis[0] = xValue;
      axis[1] = yValue;
      axis[2] = zValue;
      AddAt(axis,Val);
    }
//UPDATE-OKAY
    void Fill(const double& xValue){
      AddAt(xValue,1);
    }
//UPDATE-OKAY
    void Fill(const double& xValue, const double& yValue){
      AddAt(xValue,yValue,1);
    }
//UPDATE-OKAY
    void Fill(const double& xValue, const double& yValue, const double& zValue){
      AddAt(xValue,yValue,zValue,1);
    }
//UPDATE-TO-TEST    
    void SetBinContent(const unsigned* WhichBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        uint64_t ei = GetTotBin(WhichBin);
        if(ei>=TotNumBins || ei==-1) return;
        SetBinContent(ei,Val);
        CumUpdated = false;
    }
//UPDATE-TO-TEST   
    void SetBinError(const unsigned* WhichBin, const Type& Val){
        if(!Initialized) {InitWarning(); return;}
        uint64_t ei = GetTotBin(WhichBin);
        if(ei>=TotNumBins || ei==-1) return;
        SetBinError(ei,Val);
    }
//UPDATE-TO-TEST
    void ComputeError(){
        if(!Initialized) {InitWarning(); return;}
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            BinError[ej] = sqrt(BinValue[ej]);
//printf("BinValue[%u]=%f+/-%f\n",uBin,BinValue[uBin],BinError[uBin]);
        }
    }
//UPDATE-OKAY
    void SetBinCenter(const unsigned short& sDim, const unsigned& WhichBin, const double& Val){
        if(!Initialized) {InitWarning(); return;}
        if(sDim>=Dim)return;
        if(WhichBin>=NumBins[sDim])return;
        if(!BinCenter)return;
        if(!BinCenter[sDim])return;
        BinCenter[sDim][WhichBin] = Val;
    }
//UPDATE-OKAY
    void SetBinCenter(const unsigned short& sDim, const double* bincenter){
        if(!Initialized) {InitWarning(); return;}
        if(sDim>=Dim)return;
        if(!BinCenter)return;
        if(!BinCenter[sDim])return;
        for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
          BinCenter[sDim][uBin] = bincenter[uBin];
        }
    }
//UPDATE-OKAY
    double GetBinCenter(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim)return 0;
        if(WhichBin>=NumBins[sDim])return 0;
        if(!BinCenter)return 0;
        if(!BinCenter[sDim])return 0;
        return BinCenter[sDim][WhichBin];
    }
//UPDATE-OKAY
    Type GetBinContent(const uint64_t& WhichTotBin) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(WhichTotBin>=TotNumBins+2) return NULL_ELEM;
        unsigned ej = smid(WhichTotBin);
        return ej==-1?NULL_ELEM:BinValue[ej];
    }
//UPDATE-OKAY
    Type GetBinContent(const unsigned& WhichX, const unsigned& WhichY) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(Dim!=2)return NULL_ELEM;
        unsigned WhichBin[2];
        WhichBin[0] = WhichX;
        WhichBin[1] = WhichY;
        return GetBinContent(GetTotBin(WhichBin));
    }
//UPDATE-OKAY
    Type GetBinContent(const unsigned& WhichX, const unsigned& WhichY, const unsigned& WhichZ) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(Dim!=3)return NULL_ELEM;
        unsigned WhichBin[3];
        WhichBin[0] = WhichX;
        WhichBin[1] = WhichY;
        WhichBin[2] = WhichZ;
        return GetBinContent(GetTotBin(WhichBin));
    }
//UPDATE-OKAY
    Type GetBinContent(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        return GetBinContent(GetTotBin(WhichBin));
    }
//UPDATE-OKAY
    Type GetBinError(const uint64_t& WhichTotBin) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(WhichTotBin>=TotNumBins) return NULL_ELEM;
        unsigned ej = smid(WhichTotBin);
        return ej==-1?NULL_ELEM:BinError[ej];        
    }
//UPDATE-OKAY
    Type GetBinError(const unsigned* WhichBin) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        return GetBinError(GetTotBin(WhichBin));
    }
//UPDATE-OKAY
    Type GetBinError(const unsigned& WhichX, const unsigned& WhichY) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(Dim!=2)return NULL_ELEM;
        unsigned WhichBin[2];
        WhichBin[0] = WhichX;
        WhichBin[1] = WhichY;
        return GetBinError(GetTotBin(WhichBin));
    }
//UPDATE-OKAY
    Type GetBinError(const unsigned& WhichX, const unsigned& WhichY, const unsigned& WhichZ) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
        if(Dim!=3)return NULL_ELEM;
        unsigned WhichBin[3];
        WhichBin[0] = WhichX;
        WhichBin[1] = WhichY;
        WhichBin[2] = WhichZ;
        return GetBinError(GetTotBin(WhichBin));
    }
    //Type GetBinError(const unsigned short& sDim, const unsigned& WhichBin) const{
    //    return GetBinError(sDim,GetBin(sDim,WhichBin));
    //}

//UPDATE-OKAY
    double GetLowEdge(const unsigned short& sDim) const{
        return GetBinLowEdge(sDim,0);
    }
//UPDATE-OKAY
    double GetUpEdge(const unsigned short& sDim) const{
        return GetBinUpEdge(sDim,NumBins[sDim]-1);
    }
//UPDATE-OKAY
    double GetBinLowEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        //i.e. GetBinLowEdge(NumBins) will return the most upper limit
        if(WhichBin>=NumBins[sDim]+1) return 0;
        return BinRange[sDim][WhichBin];
    }
//UPDATE-OKAY
    double GetBinUpEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
        return BinRange[sDim][WhichBin+1];
    }
//UPDATE-OKAY
    double GetBinSize(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
//printf("sDim=%u; WhichBin=%u; BR1=%f; BR0=%f\n",sDim,WhichBin,BinRange[sDim][WhichBin+1],BinRange[sDim][WhichBin]);
        return (BinRange[sDim][WhichBin+1]-BinRange[sDim][WhichBin]);
    }
//UPDATE-OKAY
    double GetBinSize(const uint64_t& WhichTotBin) const{
        if(WhichTotBin>=TotNumBins) return 0;
        double BinSize = 1;
        unsigned* WhichBin = new unsigned [Dim];
        GetBinCoordinates(WhichTotBin,WhichBin);
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            BinSize *= GetBinSize(sDim,WhichBin[sDim]);
        }
        delete [] WhichBin;
        return BinSize;
    }
//UPDATE-OKAY
    unsigned FindBin(const unsigned short& sDim, const double& xVal) const{
        return GetBin(sDim, xVal);
    }
//UPDATE-OKAY
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
/*
    unsigned FindBin(const double* xVal) const{
        return GetBin(xVal);
    }
    unsigned GetBin(const double* xVal) const{
        if(!Initialized) {InitWarning(); return 0;}
        unsigned* WhichBin = new unsigned [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            WhichBin[sDim] = GetBin(sDim, xVal[sDim]);
        }
        delete [] WhichBin;
    }
*/
//UPDATE-OKAY
    double* GetBinRange(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        double* BINRANGE = new double [NumBins[sDim]+1];
        for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
            BINRANGE[uBin] = BinRange[sDim][uBin];
        }
        return BINRANGE;
    }
//UPDATE-OKAY
    double* GetBinCenters(const unsigned short& sDim) const{
        if(!Initialized) {InitWarning(); return 0;}
        if(sDim>=Dim) return 0;
        double* BINCENTER = new double [NumBins[sDim]+1];
        for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
            BINCENTER[uBin] = BinCenter[sDim][uBin];
        }
        return BINCENTER;
    }
    //returns the number of minima
//UPDATE-TO-TEST
    unsigned FindMinima(Type& Value) const{
        uint64_t* Dummy;
        unsigned NumMin=FindMinima(Value,Dummy);
        delete[]Dummy;
        return NumMin;
    }
//UPDATE-TO-TEST
    unsigned FindMinima(Type& Value, uint64_t*& BinId) const{
        bool MinFound=false;
        Type MIN_VAL;
        unsigned NumMin;
        unsigned SizeOfBuffer=8;
        uint64_t* WhichTotBins = new uint64_t[SizeOfBuffer];
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            uint64_t ei = ismid(ej);
            if(ei==-1) continue;
            if(!ej){
                NumMin=0;
                MIN_VAL=BinValue[ej];
                WhichTotBins[NumMin]=ei;
                NumMin=1;
            }
            else if(BinValue[ej]<MIN_VAL){
                NumMin=0;
                MIN_VAL=BinValue[ej];
                WhichTotBins[NumMin]=ei;
                NumMin=1;
            }
            else if(BinValue[ej]==MIN_VAL){
                if(SizeOfBuffer<=NumMin){
                    SizeOfBuffer*=2;
                    uint64_t* TempBins = new uint64_t [SizeOfBuffer];
                    for(unsigned uBuff=0; uBuff<SizeOfBuffer/2; uBuff++)
                        TempBins[uBuff]=WhichTotBins[uBuff];
                    delete[]WhichTotBins;
                    WhichTotBins=TempBins;
                }
                WhichTotBins[NumMin]=ei;
                NumMin++;
            }
        }
        //printf(" WhichTotBins=%p\n",WhichTotBins);
        //printf(" WhichTotBins[0]=%u\n",WhichTotBins[0]);
        BinId=WhichTotBins;
        Value = MIN_VAL;
        return NumMin;
    }
    //returns the number of maxima
//UPDATE-TO-TEST
    unsigned FindMixima(Type& Value) const{
        unsigned* Dummy;
        unsigned NumMax=FindMixima(Value,Dummy);
        delete[]Dummy;
        return NumMax;
    }
//UPDATE-TO-TEST
    Type FindMixima(Type& Value, uint64_t*& BinId) const{
        bool MaxFound=false;
        Type MAX_VAL;
        unsigned NumMax;
        unsigned SizeOfBuffer=8;
        uint64_t* WhichTotBins = new uint64_t[SizeOfBuffer];
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            uint64_t ei = ismid(ej);
            if(ei==-1) continue;
            if(!ej){
                NumMax=0;
                MAX_VAL=BinValue[ej];
                WhichTotBins[NumMax]=ei;
                NumMax=1;
            }
            else if(BinValue[ej]>MAX_VAL){
                NumMax=0;
                MAX_VAL=BinValue[ej];
                WhichTotBins[NumMax]=ei;
                NumMax=1;
            }
            else if(BinValue[ej]==MAX_VAL){
                if(SizeOfBuffer<=NumMax){
                    SizeOfBuffer*=2;
                    uint64_t* TempBins = new uint64_t [SizeOfBuffer];
                    for(unsigned uBuff=0; uBuff<SizeOfBuffer/2; uBuff++)
                        TempBins[uBuff]=WhichTotBins[uBuff];
                    delete[]WhichTotBins;
                    WhichTotBins=TempBins;
                }
                WhichTotBins[NumMax]=ei;
                NumMax++;
            }
        }
        BinId=WhichTotBins;
        Value = MAX_VAL;
        return NumMax;
    }
//UPDATE-TO-TEST (not changed, but I think its okay)
    Type Eval(const double* xVal, const bool& EvalTheError=false) const{
        if(!Initialized) {InitWarning(); return NULL_ELEM;}
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
          if(xVal[sDim]>GetUpEdge(sDim) || xVal[sDim]<GetLowEdge(sDim)){
            delete [] xValue1;
            delete [] xValue2;
            delete [] xBin1;
            delete [] xBin2;
            delete [] DeltaX1;
            delete [] DeltaX2;
            return NULL_ELEM;
          }
        }

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
        double Weight=1;
        double Norm=0;
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
            if(EvalTheError){Result += (GetBinError(BinArray)*Weight);}
            else{Result += (GetBinContent(BinArray)*Weight);}
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
//UPDATE-OKAY
    Type EvalError(const double* xVal) const{
        return Eval(xVal,true);
    }
//UPDATE-OKAY
    Type Eval(const double xVal, const bool& EvalTheError=false) const{
      if(Dim!=1) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Eval(xVal) function failed, this set up works only for Dim=1!\n"); return NULL_ELEM;}
      return Eval(&xVal,EvalTheError);
    }
//UPDATE-OKAY
    Type EvalError(const double xVal) const{
      if(Dim!=1) {printf("\033[1;33mWARNING:\033[0m DLM_Histo EvalError(xVal) function failed, this set up works only for Dim=1!\n"); return NULL_ELEM;}
      return EvalError(&xVal);
    }
//UPDATE-OKAY
    Type Eval2D(const double xVal, const double yVal, const bool& EvalTheError=false) const{
      if(Dim!=2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Eval2D(xVal) function failed, this set up works only for Dim=2!\n"); return NULL_ELEM;}
      double xPars[2];
      xPars[0] = xVal;
      xPars[1] = yVal;
      return Eval(xPars,EvalTheError);
    }
//UPDATE-OKAY
    Type EvalError2D(const double xVal, const double yVal) const{
      if(Dim!=2) {printf("\033[1;33mWARNING:\033[0m DLM_Histo EvalError2D(xVal) function failed, this set up works only for Dim=2!\n"); return NULL_ELEM;}
      double xPars[2];
      xPars[0] = xVal;
      xPars[1] = yVal;
      return EvalError(xPars);
    }    
//UPDATE-OKAY
    Type Eval3D(const double xVal, const double yVal, const double zVal, const bool& EvalTheError=false) const{
      if(Dim!=3) {printf("\033[1;33mWARNING:\033[0m DLM_Histo Eval3D(xVal) function failed, this set up works only for Dim=3!\n"); return NULL_ELEM;}
      double xPars[3];
      xPars[0] = xVal;
      xPars[1] = yVal;
      xPars[2] = zVal;
      return Eval(xPars,EvalTheError);
    }
//UPDATE-OKAY
    Type EvalError3D(const double xVal, const double yVal, const double zVal) const{
      if(Dim!=3) {printf("\033[1;33mWARNING:\033[0m DLM_Histo EvalError3D(xVal) function failed, this set up works only for Dim=3!\n"); return NULL_ELEM;}
      double xPars[3];
      xPars[0] = xVal;
      xPars[1] = yVal;
      xPars[2] = zVal;
      return EvalError(xPars);
    }   


//UPDATE-TO-TEST
    //makes this histogram equal in values to other, without re-initialization.
    //This only works if we have the same Dim.
    //Only overlapping bin regions are set equal (using Eval), the rest are zero
    bool SetEqualTo(const DLM_Histo& other){
        if(Dim!=other.Dim){
            printf("\033[1;33mWARNING:\033[0m DLM_Histo::SetEqualTo function failed, this set up works only for histos with the same Dim.\n");
            return false;
        }

        if(SameStructure(other)){
            //we set all bins to zero
            for(unsigned ej=0; ej<MaxNumBins; ej++){
                if(ismid(ej)){
                    BinValue[ej] = NULL_ELEM;
                    BinError[ej] = NULL_ELEM;
                }
            }

            //we iterate over the other bins, and set them back in our histogram
            for(unsigned oej=0; oej<other.MaxNumBins; oej++){
                //due to same structure
                uint64_t ei = other.ismid(oej);
                if(ei==-1) continue;
                SetBinContent(ei, other.GetBinContent(ei));
                SetBinError(ei, other.GetBinError(ei));
            }
            
            return true;
        }

        //different structure

        //this loop essentially makes sure that we update the bins that already exist in THIS
        for(unsigned ej=0; ej<MaxNumBins; ej++){
            //printf("uBin %u\n",uBin);
            unsigned* WhichTotBin = new unsigned[Dim];
            uint64_t ei = ismid(ej);
            if(ei==-1){
                continue;
            }
            GetBinCoordinates(ei, WhichTotBin);
            double* xVal = new double [Dim];
            bool InRange = true;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][WhichTotBin[sDim]];
                if(xVal[sDim]<other.GetLowEdge(sDim) || xVal[sDim]>other.GetUpEdge(sDim)){
                    InRange = false;
                    break;
                }
            }
            //printf("uBin=%u(%u) %i; xVal[0]=%.3f\n",uBin,TotNumBins,uBin>=TotNumBins,xVal[0]);
            if(InRange){
                BinValue[ej] = other.Eval(xVal);
                BinError[ej] = other.EvalError(xVal);
            }
            else{
                BinValue[ej] = NULL_ELEM;
                BinError[ej] = NULL_ELEM;  
            }

            delete [] WhichTotBin;
            delete [] xVal;
        }

        //this loop adds all elements that did not exist in THIS
        for(unsigned oej=0; oej<other.MaxNumBins; oej++){
            unsigned* WhichTotBin = new unsigned[Dim];
            uint64_t oei = ismid(oej);
            if(oei==-1){
                continue;
            }
            GetBinCoordinates(oei, WhichTotBin);
            double* xVal = new double [Dim];
            bool InRange = true;
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = other.BinCenter[sDim][WhichTotBin[sDim]];
                if(xVal[sDim]<GetLowEdge(sDim) || xVal[sDim]>GetUpEdge(sDim)){
                    InRange = false;
                    break;
                }
            }
            if(!InRange) continue;

            uint64_t ei = FindTotBin(xVal);
            if(ei==-1) continue;
            unsigned ej = smid(ei,true);
            if(ej==-1){
                delete [] WhichTotBin;
                delete [] xVal;                
                return false;
            }
            GetBinCoordinates(ei, WhichTotBin);
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                xVal[sDim] = BinCenter[sDim][WhichTotBin[sDim]];
            }
            SetBinContent(ei, other.Eval(xVal));
            SetBinError(ei, other.EvalError(xVal));
            //BinValue[ej] = other.Eval(xVal);
            //BinError[ej] = other.EvalError(xVal);

            delete [] WhichTotBin;
            delete [] xVal;
        }

        return true;
    }


//UPDATE-NEXT-TO-UPDATE
    bool operator=(const DLM_Histo& other){
        return Copy(other, false);
    }
    bool operator==(const DLM_Histo& other){
        if(!SameStructure(other)) return false;
        for(uint64_t uBin=0; uBin<TotNumBins+2; uBin++){
            if(BinValue[uBin]!=other.BinValue[uBin]) return false;
            if(uBin<TotNumBins){
                if(BinError[uBin]!=other.BinError[uBin]) return false;
                if(BinCenter[uBin]!=other.BinCenter[uBin]) return false;
            }
            if(uBin<TotNumBins+1){
                if(BinRange[uBin]!=other.BinRange[uBin]) return false;
            }
        }
        return true;
    }
//! FOR DIVISION/MULT:
//http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
//https://apps.dtic.mil/dtic/tr/fulltext/u2/785623.pdf
//UPDATE-OKAY
    bool operator+=(const DLM_Histo& other){
        return AddHisto(other,true,1);
    }
//UPDATE-OKAY
    bool operator-=(const DLM_Histo& other){
        return AddHisto(other,true,-1);
    }
//UPDATE-OKAY
    bool operator*=(const DLM_Histo& other){
        return MultiplyHisto(other,true);
    }
//UPDATE-OKAY
    bool operator/=(const DLM_Histo& other){
        return DivideHisto(other,true);
    }

    //note that with this type of construct, DLM_Histo<Type> histo=histo1+histo2 invokes the constructor only ones!
//UPDATE-TO-TEST
    DLM_Histo<Type> operator+(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result+=other;
//printf("Result=%p\n",&Result);
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator-(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result-=other;
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator*(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result*=other;
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator/(const DLM_Histo& other){
        DLM_Histo<Type> Result(*this);
        Result/=other;
        return Result;
    }
//UPDATE-TO-TEST
    bool operator+=(const Type& Value){
        AddToAll(Value);
        return true;
    }
//UPDATE-TO-TEST
    bool operator-=(const Type& Value){
        AddToAll(-Value);
        return true;
    }
//UPDATE-TO-TEST
    bool operator*=(const Type& Value){
        Scale(Value);
        return true;
    }
//UPDATE-TO-TEST
    bool operator/=(const Type& Value){
        Scale(1./Value);
        return true;
    }

//UPDATE-TO-TEST
    DLM_Histo<Type> operator+(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result+=Value;
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator-(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result-=Value;
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator*(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result*=Value;
        return Result;
    }
//UPDATE-TO-TEST
    DLM_Histo<Type> operator/(const Type& Value){
        DLM_Histo<Type> Result(*this);
        Result/=Value;
        return Result;
    }


//UPDATE-TO-TEST
    //write to a binary file of format
    bool QuickWrite(const char* FileName, bool Overwrite=false){
      if(!Initialized){
        printf("\033[1;31mERROR:\033[0m The histogram must be initialized before writing it to a file.\n");
        return false;
      }
      ifstream myFileIN(FileName);
      if(myFileIN.fail()==false && Overwrite==false){
        printf("\033[1;31mERROR:\033[0m The file %s exists. Change name or use QuickWrite(FileName, true) to overwrite.\n",FileName);
        return false;
      }
      myFileIN.close();

      ofstream myFileOUT(FileName, ios::out | ios::binary);
      if(!myFileOUT) {
         printf("\033[1;31mERROR:\033[0m Cannot open file %s.\n",FileName);
         return false;
      }

      //a silly check to make sure we have the correct file fromat
      //an older version which cannot handle sparse matrix. 2783 can still read it (but of course treating it as a non-sparse)
      //short Watermark = 1331;
      //the version with the sparse.
      short Watermark = 2783;
      myFileOUT.write((char *) &Watermark, sizeof(short));

      //WriteVersion:
      // -1 : Dim, TotNumBins, NumBins, BinRange, BinValue, BinError, BinCenter
      //      no overflow bins yet
      WriteVersion = -1;
      myFileOUT.write((char *) &WriteVersion, sizeof(int));
      myFileOUT.write((char *) &Dim, sizeof(unsigned short));
      //myFileOUT.write((char *) &TotNumBins, sizeof(unsigned));//we dont need it, as we can compute it
      for(unsigned short sDim=0; sDim<Dim; sDim++){
        myFileOUT.write((char *) &NumBins[sDim], sizeof(unsigned));
        for(unsigned uBin=0; uBin<NumBins[sDim]+1; uBin++){
          myFileOUT.write((char *) &BinRange[sDim][uBin], sizeof(double));
          if(uBin!=NumBins[sDim])
            myFileOUT.write((char *) &BinCenter[sDim][uBin], sizeof(double));
        }
      }


      myFileOUT.write((char *) &MaxNumBins, sizeof(unsigned));

      //sparse //NEW
      if(TotNumBins>MaxNumBins){
        myFileOUT.write((char *) &NumElPerSparseBucket, sizeof(unsigned));
        auto it = BucketId.begin();
        unsigned NumSavedBuckets = BucketOrder.size();
        myFileOUT.write((char *) &NumSavedBuckets, sizeof(unsigned));
        for(unsigned ub=0; ub<NumSavedBuckets; ub++, ++it){
            uint64_t it_first = it->first;
            unsigned it_second = it->second;
            uint64_t bucket_order = BucketOrder.at(ub);
            myFileOUT.write((char *) &it_first, sizeof(uint64_t));
            myFileOUT.write((char *) &it_second, sizeof(unsigned));
            myFileOUT.write((char *) &bucket_order, sizeof(uint64_t));
        }
      }

      //myFileOUT.write((char *) &MaxNumBins, sizeof(unsigned));//NEW
      for(unsigned ej=0; ej<MaxNumBins+2; ej++){
        myFileOUT.write((char *) &BinValue[ej], sizeof(Type));
        //if(TotNumBins>MaxNumBins){
        //    if( ej<100 ){
        //        printf("%e -- %e\n",BinValue[ej], GetBinContent(ej));
        //    }
        //}
      }
      for(unsigned ej=0; ej<MaxNumBins; ej++){
        myFileOUT.write((char *) &BinError[ej], sizeof(Type));
      }

      myFileOUT.close();

      return true;
    }
//UPDATE-TO-TEST
    bool QuickLoad(const char* FileName, const int Version = -1){
      if(Initialized){
        printf("\033[1;31mERROR:\033[0m Cannot load from file, the target object is already initialized (it should be blank).\n");
        return false;
      }
      if(Version!=-1){
        printf("\033[1;31mERROR:\033[0m Unknown file version (%i).\n",Version);
        return false;
      }

      ifstream myFileIN(FileName);
      if(myFileIN.fail()){
        printf("\033[1;31mERROR:\033[0m The file %s cannot be opened.\n",FileName);
        return false;
      }

      short Watermark;
      myFileIN.read ((char*) &Watermark,sizeof(short));
      if(Watermark!=1331 && Watermark!=2783){
        printf("\033[1;31mERROR:\033[0m Issue with the file format of %s\n",FileName);
        myFileIN.close();
        return false;
      }

      int WriteVersion;
      myFileIN.read((char *) &WriteVersion, sizeof(int));

      unsigned short dim;
      myFileIN.read((char *) &dim, sizeof(unsigned short));
      SetUp(dim);

      unsigned numbins;
      uint64_t totnumbins=1;
      double* binrange = NULL;
      double* bincenter = NULL;
      for(unsigned short sDim=0; sDim<Dim; sDim++){
        myFileIN.read((char *) &numbins, sizeof(unsigned));
        binrange = new double [numbins+1];
        bincenter = new double [numbins];
        for(unsigned uBin=0; uBin<numbins+1; uBin++){
          myFileIN.read((char *) &binrange[uBin], sizeof(double));
          if(uBin!=numbins){
            myFileIN.read((char *) &bincenter[uBin], sizeof(double));
          }
        }
        SetUp(sDim,numbins,binrange,bincenter);
        totnumbins *= numbins;
        delete [] binrange; binrange = NULL;
        delete [] bincenter; bincenter = NULL;
      }
      if(Watermark==2783){
        myFileIN.read((char *) &MaxNumBins, sizeof(unsigned));
      }
      else{
        MaxNumBins = DefaultMaxNumBins;
      }
      //can happen in older versions, should NOT be the case though
      if(totnumbins>DefaultMaxNumBins){
        printf("\033[1;31mERROR:\033[0m DLM_Histo::QuickLoad The histogram in %s cannot be initizlized. Most likely you are loading an old histogram with a CATS version after October 2026, which has a hard limit on the maximum number of bins allowed.\n",FileName);
        myFileIN.close();
        return false;
      }

      //sparse //NEW
      if(totnumbins>MaxNumBins && Watermark==2783){
        myFileIN.read((char *) &NumElPerSparseBucket, sizeof(unsigned));
        uint64_t it_first;
        unsigned it_second;
        uint64_t bucket_order;

        unsigned NumSavedBuckets;
        myFileIN.read((char *) &NumSavedBuckets, sizeof(unsigned));

        Initialize();

        for(unsigned ub=0; ub<NumSavedBuckets; ub++){
            myFileIN.read((char *) &it_first, sizeof(uint64_t));
            myFileIN.read((char *) &it_second, sizeof(unsigned));
            myFileIN.read((char *) &bucket_order, sizeof(uint64_t));
            BucketId.try_emplace(it_first,it_second);
            BucketOrder.push_back(bucket_order);
        }
      }
      else{
        Initialize();  
      }

      for(unsigned ej=0; ej<MaxNumBins+2; ej++){
        myFileIN.read((char *) &BinValue[ej], sizeof(Type));
        //if(ej<20) printf("%e\n",BinValue[ej]);
      }
      for(uint64_t ej=0; ej<MaxNumBins; ej++){
        myFileIN.read((char *) &BinError[ej], sizeof(Type));
      }

      myFileIN.close();
      return true;
    }


protected:
    void ConstructorState(){
//printf("ConstructorState\n");
        BinRange = NULL;
        BinValue = NULL;
        BinError = NULL;
        BinCenter = NULL;
        CumulativeValue = NULL;
        NumBins = NULL;
        TotNumBins=0;
        PER = NULL;
        Initialized=false;
        Dim = 0;
        rangen = NULL;
        NULL_ELEM = Type(0);
        if constexpr (std::numeric_limits<Type>::has_quiet_NaN){
            NAN_ELEM = std::numeric_limits<Type>::quiet_NaN();
        }
        if constexpr (std::numeric_limits<Type>::has_infinity) {
            INF_ELEM = std::numeric_limits<Type>::infinity();
        }
        else{
            NAN_ELEM = Type(0);
        }
        MaxNumBins = DefaultMaxNumBins;
    }
    void CleanUp(){
//printf("CleanUp %u dimensions\n",Dim);
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
        if(rangen){
          delete rangen;
          rangen = NULL;
        }
        if(CumulativeValue){delete [] CumulativeValue; CumulativeValue=NULL;}
        if(NumBins){delete [] NumBins; NumBins=NULL;}
        BucketId.clear();
        BucketOrder.clear();
        ConstructorState();
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
    //does NOT check if both are setup as sparse matrix or not
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
//UPDATE-TO-TEST
    //updates a single cumulative bin, assuming that WhichTotBin-1 is updated, as well as BinValue[WhichTotBin]
    void UpdateCum(const bool& yield){
      if(!Initialized) {InitWarning(); return;}
      if(!CumulativeValue) CumulativeValue = new Type[MaxNumBins+2];
      for(unsigned ej=0; ej<MaxNumBins+2; ej++){
        double Renorm = 1;
        uint64_t ei = ismid(ej);
        //if the bin is not set, we consider it zero (sparse)
        //N.B. ALL possible ej values will be set for Cumulative disto, even if the corresponding bins are missing
        if(ei==-1){
            if(ej>0) CumulativeValue[ej] = CumulativeValue[ej-1];
            else CumulativeValue[ej] = NULL_ELEM;
            continue;
        }
        //if this histo is a pdf, we need to multiply by the bin width in order to sample properly
        if(ej<MaxNumBins && yield==false){
          Renorm = GetBinSize(ei);//this is the bin size (N-dim volume, width for 1D)
        }
        if(ej==0) CumulativeValue[ej]=BinValue[ej]*Renorm;
        else CumulativeValue[ej]=CumulativeValue[ej-1]+BinValue[ej]*Renorm;
        //printf("b%u %e %e\n",uBin,CumulativeValue[uBin],BinValue[uBin]);
        if(BinValue[ej]<0){
          printf("\033[1;33mWARNING:\033[0m DLM_Histo cannot be used for sampling, due to negative entries!\n");
          CumUpdated = false;
          return;
        }
      }
      CumUpdated = true;
    }

    //the UnderOverFlow is NOT tested in the latest version, so be mindful for bugs
    bool Sample(double* axisValues, const bool& UnderOverFlow, const bool& yield, DLM_Random* RanGen){
        if(!RanGen) RanGen = rangen;
        if(!RanGen) RanGen = new DLM_Random(0);
        unsigned TopBin_Glob = UnderOverFlow?MaxNumBins+1:MaxNumBins-1;
        if(!CumulativeValue) UpdateCum(yield);
        //get a random number for the value of the function (histogram)
        double rannum = RanGen->Uniform(0,double(CumulativeValue[TopBin_Glob]));

        unsigned BinCandidate = TopBin_Glob / 2;
        unsigned ej = BinCandidate;
        unsigned LowBin = 0;
        unsigned TopBin = MaxNumBins-1;

        while(LowBin <= TopBin){
            BinCandidate = LowBin + (TopBin - LowBin) / 2;
            if(CumulativeValue[BinCandidate] >= rannum){
                // Check if this is the zeroth element OR if the previous element is < target
                if(BinCandidate == 0 || CumulativeValue[BinCandidate] < rannum){
                    break; // Exact transition point found
                }
                TopBin = BinCandidate - 1; // Value matches, but an earlier valid element exists to the left
            }
            else{
                LowBin = BinCandidate + 1; // arr[mid] < target, so search the right half
            }
        }

        unsigned ei = ismid(BinCandidate);
        //if we are at larger bins, we should in principle always get the first bin that works
        //so in theory, this can happen only if we are underflowing, but than we have literally no good 
        //starting point, as in N-dims we do not know what is really the zero-th bin
        if(ei==-1){
            return false;
        }

        unsigned* BinId = new unsigned [Dim];
        GetBinCoordinates(ei,BinId);
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            //underflow -> return the lowest edge
            if(BinId[sDim]==NumBins[sDim]){
                axisValues[sDim] = BinRange[sDim][0];
            }
            //overflow -> return the upper edge
            else if(BinId[sDim]==NumBins[sDim]+1){
                axisValues[sDim] = BinRange[sDim][NumBins[sDim]];
            }
            else{
                //axisValues[sDim] = BinCenter[sDim][BinId[sDim]];
                axisValues[sDim] = RanGen->Uniform(BinRange[sDim][BinId[sDim]],BinRange[sDim][BinId[sDim]+1]);
            }
        }
        delete [] BinId;
        return true;
    }

    //this is done based on the absolute value of the bin, NOT normalized to bin width
    void Sample_LEGACY(double* axisValues, const bool& UnderOverFlow, const bool& yield, DLM_Random* RanGen){
      if(!RanGen) RanGen = rangen;
      if(!RanGen) RanGen = new DLM_Random(0);
      uint64_t TopBin = UnderOverFlow?TotNumBins+1:TotNumBins-1;
      //unsigned TopBin = UnderOverFlow?MaxNumBins+1:MaxNumBins-1;
      unsigned TopBin_Glob = ismid(TopBin);
      if(!CumulativeValue) UpdateCum(yield);
      double rannum = RanGen->Uniform(0,double(CumulativeValue[TopBin]));
      //printf("rn %f cvtp %f\n",rannum,CumulativeValue[TopBin]);
      double value_down;
      double value_up;
      unsigned Fraction = 2;
      unsigned Step = TopBin/Fraction;
      unsigned BinCandidate = Step;
      while(true){
        //printf("BinCandidate = %u\n",BinCandidate);
        if(BinCandidate==0) value_down = double(CumulativeValue[BinCandidate]);
        else value_down = double(CumulativeValue[BinCandidate]+CumulativeValue[BinCandidate-1])*0.5;
        if(BinCandidate==TopBin) value_up = double(CumulativeValue[BinCandidate]);
        else value_up = double(CumulativeValue[BinCandidate]+CumulativeValue[BinCandidate+1])*0.5;
        //printf(" value_down = %f\n",value_down);
        //printf(" value_up = %f\n",value_up);
        //printf(" rannum = %f\n",rannum);
        //we found our bin
        if(value_up>=rannum && value_down<=rannum){
          //special treatment for under and overflow bins
          if(BinCandidate>=MaxNumBins){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                axisValues[sDim] = std::numeric_limits<double>::quiet_NaN();
            }
            return;
          }

          uint64_t BinCandidate_Glob = ismid(BinCandidate);
          //a zero bin in the sparse... should not really happen
          if(BinCandidate_Glob==-1){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                axisValues[sDim] = std::numeric_limits<double>::quiet_NaN();
            }     
            printf("\033[1;31mERROR:\033[0m DLM_Histo::Sample reports a bug, please contact the developer.\n");       
            return;
          }

          unsigned* BinId = new unsigned [Dim];
          GetBinCoordinates(BinCandidate_Glob,BinId);
          for(unsigned short sDim=0; sDim<Dim; sDim++){
            //axisValues[sDim] = BinCenter[sDim][BinId[sDim]];
            axisValues[sDim] = RanGen->Uniform(BinRange[sDim][BinId[sDim]],BinRange[sDim][BinId[sDim]+1]);
          }
          delete [] BinId;
          return;
        }
        else{
          Fraction*=2;
          Step=TopBin/Fraction;
          if(Step==0)Step=1;
        }

        if(value_up<rannum) BinCandidate+=Step;
        else BinCandidate-=Step;
        
        //we have to return the last (top) bin
        //this may either be the overflow bin, or the last normal bin (as saved in the array)
        if(BinCandidate==TopBin+1){
            if(UnderOverFlow){
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    axisValues[sDim] = std::numeric_limits<double>::quiet_NaN();
                }
                return;                
            }
            else{
                unsigned* BinId = new unsigned [Dim];
                GetBinCoordinates(TopBin_Glob,BinId);
                for(unsigned short sDim=0; sDim<Dim; sDim++){
                    //axisValues[sDim] = BinCenter[sDim][BinId[sDim]];
                    axisValues[sDim] = RanGen->Uniform(BinRange[sDim][BinId[sDim]],BinRange[sDim][BinId[sDim]+1]);
                }
                delete [] BinId;
                return;
            }
        }
        //we have to return the first normal bin (as saved in the array)
        if(BinCandidate==-1){
          unsigned* BinId = new unsigned [Dim];
          uint64_t BinCandidate_Glob = ismid(0);
          GetBinCoordinates(BinCandidate_Glob,BinId);
          for(unsigned short sDim=0; sDim<Dim; sDim++){
            //axisValues[sDim] = BinCenter[sDim][BinId[sDim]];
            axisValues[sDim] = RanGen->Uniform(BinRange[sDim][BinId[sDim]],BinRange[sDim][BinId[sDim]+1]);
          }
          delete [] BinId;
          return;
        }
      }
    }
    double Sample(const bool& UnderOverFlow, const bool& yield) const{
      if(Dim!=1) return 0;
      double result;
      Sample(&result,UnderOverFlow,yield,NULL);
      return result;
    }




    //sparse matrix id
    //the return value is the correponding position in the global array to find the 
    //values associated with WhichTotBin
    //N.B. this will return -1 (18446744073709551615) if it fails to find a bin
    unsigned smid(const uint64_t& WhichTotBin) const{
        if(WhichTotBin>=TotNumBins+2) return -1;
        if(TotNumBins>MaxNumBins){
            //underflow
            if(WhichTotBin==TotNumBins){
                return MaxNumBins;
            }
            //overflow
            else if(WhichTotBin==TotNumBins+1){
                return MaxNumBins+1;
            }
            uint64_t bi = WhichTotBin/NumElPerSparseBucket;
            auto it = BucketId.find(bi);
            if (it != BucketId.end()) {
                //j-th saved bucket
                unsigned bj = it->second;
                return bj*NumElPerSparseBucket + WhichTotBin%NumElPerSparseBucket;
            }
            else{
                return -1;
            }
        }
        else{
            return WhichTotBin;
        }
    }

    unsigned smid(const uint64_t& WhichTotBin, bool create_bucket){
        unsigned smid_val = smid(WhichTotBin);
        //over or underflow
        if(smid_val==MaxNumBins || WhichTotBin==TotNumBins+1){
            return smid_val;
        }
        //sparse
        if(TotNumBins>MaxNumBins){
            if(!create_bucket) return smid_val;
            else if(smid_val!=-1) return smid_val;
            else{
                uint64_t bi = WhichTotBin/NumElPerSparseBucket;
                unsigned bj = BucketId.size();
                smid_val = bj*NumElPerSparseBucket + WhichTotBin%NumElPerSparseBucket;
                if(smid_val>=MaxNumBins){
                    bool static err_sparse_mem = false;
                    if(!err_sparse_mem) printf("\033[1;31mERROR:\033[0m The sparse DLM_Histo ran out of memory. No bins will be set from now on. Try to change the SetUpSparse() settings.\n"); 
                    err_sparse_mem = true;
                    return -1;
                }
                BucketId.try_emplace(bi,bj);
                BucketOrder.push_back(bi);
                return smid_val;            
            }
        }
        else{
            return smid_val;
        }

/*
printf("smid %lu\n",WhichTotBin);
        if(WhichTotBin>=TotNumBins+2) return TotNumBins+1;
        if(TotNumBins>MaxNumBins){
            //BucketId.emplace(bn, j);
            //global bucket id
            uint64_t bi = WhichTotBin/NumElPerSparseBucket;
printf(" searching for %lu bucket (tot buckets %lu)\n",bi,BucketId.size());
            auto it = BucketId.find(bi);
            if (it != BucketId.end()) {
                //j-th saved bucket
                unsigned bj = it->second;
printf(" existing bucket: %u at position %lu*%lu+%lu = %lu\n", bj, bj, NumElPerSparseBucket, WhichTotBin%NumElPerSparseBucket, bj*NumElPerSparseBucket + WhichTotBin%NumElPerSparseBucket);
//usleep(0.1e3);
                return bj*NumElPerSparseBucket + WhichTotBin%NumElPerSparseBucket;
            }
            else if(create_bucket){
                //printf("\033[1;31mERROR:\033[0m DLM_Histo::smid Congrats, you won a bug, please contact the developer.\n");
//printf("smid NONE\n");
printf(" NO bucket\n");
usleep(100e3);
                //bi,bj (key,value)
                unsigned bj = BucketId.size();
                BucketId.try_emplace(bi,bj);
                return bj*NumElPerSparseBucket + WhichTotBin%NumElPerSparseBucket;
            }
            else{
                return -1;
            }
        }
        else{
            return WhichTotBin;
        }
*/

    }

    //the "inverse" - based on local bin id, we get the total (global) bin id
    uint64_t ismid(const unsigned& WhichLocalBin) const{
        if(TotNumBins>MaxNumBins){
            //underflow
            if(WhichLocalBin==MaxNumBins){
                return TotNumBins;
            }
            else if(WhichLocalBin==MaxNumBins+1){
                return TotNumBins+1;
            }
            //the bucket id
            unsigned bj = WhichLocalBin/NumElPerSparseBucket;
            //this is the true id of the bucket
            uint64_t bi; 
            if(bj<BucketOrder.size()){
                bi = BucketOrder.at(bj);
                //this the true id of the element
                return bi*NumElPerSparseBucket + WhichLocalBin%NumElPerSparseBucket;
            }
            else{
                return -1;
            }
        }
        else{
            return WhichLocalBin;
        }
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

//BELOW IS THE WHOLE DATA OF THE HISTO!!
//if we read/write to files, we have two options:
//  1) QuickRead/Write: assumes that Type contains only basic data types,
//      that are NOT dynamically initialized or so. I.e. we can just use
//      sizeof(Type) and cast it to get the full info of this object
//  2) NOT IMPLEMENTED YET
//     Use a read/write funciton within the object to make more complex initializations

//those are always there when the histo is initialized. To be saved in the file
//The QuickWrite ONLY saves these values, nothing related to integraion etc (below)
    //this is info about the file format.
    //current values:
    // 0 : error
    // Negative: QuickWrite, Positive : FullWrite
    // more info in the functions
    int WriteVersion=0;
    unsigned short Dim;
    unsigned* NumBins;
    double** BinRange;
    //the last two bins are under/overflow
    Type* BinValue;
    //used by the sparse matrix
    //[global_bucket_number][bucket_number_in_the_array]
    std::map<uint64_t, unsigned> BucketId;
    //[bucket_number_in_the_array] and the value is the global_bucket_number
    std::vector<uint64_t> BucketOrder;
    Type* BinError;
    double** BinCenter;
////////////////////////////////////////////////////
    uint64_t TotNumBins;
    uint64_t bid;//will be often used, init here one time only, reuse
    //if TotNumBins>MaxNumBins, we use the sparse mechanism to get around.
    //by default this is 23M
    //if the sparse if off, there is a limit of 4200M
    const unsigned DefaultMaxNumBins;
    unsigned MaxNumBins;
    //unsigned NumDesiredSparseBuckets;
    unsigned NumElPerSparseBucket;


    //used for the extrapolation algorithm.
    //this is all fixed based on the dimensions of the histo
    //NO NEED TO BE SAVED IN THE OUTPUT FILE, but needs to be reinited when reading
    unsigned NumPermutations;
    char** PER;

    Type* CumulativeValue;
    DLM_Random* rangen;

    bool Initialized;
    bool CumUpdated;
    Type INT_ERROR;
    Type NULL_ELEM;
    Type NAN_ELEM;
    Type INF_ELEM;

};
#endif
