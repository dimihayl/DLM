
#ifndef DLM_HISTO_H
#define DLM_HISTO_H

#include <stdio.h>
#include <math.h>
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
        if(WhichBin>=NumBins) return 0;
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
        BinRange = NULL;
        BinRange = NULL;
        BinRange = NULL;
        NumBins = NULL;
        TotNumBins=0;
    }
    DLM_Histo(const DLM_Histo& other){
        operator=(other);
    }
    ~DLM_Histo(){
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
        BinRange = new Type* [Dim];
        BinValue = new Type [Dim];
        BinError = new Type [Dim];
        BinCenter = new Type* [Dim];
        NumBins = new unsigned [Dim];
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            BinRange[sDim]=NULL;
            BinValue[sDim]=0;
            BinError[sDim]=0;
            BinCenter[sDim]=NULL;
            NumBins[sDim]=0;
        }
    }
    void SetUp(const unsigned short& sDim, const unsigned& numbins, const Type* bins){
        if(sDim>=Dim) return;
        if(NumBins[sDim]!=numbins){
            CleanUp(sDim);
        }
        if(!BinRange[sDim]) BinRange[sDim] = new Type[numbins+1];
        if(!BinCenter[sDim]) BinCenter[sDim] = new Type[numbins];
        NumBins[sDim] = numbins;
        for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
            BinRange[sDim][uBin] = bins[uBin];
            if(uBin) BinCenter[sDim][uBin] = (bins[uBin-1]+bins[uBin])*0.5;
        }
    }
    void SetUp(const unsigned short& sDim, const unsigned& numbins, const Type& xmin, const Type& xmax){
        if(sDim>=Dim) return;
        if(NumBins[sDim]!=numbins){
            CleanUp(sDim);
        }
        if(!BinRange[sDim]) BinRange[sDim] = new Type[numbins+1];
        if(!BinCenter[sDim]) BinCenter[sDim] = new Type[numbins];
        NumBins[sDim] = numbins;
        if(NumBins[sDim]==1){
            BinRange[sDim][0] = xmin;
            BinRange[sDim][1] = xmax;
            BinCenter[sDim][0] = (xmin+xmax)*0.5;
        }
        else{
            Type BinWidth = (xmax-xmin)/Type(NumBins[sDim]);
            for(unsigned uBin=0; uBin<=NumBins[sDim]; uBin++){
                BinRange[sDim][uBin] = xmin + Type(uBin)*BinWidth;
                if(uBin) BinCenter[sDim][uBin] = (BinRange[sDim][uBin-1]+BinRange[sDim][uBin])*0.5;
            }
        }
    }
    bool Initialize(const bool& ZeroElements=true){
        if(!Dim||!NumBins||!BinRange||!BinCenter) return false;
        TotNumBins=1;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
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
    }
    bool Copy(const DLM_Histo& other){
        return operator=(other);
    }
    bool Add(const DLM_Histo& other, const Type& scale=1){
        if(!(operator+=(other))) return false;
        Scale(scale);
        return true;
    }
    void Add(const Type& value, const Type& error=0){
        for(unsigned uBin=0; uBin<TotNumBins; uBin++){
            BinValue[uBin]+=value;
            BinError[uBin]=Type(sqrt(double(BinError[uBin]*BinError[uBin])+double(error*error)));
        }
    }
    void Scale(const Type& scale){
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                BinValue[sDim][uBin] *= scale;
                BinError[sDim][uBin] *= scale;
            }
            BinValue[sDim][TotNumBins] *= scale;
            BinValue[sDim][TotNumBins+1] *= scale;
        }
    }
    void ScaleToBinWidth(){
        //Type BinWidth;
//!STILL TO DO!
    }
    unsigned GetNbins() const{
        return TotNumBins;
    }
    unsigned GetNbins(const unsigned short& sDim) const{
        return NumBins?NumBins[sDim]:0;
    }

    void SetBinContent(const unsigned& WhichTotBin, const Type& Val){
        if(WhichTotBin>=TotNumBins) return;
        BinValue[WhichTotBin]=Val;
    }
     void SetBinError(const unsigned& WhichTotBin, const Type& Val){
        if(WhichTotBin>=TotNumBins) return;
        BinError[WhichTotBin]=Val;
    }
    void SetBinCenter(const unsigned& WhichTotBin, const Type& Val){
        if(WhichTotBin>=TotNumBins) return;
        BinCenter[WhichTotBin]=Val;
    }
    void Add(const unsigned& WhichTotBin, const Type& Val, const Type& Err=0){
        if(WhichTotBin>=TotNumBins) return;
        BinValue[WhichTotBin]+=Val;
        BinError[WhichTotBin]=Type(sqrt(double(BinError[WhichTotBin]*BinError[WhichTotBin])+double(Err*Err)));
    }

    void SetBinContent(const unsigned short& sDim, const unsigned& WhichBin, const Type& Val){
        SetBinContent(GetTotBin(sDim,WhichBin),Val);
    }
    void SetBinError(const unsigned short& sDim, const unsigned& WhichBin, const Type& Val){
        SetBinError(GetTotBin(sDim,WhichBin),Val);
    }
    void SetBinCenter(const unsigned short& sDim, const unsigned& WhichBin, const Type& Val){
        SetBinCenter(GetTotBin(sDim,WhichBin),Val);
    }
    void Add(const unsigned short& sDim, const unsigned& WhichBin, const Type& Val, const Type& Err=0){
        Add(GetTotBin(sDim,WhichBin),Val,Err);
    }

    void SetBinAt(const unsigned short& sDim, const Type& xVal, const Type& Val){
        SetBinContent(sDim,GetBin(sDim,xVal),Val);
    }
    void AddAt(const unsigned short& sDim, const Type& xVal, const Type& Val, const Type& Err=0){
        Add(sDim,GetBin(sDim,xVal),Val,Err);
    }

    Type GetBinCenter(const unsigned& WhichTotBin) const{
        if(WhichTotBin>=TotNumBins) return 0;
        return BinCenter[WhichTotBin];
    }
    Type GetBinCenter(const unsigned short& sDim, const unsigned& WhichBin) const{
        return GetBinCenter(sDim,GetTotBin(sDim,WhichBin));
    }

    Type GetBinContent(const unsigned& WhichTotBin) const{
        if(WhichTotBin>=TotNumBins+2) return 0;
        return BinValue[WhichTotBin];
    }
    Type GetBinContent(const unsigned short& sDim, const unsigned& WhichBin) const{
        return GetBinContent(sDim,GetBin(sDim,WhichBin));
    }

    Type GetBinError(const unsigned& WhichTotBin) const{
        if(WhichTotBin>=TotNumBins) return 0;
        return BinValue[WhichTotBin];
    }
    Type GetBinError(const unsigned short& sDim, const unsigned& WhichBin) const{
        return GetBinError(sDim,GetBin(sDim,WhichBin));
    }


    Type GetBinLowEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
        return BinRange[WhichBin];
    }
    Type GetBinUpEdge(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
        return BinRange[WhichBin+1];
    }
    Type GetBinWidth(const unsigned short& sDim, const unsigned& WhichBin) const{
        if(sDim>=Dim) return 0;
        if(WhichBin>=NumBins[sDim]) return 0;
        return (BinRange[sDim][WhichBin+1]-BinRange[sDim][WhichBin]);
    }

    unsigned GetBin(const unsigned short& sDim, const Type& xVal) const{
        if(sDim>=Dim) return TotNumBins;
        unsigned WhichTotBin=GetTotBin(sDim,0);
        if(NumBins[sDim]==0) return TotNumBins;
        if(NumBins[sDim]==1) return WhichTotBin;
        unsigned WhichBin=(NumBins[sDim]+1)/2;
        unsigned BinMod=4;
        unsigned BinStep;
        //makes sure that the xVal is in BinRange. If not, the returned value is either
        //NumBins or NumBins+1, depending on if we have an underflow or overflow
        if(xVal<BinRange[sDim][0]) return WhichTotBin+NumBins[sDim];
        if(xVal>BinRange[sDim][NumBins[sDim]]) return WhichTotBin+NumBins[sDim]+1;
        while(true){
            if(BinRange[sDim][WhichBin]<=xVal && BinRange[sDim][WhichBin+1]>=xVal){
                return WhichTotBin+WhichBin;
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
*/
    bool operator=(const DLM_Histo& other){
        SetUp(other.Dim);
        for(unsigned short sDim=0; sDim<Dim; sDim++) SetUp(sDim,other.NumBins[sDim],BinRange[sDim]);
        if(!Initialize(false)) return false;
        for(unsigned short sDim=0; sDim<Dim; sDim++){
            for(unsigned uBin=0; uBin<NumBins[sDim]; uBin++){
                BinRange[sDim][uBin] = other.BinRange[sDim][uBin];
                BinCenter[sDim][uBin] = other.BinCenter[sDim][uBin];
            }
            BinRange[sDim][NumBins[sDim]] = other.BinRange[sDim][NumBins[sDim]];

            for(unsigned uBin=0; uBin<TotNumBins; uBin++){
                BinValue[sDim][uBin] = other.BinValue[sDim][uBin];
                BinError[sDim][uBin] = other.BinError[sDim][uBin];
            }
            BinValue[sDim][TotNumBins] = other.BinValue[sDim][TotNumBins];
            BinValue[sDim][TotNumBins+1] = other.BinValue[sDim][TotNumBins+1];
        }
    }
    bool operator+=(const DLM_Histo& other){
        if(!SameStructure(other)) return false;
        for(unsigned uBin=0; uBin<TotNumBins; uBin++) BinValue[uBin] += other.BinValue[uBin];
        return true;
    }

protected:
    void CleanUp(){
        if(BinRange){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                delete [] BinRange[sDim]; BinRange[sDim]=NULL;
            }
            delete [] BinRange; BinRange=NULL;
        }
        if(BinValue){
            delete [] BinValue; BinValue=NULL;
        }
        if(BinError){
            delete [] BinError; BinError=NULL;
        }
        if(BinCenter){
            for(unsigned short sDim=0; sDim<Dim; sDim++){
                delete [] BinCenter[sDim]; BinCenter[sDim]=NULL;
            }
            delete [] BinCenter; BinCenter=NULL;
        }
        delete [] NumBins; NumBins=NULL;
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
    bool SameStructure(DLM_Histo& other) const{
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
    unsigned GetTotBin(const unsigned short& sDim, const unsigned& WhichBin) const{
        unsigned WhichTotBin=WhichBin;
        for(unsigned short sDim2=0; sDim2<sDim; sDim2++){
            WhichTotBin += NumBins[sDim2];
        }
        return WhichTotBin;
    }
    unsigned short Dim;
    unsigned TotNumBins;
    unsigned* NumBins;
    Type** BinRange;
    //the last two bins are under/overflow
    Type* BinValue;
    Type* BinError;
    Type** BinCenter;
};
#endif
