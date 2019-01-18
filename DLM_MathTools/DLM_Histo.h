
#ifndef DLM_HISTO_H
#define DLM_HISTO_H

#include <stdio.h>
//#include <stdint.h>
//#include <complex>

template <class Type> class DLM_Histo{
public:
    DLM_Histo(const unsigned& numbin, const Type* bins):NumBins(numbin){
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
    DLM_Histo(const unsigned& numbin, const Type& xmin, const Type& xmax):NumBins(numbin){
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
    DLM_Histo(const DLM_Histo& other):NumBins(other.NumBins){
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
    ~DLM_Histo(){
        if(BinRange) {delete [] BinRange; BinRange=NULL;}
        if(BinValue) {delete [] BinValue; BinValue=NULL;}
        if(BinCenter) {delete [] BinCenter; BinCenter=NULL;}
    }

    void SetBinContent(const unsigned& WhichBin, const double& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]=Val;
    }
    void SetBinCenter(const unsigned& WhichBin, const double& Val){
        if(WhichBin>=NumBins) return;
        BinCenter[WhichBin]=Val;
    }
    void Add(const unsigned& WhichBin, const double& Val){
        if(WhichBin>=NumBins) return;
        BinValue[WhichBin]+=Val;
    }
    void SetBinAt(const double& xVal, const double& Val){
        unsigned WhichBin = GetBin(xVal);
        SetBinContent(WhichBin, Val);
    }
    void AddAt(const double& xVal, const double& Val){
        unsigned WhichBin = GetBin(xVal);
        Add(WhichBin, Val);
    }

    unsigned GetBin(const double& xVal) const{
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
    Type Eval(const double& xVal) const{
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

    bool Copy(const DLM_Histo& other, const Type& Scale){
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
    bool Add(const DLM_Histo& other, const Type& Scale){
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

    bool operator=(const DLM_Histo& other){
        return Copy(other, 1);
    }
    bool operator+=(const DLM_Histo& other){
        return Add(other, 1);
    }

protected:

    DLM_Histo(const unsigned& numbin):NumBins(numbin){
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

#endif
