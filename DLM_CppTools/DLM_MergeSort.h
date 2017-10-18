#ifndef DLM_MergeSortH
#define DLM_MergeSortH

#include "stdio.h"
#include <stdlib.h>

#include "DLM_CppTools.h"

template <class Element, typename Num>
void MergeBoxes(Element* b1, Num& n1, Element* b2, Num& n2, Element* Result){
    Num Ntot = n1 + n2;
    Num pragmasplit = Ntot/2;

    Num Pos1 = 0;
    Num Pos2 = 0;
    for(Num i=0; i<pragmasplit; i++){
        //we have no elements left to compare
        if(Pos1==n1){
            Result[i] = b2[Pos2];
            Pos2++;
        }
        //we have no elements left to compare
        else if(Pos2==n2){
            Result[i] = b1[Pos1];
            Pos1++;
        }
        else if(b1[Pos1] > b2[Pos2]){
            Result[i] = b2[Pos2];
            Pos2++;
        }
        else{
            Result[i] = b1[Pos1];
            Pos1++;
        }
    }

    Pos1 = n1-1;
    Pos2 = n2-1;
    for(Num i=Ntot-1; i>=pragmasplit; i--){
        //we have no elements left to compare
        if(Pos1+1==0){
            Result[i] = b2[Pos2];
            Pos2--;
        }
        //we have no elements left to compare
        else if(Pos2+1==0){
            Result[i] = b1[Pos1];
            Pos1--;
        }
        else if(b1[Pos1] < b2[Pos2]){
            Result[i] = b2[Pos2];
            Pos2--;
        }
        else{
            Result[i] = b1[Pos1];
            Pos1--;
        }
    }
}

template <class Element, typename Num>
class DLM_MergeSort{
private:
    Num* Key;
    Num* KeyTemp;
    Num NumOfEl;
    const Element* Input;
    Element* Output;

    bool InitMem(){
        ClearMem();
        if(!Input){
            return false;
        }
        if(NumOfEl<=0){
            Input = NULL;
            return false;
        }
        Key = new Num [NumOfEl];
        for(Num i=0; i<NumOfEl; i++){
            Key[i] = i;
        }
        KeyTemp = new Num [NumOfEl];
        return true;
    }

   //imagine two boxes: [a,b,c,d] [e,f,g]
   //a = Start; e = Mid; g = End, where a,b,c... are the array number of KeyTemp

    void MergeSortBox(const Num& Start, const Num& Mid, const Num& End){
        Num pragmasplit = Start + (End-Start)/2;

        Num PosFirstBox = Start;
        Num PosSecondBox = Mid;
        for(Num i=Start; i<=pragmasplit; i++){
            //we have no elements left to compair
            if(PosFirstBox==Mid){
                KeyTemp[i] = Key[PosSecondBox];
                PosSecondBox++;
            }
            //we have no elements left to compair
            else if(PosSecondBox==End+1){
                KeyTemp[i] = Key[PosFirstBox];
                PosFirstBox++;
            }
            else if(Input[Key[PosFirstBox]] > Input[Key[PosSecondBox]]){
                KeyTemp[i] = Key[PosSecondBox];
                PosSecondBox++;
            }
            else{
                KeyTemp[i] = Key[PosFirstBox];
                PosFirstBox++;
            }
        }

        PosFirstBox = Mid-1;
        PosSecondBox = End;
        for(Num i=End; i>pragmasplit && i<=End; i--){
            //we have no elements left to compair
            if(PosSecondBox==Mid-1){
                KeyTemp[i] = Key[PosFirstBox];
                PosFirstBox--;
            }
            //we have no elements left to compair
            else if(PosFirstBox==Start-1){
                KeyTemp[i] = Key[PosSecondBox];
                PosSecondBox--;
            }
            else if(Input[Key[PosSecondBox]] < Input[Key[PosFirstBox]]){
                KeyTemp[i] = Key[PosFirstBox];
                PosFirstBox--;
            }
            else{
                KeyTemp[i] = Key[PosSecondBox];
                PosSecondBox--;
            }
        }

        for(Num i=Start; i<=End; i++){
            Key[i] = KeyTemp[i];
        }

    }

public:
    long ** TimeInfo1;
    long * TimeInfo2;
    int * NumBoxesPairs;
    int NumLevels;
    DLM_MergeSort(){
        //we want Num to be an integer index (or char, short, unsigned etc.)
        Num inttest = 0.4;
        Key = NULL;
        KeyTemp = NULL;
        NumOfEl = 0;
        Input = NULL;
        Output = NULL;
        TimeInfo1 = new long * [30];
        for(long i=0; i<30; i++){
            TimeInfo1[i] = new long [4000000];
        }
        TimeInfo2 = new long [30];
        if(inttest){
            printf("DLM_MergeSort says: ERROR! You cannot have float/double as a typename for this class!");
            return;
        }
        NumBoxesPairs = new int [4000000];
    }
    ~DLM_MergeSort(){
        ClearMem();
        for(long i=0; i<30; i++){
            delete [] TimeInfo1[i];
        }
        delete [] TimeInfo2;
        delete [] TimeInfo1;
        delete [] NumBoxesPairs;
    }

    void SetData(const Element* input, const Num& N){
        Input = input;
        NumOfEl = N;
        InitMem();
    }
    Num* GetKey(){
        return Key;
    }
    void GetSortedData(Element* input, Element* output){
        if(!Key) return;
        if(!output) return;
        Element* Temp;
        if(input==output){
            Temp = new Element[NumOfEl];
        }
        else{
            Temp = output;
        }

        for(Num i=0; i<NumOfEl; i++){
            Temp[i] = input[Key[i]];
        }
        if(input==output){
            for(Num i=0; i<NumOfEl; i++){
                output[i] = Temp[i];
            }
            delete [] Temp;
        }
    }
    Num GetNumOfEl(){
        return NumOfEl;
    }

    void CopyKey(Num* key){
        for(Num i=0; i<NumOfEl; i++){
            key[i] = Key[i];
        }
    }
    void CopySortedData(Element* output){
        if(!Key) {output=NULL; return;}
        for(Num i=0; i<NumOfEl; i++){
            output[i] = Input[Key[i]];
        }
    }

    //The merge-sort works on different "levels". On the zeroth level, we divide the n elements into n boxes and
    //merging and sorting each 2 neighbouring boxes. Thus on the first level we have 2 elements in all boxes. Do note,
    //that in the last box we might have less (i.e. one for the first level) elements. In any case, we again merge-sort
    //each 2 neighbouring boxes and continue until we have all elements in a single box.
    bool MergeSort(){
        if(!Input){
            printf("DLM_MergeSort says: ERROR! No data is loaded for sorting! Please use SetData(Element* input, Num N).");
            return false;
        }

        //the numver of levels is approx. log_2(NumOfEl) + 1, which means that even for 1e30 elements
        //we will have around 100 levels only
        //unsigned char level = 0;
        //# of elements in all the boxes before the last
        Num ElementsPerBox = 1;
        Num NumOfBoxes = NumOfEl;

        int level = 0;
        //iteration over all levels
        while(NumOfBoxes!=1){
            //merge-sorting of the neighboring boxes
            NumBoxesPairs[level] = (NumOfBoxes-1)/2 + (NumOfBoxes-1)%2;
            for(Num i=0; i<NumOfBoxes-1; i=i+2){
                DLM_Timer tim2;
                Num Start = i*ElementsPerBox;
                Num Mid = Start + ElementsPerBox;
                Num End = Mid + ElementsPerBox - 1;
                if(End>=NumOfEl) End = NumOfEl-1;
                MergeSortBox(Start, Mid, End);
            }
            ElementsPerBox *= 2;
            lldiv_t divresult;
            divresult = lldiv (NumOfEl,ElementsPerBox);
            NumOfBoxes = divresult.quot + !!divresult.rem;
            level++;
        }
        NumLevels = level;

        return true;
    }

    void ClearMem(){
        if(Key){
            delete [] Key;
            Key = NULL;
        }
        if(KeyTemp){
            delete [] KeyTemp;
            KeyTemp = NULL;
        }
        if(Output){
            delete [] Output;
            Output = NULL;
        }
    }

};

#endif



