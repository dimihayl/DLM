#ifndef DLM_SortH
#define DLM_SortH

#include "stdio.h"
#include <stdlib.h>
#include "math.h"

#include "DLM_CppTools.h"
#include "DLM_Histo.h"

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
class DLM_Sort{
private:
	//Key[i] says where in the input file is the i-th highest/lowest element located
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

    void MergeSortBox(const Num& Start, const Num& Mid, const Num& End, const bool& descending){
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
            else if(Input[Key[PosFirstBox]] > Input[Key[PosSecondBox]] && !descending){
                KeyTemp[i] = Key[PosSecondBox];
                PosSecondBox++;
            }
            else if(Input[Key[PosFirstBox]] < Input[Key[PosSecondBox]] && descending){
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
            else if(Input[Key[PosSecondBox]] < Input[Key[PosFirstBox]] && !descending){
                KeyTemp[i] = Key[PosFirstBox];
                PosFirstBox--;
            }
            else if(Input[Key[PosSecondBox]] > Input[Key[PosFirstBox]] && descending){
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
    DLM_Sort(){
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
            printf("DLM_Sort says: ERROR! You cannot have float/double as a typename for this class!");
            return;
        }
        NumBoxesPairs = new int [4000000];
    }
    ~DLM_Sort(){
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
    //copies the key
    void SetKey(Num* InputKey){
      if(Key){delete [] Key;}
      Key = new Num [NumOfEl];
      for(Num i=0; i<NumOfEl; i++){
          Key[i] = InputKey[i];
      }
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
    bool MergeSort(const bool& descending=false){
        if(!Input){
            printf("DLM_Sort says: ERROR! No data is loaded for sorting! Please use SetData(Element* input, Num N).");
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
                MergeSortBox(Start, Mid, End, descending);
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

    bool BucketSort(const bool& descending=false){
        if(!Input){
            printf("DLM_Sort says: ERROR! No data is loaded for sorting! Please use SetData(Element* input, Num N).");
            return false;

            Num NumOfBuckets = round(pow(2.*double(NumOfEl)*double(NumOfEl),1./3.));
            Num AvgNumElPerBucket = NumOfEl/NumOfBuckets;
            Num NumSampleElements = round(sqrt(double(NumOfEl)));
            Element* BucketRanges = new Element[NumOfBuckets+1];
            if(!NumOfBuckets) NumOfBuckets=1;
            if(!AvgNumElPerBucket) AvgNumElPerBucket=1;
            Num* MemoryPerBucket = new Num [NumOfBuckets];
            for(Num nBuck=0; nBuck<NumOfBuckets; nBuck++){
				MemoryPerBucket[nBuck] = AvgNumElPerBucket+AvgNumElPerBucket/10+1;
			}

            Num RealNumEl=NumOfEl;
            NumOfEl=NumSampleElements;
            //sort only small fraction of the data
            //this gives you an idea how the pdf of the data looks like, which is than used
            //to determine a meaningful choice of the bucket ranges, such that all bins have
            //similar number of elements in them, which leads to the best performance
            MergeSort(descending);
            NumOfEl=RealNumEl;
            BucketRanges[0] = Output[0];
            Num CurrentBucket=0;
            Num CurrentBucketElemets=1;
            //iterate over the small fraction of data, counting the number of elements
            //that are put in each bin. The moment we have more than the avg number of elements,
            //the bin is considered full and the corresponding value of the element saved
            for(Num nEl=0; nEl<NumSampleElements; nEl++){
                if(CurrentBucketElemets>AvgNumElPerBucket){
                    BucketRanges[CurrentBucket+1]=Output[nEl];
                    CurrentBucket++;
                    CurrentBucketElemets=0;
                }
                else{
                    CurrentBucketElemets++;
                }
            }
            BucketRanges[NumOfBuckets] = Output[NumSampleElements-1];
            Element** Bucket = new Element[NumOfBuckets];
            Num* NumElInBucket = new Num[NumOfBuckets];
            for(Num nBuck=0; nBuck<NumOfBuckets; nBuck++){
                Bucket[nBuck] = new Element [MemoryPerBucket[nBuck]];
                NumElInBucket[nBuck] = 0;
            }
            DLM_Histo1D<Element> Histo(NumOfBuckets,BucketRanges);
            //assign each element to a bucket
            for(Num nEl=0; nEl<NumOfEl; nEl++){
				unsigned WhichBucket = Histo.GetBin(Input[nEl]);
				if(NumElInBucket[WhichBucket]>=MemoryPerBucket[WhichBucket]){
					Num OldMem = MemoryPerBucket[WhichBucket];
					MemoryPerBucket[WhichBucket] *= 2;
					Element* Temp = new Element [MemoryPerBucket[WhichBucket]];
					for(Num nEl=0; nEl<OldMem; nEl++){
						Temp[nEl] = Bucket[WhichBucket];
					}
					delete [] Bucket[WhichBucket];
					Bucket[WhichBucket] = Temp[nEl];
				}
				Bucket[WhichBucket][NumElInBucket[WhichBucket]] = Input[nEl];
				NumElInBucket[WhichBucket]++;
            }
//the problem is the fucking key, and the fact, that Key[0] should give me the Id if the smallest element
//but how do I keep track of the Id of each element?
            //sort each bucket
            Num** BucketKey = new Num*[NumOfBuckets];
            #pragma omp parallel for
            for(Num nBuck=0; nBuck<NumOfBuckets; nBuck++){
				BucketKey[nBuck] = new Num [NumElInBucket[nBuck]];
				DLM_Sort<Element,Num> BucketSorter;
				BucketSorter.SetData(Bucket[nBuck],NumElInBucket[nBuck]);
				BucketSorter.MergeSort();
				BucketSorter.GetSortedData(Bucket[nBuck],Bucket[nBuck]);
				BucketSorter.CopyKey(BucketKey);
			}

			//combine all buckets
			Num CurrentElement=0;
			for(Num nBuck=0; nBuck<NumOfBuckets; nBuck++){
				for(Num nEl=0; nEl<NumElInBucket[nBuck]; nEl++){
					Key[CurrentElement] = BucketKey[nBuck]+CurrentElement;
					CurrentElement++;
				}
			}

            for(Num nBuck=0; nBuck<NumOfBuckets; nBuck++){
                delete [] Bucket[nBuck];
                delete [] BucketKey[nBuck];
            }
            delete [] Bucket;
            delete [] NumElInBucket;
            delete [] BucketRanges;
            delete [] MemoryPerBucket;
            delete [] BucketKey;
        }
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
