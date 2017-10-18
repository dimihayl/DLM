
#include <cstdlib>

#include "DLM_SubPads.h"
#include "TROOT.h"


DLM_SubPads::DLM_SubPads(const UInt_t PW, const UInt_t PH)
    :CanvasPadWidth(PW),CanvasPadHeight(PH),MaxSubPadsStep(64){
    if(!MaxSubPadsStep){
        Printf("FATAL ERROR in DLM_SubPads()! The constant MaxSubPadsStep should be >0 !!! ABORTING!");
        abort();
    }
    InitialValues();
    Initialize();
}

DLM_SubPads::~DLM_SubPads(){
    Destroy();
}

void DLM_SubPads::InitialValues(){
    NumSubPads = 0;
    SubPadPostion = NULL;
    ScaleFactor = NULL;
    SubPad = NULL;
    MainCanvas = NULL;
}

void DLM_SubPads::Initialize(){
    MainCanvas = new TCanvas(TString::Format("MainCanvas%p", this));
    MainCanvas->cd(0);
    MainCanvas->SetCanvasSize(CanvasPadWidth, CanvasPadHeight);
    //MainPad = new TPad("MainPad", "MainPad", 0,0,1,1);//lbrt
    //MainPad->Draw();
    MaxSubPads = MaxSubPadsStep;
    SubPad = new TPad* [MaxSubPadsStep];
    SubPadPostion = new Float_t* [MaxSubPadsStep];
    ScaleFactor = new Float_t* [MaxSubPadsStep];
}


//1) Backup the old arrays
//2) Delete the old arrays
//3) Reallocate the relevant arrays
//4) Restore the old values in the new arrays
//5) Clean up the backup
void DLM_SubPads::AllocMoreMemForSubPad(){
//1) Backup the old arrays------------------------------------------------------
    UShort_t tempSize = MaxSubPads;
    TPad** tempPad = new TPad* [tempSize];
    Float_t** tempMargin = new Float_t* [tempSize];
    Float_t** tempScaleFactor = new Float_t* [tempSize];

    for(UShort_t iPad=0; iPad<tempSize; iPad++){
        tempMargin[iPad] = new Float_t [4];
        tempScaleFactor[iPad] = new Float_t [2];

        tempPad[iPad] = SubPad[iPad];

        tempMargin[iPad][0] = SubPadPostion[iPad][0];
        tempMargin[iPad][1] = SubPadPostion[iPad][1];
        tempMargin[iPad][2] = SubPadPostion[iPad][2];
        tempMargin[iPad][3] = SubPadPostion[iPad][3];

        tempScaleFactor[iPad][0] = ScaleFactor[iPad][0];
        tempScaleFactor[iPad][1] = ScaleFactor[iPad][1];
    }
//2) Delete the old arrays------------------------------------------------------
    delete [] SubPad;
    for(UShort_t iPad=0; iPad<tempSize; iPad++){
        delete [] SubPadPostion[iPad];
        delete [] ScaleFactor[iPad];
    }
    delete [] SubPadPostion;
    delete [] ScaleFactor;
//3) Reallocate the relevant arrays------------------------------------------------------
    MaxSubPads += MaxSubPadsStep;
    SubPad = new TPad* [MaxSubPads];
    SubPadPostion = new Float_t* [MaxSubPads];
    ScaleFactor = new Float_t* [MaxSubPads];

//4) Restore the old values in the new arrays------------------------------------------------------
    for(UShort_t iPad=0; iPad<tempSize; iPad++){
        SubPadPostion[iPad] = new Float_t [4];
        ScaleFactor[iPad] = new Float_t [2];
        SubPad[iPad] = tempPad[iPad];

        SubPadPostion[iPad][0] = tempMargin[iPad][0];
        SubPadPostion[iPad][1] = tempMargin[iPad][1];
        SubPadPostion[iPad][2] = tempMargin[iPad][2];
        SubPadPostion[iPad][3] = tempMargin[iPad][3];

        ScaleFactor[iPad][0] = tempScaleFactor[iPad][0];
        ScaleFactor[iPad][1] = tempScaleFactor[iPad][1];
    }

//5) Clean up the backup------------------------------------------------------
    for(UShort_t iPad=0; iPad<tempSize; iPad++){
        delete [] tempMargin[iPad];
        delete [] tempScaleFactor[iPad];
    }
    delete [] tempMargin;
    delete [] tempScaleFactor;

    delete [] tempPad;
}

void DLM_SubPads::Destroy(){
    if(SubPad){
        for(UShort_t iPad=0; iPad<NumSubPads; iPad++){
            delete SubPad[iPad];
            SubPad[iPad] = NULL;
        }
        delete [] SubPad;
        SubPad = NULL;
    }
    if(SubPadPostion){
        for(UShort_t iPad=0; iPad<NumSubPads; iPad++){
            delete [] SubPadPostion[iPad];
            SubPadPostion[iPad] = NULL;
        }
        delete [] SubPadPostion;
        SubPadPostion = NULL;
    }
    if(ScaleFactor){
        for(UShort_t iPad=0; iPad<NumSubPads; iPad++){
            delete [] ScaleFactor[iPad];
            ScaleFactor[iPad] = NULL;
        }
        delete [] ScaleFactor;
        ScaleFactor = NULL;
    }
    InitialValues();
    if(MainCanvas){
        delete MainCanvas;
        MainCanvas = NULL;
    }
}

Float_t DLM_SubPads::ConvertYPadPos(Float_t ypos){
    return 1 - ypos;
}

Bool_t DLM_SubPads::cd(UShort_t WhichPad){
    if(WhichPad<NumSubPads) SubPad[WhichPad]->cd();
    else {Printf("WARNING! Bad input in DLM_SubPads::cd"); return kFALSE;}
    return kTRUE;
}

void DLM_SubPads::SetLogx(UShort_t WhichPad, Bool_t lg){
    if(WhichPad>=NumSubPads) {Printf("WARNING! Bad input in DLM_SubPads::SetLogx"); return;}
    SubPad[WhichPad]->SetLogx(lg);
}

void DLM_SubPads::SetLogy(UShort_t WhichPad, Bool_t lg){
    if(WhichPad>=NumSubPads) {Printf("WARNING! Bad input in DLM_SubPads::SetLogy"); return;}
    SubPad[WhichPad]->SetLogy(lg);
}

void DLM_SubPads::SetGridx(UShort_t WhichPad, Bool_t val){
    if(WhichPad>=NumSubPads) {Printf("WARNING! Bad input in DLM_SubPads::SetGridx"); return;}
    SubPad[WhichPad]->SetGridx(val);
}

void DLM_SubPads::SetGridy(UShort_t WhichPad, Bool_t val){
    if(WhichPad>=NumSubPads) {Printf("WARNING! Bad input in DLM_SubPads::SetGridy"); return;}
    SubPad[WhichPad]->SetGridy(val);
}

void DLM_SubPads::SetGrid(UShort_t WhichPad, Bool_t val){
    if(WhichPad>=NumSubPads) {Printf("WARNING! Bad input in DLM_SubPads::SetGrid"); return;}
    SubPad[WhichPad]->SetGrid(val,val);
}

Bool_t DLM_SubPads::AddSubPad(Float_t Left, Float_t Right, Float_t Bottom, Float_t Top){
    if(Left<0 || Left>1) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    if(Right<0 || Right>1) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    if(Bottom<0 || Bottom>1) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    if(Top<0 || Top>1) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    if(Left>=Right) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    if(Bottom>=Top) {Printf("WARNING! Bad input in DLM_SubPads::AddSubPad"); return kFALSE;}
    for(UShort_t iPad=0; iPad<NumSubPads; iPad++){
        //check if there are overlapping pads. If so => return false
        if (Left < SubPadPostion[iPad][kMarRight] && Right > SubPadPostion[iPad][kMarLeft] &&
            Bottom < SubPadPostion[iPad][kMarTop] && Top > SubPadPostion[iPad][kMarBottom] ) return kFALSE;
    }

    if(NumSubPads>=MaxSubPads) AllocMoreMemForSubPad();
    SubPad[NumSubPads] = new TPad(TString::Format("SubPad%i", NumSubPads),
                                  TString::Format("SubPad%i", NumSubPads),
                                  Left, Top, Right, Bottom);//ltrb

    SubPadPostion[NumSubPads] = new Float_t [4];
    ScaleFactor[NumSubPads] = new Float_t [2];

    SubPadPostion[NumSubPads][kMarLeft] = Left;
    SubPadPostion[NumSubPads][kMarRight] = Right;
    SubPadPostion[NumSubPads][kMarTop] = Top;
    SubPadPostion[NumSubPads][kMarBottom] = Bottom;

    ScaleFactor[NumSubPads][kWidth] = 1./(Right - Left);
    ScaleFactor[NumSubPads][kHeight] = 1./(Top - Bottom);

    MainCanvas->cd(0);
    SubPad[NumSubPads]->Draw();
    SubPad[NumSubPads]->cd();

    NumSubPads++;
    return kTRUE;
}

Bool_t DLM_SubPads::AddSubPadTL(Float_t Left, Float_t Right, Float_t Top, Float_t Bottom){
    return AddSubPad(Left, Right, 1 - Bottom, 1 - Top);
}

void DLM_SubPads::Reset(){
    Destroy();
    InitialValues();
    Initialize();
}

Bool_t DLM_SubPads::SetMargin(UShort_t PadNr, Float_t Left, Float_t Right, Float_t Bottom, Float_t Top,
                              Bool_t MainCoordinates){
    if(PadNr>=NumSubPads) {Printf("ERROR in DLM_SubPads::SetMargin! Accessing non-existing TPad!"); return kFALSE;}
    if(MainCoordinates){
        Left = ScaleFactor[PadNr][kWidth]*Left;
        Right = ScaleFactor[PadNr][kWidth]*Right;
        Bottom = ScaleFactor[PadNr][kHeight]*Bottom;
        Top = ScaleFactor[PadNr][kHeight]*Top;
    }
    if(Left<0 || Left>1) {Printf("WARNING! Bad input in DLM_SubPads::SetMargin"); return kFALSE;}
    if(Right<0 || Right>1) {Printf("WARNING! Bad input in DLM_SubPads::SetMargin"); return kFALSE;}
    if(Bottom<0 || Bottom>1) {Printf("WARNING! Bad input in DLM_SubPads::SetMargin"); return kFALSE;}
    if(Top<0 || Top>1) {Printf("WARNING! Bad input in DLM_SubPads::SetMargin"); return kFALSE;}

    SubPad[PadNr]->SetMargin(Left, Right, Bottom, Top);

    return kTRUE;
}

Float_t DLM_SubPads::GetScaleWidth(const UShort_t& PadNr){
    return PadNr<NumSubPads?ScaleFactor[PadNr][0]:0;
}
Float_t DLM_SubPads::GetScaleHeight(const UShort_t& PadNr){
    return PadNr<NumSubPads?ScaleFactor[PadNr][1]:0;
}

TCanvas* DLM_SubPads::GetCanvas(){
    return MainCanvas;
}

//! THE MOST IDIOTIC THING I HAVE FOUND ABOUT ROOT:
//So... if the width of the Pad (in pixels) is smaller then the height, then the size of the labels
//scales with the width. BUT if it is the other way around... surprise surprise, the size scales
//with the height.
void DLM_SubPads::SetLabelSize(const UShort_t& PadNr, TAxis* Axis, const Float_t& Size){

    //those coefficients can be used for a direct comparison between width and height,
    //since they take into account the difference in relative size compared to the TCanvas and
    //at the same time the difference in the aspect ratio (Width vs Height).
    double RealWidth = CanvasPadWidth/ScaleFactor[PadNr][kWidth];
    double RealHeight = CanvasPadHeight/ScaleFactor[PadNr][kHeight];
    //this factor is chosen such that the Label size has typical values between 1 and 100
    double OrderOfMag = 0.001*(CanvasPadWidth+CanvasPadHeight);

    if(RealHeight>=RealWidth)
        Axis->SetLabelSize(Size/RealWidth*OrderOfMag);
    else{
        Axis->SetLabelSize(Size/RealHeight*OrderOfMag);
    }
}

void DLM_SubPads::SetTitleSize(const UShort_t& PadNr, TAxis* Axis, const Float_t& Size){

    //those coefficients can be used for a direct comparison between width and height,
    //since they take into account the difference in relative size compared to the TCanvas and
    //at the same time the difference in the aspect ratio (Width vs Height).
    double RealWidth = CanvasPadWidth/ScaleFactor[PadNr][kWidth];
    double RealHeight = CanvasPadHeight/ScaleFactor[PadNr][kHeight];
    //this factor is chosen such that the Label size has typical values between 1 and 100
    double OrderOfMag = 0.001*(CanvasPadWidth+CanvasPadHeight);

    if(RealHeight>=RealWidth)
        Axis->SetTitleSize(Size/RealWidth*OrderOfMag);
    else{
        Axis->SetTitleSize(Size/RealHeight*OrderOfMag);
    }
}
