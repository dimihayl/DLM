#ifndef DLM_HistoDrawH
#define DLM_HistoDrawH

#include "TH1.h"
#include "TCanvas.h"

class DLM_HistoDraw{

private:

    Int_t NumberOfHistosX;
    Int_t NumberOfHistosY;
    TH1* ListOfHistos;

    TCanvas* MainCanvas;

    Bool_t* HideHistoXaxis;
    Bool_t* HideHistoYaxis;
    Bool_t* HideHistoZaxis;

    Bool_t SameSizeLabels;

    //Those are in relative coordinates with respect to
    //the main canvas
    Float_t* PadPosition;
    Float_t* PadWidth;
    Float_t* PadHeight;

    //[0] are in relative coordinates with respect to
    //the main canvas
    //[1] are in relative coordinates with respect to
    //the corresponding Pad
    Float_t** HistoPosition;
    Float_t** HistoWidth;
    Float_t** HistoHeight;

    //the size of the labels relative to the main canvas
    //(used only if SameSizeLabels==true)
    Float_t GlobalLabelSize;


    void Reset();
    void InitialValues();

public:

    DLM_HistoDraw();
    ~DLM_HistoDraw();

    Bool_t SetHisto(Int_t iHist, TH1* histo);
    TH1* GetHisto(Int_t iHist);

    Bool_t SetHideHistoXaxis(Int_t iHist, Bool_t hide);
    Bool_t GetHideHistoXaxis(Int_t iHist, Bool_t hide);

    Bool_t SetHideHistoYaxis(Int_t iHist, Bool_t hide);
    Bool_t GetHideHistoYaxis(Int_t iHist, Bool_t hide);

    Bool_t SetHideHistoZaxis(Int_t iHist, Bool_t hide);
    Bool_t GetHideHistoZaxis(Int_t iHist, Bool_t hide);

    TCanvas* GetCanvas();//calc all here, return NULL if some histo is missing?
};


#endif

