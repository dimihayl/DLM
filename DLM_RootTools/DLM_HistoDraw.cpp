#include "DLM_HistoDraw.h"
#include "TROOT.h"


DLM_HistoDraw::DLM_HistoDraw(){
    InitialValues();
}

DLM_HistoDraw::~DLM_HistoDraw(){
    Reset();
}


void DLM_HistoDraw::InitialValues(){
    NumberOfHistosX = 0;
    NumberOfHistosY = 0;
    SameSizeLabels = kTRUE;
    GlobalLabelSize = 0.05;

    ListOfHistos = NULL;

    HideHistoXaxis = NULL;
    HideHistoYaxis = NULL;
    HideHistoZaxis = NULL;

    PadPosition = NULL;
    PadWidth = NULL;
    PadHeight = NULL;

    HistoPosition = NULL;
    HistoWidth = NULL;
    HistoHeight = NULL;
}


void DLM_HistoDraw::Reset(){
    if(ListOfHistos){
        delete [] ListOfHistos;
        ListOfHistos = NULL;
    }
    if(MainCanvas){
        delete MainCanvas;
        MainCanvas = NULL;
    }
    if(HideHistoXaxis){
        delete [] HideHistoXaxis;
        HideHistoXaxis = NULL;
    }
    if(HideHistoYaxis){
        delete [] HideHistoYaxis;
        HideHistoYaxis = NULL;
    }
    if(HideHistoZaxis){
        delete [] HideHistoZaxis;
        HideHistoZaxis = NULL;
    }
    if(PadPosition){
        delete [] PadPosition;
        PadPosition = NULL;
    }
    if(PadWidth){
        delete [] PadWidth;
        PadWidth = NULL;
    }
    if(PadHeight){
        delete [] PadHeight;
        PadHeight = NULL;
    }
    if(HistoPosition){
        delete [] HistoPosition;
        HistoPosition = NULL;
    }
    if(HistoWidth){
        delete [] HistoWidth;
        HistoWidth = NULL;
    }
    if(HistoHeight){
        delete [] HistoHeight;
        HistoHeight = NULL;
    }
}
