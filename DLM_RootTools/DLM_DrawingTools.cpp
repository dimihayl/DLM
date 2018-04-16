#include "DLM_DrawingTools.h"

#include "TColor.h"

DLM_DtColor::DLM_DtColor():NumColors(18){

}

DLM_DtColor::~DLM_DtColor(){

}

int DLM_DtColor::GetColor(const int& WhichColor){
    switch(WhichColor){
    case 0 : return kRed;
    case 1 : return kBlue;
    case 2 : return kGreen+1;
    case 3 : return kMagenta;
    case 4 : return kCyan+1;
    case 5 : return kOrange;
    case 6 : return kRed+2;
    case 7 : return kBlue+2;
    case 8 : return kGreen-1;
    case 9 : return kMagenta+2;
    case 10: return kCyan+3;
    case 11: return kOrange-1;
    case 12: return kPink-8;
    case 13: return kViolet-9;
    case 14: return kAzure+10;
    case 15: return kTeal+10;
    case 16: return kSpring+9;
    case 17: return kYellow-6;
    default: return kBlack;
    }
}
