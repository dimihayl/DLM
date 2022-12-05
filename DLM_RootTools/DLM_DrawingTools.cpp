#include "DLM_DrawingTools.h"
#include "TColor.h"

DLM_DtColor::DLM_DtColor():NumColors(18),PresetForCatsPaper(true){
  MyColors=NULL;
}
DLM_DtColor::DLM_DtColor(const unsigned& NumCol):NumColors(NumCol),PresetForCatsPaper(false){
  MyColors=NULL;
  if(NumColors<=1) return;
  MyColors = new TColor* [NumColors];
  float RgbRed;
  float RgbGreen;
  float RgbBlue;
  float Progress;
  for(unsigned uCol=0; uCol<NumCol; uCol++){
    Progress = float(uCol)/float(NumColors-1);
    if(Progress<0.222){RgbRed=1; RgbGreen=Progress*4.5; RgbBlue=0;}
    else if(Progress<0.444){RgbRed=1-(Progress-0.222)*4.5; RgbGreen=1; RgbBlue=0;}
    else if(Progress<0.667){RgbRed=0; RgbGreen=1; RgbBlue=(Progress-0.444)*4.5;}
    else if(Progress<0.889) {RgbRed=0; RgbGreen=1-(Progress-0.667)*4.5; RgbBlue=1;}
    else {RgbRed=(Progress-0.889)*4.5; RgbGreen=0; RgbBlue=1;}
    //Printf("uCol = %i: RgbRed=%f, RgbGreen=%f, RgbBlue=%f", uCol, RgbRed, RgbGreen, RgbBlue);
    MyColors[uCol] = new TColor(float(5000+uCol), RgbRed, RgbGreen, RgbBlue);
  }
}


DLM_DtColor::~DLM_DtColor(){
  if(MyColors){
    for(unsigned uCol=0;  uCol<NumColors; uCol++){
      delete MyColors[uCol];
    }
    delete [] MyColors;
    MyColors = NULL;
  }
}

int DLM_DtColor::GetColor(const int& WhichColor){
  if(PresetForCatsPaper){
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
  else{
    if(NumColors<=1){
      return 0;
    }
    else{
      return MyColors[WhichColor]->GetNumber();
    }
  }
}
