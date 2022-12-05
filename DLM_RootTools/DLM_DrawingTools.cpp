#include "DLM_DrawingTools.h"
#include "TColor.h"
#include "TArrayI.h"
#include "TMath.h"
#include "TStyle.h"

DLM_DtColor::DLM_DtColor():NumColors(18),PresetForCatsPaper(true){
  MyColors=NULL;
  RainbowArray=NULL;
  //ColT=NULL;
  //ColR=NULL;
  //ColG=NULL;
  //ColB=NULL;
}
DLM_DtColor::DLM_DtColor(const unsigned& NumCol):NumColors(NumCol),PresetForCatsPaper(false){
  MyColors=NULL;
  RainbowArray=NULL;
  //ColT=NULL;
  //ColR=NULL;
  //ColG=NULL;
  //ColB=NULL;
  if(NumColors<=1) return;
  MyColors = new TColor* [NumColors];
  RainbowArray = new TArrayI();
  //auto OldPal = gStyle->GetPalette();
  gStyle->SetPalette(kRainBow);
  *RainbowArray = TColor::GetPalette();
  //gStyle->SetPalette(OldPal);
  //ColT = new Color_t [NumColors];
  //ColR=new float[NumColors];
  //ColG=new float[NumColors];
  //ColB=new float[NumColors];
  float RgbRed;
  float RgbGreen;
  float RgbBlue;
  float Progress;
  const float Frac = 1./4.5;
  for(unsigned uCol=0; uCol<NumCol; uCol++){
    Progress = float(uCol)/float(NumColors-1);
    if(Progress<Frac){RgbRed=1; RgbGreen=Progress*4.5; RgbBlue=0;}
    else if(Progress<2.*Frac){RgbRed=1-(Progress-Frac)*4.5; RgbGreen=1; RgbBlue=0;}
    else if(Progress<3.*Frac){RgbRed=0; RgbGreen=1; RgbBlue=(Progress-2.*Frac)*4.5;}
    else if(Progress<4.*Frac) {RgbRed=0; RgbGreen=1-(Progress-3.*Frac)*4.5; RgbBlue=1;}
    else {RgbRed=(Progress-0.889)*4.5; RgbGreen=0; RgbBlue=1;}

    if(RgbRed<=0) RgbRed=0;
    if(RgbRed>=1) RgbRed=1;
    if(RgbGreen<=0) RgbGreen=0;
    if(RgbGreen>=1) RgbGreen=1;
    if(RgbBlue<=0) RgbBlue=0;
    if(RgbBlue>=1) RgbBlue=1;

    //RgbRed = 0.99;
    //RgbBlue = 0.01;
    //RgbGreen = 0.01;

    //Printf("uCol = %i: RgbRed=%f, RgbGreen=%f, RgbBlue=%f", uCol, RgbRed, RgbGreen, RgbBlue);

    MyColors[uCol] = new TColor(RgbRed, RgbGreen, RgbBlue);
    //ColT[uCol].red = 1;
    //ColT[uCol].green = 1;
    //ColT[uCol].blue = 1;
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
  if(RainbowArray){
    delete RainbowArray;
    RainbowArray = NULL;
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

int DLM_DtColor::GetRainbow(const int& WhichColor){
  if(WhichColor>255) return 0;

  int RBA = TMath::Nint(float(WhichColor)/float(NumColors)*256.);
  return RainbowArray->At(RBA);
}
