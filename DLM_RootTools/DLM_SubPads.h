
#ifndef DLM_SubPadsH
#define DLM_SubPadsH

#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"

enum DLM_PadPos { kPadLeft, kPadTop, kPadRight, kPadBottom };
enum DLM_TextPos { kTextLeft, kTextBottom, kTextRight, kTextTop };
enum DLM_MarPos { kMarLeft, kMarRight, kMarBottom, kMarTop };

enum DLM_WH { kWidth, kHeight };

class DLM_SubPads{

private:
    TCanvas* MainCanvas;
    const UInt_t CanvasPadWidth;
    const UInt_t CanvasPadHeight;

    //the mem. allocated for SubPads
    UShort_t MaxSubPads;
    //the default steps in which the MaxSubPads is increased, in case
    //the memory runs out.
    const UShort_t MaxSubPadsStep;

    UShort_t NumSubPads;

    TPad** SubPad;
    //[PadNr][left/right/bottom/top]
    Float_t** SubPadPostion;
    //Ratio of the SubPad size with respect to the MainPad,
    //i.e. how many times smaller is the SubPad
    //[PadNr][Width/Height]
    Float_t** ScaleFactor;

    void InitialValues();
    void Initialize();
    void AllocMoreMemForSubPad();

public:

    //width and height of the canvas
    DLM_SubPads(const UInt_t PW, const UInt_t PH);
    ~DLM_SubPads();

    //adds a sub-pad. Returns false if the margins conflict
    //with an existing sub-pad. The coordinates are given, as per default,
    //with respect to the bottom-left corner of the main pad.
    //After executing this function the current directory will be set to the new TPad
    Bool_t AddSubPad(Float_t Left, Float_t Right, Float_t Bottom, Float_t Top);

    //the same as AddSubPad, but the coordinate is with reversed y-axis and
    //the 0 point is at the top-left corner
    Bool_t AddSubPadTL(Float_t Left, Float_t Right, Float_t Top, Float_t Bottom);//lrtb

    //convert the coordinates for the position on the canvas given
    //with respect to a centre in the Bottom left side of the canvas (as by default)
    //to a coordinate system with a centre at the top left side of the canvas and reversed y-axis.
    //the same function can be used for the reversed conversion as well.
    Float_t ConvertYPadPos(Float_t ypos);

    //selects the WhichPad -th SubPad. If WhichPad==65535 (==-1 for unsigned), than the MainPad as selected
    Bool_t cd(UShort_t WhichPad);
    void SetLogx(UShort_t WhichPad, Bool_t lg=true);
    void SetLogy(UShort_t WhichPad, Bool_t lg=true);
    void SetGridx(UShort_t WhichPad, Bool_t val=true);
    void SetGridy(UShort_t WhichPad, Bool_t val=true);
    void SetGrid(UShort_t WhichPad, Bool_t val=true);
    //Frees up all memory.
    //! This is necessary to use, if one wishes to delete MainCanvas from the same function in which
    //the DLM_SubPads was created. The reason for that is that all Pads etc. need to be deleted BEFORE
    //the TCanvas. Therefore use Destroy() before delete MainCanvas
    void Destroy();
    //Destroys + initializes again
    void Reset();

    //Set the margins of the pad. if MainCoordinates==true (default) the margins are set
    //in the coordinate system of the MainPad, not the local once!
    Bool_t SetMargin(UShort_t PadNr, Float_t Left, Float_t Right, Float_t Bottom, Float_t Top, Bool_t MainCoordinates=true);//lrbt

    Float_t GetScaleWidth(const UShort_t& PadNr);
    Float_t GetScaleHeight(const UShort_t& PadNr);
    TCanvas* GetCanvas();

    void SetLabelSize(const UShort_t& PadNr, TAxis* Axis, const Float_t& Size);
    void SetTitleSize(const UShort_t& PadNr, TAxis* Axis, const Float_t& Size);

    enum WhichPad { kMainPad=65535 };
};

#endif

