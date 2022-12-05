
#ifndef DLM_DrawingToolsH
#define DLM_DrawingToolsH

//#include "RtypesCore.h"

class TColor;
class TArrayI;

class DLM_DtColor{

private:

    const unsigned NumColors;
    const bool PresetForCatsPaper;
    TColor** MyColors;
    TArrayI* RainbowArray;

    //float* ColR,ColG,ColB;
    //Color_t* ColT;

public:

    //width and height of the canvas
    DLM_DtColor();
    DLM_DtColor(const unsigned& NumCol);
    ~DLM_DtColor();

    //this may fail to save the new colors in the output file...
    //to get around this, we can also take a color from the predefined 255 colors of root (rainbow)
    //the latter will only work if we have less than 255 number of colors
    int GetColor(const int& WhichColor);
    int GetRainbow(const int& WhichColor);


    //Color_t GetColorT(const int& WhichColor);

};

#endif
