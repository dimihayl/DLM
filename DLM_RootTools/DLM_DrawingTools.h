


#ifndef DLM_DrawingToolsH
#define DLM_DrawingToolsH

#include "TColor.h"


class DLM_DtColor{

private:

    const unsigned NumColors;

public:

    //width and height of the canvas
    DLM_DtColor();
    ~DLM_DtColor();

    int GetColor(const int& WhichColor);

};

#endif
