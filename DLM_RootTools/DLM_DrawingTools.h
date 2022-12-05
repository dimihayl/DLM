
#ifndef DLM_DrawingToolsH
#define DLM_DrawingToolsH

class TColor;

class DLM_DtColor{

private:

    const unsigned NumColors;
    const bool PresetForCatsPaper;
    TColor** MyColors;

public:

    //width and height of the canvas
    DLM_DtColor();
    DLM_DtColor(const unsigned& NumCol);
    ~DLM_DtColor();

    int GetColor(const int& WhichColor);

};

#endif
