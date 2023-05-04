
class TString;
class DLM_Ck;
class DLM_CkDecomposition;

void Ck_pL_Ledni_Usmani();
//save in a file the correlation function from DLM_Ck
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_Ck* CkToPlot);
//save in a file the correlation function from DLM_CkDecomposition
void RootFile_DlmCk(const TString& RootFileName, const TString& GraphName, DLM_CkDecomposition* CkToPlot, bool PlotIndividualContributions=false);

void Ck_pp_Decomposition(const TString& SourceType);
