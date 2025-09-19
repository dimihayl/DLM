
#ifndef ANA_PI_D_RUN2_H
#define ANA_PI_D_RUN2_H

#include <fstream>
#include <vector>

void SetAnalysisFolderPath(char* folder_path);
void SetOutputFolderPath(char* folder_path);

void MainAnalysisFunction(char* Description, std::vector<int> fit_types, int mt_bin, int NumIter, bool Bootstrap=true, bool DataVar=true, bool FitVar=true, int SEED=0);


#endif
