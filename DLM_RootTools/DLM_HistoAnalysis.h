

#ifndef DLM_HISTOANALYSIS_H
#define DLM_HISTOANALYSIS_H

#include "TMath.h"
#include "TH1D.h"

TH1D * MakeHistoFromTF1(TF1 * f1, double min, double max, unsigned int nbins = 100);
//gets the variance of a function in the interval [min, max]
//saves info about mean and mean2
double FunctionVar(TF1 * f1, double min, double max, double* mean=NULL, double* mean2=NULL, unsigned int steps = 100);

//Gets the central inteval of a functions, examines the function between min and max,
//alfa is the "amount" of integral that is desired, steps is the num of bins for the histogram,
//result is saved in from and to
//returns the median
void GetCentralInterval(TF1 * f1, const double& min, const double& max, double alpha, double& from, double& to, const unsigned int& steps = 100);
double GetCentralInterval(const TH1& h1, const double alpha, double& from, double& to,  bool extrapolate=false);
void DrawCentralInterval(TH1F*& h1, double Median, double LowReach, double UpReach, int NumNewBins,
                         TH1F*& h1main, TH1F*& h1bLow, TH1F*& h1bLow2, TH1F*& h1bUp, TH1F*& h1bUp2, TH1F*& h1median);
#endif



