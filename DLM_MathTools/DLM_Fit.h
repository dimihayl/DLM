#ifndef DLMFIT_H
#define DLMFIT_H

#include <vector>

template <class Type> class DLM_Histo;


//general idea:
//* perform multiple fit steps
//* what is a single fit step:
//  - random sample the paramters (uniformly) within the limits
//  - evaluate the chi2 for each sample, order the solutions based on their chi2
//  - for each step, evaluate an "intermediate" solution, where the accepted range of parameters
//    are those up to chi2 + numfreepars^2.
//    N.B. take special care to define allowed regions, e.g. in case there are values a,b,c of some par,
//    where a and c are good, but b is bad, than we cant claim (a,b) is the solution
//  - pick some fraction of best solutions
//  - have a function, that can create a secondary DLM_Fit object, tuned to a best solution
//    The parameters are set to the best solution, the limits are reduced according to some formula
//* now we need to be performing multiple steps. We can do that using two strategies, both based on the following logic:
//  create our secondary DLM_Fit objects and execute them. When all of them finish executing,
//  we collect the solutions, again do the whole procedure of finding the best solutions based on chi2
//  This is repeated until we have reached a certain desired precision in the par limits of some steps.
//  Dynimically, a good condition to stop is when the limits are smaller than the uncertainties
//  However, we separate two strategies of solving all of this:
//  1)  LOCAL
//      Everyting is evaluate locally (i.e. within the same process) and saved within this class
//  2)  INTER-PROCESS
//      Each solution can be executed as a separate process. To do this, we need to have a mother process
//      and an output folder with temp files. Each function call of Eval than creates a small temp file which
//      is to be read out by the the mother. Ones the mother creates all required funciton calls, it deletes the temp files
//      and continues the exectution.
class DLM_Fit{
public:
  DLM_Fit();
  ~DLM_Fit();

  void SetData(DLM_Histo<float>& data);

  //sets the evaluatuation function and number of parameters
  //resets all previous information regarding the fit function
  void SetEval(double (*evalfnc)(const std::vector<float>&, const std::vector<float>&, std::vector<DLM_Histo<float>*>&),
                const unsigned numpars);
  DLM_Histo<float>* Eval(const std::vector<float>& vars);
  DLM_Histo<float>* Eval(const float& var);

  void SetParLimits(const unsigned& WhichPar, const float& min, const float& max);
  void SetParameter(const unsigned& WhichPar, const float& val);
  void FixParameter(const unsigned& WhichPar, const float& val);

  //number of iterations per fit step
  void SetNumIterPerStep();
private:
  //DLM_Histo<float>* Data;
  //DLM_Histo<float>* Model;

  std::vector<float> Par;
  std::vector<float> Var;
  std::vector<DLM_Histo<float>*> Data;
  std::vector<DLM_Histo<float>*> Model;


};

#endif // DLMFIT_H
