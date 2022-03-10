#ifndef DLMFIT_H
#define DLMFIT_H

#include <vector>

template <class Type> class DLM_Histo;
class DLM_Random;

struct DLM_FitSolution{
  DLM_FitSolution(){Par = new std::vector<float> (); Chi2=-1; ID=0;}
  ~DLM_FitSolution(){delete Par;}
  std::vector<float>* Par;
  //if negative: not evaluated yet
  float Chi2;
  unsigned ID;
  bool operator>(const DLM_FitSolution& other) const{return Chi2>other.Chi2;}
  bool operator<(const DLM_FitSolution& other) const{return Chi2<other.Chi2;}
  bool operator==(const DLM_FitSolution& other) const{return Chi2==other.Chi2;}
  void operator=(const DLM_FitSolution& other) {
    Par=other.Par;
    Chi2=other.Chi2;
    ID=other.ID;
  }
  //the final name is basename + ID.dfs
  void SaveToFile(const char* basename);
  //execute a scipt that can evaluate this process
  void InterProcessEval(const char* exefile);
  void ReadFromFile(const char* basename);
};

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

  //resets all about the function and requires later a specific number of data sets and parameters
  bool SetUp(const unsigned& numdata, const unsigned& numpars);

  //we can add multiple histograms to fit
  //N.B. the FitFnct should be able to evaluate ALL of them
  bool SetData(const unsigned& WhichSet, DLM_Histo<float>& data);

  //sets the evaluation function. The input arguments of this function should be the parameters.
  //the output is saved in in vector of DLM_Histos. The histos are supposed to exist already.
  //resets all previous information regarding the fit function
  void SetFitFnct(void (*fitfnc)(const std::vector<float>&, std::vector<DLM_Histo<float>*>&));
  //evaluates the function for all variables, in the process computing Chi2, Err etc.
  //DLM_Histo<float>* Eval(const std::vector<float>& vars);
  //DLM_Histo<float>* Eval(const float& var);

  //evaluates the function for a certain set of parameters
  //it is the function used internally for fitting
  //returns the Model DLM_Histo within this class, which contains the output
  std::vector<DLM_Histo<float>*> Eval(const std::vector<float>& pars);



  //a single "fit" step, where we look at the current BestSols and random sample
  //some parameters to test next. At the end, we add the new solutions to Solution
  //and reevaluate the BestSols
  //void WalkAround(std::vector<float>& pars);


  //do multiple walk arounds, until we have found our final solutions.
  //the limit to break is either that our first 10 solutions are identical (i.e. no dependence on parameters)
  //or we are in a situation where within the first 100 best solutions we have at least 68 that end up within
  //the 1 sigma band with respect to the best solution
  void Fit();

  void SetParLimits(const unsigned& WhichPar, const float& min, const float& max);
  void SetParameter(const unsigned& WhichPar, const float& val);
  void FixParameter(const unsigned& WhichPar, const float& val);

  void SetFitRange(const unsigned& WhichData, const unsigned WhichAxis,
                    const float& lower, const float& upper);

  //number of best solutions selected per iteratopms
  void SetNumBestSols(const unsigned& num);
  //how many of the best solutions should we study
  //minimum of 1, maximum of NumBestSols. Allows to better disentangle multiple minima.
  //N.B. the CPU time goes as BestSols*WildCards, i.e. set to 1 in case
  //you know there is a single mimimum within the parameter space
  void SetNumWildCards(const unsigned& num);
  //if false, the evaluation is perfomed by this process.
  void SetInterProcess(const bool& yesno=true);
  void SetNumThreads(const unsigned& num);
  void SetSeed(const unsigned& thread, const unsigned& seed);

  float Chi2() const;
  float Chi2Ndf() const;
  float Pval() const;
  float Nsigma() const;
  unsigned Ndf() const;
  unsigned Npts() const;
  unsigned Npar(const bool& free_pars=true) const;
  float GetParameter() const;
  float GetParError() const;
  float GetParLowLimit() const;
  float GetParUpLimit() const;
  float GetParRange() const;
private:
  //DLM_Histo<float>* Data;
  //DLM_Histo<float>* Model;

  void (*FitFnct)(const std::vector<float>&, std::vector<DLM_Histo<float>*>&);
  unsigned NumPars;
  unsigned NumData;

  void OrderSolutions();

//std::vector<DLM_FitSolution>& desired_walks
  //creates a list for the set of parameters that are to be tested next,
  //and than exectutes each one of them in parallel, either on a single machine or as a job
  void WalkAround();


/*
  //parameters
  std::vector<float> Par;
  std::vector<float> ParErr;
  //Lower/Upper limit
  std::vector<float> ParL;
  std::vector<float> ParU;
  //range. The moment ParR<ParErr, it means we have converged for this parameter.
  std::vector<float> ParR;

  //variables
  std::vector<float> Var;
*/

  std::vector<float>* DataLow;
  std::vector<float>* DataUp;


  float chi2;
  unsigned NumDataPts;

  //should have the dim of the variables
  std::vector<DLM_Histo<float>*>* Data;
  std::vector<DLM_Histo<float>*>* Model;

  //saves the paramters and chi2 of each function call. This vector is ordered
  //according to the chi2 only after we find the best solution and will be used
  //to determine the final uncertainties
  std::vector<DLM_FitSolution>* Solution;
  //keeps track of the best solutions found so far
  //this vector contains the parameters and chi2, and is ordered
  //such that the last element is has the worst chi2
  std::vector<DLM_FitSolution>* BestSols;
  unsigned NumWildCards;
  unsigned NumBestSols;

  DLM_Random** RanGen;
  unsigned NumThreads;
  const unsigned MaxThreads;

};

#endif // DLMFIT_H
