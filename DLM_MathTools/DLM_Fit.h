#ifndef DLMFIT_H
#define DLMFIT_H

#include <vector>

template <class Type> class DLM_Histo;
class DLM_Random;

struct DLM_FitSolution{
  DLM_FitSolution();
  DLM_FitSolution(const DLM_FitSolution& other);
  ~DLM_FitSolution();
  std::vector<float>* Par;
  //if negative: not evaluated yet
  float Chi2;
  unsigned ID;
  bool operator>(const DLM_FitSolution& other) const{return Chi2>other.Chi2;}
  bool operator<(const DLM_FitSolution& other) const{return Chi2<other.Chi2;}
  bool operator==(const DLM_FitSolution& other) const{return Chi2==other.Chi2;}
  DLM_FitSolution operator-(const DLM_FitSolution& other) const;
  DLM_FitSolution operator+(const DLM_FitSolution& other) const;
  void Add(const DLM_FitSolution& other, const float& scale);
  void operator+=(const DLM_FitSolution& other);
  void operator-=(const DLM_FitSolution& other);
  void operator*=(const float& value);
  //void operator=(const DLM_FitSolution& other);
  DLM_FitSolution& operator=(const DLM_FitSolution& other);
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
  DLM_Fit(const unsigned num_threads=0);
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
  //!!!! THIS FUNCTION IS NOT THREAD-SAFE
  std::vector<DLM_Histo<float>*> Eval(const std::vector<float>& pars);



  //a single "fit" step, where we look at the current BestSols and random sample
  //some parameters to test next. At the end, we add the new solutions to Solution
  //and reevaluate the BestSols
  //void WalkAround(std::vector<float>& pars);


  //do multiple walk arounds, until we have found our final solutions.
  //the limit to break is either that our first 10 solutions are identical (i.e. no dependence on parameters)
  //or we are in a situation where within the first 100 best solutions we have at least 68 that end up within
  //the 1 sigma band with respect to the best solution
  //This function is automatically using OMP for multi-threading, however it is by itself NOT threadsafe,
  //i.e. it should be called only once!!
  void Fit();

  void SetParLimits(const unsigned& WhichPar, const float& min, const float& max);
  //only these orders of magnitudes will be scanned
  //N.B. This function can conflict with the SetParLimits, be careful in usage!!!
  //this function is used to navigate the random sampling in the beginning better
  //N.B. we still allow the main algorithm to test any values within the ParLimits
  void SetParOrdMag(const unsigned& WhichPar, const int& min, const int& max);
  void SetParMinMag(const unsigned& WhichPar, const int& min);
  void SetParMaxMag(const unsigned& WhichPar, const int& max);
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
  unsigned Npar(const bool& free_pars=true) const;//done
  float GetParameter() const;
  float GetParError() const;
  float GetParLowLimit() const;
  float GetParUpLimit() const;
  float GetParRange() const;
  bool ParIsFree(const unsigned& uPar) const;
  bool ParIsFixed(const unsigned& uPar) const;

  bool DEBUG_PrepareForWalk();
  void DEBUG_WanderAround(){WanderAround();}
  std::vector<DLM_FitSolution>* DEBUG_GetSolution() {return Solution;}

private:
  //DLM_Histo<float>* Data;
  //DLM_Histo<float>* Model;

  void (*FitFnct)(const std::vector<float>&, std::vector<DLM_Histo<float>*>&);
  unsigned NumPars;
  unsigned NumData;
  bool InterProcess;

  void OrderSolutions();

//std::vector<DLM_FitSolution>& desired_walks
  //creates a list for the set of parameters that are to be tested next
  bool PrepareForWalk();
  //starts the evaluation of the parameters from above
  //this is done either in parallel, on a single machine, or as a job
  void WanderAround();
  //the basis function of Eval, that is threadsafe
  bool ProbePars(const std::vector<float>* pars, const unsigned& ThId);
  //collects all the data from the 'wandering'
  //trivial if the wandering was done locally, otherwise the output of each
  //job has to be read out first.
  void MemoryLane();
  //the thread id of the currently best solution saved
  unsigned ThIdBestSol();

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
  std::vector<float>* Par;
  //Lower/Upper limit. Those are set to zero if the paremter is fixed
  std::vector<float>* ParL;
  std::vector<float>* ParU;
  std::vector<float>* ParOL;
  std::vector<float>* ParOU;



  std::vector<float>* DataLow;
  std::vector<float>* DataUp;

  unsigned ThIdSol;
  float* chi2;
  unsigned* NumDataPts;

  //should have the dim of the variables
  std::vector<DLM_Histo<float>*>* Data;
  //std::vector<DLM_Histo<float>*>* Model;
  std::vector<DLM_Histo<float>*>** Model;

  //saves the paramters and chi2 of each function call. This vector is ordered
  //according to the chi2 only after we find the best solution and will be used
  //to determine the final uncertainties
  std::vector<DLM_FitSolution>* Solution;
  //keeps track of the best solutions found so far
  //this vector contains the parameters and chi2, and is ordered
  //such that the last element is has the worst chi2
  //std::vector<DLM_FitSolution>* BestSols;
  unsigned NumWildCards;
  unsigned NumBestSols;

  DLM_Random** RanGen;
  unsigned NumThreads;
  const unsigned MaxThreads;

};




class DLM_SA_Fit{

public:
  DLM_SA_Fit(const unsigned& ndata, const unsigned& npars);
  ~DLM_SA_Fit();

  typedef double (*SA_FitFunction)(double*, double*);

  void SetParameter(const unsigned ipar, const float value);
  void SetParLimits(const unsigned ipar, const float value_min, const float value_max);
  //void SetParLogScale(const bool yesno);
  //how many evaluation step we will do +/- the default.
  //by default this is 5. The total number of iterations is 1 + par_st*2
  //the minimum value is 2, maximum 200
  void SetParSteps(const unsigned ipar, const unsigned par_st);
  //while we fit, is the fitter allowed to make the binning finer to evaluate the parameter 
  //to better precision. The refinement is done by reducing the parameter space for each parameter by a fixed factor.
  //The factor is by default 1.5, but it can be changed to any value between 1.2 and 10.
  //We can also control the number of refinements to do. Allowed range: 0 to 100 (default 30)
  void SetParRefinements(const unsigned ref_num, const float ref_fac);
  //default 0.25. This parameter gives us a cut-off at which we stop refining the parameter space, given in units 
  //of chi2. I.e. if neighbouring grid pts do not result in larger chi2 differences, we stop the whole fit
  //lower value will give more precision, but will increase computational time.
  void SetEpsChi2(const double eps_chi2);
  //fit range
  void SetDataRange(const unsigned idata, const float value_min, const float value_max);
  void SetData(const unsigned idata, DLM_Histo<float>& hdata);
  //how are zero bins with zero error treated?
  //by default the minimal error is 1e-37, however this can be modified here.
  //if set to zero, the bins will simply be rejected. If some other number, this number will be set to 
  //the error of all bins of zero error
  void SetZeroBinError(const double& zbe);
  //the parameter_ids should be a vector, containing the ID of the parameter (as defined within this object)
  //which should be used by the function fit_fun. E.g. parameter_ids = {0,5,6} would imply that fit_fun needs 3 parameters,
  //end they correspond ot the 0-th, 5-th and 6-th parameters set by SetParameter(const unsigned ipar, const float value)
  void SetFitFunction(const unsigned idata, SA_FitFunction fit_fun, std::vector<int>& parameter_ids);
  //you can also use a single fit fucntion, in which case double* x of the function is assumed to cover the data sets,
  //i.e. x[i] is the x-axis if the i-th data set. The parameters are do not hava ids any more
  //void SetFitFunction(SA_FitFunction fit_fun);

  float GetParameter(const unsigned& ipar);
  float GetParError(const unsigned& ipar);
  void GetParLimits(unsigned ipar, float& parmin, float& parmax) const;
  bool GetParLogScale();
  unsigned GetParSteps(const unsigned ipar);
  unsigned GetParRefinements();
  double GetEpsChi2();
  void GetDataRange(unsigned idata, float& datamin, float& datamax) const;

  float GetChisquare();
  int GetNDF();
  unsigned GetNumberFitPoints();
  unsigned GetNumberFreeParameters();
  float GetProb();

  void SetRandomSeed(unsigned rnd_seed);

  bool GoBabyGo();

      //double Chi2(const DLM_Histo& other){
       // if(!SameStructure(other)) return -1;
       // double chi2_val = 0;
      // for(unsigned uBin=0; uBin<TotNumBins+2; uBin++){
      //      chi2_val += pow((BinValue[uBin]-other.BinValue[uBin])/sqrt(BinError[uBin]*BinError[uBin]+other.BinError[uBin]*other.BinError[uBin]),2.);
      //  }
      //  return chi2_val;
   // }

  //bool operator+=(const DLM_SA_Fit& other);
  //bool operator/=(const float& value);

  //bool operator=(const DLM_SA_Fit& other);
  //bool operator=(const float& value);

private:

  const unsigned N_Data;
  const unsigned N_Pars;
  const unsigned min_n_grid_pts;
  const unsigned max_n_grid_pts;
  const unsigned max_n_solutions;

  std::vector<float>* Par;
  //the general min/max for the fit
  std::vector<float>* ParMin;
  std::vector<float>* ParMax;
  std::vector<unsigned>* ParSteps;

  //solutions
  std::vector<float>** ParSol;
  std::vector<float>** ParErr;

  //check if 1D
  DLM_Histo<float>** input_data;
  std::vector<float>* fit_min;
  std::vector<float>* fit_max;
  SA_FitFunction* theory_model;
  std::vector<int>** par_id;
  //-1 = not defined
  //0 = 1 big global function
  //1 = we have 1 function per data set
  //WILL NOT WORK, THE OPTION 1 IS NOW OFF
  int function_type;

  //an array of grids, 1 for solution node found
  DLM_Histo<float>** chi2_grid;

  double epsChi2;
  double zero_bin_err;
  unsigned num_refinements;
  double refinment_factor;

  unsigned num_solutions;

  double currentChi2;
  double currentBestChi2;
  unsigned num_data_pts;
  unsigned num_free_fit_pars;

  //to be used often in the code, defined here for convenience to avoid repeated mem allocation
  unsigned* which_bin;
  unsigned* which_bin_best;
  std::vector<unsigned>* tot_best_bin;



  void del_chi2_grid();
  //std::vector<unsigned> unflatten_index(unsigned linear_index);
  //unsigned flatten_index(const std::vector<unsigned>& coords);
  float EvalChisquare(double* pars);

  unsigned ExploreAndNavigate();

  std::vector<unsigned> get_all_neighbours(unsigned isol, unsigned tot_bin_id);
  std::vector<unsigned> get_the_neighborhood(unsigned isol, unsigned tot_bin_id);
  

  DLM_Random* rangen;

};

/*
//we need from the start the number of fit parameters
//we need from the start the dimension of the data (e.g. a fit in k* and mt is 2D)
class DLM_SaGd_Fit{

public:
  DLM_SaGd_Fit(const unsigned& ndata, const unsigned& npars);
  ~DLM_SaGd_Fit();

  void SetUpVar(const unsigned& ivar, const double& var_min, const double& var_max);
  //void SetUpPar(const unsigned& ipar, const double& par_min, const double& par_max);



private:

  DLM_FitParameters* fit_pars;
  //std::vector<float>* Var;

  DLM_Histo<float>* data_histo;
  DLM_Histo<float>* theo_histo;

};
*/



#endif // DLMFIT_H
