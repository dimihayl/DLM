#ifndef SCATTERINGPARAMETERS_H
#define SCATTERINGPARAMETERS_H

#include <string>
#include <vector>

#include "DLM_Histo.h"

class CATS;
class DLM_Random;

struct DLM_PotSp{
    double p0=0;
    double p1=0;
    double f=0;
    double d=0;
    long num_entries=0;

    double p0_min = 1e64;
    double p0_max = -1e64;
    double p1_min = 1e64;
    double p1_max = -1e64;

    DLM_PotSp(){DLM_PotSp(0);}
    DLM_PotSp(const double value){*this = 0;}

    void set(double p0_, double p1_, double f_, double d_, long entries_ = 1);
    bool operator+=(const DLM_PotSp& other);
    bool operator/=(const double& value);

    bool operator=(const DLM_PotSp& other);
    bool operator=(const double& value);
    DLM_PotSp operator*(const double& value);
    void compute_average();
};

class DLM_ScatteringPars{

public:
    typedef double (*PotFunction)(double*);

    DLM_ScatteringPars();
    ~DLM_ScatteringPars();

    //f is the potential function, assumed to follow the same definitions as in DLM_Potentials (used for CATS)
    //in these potentials [0] and [1] from Pars are reserved for r*,k*
    //The Nth part though is assumed to be at [2+N].
    //The p0 and p1 here define at which N are the parameters which we want to scan over
    //IMPORTANT!!! Ones p0 and p1 are set, within this class they are ALWAYS accessed with flags 0 and 1
    void SetPotFun(PotFunction f, unsigned p0=0, unsigned p1=1);
    //needed if we have more than two parameters, and some have to be fixed.
    //it is assumed that the vector contains the exact number of pars used by the potential
    void SetPotPar(std::vector<double> pot_pars_);
    void SetParLimits(unsigned WhichPar, double par_min, double par_max);
    void SetParGrid(unsigned WhichPar, unsigned num_bins);
    //do we want to sample with a log law. This may be a problem if zero is within the sampling range.
    //Thus min_abs defines a small region around zero par<|min_abs| that is EXCLUDED from the sampling
    //void SetParLogScale(bool yesno, double min_abs = 1e-32);

    //the desired scattering length f, and effective range d
    void SetTarget_f(double min_f, double max_f);
    void SetSlGrid(unsigned num_bins);
    void SetTarget_d(double min_d, double max_d);
    void SetErGrid(unsigned num_bins);

    void SetRedMass(double red_mass_);

    //void SetNumThreads(unsigned num_thr);
    void SetRandomSeed(unsigned seed);
    
    //completely random sampling within the limits of the potential parameters
    void RandomScan(unsigned NumSamples);
    //we sample target_fraction of the time by trying to sample around pot pars that already 
    //provided a (f,d) combo within the desired limits. The remaining 1-target_fraction we sample with RandomScan
    void TargetedScan(unsigned NumSamples, float target_fraction = 0.8);

    std::vector<double> GetPotPars(double f0, double d0);

    //should be one for the class to work
    double Occupancy();
    unsigned GetNumEntries();

    //saves the dlm_PotSp_Map into a file. If loaded, one can continue working with it
    void Save(std::string file_name, bool Overwrite=false);
    //dumps into an ascii file the settings used for the simulation.
    //no automatic loading available yet
    void SaveSettings(std::string file_name, bool Overwrite=false);
    void Load(std::string file_name, bool reset = false);

private:

    const double kMin;
    const double kMax;
    const unsigned kSteps;
    double EPS;

    //the level of agreement (fractional) our 3 evaluation points of the 
    //scattering parameters need to have to accept the solution
    const double eps_f;
    const double eps_d;

    unsigned id_p0;
    unsigned id_p1;

    unsigned bins_par0;
    double pot_par0[2];
    unsigned bins_par1;
    double pot_par1[2];

    unsigned bins_f;
    double target_f[2];
    unsigned bins_d;
    double target_d[2];

    double red_mass;


    //2D in f0 and d0
    DLM_Histo<DLM_PotSp>* dlm_PotSp_Map;
    DLM_Histo<DLM_PotSp>* dlm_PotSp_AvgMap;

    DLM_Random* rangen;

    PotFunction potential;
    
    CATS* Kitty;
    double* pot_pars;
    unsigned num_pot_pars;

    bool GetScatteringParameters(double& f, double& fe, double& d, double& de);
    void Reset();
    void SampleSomeStuff(unsigned NumSamples, double min_p0, double max_p0, double min_p1, double max_p1);
    


};


#endif
