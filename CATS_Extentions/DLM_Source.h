
#ifndef DLM_SOURCE_H
#define DLM_SOURCE_H

double GaussSource(double* Pars);
double GaussSourceTF1(double* x, double* Pars);
double GaussSourceTheta(double* Pars);
double CauchySource(double* Pars);
double CauchySourceTheta(double* Pars);

double DoubleGaussSource(double* Pars);
double GaussCauchySource(double* Pars);

//a monte-carlo out-side-long Gaussian source. Works very slowly!
double GaussOSL_MC(double* Pars);

double Gauss_Exp_Approx(double* Pars);
double Gauss_Exp(double* Pars);
double GaussExpSimple(double* Pars);
double GaussExpTotSimple(double* Pars);
double GaussExpTotIdenticalSimple(double* Pars);
double GaussExpTotIdenticalSimple_2body(double* Pars);
double GaussExpTotSimple_2body(double* Pars);

double MemberSourceForwarder(void* context, double* Pars);
class MemberSource{
public:
    virtual double Eval(double* Pars);
    double Eval(const double& Momentum, const double Radius, const double& Angle);
private:
    double PARS[3];
};

class MS_Gauss:public MemberSource{
public:
    double Size;
    double Eval(double* Pars);
};

//A source taking into account resonances and mT scaling. The Simple part is that resonances are back to back
//and in case of two resonances, the t*p/m are just added up. Also we use the approximate relation for small t*p/m
class MS_GaussExp_mT_Simple:public MemberSource{
public:
    MS_GaussExp_mT_Simple();
    ~MS_GaussExp_mT_Simple();

    void SetNum_mT(const unsigned& nmt);
    void SetMean_mT(const unsigned& umt, const double& mmt);
    void SetWeight_mT(const unsigned& umt, const double& wmt);
    void SetLinear_mT(const double& lin);
    void SetSlope_mT(const double& slope);
    void SetMass(const unsigned short& particle, const double& mass);
    void SetMassR(const unsigned short& particle, const double& mass);
    void SetMassD(const unsigned short& particle, const double& mass);
    void SetTau(const unsigned short& particle, const double& tau);
    void SetResonanceWeight(const unsigned short& particle, const double& weight);

    double Eval(double* Pars);
private:
    unsigned Num_mT;//number of mT bins (common for the particle pair)
    //[mT]
    double* Mean_mT;//the mean mT in each bin
    double* Weight_mT;//the weight of each mT bin
    double Linear_mT;//a+b*mT = r, linear is a, slope is b
    double Slope_mT;
    //[particle 1,2]
    double* Mass;//mass of the main particles we investigate
    double* MassR;//mass of the resonance from which each particles comes from
    double* MassD;//mass of the second daughter in the decay of the resonance
    double* Tau;//the mean lifetime of the resonance
    double* Weight_R;//the amount of resonances

    double* Parameters;
};
#endif
