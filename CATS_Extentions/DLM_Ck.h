#ifndef DLM_CK_H
#define DLM_CK_H

#include "CATS.h"
#include "CATStools.h"

class DLM_Ck : public DLM_Histo<double>{

public:
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double* bins);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar, CATS& cat, const unsigned& numbin, const double& minMom, const double& maxMom);
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
           const unsigned& numbin, const double* bins, double (*CorrFun)(const double&, const double*, const double*));
    DLM_Ck(const unsigned& nSourcePar, const unsigned& nPotPar,
           const unsigned& numbin, const double& minMom, const double& maxMom, double (*CorrFun)(const double&, const double*, const double*));
    ~DLM_Ck();

    //change UpToData when you do that
    void SetSourcePar(const unsigned& WhichPar, const double& Value);
    double GetSourcePar(const unsigned& WhichPar);
    unsigned GetNumSourcePar();
    void SetPotPar(const unsigned& WhichPar, const double& Value);
    double GetPotPar(const unsigned& WhichPar);
    //double GetPotPar(const unsigned& WhichPar);
    unsigned GetNumPotPar();
    //the momentum value after which the correlation will be interpolated
    //N.B. if we go outside the maximum momentum, the same is assumed
    //the interpolation is done so, as the C(k) is linearly propagated between the final computed point and C(kc) = 1.
    //If kc<maximum momentum C(k>maximum momentum) = C(maximum momentum)
    //If kc<0, than outside the maximum momentum C(k) = |kc|
    void SetCutOff(const double& Momentum=1e6, const double& kc=-1);
    double GetCutOff() const;
    double GetCutOff_kc() const;

    bool Status();
    void Update(const bool& FORCE=false);
    double Eval(const double& Momentum);

private:
    const unsigned NumSourcePar;
    const unsigned NumPotPar;
    //double* MomBinCopy;
    CATS* Kitty;
    double (*CkFunction)(const double&, const double*, const double*);
    double* SourcePar;
    bool SourceUpToDate;
    double* PotPar;
    bool PotUpToDate;
    double CutOff;
    double CutOff_kc;

    void DefaultConstructor();
};

#endif
