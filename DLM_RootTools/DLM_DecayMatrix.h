#ifndef DLM_DECAYMATRIX_H
#define DLM_DECAYMATRIX_H

class TLorentzVector;

class DLM_DecayMatrix{
public:
    DLM_DecayMatrix();
    ~DLM_DecayMatrix();

    void SetFileName(const char* name);
    void SetHistoName(const char* name);
    void SetBins(const unsigned& numbins, const unsigned& kmin, const unsigned& kmax);
    void SetNumDaughters1(const unsigned& numdaughters);
    void SetDaughterMass1(const unsigned& daughter, const double& mass);
    void SetMotherMass1(const double& mass);
    void SetNumDaughters2(const unsigned& numdaughters);
    void SetDaughterMass2(const unsigned& daughter, const double& mass);
    void SetMotherMass2(const double& mass);
    void SetMeanMomentum(const double& mean);
    void SetMomentumSpread(const double& spread);
    //void SetUnitsMeV();
    //void SetUnitsGeV();
    void Run(const int& SEED, const unsigned& NumIter);

private:

    //bool UseGeV;
    char* FileName;
    char* HistoName;
    unsigned NumBins;
    double kMin;
    double kMax;
    unsigned NumDaughters1;
    unsigned NumDaughters2;
    double* Mass_Daughter1;
    double* Mass_Daughter2;
    double Mass1;
    double Mass2;
    double MomMean;
    double MomSpread;

    double relKcalc(const TLorentzVector& track1, const TLorentzVector& track2);
};
#endif // DLM_DECAYMATRIX_H
