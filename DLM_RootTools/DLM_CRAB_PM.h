#ifndef DLM_CRAB_PM_H
#define DLM_CRAB_PM_H

#include "TRandom3.h"
#include "TH1F.h"

using namespace std;

class DLM_CRAB_PM{

private:
    string OutputFileName;
    //[NumSpecies][x,y,z,tau]
    double** SpaceCoord;
    //double* rx;
    //double* ry;
    //double* rz;
    //double* tau;

    unsigned NumEvents;
    unsigned NumPerEvent;

    unsigned NumSpecies;
    double Temperature;
    int* ParticleSpecies;
    double* ParticleMass;

    void MemInit(const unsigned& num);
    void MemCleanup();

    //[NumSpecies]
    //info about the shape of the source for each particle species
    //if negative than the shape corresponds to some preset function.
    //if non-negative, than it corresponds to a specific distribution, that can
    //be set via its CDF (see below).
    //! so far I have ONLY implemented the histo case (non-negative)
    //[NumSpecies][x,y,z,tau]
    char** Shape;
    //char* ShapeRx;
    //char* ShapeRy;
    //char* ShapeRz;
    //char* ShapeTau;

    //List of pointers to CDF, that should be set by hand.
    //Do note, that if a CDF is NULL and a there is a request to sample from it,
    //a Shape of type 0 will be assumed (some default gauss). A warning message should appear.
    TH1F*** CDF;

    TH1F** IntegratedSource;

    const double Pi;
public:
    DLM_CRAB_PM();
    ~DLM_CRAB_PM();
double GetRandSpace(TRandom3& rangen, unsigned Species, int Coord);
    void SetNumEvents(unsigned num);
    void SetNumPerEvent(unsigned num);
    void SetOutputFile(string OutputName);
    void SetNumSpecies(unsigned num);
    void SetParticleMass(unsigned np, double mass);
    void SetParticleMass(const double* mass);
    void SetParticleSpecies(unsigned np, int species);
    void SetParticleSpecies(const int* species);
    void SetTemperature(double temp);

    void SetRx(unsigned np, double value);
    void SetRx(const double* value);

    void SetRy(unsigned np, double value);
    void SetRy(const double* value);

    void SetRz(unsigned np, double value);
    void SetRz(const double* value);

    void SetRxyz(unsigned np, double value);
    void SetRxyz(const double* value);

    void SetTau(unsigned np, double value);
    void SetTau(const double* value);

    void SetShapeRx(unsigned np, char value);
    void SetShapeRx(const char* value);
    void SetShapeRx(unsigned np, TH1F* cdf);

    void SetShapeRy(unsigned np, char value);
    void SetShapeRy(const char* value);
    void SetShapeRy(unsigned np, TH1F* cdf);

    void SetShapeRz(unsigned np, char value);
    void SetShapeRz(const char* value);
    void SetShapeRz(unsigned np, TH1F* cdf);

    void SetShapeRxyz(unsigned np, char value);
    void SetShapeRxyz(const char* value);
    void SetShapeRxyz(unsigned np, TH1F* cdf);

    void SetShapeTau(unsigned np, char value);
    void SetShapeTau(const char* value);
    void SetShapeTau(unsigned np, TH1F* cdf);

    void RunPhasemaker(unsigned Seed=0);

    TH1F* GetIntegratedSource(const unsigned& Species);

    //enum ParticlePair { pp, pLambda};
    enum PM_COORD_POS { tau,rx,ry,rz };
    enum PM_COORD_ENE { energy,kx,ky,kz };
};




#endif // DLM_CRAB_H

