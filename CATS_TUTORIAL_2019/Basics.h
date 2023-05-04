

const double Mass_Pi0 = 134.98;
const double Mass_Pic = 139.57;
const double Mass_p = 938.272;
const double Mass_L = 1115.683;
const double Mass_Xim = 1321.7;

class TGraph;

double Basics_Potential_Usmani(double* pars);
double Basics_Source_Gauss(double* Pars);
//double CATS_FIT_PL(double* x, double* pars);

TGraph* Basics_PiPiTheory(const bool& Identical, const bool& WithCoulomb);
TGraph* Basics_PiPiCATS(const bool& Identical, const bool& WithCoulomb);
void Basics_ProtonLambda();
