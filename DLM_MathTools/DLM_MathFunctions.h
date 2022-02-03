#ifndef DLM_MATHFUNCTIONS_H
#define DLM_MATHFUNCTIONS_H

#include <vector>

double gammln(const double xx);
double factrl(const unsigned n);


//computes the value of phi, given the convention that it is in 0 to 2pi
//n.b. the atan2 gives values that are -pi to pi, and if you simply add pi to it
//the geometrical quadrants get swaped, thus we need this atanPhi
//N.B. Use these functions with care, when comparing to ROOT output, as in
//ROOT the standard atan2 function is used, which gives values between -pi and pi
double atanPhi(const double& y, const double& x);
double AnglePhi(const double& py, const double& px);
double AngleTheta(const double& pt, const double& pz);
double Pseudorapidity(const double& pt, const double& pz);
double Transverse(const double& py, const double& px);

std::vector<std::vector<unsigned>> BinomialPermutations(const unsigned& N, const unsigned& k);

#endif
