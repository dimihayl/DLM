#ifndef DLM_MATHFUNCTIONS_H
#define DLM_MATHFUNCTIONS_H

#include <vector>

//This function calculates the natural log of the gamma function for xx > 0
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

double GetPval(const double& chi2, const double& ndf);
double GetNsigma(const double& pval);
double GetNsigma(const double& chi2, const double& ndf);
//pars: [0] = goal in nsigma, [1] = ndf
//double FindDeltaChi2_GivenNsigma(const double& chi2, const double* pars);
//evaluate the difference that is needed in Chi2, that builds up a
//significance of nsigma. This is done based on the number of fit parameters.
double GetDeltaChi2(const double& nsigma, const unsigned& nfreepars);
double GetPvalFromNsig(const double& nsigma);

std::vector<std::vector<unsigned>> BinomialPermutations(const unsigned& N, const unsigned& k);

#endif
