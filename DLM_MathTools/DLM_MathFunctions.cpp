#include "DLM_MathFunctions.h"
#include "DLM_RootFinder.h"
#include "math.h"
#include "omp.h"
#include <iostream>
#include "gsl_sf_gamma.h"

using namespace std;

const double HalfPi(1.570796326794897);
const double Pi(3.141592653589793);
const double TwoPi(6.283185307179586);

double** factrl_array=NULL;

double gammln(const double xx) {
	int j;
	double x,tmp,y,ser;
	static const double cof[14]={57.1562356658629235,-59.5979603554754912,
	14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
	.465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
	-.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
	.844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
	if (xx <= 0) throw("bad arg in gammln");
	y=x=xx;
	tmp = x+5.24218750000000000;
	tmp = (x+0.5)*log(tmp)-tmp;
	ser = 0.999999999999997092;
	for (j=0;j<14;j++) ser += cof[j]/++y;
	return tmp+log(2.5066282746310005*ser/x);
}
double factrl(const unsigned n) {
    if(n>170) return 0;
    if(!factrl_array){
        factrl_array = new double* [omp_get_num_procs()];
        for(int uThr=0; uThr<omp_get_num_procs(); uThr++){
            factrl_array[uThr] = NULL;
        }
    }
    const int WhichThread = omp_get_thread_num();
    if(!factrl_array[WhichThread]){
        factrl_array[WhichThread] = new double [171];
        factrl_array[WhichThread][0] = 1.;
        for (int i=1;i<171;i++) factrl_array[WhichThread][i] = i*factrl_array[WhichThread][i-1];
    }
	return factrl_array[WhichThread][n];
}

//copy-pasted from ROOT
double DLM_Poisson(double x, double par){

  // compute the Poisson distribution function for (x,par)
  // The Poisson PDF is implemented by means of Euler's Gamma-function
  // (for the factorial), so for all integer arguments it is correct.
  // BUT for non-integer values it IS NOT equal to the Poisson distribution.
   if (x<0)
      return 0;
   else if (x == 0.0)
      return 1./exp(par);
   else {
      double lnpoisson = x*log(par)-par-gammln(x+1.);
      return exp(lnpoisson);
   }
}

bool ApproxEqualF(const float val1, const float val2){
	if(fabs(val1-val2)/fabs(val1*val2)<5e-7) return true;
	return false;
}
bool ApproxEqualD(const double val1, const double val2){
	if(fabs(val1-val2)/fabs(val1*val2)<5e-15) return true;
	return false;	
}

double atanPhi(const double& y, const double& x){
	if(x==0){
		if(y>0) return HalfPi;
		else if(y<0) return -HalfPi;
		else return Pi;//actually it can be ANY number between 0 and TwoPi, so we take the average
	}
	else if(x>0&&y>0){
		return atan(fabs(y/x));
	}
	else if(x<0&&y>0){
		return Pi-atan(fabs(y/x));
	}
	else if(x<0&&y<0){
		return Pi+atan(fabs(y/x));
	}
	else{
		return TwoPi-atan(fabs(y/x));
	}
}
double AnglePhi(const double& py, const double& px){
	return atanPhi(py,px);
}
double AngleTheta(const double& pt, const double& pz){
	if(pt<0) return 0;
	if(!pz) return HalfPi;
	else if(pz>0) return atan(pt/pz);
	else return atan(pt/pz)+Pi;
}
double Pseudorapidity(const double& pt, const double& pz){
	if(pt<0) return 0;
	else if(pt==0) return 1000;
	else return asinh(pz/pt);
}
double Transverse(const double& py, const double& px){
	return sqrt(py*py+px*px);
}

double GetPval(const double& chi2, const double& ndf){
	return gsl_sf_gamma_inc_Q(ndf*0.5,chi2*0.5);
}

double GetNsigma(const double& chi2, const double& ndf){
//printf("chi2/ndf = %.2f/%.0f\n",chi2,ndf);
	return GetNsigma(GetPval(chi2,ndf));
}

//based on TMath::ErfcInverse (ROOT)
//return sqrt(2)*TMath::ErfcInverse(pval)
double GetNsigma(const double& pval){
	// Computes quantiles for standard normal distribution N(0, 1)
	// at probability p
	// ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3, 477-484.
//printf("pval = %f\n",pval);
	double p = 0.5*pval;

	if ((p<=0)||(p>=1)) {
		return 0;
	}

	const double a0 = 3.3871328727963666080e0;
	const double a1 = 1.3314166789178437745e+2;
	const double a2 = 1.9715909503065514427e+3;
	const double a3 = 1.3731693765509461125e+4;
	const double a4 = 4.5921953931549871457e+4;
	const double a5 = 6.7265770927008700853e+4;
	const double a6 = 3.3430575583588128105e+4;
	const double a7 = 2.5090809287301226727e+3;
	const double b1 = 4.2313330701600911252e+1;
	const double b2 = 6.8718700749205790830e+2;
	const double b3 = 5.3941960214247511077e+3;
	const double b4 = 2.1213794301586595867e+4;
	const double b5 = 3.9307895800092710610e+4;
	const double b6 = 2.8729085735721942674e+4;
	const double b7 = 5.2264952788528545610e+3;
	const double c0 = 1.42343711074968357734e0;
	const double c1 = 4.63033784615654529590e0;
	const double c2 = 5.76949722146069140550e0;
	const double c3 = 3.64784832476320460504e0;
	const double c4 = 1.27045825245236838258e0;
	const double c5 = 2.41780725177450611770e-1;
	const double c6 = 2.27238449892691845833e-2;
	const double c7 = 7.74545014278341407640e-4;
	const double d1 = 2.05319162663775882187e0;
	const double d2 = 1.67638483018380384940e0;
	const double d3 = 6.89767334985100004550e-1;
	const double d4 = 1.48103976427480074590e-1;
	const double d5 = 1.51986665636164571966e-2;
	const double d6 = 5.47593808499534494600e-4;
	const double d7 = 1.05075007164441684324e-9;
	const double e0 = 6.65790464350110377720e0;
	const double e1 = 5.46378491116411436990e0;
	const double e2 = 1.78482653991729133580e0;
	const double e3 = 2.96560571828504891230e-1;
	const double e4 = 2.65321895265761230930e-2;
	const double e5 = 1.24266094738807843860e-3;
	const double e6 = 2.71155556874348757815e-5;
	const double e7 = 2.01033439929228813265e-7;
	const double f1 = 5.99832206555887937690e-1;
	const double f2 = 1.36929880922735805310e-1;
	const double f3 = 1.48753612908506148525e-2;
	const double f4 = 7.86869131145613259100e-4;
	const double f5 = 1.84631831751005468180e-5;
	const double f6 = 1.42151175831644588870e-7;
	const double f7 = 2.04426310338993978564e-15;

	const double split1 = 0.425;
	const double split2 = 5.;
	const double konst1 = 0.180625;
	const double konst2 = 1.6;

	double q, r, quantile;
	q=p-0.5;
	if (fabs(q)<split1) {
		 r=konst1-q*q;
		 quantile = q* (((((((a7 * r + a6) * r + a5) * r + a4) * r + a3)
								* r + a2) * r + a1) * r + a0) /
								(((((((b7 * r + b6) * r + b5) * r + b4) * r + b3)
								* r + b2) * r + b1) * r + 1.);
	} else {
		 if(q<0) r=p;
		 else    r=1-p;
		 //error case
		 if (r<=0)
				quantile=0;
		 else {
				r=sqrt(-log(r));
				if (r<=split2) {
					 r=r-konst2;
					 quantile=(((((((c7 * r + c6) * r + c5) * r + c4) * r + c3)
										* r + c2) * r + c1) * r + c0) /
										(((((((d7 * r + d6) * r + d5) * r + d4) * r + d3)
										* r + d2) * r + d1) * r + 1);
				} else{
					 r=r-split2;
					 quantile=(((((((e7 * r + e6) * r + e5) * r + e4) * r + e3)
										* r + e2) * r + e1) * r + e0) /
										(((((((f7 * r + f6) * r + f5) * r + f4) * r + f3)
										* r + f2) * r + f1) * r + 1);
				}
				if (q<0) quantile=-quantile;
		 }
	}
	return -quantile;
}

//pars: [0] = goal in nsigma, [1] = ndf
double FindDeltaChi2_GivenNsigma(const double& chi2, const double* pars){
  return GetNsigma(chi2,pars[1])-pars[0];
}

double GetDeltaChi2(const double& nsigma, const unsigned& nfreepars){
  if(!nfreepars) return 1e128;
  double pars[2];
  pars[0] = nsigma;
  pars[1] = double(nfreepars);
  if(nsigma<1)  return NewtonRapson(FindDeltaChi2_GivenNsigma,0,double(nfreepars)*2.,pars);
  if(nsigma<3)  return NewtonRapson(FindDeltaChi2_GivenNsigma,0,double(nfreepars)*8.,pars);
  if(nsigma<5)  return NewtonRapson(FindDeltaChi2_GivenNsigma,0,double(nfreepars)*32.,pars);
  if(nsigma<10) return NewtonRapson(FindDeltaChi2_GivenNsigma,0,double(nfreepars)*128.,pars);
  return 1e128;
}

double GetPvalFromNsig(const double& nsigma){
	return GetPval(GetDeltaChi2(nsigma,1),1);
}

/*
void BinomialPermutations(const unsigned offset, const unsigned& k,
													std::vector<std::vector<unsigned>>& permutations, std::vector<unsigned>& elements){
	if(k==0){
		std::vector<unsigned> dummy;
		permutations.push_back(dummy);
		//printf("End of the line\n");
		return;
	}
	if(permutations.size()==0){
		std::vector<unsigned> dummy;
		permutations.push_back(dummy);
	}
	printf("Go on\n");
  for(unsigned ui=offset; ui<=elements.size()-k; ui++){
    permutations.back().push_back(ui);
		printf("pb %u\n", ui);
    BinomialPermutations(ui+1,k-1,permutations,elements);
    permutations.back().pop_back();
  }
}
std::vector<std::vector<unsigned>> BinomialPermutations(const unsigned& N, const unsigned& k){
	std::vector<std::vector<unsigned>> permutations;
	if(k==0||k>N){return permutations;}

	std::vector<unsigned> elements;
	for(unsigned uN=0; uN<N; uN++){
		elements.push_back(uN);
	}
	printf("Main body\n");
	BinomialPermutations(0,k,permutations,elements);
	return permutations;
}
*/

///////////////////////

void BinomialPermutations(unsigned offset, unsigned k, vector<unsigned>& elements, vector<vector<unsigned>>& permutations){
  if (k==0) {
    permutations.push_back(permutations.back());
    return;
  }
  if(!permutations.size()){
    permutations.emplace_back();
  }
  for (unsigned ui=offset; ui<=elements.size()-k; ui++){
    permutations.back().push_back(elements[ui]);
    BinomialPermutations(ui+1, k-1, elements, permutations);
    permutations.back().pop_back();
  }
}

std::vector<std::vector<unsigned>> BinomialPermutations(const unsigned& N, const unsigned& k){
  vector<unsigned> elements;
  vector<vector<unsigned>> permutations;
	if(k==0||k>N){return permutations;}

  for (unsigned ui = 0; ui < N; ui++) { elements.push_back(ui); }
  BinomialPermutations(0,k,elements,permutations);
	permutations.pop_back();
	return permutations;
}


//par[0] is an overall normalization
//than we have a pol4 = p0*(1+p1*k+p2*k^2+p3*k^3+p4*k^4), which has 3 free arguments and the following properties
//par4!=0 (pol4 flat at 0)
//	par1,par2 the two extrema, par3 is the p4, par4 is dummy
//par4==0&&par3!=0 (pol3)
//	par1,par2 the two extrema, par3 is p3
//par4==0&&par3==0&&par2!=0 (pol2)
//	par1 is the extrema, par2 is p2
//par4==0&&par3==0&&par2==0&&par1!=0 (pol1)
//	par1 is p1
//to avoid problems with a starting parameter of zero, to switch the order of the par we use -1e6 as a value
//a Mathematica computation of the equations is in your Femto folder
double DLM_Baseline(double* xval, double* par){

    double& k = *xval;
    double& p0 = par[0];
    //constrained polynomials

    double p1;
    double p2;
    double p3;
    double p4;
    if(par[4]!=-1e6){
        p4 = par[3];
        p3 = -4./3.*(par[1]+par[2])*p4;
        p2 = 2.*par[1]*par[2]*p4;
        p1 = 0;
    }
    else if(par[3]!=-1e6){
        p4 = 0;
        p3 = par[3];
        p2 = -1.5*(par[1]+par[2])*p3;
        p1 = 3.*par[1]*par[2]*p3;
    }
    else if(par[2]!=-1e6){
        p4 = 0;
        p3 = 0;
        p2 = par[2];
        p1 = -2.*par[1]*p2;
    }
    else{
        p4 = 0;
        p3 = 0;
        p2 = 0;
        p1 = par[1];
    }
    return p0*(1.+p1*k+p2*pow(k,2)+p3*pow(k,3)+p4*pow(k,4));

}
