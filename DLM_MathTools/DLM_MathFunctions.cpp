#include "DLM_MathFunctions.h"
#include "math.h"
#include "omp.h"
#include <iostream>

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
