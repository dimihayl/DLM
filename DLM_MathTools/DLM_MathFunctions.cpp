#include "DLM_MathFunctions.h"
#include "math.h"
#include "omp.h"
#include <iostream>

using namespace std;

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
