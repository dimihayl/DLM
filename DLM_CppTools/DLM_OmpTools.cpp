#include "DLM_OmpTools.h"
#include "DLM_CppTools.h"

#include "omp.h"

#include "math.h"

unsigned GetNumThreads(){
  unsigned Nth;
  #pragma omp parallel
  {
    #pragma omp master
    {
      Nth=omp_get_num_threads();
    }
  }
  return Nth;
}

unsigned GetOptimalThreads(){
  unsigned Nth;
  #pragma omp parallel
  {
    #pragma omp master
    {
      Nth = omp_get_num_threads();
    }
  }
  unsigned Cth=Nth;
  unsigned Oth=Nth;;

  const unsigned IterStep = 1000000;
  DLM_Timer Timer;
  long long ltime_mus;
  double dtime_ms;
  double baseline_ms;
  double temp;
  double reduction_factor;
  bool BreakOut = false;
  do{
    //bool TestOngoing = true;
    Timer.Start();
    omp_set_dynamic(0);
    omp_set_num_threads(Cth);
    while(true){
      double DUMMY = 0;
      #pragma omp parallel for private(DUMMY)
      for(unsigned uIter=0; uIter<IterStep; uIter++){
        if(uIter==0){
          printf("%u (%u)\n",omp_get_num_threads(),Nth);
        }
        for(unsigned u=0; u<10000; u++){
          DUMMY += log(double(uIter));
        }
      }
      ltime_mus = Timer.Stop();
      dtime_ms = double(ltime_mus)*0.001;
      //if we are too fast, we continue iterating
      //if(dtime_ms<time_ms) continue;

      printf(" -> Cth = %u, t = %f (%f)\n",Cth,dtime_ms,baseline_ms);
      if(Cth==Nth){
        baseline_ms = dtime_ms;
        if(Cth<2){BreakOut=true;break;}
        Cth /= 2;
        break;
      }
      else{
        //if we have essentially the same time, we take that as a good solution
        //if( fabs(dtime_ms-baseline_ms)/baseline_ms < 0.1 || dtime_ms<baseline_ms){
        if(true){
          Oth = Cth;
          if(Cth<2){BreakOut=true;break;}
          Cth /= 2;
          break;
        }
        //if we have a worst solution
        else{
          BreakOut=true;
          break;
        }
      }
      Timer.Start();
      break;
    }//while
  }
  while(!BreakOut);
  return Oth?Oth:1;
}
