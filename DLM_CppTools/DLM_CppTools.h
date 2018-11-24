
#ifndef DLM_CppToolsH
#define DLM_CppToolsH

#include <sys/time.h>

extern long DLM_CppHelp_TimeS;
extern long DLM_CppHelp_TimeUS;

int ipow(int base, unsigned char exp);
unsigned uipow(unsigned base, unsigned char exp);

void ShowTime(long long T, char * str, short show=0, bool compact=true, short NumUnits=7);

class DLM_Timer{

protected:
    struct timeval start;
    struct timeval end;
public:
    DLM_Timer();
    ~DLM_Timer();
    void Start();
    long long Stop();
};

#endif

