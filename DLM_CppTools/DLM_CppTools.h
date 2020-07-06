
#ifndef DLM_CppToolsH
#define DLM_CppToolsH

#include <sys/time.h>

extern long DLM_CppHelp_TimeS;
extern long DLM_CppHelp_TimeUS;

int ipow(int base, unsigned char exp);
unsigned uipow(unsigned base, unsigned char exp);

void ShowTime(long long T, char * str, short show=0, bool compact=true, short NumUnits=7);

//0 if file does not exist
//1 if file exists and is not opened
//-1 if the file exists and is opened
//the latter does not always function, it has to be opened by a process in a linux terminal
//it was proven to always work with .root files though!
int file_status(const char* name);

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

