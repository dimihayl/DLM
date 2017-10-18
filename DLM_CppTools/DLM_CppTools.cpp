
#include "DLM_CppTools.h"

#include <stdio.h>
#include <string.h>

int ipow(int base, unsigned char exp){
    int result = 1;
    while (exp){
        if (exp & 1) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

unsigned uipow(unsigned base, unsigned char exp){
    unsigned result = 1;
    while (exp){
        if (exp & 1) result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

void ShowTime(long long Time, char * str, short show, bool compact){
    //[0] = years (1 y = 365 d), [1] = days, [2] = hours, [3] = minutes, [4] = seconds
    if(show<1 || show>5)
        show = 5;

    long long max_value[5];     //how much units of this unit constitute the "higher" unit (i.e. 365d = 1y, 24h = 1d etc.)
    long long value_in_sec[5];  //how much seconds is one unit
    long long values[5];        //values of each unit during the calculation
    char Unit_name[5][2][2][10];    //name of the units [][0][] -> singular; [][1][] -> plural; [][][0/1] -> full/compact form

    for(short i=0; i<5; i++)
        values[i] = 0;

    strcpy(Unit_name[0][0][0], "year");
    strcpy(Unit_name[0][1][0], "years");
    strcpy(Unit_name[0][0][1], "y");
    strcpy(Unit_name[0][1][1], "y");

    strcpy(Unit_name[1][0][0], "day");
    strcpy(Unit_name[1][1][0], "days");
    strcpy(Unit_name[1][0][1], "d");
    strcpy(Unit_name[1][1][1], "d");

    strcpy(Unit_name[2][0][0], "hour");
    strcpy(Unit_name[2][1][0], "hours");
    strcpy(Unit_name[2][0][1], "h");
    strcpy(Unit_name[2][1][1], "h");

    strcpy(Unit_name[3][0][0], "minute");
    strcpy(Unit_name[3][1][0], "minutes");
    strcpy(Unit_name[3][0][1], "min");
    strcpy(Unit_name[3][1][1], "min");

    strcpy(Unit_name[4][0][0], "second");
    strcpy(Unit_name[4][1][0], "seconds");
    strcpy(Unit_name[4][0][1], "s");
    strcpy(Unit_name[4][1][1], "s");

    max_value[0] = 3000000000*3000000000;
    max_value[1] = 365;
    max_value[2] = 24;
    max_value[3] = 60;
    max_value[4] = 60;

    value_in_sec[0] = 3600*24*365;
    value_in_sec[1] = 3600*24;
    value_in_sec[2] = 3600;
    value_in_sec[3] = 60;
    value_in_sec[4] = 1;

    //the first shown unit
    short first_unit = 0;

    //the last shown unit
    short last_unit;

    long long T_temp = Time;
    //sets values for y, d, h, min, s
    for(short i=0; i<5; i++){
        values[i] += T_temp/value_in_sec[i];
        T_temp -= values[i]*value_in_sec[i];
        if(i>=1){
            if(values[i] && !values[i-1] && !first_unit){
                first_unit = i;
            }
        }
    }

    last_unit = first_unit + show - 1;
    if(last_unit>4) last_unit = 4;
    if(!Time){
        first_unit = 4;
        last_unit = 4;
    }

    if(last_unit<4){
        if(values[last_unit+1] >= max_value[last_unit+1]/2){
            values[last_unit]++;
        }
    }

    for(short i=4; i>0; i--){
        if(values[i]>=max_value[i]){
            values[i] -= max_value[i];
            values[i-1]++;
        }
        if(values[i] && !values[i-1])
            first_unit = i;
    }
    if(values[0])
        first_unit = 0;

    last_unit = first_unit + show - 1;
    if(last_unit>4) last_unit = 4;
    if(!Time){
        first_unit = 4;
        last_unit = 4;
    }

    strcpy(str, "");
    char TempStr[24];

    for(short i=first_unit; i<=last_unit; i++){

        sprintf(TempStr, "%lld %s", values[i], Unit_name[i][values[i]!=1][compact]);
        strcat(str, TempStr);

        if(i<last_unit-1){
            if(compact)
                strcat(str, " ");
            else
                strcat(str, ", ");
        }
        else if(i==last_unit-1){
            if(compact)
                strcat(str, " ");
            else
                strcat(str, " and ");
        }
    }
}

DLM_Timer::DLM_Timer(){
    Start();
}
DLM_Timer::~DLM_Timer(){

}

void DLM_Timer::Start(){
    gettimeofday(&start, NULL);
}

long long DLM_Timer::Stop(){
    gettimeofday(&end, NULL);
    return (((long long)(end.tv_sec) - (long long)(start.tv_sec))*1000000 +
            ((long long)(end.tv_usec) - (long long)(start.tv_usec)));
}

