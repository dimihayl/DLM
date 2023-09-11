
#include "DLM_CppTools.h"

#include <stdio.h>
//#include <string.h>

#include <unistd.h>
#include <fcntl.h>
#include <cerrno>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <stdexcept>
#include <math.h>


bool isNumber(const char* str) {
    std::istringstream iss(str);
    double number;

    // Try to read the value as a double
    if (!(iss >> number)) {
        // Conversion failed or input is not a valid number
        return false;
    }

    // Check if there are any remaining characters (e.g., scientific notation)
    char remaining;
    if (iss >> remaining) {
        // There are remaining characters after parsing the number
        return false;
    }

    // The input is a valid number
    return true;
}

bool isInteger(const char* str){
  if(!isNumber(str)){
    return false;
  }
  double value = atof(str);
  if(value != round(value)) return false;
  else return true;
}

std::vector<std::string> ParseString(const std::string& text, const std::string& delim){
  std::vector<std::string> tokens;
  size_t start = 0U;
  size_t end = text.find(delim);
  while (end != std::string::npos){
    tokens.push_back(text.substr(start, end - start));
    start = end + delim.length();
    end = text.find(delim, start);
  }
  //std::string last = text.substr(start, end - start);
  //if(last!="" && last !="\n") tokens.push_back(last);
  tokens.push_back(text.substr(start, end - start));
  return tokens;
}


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

void ShowTime(long long Time, char * str, short show, bool compact, short NumUnits){
    //const unsigned NumUnits=6;
    if(NumUnits<1 || NumUnits>7) NumUnits=7;
    //[0] = years (1 y = 365 d), [1] = days, [2] = hours, [3] = minutes, [4] = seconds, [5] = ms
    if(show<1 || show>NumUnits)
        show = NumUnits;

    long long max_value[NumUnits];     //how much units of this unit constitute the "higher" unit (i.e. 365d = 1y, 24h = 1d etc.)
    long long value_in_sec[NumUnits];  //how much milliseconds is one unit
    long long values[NumUnits];        //values of each unit during the calculation
    char Unit_name[NumUnits][2][2][16];    //name of the units [][0][] -> singular; [][1][] -> plural; [][][0/1] -> full/compact form

    for(short sUnit=0; sUnit<NumUnits; sUnit++){
        values[sUnit] = 0;
        switch(sUnit){
        case 0 :
            strcpy(Unit_name[0][0][0], "year");
            strcpy(Unit_name[0][1][0], "years");
            strcpy(Unit_name[0][0][1], "y");
            strcpy(Unit_name[0][1][1], "y");
            break;
        case 1 :
            strcpy(Unit_name[1][0][0], "day");
            strcpy(Unit_name[1][1][0], "days");
            strcpy(Unit_name[1][0][1], "d");
            strcpy(Unit_name[1][1][1], "d");
            break;
        case 2 :
            strcpy(Unit_name[2][0][0], "hour");
            strcpy(Unit_name[2][1][0], "hours");
            strcpy(Unit_name[2][0][1], "h");
            strcpy(Unit_name[2][1][1], "h");
        break;
        case 3 :
            strcpy(Unit_name[3][0][0], "minute");
            strcpy(Unit_name[3][1][0], "minutes");
            strcpy(Unit_name[3][0][1], "min");
            strcpy(Unit_name[3][1][1], "min");
            break;
        case 4 :
            strcpy(Unit_name[4][0][0], "second");
            strcpy(Unit_name[4][1][0], "seconds");
            strcpy(Unit_name[4][0][1], "s");
            strcpy(Unit_name[4][1][1], "s");
            break;
        case 5 :
            strcpy(Unit_name[5][0][0], "millisecond");
            strcpy(Unit_name[5][1][0], "milliseconds");
            strcpy(Unit_name[5][0][1], "ms");
            strcpy(Unit_name[5][1][1], "ms");
            break;
        case 6 :
            strcpy(Unit_name[6][0][0], "microsecond");
            strcpy(Unit_name[6][1][0], "microseconds");
            strcpy(Unit_name[6][0][1], "μs");
            strcpy(Unit_name[6][1][1], "μs");
            break;
        default :  break;
        }
    }

    max_value[0] = 3000000000*3000000000;
    max_value[1] = 365;
    max_value[2] = 24;
    max_value[3] = 60;
    max_value[4] = 60;
    max_value[5] = 1000;
    max_value[6] = 1000;

    value_in_sec[NumUnits-1]=1;
    for(short sUnit=NumUnits-2; sUnit>=0; sUnit--){
        switch(sUnit){
        case 0 : value_in_sec[0] = value_in_sec[1]*365; break;
        case 1 : value_in_sec[1] = value_in_sec[2]*24; break;
        case 2 : value_in_sec[2] = value_in_sec[3]*60; break;
        case 3 : value_in_sec[3] = value_in_sec[4]*60; break;
        case 4 : value_in_sec[4] = value_in_sec[5]*1000; break;
        case 5 : value_in_sec[5] = value_in_sec[6]*1000; break;
        default :  break;
        }
    }

    //the first shown unit
    short first_unit = 0;

    //the last shown unit
    short last_unit;

    long long T_temp = Time;
    //sets values for y, d, h, min, s
    for(short i=0; i<NumUnits; i++){
        values[i] += T_temp/value_in_sec[i];
        T_temp -= values[i]*value_in_sec[i];
        if(i>=1){
            if(values[i] && !values[i-1] && !first_unit){
                first_unit = i;
            }
        }
    }

    last_unit = first_unit + show - 1;
    if(last_unit>NumUnits-1) last_unit = NumUnits-1;
    if(!Time){
        first_unit = NumUnits-1;
        last_unit = NumUnits-1;
    }

    if(last_unit<NumUnits-1){
        if(values[last_unit+1] >= max_value[last_unit+1]/2){
            values[last_unit]++;
        }
    }

    for(short i=NumUnits-1; i>0; i--){
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
    if(last_unit>NumUnits-1) last_unit = NumUnits-1;
    if(!Time){
        first_unit = NumUnits-1;
        last_unit = NumUnits-1;
    }

    strcpy(str, Time>=0?"":"- ");
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

//0 no such file
//1 file exists and is NOT opened
//-1 file exists and is opened
//-2 ERROR
int file_status(const char* name){
    #ifdef F_SETLEASE
    int file_stat = open(name, O_RDONLY);
    if (file_stat < 0) {
        //file does not exist
        close(file_stat);
        return 0;
    }
    if (fcntl(file_stat, F_SETLEASE, F_WRLCK) && EAGAIN == errno) {
        close(file_stat);
        return -1;
    }
    else {
        fcntl(file_stat, F_SETLEASE, F_UNLCK);
        close(file_stat);
        return 1;
    }
    #else
    printf("\033[1;33mWARNING:\033[0m The file_status function will only work with Linux.\n");
    return -2;
    #endif // F_SETLEASE
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
