
#include "DLM_Potentials_pp.h"
#include "DLM_StefanoPotentials.h"

#include <stdio.h>
#include <string.h>

int DlmPot=0;
int DlmPotFlag=0;
DLM_StefanoPotentials** ppV18pot=NULL;

void CleanUpV18Pot(){
    if(ppV18pot){
        for(unsigned iPot=0; iPot<30; iPot++){
            delete ppV18pot[iPot];
        }
        delete [] ppV18pot;
        ppV18pot = NULL;
    }
}

double ZeroPotential(double* Radius){
    return 0;
}

//I took it from Oli, he found it in some paper somewhere...
//the difference to the CRAB potential is only in the 3P1 (triplet) state
double ReidSoftCore(const int& polar, double* Radius){
    double r = Radius[0];
    double pmux,f1,f2,f3,f4,f7,vr;
    /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
     with Tensor forces neglected */
    if(polar==0){
        /* S=0 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
    }
    else if(polar>0){
        /* S=1 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f2=f1*f1;
        f3=f2*f1;
        f4=f2*f2;
        //f6=f4*f2;

        //ANNALS OF PHYSICS: 50, 411448 (1968)
        //This potential is originally included in CRAB (check), but describes the 3P2 - 3F2 mixing according to Reid paper
        //vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;

        //potential according to Reid:
        vr = 10.463 * ( (1. + 2./pmux + 2./(pmux*pmux))*f1 - (8./pmux + 2./(pmux*pmux))*f4 )/pmux - 135.25 * f2/pmux + 472.81 * f3/pmux;
    }

    return vr;//MeV
}

double ReidSoftCore1S0(double* Radius){
    return ReidSoftCore(0, Radius);
}

double ReidSoftCore3P(double* Radius){
    return ReidSoftCore(1, Radius);
}

//input in fm
double fReidMeVfm(const double& rad,const unsigned& polar){
    double r = rad;//convert in fm
    double pmux,f1,f2,f4,f6,f7,vr;
    /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
     with Tensor forces neglected */
    if(polar==0){
        /* S=0 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
        //vr=0;
    }
    else if(polar>0){
        /* S=1 */
        pmux=r*0.7;
        f1=exp(-pmux);
        f2=f1*f1;
        f4=f2*f2;
        f6=f4*f2;
        vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;
    }

    return vr;//MeV
}

double fReidMeVfm1S0(double* rad){
    return fReidMeVfm(*rad, 0);
}
double fReidMeVfm3P(double* rad){
    return fReidMeVfm(*rad, 1);
}

double fReidDlm(const double& rad,const unsigned short& s,const unsigned short& l,const unsigned short& j){
    const double mu_const = 0.7;
    const double h_const = 10.463;
    double f1,f2,f3,f4,f6,f7,vr;
    const double pmux=rad*mu_const;
    f1=exp(-pmux);
    f2=f1*f1;
    f3=f2*f1;
    f4=f2*f2;
    f6=f4*f2;
    f7=f4*f3;
    vr=0;
    if(s==0 && l==0 && j==0){
        /* 1S0 */
        f1=exp(-pmux);
        f4=(f1*f1*f1*f1);
        f7=f4*(f1*f1*f1);
        vr=-h_const*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
        //vr=0;
    }
    else if(s==1 && l==1 && j==0){
        /* 3P0 */
        vr =    -h_const*((1.+4./pmux+4./pmux/pmux)*f1 - (16./pmux+4./pmux/pmux)*f4)/pmux
                + 27.133*f2/pmux - 790.74*f4/pmux + 20662.*f7/pmux;
    }
    else if(s==1 && l==1 && j==1){
        /* 3P1 */
        vr =    h_const*((1.+2./pmux+2./pmux/pmux)*f1 - (8./pmux+2./pmux/pmux)*f4)/pmux
                - 135.25*f2/pmux + 472.81*f3/pmux;
    }
    else if(s==1 && l==1 && j==2){
        /* 3P2 - 3F2 */
        vr =    (h_const*f1/3. - 933.48*f4 + 4152.1*f6)/pmux;
    }
    return vr;//MeV
}

double fReidDlm1S0(double* rad){
    return fReidDlm(*rad,0,0,0);
}
double fReidDlm3P0(double* rad){
    return fReidDlm(*rad,1,1,0);
}
double fReidDlm3P1(double* rad){
    return fReidDlm(*rad,1,1,1);
}
double fReidDlm3P2(double* rad){
    return fReidDlm(*rad,1,1,2);
}
double fReidDlm3P(double* rad){
    return (fReidDlm(*rad,1,1,0)+fReidDlm(*rad,1,1,1)+fReidDlm(*rad,1,1,2))/3.;
}

//N.B. the potentials are never deleted, i.e. they stay there until termination
double fV18potential(const int& V18Pot, const int& Spin, const int& AngMom, const int& TotMom, double* Radius){
    if( (V18Pot<1||V18Pot>24) && (V18Pot<112||V18Pot>114) && (V18Pot<122||V18Pot>124) ){
        return 0;
    }
    if(!ppV18pot){
        ppV18pot = new DLM_StefanoPotentials* [30];
        for(unsigned uPot=0; uPot<30; uPot++){
            unsigned WhichPot = uPot+1;
            if(uPot==24) WhichPot=112;
            else if(uPot==24+1) WhichPot=113;
            else if(uPot==24+2) WhichPot=114;
            else if(uPot==24+3) WhichPot=122;
            else if(uPot==24+4) WhichPot=123;
            else if(uPot==24+5) WhichPot=124;
            ppV18pot[uPot] = new DLM_StefanoPotentials(WhichPot);
        }
    }
    unsigned StefPotId = V18Pot-1;
    if(V18Pot==112) StefPotId=24;
    else if(V18Pot==113) StefPotId=25;
    else if(V18Pot==114) StefPotId=26;
    else if(V18Pot==122) StefPotId=27;
    else if(V18Pot==123) StefPotId=28;
    else if(V18Pot==124) StefPotId=29;
//printf("StefPotId=%i --> V=%f\n",StefPotId,ppV18pot[StefPotId]->EvalCATS_v1_0(Radius[0],0));
    return ppV18pot[StefPotId]->Eval_PWprojector_pp(Radius[0],Spin,AngMom,TotMom,DlmPotFlag);
}

double fDlmPot(const int& Spin, const int& AngMom, const int& TotMom, double* Radius){
    //printf("V=%f\n",fV18potential(9,Spin,AngMom,TotMom,Radius)) ;
    switch(DlmPot){
        case pp_AV18 : return fV18potential(9,Spin,AngMom,TotMom,Radius);
        case pp_ReidV8 : return fV18potential(2,Spin,AngMom,TotMom,Radius);
        case pp_ReidSC : return fReidDlm(Radius[0],Spin,AngMom,TotMom);
        case pp_ReidOli : return ReidSoftCore(Spin,Radius);
        case pp_ReidCrab : return fReidMeVfm(Radius[0],Spin);
        default : return 0;
    }
}

double fDlmPot1S0(double* Radius){
    return fDlmPot(0,0,0,Radius);
}
double fDlmPot3P0(double* Radius){
    return fDlmPot(1,1,0,Radius);
}
double fDlmPot3P1(double* Radius){
    return fDlmPot(1,1,1,Radius);
}
double fDlmPot3P2(double* Radius){
    return fDlmPot(1,1,2,Radius);
}
double fDlmPot3P(double* Radius){
    return (fDlmPot3P0(Radius)+fDlmPot3P1(Radius)+fDlmPot3P2(Radius))/3.;
}

void GetDlmPotName(const int& potid, const int& potflag, char* name){
    switch(potid){
        case pp_AV18 :
            ppV18pot[8]->PotentialName(9, name);
            break;
        case pp_ReidV8 :
            ppV18pot[1]->PotentialName(2, name);
            break;
        case pp_ReidSC :
            strcpy(name,"Reid Soft-Core");
            break;
        case pp_ReidOli :
            strcpy(name,"Reid 3P2-3F2");
            break;
        case pp_ReidCrab :
            strcpy(name,"Reid 3P1");
            break;
        default :
            strcpy(name,"Unknown potential");
            break;
    }
    char Buffer[16];
    sprintf(Buffer, "%i", potflag);
    if(potflag!=0 && strcmp(name,"Unknown potential")){strcat(name, "^{(");strcat(name,Buffer);strcat(name,")}");}
}
