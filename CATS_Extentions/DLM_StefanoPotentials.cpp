
#include "DLM_StefanoPotentials.h"


//!TEST
#include <unistd.h>
#include <string.h>

DLM_StefanoPotentials::DLM_StefanoPotentials(const int& WhichPot,
            const int& XMN, const int& GAM, const int& RHO, const int& CHI, const int& OMG, const int& FTP):
    m_pi(138.03),m_pi_not(134.98),m_pi_plus(139.570),m_n(939.565),
    g_a(1.267),g_a_Tre(1.29),f_pi(92.4),pi(3.1415926535897),hc197(197.327053),
    g34 ( 1.22541670246517),EGM_C1(-0.00081),EGM_C3(-0.0034),EGM_C4(0.0034),
    Lam(1000.0),small(1e-4),vsmall(1e-10),xmn(XMN),gam(GAM),rho(RHO),chi(CHI),omg(OMG),ftp(FTP),
    mp(938.27231*xmn),mn(939.56563*xmn),lpot(WhichPot){
//printf("WhichPot = %i\n", WhichPot);
    vv = new double [18];
    vp = new double [12];
    vw = new double [14];

    pm0 = new double [12];
    pm1 = new double [12];
    pgva01 = new double [12];
    pgva11 = new double [12];
    pgva00 = new double [12];
    pgva10 = new double [12];
    pgvb01 = new double [12];
    pgvb11 = new double [12];
    pgvb00 = new double [12];
    pgvb10 = new double [12];
    pgvls1 = new double [12];
    pgvt1 = new double [12];
    pgvso21 = new double [12];
    pgvls0 = new double [12];
    pgvt0 = new double [12];
    pgvso20 = new double [12];

    if(lpot<100) setpot();
    NLO_C = new double [7];
    N2LO_C = new double [7];

}

DLM_StefanoPotentials::~DLM_StefanoPotentials(){
    delete [] vv;
    delete [] vp;
    delete [] vw;

    delete [] pm0;
    delete [] pm1;
    delete [] pgva01;
    delete [] pgva11;
    delete [] pgva00;
    delete [] pgva10;
    delete [] pgvb01;
    delete [] pgvb11;
    delete [] pgvb00;
    delete [] pgvb10;
    delete [] pgvls1;
    delete [] pgvt1;
    delete [] pgvso21;
    delete [] pgvls0;
    delete [] pgvt0;
    delete [] pgvso20;

    delete [] NLO_C;
    delete [] N2LO_C;
}


double DLM_StefanoPotentials::yc(const double& t, const double& x){
    return exp(-t)/x;
}
double DLM_StefanoPotentials::yt(const double& t, const double& x){
    return (1.+3./t+3./pow(t,2))*exp(-t)/x;
}
double DLM_StefanoPotentials::yls(const double& t, const double& x){
    return -(1.+t)*exp(-t)/pow(x,3);
}
double DLM_StefanoPotentials::yl2(const double& t, const double& x){
    return (1.+2./t)*exp(-t)/pow(x,3);
}
double DLM_StefanoPotentials::pc(const double& t){
    return exp(-t)/t;
}
double DLM_StefanoPotentials::pt(const double& t){
    return (1.+3./t+3./pow(t,2));
}
double DLM_StefanoPotentials::pls(const double& t){
    return (1.+1./t)/t;
}
double DLM_StefanoPotentials::pso2(const double& t){
    return (1.+3./t+3./pow(t,2))/pow(t,2);
}

/*
void DLM_StefanoPotentials::SetXmn(const double& val){
    xmn = val;
}

double DLM_StefanoPotentials::GetXmn(){
    return xmn;
}

void DLM_StefanoPotentials::SetGam(const double& val){
    gam = val;
}

double DLM_StefanoPotentials::GetGam(){
    return gam;
}


void DLM_StefanoPotentials::SetRho(const double& val){
    rho = val;
}

double DLM_StefanoPotentials::GetRho(){
    return rho;
}

void DLM_StefanoPotentials::SetChi(const double& val){
    chi = val;
}

double DLM_StefanoPotentials::GetOmg(){
    return omg;
}

void DLM_StefanoPotentials::SetOmg(const double& val){
    omg = val;
}

double DLM_StefanoPotentials::GetChi(){
    return chi;
}

void DLM_StefanoPotentials::SetFtp(const double& val){
    ftp = val;
}

double DLM_StefanoPotentials::GetFtp(){
    return ftp;
}
*/

void DLM_StefanoPotentials::setpot(){
    //double h2m = 41.47;
    //double h2mcsb=0.;

    //if(lpot>=4 && lpot<=6){
        //h2m=197.33*197.33/(938.9*xmn);
    //}
    //else if( (lpot>=7 && lpot<=17) || lpot>=21){
        //h2m=0.5*hc197*hc197*(1./mp+1./mn);
        //h2mcsb=0.5*hc197*hc197*(1./mp-1./mn);
    //}

    tnr=1.;
    if (lpot==3) tnr=exp(-gam*rho);
    if ((lpot>=4 && lpot<=17) || lpot>=21){
        pimass1=gam;
        pimass2=rho;
        pimass3=chi;
        wsrange=1.-(2./3.)*(omg-1.);
        ftpec=ftp;
    }

    pm1[0]=0.684026; pm1[1]=1.6; pm1[2]=2.3; pm1[3]=3.0;
    pm1[4]=3.7; pm1[5]=4.4; pm1[6]=5.1; pm1[7]=5.8;
    pm1[8]=6.5; pm1[9]=8.2; pm1[10]=9.9; pm1[11]=11.3;

    pm0[0]=0.699536; pm0[1]=1.6; pm0[2]=2.3; pm0[3]=3.0;
    pm0[4]=3.7; pm0[5]=4.4; pm0[6]=5.1; pm0[7]=5.8;
    pm0[8]=6.5; pm0[9]=8.2; pm0[10]=9.9; pm0[11]=11.3;

    pgva01[0]=-10.077427; pgva01[1]=-120.49564; pgva01[2]=-212.36460; pgva01[3]=-8717.4198;
    pgva01[4]=54383.377; pgva01[5]=-213421.47; pgva01[6]=494583.57; pgva01[7]=-667153.34;
    pgva01[8]=529575.98; pgva01[9]=-137034.12; pgva01[10]=-346971.94; pgva01[11]=0;

    pgva11[0]=3.3591422; pgva11[1]=-86.479568; pgva11[2]=-465.93111; pgva11[3]=1867.3085;
    pgva11[4]=3850.9213; pgva11[5]=-19674.338; pgva11[6]=123231.40; pgva11[7]=-314493.61;
    pgva11[8]=242424.40; pgva11[9]=166904.04; pgva11[10]=-485343.64; pgva11[11]=0;

    pgvb01[0]=0.0026851393; pgvb01[1]=0.051092455; pgvb01[2]=-0.84264258; pgvb01[3]=14.736312;
    pgvb01[4]=-145.21993; pgvb01[5]=841.58389; pgvb01[6]=-2786.1170; pgvb01[7]=5056.4510;
    pgvb01[8]=-3367.4205; pgvb01[9]=-1784.5529; pgvb01[10]=5354.8266; pgvb01[11]=0;

    pgvb11[0]=-0.00089504644; pgvb11[1]=0.037488481; pgvb11[2]=-0.89373089; pgvb11[3]=14.123475;
    pgvb11[4]=-146.60152; pgvb11[5]=841.91462; pgvb11[6]=-2839.4273; pgvb11[7]=5265.3427;
    pgvb11[8]=-3500.0430; pgvb11[9]=-2487.9479; pgvb11[10]=7306.8121; pgvb11[11]=0;

    pgvls1[0]=0; pgvls1[1]=-426.00359; pgvls1[2]=26279.517; pgvls1[3]=-575570.33;
    pgvls1[4]=6003393.4; pgvls1[5]=-34519443.; pgvls1[6]=113554590.; pgvls1[7]=-207292090.;
    pgvls1[8]=171315480.; pgvls1[9]=-86418222.; pgvls1[10]=0; pgvls1[11]=0;

    pgvt1[0]=3.3591422; pgvt1[1]=-0.85945824; pgvt1[2]=-104.76340; pgvt1[3]=1262.9465;
    pgvt1[4]=-18881.061; pgvt1[5]=106132.46; pgvt1[6]=-332119.10; pgvt1[7]=555857.62;
    pgvt1[8]=-349166.64; pgvt1[9]=-119450.13; pgvt1[10]=0; pgvt1[11]=0;

    pgvso21[0]=0; pgvso21[1]=-0.52218640; pgvso21[2]=186.44558; pgvso21[3]=-3709.1115;
    pgvso21[4]=55913.117; pgvso21[5]=-369985.60; pgvso21[6]=1453754.3; pgvso21[7]=-3135247.1;
    pgvso21[8]=2433908.1; pgvso21[9]=0; pgvso21[10]=0; pgvso21[11]=0;

    pgva00[0]=32.290874; pgva00[1]=-82.465631; pgva00[2]=1232.9384; pgva00[3]=-16859.879;
    pgva00[4]=172926.83; pgva00[5]=-768352.77; pgva00[6]=2189047.5; pgva00[7]=-3844728.7;
    pgva00[8]=2799055.9; pgva00[9]=502518.28; pgva00[10]=-2600612.4; pgva00[11]=0;

    pgva10[0]=-10.763625; pgva10[1]=-42.973669; pgva10[2]=-718.56844; pgva10[3]=4246.9120;
    pgva10[4]=-34574.024; pgva10[5]=126711.69; pgva10[6]=-274168.41; pgva10[7]=529607.24;
    pgva10[8]=-366067.13; pgva10[9]=-223036.73; pgva10[10]=406838.33; pgva10[11]=0;

    pgvb00[0]=-0.0085980096; pgvb00[1]=0.026814385; pgvb00[2]=-1.3280693; pgvb00[3]=10.324289;
    pgvb00[4]=-115.27067; pgvb00[5]=694.56175; pgvb00[6]=-2387.9335; pgvb00[7]=4238.8011;
    pgvb00[8]=-2452.1604; pgvb00[9]=-1951.2821; pgvb00[10]=4180.1160; pgvb00[11]=0;

    pgvb10[0]=0.0028660032; pgvb10[1]=-0.00081798046; pgvb10[2]=-0.53314560; pgvb10[3]=0.83162030;
    pgvb10[4]=-31.192395; pgvb10[5]=300.41384; pgvb10[6]=-1241.5067; pgvb10[7]=2476.2241;
    pgvb10[8]=-1304.3030; pgvb10[9]=-2149.6577; pgvb10[10]=4099.6917; pgvb10[11]=0;

    pgvls0[0]=0; pgvls0[1]=-66.176421; pgvls0[2]=2890.3688; pgvls0[3]=-62592.400;
    pgvls0[4]=691461.41; pgvls0[5]=-4096914.6; pgvls0[6]=14032093.; pgvls0[7]=-26827468.;
    pgvls0[8]=23511442.; pgvls0[9]=-14688461.; pgvls0[10]=0; pgvls0[11]=0;

    pgvt0[0]=-10.763625; pgvt0[1]=-0.46818029; pgvt0[2]=60.147739; pgvt0[3]=352.56941;
    pgvt0[4]=514.32170; pgvt0[5]=11637.302; pgvt0[6]=-44595.415; pgvt0[7]=69211.738;
    pgvt0[8]=-48127.668; pgvt0[9]=7051.4008; pgvt0[10]=0; pgvt0[11]=0;

    pgvso20[0]=0; pgvso20[1]=-0.62851020; pgvso20[2]=-76.290197; pgvso20[3]=-788.27581;
    pgvso20[4]=-6490.4798; pgvso20[5]=5473.4378; pgvso20[6]=-32941.912; pgvso20[7]=249491.32;
    pgvso20[8]=-16012.956; pgvso20[9]=0; pgvso20[10]=0; pgvso20[11]=0;

    if(lpot==0){
        double suma01=0;
        double suma11=0;
        double suma00=0;
        double suma10=0;
        double sumb01=0;
        double sumb11=0;
        double sumb00=0;
        double sumb10=0;
        for(int i=0; i<11; i++){
            suma01=suma01+pgva01[i]/pm1[i];
            suma11=suma11+pgva11[i]/pm1[i];
            suma00=suma00+pgva00[i]/pm0[i];
            suma10=suma10+pgva10[i]/pm0[i];
            sumb01=sumb01+pgvb01[i]/pm1[i];
            sumb11=sumb11+pgvb11[i]/pm1[i];
            sumb00=sumb00+pgvb00[i]/pm0[i];
            sumb10=sumb10+pgvb10[i]/pm0[i];
        }

        pgva01[11]=-suma01*pm1[11];
        pgva11[11]=-suma11*pm1[11];
        pgva00[11]=-suma00*pm0[11];
        pgva10[11]=-suma10*pm0[11];
        pgvb01[11]=-sumb01*pm1[11];
        pgvb11[11]=-sumb11*pm1[11];
        pgvb00[11]=-sumb00*pm0[11];
        pgvb10[11]=-sumb10*pm0[11];
        double sumt1=0;
        double sumt0=0;
        double sumls1=0;
        double sumls0=0;
        double tumt1=0;
        double tumt0=0;
        double tumls1=0;
        double tumls0=0;
        for(int i=0; i<10; i++){
            sumt1=sumt1+pgvt1[i]/pm1[i];
            sumt0=sumt0+pgvt0[i]/pm0[i];
            sumls1=sumls1+pgvls1[i]/pm1[i];
            sumls0=sumls0+pgvls0[i]/pm0[i];
            tumt1=tumt1+pgvt1[i]/pow(pm1[i],3);
            tumt0=tumt0+pgvt0[i]/pow(pm0[i],3);
            tumls1=tumls1+pgvls1[i]/pow(pm1[i],3);
            tumls0=tumls0+pgvls0[i]/pow(pm0[i],3);
        }

        pgvt1[10]=pow(pm1[10],3)*(pow(pm1[11],2)*tumt1-sumt1)
               /(pow(pm1[10],2)-pow(pm1[11],2));
        pgvt1[11]=pow(pm1[11],3)*(pow(pm1[10],2)*tumt1-sumt1)
               /(pow(pm1[11],2)-pow(pm1[10],2));
        pgvt0[10]=pow(pm0[10],3)*(pow(pm0[11],2)*tumt0-sumt0)
               /(pow(pm0[10],2)-pow(pm0[11],2));
        pgvt0[11]=pow(pm0[11],3)*(pow(pm0[10],2)*tumt0-sumt0)
               /(pow(pm0[11],2)-pow(pm0[10],2));
        pgvls1[10]=pow(pm1[10],3)*(pow(pm1[11],2)*tumls1-sumls1)
               /(pow(pm1[10],2)-pow(pm1[11],2));
        pgvls1[11]=pow(pm1[11],3)*(pow(pm1[10],2)*tumls1-sumls1)
               /(pow(pm1[11],2)-pow(pm1[10],2));
        pgvls0[10]=pow(pm0[10],3)*(pow(pm0[11],2)*tumls0-sumls0)
               /(pow(pm0[10],2)-pow(pm0[11],2));
        pgvls0[11]=pow(pm0[11],3)*(pow(pm0[10],2)*tumls0-sumls0)
               /(pow(pm0[11],2)-pow(pm0[10],2));


        double sumso21=0;
        double sumso20=0;
        double tumso21=0;
        double tumso20=0;
        double uumso21=0;
        double uumso20=0;

        for(int i=0; i<9; i++){
            sumso21=sumso21+pgvso21[i]/pm1[i];
            sumso20=sumso20+pgvso20[i]/pm0[i];
            tumso21=tumso21+pgvso21[i]/pow(pm1[i],3);
            tumso20=tumso20+pgvso20[i]/pow(pm0[i],3);
            uumso21=uumso21+pgvso21[i]/pow(pm1[i],5);
            uumso20=uumso20+pgvso20[i]/pow(pm0[i],5);
        }


        pgvso21[9]=pow(pm1[9],5)*(-pow(pm1[10],2)*pow(pm1[11],2)*uumso21
                             +(pow(pm1[10],2)+pow(pm1[11],2))*tumso21-sumso21)
                /((pow(pm1[11],2)-pow(pm1[9],2))*(pow(pm1[10],2)-pow(pm1[9],2)));
        pgvso21[10]=pow(pm1[10],5)*(-pow(pm1[11],2)*pow(pm1[9],2)*uumso21
                             +(pow(pm1[11],2)+pow(pm1[9],2))*tumso21-sumso21)
                /((pow(pm1[9],2)-pow(pm1[10],2))*(pow(pm1[11],2)-pow(pm1[10],2)));
        pgvso21[11]=pow(pm1[11],5)*(-pow(pm1[9],2)*pow(pm1[10],2)*uumso21
                             +(pow(pm1[9],2)+pow(pm1[10],2))*tumso21-sumso21)
                /((pow(pm1[10],2)-pow(pm1[11],2))*(pow(pm1[9],2)-pow(pm1[11],2)));
        pgvso20[9]=pow(pm0[9],5)*(-pow(pm0[10],2)*pow(pm0[11],2)*uumso20
                             +(pow(pm0[10],2)+pow(pm0[11],2))*tumso20-sumso20)
                /((pow(pm0[11],2)-pow(pm0[9],2))*(pow(pm0[10],2)-pow(pm0[9],2)));
        pgvso20[10]=pow(pm0[10],5)*(-pow(pm0[11],2)*pow(pm0[9],2)*uumso20
                             +(pow(pm0[11],2)+pow(pm0[9],2))*tumso20-sumso20)
                /((pow(pm0[9],2)-pow(pm0[10],2))*(pow(pm0[11],2)-pow(pm0[10],2)));
        pgvso20[11]=pow(pm0[11],5)*(-pow(pm0[9],2)*pow(pm0[10],2)*uumso20
                             +(pow(pm0[9],2)+pow(pm0[10],2))*tumso20-sumso20)
                /((pow(pm0[10],2)-pow(pm0[11],2))*(pow(pm0[9],2)-pow(pm0[11],2)));
    }

}

// -------------
// malfliet-tjon
// -------------
void DLM_StefanoPotentials::MalflietTjon(){
    const double& rr = *CurrentVerySmallRadius;
    p01=(7.39*exp(-3.11*rr)-2.64*exp(-1.55*rr))*hc197/rr;
    p10=(7.39*exp(-3.11*rr)-3.22*exp(-1.55*rr))*hc197/rr;
    p00=0;
    p11=0;
// ------------------
// V version of force
// ------------------

    vv[0]=0.5*(p01+p10);

}

void DLM_StefanoPotentials::ReidV8(){
    const double& rr = *CurrentVerySmallRadius;
    uVar=0.7;
    xVar=uVar*rr;
    y1Var=yc(xVar, xVar);
    y2Var=yc(2*xVar, xVar);
    y3Var=yc(3*xVar, xVar);
    y4Var=yc(4*xVar, xVar);
    y6Var=yc(6*xVar, xVar);
    y7Var=yc(7*xVar, xVar);
    yrVar=yt(xVar, xVar)-(12./xVar+3./pow(xVar,2))*y4Var;
    if (*CurrentRadius<=10*vsmall) yrVar=23.5/xVar;
    hrVar=10.463;
    vv[0]=-19.874*y2Var+135.21*y3Var-1432.3*y4Var+4196.4*y6Var+1215.8*y7Var;
    vv[1]=19.874*y2Var-135.21*y3Var+319.52*y4Var-1082.3*y6Var+405.3*y7Var;
    vv[2]=46.241*y2Var-135.21*y3Var-64.78*y4Var+1398.8*y6Var-1215.8*y7Var;
    vv[3]=(hrVar/3.)*y1Var-46.241*y2Var+135.21*y3Var+244.06*y4Var-360.76*y6Var-405.3*y7Var;
    vv[4]=-26.194*y3Var+87.943*y4Var-418.38*y6Var;
    vv[5]=(hrVar/3.)*yrVar-8.731*y3Var-87.943*y4Var+418.38*y6Var;
    vv[6]=177.23*y4Var-2233.9*y6Var;
    vv[7]=-177.23*y4Var+159.75*y6Var;
    vw[1]=(hrVar/3.)*y1Var;
    vw[2]=(hrVar/3.)*yrVar;

//if(rr>=0.9 && rr<0.91){
    //for(int i=0; i<18; i++){
    //    printf("rr=%f; vv[%i]=%f\n", rr, i, vv[i]);
    //}
//}
}

void DLM_StefanoPotentials::UrbanaV14(){
    const double& rr = *CurrentRadius;
    uVar=0.7;
    cpi=2;
    xVar=uVar*rr;
    if (rr<=small){
        ypi=cpi*rr/uVar;
        tpi=3*pow(cpi, 2)*rr/pow(uVar, 3);
    }
    else{
        rcut=1-exp(-cpi*rr*rr);
        ypi=yc(xVar, xVar)*rcut;
        tpi=yt(xVar, xVar)*pow(rcut, 2);
    }
    tpi2=tpi*tpi*tnr;

    ypi=10.463*ypi/3.;
    tpi=10.463*tpi/3.;
    ws=1/(1+exp((rr-.5)/.2));
    wp=1/(1+exp((rr-.36)/.17));

    p11=  -4.32  *tpi2+2145.*ws+  ypi;
    pt1=  -0.18  *tpi2         +  tpi;
    pls1=             -2200.*wp;
    pl211=            -  20.*ws;
    pls21=             147.5*ws;
// following is old urbana v14 (incorrect deuteron)
//     p10=  -6.8009*tpi2+2400.*ws-3*ypi
// following is new urbana v14 (correct deuteron)
    p10=  -6.7983*tpi2+2400.*ws-3*ypi;
    pt0=   0.75  *tpi2         -3*tpi;
    pls0=                80.*ws;
    pl210=              380.*ws;
    pls20=-0.2   *tpi2- 230.*ws;
    p01=  -6.255 *tpi2+2000.*ws-3*ypi;
    pl201=               49.*ws;
    p00= -13.2   *tpi2+8700.*ws+9*ypi;
    pl200= 0.6   *tpi2- 500.*ws;
    vv[0]=.0625*(9*p11+3*p10+3*p01+p00);
    vv[1]=.0625*(3*p11-3*p10  +p01-p00);
    vv[2]=.0625*(3*p11  +p10-3*p01-p00);
    vv[3]=.0625*(  p11  -p10  -p01+p00);
    vv[4]=.25*(3*pt1+pt0);
    vv[5]=.25*(  pt1-pt0);
    vv[6]=.25*(3*pls1+pls0);
    vv[7]=.25*(  pls1-pls0);
    vv[8]= .0625*(9*pl211+3*pl210+3*pl201+pl200);
    vv[9]=.0625*(3*pl211-3*pl210+  pl201-pl200);
    vv[10]=.0625*(3*pl211+  pl210-3*pl201-pl200);
    vv[11]=.0625*(  pl211-  pl210-  pl201+pl200);
    vv[12]=.25*(3*pls21+pls20);
    vv[13]=.25*(  pls21-pls20);
    vw[1]=ypi;
    vw[2]=tpi;
}

void DLM_StefanoPotentials::ArgonneV14(){
    const double& rr = *CurrentRadius;
    uVar=138.03/197.33;
    u1Var=pimass1*uVar;
    u2Var=pimass2*uVar;
    u3Var=pimass3*uVar;
    cpi=2;
    x1Var=u1Var*rr;
    x2Var=u2Var*rr;
    x3Var=u3Var*rr;
    if (rr<=small){
        ypi=cpi*rr/u1Var;
        tpi=3*pow(cpi,2)*rr/pow(u1Var,3);
        ypib=cpi*rr/u2Var;
        tpib=3*pow(cpi,2)*rr/pow(u2Var,3);
        ypid=cpi*rr/u3Var;
        tpid=3*pow(cpi,2)*rr/pow(u3Var,3);
    }
    else{
        rcut=1-exp(-cpi*rr*rr);
        ypi=rcut*exp(-x1Var)/x1Var;
        tpi=rcut*(1+3/x1Var+3/pow(x1Var,2))*ypi;
        ypib=rcut*exp(-x2Var)/x2Var;
        tpib=rcut*(1+3/x2Var+3/pow(x2Var,2))*ypib;
        ypid=rcut*exp(-x3Var)/x3Var;
        tpid=rcut*(1+3/x3Var+3/pow(x3Var,2))*ypid;
    }
    pifac=pow(pimass1,3);
    ypi=pifac*3.72681*ypi;
    tpi=pifac*3.72681*tpi;
    tpi2=ftpec*pow(pow(pimass2,3)*tpib, 2);
    tpi3=ftpec*pow(pow(pimass3,3)*tpid, 2);
    rws=wsrange*0.5;
    aws=wsrange*0.2;
    ws=1./(1+exp((rr-rws)/aws));
    p11=  -2.63  *tpi2+1179.*ws+  ypi;
    pt1=  -0.91  *tpi2+ 406.*ws+  tpi;
    pls1=  0.61  *tpi3- 879.*ws;
    pl211=-0.12  *tpi3-   2.*ws;
    pls21=-0.54  *tpi3+ 536.*ws;
    p10=  -6.5572*tpi2+2700.*ws-3*ypi;
    pt0=   2.10  *tpi2- 783.*ws-3*tpi;
    pls0=  0.42  *tpi3- 242.*ws;
    pl210= 0.48  *tpi3+ 110.*ws;
    pls20=-0.55  *tpi3-  44.*ws;
    p01=  -8.1188*tpi2+2800.*ws-3*ypi;
    pl201= 0.05  *tpi3+  63.*ws;
    p00=  -9.12  *tpi2+5874.*ws+9*ypi;
    pl200= 0.62  *tpi3- 363.*ws;
// -----------
    if (lpot!=6){
        pl2av=0.;
// fix average of D waves
//       if (lpot.eq.5) pl2av=.5*(pl210+pl201)+pls20/3
// fix 1D2 partial wave
        if (lpot==5) pl2av=2.*pl201/3.;
        p00=p00+2.*(pl200-pl2av);
        pls0=pls0-2.*(pl210-pl2av)-3.*pls20;
        p11=p11+2.*(pl211-pl2av)+4.*pls21/3.;
        pt1=pt1-5.*pls21/12.;
        pls1=pls1-0.5*pls21;
    }
// -----------
    vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
    vv[1]=0.0625*(3.*p11-3.*p10  +p01-p00);
    vv[2]=0.0625*(3.*p11  +p10-3.*p01-p00);
    vv[3]=0.0625*(  p11  -p10  -p01+p00);
    vv[4]=0.25*(3.*pt1+pt0);
    vv[5]=0.25*(  pt1-pt0);
    vv[6]=0.25*(3.*pls1+pls0);
    vv[7]=0.25*(  pls1-pls0);
    if (lpot==5){
        vv[8]=pl2av;
    }
    else{
        vv[8]= 0.0625*(9.*pl211+3.*pl210+3.*pl201+pl200);
        vv[9]=0.0625*(3.*pl211-3.*pl210+  pl201-pl200);
        vv[10]=0.0625*(3.*pl211+  pl210-3.*pl201-pl200);
        vv[11]=0.0625*(  pl211-  pl210-  pl201+pl200);
        vv[12]=0.25*(3*pls21+pls20);
        vv[13]=0.25*(  pls21-pls20);
    }

    vw[1]=ypi;
    vw[2]=tpi;
}

//---------------------------
//argonne v18 and derivatives
//lpot = 7 -> v18-csbl
//       8 -> v18-csbs
//       9 -> v18
//      10 -> v8'
//      11 -> v6'
//      12 -> v4'
//      13 -> v2'
//      14 -> v1'
//      15 -> vx'
//      21 -> v18 v1.9
//      22 -> v8' v1.9
//      23 -> v18 v1.7
//      24 -> v8' v1.7
//---------------------------
void DLM_StefanoPotentials::ArgonneV18(){
    const double& rr = *CurrentRadius;
    mpi0=134.9739;
    mpic=139.5675;
    mpi=(mpi0+2.*mpic)/3.;
    mpi0=pimass1*mpi0;
    mpic=pimass1*mpic;
    mpis=pimass2*mpi;
    mpiq=pimass3*mpi;
    mu=mpi/hc197;
    mu0=mpi0/hc197;
    muc=mpic/hc197;
    mus=mpis/hc197;
    muq=mpiq/hc197;
    fsq=0.075;
    if (lpot<=15){
        cpi=2.1;
        rws=wsrange*.5;
        aws=0.2*wsrange;
//      aiws=5./wsrange
    }
    else if (lpot==21 || lpot==22){
        cpi=1.9;
        rws=wsrange*0.525;
        aws=0.21*wsrange;
    }
    else if (lpot==23 || lpot==24){
        cpi=1.7;
        rws=wsrange*0.55;
        aws=0.22*wsrange;
    }
    aiws=1./aws;
    xVar=mu*rr;
    xsVar=mus*rr;
    xqVar=muq*rr;
    x0Var=mu0*rr;
    xcVar=muc*rr;
    if (rr<=small){
        tpis=3*pow(cpi,2)*rr/pow(mus,3);
        tpiq=3*pow(cpi,2)*rr/pow(muq,3);
        ypi0=pow(mpi0/mpic,2)*(mpi0/3)*cpi*rr/mu0;
        ypic=(mpic/3)*cpi*rr/muc;
        tpi0=3*cpi*ypi0/pow(mu0,2);
        tpic=3*cpi*ypic/pow(muc,2);
    }
    else{
        rcut=1-exp(-cpi*rr*rr);
        tpis=pow(rcut,2)*(1+3/xsVar+3/pow(xsVar,2))*exp(-xsVar)/xsVar;
        tpiq=pow(rcut,2)*(1+3/xqVar+3/pow(xqVar,2))*exp(-xqVar)/xqVar;
        ypi0=pow(mpi0/mpic,2)*(mpi0/3)*exp(-x0Var)*rcut/x0Var;
        ypic=(mpic/3)*exp(-xcVar)*rcut/xcVar;
        tpi0=(1+(3+3/x0Var)/x0Var)*ypi0*rcut;
        tpic=(1+(3+3/xcVar)/xcVar)*ypic*rcut;
    }
    pifac=pow(pimass1,2);
    ypi0=pifac*fsq*ypi0;
    ypic=pifac*fsq*ypic;
    tpi0=pifac*fsq*tpi0;
    tpic=pifac*fsq*tpic;
    tpi2=ftpec*pow(pow(pimass2,3)*tpis,2);
    tpi3=ftpec*pow(pow(pimass3,3)*tpiq,2);
    ews=exp((rr-rws)*aiws);
    ws=1./(1.+ews);
    ews0=exp(-rws*aiws);
    ws0=1./(1.+ews0);
    fws=aiws*ews0*ws0;
    wsp=(1.+fws*rr)*ws;
    wsz=-aiws*ews*pow(ws,2);
    wspp=fws*ws+(1.+fws*rr)*wsz;
    wszz=-2.*aiws*ews*ws*wsz-pow(aiws,2)*ews*pow(ws,2);
    wspdp=2.*fws*wsz+(1.+fws*rr)*wszz;
    wspds=wspdp+2.*wspp/rr;
    wsx=ws*xVar;
    wsx2=wsx*xVar;
    dypi00=pow(mpi0/mpic,2)*(mpi0/3.)*cpi/mu0;
    dypic0=(mpic/3.)*cpi/muc;
    ypi0p=ypi0-fsq*dypi00*ws*rr/ws0;
    ypicp=ypic-fsq*dypic0*ws*rr/ws0;
    ypi=(ypi0+2.*ypic)/3.;
    tpi=(tpi0+2.*tpic)/3.;
    ypibar=(ypi0-ypic)/3.;
    tpibar=(tpi0-tpic)/3.;
// final version 11/1/93
// nn potential added 11/15/93
// absolutely final version 3/29/94
// totally absolutely final version 6/28/94
    if (lpot<=15){
        p11pp=  -7.62701*tpi2+1815.4920*wsp+1847.8059*wsx2+ypi0p;
        p11np=  -7.62701*tpi2+1813.5315*wsp+1847.8059*wsx2-ypi0p+2*ypicp;
        p11nn=  -7.62701*tpi2+1811.5710*wsp+1847.8059*wsx2+ypi0p;
        pt1pp=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0;
        pt1np=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2-tpi0+2*tpic;
        pt1nn=   1.07985*tpi2 -190.0949*wsx -811.2040*wsx2+tpi0;
        p01pp= -11.27028*tpi2+3346.6874*wsp-3*ypi0p;
        p01np= -10.66788*tpi2+3126.5542*wsp-3*(-ypi0p+2*ypicp);
        p01nn= -11.27028*tpi2+3342.7664*wsp-3*ypi0p;
// location for CD-CSB test version :
// keeps p11, pt1, pt1cs, p01, p01cd the same
// makes p11cd = p01cd/9,
//       pt1cd = 5*p11cd/7 for tpi2 ; wsx adjusted to preserve phase shift
// while p11cs, p01cs are set with extra class iii-iv terms below
        if (lpot==7 || lpot==8){
            p11pp=  -7.64930*tpi2+1821.6119*wsp+1847.8059*wsx2+ypi0p;
            p11np=  -7.58243*tpi2+1797.3707*wsp+1847.8059*wsx2-ypi0p+2.*ypicp;
            p11nn=  -7.64930*tpi2+1821.6119*wsp+1847.8059*wsx2+ypi0p;
            pt1pp=   1.06391*tpi2 -178.6199*wsx -811.2040*wsx2+tpi0;
            pt1np=   1.11173*tpi2 -213.0449*wsx -811.2040*wsx2-tpi0+2.*tpic;
            pt1nn=   1.06391*tpi2 -178.6199*wsx -811.2040*wsx2+tpi0;
            p01pp= -11.27028*tpi2+3344.7269*wsp-3.*ypi0p;
            p01np= -10.66788*tpi2+3126.5542*wsp-3.*(-ypi0p+2.*ypicp);
            p01nn= -11.27028*tpi2+3344.7269*wsp-3.*ypi0p;
        }
// CD test complete
        pls1=    -0.62697*tpi3 -570.5571*wsp +819.1222*wsx2;
        pl211=    0.06709*tpi3 +342.0669*wsp -615.2339*wsx2;
        pls21=    0.74129*tpi3   +9.3418*wsp -376.4384*wsx2;
        p10=    -8.62770*tpi2+2605.2682*wsp +441.9733*wsx2-ypi0p-2.*ypicp;
        pt0=    1.485601*tpi2-1126.8359*wsx +370.1324*wsx2-tpi0-2.*tpic;
        pls0=     0.10180*tpi3  +86.0658*wsp -356.5175*wsx2;
        pl210=   -0.13201*tpi3 +253.4350*wsp   -1.0076*wsx2;
        pls20=    0.07357*tpi3 -217.5791*wsp  +18.3935*wsx2;
        pl201=    0.12472*tpi3  +16.7780*wsp;
        p00=    -2.09971*tpi2+1204.4301*wsp-3*(-ypi0p-2.*ypicp);
        pl200=   -0.31452*tpi3 +217.4559*wsp;
// location for v1.9
    }
    else if (lpot==21 || lpot==22){
        p11pp=  -8.25333*tpi2+1629.4894*wsp+1007.3079*wsx2+ypi0p;
        p11np=  -8.25333*tpi2+1627.8623*wsp+1007.3079*wsx2-ypi0p+2.*ypicp;
        p11nn=  -8.25333*tpi2+1626.2352*wsp+1007.3079*wsx2+ypi0p;
        pt1pp=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2+tpi0;
        pt1np=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2-tpi0+2.*tpic;
        pt1nn=   1.22738*tpi2 -331.5020*wsx -415.4240*wsx2+tpi0;
        pls1=   -1.24596*tpi2 -438.1866*wsp +881.8829*wsx2;
        pl211=   0.17268*tpi2 +210.3707*wsp -418.9703*wsx2;
        pls21=   0.68968*tpi2  -44.1763*wsp -154.7568*wsx2;
        p10=   -10.62968*tpi2+2297.1952*wsp +503.6560*wsx2-ypi0p-2.*ypicp;
        pt0=     1.43163*tpi2 -932.7628*wsx +415.7518*wsx2-tpi0-2.*tpic;
        pls0=    0.31692*tpi2   +5.9540*wsp -261.4438*wsx2;
        pl210=   0.20369*tpi2 +164.8268*wsp -133.2324*wsx2;
        pls20=  -0.08370*tpi2 -162.9074*wsp  +92.1321*wsx2;
        p01pp= -11.24918*tpi2+2446.4156*wsp -278.5780*wsx2-3.*ypi0p;
        p01np= -10.57598*tpi2+2273.0877*wsp -289.7548*wsx2
                                               -3.*(-ypi0p+2.*ypicp);
        p01nn= -11.24918*tpi2+2443.1614*wsp -278.5780*wsx2-3.*ypi0p;
        pl201=   0.12622*tpi2  +23.4445*wsp  -13.5987*wsx2;
        p00=    -2.14060*tpi2+1000.2218*wsp -167.8362*wsx2
                                               -3.*(-ypi0p-2.*ypicp);
        pl200=  -0.30287*tpi2 +178.4343*wsp  -40.3099*wsx2;
// v1.9 complete
// location for v1.7
    }
    else if (lpot==23 || lpot==24){
        p11pp=  -9.10461*tpi2+1383.9447*wsp +635.3480*wsx2+ypi0p;
        p11np=  -9.10461*tpi2+1382.5743*wsp +635.3480*wsx2-ypi0p+2*ypicp;
        p11nn=  -9.10461*tpi2+1381.2039*wsp +635.3480*wsx2+ypi0p;
        pt1pp=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2+tpi0;
        pt1np=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2-tpi0+2*tpic;
        pt1nn=   1.41302*tpi2 -353.5571*wsx -220.1474*wsx2+tpi0;
        pls1=   -2.32067*tpi2 -294.5597*wsp +905.8846*wsx2;
        pl211=   0.31011*tpi2 +124.9726*wsp -295.2155*wsx2;
        pls21=   0.56878*tpi2  -53.4883*wsp  -35.2622*wsx2;
        p10=   -12.16378*tpi2+1841.4730*wsp +490.2208*wsx2-ypi0p-2*ypicp;
        pt0=     1.38999*tpi2 -748.1056*wsx +373.6946*wsx2-tpi0-2*tpic;
        pls0=    0.38755*tpi2  -17.6948*wsp -168.8546*wsx2;
        pl210=   0.61367*tpi2 +108.7064*wsp -206.9220*wsx2;
        pls20=  -0.35600*tpi2 -118.0195*wsp +148.1324*wsx2;
        p01pp= -11.28159*tpi2+1735.7847*wsp -290.3768*wsx2-3*ypi0p;
        p01np= -11.62496*tpi2+1724.8767*wsp -116.8103*wsx2
                                               -3*(-ypi0p+2*ypicp);
        p01nn= -11.28159*tpi2+1733.0438*wsp -290.3768*wsx2-3*ypi0p;
        pl201=   0.12387*tpi2  +24.9311*wsp  -14.7805*wsx2;
        p00=    -1.94411*tpi2 +808.3897*wsp -277.7773*wsx2
                                               -3*(-ypi0p-2*ypicp);
        pl200=  -0.31130*tpi2 +139.5016*wsp  -40.4310*wsx2;
// v1.7 complete
    }
    p11=(p11pp+p11nn+p11np)/3.;
    p11cd=(0.5*(p11pp+p11nn)-p11np)/6.;
    p11cs=(p11pp-p11nn)/4.;
    pt1=(pt1pp+pt1nn+pt1np)/3.;
    pt1cd=(0.5*(pt1pp+pt1nn)-pt1np)/6.;
    pt1cs=(pt1pp-pt1nn)/4.;
    p01=(p01pp+p01nn+p01np)/3.;
    p01cd=(0.5*(p01pp+p01nn)-p01np)/6.;
    p01cs=(p01pp-p01nn)/4.;
// -----------
    if (lpot>=10 && lpot!=21 && lpot!=23){
        pl2av=0.;
        p00=p00+2*(pl200-pl2av);
        pls0=pls0-2*(pl210-pl2av)-3*pls20;
        p11=p11+2*(pl211-pl2av)+4*pls21/3;
        pt1=pt1-5*pls21/12;
        pls1=pls1-.5*pls21;
// fix deuteron in v6' case
        if (lpot>=11 && lpot<=15) p10=p10-0.3*pls0;
// fix deuteron in v4' case
        if (lpot>=12 && lpot<=15) p10=p10+0.8735*pt0;
// project only vc and vt in v2' case
        if (lpot==13){
            vv[0]=0.25*(3.*p01+p10);
            vv[1]=0.25*(  p01-p10);
            return;
        }
// average 1s & 3s in v1' case
        else if (lpot==14){
            vv[0]=0.5*(p01+p10);
            return;
        }
// combination for vx' case
        else if (lpot==15){
            vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
            vv[1]=0.0125*(9.*p11-5.*p10-5.*p01+p00);
            vv[2]=vv[1];
            vv[3]=vv[1];
            return;
        }
    }
// -----------
    vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
    vv[1]=0.0625*(3.*p11-3.*p10  +p01-p00);
    vv[2]=0.0625*(3.*p11  +p10-3.*p01-p00);
    vv[3]=0.0625*(  p11  -p10  -p01+p00);
    if (lpot==12) return;
    vv[4]=0.25*(3.*pt1+pt0);
    vv[5]=0.25*(  pt1-pt0);
    if(lpot==11) return;
    vv[6]=0.25*(3.*pls1+pls0);
    vv[7]=0.25*(  pls1-pls0);
    if (lpot!=10 && lpot!=22 && lpot!=24){
        vv[8]= .0625*(9*pl211+3*pl210+3*pl201+pl200);
        vv[9]=.0625*(3*pl211-3*pl210+  pl201-pl200);
        vv[10]=.0625*(3*pl211+  pl210-3*pl201-pl200);
        vv[11]=.0625*(  pl211-  pl210-  pl201+pl200);
        vv[12]=.25*(3*pls21+pls20);
        vv[13]=.25*(  pls21-pls20);
        vv[14]=.25*(3*p11cd+p01cd);
        vv[15]=.25*(  p11cd-p01cd);
        vv[16]=pt1cd;
        vv[17]=p01cs;
    }
//for(int i=0; i<18; i++){
//    printf("rr=%f; vv[%i]=%f\n", rr, i, vv[i]);
//}
    vw[1]=ypi;
    vw[2]=tpi;
// delta function added to ypi!
//     xlm=900./mpi
//     xlm2=xlm*xlm
//     rmu=mu*rr
//     rlm=xlm*rmu
//     ermu=exp(-rmu)/rmu
//     erlm=exp(-rlm)/rmu
//     vw(2)=fsq*(mpi/3)
//    &     *(ermu-xlm2*erlm-.5*xlm*(xlm2-1.)*(1.-2./rlm)*erlm*rmu)
//     vw(3)=fsq*(mpi/3)
//    &     *((1.+3./rmu+3./rmu**2)*ermu-xlm2*(1.+3./rlm+3./rlm**2)*erlm
//    &       -.5*xlm*(xlm2-1.)*(1.+1./rlm)*erlm*rmu)
//     b=4.27
//     br=b*rr
//     fdelta=b**3*(1+br+br**2/3)*exp(-br)/16
//     vw(2)=ypi-fsq*(mpi/3)*fdelta/mu**3
// delta function added to ypi!
    vw[3]=ypibar;
    vw[4]=tpibar;
// extra class iii-iv csb
    if (lpot==7 || lpot==8){
// -----------------------
// following is for case 1
        if (lpot==7){
            vcsb=1.50875;
            krho=6.1;
            komg=0.14;
// following is for case 2
        }
        else if (lpot==8){
            vcsb=1.11160;
            krho=3.7;
            komg=-0.12;
        }
// -----------------------
        vv[17]=vcsb*wsp;
        if (rr>=vsmall){
            deltam=(mn-mp)/(mn+mp);
            crho=2.44;
            grho=0.55;
            fmsqi=vcsb*pow(hc197/(mn+mp),2);
            vw[5]=fmsqi*(2+komg+krho)*wspds;
            vw[6]=2.*fmsqi*(1.+komg)*(1.+krho)*wspds/3.;
            vw[7]=-fmsqi*(1.+komg)*(1.+krho)*(wspds-3.*wspp/rr)/3.;
            vw[8]=4.*fmsqi*(2.+komg+krho)*wspp/rr;
            vw[9]=4.*fmsqi*(krho-komg)*wspp/rr;
            mrho=770.;
            murho=mrho/hc197;
            xrho=murho*rr;
            rcut=1-exp(-cpi*rr*rr);
            zpic=mpic*((1+1/xcVar)/xcVar)*exp(-xcVar)*pow(rcut,1.5)/xcVar;
            zcut=1-exp(-crho*rr*rr);
            zrho=mrho*((1+1/xrho)/xrho)*exp(-xrho)*pow(zcut,1.5)/xrho;
            vw[10]=2.*fsq*deltam*zpic;
            vw[11]=2.*grho*deltam*pow(mrho/(mn+mp),2)*zrho;
            vw[12]=2.*pow(1.+krho,2)*vw[11];
        }
    }
}

// ----------------------------------------
// argonne v18p  (p**2 terms)
// argonne v18pq (p**2 terms & l**2 tensor)
// lpot = 16 -> v18p
//        17 -> v18pq
// ----------------------------------------
void DLM_StefanoPotentials::ArgonneV18p(){
    const double& rr = *CurrentRadius;
    mpi0=134.9739;
    mpic=139.5675;
    mpi=(mpi0+2*mpic)/3.;
    mu0=mpi0/hc197;
    muc=mpic/hc197;
    mu=mpi/hc197;
    fsq=0.075;
    cpi=2.1;
    rws=0.5;
    aiws=5.;
    xVar=mu*rr;
    x0Var=mu0*rr;
    xcVar=muc*rr;
    if (rr<=small){
        tpi=3.*pow(cpi,2)*rr/pow(mu,3);
        ypi0=pow(mpi0/mpic,2)*(mpi0/3.)*cpi*rr/mu0;
        ypic=(mpic/3.)*cpi*rr/muc;
        tpi0=3.*cpi*ypi0/pow(mu0,2);
        tpic=3.*cpi*ypic/pow(muc,2);
    }
    else{
        expcut=exp(-cpi*rr*rr);
        rcut=1.-expcut;
        rcutp=2.*cpi*rr*expcut;
        rcutdp=2.*cpi*expcut*(1.-2.*cpi*pow(rr,2));
        ypi0=pow(mpi0/mpic,2)*(mpi0/3.)*exp(-x0Var)*rcut/x0Var;
        ypic=(mpic/3)*exp(-xcVar)*rcut/xcVar;
        tpi0=(1.+(3.+3./x0Var)/x0Var)*ypi0*rcut;
        tpic=(1.+(3.+3./xcVar)/xcVar)*ypic*rcut;
        yVar=exp(-xVar)/xVar;
        ypVar=-mu*(1.+1./xVar)*yVar;
        ydpVar=pow(mu,2)*(1.+2.*(1.+1./xVar)/xVar)*yVar;
        tVar=1.+3.*(1+1/xVar)/xVar;
        tpVar=-mu*3.*(1+2/xVar)/pow(xVar,2);
        tdpVar=pow(mu,2)*6.*(1+3/xVar)/pow(xVar,3);
        tpi=tVar*yVar*pow(rcut,2);
        tpip=(tpVar*yVar+tVar*ypVar)*pow(rcut,2)+2.*tVar*yVar*rcut*rcutp;
        tpidp=(tdpVar*yVar+2.*tpVar*ypVar+tVar*ydpVar)*pow(rcut,2)+4.*(tpVar*yVar+tVar*ypVar)*rcut*rcutp
            +2.*tVar*yVar*(rcut*rcutdp+pow(rcutp,2));
    }
    ypi0=fsq*ypi0;
    ypic=fsq*ypic;
    tpi0=fsq*tpi0;
    tpic=fsq*tpic;
    tpi2=tpi*tpi;
    tpi2p=2.*tpi*tpip;
    tpi2dp=2.*(tpi*tpidp+pow(tpip,2));
    expws=exp((rr-rws)*aiws);
    ws=1./(1.+expws);
    wsp=-aiws*expws*pow(ws,2);
    wsdp=aiws*wsp+2.*pow(wsp,2)/ws;
    ws0=1./(1.+exp(-rws*aiws));
    wsp0=-aiws*exp(-rws*aiws)*pow(ws0,2);
    wsm=ws*(1.-wsp0*xVar/(mu*ws0));
    wsmp=wsp-wsp0*(xVar*wsp+mu*ws)/(mu*ws0);
    wsmdp=wsdp-wsp0*(xVar*wsdp+2.*mu*wsp)/(mu*ws0);
    wsx=ws*xVar;
    wsx2=wsx*xVar;
    wsx2p=pow(xVar,2)*wsp+2.*mu*xVar*ws;
    wsx2dp=pow(xVar,2)*wsdp+4.*mu*xVar*wsp+2.*pow(mu,2)*ws;
    dypi00=pow(mpi0/mpic,2)*(mpi0/3.)*cpi/mu0;
    dypic0=(mpic/3.)*cpi/muc;
    ypi0p=ypi0-fsq*dypi00*ws*rr/ws0;
    ypicp=ypic-fsq*dypic0*ws*rr/ws0;
    ypi=(ypi0+2.*ypic)/3.;
    tpi=(tpi0+2.*tpic)/3.;
    ypibar=(ypi0-ypic)/3.;
    tpibar=(tpi0-tpic)/3.;
// ----------------------------------------
// first version 2002.07.18
// second version 2003.01.25 (deuteron bug fixed)
// plain p**2 version 2003.02.13
// inconsistent pl2p0x and pl2dp0x corrected 2003.08.04
// ----------------------------------------
    if (lpot==16){
        p11pp= -11.059567*tpi2+2828.6889*wsm+3178.7405*wsx2+ypi0p;
        p11np= -11.059567*tpi2+2826.8983*wsm+3178.7405*wsx2-ypi0p+2.*ypicp;
        p11nn= -11.059567*tpi2+2825.1077*wsm+3178.7405*wsx2+ypi0p;
        pt1pp=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2+tpi0;
        pt1np=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2-tpi0+2.*tpic;
        pt1nn=   1.514600*tpi2 -376.1341*wsx -955.2246*wsx2+tpi0;
        pls1=   -0.732326*tpi2 -576.5707*wsm +925.8808*wsx2;
        pl211=  -0.420808*tpi2 +209.2711*wsm  -51.6679*wsx2;
        pls21=   1.023793*tpi2  -18.8040*wsm -544.7574*wsx2;
        p10=    -8.061451*tpi2+1938.3975*wsm+1895.2537*wsx2-ypi0p-2.*ypicp;
        pt0=     1.072145*tpi2-1220.5266*wsx +872.6085*wsx2-tpi0-2.*tpic;
        pls0=    0.328041*tpi2  -74.7348*wsm -271.1134*wsx2;
        pl210=   0.226986*tpi2  -19.8530*wsm +104.6043*wsx2;
        pls20=  -0.091272*tpi2  -67.3591*wsm  -81.6822*wsx2;
    }
// ----------------------------------------
    else if (lpot==17){
        p11pp=  -9.882847*tpi2+2589.7742*wsm+2952.3910*wsx2+ypi0p;
        p11np=  -9.882847*tpi2+2587.9836*wsm+2952.3910*wsx2-ypi0p+2.*ypicp;
        p11nn=  -9.882847*tpi2+2586.1930*wsm+2952.3910*wsx2+ypi0p;
        pt1pp=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2+tpi0;
        pt1np=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2-tpi0+2.*tpic;
        pt1nn=   1.420069*tpi2 -453.5357*wsx -837.3820*wsx2+tpi0;
        pls1=   -1.749197*tpi2 -493.8470*wsm+1533.0637*wsx2;
        pl211=  -0.008159*tpi2 +132.7694*wsm -169.8510*wsx2;
        pls21=   0.135181*tpi2  -17.7975*wsm  -46.2542*wsx2;
        p10=    -8.351808*tpi2+2325.5929*wsm +957.8091*wsx2-ypi0p-2.*ypicp;
        pt0=     1.327862*tpi2-1170.8528*wsx +580.5596*wsx2-tpi0-2.*tpic;
        pls0=    0.060223*tpi2  +58.3208*wsm -126.0235*wsx2;
        pl210=  -0.023577*tpi2   +1.8164*wsm +127.2921*wsx2;
        pls20=   0.000589*tpi2  -25.1123*wsm   -4.6897*wsx2;
    }
// ----------------------------------------
    p01pp= -10.518030*tpi2+2836.0715*wsm +651.1945*wsx2-3*ypi0p;
    p01np= -10.812190*tpi2+2816.4190*wsm+1002.5300*wsx2
                                                    -3*(-ypi0p+2*ypicp);
    p01nn= -10.518030*tpi2+2832.4903*wsm +651.1945*wsx2-3*ypi0p;
    pl201=   0.134747*tpi2   -9.4691*wsm;
    p00=    -4.739629*tpi2+1121.2225*wsm+2764.3395*wsx2
                                                    -3*(-ypi0p-2*ypicp);
    pl200=  -0.227084*tpi2 +166.5629*wsm;
// ----------------------------------------
    p11=(p11pp+p11nn+p11np)/3.;
    p11cd=(0.5*(p11pp+p11nn)-p11np)/6.;
    p11cs=(p11pp-p11nn)/4.;
    pt1=(pt1pp+pt1nn+pt1np)/3.;
    pt1cd=(0.5*(pt1pp+pt1nn)-pt1np)/6.;
    pt1cs=(pt1pp-pt1nn)/4.;
    p01=(p01pp+p01nn+p01np)/3.;
    p01cd=(0.5*(p01pp+p01nn)-p01np)/6.;
    p01cs=(p01pp-p01nn)/4.;
    vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
    vv[1]=0.0625*(3.*p11-3.*p10  +p01-p00);
    vv[2]=0.0625*(3.*p11  +p10-3.*p01-p00);
    vv[3]=0.0625*(  p11  -p10  -p01+p00);
    vv[4]=0.25*(3.*pt1+pt0);
    vv[5]=0.25*(  pt1-pt0);
    vv[6]=0.25*(3.*pls1+pls0);
    vv[7]=0.25*(  pls1-pls0);
    vv[8]= 0.0625*(9.*pl211+3.*pl210+3.*pl201+pl200);
    vv[9]=0.0625*(3.*pl211-3.*pl210+  pl201-pl200);
    vv[10]=0.0625*(3.*pl211+  pl210-3.*pl201-pl200);
    vv[11]=0.0625*(  pl211-  pl210-  pl201+pl200);
    vv[12]=0.25*(3.*pls21+pls20);
    vv[13]=0.25*(  pls21-pls20);
    vv[14]=0.25*(3.*p11cd+p01cd);
    vv[15]=0.25*(  p11cd-p01cd);
    vv[16]=pt1cd;
    vv[17]=p01cs;
    vw[1]=ypi;
    vw[2]=tpi;
    vw[3]=ypibar;
    vw[4]=tpibar;
// -----------------------
// vpsq & derivative terms
// -----------------------
    if (lpot==16){
        pl2p11=   -0.420808*tpi2p  +209.2711*wsmp   -51.6679*wsx2p;
        pl2p10=    0.226986*tpi2p   -19.8530*wsmp  +104.6043*wsx2p;
        pl2dp11=  -0.420808*tpi2dp +209.2711*wsmdp  -51.6679*wsx2dp;
        pl2dp10=   0.226986*tpi2dp  -19.8530*wsmdp +104.6043*wsx2dp;
    }
// -----------------------
    else if (lpot==17){
        pl2p11=   -0.008159*tpi2p  +132.7694*wsmp  -169.8510*wsx2p;
        pl2p10=   -0.023577*tpi2p    +1.8164*wsmp  +127.2921*wsx2p;
        pl2dp11=  -0.008159*tpi2dp +132.7694*wsmdp -169.8510*wsx2dp;
        pl2dp10=  -0.023577*tpi2dp   +1.8164*wsmdp +127.2921*wsx2dp;
    }
// -----------------------
    pl2p01=    0.134747*tpi2p    -9.4691*wsmp;
    pl2p00=   -0.227084*tpi2p  +166.5629*wsmp;
    pl2dp01=   0.134747*tpi2dp   -9.4691*wsmdp;
    pl2dp00=  -0.227084*tpi2dp +166.5629*wsmdp;
// ----------------------------------------
    psq11=0.5*pow(rr,2)*pl211;
    psq10=0.5*pow(rr,2)*pl210;
    psq01=0.5*pow(rr,2)*pl201;
    psq00=0.5*pow(rr,2)*pl200;
    psqp11=rr*pl211+0.5*pow(rr,2)*pl2p11;
    psqp10=rr*pl210+0.5*pow(rr,2)*pl2p10;
    psqp01=rr*pl201+0.5*pow(rr,2)*pl2p01;
    psqp00=rr*pl200+0.5*pow(rr,2)*pl2p00;
    psqdp11=pl211+2.*rr*pl2p11+0.5*pow(rr,2)*pl2dp11;
    psqdp10=pl210+2.*rr*pl2p10+0.5*pow(rr,2)*pl2dp10;
    psqdp01=pl201+2.*rr*pl2p01+0.5*pow(rr,2)*pl2dp01;
    psqdp00=pl200+2.*rr*pl2p00+0.5*pow(rr,2)*pl2dp00;
    vp[0]=0.0625*(9.*psq11+3.*psq10+3.*psq01+psq00);
    vp[1]=0.0625*(3.*psq11-3.*psq10  +psq01-psq00);
    vp[2]=0.0625*(3.*psq11  +psq10-3.*psq01-psq00);
    vp[3]=0.0625*(  psq11  -psq10  -psq01+psq00);
    vp[4]=0.0625*(9.*psqp11+3.*psqp10+3.*psqp01+psqp00);
    vp[5]=0.0625*(3.*psqp11-3.*psqp10+  psqp01-psqp00);
    vp[6]=0.0625*(3.*psqp11+  psqp10-3.*psqp01-psqp00);
    vp[7]=0.0625*(  psqp11-  psqp10-  psqp01+psqp00);
    vp[8]= 0.0625*(9.*psqdp11+3.*psqdp10+3.*psqdp01+psqdp00);
    vp[9]=0.0625*(3.*psqdp11-3.*psqdp10+  psqdp01-psqdp00);
    vp[10]=0.0625*(3.*psqdp11+  psqdp10-3.*psqdp01-psqdp00);
    vp[11]=0.0625*(  psqdp11-  psqdp10-  psqdp01+psqdp00);
}


// ----------------------
// super-soft core(c) v14
// lpot = 18 -> v14
//        19 -> mod v8'
// ----------------------
void DLM_StefanoPotentials::SuperSoftCoreV14(){
    const double& rr = *CurrentVerySmallRadius;
    xVar=0.7*rr;
    rr4=pow(rr,4.);
    rc4=1.-exp(-rr4);
    rc6=1.-exp(-pow(rr,6));
    hrVar=10.463;
    p11=144.83*exp(-rr4/pow(0.88787,2))
        +(-241.34*yc(3.3788*xVar, xVar)+(hrVar/3)*yc(xVar, xVar))*rc4;
    p10=215.32*exp(-rr4/pow(0.85807,2))
        +(-883.6*yc(3.5042*xVar, xVar)-hrVar*yc(xVar, xVar))*rc4;
    p01=375.*exp(-rr4/pow(0.47552,2))
        +(-1001.6*yc(3.6071*xVar, xVar)-hrVar*yc(xVar, xVar))*rc4;
    p00=75.653*exp(-rr4/pow(3.0000,2))
        +(-286.26*yc(2.0254*xVar, xVar)+3*hrVar*yc(xVar, xVar))*rc4;
    pt1=36.*exp(-rr4/pow(1.0805,2))
        +(-110.*yt(3.9529*xVar, xVar)+(hrVar/3)*yt(xVar, xVar))*rc6;
    pt0=-58.951*exp(-rr4/pow(1.3171,2))
        +(395.18*yt(4.3098*xVar, xVar)-hrVar*yt(xVar, xVar))*rc6;
    pls1=(520.*yls(5.661*xVar, xVar)-54.85*yls(4.0141*xVar, xVar))*rc6;
    pls0=(-40.466*yls(5.768*xVar, xVar)-40.408*yls(4.0676*xVar, xVar))*rc6;
    pl211=(6.65*yl2(1.965*xVar, xVar)-0.959*yl2(xVar, xVar))*rc6;
    pl210=(17.626*yl2(2.6463*xVar, xVar)-0.35261*yl2(xVar, xVar))*rc6;
    pl201=(14.*yl2(2.5*xVar, xVar)-0.35*yl2(xVar, xVar))*rc6;
    pl200=(15.633*yl2(2.01*xVar, xVar)+0.72581*yl2(xVar, xVar))*rc6;
    pq0=-3.9904*yl2(2.4583*xVar, xVar)*rc6;
// fix v8'
    if (lpot==19){
        p00=p00+2.*pl200;
        pls0=pls0-2.*pl210-10.*pq0;
        p11=p11+2.*pl211;
// ------------------------
// option for modified v8'
// ------------------------
         p11=p11-111.*yc(3.3788*xVar, xVar)*rc4;
// ------------------------
    }
    vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
    vv[1]=0.0625*(3.*p11-3.*p10+  p01-p00);
    vv[2]=0.0625*(3.*p11+  p10-3.*p01-p00);
    vv[3]=0.0625*(  p11-  p10-  p01+p00);
    vv[4]=0.25*(3.*pt1+pt0);
    vv[5]=0.25*(  pt1-pt0);
    vv[6]=0.25*(3.*pls1+pls0)+0.75*pq0;
    vv[7]=0.25*(  pls1-pls0)-0.75*pq0;
    if (lpot==19) return;
    vv[8]= 0.0625*(9.*pl211+3.*pl210+3.*pl201+pl200)-0.75*pq0;
    vv[9]=0.0625*(3.*pl211-3.*pl210+  pl201-pl200)+0.75*pq0;
    vv[10]=0.0625*(3.*pl211+  pl210-3.*pl201-pl200)-0.25*pq0;
    vv[11]=0.0625*(  pl211-  pl210-  pl201+pl200)+0.25*pq0;
    vv[12]=1.5*pq0;
    vv[13]=-1.5*pq0;
    vw[1]=(hrVar/3.)*yc(xVar, xVar);
    vw[2]=(hrVar/3.)*yt(xVar, xVar);
    //rr=rrsave
}

void DLM_StefanoPotentials::Paris(){
    const double& rr = *CurrentRadius;
    p01=0;
    p11=0;
    p00=0;
    p10=0;
    pt1=0;
    pt0=0;
    pls1=0;
    pls0=0;
    pl201=0;
    pl211=0;
    pl200=0;
    pl210=0;
    pls21=0;
    pls20=0;
    for(int i=0; i<12; i++){
        if (rr<=small){
            rm1=rr*pm1[i];
            rm0=rr*pm0[i];
            vc1=-1.+0.5*rm1;
            vc0=-1.+0.5*rm0;
            vt1=0.125*rm1;
            vt0=0.125*rm0;
            vls1=1./3.-0.125*rm1;
            vls0=1./3.-0.125*rm0;
            vso21=-1./15.+rm1/48.;
            vso20=-1./15.+rm0/48.;
        }
        else{
            vc1=pc(pm1[i]*rr);
            vc0=pc(pm0[i]*rr);
            vt1=pt(pm1[i]*rr)*vc1;
            vt0=pt(pm0[i]*rr)*vc0;
            vls1=pls(pm1[i]*rr)*vc1;
            vls0=pls(pm0[i]*rr)*vc0;
            vso21=pso2(pm1[i]*rr)*vc1;
            vso20=pso2(pm0[i]*rr)*vc0;
        }
        p01=p01+pgva01[i]*vc1;
        p11=p11+pgva11[i]*vc1;
        p00=p00+pgva00[i]*vc0;
        p10=p10+pgva10[i]*vc0;
        pt1=pt1+pgvt1[i]*vt1;
        pt0=pt0+pgvt0[i]*vt0;
        pls1=pls1+pgvls1[i]*vls1;
        pls0=pls0+pgvls0[i]*vls0;
        pl201=pl201+pgvb01[i]*vc1;
        pl211=pl211+pgvb11[i]*vc1;
        pl200=pl200+pgvb00[i]*vc0;
        pl210=pl210+pgvb10[i]*vc0;
        pls21=pls21+pgvso21[i]*vso21;
        pls20=pls20+pgvso20[i]*vso20;
    }

    vv[0]=0.0625*(9.*p11+3.*p10+3.*p01+p00);
    vv[1]=0.0625*(3.*p11-3.*p10+  p01-p00);
    vv[2]=0.0625*(3.*p11+  p10-3.*p01-p00);
    vv[3]=0.0625*(  p11-  p10-  p01+p00);
    vv[4]=0.25*(3.*pt1+pt0);
    vv[5]=0.25*(  pt1-pt0);
    vv[6]=0.25*(3.*pls1+pls0);
    vv[7]=0.25*(  pls1-pls0);
    vv[8]= 0.0625*(9.*pl211+3.*pl210+3.*pl201+pl200);
    vv[9]=0.0625*(3.*pl211-3.*pl210+  pl201-pl200);
    vv[10]=0.0625*(3.*pl211+  pl210-3.*pl201-pl200);
    vv[11]=0.0625*(  pl211-  pl210-  pl201+pl200);
    vv[12]=0.25*(3.*pls21+pls20);
    vv[13]=0.25*(  pls21-pls20);
}

// lr:   potential returned is v(r)*r**lr
// rr:   relative position r in fm
void DLM_StefanoPotentials::pot(const double& rr, const double& lr){
    const double rrsmall = rr<=small?small:rr;
    const double rrvsmall = rr<=vsmall?vsmall:rr;
    CurrentRadius = &rr;
    CurrentSmallRadius = &rrsmall;
    CurrentVerySmallRadius = &rrvsmall;

    for(int i=0; i<18; i++) vv[i] = 0;
    for(int i=0; i<12; i++) vp[i] = 0;
    for(int i=0; i<14; i++) vw[i] = 0;

    const double rr427 = 4.27*rrsmall;

    double ff=1-exp(-rr427)*(1.+rr427*(33.+rr427*(9.+rr427))/48.);
    vw[0]=1.4399652*ff/rr;

    switch(lpot){
        case 1 : MalflietTjon(); break;
        case 2 : ReidV8(); break;
        case 3 : UrbanaV14(); break;
        case 4 : ArgonneV14(); break;
        case 5 : ArgonneV14(); break;
        case 6 : ArgonneV14(); break;
        case 7 : ArgonneV18(); break;
        case 8 : ArgonneV18(); break;
        case 9 : ArgonneV18(); break;
        case 10 : ArgonneV18(); break;
        case 11 : ArgonneV18(); break;
        case 12 : ArgonneV18(); break;
        case 13 : ArgonneV18(); break;
        case 14 : ArgonneV18(); break;
        case 15 : ArgonneV18(); break;
        case 16 : ArgonneV18p(); break;
        case 17 : ArgonneV18p(); break;
        case 18 : SuperSoftCoreV14(); break;
        case 19 : SuperSoftCoreV14(); break;
        case 20 : Paris(); break;
        case 21 : ArgonneV18(); break;
        case 22 : ArgonneV18(); break;
        case 23 : ArgonneV18(); break;
        case 24 : ArgonneV18(); break;
        case 112: cheft(rr); break;
        case 113: cheft(rr); break;
        case 114: cheft(rr); break;
        case 122: cheft(rr); break;
        case 123: cheft(rr); break;
        case 124: cheft(rr); break;
        default : break;
    }

    if(lr>=1){
        for(int l=0; l<18; l++){
            //printf("rr=%f; lr=%f; vv[%i]=%f; ",rr,lr,l,vv[l]);
            //printf("pow(%f,%f)=%f; ",rr,lr,pow(rr,lr));
            vv[l] = vv[l]*pow(rr,lr);
            //printf("NEWvv[%i]=%f\n",l,vv[l]);
            //usleep(100e3);
        }
        vw[0]=vw[0]*pow(rr,lr);
        vw[1]=vw[1]*pow(rr,lr);
        vw[2]=vw[2]*pow(rr,lr);
        vw[3]=vw[3]*pow(rr,lr);
        vw[4]=vw[4]*pow(rr,lr);
    }

}

void DLM_StefanoPotentials::gauleg(double& x1, double& x2, double* x, double* w, const int& n){
    const double EPS = 1e-14;
    double p1,p2,p3,pp,xl,xm,z,z1;
    int m;
    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    //! keep in mind that going from FORTRAN to C++ we change i -> i+1 with limits starting from 0 instead of 1
    //this is also true for j. As a result some of the expressions below change!
    for(int i=0; i<m; i++){
        z=cos(3.141592654*(double(i)+0.75)/(double(n)+0.5));
        GoToSucks1:
        p1=1;
        p2=0;
        for(int j=0; j<n; j++){
            p3=p2;
            p2=p1;
            p1=((2.*double(j)+1.)*z*p2-double(j)*p3)/double(j+1);
        }
        pp=double(n)*(z*p1-p2)/(z*z-1.);
        z1=z;
        z=z1-p1/pp;
        if(fabs(z-z1)>EPS) goto GoToSucks1;
        x[i]=xm-xl*z;
        x[n-1-i]=xm+xl*z;//!error prone by translation
        w[i]=2.*xl/((1.-z*z)*pp*pp);
        w[n-1-i]=w[i];//!error prone by translation
    }
}

void DLM_StefanoPotentials::spectralNLO(const double& rad){

    const int n_int=50;

    double* Xgauleg = new double [n_int];
    double* Wgauleg = new double [n_int];

    double f_pi_o_fm, m_pi_o_fm, L_cut_o_fm;
    double fmpisq, overall, co1, co2, co3, res;
    double a, b, x;

    f_pi_o_fm=f_pi/hc197;
    m_pi_o_fm=m_pi/hc197;
    fmpisq = 4.0*(m_pi_o_fm*m_pi_o_fm);
    L_cut_o_fm=Lam/hc197;

    a = 2.0*m_pi_o_fm;
    b = L_cut_o_fm;

    gauleg(a,b,Xgauleg,Wgauleg,n_int);

    W_C = 0.0;
    V_S = 0.0;
    V_T = 0.0;

    for(int i=0; i<n_int; i++){
        x = Xgauleg[i];

        overall = 1.0/(768.0*pi*pow(f_pi_o_fm,4));
        co1 = 5.0*pow(g_a,4) - 4.0*pow(g_a,2) - 1.0;
        co2 = 23.0*pow(g_a,4) - 10.0*pow(g_a,2) - 1.0;
        co3 = 48.0*pow(g_a,4)*pow(m_pi_o_fm,4);
        res = exp(-x*rad)*sqrt(pow(x,2) - fmpisq);
        res = res*(fmpisq*co1 - pow(x,2)*co2 + co3/(fmpisq - pow(x,2)));
        res = res*overall;
        res = (hc197*res)/(2.0*pow(pi,2)*rad);
        W_C = W_C + Wgauleg[i]*res;

        overall = 3.0*pow(g_a,4)/(128.0*pi*pow(f_pi_o_fm,4));
        res = exp(-x*rad)*sqrt(x*x - fmpisq);
        res = res*x*x*overall;
        res = (hc197*res)/(3.0*pow(pi,2)*rad);
        V_S = V_S + Wgauleg[i]*res;

        overall = 3.0*pow(g_a,4)/(128.0*pi*pow(f_pi_o_fm,4));
        res = exp(-x*rad)*sqrt(x*x - fmpisq);
        res = res*(3.0 + 3.0*x*rad + (x*x)*(rad*rad))*overall;
        res = -(hc197*res)/(6.0*pow(pi,2)*pow(rad,3));
        V_T = V_T + Wgauleg[i]*res;
    }

    delete [] Xgauleg;
    delete [] Wgauleg;
}

void DLM_StefanoPotentials::spectralN2LO(const double& rad){

    double f_pi_o_fm, m_pi_o_fm, L_cut_o_fm;
    double fpifour, overall, co1, co2;
    double x, y;

    f_pi_o_fm=f_pi/hc197;
    m_pi_o_fm=m_pi/hc197;
    fpifour = pow(f_pi_o_fm,4);
    L_cut_o_fm=Lam/hc197;

    x = m_pi_o_fm*rad;
    y = L_cut_o_fm*rad;

    overall = (3.0*pow(g_a,2))/(32.0*pow(pi,2)*fpifour);
    co1 = 2.0*EGM_C1*pow(x,2)*pow(1.0+x,2)
          + EGM_C3*(6.0+12.0*x+10.0*pow(x,2)+4.0*pow(x,3)+pow(x,4));
    co1 = co1*exp(-2.0*x)/pow(rad,6);
    co2 = 4.0*EGM_C1*pow(x,2)*(2.0+y*(2.0+y)-2.0*pow(x,2))
          + EGM_C3*(24.0 + y*(24.0+12.0*y+4.0*pow(y,2)+pow(y,3))
                - 4.0*pow(x,2)*(2.0+2.0*y+pow(y,2)) + 4.0*pow(x,4));
    co2 = -co2*exp(-y)/(4.0*pow(rad,6));
    V_C = overall*(co1 + co2);
    V_C = V_C*pow(hc197,2);

    overall = pow(g_a,2)/(48.0*pow(pi,2)*fpifour);
    co1 = EGM_C4*(1.0+x)*(3.0+3.0*x+2.0*pow(x,2));
    co1 = co1*exp(-2.0*x)/pow(rad,6);
    co2 = EGM_C4*(24.0+24.0*y+12.0*pow(y,2)+4.0*pow(y,3)+pow(y,4)
          - 4.0*pow(x,2)*(2.0+2.0*y+pow(y,2)));
    co2 = -co2*exp(-y)/(8.0*pow(rad,6));
    W_S = overall*(co1 + co2);
    W_S = W_S*pow(hc197,2);

    overall = pow(g_a,2)/(48.0*pow(pi,2)*fpifour);
    co1 = EGM_C4*(1.0+x)*(3.0+3.0*x+pow(x,2));
    co1 = -co1*exp(-2.0*x)/pow(rad,6);
    co2 = EGM_C4*(48.0+48.0*y+24.0*pow(y,2)+7.0*pow(y,3)+pow(y,4)
          - 4.0*pow(x,2)*(8.0+5.0*y+pow(y,2)));
    co2 = co2*exp(-y)/(16.0*pow(rad,6));
    W_T = overall*(co1 + co2);
    W_T = W_T*pow(hc197,2);

}

void DLM_StefanoPotentials::cheft(const double& rr){
    if(lpot>=112 && lpot<=114){
        r_0p=1.0;
        LO_C_s=-0.751120;
        LO_C_t=0.374089;
        LO_CIB=-0.023605;
        LO_CSB=-0.019876;
        NLO_C_s=3.168026;
        NLO_C_t=1.413952;
        NLO_C[0]=0.314202;
        NLO_C[1]=0.257857;
        NLO_C[2]=-0.131344;
        NLO_C[3]=0.118613;
        NLO_C[4]=-2.385514;
        NLO_C[5]=0.373188;
        NLO_C[6]=-0.356684;
        NLO_CIB=0.050944;
        NLO_CSB=0.008231;
        N2LO_C_s=5.438495;
        N2LO_C_t=0.276718;
        N2LO_C[0]=-0.140842;
        N2LO_C[1]=0.042427;
        N2LO_C[2]=-0.123378;
        N2LO_C[3]=0.110184;
        N2LO_C[4]=-2.112533;
        N2LO_C[5]=0.158979;
        N2LO_C[6]=-0.269935;
        N2LO_CIB=0.053204;
        N2LO_CSB=0.009759;
    }
    else if(lpot>=122 && lpot<=124){
        r_0p=1.2;
        LO_C_s=-1.796928;
        LO_C_t=0.154414;
        LO_CIB=-0.013348;
        LO_CSB=-0.019589;
        NLO_C_s=0.035507;
        NLO_C_t=0.717286;
        NLO_C[0]=0.222888;
        NLO_C[1]=0.228779;
        NLO_C[2]=-0.150425;
        NLO_C[3]=0.089285;
        NLO_C[4]=-2.029313;
        NLO_C[5]=0.340112;
        NLO_C[6]=-0.362474;
        NLO_CIB=0.054768;
        NLO_CSB=0.006596;
        N2LO_C_s=2.687641;
        N2LO_C_t=0.233817;
        N2LO_C[0]=-0.079513;
        N2LO_C[1]=0.076102;
        N2LO_C[2]=-0.169261;
        N2LO_C[3]=0.123588;
        N2LO_C[4]=-1.942800;
        N2LO_C[5]=0.214206;
        N2LO_C[6]=-0.341926;
        N2LO_CIB=0.056482;
        N2LO_CSB=0.007710;
    }

    double rad = rr>0.001?rr:0.001;

    for(int i=0; i<18; i++) vv[i] = 0;
    g_r0p = 1.0/(pi*g34*pow(r_0p,3));

    //Note that the following is used both for contacts and to regulate OPE
    del_r = exp(-pow(rad/r_0p,4));

    //Also note that contact parameters are always multiplied by hc197

    //We take OPE as Goldberger-Treiman modified, i.e. the same, at all orders
    //x = r*mu = r*m_pi*c/hbar
    x_not = m_pi_not*rad/hc197;
    pre_not = (pow(m_pi_not,3)/(12.0*pi))*pow(g_a_Tre/(2*f_pi),2)*
                  exp(-x_not)/x_not;
    x_plus = m_pi_plus*rad/hc197;
    pre_plus = (pow(m_pi_plus,3)/(12.0*pi))*pow(g_a_Tre/(2*f_pi),2)*
                   exp(-x_plus)/x_plus;

    //The pion masses are according to Eq. (2.7) in the 1997 Pudliner PRC
    sp_not = (1.0-del_r)*pre_not;
    sp_plus = (1.0-del_r)*pre_plus;
    vv[3] = (sp_not + 2.0*sp_plus)/3.0;
    ten_not = (1.0-del_r)*pre_not*(1.0 + 3.0/x_not +
              3.0/pow(x_not,2));
    ten_plus = (1.0-del_r)*pre_plus*(1.0 + 3.0/x_plus +
               3.0/pow(x_plus,2));
    vv[5] = (ten_not + 2.0*ten_plus)/3.0;
    vv[15] = (sp_not - sp_plus)/3.0;
    vv[16] = (ten_not - ten_plus)/3.0;

    if (lpot==112 || lpot==122){
        //We are including charge-breaking at each order
        vv[0] += hc197*(LO_CIB/2.0)*del_r*g_r0p;
        vv[1] += hc197*(LO_CIB/6.0)*del_r*g_r0p;
        vv[14] += hc197*(LO_CIB/6.0)*del_r*g_r0p;
        vv[17] += hc197*(LO_CSB/2.0)*del_r*g_r0p;

        //LO has 2 contacts
        vv[0] += hc197*LO_C_s*del_r*g_r0p;
        vv[2] += hc197*LO_C_t*del_r*g_r0p;
    }
    else if(lpot==113 || lpot==123){
        //We are including charge-breaking at each order
        vv[0] += hc197*(NLO_CIB/2.0)*del_r*g_r0p;
        vv[1] += hc197*(NLO_CIB/6.0)*del_r*g_r0p;
        vv[14] += hc197*(NLO_CIB/6.0)*del_r*g_r0p;
        vv[17] += hc197*(NLO_CSB/2.0)*del_r*g_r0p;

        //At NLO the LO 2 contacts change
        vv[0] += hc197*NLO_C_s*del_r*g_r0p;
        vv[2] += hc197*NLO_C_t*del_r*g_r0p;

        //NLO also has a TPE
        spectralNLO(rad);
        vv[1] += (1.0-del_r)*W_C;
        vv[2] += (1.0-del_r)*V_S;
        vv[4] += (1.0-del_r)*V_T;

        //NLO also has 7 new contacts
        lap_res = 20.0*(pow(rad,2)/pow(r_0p,4)) - 16.0*pow(rad,6)/pow(r_0p,8);
        vv[0] += hc197*NLO_C[0]*lap_res*del_r*g_r0p;
        vv[1] += hc197*NLO_C[1]*lap_res*del_r*g_r0p;
        vv[2] += hc197*NLO_C[2]*lap_res*del_r*g_r0p;
        vv[3] += hc197*NLO_C[3]*lap_res*del_r*g_r0p;
        vv[2] += hc197*NLO_C[5]*(lap_res/3.0)*del_r*g_r0p;
        vv[3] += hc197*NLO_C[6]*(lap_res/3.0)*del_r*g_r0p;
        other_res = (8.0*(pow(rad,2)/pow(r_0p,4)) - 16.0*pow(rad,6)/pow(r_0p,8))/3.0;
        vv[4] += hc197*NLO_C[5]*other_res*del_r*g_r0p;
        vv[5] += hc197*NLO_C[6]*other_res*del_r*g_r0p;
        vv[6] = hc197*NLO_C[4]*2.0*(pow(rad,2)/pow(r_0p,4))*del_r*g_r0p;
    }
    else if (lpot==114 || lpot==124){
        //We are including charge-breaking at each order
        vv[0] += hc197*(N2LO_CIB/2.0)*del_r*g_r0p;
        vv[1] += hc197*(N2LO_CIB/6.0)*del_r*g_r0p;
        vv[14] += hc197*(N2LO_CIB/6.0)*del_r*g_r0p;
        vv[17] += hc197*(N2LO_CSB/2.0)*del_r*g_r0p;

        //At N2LO the LO 2 contacts change
        vv[0] += hc197*N2LO_C_s*del_r*g_r0p;
        vv[2] += hc197*N2LO_C_t*del_r*g_r0p;

        //N2LO also has NLO's TPE
        spectralNLO(rad);
        vv[1] += (1.0-del_r)*W_C;
        vv[2] += (1.0-del_r)*V_S;
        vv[4] += (1.0-del_r)*V_T;

        //N2LO also has a different TPE
        spectralN2LO(rad);
        vv[0] += (1.0-del_r)*V_C;
        vv[3] += (1.0-del_r)*W_S;
        vv[5] += (1.0-del_r)*W_T;

        //N2LO also has new parameters for the 7 NLO contacts
        lap_res = 20.0*pow(rad,2)/pow(r_0p,4) - 16.0*pow(rad,6)/pow(r_0p,8);
        vv[0] += hc197*N2LO_C[0]*lap_res*del_r*g_r0p;
        vv[1] += hc197*N2LO_C[1]*lap_res*del_r*g_r0p;
        vv[2] += hc197*N2LO_C[2]*lap_res*del_r*g_r0p;
        vv[3] += hc197*N2LO_C[3]*lap_res*del_r*g_r0p;
        vv[2] += hc197*N2LO_C[5]*(lap_res/3.0)*del_r*g_r0p;
        vv[3] += hc197*N2LO_C[6]*(lap_res/3.0)*del_r*g_r0p;
        other_res = (8.0*pow(rad,2)/pow(r_0p,4) - 16.0*pow(rad,6)/pow(r_0p,8))/3.0;
        vv[4] += hc197*N2LO_C[5]*other_res*del_r*g_r0p;
        vv[5] += hc197*N2LO_C[6]*other_res*del_r*g_r0p;
        vv[6] = hc197*N2LO_C[4]*2.0*pow(rad,2)/pow(r_0p,4)*del_r*g_r0p;
    }
}


double DLM_StefanoPotentials::Eval(const double& Radius, const int& Spin, const int& Spin3, const int& IsoSpin, const int& IsoSpin3){
    EvalValue = 0;
    pot(Radius, 0);
    //pp
    if(IsoSpin==1 && IsoSpin3==1 && Spin==0){
        EvalValue += vv[0];
        EvalValue += vv[1];
        EvalValue += -3.*vv[2];
        EvalValue += -3.*vv[3];
        ////EvalValue += vv[4];
        ////EvalValue += vv[5];
        ////EvalValue += vv[6];
        ////EvalValue += vv[7];
        ////EvalValue += vv[8];
        ////EvalValue += vv[9];
        ////EvalValue += vv[10];
        ////EvalValue += vv[11];
        ////EvalValue += vv[12];
        ////EvalValue += vv[13];
        EvalValue += 2.*vv[14];//*2.;
        EvalValue += -6.*vv[15];//*2.;
        ////EvalValue += vv[16];//*2.;
        EvalValue += 2.*vv[17];//*2.;
    }
    else if(IsoSpin==1 && IsoSpin3==1 && Spin==1){
        EvalValue += vv[0];
        EvalValue += vv[1];
        EvalValue += vv[2];
        EvalValue += vv[3];
        ////EvalValue += 2.*vv[4];
        ////EvalValue += 2.*vv[5];
        ////EvalValue += -2.*vv[6];
        ////EvalValue += -2.*vv[7];
        ////EvalValue += vv[8];
        ////EvalValue += vv[9];
        ////EvalValue += vv[10];
        ////EvalValue += vv[11];
        ////EvalValue += vv[12];
        ////EvalValue += vv[13];
        EvalValue += 2.*vv[14];//*2.;
        EvalValue += 2.*vv[15];//*2.;
        ////EvalValue += vv[16];//*2.;
        EvalValue += 2.*vv[17];//*2.;
    }
    else if(IsoSpin==1 && IsoSpin3==1){
        //! I STILL DO NOT KNOW HOW TO TAKE THE SPIN INTO ACCOUNT!
        EvalValue += vv[0];
        EvalValue += vv[2];
        EvalValue += vv[4];
        EvalValue += vv[6];
        EvalValue += vv[8];
        EvalValue += vv[10];
        EvalValue += vv[12];
        EvalValue += vv[1];
        EvalValue += vv[3];
        EvalValue += vv[5];
        EvalValue += vv[7];
        EvalValue += vv[9];
        EvalValue += vv[11];
        EvalValue += vv[13];
        EvalValue += vv[14];//*2.;
        EvalValue += vv[15];//*2.;
        EvalValue += vv[16];//*2.;
        EvalValue += vv[17];//*2.;
    }
    else{
        printf("WARNING -> this mode is not implemented yet!\n");
    }
    //printf("EvalValue(%f)=%f\n",Radius, EvalValue);
    //usleep(100e3);
    return EvalValue;
}

double DLM_StefanoPotentials::EvalCATS_v1_0(const double& Radius, const int& Spin){
    return Eval(Radius, Spin, Spin, 1, 1);
}

//DlmFlag : 1 => modified s12 for l==j-1 (dlm magic)
//DlmFlag : 2 => take the coupled channel, but compute only the first diagonal term
//DlmFlag : 3 => take the coupled channel, but compute only the off-diagonal term
//DlmFlag : 4 => take the coupled channel, but compute only the second diagonal term
//IsoSpin is for the particle pair, PartType1/2 should be 1 for proton and -1 for neutron
double DLM_StefanoPotentials::Eval_PWprojector(const double& Radius,
                                               const int& IsoSpin, const int& PartType1, const int& PartType2,
                                               const int& Spin, const int& AngMom, const int& TotMom, const int& DlmFlag){
    pot(Radius, 0);

    const int& t1z=PartType1;//isospin3 x 2 particle1
    const int& t2z=PartType2;//isospin3 x 2 particle2
    const int s1ds2=4*Spin-3;
    //this term should be zero in case one of the particles has IsoSpin 0
    const int t1dt2=(PartType1*PartType2==0)?0:4*IsoSpin-3;
    const int t12=3*t1z*t2z-t1dt2;
    double vc = vv[0]+t1dt2*vv[1]+s1ds2*vv[2]+s1ds2*t1dt2*vv[3]+t12*vv[14]+s1ds2*t12*vv[15]+(t1z+t2z)*vv[17];
    double vt = vv[4]+t1dt2*vv[5]+t12*vv[16];
    double vls = vv[6]+t1dt2*vv[7];
    double vl2=vv[8]+t1dt2*vv[9]+s1ds2*vv[10]+s1ds2*t1dt2*vv[11];
    double vls2=vv[12]+t1dt2*vv[13];

    const int ls=(TotMom*(TotMom+1)-AngMom*(AngMom+1)-Spin*(Spin+1))/2;
    double s12=0;
    double s12m;
    double s12p;
    const double lsm=TotMom-1;
    const double lsp=-(TotMom+2);

    double RETURN_VAL=0;
//printf("r=%.2f, s=%i, l=%i, j=%i, FLAG=%i\n",Radius,Spin,AngMom,TotMom,DlmFlag);
    if(Spin==1 && AngMom==TotMom) s12=2;
    else if(AngMom==(TotMom+1)) {s12=-2.*double(TotMom+2)/double(2*TotMom+1);}
    //this I made up to make it work... it is -0.6 (-3/5) for s=1, l=1, j=2
    else if(AngMom==(TotMom-1) && DlmFlag==1) {s12=-1.*double(TotMom+1)/double(2*TotMom+1);}
    else if(Spin==1 && AngMom==(TotMom-1) && (DlmFlag==2 || DlmFlag==3 || DlmFlag==4 || DlmFlag==5)){
        s12m=-2.*double(TotMom-1.)/double(2.*TotMom+1.);
        s12=sqrt(double(36.*TotMom*(TotMom+1)))/double(2.*TotMom+1.);
        s12p=-2.*double(TotMom+2.)/double(2.*TotMom+1.);
    }

    RETURN_VAL = vc+s12*vt+ls*vls+AngMom*(AngMom+1)*vl2+ls*ls*vls2;
    if(Spin==1 && AngMom==(TotMom-1)){
        if(DlmFlag==2) RETURN_VAL = vc+s12m*vt+lsm*vls+AngMom*(AngMom+1)*vl2+lsm*lsm*vls2;
        else if(DlmFlag==3) RETURN_VAL = s12*vt;
        else if(DlmFlag==4) RETURN_VAL = vc+s12p*vt+lsp*vls+(AngMom+2.)*(AngMom+3.)*vl2+lsp*lsp*vls2;
    }

    return RETURN_VAL;
}

//DlmFlag : 1 => modified s12 for l==j-1
//DlmFlag : 2 => take the coupled channel, but compute only the first diagonal term
//DlmFlag : 3 => take the coupled channel, but compute only the off-diagonal term
//DlmFlag : 4 => take the coupled channel, but compute only the second diagonal term
double DLM_StefanoPotentials::Eval_PWprojector_pp(const double& Radius, const int& Spin,
                                               const int& AngMom, const int& TotMom, const int& DlmFlag){
    pot(Radius, 0);

    const int t1z=1;
    const int t2z=1;
    const int IsoSpin = 1;
    const int s1ds2=4*Spin-3;
    const int t1dt2=4*IsoSpin-3;
    const int t12=3*t1z*t2z-t1dt2;
    double vc = vv[0]+t1dt2*vv[1]+s1ds2*vv[2]+s1ds2*t1dt2*vv[3]+t12*vv[14]+s1ds2*t12*vv[15]+(t1z+t2z)*vv[17];
    double vt = vv[4]+t1dt2*vv[5]+t12*vv[16];
    double vls = vv[6]+t1dt2*vv[7];
    double vl2=vv[8]+t1dt2*vv[9]+s1ds2*vv[10]+s1ds2*t1dt2*vv[11];
    double vls2=vv[12]+t1dt2*vv[13];

    const int ls=(TotMom*(TotMom+1)-AngMom*(AngMom+1)-Spin*(Spin+1))/2;
    double s12=0;
    double s12m;
    double s12p;
    const double lsm=TotMom-1;
    const double lsp=-(TotMom+2);

    double RETURN_VAL=0;
//printf("r=%.2f, s=%i, l=%i, j=%i, FLAG=%i\n",Radius,Spin,AngMom,TotMom,DlmFlag);
    if(Spin==1 && AngMom==TotMom) s12=2;
    else if(AngMom==(TotMom+1)) {s12=-2.*double(TotMom+2)/double(2*TotMom+1);}
    //this I made up to make it work... it is -0.6 (-3/5) for s=1, l=1, j=2
    else if(AngMom==(TotMom-1) && DlmFlag==1) {s12=-1.*double(TotMom+1)/double(2*TotMom+1);}
    else if(Spin==1 && AngMom==(TotMom-1) && (DlmFlag==2 || DlmFlag==3 || DlmFlag==4 || DlmFlag==5)){
        s12m=-2.*double(TotMom-1.)/double(2.*TotMom+1.);
        s12=sqrt(double(36.*TotMom*(TotMom+1)))/double(2.*TotMom+1.);
        s12p=-2.*double(TotMom+2.)/double(2.*TotMom+1.);
    }

    RETURN_VAL = vc+s12*vt+ls*vls+AngMom*(AngMom+1)*vl2+ls*ls*vls2;
    if(Spin==1 && AngMom==(TotMom-1)){
        if(DlmFlag==2) RETURN_VAL = vc+s12m*vt+lsm*vls+AngMom*(AngMom+1)*vl2+lsm*lsm*vls2;
        else if(DlmFlag==3) RETURN_VAL = s12*vt;
        else if(DlmFlag==4) RETURN_VAL = vc+s12p*vt+lsp*vls+(AngMom+2.)*(AngMom+3.)*vl2+lsp*lsp*vls2;
    }

    return RETURN_VAL;
}

void DLM_StefanoPotentials::PotentialName(const int& WhichPot, char* Name){
    switch(WhichPot){
    case 1 : strcpy(Name, "Malfliet-Tjon"); break;
    case 2 : strcpy(Name, "Reid V8"); break;
    case 3 : strcpy(Name, "Urbana V14"); break;
    case 4 : strcpy(Name, "Argonne V8*"); break;
    case 5 : strcpy(Name, "Argonne V9"); break;
    case 6 : strcpy(Name, "Argonne V14"); break;
    case 7 : strcpy(Name, "Argonne V18-csbl"); break;
    case 8 : strcpy(Name, "Argonne V18-csbs"); break;
    case 9 : strcpy(Name, "Argonne V18"); break;
    case 10: strcpy(Name, "Argonne V8"); break;
    case 11: strcpy(Name, "Argonne V6"); break;
    case 12: strcpy(Name, "Argonne V4"); break;
    case 13: strcpy(Name, "Argonne V2"); break;
    case 14: strcpy(Name, "Argonne V1"); break;
    case 15: strcpy(Name, "Argonne VX"); break;
    case 16: strcpy(Name, "Argonne V18p"); break;
    case 17: strcpy(Name, "Argonne V18pq"); break;
    case 18: strcpy(Name, "Super-soft core V14"); break;
    case 19: strcpy(Name, "Super-soft core V8*"); break;
    case 20: strcpy(Name, "Paris"); break;
    case 21: strcpy(Name, "Argonne V18 1.9"); break;
    case 22: strcpy(Name, "Argonne V8 1.9"); break;
    case 23: strcpy(Name, "Argonne V18 1.7"); break;
    case 24: strcpy(Name, "Argonne V8 1.7"); break;
    case 112: strcpy(Name, "LO, R0=1.0"); break;
    case 113: strcpy(Name, "NLO, R0=1.0"); break;
    case 114: strcpy(Name, "N2LO, R0=1.0"); break;
    case 122: strcpy(Name, "LO, R0=1.2"); break;
    case 123: strcpy(Name, "NLO, R0=1.2"); break;
    case 124: strcpy(Name, "N2LO, R0=1.2"); break;
    default : strcpy(Name, "Unknown potential"); break;
    }
}
