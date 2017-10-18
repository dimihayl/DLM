

#ifndef DLM_STEFANOPOTENTIALS_H
#define DLM_STEFANOPOTENTIALS_H

#include "math.h"
#include <stdio.h>

class DLM_StefanoPotentials{

private:

    const double m_pi;
    const double m_pi_not;
    const double m_pi_plus;
    const double m_n;
    const double g_a;
    const double g_a_Tre;
    const double f_pi;
    const double pi;
    const double hc197;
    const double g34;
    const double EGM_C1;
    const double EGM_C3;
    const double EGM_C4;
    const double Lam;
    const double small;
    const double vsmall;

    const double xmn, gam, rho, chi, omg, ftp;
    const double mp,mn;


    double* vv;//size 18
    double* vp;//size 12
    double* vw;//size 14
    double mpi0,mpic,mpi,mpis,mpiq,mu0,muc,mu,mus,muq;
    double krho,komg,lam,lamp,lamw,mrho,momg,murho,muomg;
    //double* pm0, pm1,pgva01,pgva11,pgva00,pgva10,pgvb01,pgvb11,pgvb00,pgvb10,pgvls1,pgvt1,pgvso21,pgvls0,pgvt0,pgvso20;//size 12
    const int lpot;//! THIS SPECIFIES THE POTENTIAL
//        1=malfliet-tjon V
//        2=reid v8
//        3=urbana v14
//        4=argonne v8 (v8' reduction of v14)
//        5=argonne v9 (less simplified v14)
//        6=argonne v14
//        7=argonne v18-csbl (experimental cd/csb-large)
//        8=argonne v18-csbs (experimental cd/csb-small)
//        9=argonne v18
//       10=argonne v8'
//       11=argonne v6'
//       12=argonne v4'
//       13=argonne v2'
//       14=argonne v1'
//       15=argonne vx'
//       16=argonne v18p  (p**2 terms)
//       17=argonne v18pq (p**2 terms & l**2 tensor)
//       18=super-soft core(c) v14
//       19=super-soft core(c) v8' modified
//       20=paris
//       21=argonne v18 v1.9
//       22=argonne v8' v1.9
//       23=argonne v18 v1.7
//       24=argonne v8' v1.7
//       lpot=112 means LO, R0=1.0
//       lpot=113 means NLO, R0=1,0
//       lpot=114 means N2LO, R0=1.0
//       lpot=122 means LO, R0=1.2
//       lpot=123 means NLO, R0=1,2
//       lpot=124 means N2LO, R0=1.2

    double ftpec,pimass1,pimass2,pimass3,wsrange,tnr;
    bool v7cut;

    double yc(const double& t, const double& x);
    double yt(const double& t, const double& x);
    double yls(const double& t, const double& x);
    double yl2(const double& t, const double& x);
    double pc(const double& t);
    double pt(const double& t);
    double pls(const double& t);
    double pso2(const double& t);



    //parameter for nucleon mass variation
    //double xmn;
    //void SetXmn(const double& val);
    //double GetXmn();

    //parameter for pion mass variation in OPE (Argonne only)
    //double gam;
    //void SetGam(const double& val);
    //double GetGam();

    //parameter for pion mass variation in TPE-s (Argonne only)
    //double rho;
    //void SetRho(const double& val);
    //double GetRho();

    //parameter for pion mass variation in TPE-L (Argonne only)
    //double chi;
    //void SetChi(const double& val);
    //double GetChi();

    //parameter for heavy-meson mass variation (Argonne only)
    //double omg;
    //void SetOmg(const double& val);
    //double GetOmg();

    //parameter for intermediate coupling variation (Argonne only)
    //double ftp;
    //void SetFtp(const double& val);
    //double GetFtp();


    //double h2m;
    //double h2mcsb;
    //double tnr;
    //double hc;
    //double mp;
    //double mn;

    double* pm1;
    double* pm0;
    double* pgva01;
    double* pgva11;
    double* pgvb01;
    double* pgvb11;
    double* pgvls1;
    double* pgvt1;
    double* pgvso21;
    double* pgva00;
    double* pgva10;
    double* pgvb00;
    double* pgvb10;
    double* pgvls0;
    double* pgvt0;
    double* pgvso20;

    const double* CurrentRadius;
    const double* CurrentSmallRadius;
    const double* CurrentVerySmallRadius;

    double p01;
    double p10;
    double p00;
    double p11;
    double pt1;
    double pls1;
    double pl211;
    double pls21;
    double pt0;
    double pls0;
    double pl210;
    double pls20;
    double pl201;
    double pl200;
    double pq0;
    double rm0;
    double rm1;
    double vc0;
    double vc1;
    double vt0;
    double vt1;
    double vls0;
    double vls1;
    double vso20;
    double vso21;

    double uVar;
    double u1Var;
    double u2Var;
    double u3Var;
    double xVar;
    double x0Var;
    double x1Var;
    double x2Var;
    double x3Var;
    double xsVar;
    double xqVar;
    double xcVar;
    double yVar;
    double y1Var;
    double y2Var;
    double y3Var;
    double y4Var;
    double y6Var;
    double y7Var;
    double yrVar;
    double ypVar;
    double ydpVar;
    double hrVar;
    double tVar;
    double tpVar;
    double tdpVar;
    double tpip;
    double tpidp;
    double tpi2p;
    double tpi2dp;
    double expws;
    double wsdp;
    double wsp0;
    double wsm;
    double wsmp;
    double wsmdp;
    double wsx2p;
    double wsx2dp;

    double cpi;
    double ypi;
    double tpi;
    double tpi2;
    double tpi3;
    double rcut;
    double rcutp;
    double rcutdp;
    double ws;
    double wp;
    double ypib;
    double tpib;
    double ypid;
    double tpid;
    double pifac;
    double rws;
    double aws;
    double pl2av;
    double fsq;
    double aiws;
    double tpis;
    double tpiq;
    double ypi0;
    double ypic;
    double tpi0;
    double tpic;
    double ews;
    double ews0;
    double ws0;
    double fws;
    double wsp;
    double wsz;
    double wspp;
    double wszz;
    double wspdp;
    double wspds;
    double wsx;
    double wsx2;
    double dypi00;
    double dypic0;
    double ypi0p;
    double ypicp;
    double ypibar;
    double tpibar;
    double p11pp;
    double p11np;
    double p11nn;
    double pt1pp;
    double pt1np;
    double pt1nn;
    double p01pp;
    double p01np;
    double p01nn;
    double p11cd;
    double p11cs;
    double pt1cd;
    double pt1cs;
    double p01cd;
    double p01cs;
    double pl2p11;
    double pl2p10;
    double pl2dp11;
    double pl2dp10;
    double pl2p01;
    double pl2p00;
    double pl2dp01;
    double pl2dp00;
    double psq11;
    double psq10;
    double psq01;
    double psq00;
    double psqp11;
    double psqp10;
    double psqp01;
    double psqp00;
    double psqdp11;
    double psqdp10;
    double psqdp01;
    double psqdp00;
    double vcsb;
    double deltam;
    double crho;
    double grho;
    double fmsqi;
    double xrho;
    double zpic;
    double zcut;
    double zrho;
    double expcut;
    double rr4;
    double rc4;
    double rc6;
    double EvalValue;

//for the chiral potentials
    double x_not, x_plus, g_r0p, del_r;
    double pre_not, pre_plus;
    double W_C, V_S, V_T;
    double V_C, W_S, W_T;
    double lap_res, other_res;
    double sp_not, sp_plus;
    double ten_not, ten_plus;

    double r_0p;
    double LO_C_s;
    double LO_C_t;
    double LO_CIB;
    double LO_CSB;
    double NLO_C_s;
    double NLO_C_t;
    double* NLO_C;
    double NLO_CIB;
    double NLO_CSB;
    double N2LO_C_s;
    double N2LO_C_t;
    double* N2LO_C;
    double N2LO_CIB;
    double N2LO_CSB;




    void MalflietTjon();
    void ReidV8();
    void UrbanaV14();
    void ArgonneV14();
    void ArgonneV18();
    void ArgonneV18p();
    void SuperSoftCoreV14();
    void Paris();

    void setpot();
    void pot(const double& rr, const double& lr=0);

    void gauleg(double& x1, double& x2, double* x, double* w, const int& n);
    void spectralNLO(const double& rad);
    void spectralN2LO(const double& rad);
    void cheft(const double& rr);



public:
    DLM_StefanoPotentials(const int& WhichPot,
            const int& XMN=1, const int& GAM=1, const int& RHO=1, const int& CHI=1, const int& OMG=1, const int& FTP=1);
    ~DLM_StefanoPotentials();

    //for pp -> I==1; I3==1
    //for pn -> I==1,0; I3==0
    //for nn -> I==1; I3==-1
    double Eval(const double& Radius, const int& Spin, const int& Spin3, const int& IsoSpin, const int& IsoSpin3);
    double EvalCATS_v1_0(const double& Radius, const int& Spin);
    void PotentialName(const int& WhichPot, char* Name);


};

#endif

