#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "DLM_ScatteringParameters.h"
#include "DLM_Source.h"
#include "DLM_Random.h"
#include "DLM_Potentials.h"

#include "CATS.h" 
#include "CATStools.h"
#include "CATSconstants.h"


void DLM_PotSp::set(double p0_, double p1_, double f_, double d_, long entries_){
    p0 = p0_;
    p1 = p1_;
    f = f_;
    d = d_;
    num_entries = entries_;
    if(p0>p0_max){p0_max = p0;}
    if(p0<p0_min){p0_min = p0;}
    if(p1>p1_max){p1_max = p1;}
    if(p1<p1_min){p1_min = p1;}
}

bool DLM_PotSp::operator+=(const DLM_PotSp& other){
    p0 += other.p0;
    p1 += other.p1;
    f += other.f;
    d += other.d;
    num_entries += other.num_entries;

    if(other.p0_max > p0_max) p0_max = other.p0_max;
    if(other.p0_min < p0_min) p0_min = other.p0_min;
    if(other.p1_max > p1_max) p1_max = other.p1_max;
    if(other.p1_min < p1_min) p1_min = other.p1_min;

    return true;
}

bool DLM_PotSp::operator/=(const double& value){
    p0 /= value;
    p1 /= value;
    f /= value;
    d /= value;
    num_entries = lround(double(num_entries)/value);
    return true;
}

bool DLM_PotSp::operator=(const DLM_PotSp& other){
    p0 = other.p0;
    p1 = other.p1;
    f = other.f;
    d = other.d;
    num_entries = other.num_entries;
    p0_max = other.p0_max;
    p0_min = other.p0_min;
    p1_max = other.p1_max;
    p1_min = other.p1_min;
    return true;
}
bool DLM_PotSp::operator=(const double& value){
    p0 = value;
    p1 = value;
    f = value;
    d = value;
    num_entries = 0;
    p0_min = 1e64;
    p0_max = -1e64;
    p1_min = 1e64;
    p1_max = -1e64;
    return true;
}
DLM_PotSp DLM_PotSp::operator*(const double& value){
    DLM_PotSp Result;
    Result.p0 = p0 * value;
    Result.p1 = p1 * value;
    Result.f = f * value;
    Result.d = d * value;
    Result.num_entries = num_entries;  
    return Result;
}
void DLM_PotSp::compute_average(){
    *this /= num_entries;
}

DLM_ScatteringPars::DLM_ScatteringPars():kMin(0),kMax(80),kSteps(5),EPS(5e-10),eps_f(0.01),eps_d(0.1){

    id_p0 = 0;
    id_p1 = 1;

    target_f[0] = 0;
    target_f[1] = 0;
    target_d[0] = 0;
    target_d[1] = 0;    
    red_mass = 0;
    //NumThreads = 0;

    Kitty = NULL;
    pot_pars = NULL;
    num_pot_pars = 0;
    dlm_PotSp_Map = NULL;
    dlm_PotSp_AvgMap = NULL;

    rangen = NULL;

}

DLM_ScatteringPars::~DLM_ScatteringPars(){
    if(pot_pars){
        delete [] pot_pars;
        pot_pars = NULL;
    }
    if(Kitty){
        delete Kitty;
        Kitty = NULL;
    }
    if(dlm_PotSp_Map){
        delete dlm_PotSp_Map;
        dlm_PotSp_Map = NULL;
    }
    if(dlm_PotSp_AvgMap){
        delete dlm_PotSp_AvgMap;
        dlm_PotSp_AvgMap = NULL;
    }
}
void DLM_ScatteringPars::Reset(){
    if(Kitty){
        delete Kitty;
        Kitty = NULL;
    }
    if(dlm_PotSp_Map){
        delete dlm_PotSp_Map;
        dlm_PotSp_Map = NULL;
    }
    if(dlm_PotSp_AvgMap){
        delete dlm_PotSp_AvgMap;
        dlm_PotSp_AvgMap = NULL;
    }
}

void DLM_ScatteringPars::SetPotFun(PotFunction f, unsigned p0, unsigned p1){
    Reset();
    potential = f;
    id_p0 = p0;
    id_p1 = p1;
    if(!pot_pars){
        pot_pars = new double [2];
        num_pot_pars = 2;
    }
}
void DLM_ScatteringPars::SetPotPar(std::vector<double> pot_pars_){
    Reset();
    if(pot_pars){
        delete [] pot_pars;
    }
    pot_pars = new double [pot_pars_.size()];
    num_pot_pars = pot_pars_.size();
    for(int iPar=0; iPar<pot_pars_.size(); iPar++){
        pot_pars[iPar] = pot_pars_.at(iPar);
    }
}

void DLM_ScatteringPars::SetParLimits(unsigned WhichPar, double par_min, double par_max){
    Reset();
    if(WhichPar==id_p0){
        pot_par0[0] = par_min;
        pot_par0[1] = par_max;       
    }
    else if(WhichPar==id_p1){
        pot_par1[0] = par_min;
        pot_par1[1] = par_max;       
    }
}
void DLM_ScatteringPars::SetParGrid(unsigned WhichPar, unsigned num_bins){
    Reset();
    if(WhichPar==id_p0){
        bins_par0 = num_bins;
    }
    else if(WhichPar==id_p1){
        bins_par1 = num_bins;    
    }    
}

//the desired scattering length f, and effective range d
void DLM_ScatteringPars::SetTarget_f(double min_f, double max_f){
    Reset();
    target_f[0] = min_f;
    target_f[1] = max_f;
}
void DLM_ScatteringPars::SetSlGrid(unsigned num_bins){
    Reset();
    bins_f = num_bins;
}
void DLM_ScatteringPars::SetTarget_d(double min_d, double max_d){
    Reset();
    target_d[0] = min_d;
    target_d[1] = max_d;    
}
void DLM_ScatteringPars::SetErGrid(unsigned num_bins){
    Reset();
    bins_d = num_bins;
}

void DLM_ScatteringPars::SetRedMass(double red_mass_){
    Reset();
    red_mass = red_mass_;
}

void DLM_ScatteringPars::SetRandomSeed(unsigned seed){
    if(rangen){delete rangen;}
    rangen = new DLM_Random(seed);
}

bool DLM_ScatteringPars::GetScatteringParameters(double& f, double& fe, double& d, double& de){
    if(!Kitty) return false;
    if(!Kitty->CkStatus()) return false;

    double f0 = 0;
    double fe0 = 0;
    double d0 = 0;
    double de0 = 0;
    double f_val;
    double d_val;
    unsigned prm = 0;
    for(unsigned iKstar=0; iKstar<kSteps; iKstar++){
        double ki = Kitty->GetMomentum(iKstar);
        double gi = ki/(tan(Kitty->GetPhaseShift(iKstar,0,0)));
        for(unsigned jKstar=iKstar+1; jKstar<kSteps; jKstar++){
            double kj = Kitty->GetMomentum(jKstar);
            double gj = kj/(tan(Kitty->GetPhaseShift(jKstar,0,0)));
            f_val = (kj*kj-ki*ki)/(gi*kj*kj-gj*ki*ki) * hbarc;
            d_val = 2.*(gj-gi)/(kj*kj-ki*ki) * hbarc;
            f0 += f_val;
            fe0 += f_val*f_val;
            d0 += d_val;
            de0 += d_val*d_val;
            prm++;
        }
    }
    f0 /= double(prm);
    d0 /= double(prm);
    fe0 /= double(prm); 
    de0 /= double(prm); 
    if(fe0-f0*f0<=0) {fe0 = 0;}
    else {fe0 = sqrt(fe0-f0*f0);}
    if(de0-d0*d0<=0) {de0 = 0;}
    else {de0 = sqrt(de0-d0*d0);}


    f = 0;
    fe = 0;
    d = 0;
    de = 0;
    prm = 0;
    for(unsigned iKstar=0; iKstar<kSteps; iKstar++){
        double ki = Kitty->GetMomentum(iKstar);
        double gi = ki/(tan(Kitty->GetPhaseShift(iKstar,0,0)));
        //printf("ki, gi = %f %f\n",ki,gi);
        for(unsigned jKstar=iKstar+1; jKstar<kSteps; jKstar++){
            double kj = Kitty->GetMomentum(jKstar);
            double gj = kj/(tan(Kitty->GetPhaseShift(jKstar,0,0)));
            f_val = (kj*kj-ki*ki)/(gi*kj*kj-gj*ki*ki) * hbarc;
            d_val = 2.*(gj-gi)/(kj*kj-ki*ki) * hbarc;
            //printf("  errors %.3e %.3e\n",fabs(f_val-f0)/fe0, fabs(d_val-d0)/de0);
            if(fabs(f_val-f0)/fe0<1 && fabs(d_val-d0)/de0<1){
                //printf(" !!!! \n");
                f += f_val;
                fe += f_val*f_val;
                d += d_val;
                de += d_val*d_val;
                prm++;
            }
        }
    }


    f /= double(prm);
    d /= double(prm);
    fe /= double(prm); 
    de /= double(prm); 
    if(fe-f*f<=0) {fe = 0;}
    else {fe = sqrt(fe-f*f);}
    if(de-d*d<=0) {de = 0;}
    else {de = sqrt(de-d*d);}

    if(fabs(f)>=1 && fabs(fe/f)>eps_f) return false;    
    if(fabs(f)<1 && fabs(fe)>eps_f) return false;

    if(fabs(d)>=1 && fabs(de/d)>eps_d) return false;    
    if(fabs(d)<1 && fabs(de)>eps_d) return false;

    if(f<target_f[0] || f>target_f[1]) return false;
    if(d<target_d[0] || d>target_d[1]) return false;

    return true;
}


//completely random sampling within the limits of the potential parameters
void DLM_ScatteringPars::SampleSomeStuff(unsigned NumSamples, double min_p0, double max_p0, double min_p1, double max_p1){
    if(!Kitty){
        Kitty = new CATS();
        Kitty->SetMomBins(kSteps,kMin,kMax);
        Kitty->SetThetaDependentSource(false);
        CATSparameters cPars(CATSparameters::tSource, 1, true);
        cPars.SetParameter(0, 1.0);
        Kitty->SetAnaSource(GaussSource, cPars);
        Kitty->SetUseAnalyticSource(true);
        Kitty->SetMomentumDependentSource(false);
        Kitty->SetExcludeFailedBins(false);
        Kitty->SetQ1Q2(0);
        Kitty->SetQuantumStatistics(false);
        Kitty->SetNumChannels(1);
        Kitty->SetRedMass(red_mass);
        Kitty->SetNumPW(0, 1);
        Kitty->SetSpin(0, 0);
        Kitty->SetChannelWeight(0, 1.);
        CATSparameters pPars(CATSparameters::tPotential, num_pot_pars, true);
        pPars.SetParameters(pot_pars);
        Kitty->SetShortRangePotential(0, 0, potential, pPars);
        Kitty->SetNotifications(CATS::nWarning);
        Kitty->SetEpsilonProp(EPS);
    }
    if(!rangen){
        rangen = new DLM_Random(0);
    }
    if(!dlm_PotSp_Map){
        dlm_PotSp_Map = new DLM_Histo<DLM_PotSp> ();
        dlm_PotSp_Map->SetUp(2);
        dlm_PotSp_Map->SetUp(0,bins_f,target_f[0],target_f[1]);
        dlm_PotSp_Map->SetUp(1,bins_d,target_d[0],target_d[1]);
        dlm_PotSp_Map->Initialize();
    }

    double rnd_p0;
    double rnd_p1;
    for(unsigned uIter=0; uIter<NumSamples; uIter++){
        rnd_p0 = rangen->Uniform(min_p0, max_p0);
        rnd_p1 = rangen->Uniform(min_p1, max_p1);
        //printf("r pars = %f %f (%i %i)\n",rnd_p0,rnd_p1,id_p0,id_p1);
        Kitty->SetShortRangePotential(0, 0, id_p0, rnd_p0);
        Kitty->SetShortRangePotential(0, 0, id_p1, rnd_p1);
        Kitty->KillTheCat();
        double f0,d0;
        double fe0,de0;
        if( GetScatteringParameters(f0,fe0,d0,de0) ){
            DLM_PotSp potsp;
            potsp.set(rnd_p0,rnd_p1,f0,d0);
            dlm_PotSp_Map->AddAt(f0,d0,potsp);
            //printf("%f %f\n",f0,d0);
        }
        else{
            uIter--;
        }
    }
}

void DLM_ScatteringPars::RandomScan(unsigned NumSamples){
    if(dlm_PotSp_AvgMap){
        delete dlm_PotSp_AvgMap;
        dlm_PotSp_AvgMap = NULL;
    }

    SampleSomeStuff(NumSamples, pot_par0[0], pot_par0[1], pot_par1[0], pot_par1[1]);
}

//we sample target_fraction of the time by trying to sample around pot pars that already 
//provided a (f,d) combo within the desired limits. The remaining 1-target_fraction we sample with RandomScan
void DLM_ScatteringPars::TargetedScan(unsigned NumSamples, float target_fraction){
//printf("TargetedScan\n");
    if(dlm_PotSp_AvgMap){
        delete dlm_PotSp_AvgMap;
        dlm_PotSp_AvgMap = NULL;
    }

    unsigned samples = 0;
    while(samples<NumSamples){
        if(dlm_PotSp_Map && rangen->Uniform(0,1)<target_fraction){
            for(unsigned uTry=0; uTry<10; uTry++){
                int rnd_bin = rangen->Integer(0, dlm_PotSp_Map->GetNbins());
                if(dlm_PotSp_Map->GetBinContent(rnd_bin).num_entries>=3){
                    SampleSomeStuff(1,  dlm_PotSp_Map->GetBinContent(rnd_bin).p0_min,
                                        dlm_PotSp_Map->GetBinContent(rnd_bin).p0_max,
                                        dlm_PotSp_Map->GetBinContent(rnd_bin).p1_min,
                                        dlm_PotSp_Map->GetBinContent(rnd_bin).p1_max);
                    samples++;
                    break;
                }
            }
        }
        else{
            SampleSomeStuff(1, pot_par0[0], pot_par0[1], pot_par1[0], pot_par1[1]);
            samples++;
        }    
    }
}

std::vector<double> DLM_ScatteringPars::GetPotPars(double f0, double d0){

    std::vector<double> result(num_pot_pars,0);

    if(!dlm_PotSp_AvgMap){
        if(!dlm_PotSp_Map) return result;
        dlm_PotSp_AvgMap = new DLM_Histo<DLM_PotSp> (*dlm_PotSp_Map);
        for(unsigned uBin=0; uBin<dlm_PotSp_AvgMap->GetNbins(); uBin++){
            dlm_PotSp_AvgMap->GetBinElement(uBin).compute_average();
        }
    }

    if( dlm_PotSp_AvgMap->GetBinContent(dlm_PotSp_AvgMap->FindBin(0,f0), dlm_PotSp_AvgMap->FindBin(1,d0)).num_entries==0){
        return result;
    }
    //dlm_PotSp_AvgMap->GetBinContent(dlm_PotSp_AvgMap->FindBin(0,f0),dlm_PotSp_AvgMap->FindBin(0,f0));

    for(unsigned uPar=0; uPar<num_pot_pars; uPar++){
        if(!dlm_PotSp_AvgMap) break;
        if(uPar==id_p0) result[uPar]=(dlm_PotSp_AvgMap->Eval2D(f0,d0).p0);
        else if(uPar==id_p1) result[uPar]=(dlm_PotSp_AvgMap->Eval2D(f0,d0).p1);
        else result[uPar]=(pot_pars[uPar]);
    }
    return result;

}

double DLM_ScatteringPars::Occupancy(){
    unsigned filled_bins = 0;
    for(unsigned uBin=0; uBin<dlm_PotSp_Map->GetNbins(); uBin++){
        filled_bins += bool(dlm_PotSp_Map->GetBinContent(uBin).num_entries);
        //printf("   %d \n", dlm_PotSp_Map->GetBinContent(uBin).num_entries);
    }
    return double(filled_bins)/double(dlm_PotSp_Map->GetNbins());
}
unsigned DLM_ScatteringPars::GetNumEntries(){
    if(!dlm_PotSp_Map) return 0;
    unsigned nentr=0;
    for(unsigned uBin=0; uBin<dlm_PotSp_Map->GetNbins(); uBin++){
        nentr += dlm_PotSp_Map->GetBinContent(uBin).num_entries;
    }
    return nentr;
}

void DLM_ScatteringPars::Save(std::string file_name, bool Overwrite){
    if(dlm_PotSp_Map){
        dlm_PotSp_Map->QuickWrite(file_name.c_str(), Overwrite);
    }
}
void DLM_ScatteringPars::SaveSettings(std::string file_name, bool Overwrite){

// Determine the open mode based on the Overwrite flag
    // std::ios::trunc clears the file, std::ios::app appends to the end
    std::ios_base::openmode mode = Overwrite ? std::ios::trunc : std::ios::app;
    
    std::ofstream outFile(file_name, mode);

    if (!outFile.is_open()) {
        // You might want to use your class's error reporting system here
        return;
    }

    // Set formatting to scientific notation with 6 decimal places (%.6e)
    outFile << std::scientific << std::setprecision(6);

    // --- Potential Parameters ---
    outFile << "num_pot_pars\t" << num_pot_pars << "\n";
    if (pot_pars != nullptr) {
        for (unsigned i = 0; i < num_pot_pars; ++i) {
            outFile << "pot_pars[" << i << "]\t" << pot_pars[i] << "\n";
        }
    }

    // --- Particle IDs ---
    outFile << "id_p0\t\t" << id_p0 << "\n";
    outFile << "id_p1\t\t" << id_p1 << "\n";

    // --- Parameter 0 Bins and Range ---
    outFile << "bins_par0\t" << bins_par0 << "\n";
    outFile << "pot_par0[0]\t" << pot_par0[0] << "\n";
    outFile << "pot_par0[1]\t" << pot_par0[1] << "\n";

    // --- Parameter 1 Bins and Range ---
    outFile << "bins_par1\t" << bins_par1 << "\n";
    outFile << "pot_par1[0]\t" << pot_par1[0] << "\n";
    outFile << "pot_par1[1]\t" << pot_par1[1] << "\n";

    // --- Target F Bins and Range ---
    outFile << "bins_f\t\t" << bins_f << "\n";
    outFile << "target_f[0]\t" << target_f[0] << "\n";
    outFile << "target_f[1]\t" << target_f[1] << "\n";

    // --- Target D Bins and Range ---
    outFile << "bins_d\t\t" << bins_d << "\n";
    outFile << "target_d[0]\t" << target_d[0] << "\n";
    outFile << "target_d[1]\t" << target_d[1] << "\n";

    // --- Physics Constants ---
    outFile << "red_mass\t" << red_mass << "\n";

    outFile.close();

}


void DLM_ScatteringPars::Load(std::string file_name, bool reset){
    if(reset) Reset();
    if(!dlm_PotSp_Map){
        dlm_PotSp_Map = new DLM_Histo<DLM_PotSp> ();
        dlm_PotSp_Map->SetUp(2);
        dlm_PotSp_Map->SetUp(0,bins_f,target_f[0],target_f[1]);
        dlm_PotSp_Map->SetUp(1,bins_d,target_d[0],target_d[1]);
        dlm_PotSp_Map->Initialize();

    }
    DLM_Histo<DLM_PotSp> temp_hist;
    if(!temp_hist.QuickLoad(file_name.c_str())){
        //Reset();
    }
    else{
        if(reset){
            *dlm_PotSp_Map = temp_hist;
        }
        else{
            for(unsigned uBin=0; uBin<dlm_PotSp_Map->GetNbins(); uBin++){
                dlm_PotSp_Map->GetBinElement(uBin) += temp_hist.GetBinElement(uBin);
            }
        }
    }
}