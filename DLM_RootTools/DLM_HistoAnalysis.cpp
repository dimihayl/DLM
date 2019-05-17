

#include "TMath.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1D.h"

#include "DLM_HistoAnalysis.h"

using namespace TMath;


TH1D * MakeHistoFromTF1(TF1 * f1, double min, double max, unsigned int nbins){
    TH1D * h1d = new TH1D("h1d", "h1d", nbins, min, max);

    double x,y;

    for(unsigned int i=1; i<=nbins; i++){
        x = h1d->GetBinCenter(i);
        y = f1->Eval(x);
        h1d->SetBinContent(i, y);
    }

    return h1d;
}

double FunctionVar(TF1 * f1, double min, double max, double* mean, double* mean2, unsigned int steps){

    TH1D * h1d = MakeHistoFromTF1(f1, min, max, steps);

    double Mean = 0;
    double Mean2 = 0;

    Mean = 0;
    Mean2 = 0;

    double hintegral = h1d->Integral();

    for(unsigned int i=1; i<=steps; i++){
        Mean += h1d->GetBinCenter(i)*h1d->GetBinContent(i);
        Mean2  += Power(h1d->GetBinCenter(i), 2.)*h1d->GetBinContent(i);
    }
    Mean /= hintegral;
    Mean2  /= hintegral;

    if(mean && mean2){
        *mean = Mean;
        *mean2 = Mean2;
    }

    delete h1d;
    return Mean2 - Power(Mean, 2.);
}



void GetCentralInterval(TF1 * f1, const double& min, const double& max, double alpha, double& from, double& to, const unsigned int& steps){

    from = min;
    to = max;

    TH1D * h1d = MakeHistoFromTF1(f1, min, max, steps);

    double hintegral = h1d->Integral();
    double fintegral = f1->Integral(min, max);
    double myinterval = fintegral*(0.5-alpha/2.);
    double currentintegral = 0;

    h1d->Scale(fintegral/hintegral);

    for(unsigned int i=1; i<=steps; i++){
        currentintegral += h1d->GetBinContent(i);
        if(currentintegral>=myinterval){
            from = h1d->GetBinCenter(i);
            break;
        }
    }

    currentintegral = 0;
    for(unsigned int i=steps; i>0; i--){
        currentintegral += h1d->GetBinContent(i);
        if(currentintegral>=myinterval){
            to = h1d->GetBinCenter(i);
            break;
        }
    }

    delete h1d;

}



double GetCentralInterval(const TH1& h1, const double alpha, double& from, double& to, bool extrapolate){

    double hintegral = h1.Integral();
    double myinterval = hintegral*(0.5-alpha/2.);
    double currentintegral = 0;

    //Printf("hintegral = %f", hintegral);
    //Printf("myinterval = %f", myinterval);

    //h1->Scale(1/hintegral);
    int nbins = h1.GetNbinsX();

    for(int i=1; i<=nbins; i++){

        currentintegral += h1.GetBinContent(i);
        //double x2 = currentintegral;
        if(currentintegral>=myinterval){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f1 = currentintegral - h1.GetBinContent(i);
                double f2 = currentintegral;
                double fI = myinterval;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                from = (fI-b)/a;
                //Printf("x1 = %f", x1);
                //Printf("x2 = %f", x2);
                //Printf("a = %f", a);
                //Printf("b = %f", b);
                //Printf("myinterval = %f", myinterval);
            }
            else
                from = h1.GetBinLowEdge(i);
            break;
        }
    }
    currentintegral = 0;
    for(unsigned int i=nbins; i>0; i--){
        currentintegral += h1.GetBinContent(i);
        if(currentintegral>=myinterval){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f2 = hintegral-(currentintegral - h1.GetBinContent(i));
                double f1 = hintegral-currentintegral;
                double fI = hintegral-myinterval;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                to = (fI-b)/a;
            }
            else
                to = h1.GetBinLowEdge(i) + h1.GetBinWidth(i);
            break;
        }
    }

    currentintegral = 0;
    for(int i=1; i<=nbins; i++){
        currentintegral += h1.GetBinContent(i);
        if(currentintegral>=0.5*hintegral){
            if(extrapolate){
                double x1 = h1.GetBinLowEdge(i);
                double x2 = h1.GetXaxis()->GetBinUpEdge(i);
                double f1 = currentintegral - h1.GetBinContent(i);
                double f2 = currentintegral;
                double fI = 0.5*hintegral;
                double a = (f2-f1)/(x2-x1);
                double b = f1-a*x1;
                return (fI-b)/a;
            }
            return h1.GetBinLowEdge(i);
        }
    }

    return 0;
}

void DrawCentralInterval(TH1F*& h1, double Median, double LowReach, double UpReach, int NumNewBins,
                         TH1F*& h1main, TH1F*& h1bLow, TH1F*& h1bLow2, TH1F*& h1bUp, TH1F*& h1bUp2, TH1F*& h1median){

    //double from = h1->GetBinLowEdge(h1->FindFirstBinAbove(0));
    //double to = h1->GetXaxis()->GetBinUpEdge(h1->FindLastBinAbove(0));
    double from = h1->GetBinLowEdge(1);
    double to = h1->GetXaxis()->GetBinUpEdge(h1->GetNbinsX());
    int NumBins = h1->GetNbinsX();
    double NewBinRange = (to-from)/double(NumNewBins);
    //from -= NewBinRange;
    //to += NewBinRange;

    Printf("NumBins = %i", NumBins);
    Printf("from = %f", from);
    Printf("to = %f", to);

    h1main = new TH1F(TString::Format("%s_vis", h1->GetName()),
                          TString::Format("%s_vis", h1->GetName()), NumNewBins, from, to);

    int FirstNonZeroBin = h1main->FindBin(LowReach);
    //if(LowReach > h1main->GetBinCenter(FirstNonZeroBin)) FirstNonZeroBin++;
    int LastNonZeroBin = h1main->FindBin(UpReach);
    //if(UpReach < h1main->GetBinCenter(LastNonZeroBin)) LastNonZeroBin--;
    double FNZB_low = h1main->GetBinLowEdge(FirstNonZeroBin);
    double LNZB_up = h1main->GetXaxis()->GetBinUpEdge(LastNonZeroBin);

    //Printf("FirstNonZeroBin = %i; LastNonZeroBin = %i",FirstNonZeroBin, LastNonZeroBin);
    //Printf("%f : %f : %f : %f : %f : %f",from, FNZB_low, LowReach, UpReach, LNZB_up, to);


    h1bLow = new TH1F(TString::Format("%s_vis2l", h1->GetName()),
                          TString::Format("%s_vis2l", h1->GetName()), FirstNonZeroBin-1,
                             from, FNZB_low);
    h1bUp = new TH1F(TString::Format("%s_vis2u", h1->GetName()),
                          TString::Format("%s_vis2u", h1->GetName()), NumNewBins - LastNonZeroBin,
                            LNZB_up, to);
    h1bLow2 = new TH1F(TString::Format("%s_vis2l2", h1->GetName()),
                          TString::Format("%s_vis2l2", h1->GetName()), 1,
                              FNZB_low, LowReach);
    h1bUp2 = new TH1F(TString::Format("%s_vis2u2", h1->GetName()),
                          TString::Format("%s_vis2u2", h1->GetName()), 1,
                            UpReach, LNZB_up);


    h1median = new TH1F(TString::Format("%s_vis3", h1->GetName()),
                          TString::Format("%s_vis3", h1->GetName()), 1, Median-NewBinRange/8., Median+NewBinRange/8.);

    h1main->SetFillColor(h1->GetFillColor());
    h1main->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    h1main->GetYaxis()->SetTitle(h1->GetYaxis()->GetTitle());


    h1bLow->SetFillColor(kGray+1);
    h1bUp->SetFillColor(kGray+1);
    h1bLow2->SetFillColor(kGray+1);
    h1bUp2->SetFillColor(kGray+1);

    h1median->SetFillColor(kBlack);


    for(int i=1; i<=NumBins; i++){
        h1main->Fill(h1->GetBinCenter(i), h1->GetBinContent(i));
        //h1b->Fill(h1->GetBinCenter(i+1), h1->GetBinContent(i+1));
    }

    for(int i=1; i<FirstNonZeroBin; i++)
        h1bLow->SetBinContent(i, h1main->GetBinContent(i));
    h1bLow2->SetBinContent(1, h1main->GetBinContent(FirstNonZeroBin));
    h1bUp2->SetBinContent(1, h1main->GetBinContent(LastNonZeroBin));
    for(int i=LastNonZeroBin+1; i<=NumNewBins; i++)
        h1bUp->SetBinContent(i-LastNonZeroBin, h1main->GetBinContent(i));



    h1main->SetStats(false);
    h1main->SetTitle("");
    h1main->SetLineColor(h1main->GetFillColor());
    h1main->SetLineWidth(1);
    h1main->SetMarkerStyle(0);
    h1main->SetMarkerSize(0);
    h1main->GetYaxis()->SetRange(0, 1.05*h1main->GetMaximum());
    h1main->GetYaxis()->SetRangeUser(0, 1.05*h1main->GetMaximum());
    h1main->GetYaxis()->SetLimits(0, 1.05*h1main->GetMaximum());
    h1main->GetXaxis()->SetLabelSize(0.065);
    h1main->GetXaxis()->CenterTitle();
    h1main->GetXaxis()->SetTitleOffset(1.15);
    h1main->GetXaxis()->SetLabelOffset(0.02);
    h1main->GetXaxis()->SetTitleSize(0.075);
    h1main->GetYaxis()->SetLabelSize(0.065);
    h1main->GetYaxis()->CenterTitle();
    h1main->GetYaxis()->SetTitleOffset(1);
    h1main->GetYaxis()->SetTitleSize(0.075);
    h1main->GetYaxis()->SetNdivisions(505);

    h1bLow->SetStats(false);
    h1bLow->SetTitle("");
    h1bLow->SetLineColor(h1bLow->GetFillColor());
    h1bLow->SetLineWidth(1);
    h1bLow->SetMarkerStyle(0);
    h1bLow->SetMarkerSize(0);

    h1bLow2->SetStats(false);
    h1bLow2->SetTitle("");
    h1bLow2->SetLineColor(h1bLow2->GetFillColor());
    h1bLow2->SetLineWidth(1);
    h1bLow2->SetMarkerStyle(0);
    h1bLow2->SetMarkerSize(0);

    h1bUp->SetStats(false);
    h1bUp->SetTitle("");
    h1bUp->SetLineColor(h1bUp->GetFillColor());
    h1bUp->SetLineWidth(1);
    h1bUp->SetMarkerStyle(0);
    h1bUp->SetMarkerSize(0);

    h1bUp2->SetStats(false);
    h1bUp2->SetTitle("");
    h1bUp2->SetLineColor(h1bUp2->GetFillColor());
    h1bUp2->SetLineWidth(1);
    h1bUp2->SetMarkerStyle(0);
    h1bUp2->SetMarkerSize(0);

    h1median->SetStats(false);
    h1median->SetTitle("");
    h1median->SetLineColor(h1median->GetFillColor());
    h1median->SetLineWidth(1);
    h1median->SetMarkerStyle(0);
    h1median->SetMarkerSize(0);

    //for(int i=FirstNonZeroBin; i<=LastNonZeroBin; i++){
    //    h1b->SetBinContent(i,0);
    //}

    h1median->SetBinContent(1, h1main->GetBinContent(h1main->FindBin(Median)));

}
