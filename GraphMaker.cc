//macro to draw some useful graphics
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TDatabasePDG.h>
#include <iostream>
#include <TAxis.h>
#include <vector>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
using namespace std;
void GraphMaker()
{ 
  TFile * Eff1 = new TFile("Results/MCselectorchecks.root");
  TH2D * histEff = (TH2D*) Eff1->Get("fHistDeltaCtVsCtMC");
  TFile results("Graph.root", "RECREATE");
  results.cd();
  histEff->Write();
  TFile * Eff = new TFile("FitK0s.root");
  TCanvas *cv = new TCanvas("cv","cv");
  TH1D * fC[7];
  cv->cd();
  for(int Rad=7;Rad>=1;Rad--)
  {
    fC[Rad] = (TH1D*) Eff->Get(Form("histEff%i",Rad));
    fC[Rad]->Draw("SAME PMC PLC");
  }
  cv->BuildLegend();
  results.cd();
  cv->Write();

  TCanvas *cvK0s = new TCanvas("cvK0s","cvK0s");
  for(int Rad=7;Rad>=1;Rad--)
  {
    fC[Rad] = (TH1D*) Eff->Get(Form("HistoNCtCor%s Rad: %i","K0s",Rad));
    fC[Rad]->SetMarkerStyle(20);
    fC[Rad]->Draw("SAME PMC PLC");
  }
  results.cd();
  cvK0s->Write();

  TFile * MC = new TFile("Results/MCselectorchecks.root");
  TH2D * histPrim = (TH2D*) MC->Get("fHistCtMC");
  TH2D * histSec = (TH2D*) MC->Get("fHistCtMCSecondary");

  histPrim->SetTitle("primary");
  histPrim->SetLineColor(kBlue);
  histPrim->SetLineWidth(2);
  histSec->SetLineWidth(2);
  histSec->SetTitle("secondary");
  histSec->SetLineColor(kRed);
  TCanvas *cvMC = new TCanvas("cvMc","cvMc");
  TPave* StarUpsSys = new TPave(-0.1,2,0.1,2,1,"tbrl");;
  histPrim->Draw();
  histSec->Draw("SAME");
  results.cd();
  cvMC->Write();
}