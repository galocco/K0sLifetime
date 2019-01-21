//macro to study the resolution on ct then find the right binning on ct
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
using namespace std;
void DeltaCtFit()
{
  ///
  TFile *histos = new TFile("Results/MCselectorchecks.root");
  TFile results("Results/ResolutionFit.root","RECREATE");
  //fHistDeltaCtVsCtMC is the histogram of the difference between recontructed K0s ct and 
  //the Monte Carlo truth
  TH2D * histoTot= (TH2D*) histos->Get("fHistDeltaCtVsCtMC");
  int BinMax,Max;
  float MinPeak,MaxPeak;
  TAxis* ctAxis = histoTot->GetXaxis();
  //fHistDeltaCtVsCtMC bin must have all bins in ct with the same binwidth
  TH1D* histoRMS = new TH1D("histoRMS","#Deltact RMS(ct)",histoTot->GetNbinsX(),0,histoTot->GetNbinsX());
  histoRMS->GetXaxis()->SetTitle("ct (cm)")
  histoRMS->GetYaxis()->SetTitle("standard deviation (cm)");
  //loop to make slice on ct and to get the RMS of each slice
  for(int BinCt=1;BinCt<=histoTot->GetNbinsX();BinCt++)
  {
    TH1D* fHistBin=(TH1D*) histoTot->ProjectionY(Form("#Delta#it{c}t Bin%i",BinCt),BinCt,BinCt); 
    fHistBin->SetTitle(Form("%.2f #leq #it{c}t < %.2f cm", ctAxis->GetBinLowEdge(BinCt), ctAxis->GetBinUpEdge(BinCt)));
    fHistBin->GetYaxis()->SetTitle("counts");
    BinMax=1;
    for(int BinMass=2;BinMass<=fHistBin->GetNbinsX();BinMass++)
    {
      if(fHistBin->GetBinContent(BinMass)>=fHistBin->GetBinContent(BinMax))
      {
        BinMax=BinMass;
        Max=fHistBin->GetBinContent(BinMass);
      }
    }
    histoRMS->SetBinContent(BinCt,fHistBin->GetRMS());
    histoRMS->SetBinError(BinCt,fHistBin->GetRMSError());
    results.cd();
    fHistBin->Write();
  }
  histoRMS->Write();
  results.Close();
}