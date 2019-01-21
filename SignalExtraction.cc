//macro to do the signal extraction cutting on the distance between decay vertex(V0Radius) 
//and beamline vs ct
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
/*
ParticleCode is the PDG of the particle
NSigma how much times the sigma define the extraction range
*/
void SignalExtraction(const int ParticleCode=310,const float NSigma=5,const char FileName[]="Results/MassvsLifetimeschecks.root",const char HistoName[]="Results/fHistInvMassK0sCt" ,const char ParticleName[]="K0s",const char ResultsFile[]="Results/FitK0s.root",const char EffHisto[]="fHistEff",const bool Draw=true)
{
  TDatabasePDG pdg;
  TFile * Eff = new TFile("Results/MCselectorchecks.root");
  TH2D * histEff = (TH2D*) Eff->Get(EffHisto);
  //histogram where see lifetime vs V0Radius
  TH1D * histoLifeTimes = new TH1D("histoLifeTimes","",7,0,7);
  TFile results(ResultsFile, "RECREATE"); 
  double TauFit[7],Radious[7],ErrTauFit[7],ErrRadious[7];
  //loop over the cut on distance between decay vertex and beamline
  for(int Rad=1;Rad<=7;Rad++)
  {
    TH1D  *histEffProj = (TH1D*) histEff->ProjectionX(Form("histEff%i",Rad),1,Rad);
    if(Rad!=7)
      histEffProj->SetTitle(Form("V0 radius < %.1f cm",histEff->GetYaxis()->GetBinUpEdge(Rad)));
    else
      histEffProj->SetTitle("All V0s");
    TH1D * fHistPrimary = (TH1D*) Eff->Get("fHistCtMCPrimary");

    histEffProj->GetYaxis()->SetTitle("efficiency x acceptance");
    //reconstructed K0s divided by all to obtain efficiency
    histEffProj->Divide(fHistPrimary);
    gStyle->SetOptStat(0);
    TFile *histos = new TFile(FileName);
    TH3D *HistoInvMassTot = (TH3D*) histos->Get(HistoName);
    TAxis *ctAxis = HistoInvMassTot->GetXaxis();

    double Binning[HistoInvMassTot->GetNbinsX()+1];
    Binning[0]=ctAxis->GetBinLowEdge(1);
    for(int Bin=1;Bin<=HistoInvMassTot->GetNbinsX();Bin++)
      Binning[Bin]=ctAxis->GetBinUpEdge(Bin);
    //histogram of the counts vs ct
    TH1D *HistoNCtParticle = new TH1D(Form("HistoNCt%s Rad: %i",ParticleName,Rad),Form("V0Radius < %.1f cm;#it{c}t (cm);#frac{dN}{dct} (cm^{-1})",histEff->GetYaxis()->GetBinUpEdge(Rad)),HistoInvMassTot->GetNbinsX(),Binning);
    //histogram of the counts vs ct corrected by efficiency
    TH1D *HistoNCtParticleCor = new TH1D(Form("HistoNCtCor%s Rad: %i",ParticleName,Rad),Form("V0Radius < %.1f cm;#it{c}t (cm);#frac{dN}{dct} (cm^{-1})",histEff->GetYaxis()->GetBinUpEdge(Rad)),HistoInvMassTot->GetNbinsX(),Binning);
    if(Rad==7)
    {
      HistoNCtParticle->SetTitle("All V0s");
      HistoNCtParticleCor->SetTitle("All V0s");
    }
    int LowerWindowBin;
    float LowEdge,UpEdge;
    float MinPeak,MaxPeak,BinMax;
    //signal extraction from every bin of ct
    for(int BinCt=1;BinCt<=HistoInvMassTot->GetNbinsX();BinCt++)
    {
      TF1 *peak = new TF1("peak","gaus(0)");
      TF1 *noise = new TF1("noise","expo(0)");
        noise->SetLineColor(kViolet);    
      TH2D *histoInvMass = (TH2D*) HistoInvMassTot->ProjectionY(Form("%sBin%i Rad: %i",ParticleName,BinCt,Rad),BinCt,BinCt,1,Rad);
      histoInvMass->SetTitle(Form("%.2f #leq #it{c}t < %.2f cm", ctAxis->GetBinLowEdge(BinCt), ctAxis->GetBinUpEdge(BinCt)));
      histoInvMass->GetYaxis()->SetTitle("counts");
      //What I want to do is find the extraction range fitting the peak 
      BinMax=1;
      for(int BinMass=2;BinMass<=histoInvMass->GetNbinsX();BinMass++)
      {
        if(histoInvMass->GetBinContent(BinMass)>=histoInvMass->GetBinContent(BinMax))
          BinMax=BinMass;
      }
      MinPeak=histoInvMass->GetBinLowEdge(BinMax-3);
      MaxPeak=histoInvMass->GetBinWidth(BinCt)*7+MinPeak;
      peak->SetRange(MinPeak,MaxPeak); 
      histoInvMass->Fit(peak,"QLRN");
      LowEdge=peak->GetParameter(1)-NSigma*peak->GetParameter(2);
      UpEdge=peak->GetParameter(1)+NSigma*peak->GetParameter(2);
    
      peak->SetRange(histoInvMass->GetBinLowEdge(1),histoInvMass->GetBinLowEdge(histoInvMass->GetNbinsX())+histoInvMass->GetBinWidth(histoInvMass->GetNbinsX()));
      histoInvMass->GetListOfFunctions()->Add(peak);
      
      TLine *LowerBorder = new TLine(LowEdge,0,LowEdge,histoInvMass->GetMaximum());
      TLine *UpperBorder = new TLine(UpEdge,0,UpEdge,histoInvMass->GetMaximum());
      UpperBorder->SetLineColor(kGreen);
      LowerBorder->SetLineColor(kGreen);
      histoInvMass->GetListOfFunctions()->Add(LowerBorder);
      histoInvMass->GetListOfFunctions()->Add(UpperBorder);

      //Now I remove the counts in the extraction range to fit the background
      vector<float> RemovedBins;
      for(int Binj=histoInvMass->FindBin(peak->GetParameter(1)-peak->GetParameter(2)*NSigma);Binj<=histoInvMass->FindBin(peak->GetParameter(1)+peak->GetParameter(2)*NSigma);Binj++)
      {
        RemovedBins.push_back(histoInvMass->GetBinContent(Binj));
        histoInvMass->SetBinContent(Binj,0);
      }
      //Refill of the histogram
      float min=peak->GetParameter(1)-peak->GetParameter(2)*NSigma;
      float max=peak->GetParameter(1)+peak->GetParameter(2)*NSigma;
      LowerWindowBin=histoInvMass->FindBin(peak->GetParameter(1)-peak->GetParameter(2)*NSigma);
      if((max-min)>=(histoInvMass->GetNbinsX()*histoInvMass->GetBinWidth(1)))
      {
        for(auto Bin : RemovedBins)
          histoInvMass->SetBinContent(LowerWindowBin++,Bin);
        results.cd();
        histoInvMass->Write();
        continue;
      }
      histoInvMass->Fit(noise,"Q+");

      for(auto Bin : RemovedBins)
        histoInvMass->SetBinContent(LowerWindowBin++,Bin);
      results.cd();
      histoInvMass->Write();
      //Calculation of the number of particles divided by ct binwidth and setting of the uncertainties
      double sigPlusBkg = histoInvMass->Integral(histoInvMass->FindBin(min),histoInvMass->FindBin(max));
      double bkg = noise->Integral(min,max,0.1)/histoInvMass->GetBinWidth(1);
      double sig = sigPlusBkg - bkg; HistoNCtParticle->SetBinContent(BinCt,sig/HistoNCtParticle->GetBinWidth(BinCt));
      HistoNCtParticle->SetBinError(BinCt,TMath::Sqrt(sigPlusBkg)/HistoNCtParticle->GetBinWidth(BinCt));
      HistoNCtParticle->SetBinContent(BinCt,sig/HistoNCtParticle->GetBinWidth(BinCt));
      histEffProj->SetBinError(BinCt,TMath::Sqrt(histEffProj->GetBinContent(BinCt)*(1-histEffProj->GetBinContent(BinCt))/fHistPrimary->GetBinContent(BinCt)));
    
      if(histEffProj->GetBinContent(BinCt)!=0)
      {
        HistoNCtParticleCor->SetBinContent(BinCt,(sig/histEffProj->GetBinContent(BinCt))/HistoNCtParticle->GetBinWidth(BinCt));
        HistoNCtParticleCor->SetBinError(BinCt,TMath::Sqrt(sigPlusBkg + bkg)/histEffProj->GetBinContent(BinCt)/HistoNCtParticle->GetBinWidth(BinCt));
      }
      else
        HistoNCtParticleCor->SetBinError(BinCt,0);
    }
    //fit of the counts vs ct to get the lifetime from the slope of the exponential
    TF1 * expo = new TF1("expo","[0]*TMath::Exp(-x /(100*TMath::C() * [1]))",0.75,7.5);
    expo->SetParameter(1,8.954e-11);
    expo->SetParameter(0,1.e5);
    HistoNCtParticleCor->Fit(expo,"R");
    float tau=expo->GetParameter(1)*pow(10,11);
    float Errtau=expo->GetParError(1)*pow(10,11);
    TauFit[Rad-1]=tau;
    Radious[Rad-1]=histEff->GetYaxis()->GetBinUpEdge(Rad);
    ErrTauFit[Rad-1]=Errtau;
    ErrRadious[Rad-1]=0;
    
    float ChiNorm=expo->GetChisquare()/expo->GetNDF();
    if(Rad==7)
    {
      TLegend* legend = new TLegend(0.2,0.9,0.5,1);
      legend->SetTextSize(0.05);
      legend->SetLineColor(kWhite);
      
      legend->AddEntry(HistoNCtParticleCor,Form("#tau = (%.2f #pm %.2f)#upoint10^{-11} s",tau,Errtau),"");
      HistoNCtParticleCor->GetListOfFunctions()->Add(legend);
    }

    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerColor(kBlue);
    histEffProj->SetMarkerColor(kBlue);
    results.cd();
    histEffProj->Write();
    HistoNCtParticle->Write();
    HistoNCtParticleCor->Write();
    histoLifeTimes->SetBinContent(Rad,tau);
    histoLifeTimes->SetBinError(Rad,Errtau);
  }
  //stuff to see the results
  if(Draw)
  {
    TCanvas *cv = new TCanvas("cv","cv");
    TGraphErrors *LifeTimes = new TGraphErrors(7,Radious,TauFit,ErrRadious,ErrTauFit);
    LifeTimes->SetMarkerStyle(20);
    LifeTimes->SetMarkerColor(kRed);
    LifeTimes->SetLineColor(kRed);
    cv->SetLogx();
    LifeTimes->GetYaxis()->SetTitle("#tau (10^{-11} s)");
    LifeTimes->GetXaxis()->SetTitle("V0 radius (cm)");
    TLine *lu= new TLine(0,8.958,440,8.958);
    TLine *PDG= new TLine(0,8.954,440,8.954);
    PDG->SetLineWidth(3);
    TLine *ld= new TLine(0,8.950,440,8.950);
    PDG->SetLineColor(kBlue);
    lu->SetLineColor(kBlue);
    ld->SetLineColor(kBlue);
    LifeTimes->Draw("AP");
    PDG->Draw("SAME");
    TLegend *leg = new TLegend(0.2,0.2,0.8,0.4);
    leg->AddEntry(LifeTimes,"#tau measured for each cut","p");
    leg->AddEntry(PDG,"PDG quoted #tau = (8.954 #pm 0.004)#upoint10^{-11} s","l");
    leg->Draw("SAME");
    results.cd();
    cv->Draw();
  }
  histoLifeTimes->Write();
}