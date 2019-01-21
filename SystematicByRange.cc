#include <iostream>
#include <TAxis.h>
#include <vector>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
using namespace std;
float NParticles(TF1* peak,TF1 * noise,TH1D * histo, float NSigma)
{
  float min=peak->GetParameter(1)-peak->GetParameter(2)*NSigma;
  float max=peak->GetParameter(1)+peak->GetParameter(2)*NSigma;
  return (peak->Integral(min,max,0.1)-noise->Integral(min,max,0.1))/histo->GetBinWidth(1);
}
void SystematicByRange(const int ParticleCode=310,const char FileName[]="MassvsLifetimeschecks.root",const char HistoName[]="fHistInvMassK0sCt" ,const char ParticleName[]="K0s",const char ResultsFile[]="Systematic.root",const char EffHisto[]="fHistEff")
{

  TDatabasePDG pdg;
  TFile results(ResultsFile, "RECREATE");

  TFile * Eff = new TFile("Results/MCselectorchecks.root");
  TH1D * histEff = (TH1D*) Eff->Get("fHistEff");
  gStyle->SetOptStat(0);
  TFile *histos = new TFile(FileName);
  TH3D *HistoInvMassTot = (TH3D*) histos->Get(HistoName);
  TAxis *ctAxis = HistoInvMassTot->GetXaxis();

  //binning
  double Binning[HistoInvMassTot->GetNbinsX()+1];
  Binning[0]=ctAxis->GetBinLowEdge(1);
  for(int Bin=1;Bin<=HistoInvMassTot->GetNbinsX();Bin++)
    Binning[Bin]=ctAxis->GetBinUpEdge(Bin);

  TH1D *HistoNCtParticle = new TH1D(Form("HistoNCt%s",ParticleName),";#it{c}t (cm);#frac{dN}{dct} (cm^{-1})",HistoInvMassTot->GetNbinsX(),Binning);

  int LowerWindowBin;
  float LowEdge,UpEdge;
  float MinPeak,MaxPeak,BinMax;
  int counter=0;
  vector<TH1D> histoSigma; 
  for(double NSigma=3;NSigma<=5;NSigma+=0.25)
  {
  TH1D HistoNCtParticle = TH1D(Form("HistoNCt%s %.2f sigma ",ParticleName,NSigma),Form("Bin counting %.2f#sigma ;#it{c}t (cm);#frac{dN}{dct} (cm^{-1})",NSigma),HistoInvMassTot->GetNbinsX(),Binning);
  counter++;
  for(int BinCt=1;BinCt<=HistoInvMassTot->GetNbinsX();BinCt++)
  {
    TF1 *peak = new TF1("peak","gaus(0)");
    TF1 *noise = new TF1("noise","expo(0)");
    noise->SetLineColor(kBlue);    

    TH1D* histoInvMass=(TH1D*) HistoInvMassTot->ProjectionY(Form("%sBin%i",ParticleName,BinCt),BinCt,BinCt,1,7);
    histoInvMass->SetTitle(Form("%.2f #leq #it{c}t < %.2f cm", ctAxis->GetBinLowEdge(BinCt), ctAxis->GetBinUpEdge(BinCt)));
    BinMax=1;
    for(int BinMass=2;BinMass<=histoInvMass->GetNbinsX();BinMass++)
    {
      if(histoInvMass->GetBinContent(BinMass)>=histoInvMass->GetBinContent(BinMax))
        BinMax=BinMass;
    }
    MinPeak=histoInvMass->GetBinLowEdge(BinMax-3);
    MaxPeak=histoInvMass->GetBinWidth(BinCt)*7+MinPeak;
    peak->SetRange(MinPeak,MaxPeak); 
    histoInvMass->Fit(peak,"QLR");

    LowEdge=peak->GetParameter(1)-NSigma*peak->GetParameter(2);
    UpEdge=peak->GetParameter(1)+NSigma*peak->GetParameter(2);
  
    peak->SetRange(histoInvMass->GetBinLowEdge(1),histoInvMass->GetBinLowEdge(histoInvMass->GetNbinsX())+histoInvMass->GetBinWidth(histoInvMass->GetNbinsX()));
    vector<float> RemovedBins;
    for(int Binj=histoInvMass->FindBin(peak->GetParameter(1)-peak->GetParameter(2)*NSigma);Binj<=histoInvMass->FindBin(peak->GetParameter(1)+peak->GetParameter(2)*NSigma);Binj++)
    {
      RemovedBins.push_back(histoInvMass->GetBinContent(Binj));
      histoInvMass->SetBinContent(Binj,0);
    }
    float min=peak->GetParameter(1)-peak->GetParameter(2)*NSigma;
    float max=peak->GetParameter(1)+peak->GetParameter(2)*NSigma;
    LowerWindowBin=histoInvMass->FindBin(peak->GetParameter(1)-peak->GetParameter(2)*NSigma);
    if((max-min)>=(histoInvMass->GetNbinsX()*histoInvMass->GetBinWidth(1)))
    {
      for(auto Bin : RemovedBins)
        histoInvMass->SetBinContent(LowerWindowBin++,Bin);
      continue;
    }
    histoInvMass->Fit(noise,"Q+");

    for(auto Bin : RemovedBins)
      histoInvMass->SetBinContent(LowerWindowBin++,Bin);

    double sigPlusBkg = histoInvMass->Integral(histoInvMass->FindBin(min),histoInvMass->FindBin(max));
    double sig = sigPlusBkg - noise->Integral(min,max,0.1)/histoInvMass->GetBinWidth(1);
    HistoNCtParticle.SetBinContent(BinCt,sig/HistoNCtParticle.GetBinWidth(BinCt));
    HistoNCtParticle.SetBinError(BinCt,TMath::Sqrt(sigPlusBkg)/HistoNCtParticle.GetBinWidth(BinCt));

  }
    histoSigma.push_back(HistoNCtParticle);
  }
  //istogramma con RMS dei conteggi per bin 
  TH1D *HistoRMSRel = new TH1D("HistoRMSRel","relative systematic uncertainty",HistoInvMassTot->GetNbinsX(),Binning);
  HistoRMSRel->GetXaxis()->SetTitle("ct (cm)");
  HistoRMSRel->GetYaxis()->SetTitle("standard deviation / counts");//controllare
  TH1D *HistoRMS = new TH1D("HistoRMS","systematic error",HistoInvMassTot->GetNbinsX(),Binning);
  HistoRMS->GetXaxis()->SetTitle("ct (cm)");
  HistoRMS->GetYaxis()->SetTitle("RMS");//controllare

  for(int BinCt=1;BinCt<=20;BinCt++)
  {
    vector <double> mis;
    for(int cont=0;cont<=8;cont++)
      mis.push_back(histoSigma[cont].GetBinContent(BinCt));
    HistoRMSRel->SetBinContent(BinCt,TMath::RMS(mis.begin(),mis.end())/histoSigma[2].GetBinContent(BinCt));
    HistoRMS->SetBinContent(BinCt,TMath::RMS(mis.begin(),mis.end()));
  }
  results.cd();
  HistoRMS->Write();
  HistoRMSRel->Write();

  ////
  TCanvas *cv= new TCanvas("cv","cv");
  for(int cont=0;cont<=8;cont++)
  {
    histoSigma[cont].SetTitle("");
    histoSigma[cont].Draw("SAME PLC PMC");
  }
  cv->BuildLegend();
  results.cd();
  cv->Write();
  
}