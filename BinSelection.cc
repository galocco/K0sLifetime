

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
//the function return the bin with the max value between the two double
int FindMax(TH1D*,double,double);
//the function return the binwidth according to our "taste"
double BinAdaptation(double,double,double);
//Times is how much times multiply the RMS to define the binwidth,binwidth is larger than
//risolution because we want to be sure that the ct misured really stay in the bin
//BinUnit is the pace of cm in ct (just to make the histogram look better)
//CtLimit is the max ct that we want to consider
void BinSelection(const double Times=6,const double BinUnit=0.25,double CtLimit)
{
  TFile * FileRMS = new TFile("Results/ResolutionFit.root");
  //HistoRMS is the histogram with the RMS
  TH1D * HistoRMS = (TH1D*) FileRMS->Get("histoRMS");
  double WidthCandidate;
  double past=0;
  int counter=0;
  vector<double> BinWidth;
  vector<double> BinWidthConverted;
  //the binwidth is N-times the larger RMS in the ct interval
  for(int Bin=1;past<CtLimit;Bin++)
  {
    counter++;
    WidthCandidate=6*HistoRMS->GetBinContent(Bin);
    if(WidthCandidate>=FindMax(HistoRMS,past,past+WidthCandidate))
    {
      past+=WidthCandidate;
      Bin=HistoRMS->GetBin(past);
      BinWidth.push_back(WidthCandidate);
      BinWidthConverted.push_back(BinAdaptation(WidthCandidate/6,Times,BinUnit));
    }
  }
  TH1D * HistoBW = new TH1D ("HistoBW",";bin;binwidth (cm)",counter,0,counter);
  for(int BinBW=1;BinBW<=counter;BinBW++)
    HistoBW->SetBinContent(BinBW,BinWidthConverted[BinBW-1]);
  TFile results("Results/BinWidth.root", "RECREATE");
  HistoBW->Draw();
  HistoBW->Write();
}

int FindMax(TH1D * Histo,double Min,double Max)
{
  int BinMax=Histo->GetBinContent(Histo->GetBin(Min));
  for(int Bin=Histo->GetBin(Min)+1;Bin<=Histo->GetBin(Max);Bin++)
    if(Histo->GetBinContent(Bin)>BinMax)
      BinMax=Histo->GetBinContent(Bin);
  return BinMax;
}

double BinAdaptation(double BinWidth,double Times,double BinUnit)
{
  if(Times*BinWidth-BinUnit*static_cast<int>(Times*BinWidth/BinUnit)>BinUnit/2)
    return (static_cast<int>(Times*BinWidth/BinUnit)+1)*BinUnit;
  else
    return (static_cast<int>(Times*BinWidth/BinUnit)*BinUnit);
}