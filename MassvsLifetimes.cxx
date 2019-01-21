//Selector of the real data
//Run with MvLrun.cc
#define MassvsLifetimes_cxx


#include "MassvsLifetimes.h"
#include <TH2.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TF1.h>
#include "MiniV0.h"
#include <vector>
using namespace Lifetimes;

MassvsLifetimes::MassvsLifetimes(TTree *) :
fChain{nullptr},
fOutputFileName{"MassvsLifetimeschecks.root"} {}
//______________________________________________________________________________
//______________________________________________________________________________



void MassvsLifetimes::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

//______________________________________________________________________________
//______________________________________________________________________________


void MassvsLifetimes::SlaveBegin(TTree * /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  ////
  //Setting of the binning parameter
  int Bin=20;
  double BinMax=20;  
  //Binwidths determined by the macro BinSelection,the array must be fill by hand
  double Binning[] ={0.,0.75,0.75,0.75,1.25,1.75,2.25,2.25,2.25,2.25,2.25,2.25,2.25};
  //distance of each layer, the last (400.) is number big enough to represent all ALICE
  double Radious[]={0.,3.9,7.6,15.,23.9,38.,43.,400.};
  int RadiousLenght=7;
  //loop to obtain an array of the binnig
  double sum=0;
  for(int Cell=0;Cell<=12;Cell++)
  {
    sum+=Binning[Cell];
    Binning[Cell]=sum;
  }
  //parameter of the Invariant Mass vs ct
  double MBin[81+8];
  MBin[0]=0.45;
  for(int Cell=1;Cell<=80+8;Cell++)
    MBin[Cell]=MBin[Cell-1]+0.10/80;
  
  fHistV0radius = new TH1D("fHistV0radius", ";V0 r (cm); Counts", 250, 0, 250);
  fHistV0pt =new TH1D("fHistV0pt", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 4.);
  fHistV0eta = new TH1D("fHistV0eta", ";V0 #eta; Counts", 80, -0.8, 0.8);
  fHistInvMassK0s = new TH2D("fHistInvMassK0s",";V0 #it{p}_{T} (GeV/#it{c}); m_{#pi#pi} (GeV/#it{c}^{2})", 20,0, 2, 80, 0.46, 0.54);
  fHistInvMassLambda = new TH2D("fHistInvMassLambda",";V0 #it{p}_{T} (GeV/#it{c}); m_{p#pi} (GeV/#it{c}^{2})", 20, 0,2, 80, 1.075, 1.155);
  fHistDistOverTotMom = new TH1D("fHistDistOverTotMom", ";V0 L/#it{p} (#it{c} cm / GeV); Counts",250, 0, 250);
  fHistV0CosPA = new TH1D("fHistV0CosPA", ";V0 cos(#theta_{P}); Counts",MiniV0::fgkV0cosPA_n, MiniV0::fgkV0cosPA_f,MiniV0::fgkV0cosPA_l);
  fHistChi2V0 =new TH1D("fHistChi2V0", ";V0 #chi^{2}; Counts", MiniV0::fgkV0chi2_n,MiniV0::fgkV0chi2_f, MiniV0::fgkV0chi2_l);
  fHistDcaNeg2PrimaryVertex =new TH1D("fHistDcaNeg2PrimaryVertex", ";Neg prong DCA (cm); Counts",MiniV0::fgkDCAProng2PV_n, MiniV0::fgkDCAProng2PV_f,MiniV0::fgkDCAProng2PV_l);
  fHistDcaPos2PrimaryVertex =new TH1D("fHistDcaPos2PrimaryVertex", ";Pos prong DCA (cm); Counts",MiniV0::fgkDCAProng2PV_n, MiniV0::fgkDCAProng2PV_f,MiniV0::fgkDCAProng2PV_l);
  fHistDcaV0daughters = new TH1D("fHistDcaV0daughters", ";Prongs DCA; Counts",MiniV0::fgkDCAProngs_n, MiniV0::fgkDCAProngs_f,MiniV0::fgkDCAProngs_l);
  fHistV0armAlpha = new TH1D("fHistV0armAlpha", ";Armenteros #alpha; Counts",MiniV0::fgkArmAlpha_n, MiniV0::fgkArmAlpha_f,MiniV0::fgkArmAlpha_l);
  fHistV0armAlpha2 = new TH1D("fHistV0armAlpha2", ";Armenteros #alpha; Counts",MiniV0::fgkArmAlpha_n, -1.1,1.1);
  fHistV0armPt = new TH1D("fHistV0armPt", ";Armenteros #it{p}_{T} (GeV/#it{c}); Counts",MiniV0::fgkArmPt_n, MiniV0::fgkArmPt_f, MiniV0::fgkArmPt_l);
  fHistLeastNxedRows = new TH1D("fHistLeastNxedRows", ";Min # of crossed rows; Counts", 256, -0.5, 255.5);
  fHistLeastXedOverFindable = new TH1D("fHistLeastXedOverFindable",";Min # of crossed rows / findable clusters; Counts",MiniV0::fgkXedOverFindable_n, MiniV0::fgkXedOverFindable_f,MiniV0::fgkXedOverFindable_l);
  fHistMaxChi2PerCluster =new TH1D("fHistMaxChi2PerCluster", ";Min #chi^{2}/TPC clusters; Counts",MiniV0::fgkChi2xCluster_n, MiniV0::fgkChi2xCluster_f,MiniV0::fgkChi2xCluster_l);
  fHistNsigmaPosPion = new TH1D("fHistNsigmaPosPion", ";n_{#sigma} TPC Pos Pion; Counts",MiniV0::fgkTPCsigma_n, MiniV0::fgkTPCsigma_f, MiniV0::fgkTPCsigma_l);
  fHistNsigmaPosProton = new TH1D("fHistNsigmaPosProton", ";n_{#sigma} TPC Pos Proton; Counts",MiniV0::fgkTPCsigma_n, MiniV0::fgkTPCsigma_f, MiniV0::fgkTPCsigma_l);
  fHistNsigmaNegPion = new TH1D("fHistNsigmaNegPion", ";n_{#sigma} TPC Neg Pion; Counts",MiniV0::fgkTPCsigma_n, MiniV0::fgkTPCsigma_f, MiniV0::fgkTPCsigma_l);
  fHistNsigmaNegProton = new TH1D("fHistNsigmaNegProton", ";n_{#sigma} TPC Neg Proton; Counts",MiniV0::fgkTPCsigma_n, MiniV0::fgkTPCsigma_f, MiniV0::fgkTPCsigma_l);
  fHistEtaPos = new TH1D("fHistEtaPos", ";Pos prong #eta; Counts",MiniV0::fgkEta_n, MiniV0::fgkEta_f, MiniV0::fgkEta_l);
  fHistEtaNeg = new TH1D("fHistEtaNeg", ";Neg prong #eta; Counts",MiniV0::fgkEta_n, MiniV0::fgkEta_f, MiniV0::fgkEta_l);
  fHistArment = new TH2D("fHistArmenteros",";Armenteros #alpha;Armenteros #it{p}_{T} (GeV/#it{c})",MiniV0::fgkArmAlpha_n, MiniV0::fgkArmAlpha_f,MiniV0::fgkArmAlpha_l,MiniV0::fgkArmPt_n, MiniV0::fgkArmPt_f, MiniV0::fgkArmPt_l);
  fHistInvMassK0sCt = new TH3D("fHistInvMassK0sCt",";#it{c}t (cm); m_{#pi#pi} (GeV/#it{c}^{2})",12,Binning, 80,MBin,RadiousLenght,Radious);
  fHistInvMassLambdaCt = new TH2D("fHistInvMassLambdaCt",";#it{c}t (cm); m_{#pi#pi} (GeV/#it{c}^{2})",60,0,30, 80, 1.075, 1.155);
  fHistAntiLambda = new TH2D("fHistAntiLambda","#it{c}t (cm);m_{#pi#pi} (GeV/#it{c}^{2}",60,0,30,80,1.075,1.555);
  fHistLambda = new TH2D("fHistLambda","#it{c}t (cm);m_{#pi#pi} (GeV/#it{c}^{2}",60,0,30,80,1.075,1.155);

  GetOutputList()->Add(fHistV0radius);
  GetOutputList()->Add(fHistV0pt);
  GetOutputList()->Add(fHistV0eta);
  GetOutputList()->Add(fHistInvMassK0s);
  GetOutputList()->Add(fHistInvMassLambda);
  GetOutputList()->Add(fHistDistOverTotMom);
  GetOutputList()->Add(fHistV0CosPA);
  GetOutputList()->Add(fHistChi2V0);
  GetOutputList()->Add(fHistDcaNeg2PrimaryVertex);
  GetOutputList()->Add(fHistDcaPos2PrimaryVertex);
  GetOutputList()->Add(fHistDcaV0daughters);
  GetOutputList()->Add(fHistV0armAlpha);
  GetOutputList()->Add(fHistV0armPt);
  GetOutputList()->Add(fHistLeastNxedRows);
  GetOutputList()->Add(fHistLeastXedOverFindable);
  GetOutputList()->Add(fHistMaxChi2PerCluster);
  GetOutputList()->Add(fHistNsigmaPosPion);
  GetOutputList()->Add(fHistNsigmaPosProton);
  GetOutputList()->Add(fHistNsigmaNegPion);
  GetOutputList()->Add(fHistNsigmaNegProton);
  GetOutputList()->Add(fHistEtaPos);
  GetOutputList()->Add(fHistEtaNeg);
  GetOutputList()->Add(fHistArment);
  GetOutputList()->Add(fHistInvMassK0sCt);
  GetOutputList()->Add(fHistInvMassLambdaCt);
  GetOutputList()->Add(fHistLambda);
  GetOutputList()->Add(fHistAntiLambda);
  ////

}

//______________________________________________________________________________
//______________________________________________________________________________


Bool_t MassvsLifetimes::Process(Long64_t entry) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  for (auto mini : V0s) {

    fHistV0radius->Fill(mini.GetV0radius());
    fHistV0pt->Fill(mini.GetV0pt());
    fHistV0eta->Fill(mini.GetV0eta());
    fHistInvMassK0s->Fill(mini.GetV0pt(), mini.GetK0sInvMass());
    fHistInvMassLambda->Fill(mini.GetV0pt(),mini.GetLambdaInvMass());
    fHistDistOverTotMom->Fill(mini.GetDistOverP());
    fHistV0CosPA->Fill(mini.GetV0CosPA());
    fHistChi2V0->Fill(mini.GetV0chi2());
    fHistDcaNeg2PrimaryVertex->Fill(mini.GetNegProngPvDCA());
    fHistDcaPos2PrimaryVertex->Fill(mini.GetPosProngPvDCA());
    fHistDcaV0daughters->Fill(mini.GetProngsDCA());//non sicuro
    fHistV0armAlpha->Fill(mini.GetArmenterosAlpha());
    fHistV0armPt->Fill(mini.GetArmenterosPt());
    fHistLeastNxedRows->Fill(mini.GetLeastNumberOfXedRows());
    fHistLeastXedOverFindable->Fill(mini.GetLeastXedRowsOverFindable());
    fHistMaxChi2PerCluster->Fill(mini.GetMaxChi2perCluster());
    fHistNsigmaPosPion->Fill(mini.GetPosProngTPCnsigmaPion());
    fHistNsigmaPosProton->Fill(mini.GetPosProngTPCnsigmaProton());
    fHistNsigmaNegPion->Fill(mini.GetNegProngTPCnsigmaPion());
    fHistNsigmaNegProton->Fill(mini.GetNegProngTPCnsigmaProton());
    fHistEtaPos->Fill(mini.GetPosProngEta());
    fHistEtaNeg->Fill(mini.GetNegProngEta());
    fHistArment->Fill(mini.GetArmenterosAlpha(),mini.GetArmenterosPt());
    fHistInvMassK0sCt->Fill(pdg.GetParticle(310)->Mass()*mini.GetDistOverP(),mini.GetK0sInvMass(),mini.GetV0radius());
    fHistInvMassLambdaCt->Fill(pdg.GetParticle(3122)->Mass()*mini.GetDistOverP(),mini.GetLambdaInvMass());
    if(mini.GetArmenterosAlpha()<0)
      fHistAntiLambda->Fill(pdg.GetParticle(3122)->Mass()*mini.GetDistOverP(),mini.GetLambdaInvMass());
    else
      fHistLambda->Fill(pdg.GetParticle(3122)->Mass()*mini.GetDistOverP(),mini.GetLambdaInvMass());  
  }
  
  return kTRUE;
}

//______________________________________________________________________________
//______________________________________________________________________________


void MassvsLifetimes::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}

//______________________________________________________________________________
//______________________________________________________________________________


void MassvsLifetimes::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TFile output(Form("Results/%s", fOutputFileName.data()),"RECREATE");
  GetOutputList()->Write();
  output.Close();

}

//______________________________________________________________________________
//_______________________________________________________________________________


void MassvsLifetimes::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);

}

//______________________________________________________________________________
//_______________________________________________________________________________


Bool_t MassvsLifetimes::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}
