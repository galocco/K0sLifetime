
#ifndef MCselector_h
#define MCselector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TDatabasePDG.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TF1.h>
#include "TRandom3.h"
#include <vector>

/// Headers needed by this particular selector
#include "MiniV0.h"
#include "MCparticle.h"
#include "TH2.h"
#include "TH3.h"


class MCselector : public TSelector {

public:
  TTreeReader     fReader;      //! the tree reader
  TTree          *fChain = 0;   //! pointer to the analyzed TTree or TChain

  /// Readers to access the data
  TTreeReaderValue<float>  fMultiplicity = {fReader, "fMultiplicity"};
  TTreeReaderArray<Lifetimes::MiniV0> V0s = {fReader, "V0s"};
  TTreeReaderArray<Lifetimes::MCparticle> MCparticles = {fReader, "MCparticles"};

  /// standard selector methods
  MCselector(TTree * /*tree*/ =0);
  virtual ~MCselector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();


  /// output file name
  string fOutputFileName;


  /// Control histograms to monitor the filtering
  TH1D* fHistV0radius;              // V0 decay vertex radius
  TH1D* fHistV0pt;                  // V0 transverse momentum
  TH1D* fHistV0eta;                 // V0 pseudorapidity
  TH2D* fHistInvMassK0s;            // Invariant mass for K0s
  TH2D* fHistInvMassLambda;         // Invariant mass for (anti-)Lambda
  TH1D* fHistDistOverTotMom;        // L/p
  TH1D* fHistV0CosPA;               // V0 cosine of pointing angle
  TH1D* fHistChi2V0;                // V0 fit chi2
  TH1D* fHistDcaNeg2PrimaryVertex;  // DCA of the negative prong to the PV
  TH1D* fHistDcaPos2PrimaryVertex;  // DCA of the positive prong to the PV
  TH1D* fHistDcaV0daughters;        // DCA between the two prongs
  TH1D* fHistV0armAlpha;            // Armenteros alpha
  TH1D* fHistV0armAlpha2;            // Armenteros alpha
  TH1D* fHistV0armPt;               // Armenteros pt
  TH1D* fHistLeastNxedRows;         // Min number of xed roads
  TH1D* fHistLeastXedOverFindable;  // Min number of xed roads/findable clusters
  TH1D* fHistMaxChi2PerCluster;     // Max chi2 per cluster in TPC
  TH1D* fHistNsigmaPosPion;         // # sigma TPC pion for the positive prong
  TH1D* fHistNsigmaPosProton;       // # sigma TPC proton for the positive prong
  TH1D* fHistNsigmaNegPion;         // # sigma TPC pion for the negative prong
  TH1D* fHistNsigmaNegProton;       // # sigma TPC proton for the negative prong
  TH1D* fHistEtaPos;                // Pseudorapidity of the positive prong
  TH1D* fHistEtaNeg; // Pseudorapidity of the negative prong
  TH2D* fHistArment; 
  TH2D* fHistInvMassK0sCt;
  TH2D* fHistInvMassLambdaCt;
  TH2D* fHistAntiLambda;
  TH2D* fHistLambda;
  
  TH1D* fHistCtV0;//dist over top
  TH1D* fHistCtMC;
  TH2D* fHistDeltaCtVsCtV0;//diff su ct
  TH2D* fHistDeltaCtVsCtMC;
  TH2D* fHistDeltaXBinCenter;
  TH1D* fHistCtMCTot;
  TH1D* fHistMother;
  TH1D* fHistSecWeakDecay;
  TH1D* fHistSecMaterial;
  TH2D* fHistEff;
  TH1D* fHistCtMCPrimary;// istogrammi per vedere se cambia l'andamento
  TH1D* fHistCtMCSecondary;
  TH2D* fHistCtRecVsCtGen;

  TDatabasePDG pdg;
  ClassDef(MCselector,0);
};

#endif
