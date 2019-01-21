//macro just to run all easily 
#include <TProof.h>
#include <TChain.h>
void MCrun() {
  TProof::Open("");
  TChain V0tree("fTreeV0");
  V0tree.Add("V0treeMC.root");
  V0tree.SetProof();
  V0tree.Process("MCselector.cxx+");
}