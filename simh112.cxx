#include "simh112.h"
ClassImp(simh112)

simh112::simh112(char*  files)
{
  Input=new TChain("h112");
  Input->Add(files);
  Init(Input);
}

simh112::~simh112()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t simh112::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t simh112::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void simh112::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("sume", sume, &b_sume);
   fChain->SetBranchAddress("sump", sump, &b_sump);
   fChain->SetBranchAddress("sums", sums, &b_sums);
   fChain->SetBranchAddress("Nhit", &Nhit, &b_Nhit);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("hit", hit, &b_hit);
   fChain->SetBranchAddress("trk", trk, &b_trk);
   fChain->SetBranchAddress("Ntrack", &Ntrack, &b_Ntrack);
   fChain->SetBranchAddress("ptrk", ptrk, &b_ptrk);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("pdg", pdg, &b_pdg);
   fChain->SetBranchAddress("mo", mo, &b_mo);
   fChain->SetBranchAddress("da", da, &b_da);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("m", m, &b_m);
   fChain->SetBranchAddress("v", v, &b_v);
   fChain->SetBranchAddress("gpid", gpid, &b_gpid);
   fChain->SetBranchAddress("gtrk", gtrk, &b_gtrk);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("proc", &proc, &b_proc);
   fChain->SetBranchAddress("shat", &shat, &b_shat);
   fChain->SetBranchAddress("that", &that, &b_that);
   fChain->SetBranchAddress("uhat", &uhat, &b_uhat);
   fChain->SetBranchAddress("cost", &cost, &b_cost);
   fChain->SetBranchAddress("q2", &q2, &b_q2);
   fChain->SetBranchAddress("in", in, &b_in);
   fChain->SetBranchAddress("xin", xin, &b_xin);
   fChain->SetBranchAddress("x", x, &b_x);
   Notify();
}

Bool_t simh112::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void simh112::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
