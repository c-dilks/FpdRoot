//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug 12 22:26:13 2007 by ROOT version 5.12/00f
// from TTree h111/Title
// found on file: ../tmp/run6173032.1.root
//////////////////////////////////////////////////////////

#ifndef h111_h
#define h111_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class h111 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           event;
   UChar_t         ipre;
   UChar_t         ipost;
   UInt_t          Bchi;
   UInt_t          Bclo;
   UInt_t          Awdbm;
   UChar_t         Awtrg;
   UChar_t         Awdaq;
   UShort_t        Token;
   UInt_t          addbit;
   UInt_t          Dsm;
   UInt_t          extbusy;
   UInt_t          modbusy;
   UInt_t          physwd;
   UInt_t          Trgwd;
   UInt_t          Dsmad;
   UInt_t          cotamibusy;
   UChar_t         npre;
   UChar_t         npost;
   UInt_t          Cpa[32];
   UInt_t          qdsm[8];
   UInt_t          lastdsm[8];
   UInt_t          Bcd[16];
   UInt_t          L1[32];
   UInt_t          Vtxdsm[8];
   UInt_t          Fpddsm[8];
   UInt_t          Emcdsm[8];
   UInt_t          Bbcl1[16];
   UChar_t         Zdc[16];
   UChar_t         Zdcsmd[32];
   UChar_t         Ctb[256];
   UChar_t         Mwc[128];
   UChar_t         Bemc[128];
   UChar_t         Bbc[96];
   UChar_t         Fpdens[112];
   UChar_t         Fpdetb[64];
   UChar_t         Fpdwns[112];
   UChar_t         Fpdwtb[64];
   UInt_t          Fpdensl1[8];
   UInt_t          Fpdetbl1[8];
   UInt_t          Fpdwnsl1[8];
   UInt_t          Fpdwtbl1[8];
   UInt_t          bunchid;
   UInt_t          bunchid7bit;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_ipre;   //!
   TBranch        *b_ipost;   //!
   TBranch        *b_Bchi;   //!
   TBranch        *b_Bclo;   //!
   TBranch        *b_Awdbm;   //!
   TBranch        *b_Awtrg;   //!
   TBranch        *b_Awdaq;   //!
   TBranch        *b_Token;   //!
   TBranch        *b_addbit;   //!
   TBranch        *b_Dsm;   //!
   TBranch        *b_extbusy;   //!
   TBranch        *b_modbusy;   //!
   TBranch        *b_physwd;   //!
   TBranch        *b_Trgwd;   //!
   TBranch        *b_Dsmad;   //!
   TBranch        *b_cotamibusy;   //!
   TBranch        *b_npre;   //!
   TBranch        *b_npost;   //!
   TBranch        *b_Cpa;   //!
   TBranch        *b_qdsm;   //!
   TBranch        *b_lastdsm;   //!
   TBranch        *b_Bcd;   //!
   TBranch        *b_L1;   //!
   TBranch        *b_Vtxdsm;   //!
   TBranch        *b_Fpddsm;   //!
   TBranch        *b_Emcdsm;   //!
   TBranch        *b_Bbcl1;   //!
   TBranch        *b_Zdc;   //!
   TBranch        *b_Zdcsmd;   //!
   TBranch        *b_Ctb;   //!
   TBranch        *b_Mwc;   //!
   TBranch        *b_Bemc;   //!
   TBranch        *b_Bbc;   //!
   TBranch        *b_Fpdens;   //!
   TBranch        *b_Fpdetb;   //!
   TBranch        *b_Fpdwns;   //!
   TBranch        *b_Fpdwtb;   //!
   TBranch        *b_Fpdensl1;   //!
   TBranch        *b_Fpdetbl1;   //!
   TBranch        *b_Fpdwnsl1;   //!
   TBranch        *b_Fpdwtbl1;   //!
   TBranch        *b_bunchid;   //!
   TBranch        *b_bunchid7bit;   //!

   h111(TTree *tree=0);
   virtual ~h111();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef h111_cxx
h111::h111(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../tmp/run6173032.1.root");
      if (!f) {
         f = new TFile("../tmp/run6173032.1.root");
      }
      tree = (TTree*)gDirectory->Get("h111");

   }
   Init(tree);
}

h111::~h111()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t h111::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t h111::LoadTree(Long64_t entry)
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

void h111::Init(TTree *tree)
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
   fChain->SetBranchAddress("ipre", &ipre, &b_ipre);
   fChain->SetBranchAddress("ipost", &ipost, &b_ipost);
   fChain->SetBranchAddress("Bchi", &Bchi, &b_Bchi);
   fChain->SetBranchAddress("Bclo", &Bclo, &b_Bclo);
   fChain->SetBranchAddress("Awdbm", &Awdbm, &b_Awdbm);
   fChain->SetBranchAddress("Awtrg", &Awtrg, &b_Awtrg);
   fChain->SetBranchAddress("Awdaq", &Awdaq, &b_Awdaq);
   fChain->SetBranchAddress("Token", &Token, &b_Token);
   fChain->SetBranchAddress("addbit", &addbit, &b_addbit);
   fChain->SetBranchAddress("Dsm", &Dsm, &b_Dsm);
   fChain->SetBranchAddress("extbusy", &extbusy, &b_extbusy);
   fChain->SetBranchAddress("modbusy", &modbusy, &b_modbusy);
   fChain->SetBranchAddress("physwd", &physwd, &b_physwd);
   fChain->SetBranchAddress("Trgwd", &Trgwd, &b_Trgwd);
   fChain->SetBranchAddress("Dsmad", &Dsmad, &b_Dsmad);
   fChain->SetBranchAddress("cotamibusy", &cotamibusy, &b_cotamibusy);
   fChain->SetBranchAddress("npre", &npre, &b_npre);
   fChain->SetBranchAddress("npost", &npost, &b_npost);
   fChain->SetBranchAddress("Cpa", Cpa, &b_Cpa);
   fChain->SetBranchAddress("qdsm", qdsm, &b_qdsm);
   fChain->SetBranchAddress("lastdsm", lastdsm, &b_lastdsm);
   fChain->SetBranchAddress("Bcd", Bcd, &b_Bcd);
   fChain->SetBranchAddress("L1", L1, &b_L1);
   fChain->SetBranchAddress("Vtxdsm", Vtxdsm, &b_Vtxdsm);
   fChain->SetBranchAddress("Fpddsm", Fpddsm, &b_Fpddsm);
   fChain->SetBranchAddress("Emcdsm", Emcdsm, &b_Emcdsm);
   fChain->SetBranchAddress("Bbcl1", Bbcl1, &b_Bbcl1);
   fChain->SetBranchAddress("Zdc", Zdc, &b_Zdc);
   fChain->SetBranchAddress("Zdcsmd", Zdcsmd, &b_Zdcsmd);
   fChain->SetBranchAddress("Ctb", Ctb, &b_Ctb);
   fChain->SetBranchAddress("Mwc", Mwc, &b_Mwc);
   fChain->SetBranchAddress("Bemc", Bemc, &b_Bemc);
   fChain->SetBranchAddress("Bbc", Bbc, &b_Bbc);
   fChain->SetBranchAddress("Fpdens", Fpdens, &b_Fpdens);
   fChain->SetBranchAddress("Fpdetb", Fpdetb, &b_Fpdetb);
   fChain->SetBranchAddress("Fpdwns", Fpdwns, &b_Fpdwns);
   fChain->SetBranchAddress("Fpdwtb", Fpdwtb, &b_Fpdwtb);
   fChain->SetBranchAddress("Fpdensl1", Fpdensl1, &b_Fpdensl1);
   fChain->SetBranchAddress("Fpdetbl1", Fpdetbl1, &b_Fpdetbl1);
   fChain->SetBranchAddress("Fpdwnsl1", Fpdwnsl1, &b_Fpdwnsl1);
   fChain->SetBranchAddress("Fpdwtbl1", Fpdwtbl1, &b_Fpdwtbl1);
   fChain->SetBranchAddress("bunchid", &bunchid, &b_bunchid);
   fChain->SetBranchAddress("bunchid7bit", &bunchid7bit, &b_bunchid7bit);
   Notify();
}

Bool_t h111::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void h111::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t h111::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef h111_cxx
