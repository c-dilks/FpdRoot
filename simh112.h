//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Dec 28 12:03:35 2007 by ROOT version 5.12/00h-rc2
// from TTree h112/SIMNTP
// found on file: fpd++sim_pythia10014.root
//////////////////////////////////////////////////////////

#ifndef simh112_h
#define simh112_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class simh112 {
public :
   TChain *Input;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   UInt_t          event;
   Float_t         sume[2][6];
   Float_t         sump[2][6];
   Float_t         sums[2][2][2];
   Int_t           Nhit;
   UInt_t          id[10000];   //[Nhit]
   Float_t         hit[10000];   //[Nhit]
   UInt_t          trk[10000];   //[Nhit]
   Int_t           Ntrack;
   UInt_t          ptrk[500];   //[Ntrack]
   UInt_t          stat[500];   //[Ntrack]
   UInt_t          pdg[500];   //[Ntrack]
   UInt_t          mo[500][2];   //[Ntrack]
   UInt_t          da[500][2];   //[Ntrack]
   Float_t         p[500][3];   //[Ntrack]
   Float_t         e[500];   //[Ntrack]
   Float_t         m[500];   //[Ntrack]
   Float_t         v[500][4];   //[Ntrack]
   UInt_t          gpid[500];   //[Ntrack]
   UInt_t          gtrk[500];   //[Ntrack]
   Float_t         pt[500];   //[Ntrack]
   Float_t         eta[500];   //[Ntrack]
   Float_t         phi[500];   //[Ntrack]
   UInt_t          proc;
   Float_t         shat;
   Float_t         that;
   Float_t         uhat;
   Float_t         cost;
   Float_t         q2;
   UInt_t          in[2];
   Float_t         xin[2];
   Float_t         x[2];

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_sume;   //!
   TBranch        *b_sump;   //!
   TBranch        *b_sums;   //!
   TBranch        *b_Nhit;   //!
   TBranch        *b_id;   //!
   TBranch        *b_hit;   //!
   TBranch        *b_trk;   //!
   TBranch        *b_Ntrack;   //!
   TBranch        *b_ptrk;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_pdg;   //!
   TBranch        *b_mo;   //!
   TBranch        *b_da;   //!
   TBranch        *b_p;   //!
   TBranch        *b_e;   //!
   TBranch        *b_m;   //!
   TBranch        *b_v;   //!
   TBranch        *b_gpid;   //!
   TBranch        *b_gtrk;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_proc;   //!
   TBranch        *b_shat;   //!
   TBranch        *b_that;   //!
   TBranch        *b_uhat;   //!
   TBranch        *b_cost;   //!
   TBranch        *b_q2;   //!
   TBranch        *b_in;   //!
   TBranch        *b_xin;   //!
   TBranch        *b_x;   //!

   simh112(char* files);
   virtual ~simh112();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
 private:
    ClassDef(simh112,1);

};


#endif
