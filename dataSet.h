#ifndef dataSet_h
#define dataSet_h
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TObject.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TH2F.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include "Fill.h"
#include "RunData.h"
#include "FilesSet.h"
#include "Geom.h"
#include "Qt.h"

class dataSet {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TString NtpName;
   dataSet(char*,FilesSet*,char*);
   dataSet(char*,char*,char*,char*);
   dataSet(char*,char*,char*,char*,char*);
   dataSet();
   Int_t dataSetBuild(char*,char*,char*,char*);
   TString LastError;
   FilesSet* pFiles;
   Geom* pGeom;
   Bool_t UseQt;
   Qt* pQt;
   Bool_t decodeQT();
   Bool_t decodeDSM();
   // Declaration of leave types
   Int_t           event;
   UChar_t         ipre;
   UChar_t         ipost;
   UInt_t          Bchi;
   UInt_t          Bclo;
   UShort_t        Token;
   UInt_t          addbit;
   UInt_t          Dsm;
   UChar_t         npre;
   UChar_t         npost;
   UInt_t          lastdsm[8];
   UInt_t          Vtxdsm[8];
   UInt_t          Emcdsm[8];
   UInt_t          L1[32];
   UChar_t         Ctb[256];
   UChar_t         Fpdens[112];
   UInt_t          Fpdensl1[8];
   UInt_t          Bbcl1[16];
   UChar_t         Zdc[16];
   UChar_t         Zdcsmd[32];

   UInt_t          Awdbm;
   UChar_t         Awtrg;
   UChar_t         Awdaq;
   UInt_t          extbusy;
   UInt_t          modbusy;
   UInt_t          physwd;
   UInt_t          Trgwd[3];
   UInt_t          Dsmad;
   UInt_t          cotamibusy;
   UInt_t          Cpa[32];
   UInt_t          qdsm[8];
   UInt_t          Bcd[16];
   UInt_t          bunchid;
   UInt_t          bunchid7bit;

   UInt_t          Fpddsm[8];
   UChar_t         Mwc[128];
   UChar_t         Bemc[128];
   UChar_t         Bbc[96];
   UChar_t         Fpdetb[64];
   UChar_t         Fpdwns[112];
   UChar_t         Fpdwtb[64];
   UInt_t          Fpdetbl1[8];
   UInt_t          Fpdwnsl1[8];
   UInt_t          Fpdwtbl1[8];
   
   //h112 ntuples

   Int_t           num10;
   Int_t           num50;
   UInt_t          rega;
   UInt_t          cwnentry;
   UInt_t          fpdwadc[352];
   UInt_t          tdc[8];

   //h111_08 ntuples
   UInt_t          Awd[3];
   UInt_t          busystatus[3];
   UInt_t          Mtddsm[8];
   UInt_t          Vpddsm[8];
   UInt_t          Ctbdsm[8];
   UInt_t          specialtriggers[8];
   UInt_t          L2[64];


   UChar_t         Mtd[32];
   UChar_t         Vpd[64];
   UChar_t         Tof[16];
   UChar_t         Bemceast[240];
   UChar_t         Bemcwest[240];
   UInt_t          Bemclayer1[48];
   UInt_t          Eemclayer1[16];
   UChar_t         Eemc[144];
   UChar_t         Fpdw[256];
   UInt_t          Zdcl1[8];
   Int_t           nqtdata;
   UInt_t          Qtdata[3700];   //[nqtdata]
   Int_t           nemcdata;
   UInt_t          Emcdata[100];   //[nemcdata]

   //h111_09 ntuples
   Int_t           evt;
   UChar_t         iprepost;
   UInt_t          Bc[2];
   UChar_t         nprepost[2];
   UInt_t          Tofdsm[8];
   UInt_t          Vpd_09[8];
   UInt_t          Fpde[8];
   UChar_t         Mtd_09[16];
   UInt_t          Tofl1[8];
   UInt_t          Tof_09[48];
   UInt_t          L2sum[2];
 

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

   //h111.h

   TBranch        *b_num10;   //!
   TBranch        *b_num50;   //!
   TBranch        *b_rega;   //!
   TBranch        *b_cwnentry;   //!
   TBranch        *b_fpdwadc;   //!
   TBranch        *b_tdc;   //!

   //h111_08

   TBranch        *b_Awd;   //!
   TBranch        *b_busystatus;   //!
   TBranch        *b_Mtddsm;   //!
   TBranch        *b_Vpddsm;   //!
   TBranch        *b_Ctbdsm;   //!
   TBranch        *b_specialtriggers;   //!
   TBranch        *b_L2;   //!

   TBranch        *b_Mtd;   //!
   TBranch        *b_Vpd;   //!
   TBranch        *b_Tof;   //!
   TBranch        *b_Bemceast;   //!
   TBranch        *b_Bemcwest;   //!
   TBranch        *b_Bemclayer1;   //!
   TBranch        *b_Eemclayer1;   //!
   TBranch        *b_Eemc;   //!
   TBranch        *b_Fpdw;   //!
   TBranch        *b_Zdcl1;   //!
   TBranch        *b_nqtdata;   //!
   TBranch        *b_Qtdata;   //!
   TBranch        *b_nemcdata;   //!
   TBranch        *b_Emcdata;   //!

  //h111_09

   TBranch        *b_evt;   //!
   TBranch        *b_iprepost;   //!
   TBranch        *b_Bc;   //!
   TBranch        *b_nprepost;   //!
   TBranch        *b_Tofdsm;   //!
   TBranch        *b_Tofl1;   //!
   TBranch        *b_Fpde;   //!
   TBranch        *b_L2sum;   //!

  ~dataSet();
  TMatrix dMatrix(TMatrix* ,UChar_t* );
  TMatrix dMatrix(TMatrix* ,UInt_t* );
  TMatrix OutProd(TVector&, TVector&);
  TMatrix dMatrixQt(Int_t);
  TChain* Input;
  Int_t error;
  TMatrix ZeroNegativePed2(TMatrix*,TMatrix*,TMatrix*);
  TMatrix ZeroNegativePedGain2(TMatrix*,TMatrix*,TMatrix*,TMatrix*);
  TMatrix ZeroNegativePed2Gain2(TMatrix*,TMatrix*,TMatrix*,TMatrix*,TMatrix*);
  CalibStr* Rped;  
  CalibStr* Rgain;  
  CalibStr* Rgaincorr;
  CalibStr* RunPed;
  Fill* RFill;
  TObjArray* RunsArray;
  RunData* thisRunDat;
  Int_t    Cut(Long64_t entry);
  Int_t    GetEntry(Long64_t entry);
  Long64_t LoadTree(Long64_t entry);
  void     Init(TTree *tree);
  void     Loop();
  Bool_t   Notify();
  void     Show(Long64_t entry = -1);
  TString FileName;
  Int_t    CurrentRunNumber;
  Int_t    CurrentSegNumber;
  Int_t    CurrentEventNumber;
  Float_t R1Lum;
  Float_t R2Lum;
  Float_t ErrR1Lum;
  Float_t ErrR2Lum;
  Int_t   BlueSpin;
  Int_t   YellowSpin;
  Bool_t  kicked;
  TMatrixF Em(TMatrix* t0,Int_t, Int_t);
 private:
  Int_t currentTreeNumber;
  ClassDef(dataSet,1);
};

#endif

