#ifndef TrigQt_
#define TrigQt_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "Qt.h"
#include "ClusterQt.h"
#include "TObjArray.h"
#include "TString.h"
#include <iostream>

class TrigQt : public TObject
{
 public:
  TrigQt(Qt* p_Qt);
  TrigQt(Qt* p_Qt,Bool_t qt8trunc);//truncate qt8 sums before forming clusters to emulate L0dsm
  void Init();
  ~TrigQt();
  Qt* pQt;
  TString Loc(Int_t chan, Int_t slot, Int_t crate);
  TObjArray* TrigClusters;
  Int_t L0_clsum[4][3]; //[L0dsm index in L1dsm][L1dsm]
  Int_t L0_htadr[4][3]; //[L0dsm index in L1dsm][L1dsm] 
  Int_t L0_ht[4][3]; //[L0dsm index in L1dsm][L1dsm]
  Int_t L0_ht_thval[4][3]; //[L0dsm index in L1dsm][L1dsm]
  Bool_t L0_ht_thbit[4][3]; //[L0dsm index in L1dsm][L1dsm]

  UInt_t L1_TrigBits[3]; //[L1 dsms]

  UInt_t L2_TrigBits;

  Int_t thS[3]; //[L1 dsms]
  Int_t thL[3]; //[L1 dsms]
  UInt_t ntrig[2][3]; //[SL cells][3 threshold]

  UInt_t clbitsSmall;
  UInt_t mclbitsSmall;
  UInt_t htbitsSmall;
  UInt_t clbitsLarge;
  UInt_t mclbitsLarge;
  UInt_t htbitsLarge;

  void AddRealEdgeClusters(); //add qt8sum to edge high towers
  void AddEdgeClusters(); //add qt8sum to edge high towers
  void L0dsm(); //combine qt slot to form L0 dsm data
  void L12dsm(); //apply L1 and L2 trigger algo.
  void printDSM();
  
 private:
  ClassDef(TrigQt,1);
};
#endif
