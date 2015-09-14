#ifndef ClusterQt_
#define ClusterQt_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "Qt.h"
#include "TObjArray.h"
#include "TString.h"
#include <iostream>


class ClusterQt: public TObject
{
 public:
  ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_);//everything counts from 0 here
  ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Bool_t qt8trunc_);//truncate qt8 sums before forming clusters to emulate L0dsm
  ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Int_t ch_adr );
  ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Int_t ch_adr,Bool_t qt8trunc_);//truncate qt8 sums before forming clusters to emulate L0dsm
  ~ClusterQt();

  Qt* pQt;
  Int_t sum;
  Int_t slot;
  Int_t crate;
  Int_t qthtadr;//0~31
  Int_t qtht;
  Int_t detnstb;
  Int_t detchan;
  Int_t Nqt8;
  Int_t qt8sum[4];
  Int_t Cl_crate[4];
  Int_t Cl_slot[4];
  Int_t Cl_qt8[4];
  Bool_t Defined;
  Bool_t L0select;
  Bool_t qt8trunc;
   
  void GetClusterSum();
  Bool_t AddQt8(Int_t crate_,Int_t slot_,Int_t qt8_);
  void Print();
  void Init();

 private:
  ClassDef(ClusterQt,1);
};
#endif
