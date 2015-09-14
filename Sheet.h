#ifndef Sheet_
#define Sheet_
#include "TObject.h"
#include "TObjArray.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TKey.h"
#include "TF1.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "TString.h"
#include "TObjString.h"
#include "TAxis.h"
#include "TDatime.h"

class SheetHistory : public TObject
{ 
 public:
  SheetHistory(Float_t oldhv,Float_t newhv,Float_t oldgcorr,
	       Float_t newgcorr,Int_t oldshift,Int_t newshift,
	       Float_t gain=.038);
  SheetHistory(){};
  ~SheetHistory(){};
  Float_t OldHV;
  Float_t NewHV;
  Float_t OldGcorr;
  Float_t NewGcorr;
  Int_t OldShift;
  Int_t NewShift;
  Float_t Gain;
  TString CreateDate;
  ClassDef(SheetHistory,1);   
};

class Sheet : public TObject
{
 public:
  Sheet();
  ~Sheet(){};
  Sheet(char* cvsHead,char* cvsLine);
  const char*  GetName()const {return (const char*) Name;};
  void SetName(char* nm){Name=nm;};
  void Print();
  void SetShift(Int_t shft){Parsecol=shft;};
  Int_t GetShift(){return Parsecol;};
  TString DoString(char* fname);
  void Parse(TString Str,Int_t Iport,Int_t NSTB1);
  Int_t Parsecol;
  Int_t NumberElem;
  Int_t nstb1;
  Int_t chan1;
  Int_t myrow0;
  Int_t mycol0;
  Int_t PatchPanel;
  Int_t PatchRow;
  Int_t PatchCol;
  Int_t QTCrate;
  Int_t QTSlot;
  Int_t QTCard4;
  Int_t QTchannel8;
  Int_t id1;
  Int_t id2;
  Int_t id3;
  Int_t id4;
  Int_t HVsetValue;
  Float_t Run9Gain;
  TString Vname[20];
  TString TheLine;
  TString HeadLine;
  TString SetCommand(Int_t bitmask,Bool_t nw=true);
  Float_t HVfactor;
  TGraphErrors* GetGr(){return gr;};
  TGraphErrors* gr;
  TF1* GetGrFun();
  TF1* fun;
  void AddGr(TFile* HvGr);
  TObjArray* FillResets(char* NewSetFile);
  TObjArray* Resets;  
  TObjArray* Extras;
  Int_t DoReset(TObjArray* resets);
  Int_t HvScale(Float_t gcorr,Int_t shift=0);
  void AddHistory(Float_t oldhv,Float_t newhv,Float_t oldgcorr,
		  Float_t newgcorr,Int_t oldshift,Int_t newshift,
		  Float_t gain=.038);
  SheetHistory* GetHistory(Int_t depth=0);
  TString StringHistory(Int_t depth=0);
  TString QTline(Bool_t nw=true);
  Int_t Getshift(Bool_t nw=true,Int_t depth=0);  
  Int_t GetHV(Bool_t nw=true,Int_t depth=0);  
  Float_t GetGcorr(Bool_t nw=true,Int_t depth=0);
  void BackupHistory();
  TString MakecsvLine();
  TString MakecsvHead();
  static Sheet* FindSheetid(TObjArray* ar,Int_t id1,Int_t id2,Int_t id3, Int_t id4);
  static Sheet* FindSheetqt(TObjArray* ar,Int_t qtcrate,Int_t qtslot, Int_t qtchan);

 private:
  TString Name;
  ClassDef(Sheet,1);
};

#endif
