#ifndef RfEntry
#define RfEntry
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TObject.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TVector.h"
#include <iostream>
#include  "fOneD.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayI.h"
#include "TString.h"
#include "Rtypes.h"

/*
extern Int_t fpdenmap[7][7];
extern Int_t fpdesmap[7][7];
extern Int_t fpdetmap[5][5];
extern Int_t fpdebmap[5][5];
extern Int_t fpdwsmap[7][7];
extern Int_t fpdwbmap[5][5];
*/

TMatrix MakeMatI(Int_t ,Int_t* );

class fEntry :public TObject
{
public:
  fEntry(char*,int, int ,int,int,TMatrix&);
  fEntry(){};
  ~fEntry();
  Int_t Row;
  Int_t Col;
  Int_t Chan;
  Int_t Nbins;
  TString DetName;
  TString HistName;
  TString HistTitle;
  TString PHistName;
  TH1F* Hist;
  Float_t Ped;
  void SetPed(ostream &);
  TH1F* HistPedSubtract;
  Float_t Gain();
  /*
  TMatrix FpdWSmap(){return MakeMatI(7,fpdwsmap[0]);};
  TMatrix FpdESmap(){return MakeMatI(7,fpdesmap[0]);};
  TMatrix FpdENmap(){return MakeMatI(7,fpdenmap[0]);};
  TMatrix FpdETmap(){return MakeMatI(5,fpdetmap[0]);};
  TMatrix FpdEBmap(){return MakeMatI(5,fpdebmap[0]);};
  TMatrix FpdWBmap(){return MakeMatI(5,fpdwbmap[0]);};
  */
  ClassDef(fEntry,1);
};

#endif
