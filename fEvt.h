#ifndef RfEvt
#define RfEvt
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
#include "TMinuit.h"
#include "Fit49.h"
#include "photon.h"

typedef TMatrixTRow<Float_t> TMatrixRow; 

class fEvt  :public TObject
{
public:
  fEvt(Int_t,TH2F*,Int_t,Int_t,ostream&);
  ~fEvt();
  fEvt(){Ehist=0;};
  TH1F HistVec(char*,TVector&,Float_t);
  void MakeMat();
  TVector SortBins(TMatrix*);
  TH2F* Ehist;
  TH2F* Fithist;
  TH2F* ChiFit;
  TVector* Energies;
  TVector* V1x;
  TVector* V1y;
  TVector* V2x;
  TVector* V2y;
  Float_t PosOnEigenvector(TMatrixRow);
  ostream* Os;
  Int_t Event;
  Int_t Run;
  Int_t SubRun;
  Int_t MaxRow;
  Int_t MaxCol;
  Int_t NumberGammas;
  Double_t Split1pos;
  Double_t Split2pos;
  Double_t SplitE1;
  Double_t SplitE2;
  Int_t Noshape(TVector*);
  Double_t Sum1;
  Double_t Sum2;
  Double_t E1;
  Double_t E2;
  Double_t ZSum;
  Double_t Epair;
  Double_t X1;
  Double_t Y1;
  Double_t X2;
  Double_t Y2;
  Int_t NGoodTracks;
  Int_t NextRow;
  Int_t NextCol;
  Double_t CountSum;
  Double_t BinContents[7][7];
  Double_t Angle;
  Double_t Mass;
  TMatrix* mat;
  TMatrix* E1FUN;
  Float_t ADCsum;
  Double_t Dot(TVector,TVector);
  Double_t Dot1;
  Double_t Dot2;
  Double_t X1F;
  Double_t X2F;
  Double_t Y1F;
  Double_t Y2F;
  Double_t MassF;
  Double_t E1F;
  Double_t E2F;
  Double_t EpairF;
  Double_t AngleF;  
  Double_t ZSumF;
  Double_t ChisquareF;
  TVector* peval;
  void PrintInfo();
  Int_t DoFit();
  Int_t DoFit1();
  Int_t ReAnalyse();
  Int_t HiTest();
private:
  TVector* AllOnes;
  ClassDef(fEvt,1);
};
#endif
