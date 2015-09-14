#ifndef Asymmetry_
#define Asymmetry_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TGraphErrors.h>


class Asymmetry : public TObject
{
 public:
  Asymmetry(TH2F* h,TString name="h2Asm",Float_t norm=1.,TString titleString="");
  TH2F* h2dAsm; 
  Float_t Norm;
  TString TitleString;
  Int_t NbinX(){return h2dAsm->GetNbinsX();};
  Float_t LowX(){return h2dAsm->GetXaxis()->GetXmin();};
  Float_t HiX(){return h2dAsm->GetXaxis()->GetXmax();};
  void ReNorm(Float_t factor);
  TString SpinNames[8];
  TH1D* CrossAN(TString beam="yellow");
  Int_t SpinLoc[8];
  TH1D* SpinH[8];
  TString NSsumNames[4];
  TH1D* ANS[4];
  TH1D* NSSpinH[4];//sum over N & S
  TH1D* BN[2];
  TH1D* BS[2];
  TH1D* YN[2];
  TH1D* YS[2];
  TH1D* BNorth;
  TH1D* BSouth;
  TH1D* YNorth;
  TH1D* YSouth;
  TH1D* North;
  TH1D* South;
  TH1D* All;
  TH1D* AYN;
  TH1D* AYS;
  TH1D* ABN;
  TH1D* ABS;
  TH1D* AS_xx_yx;
  TH1D* AN_xx_yx;
  TH1D* AS_xx_yxp;
  TH1D* AN_xx_yxp;
  TH1D* AN_pp_mm;
  TH1D* AS_pp_mm;
  TH1D* CrossY;
  TH1D* CrossB;
  void PlotYBNS(TCanvas* cc=0,TString XAxisName="",Float_t xmin=0,Float_t xmax=250);
  TGraphErrors  ValueDist(TH1D h1,TString name="Values");
  TH1D ErrorDist(TH1D h1,TString name="eValues");
  ~Asymmetry();
 private:
    ClassDef(Asymmetry,1); 
};

#endif
