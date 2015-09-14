#ifndef Gen_
#define Gen_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "FillData.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TString.h"
#include "TComplex.h"
#include "FilesSet.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Geom.h"
#include "TH2F.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "dataSet.h"
#include "FpdMap.h"
#include "CalibStr.h"
#include "FitTower.h"
#include "Yiqun.h"
#include "TMath.h"
#include "Trigger.h"
#include "Sim.h"
#include "TKey.h"
#include "TIterator.h"


class Gen : public TObject
{
 public:
  Gen(Int_t IEW, Int_t NSTB, Int_t nev,Geom* p_Geom,CalibStr* RGain,CalibStr* RGcor,Int_t seed);
  ~Gen(){};
  TF2* p_fcnSS;
  TLorentzVector* p_pi;
  Float_t theta_s;
  Int_t Error;
  Double_t (* pGGams)(Double_t *x, Double_t *par);
  static  Double_t EPdist(Double_t *x, Double_t *par);
  static  Double_t Ptdist(Double_t *x, Double_t *par);
  static  Double_t Xdist(Double_t *x, Double_t *par);
  Float_t Fill();
  Float_t TowerWidth;
  FitTower* fitter;
  Int_t Nevents;
  Double_t NominalEta;
  Double_t NominalPhi;
  Int_t Iew;
  Int_t Instb;
  Geom* pGeom;
  Int_t Nrows;
  Int_t Ncols;
  Double_t xx,yy,eSS;
  CalibStr* Rgain;
  CalibStr* Rgaincorr; 
  TMatrix* NextEvt;
  TMatrix* Digitized;
  TMatrix* UnDigitized;
  TMatrix* DigitizeEnergy();
  Double_t BeamZ;
  Float_t PreDigitizedEnergy;
  Yiqun* rec;
  TF1* PtGen;
  TF1* XGen;
  TF2* EPGen;
  TH1F* Mass;
  TTree* g_out;
  TH2F* MassvsEnergy;
  TLorentzVector* p_gam1;
  TLorentzVector* p_gam2;
  Bool_t DoRecon;
  TH1F* Nphotons;
  TH2F* Eratch1;
  TH2F* Eratch1F;
  TH2F* EratchP;
  TH2F* Eratch2;
  TH2F* Eratch40;
  TH1F* phi1;
  TH1F* phi2;
  TH1F* phi3;
  Double_t Energy;
  Double_t Y;
  Double_t pt;
  Double_t xF;
  Float_t phi;
  Float_t phi_s;
  Float_t Z;
  TVector3 LocalGamPos[2];
  TRandom* rand;
  TH1F* hEgen;
  TH1F* hEpass;
  Bool_t piInside;
  Bool_t inFiducial;
  void CallFill(Int_t ncalls);
  TH1F* Exy[7][7];
  TH1F* Mcell[7][7];
  void WriteAll();
  TH2F* Episum;
  TH2F* Epipsum;
  TH2F* EpipbRsum;
  TH2F* EpipRsum;
  TH2F* Epigsum;
  TH1F* GenZbeam;
  TH1F* GenPZbeam;
  
  TH2F* GenEvsPt;
  TH2F* GenPEvsPt;
  TH2F* GenEvsY;
  TH2F* GenpEvsY;
  TH2F* GenThvsPh;
  TH2F* GenPThvsPh;
  TH2F* GenYvsPtWPt;
  TH2F* GenPYvsPtWPt;
  TH2F* ZvsADC[7][7];
  TH2F* PhivsADC[7][7];
  TH2F* LocalPos[7][7];
  TH1F* EphXY[7][7];
  TH2F* LocalPos125[7][7];
  TH2F* NphotvsE[7][7];
  TH2F* GenPPEvsPt;
  TH2F* GenPPThvsPh;
  TH2F* GenPPYvsPtWPt;
  TH2F* Z_Phi_40[10][10];
  void EnergyFromPhoton(TVector3* p_v3,Double_t photE ); 
  Float_t Zlimit;
  Float_t CallFillLimitE;
  Float_t EnergySum;
  Int_t GenCnt;
  Int_t NPh;
  Float_t Zrec;
  Float_t Eph1;
  Float_t Eph2;
  Float_t pimass;
  Float_t mass;
  Float_t piPt;
  Float_t piEta;
  Float_t piPhi;
 private:
  ClassDef(Gen,1);
};
#endif
