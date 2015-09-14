#ifndef Cell_
#define Cell_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "FillData.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TStyle.h"
#include "TMatrix.h"
#include "TString.h"
#include "TComplex.h"
#include "FilesSet.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Geom.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "dataSet.h"
#include "FpdMap.h"
#include "CalibStr.h"
#include "FitTower.h"
#include "TMath.h"
#include "Trigger.h"
#include "Sim.h"
#include "TKey.h"
#include "TIterator.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "Yiqun.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "RunDepMgr.h"

class tpCellDat : public TObject 
{
 public:
  tpCellDat(){};
  Float_t Mcell;
  Float_t Epair;
  Float_t Ypair;
  Float_t E1;  Float_t E2;
  Float_t X1;
  Float_t X2;
  Float_t Y1;
  Float_t Y2;
  unsigned SavedHits[100];
  Int_t nSavedHits;
  Int_t NSTB;
  Int_t ievt;
  Int_t Rnum;
 private:
  ClassDef(tpCellDat,2)
};

class tpCellRec : public TObject 
{
 public:
  tpCellRec(){};
  Float_t mr;
  Float_t er;
  Float_t e1;
  Float_t x1;
  Float_t x2;
  Float_t y1;
  Float_t y2;
  Int_t n1;
  Int_t n2;
  Int_t singles;
  Int_t NRclust;
  Float_t Sel;
  Float_t clE1;
  Float_t clE2;
  Float_t sigmax;
  Float_t sigmin;
  Float_t chi1;
  Int_t Rank1;
  Int_t Rank2;
 private:
  ClassDef(tpCellRec,2)
};

class CellHit : public TObject 
{
 public:
  CellHit(){};
  CellHit(unsigned int );
  unsigned int packed;
  Int_t EW;
  Int_t NSTB;
  Int_t Row;
  Int_t Col;
  Int_t ADC;
 private:
  ClassDef(CellHit,1)
};

class Cell : public TObject
{
 public:
  Cell();
  Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,CalibStr* rgcor);
  Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,CalibStr* rgain,CalibStr* rgcor); 
  Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,TFile* fl,CalibStr* rgcor); 
  Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,TFile* fl,CalibStr* rgcor,CalibStr* NominalMass);
  Bool_t Setup(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum); 
  ~Cell();
  tpCellDat celldat;
  TTree* CellD;
  void InitTree(const char*);
  void InitBranches();
  void CreateBranches();
  void FillTree(tpCellDat);
  
  Int_t Iew,Instb,Row1,Col1,Rnum;
  Float_t xLocal,yLocal;
  Float_t xGlobal,yGlobal;
  Bool_t FindHists(TFile* tfile);
  CalibStr* p_NomMass;//!
  TH1F* p_Mcell;
  TH2F* p_MvsZ;
  TH1F* p_EphXY;
  TH1F* p_RawAdc;
  TH1F* p_RawAdcL;
  TH1D* p_adc;
  TH1D* p_adcLed;
  TH2F* p_MvsD;
  TH2F* p_MvsE;
  TH2F* p_MvsY;
  TH3F* p_dxdyM;
  TH3F* p_dxdyE;
  TH3F* p_dxdyER;
  
  CalibStr* fpdcorr;
  CalibStr* fpdgain;
  Float_t gcorr;
  Float_t gain;
  Float_t Mass0;
  Float_t FitMcell(Float_t m1=.08,Float_t m2=.18);
  Float_t FitMcell2(Float_t m1=.08,Float_t m2=.18);
  Float_t FitMass;
  Float_t FitSlope;
  Float_t erFitMass;
  Float_t FitMassWidth;
  Float_t Massmaxcontents;
  Float_t Masspeakfraction;
  Float_t MassFitError;
  Float_t MassEntries;
  Float_t FitRawCell();
  Float_t FitRawLCell(Float_t a1=10,Float_t a2=2);
  Float_t hiRawSlope;
  Float_t hiRawIntdorig;
  Float_t hiRawSlopeL;
  Float_t ReFitQuality0to3;
  Float_t erRawSlope;
  Float_t Beyond120;
  Float_t Beyond145;
  Float_t BeyondL10;
  Float_t raw120;
  Float_t hiRawInt;
  Float_t erhiRawInt;
  Float_t erRawSlopeL;
  Float_t IntRawFitwithgcorr;
  TVector LocalGamPos[2];
  Float_t Ptp01pct;
  Float_t Energyp01pct;
  Float_t Energyp1pct;
  Float_t Eta;
  Float_t Phi;
  FilesSet* pfiles;
  Geom* p_Geom;
  Bool_t Geom_ok;
  Bool_t Setpfiles(FilesSet* _pfiles=0);
  Float_t Newcorr();
  Float_t gcorrNew;
  Float_t gcorrOrig;
  TH1F* p_Nphthis;
  TH1F* p_NphAll;
  TH1F* p_DetE;
  TH1F* p_FmsE;
  char ReasonForChange[200];
  Int_t Compare(const TObject*) const;
  void SetName(char* nm){fName=nm;};
  Bool_t AddToReason(char* reason);
  const char* GetName() const {return (const char*) fName;};
  TH1D* ReScaleHist(TH1D* p_hist,Float_t xfactor,Int_t minbin,Bool_t Normalize);  Float_t TriggerPeak();//return trigger threshold Energy based on p_DetE
  TH1D* ReScaledHist;
  Float_t LastE_Det;
  Float_t LastE_Det_Error;
  Float_t NDF_Det;
  Float_t Chi2_Det;
  Float_t Slope_Det;
  Float_t Det_Peak_Loc;
  TGraphErrors* FitMvsE(TH2F* hme,Int_t ndiv=5, char* opt="LQ0",Int_t nstart=1,Int_t nend=10000,TH1D* slicearray=0);
  //  TGraphErrors* p_gFME;
  void DrawFMSADC(Int_t NSTB,Bool_t UseEnergy=false);
  TMatrixT<float> FillFMSADC(Int_t NSTB);
  TMatrixT<float> FillFMSEnergy(Int_t NSTB);
  TObjArray* CellHitList();
  TH2F* Getgradc(Int_t NSTB){return gradc[NSTB-1];};
  TTree* ReRecon(tpCellRec* cellrec,CalibStr* Fpdcorr=0,Int_t samplemax=-1,CalibStr* Rank=0,TTreeFormula* precut=0,Bool_t UseQcut=false);
  Int_t ClNearX(Yiqun* prec,Float_t x1, Float_t y1);
  Int_t PhNearX(Yiqun* prec,Float_t x1,Float_t y1);
  Int_t RankCell(CalibStr* rank);
  Int_t Quality(CalibStr* rank);
  void SetEmin(Float_t emin){Emin=emin;};
  RunDepCor* Rdep;
  Int_t RdepBaseOverride;
  private:
  TGraph* gr[4];
  TH2F* gradc[4];
  TString fName;
  Bool_t CheckFiducial(int iew0,int instb0,TVector3* vlocal, CalibStr* rank,Float_t cellfract=.5);
  Float_t Emin;
  RunDepMgr* Mgr;
  ClassDef(Cell,10);
};

#endif
