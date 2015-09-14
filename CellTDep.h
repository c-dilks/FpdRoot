#ifndef CellTDep_
#define CellTDep_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "TObject.h"
#include "CalibStr.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TFile.h"

class CellTDep : public TObject
{
 public:
  CellTDep()
    {
      gr=new TGraphErrors();
      fit=new TF1();
    };
  CellTDep(TGraphErrors* Gr,Int_t NSTB,Int_t Row0, Int_t Col0, CalibStr* GCor,Long_t Rnum,Float_t setNominalADC=0);
  ~CellTDep();
  TF1* GetFit(){return fit;};
  Int_t GetNstb(){return nstb;};
  Int_t GetRow0(){return row0;};
  Int_t GetCol0(){return col0;};
  Int_t GetChan(){return chan;};
  TGraphErrors* GetGraph(){return gr;};
  Int_t ErrorStatus;//zero for a good fit
  Float_t GetNominalLedADC(){return NominalLedADC;};
  //  Float_t SetNominalLedADC(Int_t nom){NominalLedADC=nom;};
  Float_t SetNominalLedADC(Float_t nom){NominalLedADC=nom;};
  Float_t GetNominalGcor(){return NominalGcor;};
  Float_t GetthisLedADC(){return thisLedADC;};

  Float_t FixGcor(CalibStr* FpdCor,Int_t Segnum,Int_t Evtnum,Int_t mode); 
  Float_t FixGcor(CalibStr* FpdCor,Int_t LongEvtNum,Int_t mode)
  {return FixGcor(FpdCor,LongEvtNum/10000,LongEvtNum%10000, mode);};
  
  Float_t FixFactor(Int_t Segnum,Int_t Evtnum,Int_t mode)
  {CalibStr* FpdCor=0;return  FixGcor(FpdCor,Segnum,Evtnum,mode); };
  Float_t FixFactor(Int_t LongEvtNum,Int_t mode)
  {return FixFactor(LongEvtNum/10000,LongEvtNum%10000, mode);};
  
  Float_t GetLedError(){return thisLedError;};
  void SetName(char* nm){fName=nm;};
  const char* GetName() const {return (const char*) fName;};
  Float_t EvNumRange[2];
  Float_t GetNDFpol0(){return NDFpol0;};
  Float_t GetChipdfpol0(){return chipdfpol0;};
 private:
  TString fName;
  Int_t nstb;
  Int_t row0;
  Int_t col0;
  Int_t chan;
  Int_t N0;//number of points in graph
  Float_t Sigma;//average error on points in graph
  TGraphErrors* gr;
  TF1* fit;
  Float_t NominalLedADC;// base LED average from which changes are calculated
  Float_t thisLedADC;// LED average for this run
  Float_t NominalGcor;
  Float_t thisLedError;
  Float_t thisfunError;
  Float_t chipdfpol0;//chi2 per DOF pol0 fit
  Float_t NDFpol0;// NDF pol0 fit
  Int_t Error;
  ClassDef(CellTDep,1);
};
#endif
