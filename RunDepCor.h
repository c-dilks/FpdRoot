#ifndef RunDepCor_
#define RunDepCor_
#include "CellTDep.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TObjArray.h"
#include "TSystem.h"

/*
  An instance of the RunDepCor  class will be streamable and will be an entry in a TObjArray 
of similar instances. The name will be RunDepCor[run number].
The instance will contain.
an array of CellTDep elements.

 */

class RunDepCor : public TObject
{
 public:
  RunDepCor();
  RunDepCor(Long_t Rnum,TString ledpath,CalibStr* initgCor,RunDepCor* base=0);
  void MakeDepCor(TFile* ledfile,CalibStr* initgCor,RunDepCor* base=0);
  TFile* FillLedADC(TString adcTrLpath);
  Float_t GetCor(Int_t Segment_number, Int_t evNumber, Int_t NSTB,
		 Int_t Row, Int_t Col, CalibStr* GCor=0,Int_t mode=1);
  //  CalibStr* Initgcor;//pointer to Initial gaincorr 
  TObjArray* array;
  TString dname[4][34][17];
  CellTDep* Cdep(Int_t nstb,Int_t row0,Int_t col0);//pointer to nstb,row0,col0 time dependence element
  Int_t Location[4][34][17];
  Bool_t Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0);
  ~RunDepCor();
  void SetName(char* nm){fName=nm;};
  const char* GetName() const {return (const char*) fName;};
  Long_t RunNumber;
  Long_t BaseRunNumber;
 private:
  TString fName;
  ClassDef(RunDepCor,1);
  Int_t nrows[4];

};
#endif
