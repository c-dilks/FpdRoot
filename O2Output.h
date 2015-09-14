#ifndef O2Output_
#define O2Output_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "poutTree.h"
#include "LVec.h"
#include "FilesSet.h"
#include "TObject.h"
#include "Cell.h"
#include "Qt.h"
#include "TObjArray.h"
#include "TVector.h"
#include "CalibStr.h"
#include "TIterator.h"
#include "TMatrix.h"
#include "CalibStr.h"
#include "RunDepCor.h"


class O2Output : public TObject {
  public:
    O2Output();
    ~O2Output();
    Int_t OutputTwoTr(char* filelist, FilesSet* p_files);
    Int_t OutputMycells(char* filelist, FilesSet* p_files, Bool_t doSingle=false);
    Int_t NumberEventsToProcess;
    TString OutFileName;
    Int_t keyIndex;
    poutTree* pout;  
    Float_t maxsep; 
    Float_t minsep;
    Float_t ECutMin;
    Float_t ECutMax;
    Float_t ZMax;
  private:
    ClassDef(O2Output, 1);
};

#endif
