#ifndef RunData_
#define RunData_
#include <stdio.h>
#include <stdlib.h>
#include "TObject.h"
#include "TString.h"
#include "FilesSet.h"
#include "Fill.h"
#include "CalibStr.h"

class RunData : public TObject
{
 public:
  RunData(Int_t ,Fill*,FilesSet*);
  FilesSet* p_Files;
  Int_t RunNumber;
  Int_t FillNumber;
  const char* GetName();
  TString charname;
  Float_t R1Lum;
  Float_t ErrR1Lum;
  Float_t R2Lum;
  Float_t ErrR2Lum;
  TString RunPedFileName;
  CalibStr* RunPed;
  Int_t GetRunPed();
  Int_t error;
 private:
  ClassDef(RunData,1);
};
#endif
