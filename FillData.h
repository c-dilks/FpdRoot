#ifndef FillData_
#define FillData_
#include <stdio.h>
#include <stdlib.h>
#include "TObject.h"
#include "Filentry.h"

class FillData : public TObject
{
 public:
  FillData(Int_t);
  Int_t FillNumber;
  Int_t FirstRun;
  Int_t LastRun;
  Int_t BluePattern[120];
  Int_t YellowPattern[120];
  Filentry* BluePatternFile;
  Filentry* YellowPatternFile;
  Int_t blue_offset;
  Int_t yellow_offset;
  Int_t Lum_r1;
  Bool_t BlueGood;
  Bool_t YellowGood;
  Int_t NumberDropped;
  ULong_t Hash() const {return (ULong_t) FirstRun;};
  Bool_t  IsEqual(const TObject *obj) 
    const { return FirstRun == ((FillData*) obj)->FirstRun; }
  Bool_t  IsSortable() const { return kTRUE; }
  Int_t   Compare(const TObject *obj) 
    const { if (FirstRun > ((FillData*)obj)->FirstRun)
      return 1;
    else if (FirstRun < ((FillData*) obj)->FirstRun)
      return -1;
    else
      return 0; };
  Bool_t kicked(Int_t);
  Int_t DropLst[30];

 private:
  ClassDef(FillData,1);
};

#endif
