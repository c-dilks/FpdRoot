#ifndef Trigger_
#define Trigger_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"

class Trigger : public TObject
{
 public:
  Trigger(UInt_t* P_Bbcl1);
  Int_t strig;
  Int_t l1sume;
  Int_t l1sumw;
  Int_t l1maxtace;
  Int_t l1maxtacw;
  Int_t iadc,itac;
  Int_t l1tac[4];
  Int_t l1adc[4];
  Char_t  BBcSums[5];
 private:
  ClassDef(Trigger,1);
};

#endif

