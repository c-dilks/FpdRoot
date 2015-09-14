
#ifndef RfOneD
#define RfOneD

#include <iostream>
#include <stdio.h>     
#include  "TObjFloat_t.h"
#include  "TVector.h"
#include  "TSortedList.h"
#include "TList.h"
#include "TArrayI.h"
class fOneD 
{
 public:
  fOneD();
  fOneD(TVector&);
  ~fOneD();
  Int_t* order;
  Int_t  NumberElements;
  TArrayI* TA;
  //  ClassDef(fOneD,1);
};

#endif
