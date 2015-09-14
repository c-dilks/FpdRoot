#ifndef FpdMap_
#define FpdMap_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "TObject.h"
#include "TMatrix.h"


class FpdMap: public TObject
{
 public:
  FpdMap(Int_t );
  Int_t RunNum;
  TMatrix GetMatrix(Int_t , Int_t );
 private:
  ClassDef(FpdMap,1);
};

#endif
