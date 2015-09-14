#ifndef Filentry_
#define Filentry_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include "TObject.h"
#include "TString.h"

class Filentry : public TObject
{
 public: 

  Filentry(TString,TString,TString);
  ~Filentry();
  Filentry(){};
  TString Directory;
  TString FileName;
  const char*  Path();
  TString Type;
  TString Description;
  Int_t Error;
  Bool_t exists;
  Bool_t isOpen;
  Bool_t check();
  Bool_t open();
  Bool_t close();
  FILE* fp;//! not streamed
  Int_t Print();
  TString GetRunpedFileName(Int_t);
  TString GetRunFillBlueName(Int_t,Int_t);
  TString GetRunFillYellowName(Int_t,Int_t);
 private:
  ClassDef(Filentry,1);
 TString  path;
};

#endif
