#ifndef Vertex_
#define Vertex_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"

class Vertex : public TObject
{
 public:
  Vertex();
  ~Vertex();
  Float_t GetBbcVertex(UChar_t* p_Bbc);
  void  InitializeBbc();
  Int_t east_q_map[24];
  Int_t east_t_map[16];
  Int_t west_q_map[24];
  Int_t west_t_map[16];
  Int_t offset[16][2];
  Int_t Bbc_threshold;
  Float_t zvt;
  Int_t ewdiff;
  Int_t ewsum;
  Int_t maxqe;//maximum adc
  Int_t maxqw;
  Int_t ieqmax;//maximum adc channel
  Int_t iwqmax;
  Int_t iemax;//maximum tac channel
  Int_t iwmax;
  Int_t maxtace;//maximum tac
  Int_t maxtacw;
  Int_t qe;//maximum tac channel adc
  Int_t qw;

 private:
  ClassDef(Vertex,1);
};

#endif

