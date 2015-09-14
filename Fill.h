#ifndef Fill_
#define Fill_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "FillData.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TString.h"
#include "TComplex.h"
#include "FilesSet.h"
#include "TVector.h"
#include "TVectorT.h"

class Fill : public TObject
{
 public:
    Fill(FilesSet*,Int_t,Int_t,Int_t,Int_t); 
    ~Fill();
    Filentry* FillFile;
    FilesSet* p_Files;
    Int_t SetFillNumber(Int_t);
    Float_t* GetYellow();
    Float_t* GetBlue();
    Int_t FirstRun;
    Int_t LastRun;
    Int_t FirstFill;
    Int_t LastFill;
    TObjArray* pFillArray;
    Int_t CurrentRunNumber;
    Int_t ReadFillFile();
    Int_t GetPattern();
    TString  ErrorMessage;
    Int_t error;
    Int_t Print();
    Int_t SetFillNumberforRun(Int_t);
    Int_t GetFillNumber(){return FillNumber;};
    FillData* pfd;
    Int_t BlueCnt;
    Int_t YellowCnt;
    Int_t ListFills();
    Int_t BlueSpin(Int_t);
    Int_t YellowSpin(Int_t);

    Float_t BlueEvenOddCor();
    Float_t YellowEvenOddCor();
    Float_t* pBlueEvenOddCor(Float_t);
    Float_t* pYellowEvenOddCor(Float_t);
    Float_t ycorr[4];
    Float_t bcorr[4];
    Float_t Dot(TVectorT<float>,TVectorT<float>);
    TVectorT<float> Mult(TVectorT<float>,TVectorT<float>);
    TVectorT<float> YellowVec();
    TVectorT<float> BlueVec();
    TVectorT<float> Mask();
      
 private:
    TString Bluestring[2];
    TString Yellowstring[2];
    Float_t blueSpin[120];
    Float_t yellowSpin[120];
    Int_t FillNumber;
    ClassDef(Fill,1);
};

#endif
