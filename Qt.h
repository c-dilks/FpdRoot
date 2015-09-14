#ifndef Qt_
#define Qt_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "FilesSet.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TH2F.h"

class Qt : public TObject
{
 public:
  Qt(FilesSet* p_files);
  ~Qt();
  // Large 17x34=578
  // Small 12x24=288
  Int_t pr[578][4];// [det channels][nstb] counts from 1
  Int_t rr[578][4];// [det channels][nstb] counts from 1
  Int_t cr[578][4];// [det channels][nstb] counts from 1
  Int_t rcrate[16][20][2][2];//[cr-1],[rr-1],[pr-1],[ns] counts from 1
  Int_t rslot[16][20][2][2];//[cr-1],[rr-1],[pr-1],[ns] counts from 1
  Int_t rdcd[16][20][2][2];//[cr-1],[rr-1],[pr-1],[ns] counts from 0
  Int_t rdch[16][20][2][2];//[cr-1],[rr-1],[pr-1],[ns] counts from 0

  Int_t dcrate[578][4];// [det channels][nstb] counts from 1
  Int_t dslot[578][4];// [det channels][nstb] counts from 1
  Int_t ddcd[578][4];// [det channels][nstb] counts from 0
  Int_t ddch[578][4];// [det channels][nstb] counts from 0
  
  Int_t detchan[32][16][7];//[chan][slot][crate] counts from 1
  Int_t detnstb[16][7];//[slot][crate] counts from 1
  
  Int_t NTrigMap[32][16];//[chan][slot]
  Int_t TrigMapslot[32][16][4];//[chan][slot][num qt8];
  Int_t TrigMapqt8[32][16][4];//[chan][slot][num qt8];
  void DefineTrigMap();
  void DefineRealTrigMap();
  void AddToTrigMap(const char* slot,Int_t ch_min,Int_t ch_max, 
		    const char* s0="",const char* s1="",const char* s2="",const char* s3="");
  void AddToTrigMapSq(const char* slot,Int_t ch_min,Int_t ch_max, 
		    const char* s0="",Int_t nqt8=0);
  Float_t GetTrgAdj(Int_t ch, Int_t slot, Int_t crate); //for simulating trigger response to gain reduction

  void printDSM();
  void printADC() {for(int i=0; i<4; i++) qtMat->At(i)->Print();};
  void printQt();

  Int_t COL_NUM[4]; //[nstb]
  Int_t ROW_NUM[4];
  TObjArray* qtMat;
  TH2F* qtHist[4];           //[nstb]
  TH2F* qtsumadcmcnt;      
  Int_t qtadc[32][16][7];    //[chan][slot][crate]
  Int_t qttac[32][16][7];
  Int_t qtnht[16][7];        //[slot][crate]
  Int_t qt8sum[4][16][7];    //[qt8][slot][crate]
  Int_t qtsum[16][7];
  Int_t qtht[16][7];
  Int_t qthtadr[16][7];
  Int_t fmssum[5];
  Int_t errfff;
  Int_t QtError;
  Filentry* p_qtmap;
  Filentry* p_qtmap2pp;
  Int_t qtcrate_off[14];
  Int_t qtadd_off;
  Int_t qtmaxch[4];
  Bool_t decodeQT(Int_t nQTdata,UInt_t* QTdata, Int_t IDNEVT);
  Bool_t encodeQT(TMatrix** m_adc,UInt_t* nQTdata,UInt_t* QTdata);
  void initQT();
  TMatrix* tm(Int_t instb);
  Int_t ndata;
  Int_t NFmsAdc;
  Int_t FmsSumAdc;
  Int_t sumadc;
  Int_t mcnt;
  void EnableqtHist();
  Bool_t trgadj;

  Bool_t decodeDSM(UChar_t* DSMdata,UInt_t* dsm2,UInt_t dsm);

  //output from 40 qt boards
  Int_t qtdsm0sum[10][4];//[qt32][crate] (10/16 qt32's are in trigger)
  Int_t qtdsm0sum2[10][4];//[qt32][crate]
  Int_t qtdsm0dcsum[4][10][4];//[qt8][qt32][crate] dc qt8 sum
  Int_t qtdsm0ht[10][4];  //ht value
  Int_t qtdsm0htadr[10][4];  //ht adr

  //output from 12 L1 dms
  Int_t qtdsm1sum[2][4][3];//[2x16bit blocks][L0dsm][L1dsm]
  Int_t qtdsm1clsum[12];//[L0dsm] cluster sum
  Int_t qtdsm1clsum_apnd[12];//[L0dsm] appended with 4th dcsum
  Int_t qtdsm1cladr[12];//extended version of qtdsm0htadr
  Int_t qtdsm1dcsuma[12];//qt8 sum to L2
  Int_t qtdsm1dcsumd[12];//qt8 sum to L2
  Int_t qtdsm1ht[12];//ht threshold bit

  Int_t qtdsm2_trgbits[3]; //[L1dsms] L1 dsm output bits small, south large and north large
  Int_t qtdsm_trgbits;  //output bits from dsm2 to last dsm

 private:
  Int_t mm[4];
  void set_mm_order(Int_t m1=1,Int_t m2=3,Int_t m3=2,Int_t m4=4)
    {
      mm[0]=m1;mm[1]=m2;mm[2]=m3;mm[3]=m4;
    };
  ClassDef(Qt,1);
};
#endif
