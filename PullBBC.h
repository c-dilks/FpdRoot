#ifndef PullBBC_
#define PullBBC_
#include "Qt.h"
#include "TArray.h"
#include "TMath.h"

class QTBBCInfo : public TObject
{
 public:
  QTBBCInfo(){};
  Char_t QTNE;
  Char_t QTNW;
  Char_t QTEBBCInd[24];
  Char_t QTWBBCInd[24];
  Short_t QTEBBCTAC[24];
  Short_t QTWBBCTAC[24];
  Short_t QTEBBCADC[24];
  Short_t QTWBBCADC[24];
  Float_t vertex;
  void Print()
  {
    printf("QTNE= %d vertex=%f \n ",QTNE,vertex);
    for(int j=0;j<QTNE;j++)
      {
	printf("E%d: Ind=%d TAC=%d ADC=%d \n",j, QTEBBCInd[j],QTEBBCTAC[j],QTEBBCADC[j]);};
    printf("QTNW= %d \n",QTNE);
    for(int j=0;j<QTNW;j++)
      {
	printf("W%d: Ind=%d TAC=%d ADC=%d \n",j, QTWBBCInd[j],QTWBBCTAC[j],QTWBBCADC[j]);
      }
  }
  void SetZero()
  {
    QTNE=0;
    QTNW=0;
    for(int j=0;j<24;j++){QTEBBCInd[j]=0;QTEBBCTAC[j]=0;QTEBBCADC[j]=0;QTWBBCInd[j]=0;QTWBBCADC[j]=0;QTWBBCTAC[j]=0;vertex=0.;};
  }
 private:
  ClassDef(QTBBCInfo,1);
};

class PullBBC : public TObject
{
 public:
  PullBBC();
  void GetADCTAC(Qt* pQt);
  
  int EWbADC[2][33];//ADC...[E/W][channel] starts with channel 1
  int EWbTAC[2][33];//TAC...[E/W][channel] starts with channel 1
  int EWbADCT[2][33];//ADC...[E/W][channel]  (QT time of EWbADC)
  Qt* PQt;
  int BBCslots[17];
  int bTACchan[6][4][2][4];
  int bADCchan[6][4][2][4];
  Float_t GetVertex();
  Float_t GetVertex(Qt* pQt)
  {
    GetADCTAC(pQt);
    return GetVertex();
  }
  TArrayI sE;
  TArrayI sW;
  TArrayI SortByTime(int iew,int adcmn);
  void PrintSorted();
  Float_t ADCmin;
  Float_t MaxtacW;
  Float_t MaxtacE;
  Int_t iwmax;
  Int_t iemax;
  Float_t scale;//conversion from delta t to distance
  Float_t vtx;
  Float_t EWcal[2][25];  
  Int_t VtxRange[2];
  Int_t NAboveAdcmin[2];//number for E/W
  Float_t SumADCinTRange(int iew,float dt1=200,float t2=-10000.);// t1=t2+-dt1 and if t2=-1000 t2 is max time
  void UpdateQTBBCinfo(QTBBCInfo* qtbbc,Qt* pqt);
 private:
  ClassDef(PullBBC,1);
  
};
#endif
