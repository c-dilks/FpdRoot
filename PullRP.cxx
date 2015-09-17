#include "PullRP.h"
//using namespace std;
ClassImp(PullRP)

namespace
{
  enum ew_enum {kE,kW}; // east west
  enum io_enum {kI,kO}; // inner outer
  enum ud_enum {kU,kD}; // up down
  enum ns_enum {kN,kS}; // north south
  const Int_t MXQ_CRATE_CONF_NUM = 4; // root12fms crate num - 1 (see Qt.h)
};

PullRP::PullRP()
{
  // ADC thresholds
  ADCmin[0] = 97; // first minimum @ 111
  ADCmin[1] = 100; // first minimum @ 110
  ADCmin[2] = 100; // first minimum @ 113
  ADCmin[3] = 84; // first minimum @ ~90
  ADCmin[4] = 103; // first minimum @ 120
  ADCmin[5] = 83; // first minimum @ 90
  ADCmin[6] = 97; // first minimum @ 111
  ADCmin[7] = 97; // first minimum @ 116
  ADCmin[8] = 100; // first minimum @ 127 
  ADCmin[9] = 98; // first minimum @ 117
  ADCmin[10] = 102; // first minimum @ 124
  ADCmin[11] = 103; // first minimum @ 122
  ADCmin[12] = 131; // first minimum @ 175
  ADCmin[13] = 83; // first minimum @ 87
  ADCmin[14] = 110; // first minimum @ 130
  ADCmin[15] = 91; // first minimum @ ~102



  //calibration -- shift TAC peaks to ~1000 counts
  TACcalib[0]  = 553;
  TACcalib[1]  = 537;
  TACcalib[2]  = 528;
  TACcalib[3]  = 600;
  TACcalib[4]  = 628;
  TACcalib[5]  = 446;
  TACcalib[6]  = 491;
  TACcalib[7]  = 521;
  TACcalib[8]  = 338;
  TACcalib[9]  = 290;
  TACcalib[10] = 440;
  TACcalib[11] = 163;
  TACcalib[12] = 704;
  TACcalib[13] = 768;
  TACcalib[14] = 642;
  TACcalib[15] = 757;

  // uncomment to disable ADC and TAC calibrations
  //for(int i=0; i<16; i++) { ADCmin[i]=1; TACcalib[i]=0; };



  // PP001 QT board is in slot 8 of MXQ crate; 32 inputs (16 ADC & 16 TAC)
  // 2 16-bit outputs sent to PP101 DSM board in slot 6 (3A) of MIX crate
  // DSM inputs are on channels 0 and 1, which is mapped to FPDE[3] and FPDE[2] in 
  // the h111 ntuple (since the mapping is: DSM input ---- 0 1 2 3 4 5 6 7
  //                                        array index -- 3 2 1 0 7 6 5 4
  // cable maps:
  // - MXQ: http://www.star.bnl.gov/public/trg/TSL/Schematics/MXQ_Crate_Cable_Map.pdf
  // - MIX: http://www.star.bnl.gov/public/trg/TSL/Schematics/MIX_Crate_Cable_Map.pdf
  // crate layouts:
  // - MXQ: http://www.star.bnl.gov/public/trg/TSL/Schematics/MXQ_Crate_Layout.pdf
  // - MIX: http://www.star.bnl.gov/public/trg/TSL/Schematics/MIX_Crate_Layout.pdf
  PQt=0;
  scale=1.;//was 1.7
  //                0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6
  int RPslots0[17]={0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0};// [PP001 qt slot - 5] --> {"myslot"}
  for(int jj=0;jj<17;jj++)RPslots[jj]=RPslots0[jj];

  int bADCchan0[6][4][2][4]; // [myslot] [card] [half card] [ch4]
  memset(bADCchan0,0,sizeof(bADCchan0));
  memset(bADCchan,0,sizeof(bADCchan0));
  memset(bTACchan,0,sizeof(bADCchan0));

  bADCchan0[0][0][0][0] = EiunToIdx(kE,kO,kU,kN);
  bADCchan0[0][0][0][1] = EiunToIdx(kE,kO,kU,kS);
  bADCchan0[0][0][0][2] = EiunToIdx(kE,kO,kD,kS);
  bADCchan0[0][0][0][3] = EiunToIdx(kE,kO,kD,kN);

  bADCchan0[0][1][0][0] = EiunToIdx(kW,kO,kU,kS);
  bADCchan0[0][1][0][1] = EiunToIdx(kW,kO,kU,kN);
  bADCchan0[0][1][0][2] = EiunToIdx(kW,kO,kD,kN);
  bADCchan0[0][1][0][3] = EiunToIdx(kW,kO,kD,kS);

  bADCchan0[0][2][0][0] = EiunToIdx(kE,kI,kD,kS);
  bADCchan0[0][2][0][1] = EiunToIdx(kE,kI,kD,kN);
  bADCchan0[0][2][0][2] = EiunToIdx(kE,kI,kU,kN);
  bADCchan0[0][2][0][3] = EiunToIdx(kE,kI,kU,kS);

  bADCchan0[0][3][0][0] = EiunToIdx(kW,kI,kD,kN);
  bADCchan0[0][3][0][1] = EiunToIdx(kW,kI,kD,kS);
  bADCchan0[0][3][0][2] = EiunToIdx(kW,kI,kU,kS);
  bADCchan0[0][3][0][3] = EiunToIdx(kW,kI,kU,kN);

  //copy initialized bADCchan0 array to bADCchan and bTACchan
  for(int slot_i=0;slot_i<6;slot_i++)
  {
    for(int card_i=0;card_i<4;card_i++)
    {	
      for(int half_card_i=0;half_card_i<2;half_card_i++)
      {
        for(int ch_i=0;ch_i<4;ch_i++)
        {
          bADCchan[slot_i][card_i][half_card_i][ch_i]=bADCchan0[slot_i][card_i][half_card_i][ch_i];
          bTACchan[slot_i][card_i][1-half_card_i][ch_i]=bADCchan[slot_i][card_i][half_card_i][ch_i];
          int btac=bTACchan[slot_i][card_i][1-half_card_i][ch_i];
          /*
          printf("slot_i=%d card_i=%d half_card_i=%d ch_i=%d badc=%d btac=%d\n",
            slot_i,card_i,half_card_i,ch_i,
            bADCchan[slot_i][card_i][half_card_i][ch_i],
            bTACchan[slot_i][card_i][half_card_i][ch_i]);
          */
        };
      };
    };	
  };
  iwmax=-1;
  iemax=-1;
  VtxRange[0]=0;
  VtxRange[1]=15;
}


void PullRP::GetADCTAC(Qt* pQt)
{
  PQt=pQt;
  memset(bADC,0,sizeof(bADC));
  memset(bADCT,0,sizeof(bADC));
  memset(bTAC,0,sizeof(bTAC));
  int hcard2=0;
  int card4=0;
  int ch4;
  int slot6;

  for(int chan=0;chan<32;chan++)
  {
    ch4=chan%4;
    int  hc2=(chan/4)&1; // half card
    card4=chan/8;

    for(int slot=0;slot<7;slot++)
    {
      int ew=0;
      slot6=RPslots[slot];
      int badc=0;
      int btac=0;
      if(slot6>0 && slot6<7)
      {
        int val=PQt->qtadc[chan][slot-1][MXQ_CRATE_CONF_NUM];
        int tval=PQt->qttac[chan][slot-1][MXQ_CRATE_CONF_NUM];

        if(val<0)val=0;
        if(val>0 ||true)
        {
          badc=bADCchan[slot6-1][card4][hc2][ch4];
          btac=bTACchan[slot6-1][card4][hc2][ch4];

          ew=0;
          if(slot6==2){ew=1;}
          else
            if(slot6==3&&card4>1){ew=1;};
          if(hc2==0)
          {
            bADC[badc]=val;
            bADCT[badc]=tval;
          }
          else if(hc2==1)
          {
            bTAC[btac]=val-TACcalib[btac];
          }
        } 
      };
    }
  }
}

Float_t PullRP::GetVertex()
{
  if(PQt==0)
  {
    vtx=-1001.;
    return vtx;
  };
  MaxtacE=0;
  MaxtacW=0;
  iemax=-1;
  iwmax=-1;
  sE= SortByTime(0); // index of sorted array
  sW= SortByTime(1);
  NAboveAdcmin[0]=sE.GetSize();
  NAboveAdcmin[1]=sW.GetSize();
  if(NAboveAdcmin[0]>0){iemax=sE[0];MaxtacE=bTAC[sE[0]];};
  if(NAboveAdcmin[1]>0){iwmax=sW[0];MaxtacW=bTAC[sW[0]];};
  vtx=-997;
  if(iemax<0&&iwmax<0)vtx=-10000;
  if(iemax<0&&iwmax>0)vtx=-9999;
  if(iwmax<0&&iemax>0)vtx=-9998;
  if(vtx<-9977)return vtx;
  vtx=(MaxtacE-MaxtacW)*scale;
  return vtx; 
};

Float_t PullRP::SumADCinTRange(int iew,float dt1,float t2)
{
  float vertex=GetVertex();
  if((iew==0&&t2==-10000.)&&MaxtacE>0)t2=MaxtacE;
  if((iew==1&&t2==-10000.)&&MaxtacW>0)t2=MaxtacW;
  if(t2==-10000.)return -9999.;
  float sum=0;
  for(int j=VtxRange[0];j<=VtxRange[1];j++)
  {
    if((int)((j>>3)&0x1)==iew)
    {
      if( (bADC[j]>ADCmin[j]) && (bTAC[j]>t2-dt1) && (bTAC[j]<t2+dt1) )
        sum += bADC[j];
    };
  };
  return sum;
};

TArrayI PullRP::SortByTime(int iew)
{
  if(iew<0 || iew>1)return TArrayI(0);
  Int_t ind[16];
  Float_t sortval[16]={};
  for(int j=0;j<16;j++)
  {
    if( (int)((j>>3)&0x1)==iew && (bADC[j]>ADCmin[j]) && (j>=VtxRange[0]) && (j<=VtxRange[1]) )
      sortval[j]=bTAC[j];
    else 
      sortval[j]=-10000;
  };
  TMath::Sort(16,sortval,ind);
  int cnt=0;
  for(int j=0;j<16;j++)
  {
    if(sortval[ind[j]]>0)cnt++;
  };
  return TArrayI(cnt,ind);
};

void PullRP::PrintSorted()
{
  TString namew[2]={"East","West"};
  for (int iew=0;iew<2;iew++)
  {
    TArrayI a=SortByTime(iew);
    for(int j=0;j<a.GetSize();j++)
    {
      printf("%s:a[%d]=%d,adc=%d tdc=%d adcT=%d\n",
          namew[iew].Data(), j,a[j],
          bADC[a[j]],bTAC[a[j]],bADCT[a[j]]);
    };
  };
};

void PullRP::UpdateQTRPinfo(QTRPInfo* qtrp,Qt* pqt)
{
  memset(qtrp,0,sizeof(qtrp));
  GetVertex(pqt);
  qtrp->NE=NAboveAdcmin[0];
  qtrp->NW=NAboveAdcmin[1];
  qtrp->vertex=vtx;
  for(int j=0;j<NAboveAdcmin[0];j++)
  {
    int i=sE[j];
    if(i>=0&&i<16)
    {
      qtrp->RPE_Idx[j]=i;
      qtrp->RPE_TAC[j]=bTAC[i];
      qtrp->RPE_ADC[j]=bADC[i];
    };
  };
  for(int j=0;j<NAboveAdcmin[1];j++)
  {
    int i=sW[j];
    if(i>=0&&i<16)
    {
      qtrp->RPW_Idx[j]=i;
      qtrp->RPW_TAC[j]=bTAC[i];
      qtrp->RPW_ADC[j]=bADC[i];
    };
  };
};

void PullRP::IdxToEiun(Int_t idx0, Int_t &ew0, Int_t &io0, Int_t &ud0, Int_t &ns0)
{
  // Idx = 4 bits represinging [east/west] [inner/outer] [up/down] [north/south]
  ew0 = (idx0 >> 3) & 0x1;
  io0 = (idx0 >> 2) & 0x1;
  ud0 = (idx0 >> 1) & 0x1;
  ns0 = (idx0 >> 0) & 0x1;
};

Int_t PullRP::EiunToIdx(Int_t ew0, Int_t io0, Int_t ud0, Int_t ns0)
{
  if(abs(ew0)<2 && abs(io0)<2 && abs(ud0)<2 && abs(ns0)<2)
    return 8*ew0+4*io0+2*ud0+ns0;
};
