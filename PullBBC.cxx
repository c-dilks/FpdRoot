#include "PullBBC.h"
using namespace std;
ClassImp(PullBBC)

PullBBC::PullBBC()
{
  
  //calibration
  Float_t calibE[25]={0,1555,1489,1527,1536,1524,1550,1512,1528,1507,1520,1521,1533,1515,1511,1537,1523,1042, 984, 953,1053,1041, 981,1034,1066};

  Float_t calibW[25]={0,1452,1518,1457,1445,1454,1438,1509,1444,1415,1462,1413,1442,1460,1451,1426,1424,   0,1010,1321, 976,1014, 973,1022, 872};
  for(int jk=0;jk<25;jk++)
    {
      EWcal[0][jk]=calibE[jk]-1000;
      EWcal[1][jk]=calibW[jk]-1000;
    };

  //this assumes the following map is remapped from slots 6,8 and 10 to 1,3,5
  //http://www.star.bnl.gov/public/trg/TSL/Schematics/BBQ_Crate_Cable_Map.pdfc
  //assuming 
  PQt=0;
  scale=1.;//was 1.7
  int BBCslots0[17]={0,1,0,2,0,3,0,4,0,5,0,0,0,0,0,6,0};//6 BBC brds in these slots
  for(int jj=0;jj<17;jj++)BBCslots[jj]=BBCslots0[jj];
  int bADCchan0[6][4][2][4]={
    {
      { {1,7,2,8},  {}},  {{3,9,10,11}, {}}, //Ebrdsm 
      { {4,12,13,5}, {}},{{6,14,15,16}, {}}
    },
    {
      { {1,7,2,8},  {}},{{3,9,10,11}, {}}, //Wbrdsm
      { {4,12,13,5},{}},{{6,14,15,16}, {}}
    },
    {
      { {17,18,19,20},  {}},{{21,22,23,24}, {}}, //Ebrdlg
      { {17,18,19,20},  {}},{{21,22,23,24}, {}} //Wbrdlg
    },
    {},
    {},
    {}
  };//[myslot][card][half card][ch4]

  //copy initialized bADCchan0 array to bADCchan and bTACchan
  for(int slt=0;slt<6;slt++)
    {
      for(int st=0;st<4;st++)
	{	
	  for(int ch1=0;ch1<2;ch1++)
	    {
	      for(int ch0=0;ch0<4;ch0++)
		{
		  bADCchan[slt][st][ch1][ch0]=bADCchan0[slt][st][ch1][ch0];
		  bTACchan[slt][st][1-ch1][ch0]=bADCchan[slt][st][ch1][ch0];
		  int btac=bTACchan[slt][st][1-ch1][ch0];
		  //		  printf("slt=%d st=%d ch1=%d ch0=%d btac=%d\n",slt,st,ch1,ch0,btac);
		}
	    }
	};	
    };
  ADCmin=10;
  iwmax=-1;
  iemax=-1;
  VtxRange[0]=1;
  VtxRange[1]=16;
}
void PullBBC::GetADCTAC(Qt* pQt)
{
  PQt=pQt;
  memset(EWbADC,0,sizeof(EWbADC));
  memset(EWbADCT,0,sizeof(EWbADC));
  memset(EWbTAC,0,sizeof(EWbTAC));
  int hcard2=0;
  int card4=0;
  int ch4;
  int slot6;
  
  for(int chan=0;chan<32;chan++)
    {
      ch4=chan%4;
      int  hc2=(chan/4)&1;
      card4=chan/8;
      for(int slot=0;slot<7;slot++)
	{
	  int ew=0;
	  slot6=BBCslots[slot];
	  int badc=0;
	  int btac=0;
	  if(slot6>0 && slot6<7)
	    {
	      int val=PQt->qtadc[chan][slot-1][6];//BBQ crate
	      int tval=PQt->qttac[chan][slot-1][6];//BBQ crate

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
		      EWbADC[ew][badc]=val;
		      EWbADCT[ew][badc]=tval;
		    }
		  else if(hc2==1)
		    {
		      EWbTAC[ew][btac]=val-EWcal[ew][btac];
		    }
		} 
	    };
	}
    }
}

Float_t PullBBC::GetVertex()
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
  sE= SortByTime(0,ADCmin);
  sW= SortByTime(1,ADCmin);
  NAboveAdcmin[0]=sE.GetSize();
  NAboveAdcmin[1]=sW.GetSize();
  if(NAboveAdcmin[0]>0){iemax=sE[0];MaxtacE=EWbTAC[0][sE[0]];};
  if(NAboveAdcmin[1]>0){iwmax=sW[0];MaxtacW=EWbTAC[1][sW[0]];};
  vtx=-997;
  if(iemax<0&&iwmax<0)vtx=-10000;
  if(iemax<0&&iwmax>0)vtx=-9999;
  if(iwmax<0&&iemax>0)vtx=-9998;
  if(vtx<-9977)return vtx;
  vtx=(MaxtacE-MaxtacW)*scale;
  return vtx; 
};
Float_t PullBBC::SumADCinTRange(int iew,float dt1,float t2)
{
  float_t vertex=GetVertex();
  if((iew==0&&t2==-10000.)&&MaxtacE>0)t2=MaxtacE;
  if((iew==1&&t2==-10000.)&&MaxtacW>0)t2=MaxtacW;
  if(t2==-10000.)return -9999.;
  float sum=0;
  for(int j=VtxRange[0];j<VtxRange[1];j++)
    {
      if((EWbADC[iew][j]>ADCmin)&& (EWbTAC[iew][j]>t2-dt1)&&(EWbTAC[iew][j]<t2+dt1))sum+=EWbADC[iew][j];
    }
  return sum;
};
TArrayI PullBBC::SortByTime(int iew,int adcmn)
{

  if(iew<0 || iew>1)return TArrayI(0);
  Int_t ind[24];
  Float_t sortval[24]={};
  for(int j=0;j<24;j++)
    {
      if((EWbADC[iew][j]>adcmn)&&(j>=VtxRange[0])&&(j<=VtxRange[1]))
	{sortval[j]=EWbTAC[iew][j];}
      else {sortval[j]=-10000;};
    }
  TMath::Sort(24,sortval,ind);
  int cnt=0;
  for(int j=0;j<24;j++)
    {
      if(sortval[ind[j]]>0)cnt++;
    };
  return TArrayI(cnt,ind);
};
void PullBBC::PrintSorted()
{
  int adcmn=ADCmin;
  TString namew[2]={"East","West"};
  for (int iew=0;iew<2;iew++)
    {
      TArrayI a=SortByTime(iew,adcmn);
      for( int j=0;j<a.GetSize();j++)
	{
	  printf("%s:a[%d]=%d,adc=%d tdc=%d adcT=%d\n",
		 namew[iew].Data(), j,a[j],
		 EWbADC[iew][a[j]],EWbTAC[iew][a[j]],EWbADCT[iew][a[j]]);
	}
    }
}
void PullBBC::UpdateQTBBCinfo(QTBBCInfo* qtbbc,Qt* pqt)
{
    memset(qtbbc,0,sizeof(qtbbc));
    GetVertex(pqt);
    qtbbc->QTNE=NAboveAdcmin[0];
    qtbbc->QTNW=NAboveAdcmin[1];
    qtbbc->vertex=vtx;
    for(int j=0;j<NAboveAdcmin[0];j++)
      {
	int i=sE[j];
	if(i>0&&i<25)
	  {
	    qtbbc->QTEBBCInd[j]=i;
	    qtbbc->QTEBBCTAC[j]=EWbTAC[0][i];
	    qtbbc->QTEBBCADC[j]=EWbADC[0][i];
	  };
      };
    for(int j=0;j<NAboveAdcmin[1];j++)
      {
	int i=sW[j];
	if(i>0&&i<25)
	  {
	    qtbbc->QTWBBCInd[j]=i;
	    qtbbc->QTWBBCTAC[j]=EWbTAC[1][i];
	    qtbbc->QTWBBCADC[j]=EWbADC[1][i];
	  };
      };
};
 
