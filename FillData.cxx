#include "FillData.h"

ClassImp(FillData)
  
  FillData::FillData(Int_t Fnumber)
{
  FillNumber=Fnumber;
  NumberDropped=22;
  Int_t* dtmp;
  Int_t dtmp1[22]={0,14,54,119,  
		      25,26,27,28,29,30,31,32,33,105,106,107,108,109,110,111,112,113};
  Int_t dtmp2[20]={20,60,        
		     31,32,33,34,35,36,37,38,39,111,112,113,114,115,116,117,118,119};
  Int_t dtmp3[23]={5,20,45,53,60,
		    31,32,33,34,35,36,37,38,39,111,112,113,114,115,116,117,118,119};
  Int_t dtmp4[21]={5,20,60,       
		    31,32,33,34,35,36,37,38,39,111,112,113,114,115,116,117,118,119};
  Int_t dtmp5[18]={        
		   31,32,33,34,35,36,37,38,39,111,112,113,114,115,116,117,118,119};

  if(FillNumber<7700)
    {
      dtmp=dtmp1;
      NumberDropped=22;
    }
  else if((FillNumber<7785))
    {
      dtmp=dtmp2;
      NumberDropped=20;
    }
  else if((FillNumber>=7785)&&(FillNumber<=7788))
  {dtmp=dtmp3;
  NumberDropped=23;}
	  else if(FillNumber>7788 && FillNumber<7833)
  {
    dtmp=dtmp4;
  NumberDropped=21;}
	  else
  {
    dtmp=dtmp5;
    NumberDropped=18;};
	  
  for(Int_t i=0;i<NumberDropped;i++)DropLst[i]=dtmp[i];
  FirstRun=1000000000;
  LastRun=0;
  for(Int_t i=0;i<120;i++)
    {
      BluePattern[i]=0;
      YellowPattern[i]=0;
    };
  if(FillNumber<7700)
    {
      blue_offset=34;
      yellow_offset=114;
    }
  else if(FillNumber<7958)
    {
      blue_offset=40;
      yellow_offset=120;
    }
  else if(FillNumber<16243)// up to including Run 12
    {      
      blue_offset=120;
      yellow_offset=40;

    }
  else if(FillNumber<18536)//before Run 15
    {
      blue_offset=119;
      yellow_offset=39;
    }
  else 
    {
      blue_offset=120;
      yellow_offset=40;
    };

  BlueGood=false;
  YellowGood=false;
  BluePatternFile=0;
  YellowPatternFile=0;
  Lum_r1=0;

  /*
    pattern=0 : no valid spin
    pattern=1 : valid +
    pattern=-1: vallid -
   */
};
Bool_t FillData::kicked(Int_t bnch)
{
  for(Int_t i=0;i<NumberDropped;i++)
    {
      if(bnch==DropLst[i])return true;
    };
  return false;
};
