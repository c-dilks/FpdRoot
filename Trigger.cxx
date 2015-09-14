#include "Trigger.h"

ClassImp(Trigger)

  Trigger::Trigger(UInt_t* p_Bbcl1)
{
  strig=l1sume=l1sumw=0;
  l1maxtace=l1maxtacw=0;
  for(Int_t i=0;i<5;i++)BBcSums[i]=0;
  for(Int_t i=0;i<4;i++)
    {
      itac=6-i*2;
      iadc=7-i*2;
      l1tac[i]=(int) fmod(p_Bbcl1[itac],256);
      l1adc[i]=(int) fmod(p_Bbcl1[iadc],256);
      if(i==0 || i==2)
	{
	  l1sume+=l1adc[i];
	  if((l1tac[i] > l1maxtace) && (l1tac[i]<=250))
	    {
	      l1maxtace=l1tac[i];
	    };
	}
      else
	{
	  l1sumw+=l1adc[i];
	  if((l1tac[i]>l1maxtacw) && (l1tac[i]<=250))
	    {
	      l1maxtacw=l1tac[i];
	    };
	};
    };
  if(l1maxtace>20 && l1maxtace< 175 && l1maxtacw>20 && l1maxtacw<175 
     && l1sume>10 && l1sumw>10)
    {
      strig=1;
    };
  BBcSums[0]=(char) strig;
  BBcSums[1]=(char) l1sume;
  BBcSums[2]=(char) l1sumw;
  BBcSums[3]=(char) l1maxtace;
  BBcSums[4]=(char) l1maxtacw;
  if(l1sume>255)BBcSums[1]=255;
  if(l1sumw>255)BBcSums[2]=255;
};
