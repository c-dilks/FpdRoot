#include "ClusterQt.h"
ClassImp(ClusterQt)
using namespace std;

ClusterQt::ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_)
{
  Init();
  pQt=p_Qt;
  crate=crate_;
  slot=slot_;
  qthtadr=pQt->qthtadr[slot][crate];//get high tower for the input slot from Qt
  qtht=pQt->qtht[slot][crate];
  detnstb=pQt->detnstb[slot][crate];
  detchan=0;//counts from 1
  //cout<<"crate="<<crate<<" slot="<<slot<<" qthtadr="<<qthtadr<<" qtht="<<qtht<<endl;

  if(qtht!=0 && qthtadr!=-1){
    detchan=pQt->detchan[qthtadr][slot][crate];
    //cout<<"test crt="<<crate<<" slot="<<slot<<" qtch="<<qthtadr<<" detchan="<<detchan<<endl;
    GetClusterSum();
  };
};


ClusterQt::ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Bool_t qt8trunc_)
{
  Init();
  pQt=p_Qt;
  crate=crate_;
  slot=slot_;
  qthtadr=pQt->qthtadr[slot][crate];//get high tower for the input slot from Qt
  qtht=pQt->qtht[slot][crate]/32;
  detnstb=pQt->detnstb[slot][crate];
  detchan=0;
  qt8trunc=qt8trunc_;
  //cout<<"crate="<<crate<<" slot="<<slot<<" qthtadr="<<qthtadr<<" qtht="<<qtht<<endl;

  if(qtht>0 && qthtadr!=-1){
    detchan=pQt->detchan[qthtadr][slot][crate];
    GetClusterSum();
  }
};

ClusterQt::ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Int_t ch_adr)
{
  Init();
  pQt=p_Qt;
  crate=crate_;
  slot=slot_;
  qthtadr=ch_adr;//Build the cluster around the input channel, high tower or not
  qtht=pQt->qtadc[qthtadr][slot][crate];
  detnstb=pQt->detnstb[slot][crate];
  detchan=pQt->detchan[qthtadr][slot][crate]; 
  
  GetClusterSum();
};

ClusterQt::ClusterQt(Qt* p_Qt,Int_t crate_,Int_t slot_,Int_t ch_adr,Bool_t qt8trunc_)
{
  Init();
  pQt=p_Qt;
  crate=crate_;
  slot=slot_;
  qthtadr=ch_adr;//Build the cluster around the input channel, high tower or not
  qtht=pQt->qtadc[qthtadr][slot][crate]/32;
  detnstb=pQt->detnstb[slot][crate];
  detchan=pQt->detchan[qthtadr][slot][crate]; 
  qt8trunc=qt8trunc_;
  
  GetClusterSum();
};

ClusterQt::~ClusterQt()
{

};

void ClusterQt::GetClusterSum()
{
  int qthtadr2 = qthtadr;
  //if(crate%2>0 && slot==3 && qthtadr>=25)
  //qthtadr2 = qthtadr-1;
  if(crate==3 && slot==2 && qthtadr==16) qthtadr2 = 24;
  if(crate==3 && slot==2 && qthtadr==24) qthtadr2 = 16;
      Nqt8=pQt->NTrigMap[qthtadr2][slot];
      if(Nqt8>0)
	{
	  Defined=true;
	  sum=0;
	  for(Int_t i=0;i<Nqt8;i++)
	    {
	      Cl_crate[i]=crate;
	      Cl_slot[i]=pQt->TrigMapslot[qthtadr2][slot][i];
	      Cl_qt8[i]=pQt->TrigMapqt8[qthtadr2][slot][i];
	      qt8sum[i]=pQt->qt8sum[Cl_qt8[i]][Cl_slot[i]][crate];
	      //cout<<"i="<<i<<" crate="<<crate<<" slot="<<Cl_slot[i]<<" qt8="<<Cl_qt8[i]<<" qt8sum="<<qt8sum[i]<<endl;
	      if(qt8trunc)//emulate L0 qt8 sum truncation
		{
		  qt8sum[i]=qt8sum[i]/32;//5 bit suppression 
		  if(qt8sum[i]>31){qt8sum[i]=31;};//5 bit maximum
		  //qt8sum[i]=qt8sum[i]*32;
		};
	      sum+=qt8sum[i];
	    };
	};

  //cout<<crate<<" "<<slot<<" "<<qthtadr2<<" "<<Nqt8<<" "<<sum<<endl;
};

Bool_t ClusterQt::AddQt8(Int_t crate_,Int_t slot_,Int_t qt8_)
{
  if(Defined)
    {
      if(Nqt8>=4){printf("%d QT8 boards already included. Up to 4 QT8s allowed in a cluster \n",Nqt8);return false;};
      if(crate_<0 || crate_>3){printf("Invalid crate id %d (0~3) \n",crate_);return false;};
      if(slot_<0 || slot_>9){printf("Invalid slot id %d (0~9) \n",slot_);return false;};
      if(qt8_<0 || qt8_>3){printf("Invalid QT8 id %d (0~3) \n",qt8_);return false;};
      Cl_crate[Nqt8]=crate_;
      Cl_slot[Nqt8]=slot_;
      Cl_qt8[Nqt8]=qt8_;
      //printf("Adding Crate %d %c%d to Cluster(Crate=%d Slot=%c HT=%d) \n",
      //     Cl_crate[Nqt8],Cl_slot[Nqt8]+65,Cl_qt8[Nqt8],crate,65+slot,qtht);
      qt8sum[Nqt8]=pQt->qt8sum[qt8_][slot_][crate_];
      if(qt8trunc)//emulate L0 qt8 sum truncation
	{
	  qt8sum[Nqt8]=qt8sum[Nqt8]/32;//5 bit suppression 
	  if(qt8sum[Nqt8]>31){qt8sum[Nqt8]=31;};//5 bit maximum
	  //qt8sum[Nqt8]=qt8sum[Nqt8]*32;
	};
      sum+=qt8sum[Nqt8];;
      Nqt8++;
      //cout<<"add crate="<<crate_<<" slot="<<slot_<<" qt8="<<qt8_<<" sum="<<sum<<" qt8sum="<<qt8sum[Nqt8-1]<<endl;
      return true;
    }
  else
    {
      printf("Cluster Not Defined \n");
      return false;
    };
};

void ClusterQt::Print()
{
  printf("\n");
  printf("QT 2009 Trigger Cluster\n");
  printf("Crate(0~3) = %d    Slot(A~J) = %c    ",crate,65+slot);
  if(Defined)
    {
      printf("HT Channel(0~31) = %d    HT ADC = %d\n",qthtadr,qtht);
      printf("FMS Module(NL,SL,NS,SS,1~4) = %d    FMS Channel(1~) = %d\n",detnstb,detchan);
      printf("List of QT8s Included in The Cluster \n");
      if(Nqt8==0){printf("NONE \n");};
      for(Int_t i=0;i<Nqt8;i++)
	{
	  printf("Crate = %d    Slot = %c    QT8 = %d    QT8 ADC Sum = %d\n",Cl_crate[i],Cl_slot[i]+65,Cl_qt8[i],qt8sum[i]);
	};
      printf("----------------------------------------------------------- \n");
      printf("Cluster ADC Sum = %d \n",sum);
      if(L0select){printf("Selected by L0 DSM \n");};
    }
  else{printf("No Cluster in this QT32 \n");};
  printf("\n");
};

void ClusterQt::Init()
{
  Defined=false;
  L0select=false;
  qt8trunc=false;
  sum=0;
  Nqt8=0;
  for(Int_t i=0;i<4;i++)
    {
      qt8sum[i]=0;
      Cl_crate[i]=-1;
      Cl_slot[i]=-1;
      Cl_qt8[i]=-1;
    };
};
