#include "TrigQt.h"
ClassImp(TrigQt)
using namespace std;

TrigQt::TrigQt(Qt* p_Qt)
{
  Init();
  pQt=p_Qt;
  TrigClusters=new TObjArray(100,0);
  // create a HT cluster for every Qt32 board
  for(Int_t crt=0;crt<4;crt++)
    {
      for(Int_t sl=0;sl<10;sl++)
	{
	  ClusterQt* Clust=new ClusterQt(pQt,crt,sl);
	  TrigClusters->AddAt(Clust,crt*10+sl);
	};
    };
};

TrigQt::TrigQt(Qt* p_Qt,Bool_t qt8trunc_)
{
  Init();
  pQt=p_Qt;
  TrigClusters=new TObjArray(100,0);
  // create a HT cluster for every Qt32 board 
  for(Int_t crt=0;crt<4;crt++)
    {
      for(Int_t sl=0;sl<10;sl++)
	{
	  ClusterQt* Clust=new ClusterQt(pQt,crt,sl,qt8trunc_);
	  TrigClusters->AddAt(Clust,crt*10+sl);
	};
    };
};

void TrigQt::Init()
{
  memset(L0_clsum,0,sizeof(L0_clsum));
  memset(L0_htadr,0,sizeof(L0_htadr));
  memset(L0_ht,0,sizeof(L0_ht));
  memset(L0_ht_thval,0,sizeof(L0_ht_thval));
  memset(L0_ht_thbit,0,sizeof(L0_ht_thbit));
  L2_TrigBits=0;
  memset(L1_TrigBits,0,sizeof(L1_TrigBits));

  //high tower threshold values
  L0_ht_thval[0][0]=20;
  L0_ht_thval[0][1]=10;
  L0_ht_thval[0][2]=10;
  L0_ht_thval[1][0]=20;
  L0_ht_thval[1][1]=10;
  L0_ht_thval[1][2]=10;
  L0_ht_thval[2][0]=20;
  L0_ht_thval[2][1]=10;
  L0_ht_thval[2][2]=10;
  L0_ht_thval[3][0]=20;
  L0_ht_thval[3][1]=10;
  L0_ht_thval[3][2]=10;
  //Int_t thS[3]={2,3,4};
  //Int_t thL[3]={2,3,4};
  //cluster threshold values for small and large dets
  thS[0]=20;
  thS[1]=31;
  thS[2]=31;
  thL[0]=10;
  thL[1]=21;
  thL[2]=17;
  memset(ntrig,0,sizeof(ntrig));
  clbitsSmall = 0;
  mclbitsSmall = 0;
  htbitsSmall = 0;
  clbitsLarge = 0;
  mclbitsLarge = 0;
  htbitsLarge = 0;
}

TrigQt::~TrigQt()
{
  TrigClusters->Delete();
};

void TrigQt::AddEdgeClusters()
{
  //Finish cluster sums in the boundary region
  TIter clusters(TrigClusters);
  while( ClusterQt* Cl = (ClusterQt*)clusters() )
    {
      if(Cl->Defined);
      {
	Int_t crt0=Cl->crate;
	Int_t slot0=Cl->slot;
	Int_t ch0=Cl->qthtadr;
	Int_t neighborAE=(crt0+2)%4;
	Int_t neighborDJ=crt0+(Int_t)pow(-1,crt0);
	if(slot0==0 && ch0>0  && ch0<6 )Cl->AddQt8(neighborAE,0,0); //A0
	if(crt0%2==0 && slot0==3 && ch0>24 && ch0<30)Cl->AddQt8(neighborDJ,3,3); //top D3
	if(crt0%2==1 && slot0==3 && ch0>25 && ch0<31)Cl->AddQt8(neighborDJ,3,3); //botom D3
	if(slot0==7 && ch0>24 && ch0<31)Cl->AddQt8(crt0,8,0); //H3
	if(slot0==8 && ch0>0  && ch0<7 )Cl->AddQt8(crt0,7,3); //I0
	if(slot0==9 && ch0>24 && ch0<31)Cl->AddQt8(neighborDJ,9,3); //J3	
      };
    };
};

void TrigQt::AddRealEdgeClusters()
{
  //Finish cluster sums in the boundary region
  TIter clusters(TrigClusters);
  while( ClusterQt* Cl = (ClusterQt*)clusters() )
    {
      if(Cl->Defined);
      {
	Int_t crt0=Cl->crate;
	Int_t slot0=Cl->slot;
	Int_t ch0=Cl->qthtadr;
	Int_t neighborAE=(crt0+2)%4;
	Int_t neighborDJ=crt0+(Int_t)pow(-1,crt0);
	if(slot0==0 && ch0>=0  && ch0<=6 )Cl->AddQt8(neighborAE,0,0); //A0
	if(crt0%2==0 && slot0==3 && ch0>=24 && ch0<=30)Cl->AddQt8(neighborDJ,3,3); //top D3
	if(crt0%2==1 && slot0==3 && ch0>=25 && ch0<=31)Cl->AddQt8(neighborDJ,3,3); //botom D3
	if(slot0==7 && ch0>=24 && ch0<=31)Cl->AddQt8(crt0,8,0); //H3
	if(slot0==8 && ch0>=0  && ch0<=7 )Cl->AddQt8(crt0,7,3); //I0
	if(slot0==9 && ch0>=24 && ch0<=31)Cl->AddQt8(neighborDJ,9,3); //J3	
      };
    };
};

void TrigQt::L0dsm()
{
  for(Int_t crt=0;crt<4;crt++)
    {
      Int_t HTstore[3][4];
      memset(HTstore,-1,sizeof(HTstore));
      for(Int_t sl=0;sl<10;sl++)
	{
	  ClusterQt* Cl=(ClusterQt*)TrigClusters->At(crt*10+sl);
	  Int_t L0DSM = sl/4;
	  Int_t L0Slot = sl%4;
	  //cout<<"crt="<<crt<<" slt="<<sl<<" qthtadr="<<Cl->qthtadr<<" qtht="<<Cl->qtht<<" sum="<<Cl->sum<<endl;
	  HTstore[L0DSM][L0Slot]=Cl->qtht;
	};
      Int_t pick=-1;
      Int_t npick=0;
      for(Int_t L0=0;L0<3;L0++)
	{
	  Int_t A=HTstore[L0][0];
	  Int_t B=HTstore[L0][1];
	  Int_t C=HTstore[L0][2];
	  Int_t D=HTstore[L0][3];

	  Int_t L1 = 0;
	  Int_t L0index = 0;
	  if(L0==0) {L1=0; L0index=crt;}
	  else if(crt==0) {L1=1; L0index=L0-1;}
	  else if(crt==1) {L1=1; L0index=L0+1;}
	  else if(crt==2) {L1=2; L0index=L0-1;}
	  else if(crt==3) {L1=2; L0index=L0+1;}

	  if(  A>B  && ((!(C>D) && !(D>A)) || ( C>D &&   A>C )) ){pick=0;npick++;};
	  if(!(A>B) && ((!(C>D) &&   B>D ) || ( C>D &&   B>C )) ){pick=1;npick++;};
	  if(  C>D  && ((!(A>B) && !(B>C)) || ( A>B && !(A>C))) ){pick=2;npick++;};
	  if(!(C>D) && ((!(A>B) && !(B>D)) || ( A>B &&   D>A )) ){pick=3;npick++;};

	  if(pick==-1){printf("Cluster HT ERROR!! %d %d %d %d \n",A,B,C,D);};

	  Int_t pickSl=L0*4+pick;
	  ClusterQt* Cl=(ClusterQt*)TrigClusters->At(crt*10+pickSl);
	  Int_t qha=Cl->qthtadr;//allows cluster selection when all HTs are zero
	  //if(!Cl->Defined){qha=0;};

	  L0_htadr[L0index][L1]=100*pick+qha;
	  L0_clsum[L0index][L1]=Cl->sum;
	  L0_ht[L0index][L1]=Cl->qtht;
	  Cl->L0select=true;

	  if(L0_ht[L0index][L1]>L0_ht_thval[L0index][L1]) L0_ht_thbit[L0index][L1]=true;
	  else L0_htadr[L0index][L1]=0;
	  //cout<<"crt="<<crt<<" L0="<<L0<<" L1="<<L1<<" L0index="<<L0index<<" htadr="
	  //<<L0_htadr[L0index][L1]<<" clsum="<<L0_clsum[L0index][L1]<<" qtht="<<L0_ht[L0index][L1]
	  //<<" L0Sel="<<L0_ht_thbit[L0index][L1]<<" A="<<A<<" B="<<B<<" C="<<C<<" D="<<D<<endl;
	};
    }


}


void TrigQt::L12dsm()
{
  for(int k=0; k<3; k++)
    {
      for(int j=0; j<4; j++)
	{
	  int htadr = L0_htadr[j][k];
	  int ht = L0_ht[j][k];
	  int thval = L0_ht_thval[j][k];
	  int crt0 = 0;
	  int slt0 = 0;
	  if(k==0)
	    {
	      crt0 = j;
	      slt0 = htadr/100;
	    }
	  else if(k==1) 
	    {
	      crt0 = j/2;
	      slt0 = 4+j%2*4+htadr/100;
	    }
	  else if(k==2)
	    {
	      crt0 = j/2+2;
	      slt0 = 4+j%2*4+htadr/100;
	    }
	  ClusterQt* Cl=(ClusterQt*)TrigClusters->At(crt0*10+slt0);

	  //int clsum = L0_clsum[j][k];
	  int clsum = Cl->sum;
	  //cout<<"crt="<<j<<" dsm0="<<k<<" k="<<k<<" clsum="<<clsum<<" clsum2="<<L0_clsum[j][k]<<" ht="<<ht<<endl;
	  if(k==0) //small cells
	    {
	      for(int kk=0; kk<3; kk++)
		{
		  if(clsum>thS[kk])
		    {
		      L1_TrigBits[k] = L1_TrigBits[k] | (1 << (j*3+kk)); //add cluster bits for small det
		      //cout<<"clsum="<<clsum<<" "<<thS[kk]<<" "<<k<<" "<<j<<" "<<kk<<" "<<hex<<L1_TrigBits[k]<<dec<<endl;
		    }
		}
	      if(ht>thval)
		{
		  L1_TrigBits[k] = L1_TrigBits[k] | (1 << (12+j)); //add ht bit for small det
		}
	    }
	  else //large cells
	    {
	      int tb = j/2; //top or bottom
	      for(int kk=0; kk<3; kk++)
		{
		  if(clsum>thL[kk])
		    {
		      L1_TrigBits[k] = L1_TrigBits[k] | (1 << (tb*3+kk)); //add cluster bits for large det
		    }
		}
	      if(ht>thval)
		{
		  L1_TrigBits[k] = L1_TrigBits[k] | (1 << (6+tb)); //add ht bit for large det
		}
	    }
	}
    }
  
  UInt_t bit[4]; //tmp
  //count multi trigger for small det
  for(int kk=0; kk<4; kk++) //four quardrants
    {
      bit[kk] = ( L1_TrigBits[0] >> (kk*3) ) & 0x7;
      clbitsSmall = clbitsSmall | bit[kk];
      if(bit[kk]&0x1) ntrig[0][0]++;
      if(bit[kk]&0x2) ntrig[0][1]++;
      if(bit[kk]&0x4) ntrig[0][2]++;
      //cout<<"test "<<kk<<" "<<clbitsSmall<<endl;
    }
  //set multi cluster bits for small det
  for(int kk=0; kk<3; kk++) //three threshold
    {
      if(ntrig[0][2]>1) mclbitsSmall = 1;
      else mclbitsSmall = 0;
      if(ntrig[0][1]>1) mclbitsSmall = (mclbitsSmall<<1) + 1;
      else mclbitsSmall = (mclbitsSmall<<1) + 0;
      if(ntrig[0][0]>1) mclbitsSmall = (mclbitsSmall<<1) + 1;
      else mclbitsSmall = (mclbitsSmall<<1) + 0;
    }
  //high tower bits for three ht threshold
  htbitsSmall  = (((L1_TrigBits[0]>>12)&0x1) | ((L1_TrigBits[0]>>13)&0x1) | ((L1_TrigBits[0]>>14)&0x1) | ((L1_TrigBits[0]>>15)&0x1));
  //same thing for large det
  bit[0] = (L1_TrigBits[1]&0x7);
  bit[1] = ((L1_TrigBits[1]>>3)&0x7);
  bit[2] = (L1_TrigBits[2]&0x7);
  bit[3] = ((L1_TrigBits[2]>>3)&0x7);
  for(int kk=0; kk<4; kk++)
    {
      clbitsLarge = clbitsLarge | bit[kk];
      if(bit[kk]&0x1) ntrig[1][0]++;
      if(bit[kk]&0x2) ntrig[1][1]++;
      if(bit[kk]&0x4) ntrig[1][2]++;
    }
  for(int kk=0; kk<3; kk++)
    {
      if(ntrig[1][2]>1) mclbitsLarge = 1;
      else mclbitsLarge = 0;
      if(ntrig[1][1]>1) mclbitsLarge = (mclbitsLarge<<1) + 1;
      else mclbitsLarge = (mclbitsLarge<<1) + 0;
      if(ntrig[1][0]>1) mclbitsLarge = (mclbitsLarge<<1) + 1;
      else mclbitsLarge = (mclbitsLarge<<1) + 0;
    }
  htbitsLarge  = ((L1_TrigBits[1]>>6)&0x1) | ((L1_TrigBits[1]>>7)&0x1) | ((L1_TrigBits[2]>>6)&0x1) | ((L1_TrigBits[2]>>7)&0x1);

  //L2 trigger bits
  L2_TrigBits = (htbitsLarge<<13) + (mclbitsLarge<<10) + (clbitsLarge<<7) + (htbitsSmall<<6) + (mclbitsSmall<<3) + clbitsSmall;
  
  //for(int j=0; j<3; j++)
  //cout<<"dsm1="<<j<<" L1TrigBits="<<hex<<L1_TrigBits[j]<<dec<<endl;
  //cout<<"L2_TrigBits="<<hex<<L2_TrigBits<<dec<<endl;
}

void TrigQt::printDSM()
{
  //printf("\ncrate   slot   qtdsm0ht   qtdsm0htadr   sum0   sum1   sum2   sum3   qtdsm0sum\n");
  printf("\ncrate   slot   qtdsm0ht   qtdsm0htadr   sum0   sum1   sum2   sum3\n");
  TIter clusters(TrigClusters);
  int ii = 0;
  while( ClusterQt* Cl = (ClusterQt*)clusters() )
    {
      if(Cl->Defined);
      {
	//printf("%5d   %4d  %8d   %11d   %4d   %4d   %4d   %4d   %9d\n",ii/10,ii%10,Cl->qtht,Cl->qthtadr,
	printf("%5d   %4d   %8d   %11d   %4d   %4d   %4d   %4d\n",ii/10,ii%10,Cl->qtht,Cl->qthtadr>0?Cl->qthtadr:0,
	       pQt->qt8sum[0][ii%10][ii/10]/32<=31?pQt->qt8sum[0][ii%10][ii/10]/32:31,
	       pQt->qt8sum[1][ii%10][ii/10]/32<=31?pQt->qt8sum[1][ii%10][ii/10]/32:31,
	       pQt->qt8sum[2][ii%10][ii/10]/32<=31?pQt->qt8sum[2][ii%10][ii/10]/32:31,
	       pQt->qt8sum[3][ii%10][ii/10]/32<=31?pQt->qt8sum[3][ii%10][ii/10]/32:31
	       //,Cl->sum
	       );
      }
      ii++;
    }
  printf("dsm1   dsm0   qt   qtdsm1clsum   qtdsm1htbit   qtdsm1htadr\n");
  for(Int_t i=0; i<3; i++) //L1 dsm
    {
      for(Int_t j=0; j<4; j++) //L0index
	{
	  if(L0_clsum[j][i])printf("%4d   %4d   %2d   %11d   %8d   %11d\n",i,j,L0_htadr[j][i]/100,L0_clsum[j][i],L0_ht_thbit[j][i],L0_htadr[j][i]%100<32?L0_htadr[j][i]%100:0);
	}
    }

  //cout<<"L1 trigger bits: "<<hex<<L1_TrigBits[0]<<" "<<L1_TrigBits[1]<<" "<<L1_TrigBits[2]<<dec<<endl;
  //cout<<"L2 trigger bits: "<<hex<<L2_TrigBits<<dec<<endl;
}
