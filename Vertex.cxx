#include "Vertex.h"

ClassImp(Vertex)

Vertex::Vertex()
{
  InitializeBbc();
};

Vertex::~Vertex()
{
};

void Vertex::InitializeBbc()
{
  Bbc_threshold=10;
  Int_t east_q_map_tmp[24]={8, 5, 4, 40, 37, 36, 7, 6, 
			    3, 2, 1, 39, 38, 35, 34, 33, 
			    72, 71, 70, 69, 68, 67, 66, 65}; 
  Int_t east_t_map_tmp[16]={16, 13, 12, 48, 45, 44, 15, 14, 
			    11, 10, 9, 47, 46, 43, 42, 41}; 
  Int_t west_q_map_tmp[24]={24, 21, 20, 56, 53, 52, 23, 22, 
			    19, 18, 17, 55, 54, 51, 50, 49, 
			    80, 79, 78, 77, 76, 75, 74, 73}; 
  Int_t west_t_map_tmp[16]={32, 29, 28, 64, 61, 60, 31, 30, 
			    27, 26, 25, 63, 62, 59, 58, 57}; 
  Int_t offsetpp_east[16]={-999, -34, -31, -31, -37, -29, -19, -37, 
			  -37, -22, -999, -34, -35, -37, -35, -999};
  Int_t offsetpp_west[16]={-27, -25, -28, -30, -37, -30, -23, -28, 
			   -32, -32, -31, -24, -30, -999, -26, -22};

  for(Int_t i=0;i<24;i++)
    {
      east_q_map[i]=east_q_map_tmp[i]-1;
      west_q_map[i]=west_q_map_tmp[i]-1;
      if(i<16)
	{
	  east_t_map[i]=east_t_map_tmp[i]-1;
	  west_t_map[i]=west_t_map_tmp[i]-1;
	  offset[i][0]=offsetpp_east[i];
	  offset[i][1]=offsetpp_west[i];
	};
    };
};

Float_t Vertex::GetBbcVertex(UChar_t* p_Bbc)
{   
  Int_t chg[24][2];
  Int_t tac[16][2];
  maxqe=0;
  maxqw=0;
  ieqmax=0;
  iwqmax=0;
    
  iemax=0;
  iwmax=0;
  maxtace=0;
  maxtacw=0;
  qe=0;
  qw=0;

  zvt=999.;
  ewdiff=-9999;
  ewsum=-9999;

  for(Int_t i=0;i<24;i++)
    {
      if(east_q_map[i]>=0)
	{
	  chg[i][0]=p_Bbc[east_q_map[i]];
	}
      else
	{
	  chg[i][0]=0;
	};
      
      if(chg[i][0]>ieqmax)
	{
	  maxqe=chg[i][0];
	  ieqmax=i;
	};

      if(west_q_map[i]>=0)
	{
	  chg[i][1]=p_Bbc[west_q_map[i]];
	}
      else
	{
	  chg[i][1]=0;
	};
      
      if(chg[i][1]>iwqmax)
	{
	  maxqw=chg[i][1];
	  iwqmax=i;
	};
 
      if(i<16)
	{
	  if(east_t_map[i]>=0)
	    {
	      tac[i][0]=p_Bbc[east_t_map[i]]-offset[i][0];
	      if(tac[i][0]>maxtace && tac[i][0]<=200 && chg[i][0]>Bbc_threshold)
		{
                  maxtace=tac[i][0];
                  iemax=i;
		  qe=chg[i][0];
		};
	    }
	  else
	    {
	      tac[i][0]=0;
	    };
	  if(west_t_map[i]>=0)
	    {
	      tac[i][1]=p_Bbc[west_t_map[i]]-offset[i][1];
	      if(tac[i][1]>maxtacw && tac[i][1]<=200 && chg[i][1]>Bbc_threshold)
		{
                  maxtacw=tac[i][1];
                  iwmax=i;
		  qw=chg[i][1];
		};
	    }
	  else
	    {
	      tac[i][1]=0;
	    };
	};
    };

if(maxtace>0 && maxtacw>0)
  {
    ewdiff=maxtacw-maxtace;
    ewsum=maxtacw+maxtace;
    Float_t scale=1.7;
    zvt=scale*ewdiff;
  };
return zvt;
};
