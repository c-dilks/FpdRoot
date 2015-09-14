#include "CellTDep.h"
ClassImp(CellTDep)
CellTDep::CellTDep(TGraphErrors* Gr,Int_t NSTB,Int_t Row0, Int_t Col0, CalibStr* GCor,Long_t Rnum,Float_t setNominalADC)
{
  /*
    setNominalADC to zero to cause this cell to set its current 
    adc average to the nominal adc
   */
  Error=0;
  NominalLedADC=setNominalADC;
  EvNumRange[0]=0;
  EvNumRange[1]=0;
  char name[200];
  sprintf(name,"ctdpr%d_c%d_%d_%ld",Row0,Col0,NSTB-1,Rnum);
  printf("name=%s\n",name);
  char gname[200];
  sprintf(gname,"gctdpr%d_c%d_%d_%ld",Row0,Col0,NSTB-1,Rnum);
  char fname[200];
  sprintf(fname,"fctdpr%d_c%d_%d_%ld",Row0,Col0,NSTB-1,Rnum);
  SetName(name);
  gr=new TGraphErrors();
  fit=new TF1();
  fit->SetName(fname);
  *gr=(TGraphErrors) *Gr;
  gr->SetName(gname);
  gr->SetTitle(gname);
  nstb=NSTB;
  row0=Row0;
  col0=Col0;
  N0=0;
  if(nstb<3)
    {
      chan=row0*17+col0+1;
    }
  else
    {
      chan=row0*12+col0+1;
    };
  NominalGcor=GCor->GetValue(2,nstb,row0,col0);
  if(gr)
    {
      int npnts=0;
      if((npnts=gr->GetN())>3)
	{
	  N0=npnts;
	  EvNumRange[0]=10000000;
	  EvNumRange[1]=0;
	  Sigma=0;
	  for(int j=0;j<npnts;j++)
	    {
	      Double_t x,y,ex,ey;
	      gr->GetPoint(j,x,y);
	      Sigma+=gr->GetErrorY(j);
	      if(x<EvNumRange[0])EvNumRange[0]=x;
	      if(x>EvNumRange[1])EvNumRange[1]=x;
	    };
	  Sigma=Sigma/npnts;
	  if(EvNumRange[0]<EvNumRange[1])
	    {
	      int fitReturn=0;
	      printf("Fit with %d points \n",npnts);
	      fitReturn=gr->Fit("pol0","","",EvNumRange[0],EvNumRange[1]);
	      printf("fitReturn=%d \n",fitReturn);
	      if(fitReturn==0)
		{
		  gr->Draw("A*");
		  TF1* fun=gr->GetFunction("pol0");
		  if(fun)
		    {
		      thisLedADC=fun->GetParameter(0);
		      thisfunError=fun->GetParError(0);
		      NDFpol0=fun->GetNDF();
		      chipdfpol0=0;
		      if(NDFpol0>0)chipdfpol0=fun->GetChisquare()/NDFpol0;
		      if(NominalLedADC==0)NominalLedADC=thisLedADC;
		      *fit=*fun;
		      fit->SetName(fname);
		    }
		  else
		    {
		      if(fit)delete fit;
		      fit=0;
		      printf("no pol0 function found\n");
		    };
		}
	      else
		{
		  Error=4;  //pol0 fit failed
		}
	    }
	  else
	    {
	      Error=3; //invalid eventnumber range
	    };
	}
      else
	{
	  Error=2;// to few hits to consider
	};
    }
  else
    {
      Error=1;// no graph
    };
  
};  
 
Float_t CellTDep::FixGcor(CalibStr* FpdCor,Int_t Segnum,Int_t Evtnum,Int_t mode)
{

  Float_t Gcor=1.;
  if(FpdCor==0)
    {
      Gcor=1;
    }
  else  
    {
      Gcor=FpdCor->GetValue(2,nstb,row0,col0);
    };

  if(mode==1 || (mode==2 &&chipdfpol0>0 && chipdfpol0<2))
    {
        Gcor=(Gcor*NominalLedADC)/thisLedADC;
      thisLedError=thisfunError;
    }
  else if(mode==2&&EvNumRange[1]>1)
    {
      Long_t evnumber=(Segnum-1)*10000+Evtnum;
      Float_t evperpoint=EvNumRange[1]/N0;
      Int_t pointnumber=(evnumber*N0)/(EvNumRange[1]);

      if(pointnumber<0)pointnumber=0;
      if(pointnumber>N0-1)pointnumber=N0-1;
      
      Double_t* xe=gr->GetX();
      Float_t xev=xe[pointnumber];
      int safecnt=0;
      
      while( (xev>evnumber+evperpoint) && pointnumber>0 && safecnt<10)
	{
	  safecnt++;
	  pointnumber=pointnumber-1;
	  xev=xe[pointnumber];
	  //	  printf(" xev=%f ",xev);

	} 
      while( (xev<evnumber-evperpoint) && pointnumber<(N0-1) && safecnt<10)
	{
	  safecnt++;
	  pointnumber=pointnumber+1;
	  xev=xe[pointnumber];
	  //	  printf("    xev=%f ",xev);
	} 
      Float_t ladc=gr->GetY()[pointnumber];
      //      printf(" pointnumber=%d ladc=%f \n",pointnumber,ladc);
      if(ladc>0)
	{
	  Gcor=(Gcor*NominalLedADC)/ladc;
	  thisLedError=Sigma;
	};

    };
  return Gcor;
};
CellTDep::~CellTDep()
{
  if(fit)delete fit;
  if(gr)delete gr;
};
