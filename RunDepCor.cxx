#include "RunDepCor.h"

ClassImp(RunDepCor)
  RunDepCor::RunDepCor(Long_t Rnum,TString adcTrLpath,CalibStr* initgCor,RunDepCor* base)
{
  RunNumber=Rnum;
  if(base)
    {
      BaseRunNumber=base->RunNumber;
    }
  else
    {
      BaseRunNumber=RunNumber;
    };
  char rnm[200];
  sprintf(rnm,"RDC%ld",Rnum);
  fName=rnm;
  TFile* f=FillLedADC(adcTrLpath);
  array=new TObjArray(1300,0);
  MakeDepCor(f,initgCor,base);
  if(f)delete f;
};

void RunDepCor::MakeDepCor(TFile* ledfile,CalibStr* initgCor,RunDepCor* base)
{
  Int_t nrws[4]={34,34,24,24};
  for(int j=0;j<4;j++)nrows[j]=nrws[j];
  CalibStr* Initgcor=initgCor;
  CellTDep* ct;
  for(int d=0;d<4;d++)
    {
      for(int row=0;row<nrows[d];row++)
	{
	  for(int col=0;col<nrows[d]/2;col++)
	    {
	      Location[d][row][col]=-1;
	      if(Legal(2,d+1,row,col))
		{
		  char name[200];
		  sprintf(name,"rTadcr%d_c%d_%d",row,col,d);
		  TGraphErrors* ggr=(TGraphErrors*) ledfile->Get(name);
		  CellTDep* bct=0;
		  Float_t LedNom=0.;
		  if(base)
		    {
		      bct=base->Cdep(d+1,row,col);
		      if(bct)LedNom=bct->GetNominalLedADC();
		    };
		  array->Add(ct=new CellTDep(
					     ggr,d+1,row,col,
					     Initgcor,RunNumber,LedNom));
		  dname[d][row][col]=ct->GetName();
		  Location[d][row][col]=array->GetLast();
		}
	    };
	};
    };
};

TFile* RunDepCor::FillLedADC(TString adcTrLpath)
{
  TFile *f = new TFile(adcTrLpath,"read"); 
  if(f==0)return 0;
  TTree *Tr_adc = (TTree*)f->Get("Tr_adc");
  const char* tmpD=gSystem->Getenv("RDCtmpDir");
  if(tmpD=="")tmpD="./";
  TString TmpL=tmpD;
  TmpL=TmpL+"/tmpL.root";
  TFile* tmpL=new TFile(TmpL,"recreate");
  //Declaration of leaves types
  Int_t           ADC;
  Int_t           rowADC;
  Int_t           colADC;
  Int_t           nstbADC;
  Int_t           runnum;
  Int_t           evtnum;
  Int_t           segnum;
  Int_t           led;
  Int_t           Bunchid7bit;
  
  // Set branch addresses.
  Tr_adc->SetBranchAddress("br_ADC",&ADC);
  Tr_adc->SetBranchAddress("br_rowADC",&rowADC);
  Tr_adc->SetBranchAddress("br_colADC",&colADC);
  Tr_adc->SetBranchAddress("br_nstbADC",&nstbADC);
  Tr_adc->SetBranchAddress("br_runnum",&runnum);
  Tr_adc->SetBranchAddress("br_evtnum",&evtnum);
  Tr_adc->SetBranchAddress("br_segnum",&segnum);
  Tr_adc->SetBranchAddress("br_led",&led);
  Tr_adc->SetBranchAddress("br_Bunchid7bit",&Bunchid7bit);
  
  Long64_t nentries = Tr_adc->GetEntries();
  
  TGraph* tdep[4][34][17]; 
  TGraphErrors* tdepr[4][34][17];
  TString Tname[4][34][17];
  TString Tnamer[4][34][17];
  
  Int_t nr[4]={34,34,24,24};
  for(int d=0;d<4;d++)
    {
      for(int row=0;row<nr[d];row++)
	{
	  for(int col=0;col<nr[d]/2;col++)
	    {
	      char str[100];
	      sprintf(str,"Tadcr%d_c%d_%d",row,col,d);
	      tdep[d][row][col]=new TGraph();
	      tdep[d][row][col]->SetName(str);
	      tdep[d][row][col]->SetMarkerSize(.2);
	      tdep[d][row][col]->SetMarkerColor(3);	       
	      Tname[d][row][col]=str;
	      sprintf(str,"rTadcr%d_c%d_%d",row,col,d);
	      tdepr[d][row][col]=new TGraphErrors(0);
	      tdepr[d][row][col]->SetName(str);
	      tdepr[d][row][col]->SetMarkerSize(.2);	       
	      tdepr[d][row][col]->SetMarkerColor(2);	       
	      Tnamer[d][row][col]=str;
	    };
	};
    };
  Int_t cnt=0;
  if(nentries>100000000)nentries=100000000;
  for (Long64_t i=0; i<nentries;i++) 
    {
     Tr_adc->GetEntry(i);
      if(Legal(2,nstbADC,rowADC,colADC))
	{
	  
	  if(cnt++%100000==0)printf("cnt=%d adc hits\n",cnt);
	  if(led==1)
	    {
	      TGraph* pg=tdep[nstbADC-1][rowADC][colADC];
	      int np=pg->GetN();
	      int n0=np-40;
	      if(n0<0)n0=0;
	      Double_t x,y;
	      y=ADC;
	      if(y<=0)y=.1;
	      x=((segnum-1)*10000+evtnum);
	      pg->SetPoint(np,x,y);	     
	    };
	};
    };
  for(int d=0;d<4;d++)
    {
      for(int row=0;row<nr[d];row++)
	{
	  for(int col=0;col<nr[d]/2;col++)
	    {
	      if(Legal(2,d+1,row,col))
		{
		  TGraph* pg=tdep[d][row][col];
		  pg->Sort();
		  TGraphErrors* pgr=tdepr[d][row][col];
		  int N0= pg->GetN();
		  int NN0=pgr->GetN();
		  Double_t x,y;
		  pg->GetPoint(0,x,y);
		  Double_t xmin=x;
		  pg->GetPoint(N0-1,x,y);
		  Double_t xmax=x;
		  if(N0>40)
		    {
		      
		      for(int nn=20;nn<N0;nn=nn+40)
			{
			  Double_t x0,y0;
			  pg->GetPoint(nn,x0,y0);
			  Double_t x10=x0;
			  Double_t y10=0;
			  Double_t sig10=0;
			  Int_t count=0;
			  for(int np=nn-20;np<nn+20;np++)
			    {
			      count++;
			      pg->GetPoint(np,x0,y0);
			      y10+=y0;
			      sig10+=y0*y0;
			    };
			  if(count==0)count==1;
			  y10=y10/count;
			  sig10=sqrt((sig10-count*y10*y10))/count;
			  NN0=pgr->GetN();
			  pgr->SetPoint(NN0,x10,y10);
			  pgr->SetPointError(NN0,0,sig10);
			}
		    };
		  Float_t chi=0;
		  Float_t dof=0;
		  Float_t chipd=chi;
		  tdepr[d][row][col]->Write(Tnamer[d][row][col]);
		  
		  if(tdepr[d][row][col]->GetN()>3)
		    {
		      if(tdepr[d][row][col]->Fit("pol0")==0)
			{		
			  printf("fit 1 good\n");
			  chi=tdepr[d][row][col]->GetFunction("pol0")->GetChisquare();
			  dof=tdepr[d][row][col]->GetFunction("pol0")->GetNDF();
			  chipd=chi;
			};
		    };
		  if(dof>10)
		    {
		      chipd=chi/dof;
		      
		      if(chipd>2)
			{
			  printf("try pol3\n");
			  Int_t retval=tdepr[d][row][col]->Fit("pol3");
			  printf("Fit Retval=%d \n",retval);
			  if(retval==0)
			    {
			      chi=tdepr[d][row][col]->GetFunction("pol3")->GetChisquare();
			      dof=tdepr[d][row][col]->GetFunction("pol3")->GetNDF();
			      chipd=chi;
			      if(dof>2)chipd=chi/dof;
			    }
			  else
			    {
			      printf("failed pol3\n");
			    };
			}
		    };
		  tdepr[d][row][col]->Write(Tnamer[d][row][col]);
		};
	      if(tdepr[d][row][col])delete tdepr[d][row][col];
	      if(tdep[d][row][col])delete tdep[d][row][col];
	    };
	};
      
    };
  if(Tr_adc)delete Tr_adc;
  printf("tree Tr_adc deleted\n");
  if(f) delete f;
  printf("file f deleted\n");

  return tmpL;
  
};
RunDepCor::~RunDepCor()
{
  //  if(Initgcor)delete Initgcor;
  array->Delete();
  if(array)delete array;
};
RunDepCor::RunDepCor()
{
  array=new TObjArray();
  array->SetOwner();
};
Bool_t RunDepCor::Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0)
{
  if(iew>0 && iew<2)return false;
  if(nstb<1 || nstb>4)return false;
  if(nstb>2)
    {
      if(row0<0 || row0>23)return false;
      if(col0<0 || col0>11)return false;
      if(fabs(1.*row0-11.5)<5 && col0<5)return false;
    }
  else
    {
      if(row0<0 || row0>33)return false;
      if(col0<0 || col0>16)return false;
      if(fabs(1.*row0-16.5)<8 && col0<8)return false;
      if(row0<col0-9.5)return false;
      if(33-row0<col0-9.5)return false;
    };
  return true;
}
Float_t RunDepCor::GetCor(Int_t Segment_number, Int_t evNumber, Int_t NSTB,
			  Int_t Row, Int_t Col, CalibStr* GCor,Int_t mode)
{
  if(Legal(2,NSTB,Row,Col)&&dname[NSTB-1][Row][Col]!="")
    {
      CellTDep* ct=Cdep(NSTB,Row,Col);
      CalibStr* cal=0;
      if(GCor)
	{cal=GCor;}
      else
	{
	  //cal=Initgcor;
	  cal=0;
	};
      if(ct!=0 && cal!=0)
	{
	  return ct->FixGcor(cal,Segment_number,evNumber,mode);
	}
      else
	{
	  return 1.;
	};
    }
  else
    {
      return 1.;
    };
}

CellTDep*  RunDepCor::Cdep(Int_t nstb,Int_t row0,Int_t col0)
{

  CellTDep* ct=0;
  if(!Legal(2,nstb,row0,col0))return ct;
  if(Location[nstb-1][row0][col0]>=0)
    {
      ct=(CellTDep*) array->At(Location[nstb-1][row0][col0]);
    }	   
  else
    {
      printf("getting CellTDep by name??\n");
      ct=(CellTDep*) array->FindObject(dname[nstb-1][row0][col0]);
      if(ct)Location[nstb-1][row0][col0]=array->IndexOf(ct);
    }
  return ct;
}
