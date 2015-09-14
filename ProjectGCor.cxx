#include "ProjectGCor.h"
ClassImp(ProjectGCor)
ProjectGCor:: ProjectGCor(TString oldgcorpath,TString RunDepPath,
			  TString newgcorpath)
{ 
  kPrint=false;
  char* path=(char*) RunDepPath.Data();
  Mgr=new RunDepMgr(path);
  Mgr->SetBase(16085005);
  Cal1=new CalibStr(16000000,oldgcorpath);
  Cal1Out=new CalibStr(16000000,oldgcorpath);
  Cal2=0;
  if(newgcorpath!="")  Cal2=new CalibStr(16000000,newgcorpath);
  for(int j=0;j<4;j++)
    {
      ratio1[j]=0;
      ratio2[j]=0;
      rat1[j]=0;
      rat2[j]=0;

    };
};
ProjectGCor::~ProjectGCor()
{
  for(int j=0;j<4;j++)
    {
      if(ratio1[j])delete ratio1[j];
      if(ratio2[j])delete ratio2[j];
      if(rat1[j])delete rat1[j];
      if(rat2[j])delete rat2[j];
    }
};

CalibStr* ProjectGCor::ProjectTo(Int_t RunNumber,Int_t OldRunNumber)
{
  Mgr->SetRdep(RunNumber);
  Mgr->SetBase(OldRunNumber);

  Int_t nrows[4]={34,34,24,24};
  Int_t ncols[4]={17,17,12,12};
  for(int j=0;j<4;j++)
    {
      if(ratio1[j])delete ratio1[j];
      if(ratio2[j])delete ratio2[j];
      if(rat1[j])delete rat1[j];
      if(rat2[j])delete rat2[j];
      char titl[100];
      sprintf(titl,"ratio12_%d",j);
      ratio1[j]=new TH2F(titl,titl,ncols[j],0,ncols[j],nrows[j],0,nrows[j]);
      sprintf(titl,"rat1_%d",j);
      rat1[j]=new TH1F(titl,titl,200,.25,2.25);
      sprintf(titl,"ratioproj_%d",j);
      ratio2[j]=new TH2F(titl,titl,ncols[j],0,ncols[j],nrows[j],0,nrows[j]);
      sprintf(titl,"rat2_p_%d",j);
      rat2[j]=new TH1F(titl,titl,200,.25,2.25);
    };
  for(int nstb=1;nstb<5;nstb++)
    {
      for(int row0=0;row0<nrows[nstb-1];row0++)
	{
	  for(int col0=0;col0<ncols[nstb-1];col0++)
	    {
	      if(Mgr->Rdep->Legal(2,nstb,row0,col0))
		{
		  int nbin=ratio1[nstb-1]->FindBin(col0+.5,row0+.5);
		  ratio1[nstb-1]->SetBinContent(nbin,1.);
		  ratio2[nstb-1]->SetBinContent(nbin,1.);
		  CellTDep* c;
		  if(c=Mgr->Rdep->Cdep(nstb,row0,col0))
		    {
		      TGraphErrors* gg;
		      if(gg=c->GetGraph())
			{
			  if(gg->GetFunction("pol0"))
			    {
			      if(gg->GetN()>10)
				{
				  Float_t nom=c->GetNominalLedADC();
				  Float_t thadc=c->GetthisLedADC();
				  Float_t oldval=Cal1->GetValue(2,nstb,row0,col0);
				  Float_t cal2val=Cal2->GetValue(2,nstb,row0,col0);
				  Float_t rc12=0;
				  Float_t rp12=0;
				  if(oldval!=0)rc12=cal2val/oldval;
				  // Float_t corr= Mgr->Rdep->GetCor(20,1,nstb,row0,col0,Cal1);
				  if(thadc>0)
				    {
				      Float_t corr=oldval*nom/thadc;
				      if(oldval!=0)rp12=corr/oldval;
				      Cal1Out->SetValue(2,nstb,row0,col0,corr);
				      nbin=ratio1[nstb-1]->FindBin(col0+.5,row0+.5);
				      ratio1[nstb-1]->SetBinContent(nbin,rc12);
				      ratio2[nstb-1]->SetBinContent(nbin,rp12);
				      rat1[nstb-1]->Fill(rc12);
				      rat2[nstb-1]->Fill(rp12);
				      if(kPrint)
					{
					  printf("nomled=%3.3f thisled=%3.3f Cal1=%3.3f Cal2=%3.3f projCal2=%3.3f  \n",nom,thadc,oldval,cal2val,corr);
					};
				    };
				};
			    }
			}
		    }
		}
	    }
	}
    };

  return Cal1Out;
};

