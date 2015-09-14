#include "RunDepMgr.h"
ClassImp(RunDepMgr);
ClassImp(DAtaEv);
DAtaEv::DAtaEv(int CurrentRun,int CurrentSeg, int CurrentEvent)
{
  CurrentRunNumber=CurrentRun;
  CurrentSegNumber=CurrentSeg;
  CurrentEventNumber=CurrentEvent;
};

DAtaEv::DAtaEv()
{
  CurrentRunNumber=0;
  CurrentSegNumber=0;
  CurrentEventNumber=0;
};

RunDepMgr::RunDepMgr(char* rundepPath,Int_t BaseOverride)
{
  da=new DAtaEv();
  RunDepPath=rundepPath;
  MissingRun=0;
  if(RunDepPath=="")RunDepPath=gSystem->Getenv("RunDepPath");
  printf("RunDepPath =%s\n",(const char*) RunDepPath);
  saveday=0;
  Base=0;
  Rdep=0;
  mvsday=new TObjArray(1300);
  mvsday->SetOwner();
  RunDepBaseOverride=BaseOverride;
  if(BaseOverride>0)SetBase(BaseOverride);
  SetUnitLedFactors(); 
}
void RunDepMgr::SetUnitLedFactors()
{
  rows[0]=rows[1]=34;  rows[2]=rows[3]=24;
  for(int det=0;det<4;det++)
    {
      for(int ro=0;ro<34;ro++)
	{
	  for(int co=0;co<17;co++)
	    {
	      ledFactor[det][ro][co]=0;
	      if(Legal(2,det+1,ro,co))
		{
		  ledFactor[det][ro][co]=1;
		}
	    };
	}
    };
  
};
RunDepMgr::RunDepMgr()
{
  Rdep=new RunDepCor();
  Base=new RunDepCor();
  mvsday=new TObjArray();
  da=new DAtaEv();
}
RunDepMgr::~RunDepMgr()
{
  if(Base)delete Base;
  if(Rdep)delete Rdep;
  if(mvsday)delete mvsday;
  if(da)delete da;
}
RunDepCor* RunDepMgr::SetRdep(Int_t Rnum)
{
  /*
    replace Rdep with one corresponding to Run Rnum
    -if Rnum mismatch, delete old Rdep and creat new one from RunDep.root
    -Use predefined base id predefined
    -Setup new Rdep to reflect this base and save if changed
  */
  
  bool WriteRdep=false;
  if(Rdep)
    {
      if(Rdep->RunNumber==Rnum)return Rdep;
      delete Rdep;
    }
  
  printf("opening RunDepPath=%s\n",(const char*) RunDepPath);

  TFile* SaveDefaultFile=gROOT->GetFile();
  TFile* RdepFile=new TFile(RunDepPath);
  if(RdepFile) 
    {
      char sflb[200];
      sprintf(sflb,"RDC%d",Rnum);
      Rdep=(RunDepCor*) RdepFile->Get(sflb);
      if(!Rdep) 
	{
	  printf("Rdep for %s not in database\n",sflb);
	  delete RdepFile;
	  return 0;
	};
      delete RdepFile;
    }
  else
    {
      return 0;
    };
  if(RunDepBaseOverride>12000000)
    {
      SetBase(RunDepBaseOverride);
    };    
  if(SaveDefaultFile)SaveDefaultFile->cd();
  TIter next(Rdep->array);
  CellTDep* cp=0;
  while(cp=(CellTDep*) next())
    {
      Int_t cpnstb=cp->GetNstb();
      Int_t cprow0=cp->GetRow0();
      Int_t cpcol0=cp->GetCol0();
      if(Legal(2,cpnstb,cprow0,cpcol0))
	{
	  CellTDep* cpb=0;
	  //	  printf("nstb=%d row0=%d col0=%d \n",cpnstb,cprow0,cpcol0);
	  cpb=Base->Cdep(cpnstb,cprow0,cpcol0);
	  if(cpb)
	    { 
	      Float_t newbase=cpb->GetthisLedADC();
	      cp->SetNominalLedADC(newbase);
	    };
	  
	};
    };
  if(WriteRdep)
    {
      TFile* SaveDefaultFile=gROOT->GetFile();

      RdepFile=new TFile(RunDepPath,"Update");
      if(RdepFile)
	{
	  Rdep->Write();
	  printf("%s replaced in RunDep.root file\n",(const char*) Rdep->GetName());
	  RdepFile->Purge();
	  delete RdepFile;
	}
      if(SaveDefaultFile)SaveDefaultFile->cd();

    };
  UpdateledFactors();  
  return Rdep;
}
void RunDepMgr::UpdateledFactors()
{
  SetUnitLedFactors();
  if(!mvsday)return;
  if(!Rdep)return;
  if(!Base)return;
  printf("Updating Led Factors for %d base:%d \n",(int) Rdep->RunNumber,(int) Base->RunNumber);

  if(Day(Rdep->RunNumber)==saveday)return;
  Long_t runbase=Base->RunNumber;
  Long_t runR=Rdep->RunNumber;
  Int_t daybase=Day(Base->RunNumber);
  Int_t dayR=Day(Rdep->RunNumber);
  TIter next(mvsday);
  //  printf("Update Led Factors now\n");
  TGraphErrors* gr;
  TGraphErrors* nomgr=0;
  Bool_t ok;
  while(gr=(TGraphErrors*) next())
    {
      const char* snm=(const char*)  gr->GetName();
      int ro,co,de;
      sscanf(snm,"mday_r%d_c%d_%d",&ro,&co,&de);
      if(Legal(2,de+1,ro,co))
	{
	  TF1* f=gr->GetFunction("pol1");
	  if(f)
	    {
	      Float_t MdayR=f->Eval(dayR);
	      Float_t MdayB=f->Eval(daybase);
	      ledFactor[de][ro][co]=MdayB/MdayR;
	      CellTDep* rct=Rdep->Cdep(de+1,ro,co);
	      CellTDep* bct=Base->Cdep(de+1,ro,co);
	      ok=false;
	      nomgr=0;
	      if(rct!=0 && bct!=0)
		{
		  Float_t nomba=bct->GetNominalLedADC();
		  Float_t nomrd=rct->GetNominalLedADC();
		  nomgr=bct->GetGraph();
		  if(nomgr)
		    {
		      TF1* nomf1=nomgr->GetFunction("pol0");
		      if(nomf1)
			{
			  Float_t nompar0=nomf1->GetParameter(0);
			  Float_t factor=Rdep->Cdep(de+1,ro,co)->FixFactor(5000,2);
			  //			  printf("de=%d ro=%d co=%d ledfactor=%f factor=%f mpar0=%f NominalLedADC=%f thisadc=%f\n",de,ro,co,ledFactor[de][ro][co],factor,nompar0,rct->GetNominalLedADC(),rct->GetthisLedADC());
			  ok=true;
			};
		    };
		};
	    }
	  if(!ok)
	    {
	      printf("Not ok : de=%d ro=%d co=%d ",de,ro,co);
	      if(nomgr)printf(" npts=%d ",nomgr->GetN());
	      printf("\n");
	    };
	};
    };
  printf("Exit Update Led\n");
};

RunDepCor* RunDepMgr::SetBase(Int_t BaseOverride)
{
  RunDepBaseOverride=BaseOverride;
  if(!Rdep)
    {
      printf("Error: can't call SetBase if Rdep not defined \n"); 
      return 0;
    };
  if(Base)
    {
      if(BaseOverride!=Base->RunNumber) 
	{
	  delete Base;
	}
      else
	{
	  Rdep->BaseRunNumber=BaseOverride;
	  return Base;
	};
    };
  Rdep->BaseRunNumber=BaseOverride;  
  
  printf("form SetBase: opening RunDepPath=%s\n",(const char*) RunDepPath);
  TFile* SaveDefaultFile=gROOT->GetFile();
  TFile* RdepFile=new TFile(RunDepPath);
  if(RdepFile)
    {
      char sflbase[200];
      sprintf(sflbase,"RDC%d",BaseOverride);
      Base=(RunDepCor*) RdepFile->Get(sflbase);
      if(!Base)
	{
	  printf("could not find base: %s  in file (use self)\n",sflbase);
	  Base=Rdep;
	  Rdep->BaseRunNumber=Rdep->RunNumber;
	}
      delete RdepFile;
      if(SaveDefaultFile)SaveDefaultFile->cd();
    };
 
  UpdateledFactors();
  UpdateRdep();
  printf("return SetBase\n");
  return Base;
}

void RunDepMgr::Rcorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW,DAtaEv* dda)
  {
    da->CurrentRunNumber=dda->CurrentRunNumber;
    da->CurrentSegNumber=dda->CurrentSegNumber;
    da->CurrentEventNumber=dda->CurrentEventNumber;
    Rcorrect(padc,pEmat,det,iEW);
    return;
  };

void RunDepMgr::Rcorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW,dataSet* dda)
{
  da->CurrentRunNumber=dda->CurrentRunNumber;
  da->CurrentSegNumber=dda->CurrentSegNumber;
  da->CurrentEventNumber=dda->CurrentEventNumber;
  Rcorrect(padc,pEmat,det,iEW);
  return;
};
void RunDepMgr::Rcorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW)
{
  if(da==0)return;
  if(da->CurrentRunNumber<12000000)return;
  if(da->CurrentRunNumber==MissingRun)return;
  MissingRun=0;

  if(Rdep==0)SetRdep(da->CurrentRunNumber);
  if(da->CurrentRunNumber!=Rdep->RunNumber)SetRdep(da->CurrentRunNumber);

  if(Rdep==0) {return;};
  int nr=padc->GetNrows();
  int nc=padc->GetNcols();
  Bool_t goodpar=true;
  if(iEW!=2)goodpar=false; 
  if(det<0 || det>3)goodpar=false;
  if(det<2)
    {
      if(nr!=34)goodpar=false;
      if(nc!=17)goodpar=false;
    }
  else
    {
      if(nr!=24)goodpar=false;
      if(nc!=12)goodpar=false;
    };
  if(goodpar==false)
    {
      printf("bad parameters\n");
      return;
    }
  CellTDep* ct;
  for(int ro=0;ro<nr;ro++)
    {
      for(int co=0;co<nc;co++)
	{
	  if(((*padc)(ro,co))>1)
	    {
	      if(Legal(2,det+1,ro,co))
		{
		  ct=Rdep->Cdep(det+1,ro,co);
		  /*		  printf("ct->thisLedADC= %f   ct->NominalLedADC= %f \n",
			 ct->GetthisLedADC(),ct->GetNominalLedADC());
		  */
		  
		  Float_t ffac=ct->FixFactor(da->CurrentSegNumber,
					     da->CurrentEventNumber,2);
		  Float_t ledfac=ledFactor[det][ro][co];
		  (*pEmat)(ro,co)=((*pEmat)(ro,co))*ffac*ledfac;
		  /*		  printf("det=%d ro=%d co=%d ffac=%f ledfac=%f en=%f\n",det,ro,co,ffac,ledfac,(*pEmat)(ro,co));
		   */
		};
	    }
	};
    };
  //  printf("\nCorrected Esum=%f \n",pEmat->Sum());
}
Bool_t RunDepMgr::Legal(Int_t iew,Int_t nstb,Int_t row0,Int_t col0)
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
void RunDepMgr::Addmvsday(Int_t nstb,Int_t row,Int_t col,TGraphErrors* gr)
{
  if(gr)
    {
      TGraphErrors* ngf=new TGraphErrors(*gr);
      if(ngf)
	{
	  char  gnam[300];
	  sprintf(gnam,"mday_r%d_c%d_%d",row,col,nstb-1);
	  ngf->SetName(gnam);
	  mvsday->Add(ngf);   
	}
      else
	{
	  printf("Failed to make graph\n");
	};
    }
}
void RunDepMgr::UpdateRdep()
{
  
  if(!Rdep)return;
  if(!(Rdep->array))return;
  if(!Base)return;
  Rdep->BaseRunNumber=Base->RunNumber;
  TIter next(Rdep->array);
  CellTDep* ct;
  printf("UpdateRdep for %d with base %d \n",(int) Rdep->RunNumber,(int) Base->RunNumber);
  while(ct=(CellTDep*) next())
    {
      CellTDep* ctb=Base->Cdep(ct->GetNstb(),ct->GetRow0(),ct->GetCol0());
      if(ctb)
	{
	  TGraphErrors* gr;
	  if(gr=ctb->GetGraph())
	    {
	      TF1* f=gr->GetFunction("pol0");
	      if(f)
		{
		  Float_t par0=f->GetParameter(0);
		  if(par0>1 && par0<4096)ct->SetNominalLedADC(par0);
		};
	    }
	}
    };
};

