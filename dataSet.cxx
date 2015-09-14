#include <stdio.h>
#include "Fpdchan.h"
#include "CalibStr.h"
#include "dataSet.h"
#include "Constants.h"
ClassImp(dataSet)
ClassImp(FpdMap)

dataSet::dataSet(char* RunString,FilesSet* pfset,char* ntupleName)
{
  pFiles=pfset;
  NtpName=ntupleName;
  dataSetBuild((char*) RunString,(char*) pFiles->p_fpdped()->Path(),
				 (char*) pFiles->p_fpdgain()->Path(),
				 (char*) pFiles->p_fpdgaincorr()->Path());
};

dataSet::dataSet(char*  RunString, 
		 char*  fpdpedfile,
		 char*  fpdgainfile,
		 char*  fpdgaincorrfile)
{
  NtpName="h111";
  pFiles=0;
  dataSetBuild(RunString,fpdpedfile,fpdgainfile,fpdgaincorrfile);
};

dataSet::dataSet(char*  RunString, 
		 char*  fpdpedfile,
		 char*  fpdgainfile,
		 char*  fpdgaincorrfile,char* ntupleName)
{
  NtpName=ntupleName;
  pFiles=0;
  dataSetBuild(RunString,fpdpedfile,fpdgainfile,fpdgaincorrfile);
};

Int_t dataSet::dataSetBuild(char*  RunString, 
		  char*  fpdpedfile,
		  char*  fpdgainfile,
		  char*  fpdgaincorrfile)
{
  if(pFiles==0)
    {

      pFiles=new FilesSet("./",fpdpedfile,fpdgainfile,fpdgaincorrfile,"fill.txt","rlumi.txt","spinpat","geom.txt");
    };
  error=0;
  LastError="";
  CurrentRunNumber=0;
  currentTreeNumber=-2;
  RunPed=0;
  R1Lum=0.;
  R2Lum=0.;
  pGeom=new Geom(pFiles);
  if(NtpName=="h112")CurrentRunNumber=7094066;
  if(NtpName=="h111_08")CurrentRunNumber=9067059;
  if(NtpName=="h111_09")CurrentRunNumber=10063033;

  ErrR1Lum=0.;
  ErrR2Lum=0.;
  BlueSpin=0;
  YellowSpin=0;
  kicked=true;
  // Setup chain
  std::cout<<"enter with: "<<NtpName<<"\n";

  if(NtpName=="h111")
    {
      printf("opening h111\n");
      Input=new TChain("h111");
    }
  else if(NtpName=="h112")
    {
      Input=new TChain("h112");
      printf("opening h112\n");
    }
  else if(NtpName=="h111_08")
    {
      Input=new TChain("h111");
      printf("opening h111_08\n");
    }
  else if(NtpName=="h111_09")
    {
      Input=new TChain("h111");
      printf("opening h111_09\n");
    }
  else
    {
      error=-1000;
      LastError="Error: Neither h111 or h112 selected";
      return error;
    };
  Input->Add(RunString);
  std::cout<<"Input from "<<RunString<<"\n";
  Input->ls();
  //  TTree* tree=(TTree*)gDirectory->Get("h111");
  Init(Input);
  FileName=((TChain*) fChain)->GetListOfFiles()->At(0)->GetTitle();
  const char* fnam=(const char*) FileName;
  Int_t loc=FileName.Index("run")+3;
  Int_t runnum;
  sscanf(&(fnam[loc]),"%d",&runnum);
  CurrentRunNumber=runnum;
  // Read In calibration
  Rped=new CalibStr(CurrentRunNumber,fpdpedfile);
  std::cout<<"Rped->filename: "<<Rped->filename<<"\n";
  Rgain=new CalibStr(CurrentRunNumber,fpdgainfile);
  Rgaincorr=new CalibStr(CurrentRunNumber,fpdgaincorrfile);
  UseQt=false;
  if(pFiles->p_qtmap())
    {
      if(pFiles->p_qtmap()->FileName!="none" && CurrentRunNumber>9000000)
	{
	  printf("----------Use QT---------- \n");
	  UseQt=true;
	  pQt=new Qt(pFiles);
	};
    };
  if(Rped->error!=0)error=error+1;
  if(error!=0)LastError="Ped error";
  if(Rgain->error!=0)error=error+2;
  if(Rgaincorr->error!=0)error=error+4;
  if(error>1)LastError="Gain Error";

  RFill=0;
  RunsArray=0;
  Int_t entlast=Input->GetEntries();
  currentTreeNumber=-1;
  GetEntry(entlast);
  Int_t lastrun=CurrentRunNumber;
  GetEntry(1);
  if(error==-4000)error=0;
  Int_t firstrun=CurrentRunNumber;
  printf("Calling Fill with run range %d to %d\n",firstrun,lastrun);
  RFill=new Fill(pFiles,10000,24000,firstrun,lastrun);
  currentTreeNumber=-1;
  if(RFill->error>0)
    {
      error=error+8;
      printf("%s \n",(const char*) RFill->ErrorMessage);
      LastError="Rfill message:"+RFill->ErrorMessage;
    };
  if(error==0 )
    {
      printf("RunSArray Created\n");
      RunsArray=new TObjArray(1000,0);
      FILE* lum;

	if(pFiles->p_fpdlumi()->open())
	  {
	    lum=pFiles->p_fpdlumi()->fp;
	    char ln[100];
	    Int_t rcnt=0;
	    while(!feof(lum))
	    {
	      if(fgets(ln,90,lum))
		{
		  Int_t runn,v1,v2;
		  RunData* rdat;
		  Float_t rluma,erluma,rlumb,erlumb;
		  sscanf(ln,"%d %d %d %f %f %f %f",&runn,&v1,&v2,&rluma,&erluma,&rlumb,&erlumb);
		  //printf("Runn=%d \n",runn);
		  if( (runn>=RFill->FirstRun) && (runn<=RFill->LastRun) )
		    {
		      if(RFill->SetFillNumberforRun(runn)>0)
			{
			  printf("Adding %d  to RunsArray\n",runn);
			  RunsArray->AddAtAndExpand(rdat=new RunData(runn,RFill,pFiles),rcnt);
			  rdat->R1Lum=rluma;
			  rdat->R2Lum=rlumb;
			  rdat->ErrR1Lum=erluma;
			  rdat->ErrR2Lum=erlumb;
			  rcnt++;
			};
		    };
		};
	    };
	    
	    pFiles->p_fpdlumi()->close();
	  }
	else
	  {
	    printf("%s file not found\n",
		   (const char*) pFiles->p_fpdlumi()->Path());
	    error=error+64;
	  };
    };
};

TMatrix dataSet::OutProd(TVector& v1,TVector& v2)
{
  Int_t nr,nc;
  TMatrix outp(nr=v1.GetNrows(),nc=v2.GetNrows());
  for(Int_t i=0;i<nr;i++)
    {
      for(Int_t j=0;j<nc;j++)
	{
	  outp(i,j)=v1(i)*v2(j);
	};
    };
  return outp;
};

TMatrix dataSet::ZeroNegativePedGain2(TMatrix* t0,TMatrix* t1,TMatrix* t2,TMatrix* t3)
{
  Int_t nrows0=t0->GetNrows();
  Int_t ncols0=t0->GetNcols();
  Int_t nrows1=t1->GetNrows();
  Int_t ncols1=t1->GetNcols();
  Int_t nrows2=t2->GetNrows();
  Int_t ncols2=t2->GetNcols();
  Int_t nrows3=t3->GetNrows();
  Int_t ncols3=t3->GetNcols();
  if(nrows0!=nrows1)error=-100;
  if(nrows0!=nrows2)error=-100;
  if(nrows0!=nrows3)error=-100;
  if(error<0)return TMatrix(1,1);
  *t0=(*t0-(*t1))*(*t2)*(*t3);

  return *t0;

};

TMatrix dataSet::ZeroNegativePed2(TMatrix* t0,TMatrix* t1,TMatrix* t2)
{
  Int_t nrows0=t0->GetNrows();
  Int_t ncols0=t0->GetNcols();
  Int_t nrows1=t1->GetNrows();
  Int_t ncols1=t1->GetNcols();
  Int_t nrows2=t2->GetNrows();
  Int_t ncols2=t2->GetNcols();

  if(nrows0!=nrows1)error=-100;
  if(nrows0!=nrows2)error=-100;

  *t0=(*t0-(*t1)-(*t2));
  return *t0;
};

TMatrix dataSet::ZeroNegativePed2Gain2(TMatrix* t0,TMatrix* t1,TMatrix* t2,TMatrix* t3,TMatrix* t4)
{
  Int_t nrows0=t0->GetNrows();
  Int_t ncols0=t0->GetNcols();
  Int_t nrows1=t1->GetNrows();
  Int_t ncols1=t1->GetNcols();
  Int_t nrows2=t2->GetNrows();
  Int_t ncols2=t2->GetNcols();
  Int_t nrows3=t3->GetNrows();
  Int_t ncols3=t3->GetNcols();
  Int_t nrows4=t4->GetNrows();
  Int_t ncols4=t4->GetNcols();

  if(nrows0!=nrows1)error=-100;
  if(nrows0!=nrows2)error=-100;
  if(nrows0!=nrows3)error=-100;
  if(nrows0!=nrows4)error=-100;
  if(error<0)return TMatrix(1,1);
  *t0=(*t0-(*t1)-(*t2))*(*t3)*(*t4);
  return *t0;

};

TMatrix dataSet::Em(TMatrix* t0,Int_t ew,Int_t nstb)
{
  Int_t nrows=t0->GetNrows();
  Int_t ncols=t0->GetNcols();

  Int_t nelements= nrows*ncols;
  error=0;
  TMatrix nmat(nrows,ncols);
  nmat=0;
  TMatrix*  RP1=Rped->tm(ew,nstb);
  TMatrix*  RG1=Rgain->tm(ew,nstb);
  TMatrix*  RG2=Rgaincorr->tm(ew,nstb);
  nmat=*t0-*RP1;
  if(RunPed!=0)
    {
      TMatrix tmp(nrows,ncols);
      tmp=*(RunPed->tm(ew,nstb));

      error=-2000;
      if(tmp.GetNrows()==nrows)
	{
	  nmat=nmat-tmp;
	  error=0;
	}
      else
	{
	  error=-2000;
	};
    };
  if(RP1->GetNoElements()!=nelements && RG1->GetNoElements()!=nelements && RG2->GetNoElements()!=nelements )
    {
      error=-3000;
    };

  if(error<0)
    {
      TMatrix nmat(2,2);nmat=-1000;}
  else
    {
      for(int ir=0;ir<nrows;ir++)
	{
	  for(int ic=0;ic<ncols;ic++)
	    {
	      nmat(ir,ic)=nmat(ir,ic)*(*RG1)(ir,ic)*(*RG2)(ir,ic);
	    };
	};
    };
  return nmat;
};


TMatrix dataSet::dMatrix(TMatrix* map, UChar_t* FpdNtpBlock)
{
  Int_t nrows=map->GetNrows();
  Int_t ncols=map->GetNcols();
  TMatrix out(nrows,ncols);
  Int_t loc=0;
  for(Int_t i=0;i<ncols;i++)
    {
      for(Int_t j=0;j<nrows;j++)
	{
	  loc=(Int_t) (*map)(j,i);
	  if(loc>=0)
	    {
	      out(j,i)=FpdNtpBlock[loc];
	    }
	  else
	    {
	      out(j,i)=0.;
	    };
	};
    };
  return out;
}
TMatrix dataSet::dMatrixQt(Int_t i)
{
  return *(pQt->tm(i));
};
TMatrix dataSet::dMatrix(TMatrix* map, UInt_t* FpdNtpBlock)
{
  Int_t nrows=map->GetNrows();
  Int_t ncols=map->GetNcols();
  TMatrix out(nrows,ncols);
  Int_t loc=0;
  for(Int_t i=0;i<ncols;i++)
    {
      for(Int_t j=0;j<nrows;j++)
	{
	  loc=(Int_t) (*map)(j,i);
	  if(loc>=0)
	    {
	      out(j,i)=FpdNtpBlock[loc];
	    }
	  else
	    {
	      out(j,i)=0.;
	    };
	};
    };
  return out;
};

Bool_t dataSet::decodeQT()
{
  Bool_t decode=pQt->decodeQT(nqtdata,Qtdata,0);
  if(!decode)return false;
  return true;
};

Bool_t dataSet::decodeDSM()
{
  Bool_t decode=pQt->decodeDSM(Fpdw,Fpddsm,lastdsm[6]);
  if(!decode)return false;
  return true;
};

dataSet::~dataSet()
{
  if(Rped!=0)delete Rped;
  if(Rgain!=0)delete Rgain;
  if(Rgaincorr!=0)delete Rgaincorr;
  if(RFill!=0)delete RFill;
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t dataSet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   Int_t entr=fChain->GetEntry(entry);
   CurrentEventNumber++;
   Int_t treenum=((TChain*) fChain)->GetTreeNumber();
   if(treenum!=currentTreeNumber)
     {
       CurrentSegNumber=1;
       currentTreeNumber=treenum;
       FileName=((TChain*) fChain)->GetListOfFiles()->At(treenum)->GetTitle(); 
       const char* fnam=(const char*) FileName;
       Int_t loc=FileName.Index("run")+3;
       Int_t runnum,segnum;
       sscanf(&(fnam[loc]),"%d.%d",&runnum,&segnum);
       CurrentRunNumber=runnum;
       CurrentSegNumber=segnum;
       CurrentEventNumber=0;
       thisRunDat=0;
       if(RFill!=0)
	 {
	   RFill->pfd=0;
	   RFill->CurrentRunNumber=0;
	   RFill->SetFillNumberforRun(CurrentRunNumber);
	 };
       if(RunsArray!=0)
	 {
	   thisRunDat=0;
	   TIter next(RunsArray);
	   RunData* tmpRD;
	   while(tmpRD=(RunData*) next())
	     {
	       if(tmpRD->RunNumber==runnum)
		 {
		   thisRunDat=tmpRD;
		   RunPed=thisRunDat->RunPed;
		 };
	     };
	 };
     };
   if(RFill!=0)
     {
       if(RFill->GetFillNumber()<0)RFill->SetFillNumberforRun(CurrentRunNumber);
     };
   if(thisRunDat!=0) 
     {
       R1Lum= thisRunDat->R1Lum;
       R2Lum= thisRunDat->R2Lum;
       ErrR1Lum= thisRunDat->ErrR1Lum;
       ErrR2Lum= thisRunDat->ErrR2Lum;
       RunPed=thisRunDat->RunPed;
     };
   if(RFill!=0)
     {
       BlueSpin=RFill->BlueSpin(bunchid7bit);
       YellowSpin=RFill->YellowSpin(bunchid7bit);
       kicked=RFill->pfd->kicked(bunchid7bit);
     };
   if(thisRunDat==0)error=-4000;
   return entr;
}

Long64_t dataSet::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->IsA() != TChain::Class()) return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void dataSet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Token", &Token, &b_Token);
   fChain->SetBranchAddress("addbit", &addbit, &b_addbit);
   fChain->SetBranchAddress("Dsm", &Dsm, &b_Dsm);
   fChain->SetBranchAddress("Trgwd", Trgwd, &b_Trgwd);
   fChain->SetBranchAddress("lastdsm", lastdsm, &b_lastdsm);
   fChain->SetBranchAddress("Bcd", Bcd, &b_Bcd);
   fChain->SetBranchAddress("L1", L1, &b_L1);
   fChain->SetBranchAddress("Vtxdsm", Vtxdsm, &b_Vtxdsm);
   fChain->SetBranchAddress("Fpddsm", Fpddsm, &b_Fpddsm);
   fChain->SetBranchAddress("Emcdsm", Emcdsm, &b_Emcdsm);
   fChain->SetBranchAddress("Bbcl1", Bbcl1, &b_Bbcl1);
   fChain->SetBranchAddress("bunchid", &bunchid, &b_bunchid);
   fChain->SetBranchAddress("bunchid7bit", &bunchid7bit, &b_bunchid7bit);

   if(NtpName=="h111" || NtpName=="h112")
     {
       fChain->SetBranchAddress("event", &event, &b_event);
       fChain->SetBranchAddress("ipre", &ipre, &b_ipre);
       fChain->SetBranchAddress("ipost", &ipost, &b_ipost);
       fChain->SetBranchAddress("Bchi", &Bchi, &b_Bchi);
       fChain->SetBranchAddress("Bclo", &Bclo, &b_Bclo);
       fChain->SetBranchAddress("npre", &npre, &b_npre);
       fChain->SetBranchAddress("npost", &npost, &b_npost);
       fChain->SetBranchAddress("Cpa", Cpa, &b_Cpa);
       fChain->SetBranchAddress("Zdc", Zdc, &b_Zdc);
       fChain->SetBranchAddress("Zdcsmd", Zdcsmd, &b_Zdcsmd);
       fChain->SetBranchAddress("Ctb", Ctb, &b_Ctb);
       fChain->SetBranchAddress("Bbc", Bbc, &b_Bbc);
       fChain->SetBranchAddress("Fpdens", Fpdens, &b_Fpdens);
       fChain->SetBranchAddress("Fpdetb", Fpdetb, &b_Fpdetb);
       fChain->SetBranchAddress("Fpdensl1", Fpdensl1, &b_Fpdensl1);
       fChain->SetBranchAddress("Awdbm", &Awdbm, &b_Awdbm);
       fChain->SetBranchAddress("Awtrg", &Awtrg, &b_Awtrg);
       fChain->SetBranchAddress("Awdaq", &Awdaq, &b_Awdaq);
       fChain->SetBranchAddress("extbusy", &extbusy, &b_extbusy);
       fChain->SetBranchAddress("modbusy", &modbusy, &b_modbusy);
       fChain->SetBranchAddress("physwd", &physwd, &b_physwd);
       fChain->SetBranchAddress("Dsmad", &Dsmad, &b_Dsmad);
       fChain->SetBranchAddress("cotamibusy", &cotamibusy, &b_cotamibusy);
       fChain->SetBranchAddress("qdsm", qdsm, &b_qdsm);
       fChain->SetBranchAddress("Mwc", Mwc, &b_Mwc);
       fChain->SetBranchAddress("Bemc", Bemc, &b_Bemc);
       fChain->SetBranchAddress("Fpdwns", Fpdwns, &b_Fpdwns);
       fChain->SetBranchAddress("Fpdwtb", Fpdwtb, &b_Fpdwtb);
       fChain->SetBranchAddress("Fpdwnsl1", Fpdwnsl1, &b_Fpdwnsl1);
       fChain->SetBranchAddress("Fpdwtbl1", Fpdwtbl1, &b_Fpdwtbl1);
       fChain->SetBranchAddress("Fpdetbl1", Fpdetbl1, &b_Fpdetbl1);

     };

   if(NtpName=="h112")
     {
       fChain->SetBranchAddress("event", &event, &b_event);
       fChain->SetBranchAddress("ipre", &ipre, &b_ipre);
       fChain->SetBranchAddress("ipost", &ipost, &b_ipost);
       fChain->SetBranchAddress("Bchi", &Bchi, &b_Bchi);
       fChain->SetBranchAddress("Bclo", &Bclo, &b_Bclo);
       fChain->SetBranchAddress("npre", &npre, &b_npre);
       fChain->SetBranchAddress("npost", &npost, &b_npost);
       fChain->SetBranchAddress("Cpa", Cpa, &b_Cpa);
       fChain->SetBranchAddress("Zdc", Zdc, &b_Zdc);
       fChain->SetBranchAddress("Zdcsmd", Zdcsmd, &b_Zdcsmd);
       fChain->SetBranchAddress("Ctb", Ctb, &b_Ctb);
       fChain->SetBranchAddress("Bbc", Bbc, &b_Bbc);
       fChain->SetBranchAddress("Fpdens", Fpdens, &b_Fpdens);
       fChain->SetBranchAddress("Fpdetb", Fpdetb, &b_Fpdetb);
       fChain->SetBranchAddress("Fpdensl1", Fpdensl1, &b_Fpdensl1);
       fChain->SetBranchAddress("num10", &num10, &b_num10);
       fChain->SetBranchAddress("num50", &num50, &b_num50);
       fChain->SetBranchAddress("rega", &rega, &b_rega);
       fChain->SetBranchAddress("cwnentry", &cwnentry, &b_cwnentry);
       fChain->SetBranchAddress("fpdwadc", fpdwadc, &b_fpdwadc);
       fChain->SetBranchAddress("tdc", tdc, &b_tdc);
     };   

   if(NtpName=="h111_08")
     {
       fChain->SetBranchAddress("event", &event, &b_event);
       fChain->SetBranchAddress("ipre", &ipre, &b_ipre);
       fChain->SetBranchAddress("ipost", &ipost, &b_ipost);
       fChain->SetBranchAddress("Bchi", &Bchi, &b_Bchi);
       fChain->SetBranchAddress("Bclo", &Bclo, &b_Bclo);
       fChain->SetBranchAddress("npre", &npre, &b_npre);
       fChain->SetBranchAddress("npost", &npost, &b_npost);
       fChain->SetBranchAddress("Cpa", Cpa, &b_Cpa);
       fChain->SetBranchAddress("Zdc", Zdc, &b_Zdc);
       fChain->SetBranchAddress("Zdcsmd", Zdcsmd, &b_Zdcsmd);
       fChain->SetBranchAddress("Ctb", Ctb, &b_Ctb);
       fChain->SetBranchAddress("Bbc", Bbc, &b_Bbc);
       fChain->SetBranchAddress("Fpdens", Fpdens, &b_Fpdens);
       fChain->SetBranchAddress("Fpdetb", Fpdetb, &b_Fpdetb);
       fChain->SetBranchAddress("Fpdensl1", Fpdensl1, &b_Fpdensl1);
       fChain->SetBranchAddress("Awd", Awd, &b_Awd);
       fChain->SetBranchAddress("busystatus", busystatus, &b_busystatus);
       fChain->SetBranchAddress("Mtddsm", Mtddsm, &b_Mtddsm);
       fChain->SetBranchAddress("Vpddsm", Vpddsm, &b_Vpddsm);
       fChain->SetBranchAddress("Ctbdsm", Ctbdsm, &b_Ctbdsm);
       fChain->SetBranchAddress("specialtriggers", specialtriggers, &b_specialtriggers);
       fChain->SetBranchAddress("L2", L2, &b_L2);
       fChain->SetBranchAddress("Mtd", Mtd, &b_Mtd);
       fChain->SetBranchAddress("Vpd", Vpd, &b_Vpd);
       fChain->SetBranchAddress("Tof", Tof, &b_Tof);
       fChain->SetBranchAddress("Bemceast", Bemceast, &b_Bemceast);
       fChain->SetBranchAddress("Bemcwest", Bemcwest, &b_Bemcwest);
       fChain->SetBranchAddress("Bemclayer1", Bemclayer1, &b_Bemclayer1);
       fChain->SetBranchAddress("Eemclayer1", Eemclayer1, &b_Eemclayer1);
       fChain->SetBranchAddress("Eemc", Eemc, &b_Eemc);
       fChain->SetBranchAddress("Fpdw", Fpdw, &b_Fpdw);
       fChain->SetBranchAddress("Zdcl1", Zdcl1, &b_Zdcl1);
       fChain->SetBranchAddress("nqtdata", &nqtdata, &b_nqtdata);
       fChain->SetBranchAddress("Qtdata", Qtdata, &b_Qtdata);
       fChain->SetBranchAddress("nemcdata", &nemcdata, &b_nemcdata);
       fChain->SetBranchAddress("Emcdata", Emcdata, &b_Emcdata);
     };

   if(NtpName=="h111_09")
     {
       fChain->SetBranchAddress("evt", &evt, &b_evt); 
       fChain->SetBranchAddress("iprepost", &iprepost, &b_iprepost); 
       fChain->SetBranchAddress("Bc", Bc, &b_Bc);
       fChain->SetBranchAddress("nprepost", nprepost, &b_nprepost); 
       fChain->SetBranchAddress("Tofdsm", Tofdsm, &b_Tofdsm); 
       fChain->SetBranchAddress("Vpd", Vpd_09, &b_Vpd);
       fChain->SetBranchAddress("Mtd", Mtd_09, &b_Mtd);
       fChain->SetBranchAddress("Tof", Tof_09, &b_Tof);
       fChain->SetBranchAddress("Tofl1", Tofl1, &b_Tofl1);
       fChain->SetBranchAddress("Fpde", Fpde, &b_Fpde);
       fChain->SetBranchAddress("busystatus", busystatus, &b_busystatus);
       //fChain->SetBranchAddress("Mtddsm", Mtddsm, &b_Mtddsm);
       //fChain->SetBranchAddress("Vpddsm", Vpddsm, &b_Vpddsm);
       fChain->SetBranchAddress("specialtriggers", specialtriggers, &b_specialtriggers);
       fChain->SetBranchAddress("L2", L2, &b_L2);
       fChain->SetBranchAddress("Bemceast", Bemceast, &b_Bemceast);
       fChain->SetBranchAddress("Bemcwest", Bemcwest, &b_Bemcwest);
       fChain->SetBranchAddress("Bemclayer1", Bemclayer1, &b_Bemclayer1);
       fChain->SetBranchAddress("Eemclayer1", Eemclayer1, &b_Eemclayer1);
       fChain->SetBranchAddress("Eemc", Eemc, &b_Eemc);
       fChain->SetBranchAddress("Fpdw", Fpdw, &b_Fpdw);
       fChain->SetBranchAddress("Zdcl1", Zdcl1, &b_Zdcl1);
       fChain->SetBranchAddress("nqtdata", &nqtdata, &b_nqtdata);
       fChain->SetBranchAddress("Qtdata", Qtdata, &b_Qtdata);
       fChain->SetBranchAddress("nemcdata", &nemcdata, &b_nemcdata);
       fChain->SetBranchAddress("Emcdata", Emcdata, &b_Emcdata);
       fChain->SetBranchAddress("L2sum", L2sum, &b_L2sum);
     };

   Notify();
}

Bool_t dataSet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dataSet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t dataSet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
