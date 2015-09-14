#include "Sheet.h"
using namespace std;
#include <iostream>
#include <string>
#include <sstream>

ClassImp(Sheet);
ClassImp(SheetHistory);

Sheet::Sheet(char* cvsHead,char* cvsLine)
{
  Parsecol=0;
  Resets=new TObjArray();
  Extras=new TObjArray();
  Extras->IsOwner();
  Run9Gain=-.00001;
  gr=0;
  fun=0;
  HVfactor=1.;
  TString sh=cvsHead;
  HeadLine=sh;
  TObjArray* sho=sh.Tokenize(",");

  TString sl=cvsLine;
  TheLine=sl;
  TObjArray* slo=sl.Tokenize(",");
  Int_t Nslo=slo->GetEntries();
  Int_t Nsho=sho->GetEntries();
  if(Nslo>18)Nslo=18;
  if(Nsho>18)Nsho=18;
  NumberElem=Nslo;
  if(Nsho!=Nslo)
    {
      printf("header entry mismatch ");
      if(Nsho>Nslo)Nsho=Nslo;
      if(Nsho<Nslo)Nslo=Nsho;
    };
  Int_t val[20];
  Int_t max=Nsho;

  for(int k=0;k<max;k++)
    {
      TObjString* ss=(TObjString*) slo->At(k);
      TObjString* hs=(TObjString*) sho->At(k);
      TString sss=ss->GetString();
      TString hss=hs->GetString();
      if(k<18)Vname[k]=hss;

      const char* cl=(const char*) sss;

      Int_t newvar=0;
      if(k==17 || k==16)
	{
	  if(val[0]<3)
	    {
	      sscanf(cl,"%d",&newvar);
	    }
	  else
	    {
	      sscanf(cl,"%X",&newvar);
	    };
	}
      else
	{
	      sscanf(cl,"%d",&newvar);
	};

      val[k]=newvar;      
    };
  char SheetName[200];
  if(max>4)
    {
      nstb1=val[0];
      chan1=val[1];
      myrow0=val[2];
      mycol0=val[3];
      PatchPanel=val[6];
      PatchRow=val[7];
      PatchCol=val[8];
      QTCrate=val[9];
      QTSlot=val[10];
      QTCard4=val[11];
      QTchannel8=val[12];
      id1=val[13];
      id2=val[14];
      id3=val[15];
      id4=val[16];
      HVsetValue=val[17];

      sprintf(SheetName,"Shr%d_c%d_%d",myrow0,mycol0,nstb1-1);
      Name=SheetName;
    };
  delete sho;
  delete slo;
  
};  
TString Sheet::MakecsvHead()
{
  TString hd="Nstb,Chan1,Row0,Col0,,,PPanel,PRow,PCol,QTCrate,QTSlot,QTCard4,QTchan8,id1,id2,id3,id4,OldHV,OldGcorr,NewHV,NewGcorr,Gain,OldShift,NewShift";
  return hd;
};
TString Sheet::MakecsvLine()
{
  std::stringstream sss;
  
  sss<<nstb1<<","<<
    chan1<<","<<
    myrow0<<","<<
    mycol0<<",,,"<<
    PatchPanel<<","<<
    PatchRow<<","<<
    PatchCol<<","<<
    QTCrate<<","<<
    QTSlot<<","<<
    QTCard4<<","<<
    QTchannel8<<","<<
    id1<<","<<
    id2<<","<<
    id3<<","<<
    id4<<",";
  SheetHistory* hh;
  if(hh=GetHistory())
    {
      sss<<hh->OldHV <<","<<
	hh->OldGcorr <<","<<
	hh->NewHV <<","<<
	hh->NewGcorr <<","<<
	hh->Gain <<","<<
	hh->OldShift <<","<<
	hh->NewShift;
    };
  sss<<"\n";
  return TString(sss.str());

}
void Sheet::Print()
{
  
  std::cout<<Name<<"\n";
  std::cout<<Vname[0]<<"="<<nstb1<<"\n";
  std::cout<<Vname[2]<<"="<<myrow0<<"\n";
  std::cout<<Vname[3]<<"="<<mycol0<<"\n";

  std::cout<<Vname[9]<<"="<<QTCrate<<"\n";
  std::cout<<Vname[10]<<"="<<QTSlot<<"\n";
  std::cout<<Vname[11]<<"="<<QTCard4<<"\n";
  std::cout<<Vname[12]<<"="<<QTchannel8<<"\n";


  std::cout<<Vname[13]<<"="<<id1<<"\n";
  std::cout<<Vname[14]<<"="<<id2<<"\n";
  std::cout<<Vname[15]<<"="<<id3<<"\n";
  std::cout<<Vname[16]<<"="<<id4<<"\n";
  std::cout<<Vname[17]<<"="<<HVsetValue<<"\n";
  
};
TString Sheet::SetCommand(Int_t bitmask,Bool_t nw)
{
  /*
    make HV command string
    small cells:
    bitmask&1 for device command
    bitmask&2 for setctrl command
    bitmask&4 for rdac command
   */
  TString quote="\"";
  const char* qt=(const char* ) quote;
  if(nstb1<1 || nstb1>4)return "";
  if(nstb1==3 || nstb1==4)
    {
      char s1[200];
      char s2[200];
      char s3[200];
      sprintf(s1,"!SETdevice %d \n",id1-1);
      sprintf(s2,"!setctrl %d %d \n",id2,id3);

      Int_t HV=(Int_t) HVsetValue*HVfactor;
      if(!nw)HV=HVsetValue;
      if(HV>255)HV=255;
      if(HV<0)HV=0;
       sprintf(s3,"!rdac %2X %2X\n!sleep 100\n",id4&0xFF,HV&0xFF);
      TString rval="";
      if(bitmask&1)rval=rval+s1;
      if(bitmask&2)rval=rval+s2;
      if(bitmask&4)rval=rval+s3;
      return rval;
    }
  else
    {
      char s[200];
      Int_t HV=(Int_t) HVsetValue*HVfactor;
      if(!nw)HV=HVsetValue;
      sprintf(s,"# port %d r%dc%d NSTB=%d CHAN=%d\n;sleep 1;echo -e %swrite (%2d,%2d) %d%s\n",id1,myrow0+1,mycol0+1,nstb1,chan1,qt,id3,id4,HV,qt);
      return s;
    }; 
  return "";
}
				     
void Sheet::Parse(TString Str,Int_t Iport,Int_t NSTB1)
{
  if(!Resets)Resets=new TObjArray();
  if(NSTB1!=nstb1)return;
 const char* s=(const char*) Str;
  int off=0;
  if(NSTB1<3 && Iport==id1)
    {
      while(strlen(s)>2)
	{
	  off=TString(s).First("i");
	  if(off<0)off=strlen(s)-2;
	  s+=off;
	  if(s[1]=='t' && s[2]=='e')
	    {
	      s+=3;
	      int r=-1;
	      int c=-1;
	      int v=-1;
	      if(sscanf(s,"(%d,%d) %d ",&r,&c,&v)==3)
		{
		  char rcv[200];
		  sprintf(rcv,"NSTB1_%d_Port%d_r%d_c%d_val%d",NSTB1,Iport,r,c,v);
		  Resets->Add(new TObjString(TString(rcv)));		  
		  if(r==id3 && c==id4)
		    {
		      std::cout<<GetName()<<"HV change "<<HVsetValue<<"->"<<v<<"\n";
		      HVsetValue=v;
		      
		    };
		};
	    }
	}
    }
  else
    {
      int chip=-1;
      int chan=-1;
      int addr=-1;
      while(strlen(s)>2)
	{
	  
	  off=TString(s).First("!");
	  int chp,chn,ad,v;
	  if(sscanf(s,"setctrl %d %d",&chp,&chn)==2)
	    {
	      chip=chp;
	      chan=chn;
	    }
	  else if(sscanf(s,"rdac %x %x",&ad,&v)==2)
	    {
	      std::cout<<GetName()<<"HV change "<<HVsetValue<<"->"<<v<<"\n";
	      if(id2==chp && id3==chn && id4==ad)HVsetValue=v;
	      char rcv[200];
	      sprintf(rcv,"NSTB1_%d_chp%d_chn%2x_ad%2x_val%d",NSTB1,chp,chn,ad,v);
	      Resets->Add(new TObjString(TString(rcv)));		  
	    }
	};
    };
};
TString Sheet::DoString(char* fname)
{
  FILE* fp;
  fp=fopen(fname,"r");
  char ss[200];
  TString mys="";
  while(!feof(fp))
    {
      char* ch=fgets(ss,199,fp);
      mys=mys+TString(ss);
    };
  return mys;
};
void Sheet::AddGr(TFile* HvGr)
{
  printf("Enter AddGr\n");
  if(HvGr==0)return;
  char grname[100];
  char funname[100];
  sprintf(grname,"gr_r%d_c%d_%d",myrow0,mycol0,nstb1-1);
  sprintf(funname,"fun_r%d_c%d_%d",myrow0,mycol0,nstb1-1);
  TKey* key=HvGr->FindKey(grname);
  gr=0;
  if(key)
    {
      printf("key non-zero\n");
      TGraphErrors* tmpgr=(TGraphErrors*) key->ReadObj();
      if(tmpgr)gr=new TGraphErrors(*tmpgr);
    };
  if(gr)
    {
      printf("gr non-zero\n");
      TString gstr="";
      gstr=gstr+grname;
      gr->SetName(gstr);
      TString fstr="";
      TF1* tmpfun=0;
      TKey* key2=HvGr->FindKey(funname);
      if(key2)
	{
	  if(tmpfun=(TF1*) HvGr->FindKey(funname)->ReadObj())
	    {
	      printf("tmpfun non-zero\n");
	      fun=new TF1(*tmpfun); 
	      fstr=fstr+funname;
	      fun->SetName(fstr);
	    } else fun=0;
	}
      else fun=0;
    }
  else
    {
      fun=0;
    };
  printf("Exit AddGr \n ");
};
Int_t Sheet::DoReset(TObjArray* resets)
{
  if(resets)
    {
      TIter next(resets);
      TObjString* str;
      char rname[200];
      if(nstb1<3)
	{
	  sprintf(rname,"NSTB1_%d_Port%d_r%d_c%d_val",nstb1,id1,id3,id4);
	}
      else
	{
	  sprintf(rname,"NSTB1_%d_chp%d_chn%2x_ad%2x_val",nstb1,id2,id3,id4);
	}
      TString rstring=rname;
      while(str=(TObjString*) next())
	{
	  TString st;
	  st=str->GetString();
	  if(st.Index(rstring)==0)
	    {
	      Int_t len=rstring.Length();
	      char* cst=(char*) &st;
	      Int_t val;
	      sscanf(&(cst[len]),"%d",&val);
	      HVsetValue=val;
	    }
	}
    };
  return HVsetValue;
};
Sheet::Sheet()
{
  gr=new TGraphErrors();
  fun=new TF1();
  Resets=new TObjArray();
  Extras=new TObjArray();
};
TF1* Sheet::GetGrFun()
{
  if(!gr)return 0;
  return gr->GetFunction("fun");
};
Int_t Sheet::HvScale(Float_t gcorr,Int_t shift)
{
  /*
    This routin returns the high voltage that 
  1) reduces gcorr to unity
  2) compensats for a change of shift by "shift bits" 
   */
  Int_t retval=-1;
  Float_t savesign=1.;

  if(gr)
    {

      if(fabs(gr->GetXaxis()->GetXmax()-gr->GetXaxis()->GetXmin())<10)return -1;
      TF1* gfun=gr->GetFunction("fun");
      if(gfun)
	{
	  Float_t Sfactor=pow(2,shift);
	  if(HVsetValue<0)savesign=-1;
	  Float_t startADC=gr->Eval(fabs(HVsetValue));

	  Int_t LargeRange[2]={1100,1600};
	  Int_t SmallRange[2]={0,255};
	  Int_t* Range=LargeRange;
	  if(nstb1>2)Range=SmallRange;
	  Float_t lrat=log(gcorr)+log(Sfactor);
	  Int_t best=fabs(HVsetValue);
	  Float_t dist0=fabs(lrat);
	  retval=best;
	  printf("best=%d\n",best);
	  for(int hv=Range[0]-50;hv<Range[1]+50;hv++)
	    {
	      Float_t fval=gfun->Eval(hv);
	      Float_t diff=fabs(log(fval)-log(startADC) -lrat);
	      if(diff<dist0)		{
		  dist0=diff;
		  best=hv;
		  retval=best;
		}
	    }
	};
    }

  return retval*savesign;
}
SheetHistory::SheetHistory(Float_t oldhv,Float_t newhv,Float_t oldgcorr,
			   Float_t newgcorr,Int_t oldshift,Int_t newshift,
			   Float_t gain)
{
  OldHV=oldhv;
  NewHV=newhv;
  OldGcorr=oldgcorr;
  NewGcorr=newgcorr;
  OldShift=oldshift;
  NewShift=newshift;
  Gain=gain;
  TDatime td;
  CreateDate=td.AsString();
};
TString Sheet::QTline(Bool_t nw)
{
  SheetHistory* hh;
  if(hh=GetHistory())
    {
      Int_t mone=-1;
      Int_t shift=hh->NewShift;
      if(!nw)shift=hh->OldShift;
      char str[100];
      sprintf(str,"%8d%8d%+8d%+8d\n",QTSlot,QTCard4*8+QTchannel8,mone,shift);
      return TString(str);
    };
  return TString("");
}
void Sheet::BackupHistory()
{
  SheetHistory* h=GetHistory();
  SheetHistory* next=new SheetHistory(h->OldHV,h->NewHV,h->OldGcorr,h->NewGcorr,h->OldShift,h->NewShift,h->Gain);
  Extras->Add(next);
};
void Sheet::AddHistory(Float_t oldhv,Float_t newhv,Float_t oldgcorr,
		Float_t newgcorr,Int_t oldshift,Int_t newshift,
		Float_t gain)
{
  Extras->Add(new SheetHistory(oldhv,newhv,oldgcorr,
			       newgcorr,oldshift,newshift,
			       gain));
};
TString Sheet::StringHistory(Int_t depth)
{
  char str[400];
  Int_t nentries=Extras->GetEntries();
  sprintf(str,"none\n");
  if(nentries>depth)
    {
      SheetHistory* shh=GetHistory(depth);
      Int_t Eindex=-1;
      if(Extras)Eindex=Extras->GetEntries()-1-depth;
      void* vint=(void*) shh;
      std::stringstream sss;
      sss << vint;
      TString ss0=sss.str();

      TString ss=ss0;
      sprintf(str,"%s Eindex=%d Addr=%s \n OldHV=%f NewHV=%f OldGcorr=%f NewGcorr=%f\n OldShift=%d NewShift=%d  Gain=%f \n"
	      ,(const char*) shh->CreateDate,Eindex,(const char*) ss0,shh->OldHV,shh->NewHV,
	      shh->OldGcorr,shh->NewGcorr,shh->OldShift,shh->NewShift,shh->Gain);
	};
  return TString(str);
};
SheetHistory* Sheet::GetHistory(Int_t depth)
{  
  Int_t nentries=Extras->GetEntries();
  if(nentries>depth)
    {

      SheetHistory* shh=(SheetHistory*) (*Extras)[nentries-1-depth];
      return shh;
    };
  return 0;
};
Int_t Sheet::Getshift(Bool_t nw, Int_t depth)
{
  SheetHistory* shh= GetHistory(depth);
  if(shh)
    {
      if(nw)return shh->NewShift;
      return shh->OldShift;
    };
  return 0;
};
Float_t Sheet::GetGcorr(Bool_t nw,Int_t depth)
{
  SheetHistory* shh= GetHistory(depth);
  if(shh)
    {
      if(nw)return shh->NewGcorr;
      return shh->OldGcorr;
    };
  return 0;
};
Int_t Sheet::GetHV(Bool_t nw,Int_t depth)
{
  SheetHistory* shh= GetHistory(depth);
  if(shh)
    {
      if(nw)return shh->NewHV;
      return shh->OldHV;
    };
  return 0;
};
Sheet* Sheet::FindSheetid(TObjArray* ar,Int_t id1,Int_t id2, Int_t id3, Int_t id4)
{
  if(!ar)return 0;
  TIter next(ar);
  Sheet* shreturn=0;
  Sheet* sh;
  while(sh=(Sheet*) next())
      {
	if(sh->id1==id1)
	  {
	    if(sh->id2==id2 || id1>3)
	      {
		if(sh->id3==id3 && sh->id4==id4)return shreturn=sh;
	      }
	  };
      };
  return shreturn;
};
Sheet* Sheet::FindSheetqt(TObjArray* ar,Int_t qtcrate, Int_t qtslot, Int_t qtchan)
{
  if(!ar)return 0;
  TIter next(ar);
  Sheet* sh;
  while(sh=(Sheet*) next())
      {
	if(sh->QTCrate==qtcrate)
	  {
	    if(sh->QTSlot==qtslot)
	      {
		if(qtchan==sh->QTCard4*8+sh->QTchannel8)return sh;
	      }
	  };
      };
  return 0;
};




