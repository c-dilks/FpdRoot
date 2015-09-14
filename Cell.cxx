#include "Cell.h"
ClassImp(Cell)
ClassImp(tpCellDat); 
ClassImp(CellHit);
ClassImp(tpCellRec);

Cell::Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,CalibStr* rgcor)
{
  Rdep=0;
  Mgr=0;
  RdepBaseOverride=0;
  pfiles=0;
  p_Geom=0;
  p_NomMass =0;
  Geom_ok=false;
  fpdcorr=rgcor;
  fpdgain=0;
  Setup(iEW, iNSTB,row1,col1,rnum);
};

Cell::Cell()
{
  Mgr=new RunDepMgr();
  Rdep=new RunDepCor();
  ReScaledHist=new TH1D();
  p_Mcell=new TH1F();
  p_MvsZ=new TH2F();
  p_MvsD=new TH2F();
  p_MvsE=new TH2F();
  p_MvsY=new TH2F();
  p_EphXY=new TH1F();
  p_RawAdc=new TH1F();
  p_RawAdcL=new TH1F();
  p_adc=new TH1D();
  p_adcLed=new TH1D();
  p_dxdyM=new TH3F();
  p_dxdyE=new TH3F();
  p_dxdyER=new TH3F();
  p_Nphthis=new TH1F();
  p_NphAll=new TH1F();
  p_DetE=new TH1F();
  p_FmsE=new TH1F();
  fpdcorr=new CalibStr();
  fpdgain=new CalibStr();
  p_Geom=new Geom();
  pfiles=new FilesSet();
  CellD=new TTree();
  for(int n=0;n<4;n++)
    {
      gr[n]=new TGraph();
      gradc[n]=new TH2F();
    };
  //  if(CellD)InitBranches();
  //  p_gFME=new TGraphErrors();
};
  Cell::Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,TFile* fl,CalibStr* rgcor,CalibStr* NominalMass)
{
  Rdep=0;
  Mgr=0;
  p_NomMass =NominalMass;
  fpdcorr=rgcor;
  fpdgain=0;
  pfiles=0;
  p_Geom=PGeom;
  Geom_ok=false;
  if(p_Geom)Geom_ok=true;
  Setup(iEW, iNSTB,row1,col1,rnum);
  Setpfiles();
  if(fl!=0)FindHists(fl);
};
Cell::Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,TFile* fl,CalibStr* rgcor)
{
  Rdep=0;
  Mgr=0;
  p_NomMass =0;
  fpdcorr=rgcor;
  fpdgain=0;
  pfiles=0;
  p_Geom=PGeom;
  Geom_ok=false;
  if(p_Geom)Geom_ok=true;
  FitMass=0.;
  erFitMass=0.;
  FitMassWidth=0.;
  MassEntries=0.;
  hiRawSlope=0.;
  hiRawIntdorig=0.;
  ReFitQuality0to3=0.;
  erRawSlope=0.;
  hiRawInt=0.;
  erhiRawInt=0.;
  FitSlope=0.;
  
  Setup(iEW, iNSTB,row1,col1,rnum);
  Setpfiles();
  if(fl!=0)FindHists(fl);
};
Cell::Cell(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum,Geom* PGeom,CalibStr* rgain,CalibStr* rgcor)
{
  Mgr=0;
  Rdep=0;
  p_NomMass =0;
  fpdcorr=rgcor;
  fpdgain=rgain;
  pfiles=0;
  p_Geom=PGeom;
  Geom_ok=false;
  if(p_Geom)Geom_ok=true;
  Geom_ok=true;
  FitMass=0.;
  erFitMass=0.;
  FitMassWidth=0.;
  MassEntries=0.;
  hiRawSlope=0.;
  hiRawIntdorig=0.;
  ReFitQuality0to3=0.;
  erRawSlope=0.;
  hiRawInt=0.;
  erhiRawInt=0.;
  FitSlope=0.;
  
  Setup(iEW, iNSTB,row1,col1,rnum);
  Setpfiles();
};
Float_t Cell::TriggerPeak()
{
  if(p_DetE==0)return 0.;
  Int_t nbins=p_DetE->GetNbinsX();
  Double_t* entries= p_DetE->GetIntegral();
  Float_t lastmax=0;
  Float_t downsigma=0;
  Float_t peakloc=0.;
  LastE_Det=100.;
  LastE_Det_Error=0;
  Int_t cnt=0;
  Det_Peak_Loc=0.;
  Float_t NominalThreshold[4]={30,30,40,40};
  Float_t NRows[4]={34,34,24,24};
  Float_t NominalThresh=NominalThreshold[Instb-1];
  Float_t dnsig=0;
  Float_t Threshgcorr=NominalThresh*gcorr;

  if(Threshgcorr>=100)return (peakloc=100.);
  if(Threshgcorr<=0)return 0.;
  
  while(cnt<2 && nbins>2)
    {
      nbins--;
      Float_t val=p_DetE->GetBinContent(nbins);
      cnt+= (Int_t) val;
      if(cnt>1)
	{
	  LastE_Det=p_DetE->GetBinLowEdge(nbins);
	};
    };
  Float_t val=0;
  Int_t nb0=0;
  while(val<.5 && nb0>0)
    {
      nb0++;
      Float_t val=p_DetE->GetBinContent(nb0);
    };
  
  if(nbins<=2)LastE_Det=0;
  Float_t newpeak=Threshgcorr;
  Float_t Slp1=0;
  Float_t Icpt1=0;
  Float_t Slp2=0;
  Float_t Icpt2=0;
  Int_t ntrigbin=p_DetE->FindBin(newpeak);
  Float_t trght=p_DetE->GetBinContent(ntrigbin);
  Int_t ntcnt=ntrigbin;
  for(Int_t ib=ntrigbin-2;ib<ntrigbin+10;ib++)
    {
      if(ib>0 && ib<100)
	{
	  Int_t tmpcnt=(Int_t) p_DetE->GetBinContent(ib);
	  if(tmpcnt>trght)
	    {
	      ntcnt=ib;
	      trght=tmpcnt;
	    };
	};
    };
  ntrigbin=ntcnt;
  newpeak=p_DetE->GetBinLowEdge(ntrigbin);
  Det_Peak_Loc=newpeak;
  Int_t nendbin=p_DetE->FindBin(LastE_Det);
  Float_t nover=p_DetE->Integral(ntrigbin+1,nendbin);
  Float_t nunder=p_DetE->Integral(1,ntrigbin);
  Chi2_Det=0;
  NDF_Det=0;
  if(nunder>20)
    {
      p_DetE->Fit("expo","","",5,newpeak);
      TF1* f1=p_DetE->GetFunction("expo");
      Icpt2=f1->GetParameter(0);
      Slp2=f1->GetParameter(1);
      
      Float_t NDF2=f1->GetNDF();
      Float_t Chi2_2=f1->GetChisquare();
      if(NDF_Det>20)
	{
	  Float_t chi_per_DF=Chi2_Det/NDF_Det;
	  if(chi_per_DF>2 && newpeak>19)
	    {
	      p_DetE->Fit("expo","","",newpeak-8,newpeak);
	      TF1* f10=p_DetE->GetFunction("expo");
	      Icpt2=f10->GetParameter(0);
	      Slp2=f10->GetParameter(1);
	    };
	};
      
    };
  Slope_Det=Slp2;
  if(nover<2)
    {
      return 0;
    }
  else
    {
      p_DetE->Fit("expo","","",newpeak+1,LastE_Det);
      TF1* f1=p_DetE->GetFunction("expo");
      NDF_Det=f1->GetNDF();
      Chi2_Det=f1->GetChisquare();
      if(NDF_Det>5)
	{
	  Float_t chi_per_DF=Chi2_Det/NDF_Det;
	  Icpt1=f1->GetParameter(0);
	  Slp1=-f1->GetParameter(1);
	  Float_t par0err=f1->GetParError(0);
	  Float_t par1err=f1->GetParError(1);
	  if(Icpt1>0 && Slp1>0)
	    {
		Float_t le=Icpt1/Slp1;
	      LastE_Det_Error=LastE_Det*sqrt(pow(par0err/Icpt1,2)+
					     pow(par1err/Slp1,2));
	      if(LastE_Det_Error<10)
		{
		  LastE_Det=le;
		}
	      else
		{
		  LastE_Det_Error=0.;
		};
	    };
	};
      Slope_Det=Slp1;
      Det_Peak_Loc=newpeak;
      return newpeak;
    };
  Det_Peak_Loc=newpeak;
  return newpeak;
};

Bool_t Cell::Setup(Int_t iEW, Int_t iNSTB,Int_t row1,Int_t col1,Int_t rnum)
{
  p_dxdyM=0;
  p_dxdyE=0;
  p_dxdyER=0;
  sprintf(ReasonForChange,"None\n");
  ReScaledHist=0;
  Iew=iEW;
  Instb=iNSTB;
  Row1=row1;
  Col1=col1;
  Rnum=rnum;
  p_Mcell=0;
  p_MvsZ=0;
  p_MvsD=0;
  p_MvsE=0;
  p_MvsY=0;
  p_EphXY=0;
  p_RawAdc=0;
  p_RawAdcL=0;
  p_adc=0;
  p_adcLed=0;
  gcorr=0;
  gain=0;
  p_Nphthis=0;
  p_NphAll=0;
  p_DetE=0;
  p_FmsE=0;
  if(fpdcorr!=0)
    {
      gcorr=fpdcorr->GetValue(Iew,Instb,row1-1,col1-1);
      gcorrOrig=gcorr;
    };  
  if(fpdgain!=0)
    {
      gain=fpdgain->GetValue(Iew,Instb,row1-1,col1-1);
    };
  Mass0=.135;
  if(p_NomMass!=0)
    {
      Mass0=p_NomMass->GetValue(Iew,Instb,row1-1,col1-1);
      if(Mass0<.04)Mass0=.135;
    };
  CellD=0;
  for(int n=0;n<4;n++)
    {
      gr[n]=new TGraph();
      gradc[n]=new TH2F();
    };
  return true;
};

Cell::~Cell()
{
  if(Mgr)delete Mgr;
  if(Rdep)delete Rdep;
  if(p_Mcell)delete p_Mcell;
  if(p_MvsZ)delete p_MvsZ;
  if(p_MvsD)delete p_MvsD;
  if(p_MvsE)delete p_MvsE;
  if(p_MvsY)delete p_MvsY;
  if(p_EphXY)delete p_EphXY;
  if(p_RawAdc)delete p_RawAdc;
  if(p_RawAdcL)delete p_RawAdcL;
  if(p_adc)delete p_adc;
  if(p_adcLed)delete p_adcLed;
  if(p_Nphthis)delete p_Nphthis;
  if(p_NphAll)delete p_NphAll;
  if(p_DetE)delete p_DetE;
  if(p_FmsE)delete p_FmsE;
  if(ReScaledHist)delete ReScaledHist;
  for(int n=0;n<4;n++)
    {
      if(gr[n])delete gr[n];
      if(gradc[n])delete gradc[n];
    };
};

 Bool_t  Cell::FindHists(TFile* tfile)
{
  tfile->cd();
  char name[100];
  char* ns[2];
  ns[0]="N";
  ns[1]="S";
  Int_t n_s;
  if(Iew!=1)return false;

  if(Instb==1){n_s=0;}
  else if(Instb==2){n_s=1;}
  else return false;

  sprintf(name,"Mcell%d_%d_%s",Col1,Row1,ns[n_s]);
  p_Mcell=((TH1F*) tfile->FindKey(name)->ReadObj());
  sprintf(name,"EphXY%d_%d_%s",Col1,Row1,ns[n_s]);
  p_EphXY=((TH1F*) tfile->FindKey(name)->ReadObj());
  sprintf(name,"RawAdc%d_%d_%s",Col1,Row1,ns[n_s]);
  p_RawAdc=((TH1F*) tfile->FindKey(name)->ReadObj());
  sprintf(name,"RawAdcL%d_%d_%s",Col1,Row1,ns[n_s]);
  p_RawAdcL=((TH1F*) tfile->FindKey(name)->ReadObj());
  FitMcell();
  FitRawCell();
  FitRawLCell(10/gcorr,40/gcorr);
  Newcorr();
};

Bool_t Cell::AddToReason(char* reason)
{
  TString oldreason=ReasonForChange;
  const char* oreason=(const char*) oldreason;
  sprintf(ReasonForChange,"%s - %s \n",oreason,reason);
  return true;
};

Float_t Cell::FitMcell2(Float_t im1,Float_t im2)
{
  Float_t m1=im1;
  Float_t m2=im2;
  Int_t Im1,Im2;
  Im1=p_Mcell->FindBin(m1);
  Im2=p_Mcell->FindBin(m2);
  MassFitError=0;
  TF1 gfun("gfun","[0]* exp(-(x-[1])*(x-[1])/2/[2]/[2])");
  gfun.SetParameter(0,10.);
  gfun.SetParameter(1,.135);
  gfun.SetParameter(2,.03);
  gfun.SetParLimits(2,.025,.065);
  gfun.SetParLimits(1,.02,.35);
  Float_t intval=p_Mcell->Integral(Im1,Im2);
  if(intval<3)
    {
      gfun.SetParameter(0,.0);
      gfun.SetParLimits(0,0.,1.);      
    };
  p_Mcell->Fit("gfun","L","",im1,im2);

  FitMass=p_Mcell->GetFunction("gfun")->GetParameter(1);
  if(FitMass<.2 && im2>.22)
    {
      p_Mcell->Fit("gfun","L","",im1,.22);
    };
  FitMass=p_Mcell->GetFunction("gfun")->GetParameter(1);
  MassFitError=p_Mcell->GetFunction("gfun")->GetParError(1);
  Massmaxcontents=p_Mcell->GetFunction("gfun")->GetParameter(0);
  Int_t bmax=p_Mcell->FindBin(FitMass);
  Int_t low=p_Mcell->FindBin(im1);
  Int_t hi=p_Mcell->FindBin(im2);
  Float_t integral=p_Mcell->Integral(low,hi);
  MassEntries=p_Mcell->GetEntries();
  Masspeakfraction=0.;
  if(MassEntries>0)  Masspeakfraction=integral/MassEntries;
  Float_t mean=p_Mcell->GetMean();
  FitMassWidth=p_Mcell->GetFunction("gfun")->GetParameter(2);
  erFitMass=p_Mcell->GetFunction("gfun")->GetParError(1);

  return FitMass;
};

Float_t Cell::FitMcell(Float_t im1,Float_t im2)
{
  Float_t m1=im1;
  Float_t m2=im2;
  MassFitError=0;
  Int_t bmax=p_Mcell->GetMaximumBin();
  Float_t cmax=p_Mcell->GetBinCenter(bmax);
  Massmaxcontents=p_Mcell->GetBinContent(bmax);
  Masspeakfraction=0;
  MassEntries=p_Mcell->GetEntries();
  if(Massmaxcontents>0){
    for(Int_t j=-4;j<5;j++){Masspeakfraction+=
			      p_Mcell->GetBinContent(bmax+j)/MassEntries;};
  };
  Float_t mean=p_Mcell->GetMean();
  if(Masspeakfraction>.5){
    m1=p_Mcell->GetBinCenter(bmax-4);
    m2=p_Mcell->GetBinCenter(bmax+4);
  };
  if(Masspeakfraction< .5 && cmax<.07)
    {
      cmax=Mass0;
      m1=.08;
      m2=.19;  
      bmax=p_Mcell->FindBin(Mass0);
      Masspeakfraction=0;
      if(MassEntries>0){
	for(Int_t j=-3;j<4;j++){Masspeakfraction+=
				  p_Mcell->GetBinContent(bmax+j)/MassEntries;};
      };
      
    };
  if(Masspeakfraction< .3 && cmax>.25)
    {
      cmax=Mass0;
      m1=.08;
      m2=.19;  
      bmax=p_Mcell->FindBin(Mass0);
      Masspeakfraction=0;
      if(MassEntries>0){
	for(Int_t j=-3;j<4;j++){Masspeakfraction+=
				  p_Mcell->GetBinContent(bmax+j)/MassEntries;};
      };
      
    };

  p_Mcell->Fit("gaus","LE","",m1,m2);
  Float_t sig=p_Mcell->GetFunction("gaus")->GetParameter(2);
  p_Mcell->GetFunction("gaus")->SetParLimits(2,.015,.05);
  if(sig<.015 ||sig>.05)
    {
      p_Mcell->GetFunction("gaus")->SetParameter(2,.025);
      p_Mcell->Fit("gaus","LE","",m1,m2);
    };
  FitMass=p_Mcell->GetFunction("gaus")->GetParameter(1);
  if(FitMass<m1)  
    {
      m1=m1-.02;
      p_Mcell->Fit("gaus","LE","",m1,m2);     
      FitMass=p_Mcell->GetFunction("gaus")->GetParameter(1);  
    };   
  Int_t  loop=0;
  while(loop<3)
    { 
      loop++;
      Int_t bn0=p_Mcell->FindBin(FitMass);
      if(FitMass<m2)
	{
	  Int_t bn0=p_Mcell->FindBin(FitMass);
	  m1=p_Mcell->GetBinCenter(bn0-4);
	  m2=p_Mcell->GetBinCenter(bn0+4); 
	}
      else
	{
	  m2=m2+.01;
	};
      p_Mcell->Fit("gaus","LE","",m1,m2);
      bmax=p_Mcell->FindBin(FitMass);
      FitMass=p_Mcell->GetFunction("gaus")->GetParameter(1);  
      Float_t center=(FitMass-m1)/(m2-m1);
      
      printf("row1,col1,FitM=%d %d %f : m1,m2,center=%f %f %f \n",Row1,Col1,FitMass,m1,m2,center);
      if(.3<center && center<.7)break;

    };
  bmax=p_Mcell->FindBin(FitMass);

  Masspeakfraction=0;
  if(MassEntries>0){
    for(Int_t j=-3;j<4;j++){Masspeakfraction+=
			      p_Mcell->GetBinContent(bmax+j)/MassEntries;};
    
  };
  FitMassWidth=p_Mcell->GetFunction("gaus")->GetParameter(2);
  erFitMass=p_Mcell->GetFunction("gaus")->GetParError(1);
  return FitMass;
};

Float_t Cell::FitRawCell()
{
  if(fpdcorr->GetValue(Iew,Instb,Row1-1,Col1-1)<=0)return 0.;
  Beyond120=p_RawAdc->Integral(120,200);
  Beyond145=p_RawAdc->Integral(145,200);
  TF1 forig("forig","pol2",-4.4,-3.3);
  forig.SetParameter(0,-120.14); 
  forig.SetParameter(1,-61.24);
  forig.SetParameter(2,-6.963);
  Float_t orig=forig.Eval(Eta);
  TF1 fexp("fexp","expo",145,200.);
  //  fexp.SetParameter(0,orig);
  //  fexp.SetParLimits(0,orig-3.0001,orig+3.0001);
  p_RawAdc->Fit("fexp","LE","",145,220);
  hiRawSlope=p_RawAdc->GetFunction("fexp")->GetParameter(1);
  erRawSlope=p_RawAdc->GetFunction("fexp")->GetParError(1);
  hiRawInt=p_RawAdc->GetFunction("fexp")->GetParameter(0);

  printf("orig=%f hiRawInt=%f \n",orig,hiRawInt);
  raw120=hiRawInt+120*hiRawSlope/gcorr;
  IntRawFitwithgcorr=exp(raw120)/(-1*hiRawSlope)*gcorr;
  printf("r120=%f corrected int=%f \n",exp(raw120),IntRawFitwithgcorr);
  hiRawIntdorig=hiRawInt-orig;

  erhiRawInt=p_RawAdc->GetFunction("fexp")->GetParError(0);

  
  return hiRawSlope;
};
Float_t Cell::FitRawLCell(Float_t a1,Float_t a2)
{
  p_RawAdcL->Fit("expo","LE","",a1,a2);
  BeyondL10=p_RawAdcL->Integral((Int_t)a1,255);
  hiRawSlopeL=p_RawAdcL->GetFunction("expo")->GetParameter(1);
  erRawSlopeL=p_RawAdcL->GetFunction("expo")->GetParError(1);
  return hiRawSlopeL;
};
Bool_t Cell::Setpfiles(FilesSet* _pfiles)
{
  pfiles=_pfiles;
  if(!Geom_ok)
    {
      if(!pfiles->p_Geom()->check())
	{
	  Geom_ok=false;
	  return false;
	};
      Geom_ok=true;
      p_Geom=new Geom(pfiles);
    };
  
  Float_t towerwidth=*p_Geom->FpdTowWid(Iew,Instb);
  xLocal=towerwidth*(Col1-.5);
  yLocal=towerwidth*(Row1-.5);
  Float_t z=0.;
  TVector3 loctow(xLocal,yLocal,z);
  TVector3 vtow=p_Geom->GlobalXYZ(Iew,Instb,loctow);
  Eta=vtow.PseudoRapidity();
  Phi=vtow.Phi();
  xGlobal=vtow.X();
  yGlobal=vtow.Y();
  return true;
};

Float_t Cell::Newcorr()
{
  if(gcorr<=0)return gcorr;
  Float_t y[3];
  y[0]=hiRawSlopeL/(-.12)/gcorr-.5*(Eta+3.8);
  y[1]=pow(Mass0/FitMass,1.);
  y[2]=-(hiRawSlope/.052/gcorr)-.286*(Eta+4.05);
  Float_t Usey[3];
  Float_t Err_y[3];
  Err_y[0]=fabs(erRawSlopeL/(-.12));
  Err_y[1]=fabs(2*erFitMass/Mass0);
  Err_y[2]=fabs(erRawSlope/.052);
  Usey[2]=0.;
  Usey[1]=0.;
  Usey[0]=0.;
  Bool_t VeryGoodMass=false;
  if(Masspeakfraction>.55 && Err_y[1]<.01)VeryGoodMass=true;
  if(y[1]>.1 && y[1]<5. && Err_y[1]<.2 && Masspeakfraction>.2)
    {
      Usey[1]=.01/(fabs(Err_y[1])+.01);
      Usey[1]=Usey[1]*TMath::Min(Masspeakfraction/.8,1.);
    };

  if(y[2]>.1 && y[2]<5. && Err_y[2]<.2 && Beyond145>100)
    {
      Usey[2]=.01/(fabs(Err_y[2])+.01);
      Float_t denom=500;
      if(VeryGoodMass)denom=10000;      
      Usey[2]=Usey[2]*TMath::Min(sqrt(Beyond145/500),1.);
    };


  if(y[0]>.1 && y[0]<5. && Err_y[0]<.1 && BeyondL10>1000 && erRawSlopeL/(-.12)<.05)
    {
      Usey[0]=.002/(fabs(Err_y[0])+.005);
    };
  if(VeryGoodMass)Usey[0]=0.;
  // only mass

  Usey[0]=0;
  Usey[2]=0;
  printf("usey[0-2]: %f %f %f ",Usey[0],Usey[1],Usey[2]);
  Float_t Usedy=Usey[0]+Usey[1]+Usey[2];
  ReFitQuality0to3=Usedy;
  printf("ReRitQ=%f \n",ReFitQuality0to3);
  gcorrNew=gcorr;
  if(Usedy<=0)return gcorr;
  Usey[0]=Usey[0]/Usedy;
  Usey[1]=Usey[1]/Usedy;
  Usey[2]=Usey[2]/Usedy;
  gcorrNew=gcorr*(Usey[0]*y[0]+Usey[1]*y[1]+Usey[2]*y[2]);
  // limit to 10% change
  if(gcorrNew>1.3* gcorr)gcorrNew=1.3*gcorr;
  if(gcorrNew<.75 *gcorr)gcorrNew=.75*gcorr;

  printf("old gcorr=%f new gcorr =%f \n",gcorr,gcorrNew);
  fpdcorr->SetValue(Iew,Instb,Row1-1,Col1-1,gcorrNew);
  return gcorrNew;
};
Int_t Cell::Compare(const TObject* p_Ob) const
{
  if(strcmp(p_Ob->ClassName(),"Cell")!=0)return -1;
  Cell* pcell=(Cell*) p_Ob;
  if(pcell->Eta<Eta)return -1;
  return 1;
};

TH1D* Cell::ReScaleHist(TH1D* p_hist,Float_t xfactor,Int_t minbin,Bool_t Normalize)
{
  printf("ReScaleHist called with xfactor=%f minbin=%d \n",xfactor,minbin);
  p_hist->Print();
  if(ReScaledHist)delete ReScaledHist;
  
  Int_t nb1x=p_hist->GetNbinsX();
  Float_t lowx=p_hist->GetBinLowEdge(1);
  Float_t highx0=p_hist->GetBinLowEdge(nb1x+1);
  Float_t highx=lowx+(highx0-lowx)*xfactor;
  ReScaledHist=new TH1D("ReScaledHist","ReScaledHist",nb1x,lowx,highx);
  Float_t binw0=p_hist->GetBinWidth(1);
  Float_t binw=ReScaledHist->GetBinWidth(1);
  Float_t ReScaledBins[10000];
  Float_t Sumsq=0;
  for(Int_t j=1;j<nb1x+1;j++)
    {
      if(j>=10000)continue;
      Float_t bincontent=0;
      Float_t ledg=ReScaledHist->GetBinLowEdge(j);
      Float_t hiedg=ledg+binw;
      Int_t olow=p_hist->FindBin(ledg);
      Int_t ohi=p_hist->FindBin(ledg+binw);
      bincontent=p_hist->Integral(olow,ohi);
      Float_t lowe0=p_hist->GetBinLowEdge(olow);
      Float_t subtract=p_hist->GetBinContent(olow)*(ledg-lowe0)/binw0;
      Float_t hie0=p_hist->GetBinLowEdge(ohi)+binw0;
      subtract+=p_hist->GetBinContent(ohi)*(hie0-hiedg)/binw0;
      bincontent=bincontent-subtract;
      if(minbin>j)bincontent=0;
      ReScaledBins[j]=bincontent;
      ReScaledHist->SetBinContent(j,bincontent);
      Sumsq+=bincontent*bincontent;
    };

 
  delete ReScaledHist;
  if(nb1x>9998)
    {
      highx0=lowx+(highx0-lowx)*9998/nb1x;
      nb1x=9998;
    };
  ReScaledHist=new TH1D("ReScaledHist","ReScaledHist",nb1x,lowx,highx0);
  
  
  Float_t sqrsumsq=1;
  if(Normalize)sqrsumsq=sqrt(Sumsq);
  
  for(Int_t j=1;j<nb1x+1;j++)ReScaledHist->SetBinContent(j,ReScaledBins[j]/sqrsumsq);
  return ReScaledHist;  
};
void Cell::CreateBranches()
{
  CellD->Branch("Mcell",&celldat.Mcell,"Mcell/F");
  CellD->Branch("Epair",&celldat.Epair,"Epair/F");
  CellD->Branch("Ypair",&celldat.Ypair,"Ypair/F");
  CellD->Branch("E1",&celldat.E1,"E1/F");
  CellD->Branch("E2",&celldat.E2,"E2/F");
  CellD->Branch("X1",&celldat.X1,"X1/F");
  CellD->Branch("X2",&celldat.X2,"X2/F");
  CellD->Branch("Y1",&celldat.Y1,"Y1/F");
  CellD->Branch("Y2",&celldat.Y2,"Y2/F");
  CellD->Branch("NSTB",&celldat.NSTB,"NSTB/I");
  CellD->Branch("br_nSavedHits",&(celldat.nSavedHits),"nSavedHits/I");
  CellD->Branch("br_SavedHits",&(celldat.SavedHits),"SavedHits[100]/I");
  CellD->Branch("br_ievt",&(celldat.ievt),"ievt/I");
  CellD->Branch("br_Rnum",&(celldat.Rnum),"Rnum/I");
};
void Cell::InitBranches()
{
  CellD->SetBranchAddress("Mcell",&celldat.Mcell);
  CellD->SetBranchAddress("Epair",&celldat.Epair);
  CellD->SetBranchAddress("Ypair",&celldat.Ypair);
  CellD->SetBranchAddress("E1",&celldat.E1);
  CellD->SetBranchAddress("E2",&celldat.E2);
  CellD->SetBranchAddress("X1",&celldat.X1);
  CellD->SetBranchAddress("X2",&celldat.X2);
  CellD->SetBranchAddress("Y1",&celldat.Y1);
  CellD->SetBranchAddress("Y2",&celldat.Y2);
  CellD->SetBranchAddress("NSTB",&celldat.NSTB);
  CellD->SetBranchAddress("br_nSavedHits",&(celldat.nSavedHits));
  CellD->SetBranchAddress("br_SavedHits",&(celldat.SavedHits));
  CellD->SetBranchAddress("br_ievt",&(celldat.ievt));
  CellD->SetBranchAddress("br_Rnum",&(celldat.Rnum));

};
void Cell::InitTree(const char* name)
{
  if(CellD)delete CellD;
  printf("InitTree local file = %s to file: %s\n",name,gDirectory->GetName());
  CellD=new TTree(name,name);
  CreateBranches();
};
void Cell::FillTree(tpCellDat cd)
{
  if(CellD)
    {
      celldat=cd;
      CellD->Fill();
    };
};
TGraphErrors* Cell::FitMvsE(TH2F* hme,Int_t ndiv,char* opt, Int_t nstart,Int_t nend,TH1D* slicearray)
{
  Int_t ne=hme->GetNbinsX();
  if(nend>ne)nend=ne;
  Int_t nm=hme->GetNbinsY();
  Int_t nw=(nend-nstart+1)/ndiv;
  
  TGraphErrors* p_gFME=0;
  Float_t e[100],m[100],de[100],dm[100];
  Int_t ndat=0;
  Float_t Emax=hme->GetXaxis()->GetXmax();
  Float_t Emin=hme->GetXaxis()->GetXmin();
  Float_t dEm=(Emax-Emin)/ndiv;
  Float_t dde=(Emax-Emin)/ne;
  Float_t dn=(1.* ne)/ndiv;
  Float_t flo=nstart;
  for(int j=0;j<ndiv;j++)
  {
    char nm[100];
    sprintf(nm,"slice_%d",j);
    printf("%s\n",nm);
    Int_t lo=flo+.5;
    Int_t hi=flo+dn;
    flo=flo+dn+.5;
    e[j]=m[j]=de[j],dm[j]=0;
    TH1D* hsl=hme->ProjectionY(nm,lo,hi);
    printf("range =(%d .. %d ) ",lo,hi);
    printf("N(%d)=%f\n",j,hsl->Integral());

    if(hsl->Integral()>3)
      {
	hsl->Fit("gaus",opt,"",.09,.16);
	e[ndat]=Emin+(hi+lo)/2.*dde;
	printf("e[%d]=%f \n",ndat,e[ndat]); 
	de[ndat]=dEm/2.;
	m[ndat]=hsl->GetFunction("gaus")->GetParameter(1);
	dm[ndat]=hsl->GetFunction("gaus")->GetParError(1);
	if(m[ndat]>.15)
	  {
	    hsl->Fit("gaus",opt,"",.1,.25);
	    e[ndat]=Emin+(hi+lo)/2.*dde;
	    printf("e[%d]=%f \n",ndat,e[ndat]); 
	    de[ndat]=dEm/2.;
	    m[ndat]=hsl->GetFunction("gaus")->GetParameter(1);
	    dm[ndat]=hsl->GetFunction("gaus")->GetParError(1);
	  };
	if(m[ndat]<.125)
	  {
	    hsl->Fit("gaus",opt,"",.06,.15);
	    e[ndat]=Emin+(hi+lo)/2.*dde;
	    printf("e[%d]=%f \n",ndat,e[ndat]); 
	    de[ndat]=dEm/2.;
	    m[ndat]=hsl->GetFunction("gaus")->GetParameter(1);
	    dm[ndat]=hsl->GetFunction("gaus")->GetParError(1);
	  };
	ndat++;
      };
    if(slicearray!=0)
      {
	*slicearray=*hsl;
	slicearray++;
      }
    
    delete hsl;
  };
  if(ndat==0)ndat=1;
  p_gFME=new TGraphErrors(ndat,e,m,de,dm);
  p_gFME->GetYaxis()->SetRangeUser(0.,.2);
  return p_gFME;
};
TMatrixT<float> Cell::FillFMSEnergy(Int_t NSTB)
{  

  TFile* SaveDefaultFile=gROOT->GetFile();  
  int nrows=34;
  int ncols=17;
  if(NSTB>2)
    {
      nrows=24;
      ncols=12;
    };
  TMatrixT<float> tm(nrows,ncols);
  Bool_t RunCorrect=false;
  
  if(Mgr==0)
    {
      if(gSystem->Getenv("RunDepPath"))
	{ 
	  printf("reread RunDepPath\n");
	  TFile* RdepFile=new TFile("$RunDepPath");//set in SetFMSEnv
	  Mgr=(RunDepMgr*) RdepFile->Get("Mgr");
	  Mgr->RunDepPath="$RunDepPath";
	  delete RdepFile;
	  Rdep=0;
	};
    }
  
  if(Mgr&& Rdep &&celldat.Rnum>=12000000)
    {
      if(celldat.Rnum==Rdep->RunNumber)RunCorrect=true;
      printf("RunCorrect set to true\n");
      if(Mgr) printf ("Mgr nonzero\n");
    };

  if((! RunCorrect) && (celldat.Rnum>=12000000)&&Mgr )
    {
      Rdep= Mgr->SetRdep(celldat.Rnum);
      if(RdepBaseOverride>12000000)
	{
	  Mgr->SetBase(RdepBaseOverride);
	};
    }
  else Rdep=0;
    
  

  for(int nadc=0;nadc<celldat.nSavedHits;nadc++)
    {
      //      printf("adc %d of %d\n",nadc,celldat.nSavedHits);
      unsigned int s=(unsigned int) celldat.SavedHits[nadc];
      Int_t sew,snstb,srow,scol,sadc;
      Float_t energy;
      sew=1;
      if(s&&0x80000000)sew=2;
      if(sew!=2)continue;
      s=s&0x7FFFFFFF;
      snstb=((s/0x10000000)&7)+1;
      if(snstb!=NSTB)continue;
      srow= ((s/0x00400000)&0x3F)+1;
      scol= ((s/0x00010000)&0x3F)+1;
      sadc= s&0xFFF;
      energy=sadc;
      //      printf("row=%d col=%d sadc=%d   ",srow,scol,sadc);
      if(srow<1 || srow>nrows ||scol<1||scol>ncols)
	printf("Bad Packed Data: nstb=%d row1=%d col=%d adc=%d\n",snstb,srow,scol,sadc);      
      if(Rdep&&Mgr)
	{
	  Int_t nseg=(celldat.ievt/10000)+1;
	  Int_t nevt=celldat.ievt%10000;
	  Float_t gcr=Rdep->GetCor(nseg,nevt,snstb,srow-1,scol-1,fpdcorr,2);
	  Float_t ledfac=Mgr->ledFactor[snstb-1][srow-1][scol-1];
	  energy=energy*gcr*ledfac;
	}
      else
	{
	  energy=energy*fpdcorr->GetValue(sew,snstb,srow-1,scol-1);
	};
      
      energy=energy*fpdgain->GetValue(sew,snstb,srow-1,scol-1);
      if(energy<Emin)energy=0.;
      tm(srow-1,scol-1)=energy;
    }
    //  tm.Print();
  if(SaveDefaultFile)SaveDefaultFile->cd();
  return tm;
};
TObjArray* Cell::CellHitList()
{   
  Int_t nhits=celldat.nSavedHits;
  TObjArray* rar=new TObjArray(nhits+1,0);
  rar->SetOwner();
  for(int nadc=0;nadc<celldat.nSavedHits;nadc++)
    {
      rar->AddAt(new CellHit(celldat.SavedHits[nadc]),nadc);
    };
  return rar;
};
TMatrixT<float> Cell::FillFMSADC(Int_t NSTB)
{
   int nrows=34;
   int ncols=17;
   if(NSTB>2)
     {
       nrows=24;
       ncols=12;
     };
   TMatrixT<float> tm(nrows,ncols);
   for(int nadc=0;nadc<celldat.nSavedHits;nadc++)
     {
       unsigned int s=(unsigned int) celldat.SavedHits[nadc];
       Int_t sew,snstb,srow,scol,sadc;
       sew=1;
       if(s&&0x80000000)sew=2;
       if(sew!=2)continue;
       s=s&0x7FFFFFFF;
       snstb=((s/0x10000000)&7)+1;
       if(snstb!=NSTB)continue;
       srow= ((s/0x00400000)&0x3F)+1;
       scol= ((s/0x00010000)&0x3F)+1;
       sadc= s&0xFFF;
       tm(srow-1,scol-1)=sadc;

     };
   return tm;
};
void Cell::DrawFMSADC(Int_t instb,Bool_t UseEnergy)
{
  gStyle->SetPalette(1);
  if(Iew!=2 || instb<1 || instb>4)return;
  if(gr[instb-1])delete gr[instb-1];
  if(gradc[instb-1])delete gradc[instb-1];
  if(UseEnergy)
    {
      TMatrixT<float> tm=FillFMSEnergy(instb);
      gradc[instb-1]=new TH2F(tm);
      char str[100];
      sprintf(str,"NSTB%d_r%d_c%d_En=%4.1f\n",Instb,Row1,Col1,tm.Sum());
      gradc[instb-1]->SetTitle(str);
    }
  else
    {
      TMatrixT<float> tm=FillFMSADC(instb);
      gradc[instb-1]=new TH2F(tm);
    };
  gradc[instb-1]->SetStats(0);
  Float_t xhit[2],yhit[2];
  gradc[instb-1]->Draw("zcol");
  gradc[instb-1]->Draw("samebox");
  xhit[0]=celldat.X1;
  xhit[1]=celldat.X2;
  yhit[0]=celldat.Y1;
  yhit[1]=celldat.Y2;
  gr[instb-1]=new TGraph(2,xhit,yhit);
  gr[instb-1]->Draw("*");

};
CellHit::CellHit(unsigned int sp)
{
  unsigned int s=sp;
  Int_t sew,snstb,srow,scol,sadc;
       EW=1;
       if(s&&0x80000000)EW=2;
       if(sew!=2)
	 {
	   NSTB=0;
	   Row=0;
	   Col=0;
	   ADC=0;
	 };
       s=s&0x7FFFFFFF;
       NSTB=((s/0x10000000)&7)+1;
       Row= ((s/0x00400000)&0x3F)+1;
       Col= ((s/0x00010000)&0x3F)+1;
       ADC= s&0xFFF;
};
TTree* Cell::ReRecon(tpCellRec* cr,CalibStr* Fpdcorr,Int_t SampleMax,CalibStr* pRank,TTreeFormula* precut,Bool_t UseQcut)
{
  Int_t NskipPrint=50;
  Int_t NROWS[4]={34,34,24,24};
  //  precut->Print();
  Bool_t evplot=false;
  Int_t Nevplot=0;
  if(const char* ReRecEvDisplayCnt=gSystem->Getenv("ReRecEvDisplayCnt"))
    {
      printf("ReRecEvDisplayCnt=%s \n",ReRecEvDisplayCnt);
      if(sscanf(ReRecEvDisplayCnt,"%d",&Nevplot))
	{
	  evplot=true;
	  printf("evplot set to true\n");
	}
    };
  Int_t fillcnt=0;
  if(Fpdcorr)fpdcorr=Fpdcorr;

  TTree* trr=CellD->CloneTree(0);
  trr->SetDirectory(0);

  trr->Branch("br_mr",&(cr->mr),"mr/F");
  trr->Branch("br_er",&(cr->er),"er/F");
  trr->Branch("br_e1",&(cr->e1),"e1/F");
  trr->Branch("br_x1",&(cr->x1),"x1/F");
  trr->Branch("br_x2",&(cr->x2),"x2/F");
  trr->Branch("br_y1",&(cr->y1),"y1/F");
  trr->Branch("br_y2",&(cr->y2),"y2/F");
  trr->Branch("br_n1",&(cr->n1),"n1/I");
  trr->Branch("br_n2",&(cr->n2),"n2/I");
  trr->Branch("br_singles",&(cr->singles),"singles/I");

  trr->Branch("br_NRclust",&(cr->NRclust),"NRclust/I");
  trr->Branch("br_Sel",&(cr->Sel),"Sel/I");
  trr->Branch("br_clE1",&(cr->clE1),"clE1/F");
  trr->Branch("br_clE2",&(cr->clE2),"clE2/F");
  trr->Branch("br_sigmax",&(cr->sigmax),"sigmax/F");
  trr->Branch("br_sigmin",&(cr->sigmin),"sigmin/F");
  trr->Branch("br_chi1",&(cr->chi1),"chi1/F");
  trr->Branch("br_Rank1",&(cr->Rank1),"Rank1/I");  
  trr->Branch("br_Rank2",&(cr->Rank2),"Rank2/I");
  cr->n1=-1;
  cr->n2=-1;
  cr->x1=0.;
  cr->x2=0.;
  cr->y1=0.;
  cr->y2=0.;
  cr->e1=0.;
  cr->er=0.;
  cr->mr=0.;
  cr->clE1=0;
  cr->clE2=0;
  cr->chi1=0;
  cr->sigmax=0;
  cr->sigmin=0;	
  
  Int_t nentries=CellD->GetEntries();
  printf("nentries=%d \n",nentries);
  Int_t skip=1;
  Int_t cntpass=0;
  Float_t x[100][4],y[100][4],e[100][4];
  Yiqun* rc_[4];
  TGraph* hit_[4];
  Int_t N_phot[4];
  Int_t In2[4];
  Int_t Instb2=-1;
  TMatrixT<float>* tm_[4];
  TH2F* hm_[4];
  TMatrixT<float>* tf_[4];
  TH2F* hf_[4];
  for(int jj=0;jj<4;jj++)
    {
      tm_[jj]=0;
      hm_[jj]=0;
      tf_[jj]=0;
      hf_[jj]=0;
      hit_[jj]=0;
      tm_[jj]=new TMatrixT<float>(NROWS[jj]/2,NROWS[jj]);
      hm_[jj]=new TH2F(*tm_[jj]);
      tf_[jj]=new TMatrixT<float>(NROWS[jj]/2,NROWS[jj]);
      hf_[jj]=new TH2F(*tf_[jj]);
      hit_[jj]=new TGraph();
    };

  TLorentzVector phv[100][4];
  TCanvas* c4;
  if(evplot && Nevplot>0 ) 
    {
      printf("starting ev display\n");
      c4=new TCanvas("c4","c4",500,500);
      c4->Print("pltsel.ps(");
    };
  for(int in_stb=1;in_stb<5; in_stb++){rc_[in_stb-1]=0;};
  if(SampleMax==0)SampleMax=nentries;

  CellD->Print();
  for(Int_t j=0;j<nentries;j++)
    {

      if(evplot&&Nevplot>0)
	{
	  c4->Clear();
	  //	c4->Divide(4,2);
	  c4->Divide(2,1);
	};
      for(int in_stb=1;in_stb<5; in_stb++)
	{
	  if(hit_[in_stb-1]!=0){delete hit_[in_stb-1];};
	  if(rc_[in_stb-1]!=0){delete rc_[in_stb-1];};
	  if(tm_[in_stb-1]!=0){delete tm_[in_stb-1];};
	  if(hm_[in_stb-1]!=0){delete hm_[in_stb-1];};
	  if(tf_[in_stb-1]!=0){delete tf_[in_stb-1];};
	  if(hf_[in_stb-1]!=0){delete hf_[in_stb-1];};
	  hit_[in_stb-1]=0;
	  rc_[in_stb-1]=0;
	  tm_[in_stb-1]=0;
	  tf_[in_stb-1]=0;
	  hm_[in_stb-1]=0;
	  hf_[in_stb-1]=0;
	  N_phot[in_stb-1]=0;
	};
      if(cntpass>SampleMax)j=nentries;      
      celldat.ievt=0;
      //      printf("ievt branch address=%x\n",CellD->GetBranch("br_ievt")->GetAddress());
      CellD->GetEntry(j);
      //      printf("j=%d : celldat.ievt address=%x\n",j,&celldat.ievt);
      Int_t tmpe1,tmpe2;
      tmpe1=celldat.ievt/1000;
      tmpe2=celldat.ievt%1000;
      //      printf("Seg =%d evt=%d run %d Mcell=%f Epair=%f (%x)\n",tmpe1,tmpe2,celldat.Rnum,celldat.Mcell,celldat.Epair,celldat.ievt);
      cr->Rank1=-1;
      cr->Rank2=-1;
      if(pRank)cr->Rank1=(Int_t) pRank->GetValue(Iew,Instb,Row1-1,Col1-1);
      Bool_t skip=false;
      if(precut!=0)
	{
	  if(precut->EvalInstance()==0)
	    {
	      skip=true;
	    };
	};
      
      int r2=(int) celldat.Y2;
      int c2=abs((int) celldat.X2);
      int qval;
      if(r2>=0 && r2<NROWS[Instb-1] && c2>=0 && c2<NROWS[Instb-1]/2 && pRank)
	{
	  qval=pRank->GetValue(Iew,Instb,r2,c2);
	};

      if(celldat.nSavedHits<=0)
	{
	  continue;
	  printf("no saved hits \n");
	};
      if(skip)continue;
      cntpass++;

      char str[200];
      if((evplot&&Nevplot>0) || true)
	{
	  sprintf(str,"(%d) Y2=%f X2=%f r2=%d c2=%d q=%x E=%f E1=%f nhits=%d Mcell=%f \n",Nevplot,celldat.Y2,celldat.X2,r2,c2,qval,celldat.Epair,celldat.E1,celldat.nSavedHits,celldat.Mcell);
	  //	  printf("%s \n",str);
	  
	};
      
      //TMatrix Emat=FillFMSEnergy(Instb);
      //	Yiqun rec(&Emat,p_Geom,fpdgain,fpdcorr,2,Instb);
      for(int in_stb=1;in_stb<5; in_stb++)
	{
	  tm_[in_stb-1]=new TMatrixT<float>(FillFMSEnergy(in_stb));
	  if(evplot&&Nevplot>0)
	    {
	      
	      printf("sum(%d)=%f\n",in_stb-1,tm_[in_stb-1]->Sum());
	      //	      c4->cd(2*(in_stb-1)+1);
	      c4->cd(1);
	      hm_[in_stb-1]=new TH2F(*(tm_[in_stb-1]));
	      hm_[in_stb-1]->SetMaximum(50);
	      hm_[in_stb-1]->SetStats(0);
	      if(in_stb==Instb)
		{
		  hm_[in_stb-1]->Draw("zcol");
		  hm_[in_stb-1]->SetMarkerSize(.4);
		  hm_[in_stb-1]->Draw("textsame");
		  hm_[in_stb-1]->SetTitle(str);
		  //	      c4->GetPad(2*(in_stb-1)+1)->SetLogz();
		  c4->GetPad(1)->SetLogz();
		}
	    }; 
	  Bool_t Rec4Dets=false;
	  
	  if(tm_[in_stb-1]->Sum()>.1 && (in_stb==Instb || Rec4Dets))
	    {
	      rc_[in_stb-1] =new Yiqun(tm_[in_stb-1],p_Geom,fpdgain,fpdcorr,2,in_stb);
	      N_phot[in_stb-1]=rc_[in_stb-1]->NPh;
	      if(evplot&&Nevplot>0)
		{
		  //		  c4->cd(2*(in_stb-1)+2);
		  c4->cd(2);
		  tf_[in_stb-1]=new TMatrixT<float>(rc_[in_stb-1]->FittedMat());
		  hf_[in_stb-1]=new TH2F(*(tf_[in_stb-1]));		  
		  hf_[in_stb-1]->SetMaximum(50);
		  hf_[in_stb-1]->SetStats(0);
		  hit_[in_stb-1]=new TGraph();
		  Float_t ene[10];
		  Float_t m12=0;
		  ene[0]=ene[1]=ene[2]=0;
		  for(int hnum=0;hnum<rc_[in_stb-1]->NPh;hnum++)
		    {
		      printf("hit section \n");
		      if(hnum==1)
			{
			  m12=(rc_[in_stb-1]->mom(0)+rc_[in_stb-1]->mom(1)).Mag();
			};
		      Float_t wx=rc_[in_stb-1]->widLG[0];
		      Float_t wy=rc_[in_stb-1]->widLG[1];
		      Float_t xx=rc_[in_stb-1]->photons[hnum].xPos/wx;
		      Float_t yy=rc_[in_stb-1]->photons[hnum].yPos/wy;
		      ene[hnum]=rc_[in_stb-1]->mom(hnum).E();
		      //		      printf("xx=%f yy=%f \n",xx,yy);
		      hit_[in_stb-1]->SetPoint(hnum,xx,yy);
		      hit_[in_stb-1]->Print();
		    };
		  //		  c4->GetPad(2*(in_stb-1)+2)->SetLogz();
		  char figtit[200];
		  sprintf(figtit,"e1=%f e2=%f m12=%f \n",ene[0],ene[1],m12);
		  if(in_stb==Instb)
		    {
		      hf_[in_stb-1]->SetTitle(figtit);
		      hf_[in_stb-1]->Draw("zcol");
		      hf_[in_stb-1]->SetMarkerSize(.4);
		      hf_[in_stb-1]->Draw("textsame");
		      hit_[in_stb-1]->Draw("*");
		      c4->GetPad(2)->SetLogz();
		      printf("draw hit \n");
		      hit_[in_stb-1]->Draw("*");
		      c4->cd(1);
		      hit_[in_stb-1]->Draw("*");
		      printf("draw hit done \n");
		    }
		};
	    };
	};
      if(evplot&&Nevplot>0)
	{
	  c4->Update();
	  if(Nevplot>1)
	    c4->Print("pltsel.ps");
	  else
	    c4->Print("pltsel.ps)");

	};
      Nevplot--;


      Bool_t fiducial=true;
      Yiqun* p_rec=0;

      if(rc_[Instb-1])
	{
	  p_rec=rc_[Instb-1];
	  if(p_rec!=0)
	    {
	      Float_t wx=p_rec->widLG[0];
	      Float_t wy=p_rec->widLG[1];
	      Int_t kcnt=0;

	      for(int in_stb=1;in_stb<5; in_stb++)
		{
		  p_rec=rc_[in_stb-1];
		  if(p_rec!=0)
		    {
		      //printf("rc_[%d]->NPh=%d \n",in_stb-1,rc_[in_stb-1]->NPh);
		      for(int k=0;k<p_rec->NPh;k++)
			{
			  if(k>98)continue;
			  x[k][in_stb-1]=p_rec->photons[k].xPos/wx;
			  y[k][in_stb-1]=p_rec->photons[k].yPos/wy;
			  e[k][in_stb-1]=p_rec->mom(k).E();
			  phv[k][in_stb-1]=p_rec->mom(k);
			  kcnt++;
			  if(cntpass%NskipPrint ==0)printf("ev: %d:",j);
			  if(cntpass%NskipPrint==0)printf("E[%d]=%f)",in_stb-1,e[k][in_stb-1]);
			};
		    };
		};

	      if(kcnt==0)continue;
	      Yiqun* p_rec=rc_[Instb-1];
	      if(cntpass%NskipPrint==0) printf("\n");
	      Int_t FirstClust=-1;
	      Int_t SecondClust=-1;
	      cr->n1=-1;
	      cr->n2=-1;
	      cr->x1=0.;
	      cr->x2=0.;
	      cr->y1=0.;
	      cr->y2=0.;
	      cr->e1=0.;
	      cr->er=0.;
	      cr->mr=0.;
	      cr->clE1=0;
	      cr->clE2=0;
	      cr->chi1=0;
	      cr->sigmax=0;
	      cr->sigmin=0;	
	      Int_t FirstCluster=-1;
	      Int_t SecondCluster=-1;
	      Int_t FirstClusterInstb=-1;
	      Int_t SecondClusterInstb=-1;
	      // require first cluster to be in local Instb

	      if(p_rec->NPh>0)
		{
		  cr->n1=PhNearX(p_rec,celldat.X1*wx,celldat.Y1*wy);

		  if(cr->n1>-1)
		    {
		      cr->x1=x[cr->n1][Instb-1];
		      cr->y1=y[cr->n1][Instb-1];
		      cr->e1=e[cr->n1][Instb-1];
		      cr->er=cr->e1;
		      FirstClust=ClNearX(p_rec,cr->x1*wx,cr->y1*wy);
		      Int_t FirstClusterInstb=Instb;;
		      // printf("FirstClust (Instb=%d) =%d \n",Instb,FirstClust);
		      TVector3 vxy1(cr->x1,cr->y1,0);
		      // printf("Check1\n");
		      fiducial=CheckFiducial(2,Instb,&vxy1,pRank,.5);
		      if(!fiducial)cr->Rank1=cr->Rank1|0x0100000;
		    };  	      
		};


	      //look first in Instb for second cluster
	      Bool_t Keeplooking=true;
	      for(int clc=0;clc<4;clc++)
		{
		  Int_t current_instb=(Instb-1+clc)%4;
		  p_rec=rc_[current_instb];
		 
		  if(N_phot[current_instb]>1 && Keeplooking)
		    {
		      cr->n2=PhNearX(p_rec,celldat.X2*wx,celldat.Y2*wy);
		      Float_t distmiss1=sqrt(pow(cr->x1-celldat.X1,2)+pow(cr->y1-celldat.Y1,2));
		      //		      printf("distmiss1=%f \n",distmiss1);

		      if(cr->n1>-1 &&cr->n2>-1)
			{
			  cr->x2=x[cr->n2][Instb-1];
			  cr->y2=y[cr->n2][Instb-1];
			  cr->er=(phv[cr->n1][Instb-1]+phv[cr->n2][current_instb]).E();;
			  cr->mr=(phv[cr->n1][Instb-1]+phv[cr->n2][current_instb]).Mag();;
			  SecondClust=ClNearX(p_rec,cr->x2*wx,cr->y2*wy);
			  // printf("  n2=%d en2=%f",cr->n2,cr->er-cr->e1);
			  // printf(" M12= %f \n",cr->mr);
			  Float_t distmiss2=sqrt(pow(cr->x2-celldat.X2,2)+pow(cr->y2-celldat.Y2,2));
			  //			  printf("distmiss2=%f\n",distmiss2);
			  if(distmiss2<2.)Keeplooking=false;
			  Instb2=current_instb+1;
			};
		    };
		};
	      if(FirstClust>-1)
		{
		  p_rec=rc_[Instb-1];
		  HitCluster* cl=&(p_rec->clust[FirstClust]);
		  cr-> clE1=cl->energy;
		  cr->chi1=cl->chiSquare;
		  cr-> sigmax=cl->sigmaMax;
		  cr->sigmin=cl->sigmaMin;
		  cr->singles=0;
		};
	      if(SecondClust>-1 && (SecondClust!=FirstClust || Instb2 != Instb))
		{
		  p_rec=rc_[Instb2-1];
		  HitCluster* cl=&(p_rec->clust[SecondClust]);
		  cr->clE2=cl->energy;
		  cr->singles=Instb*10+Instb2*100;
		};
	      //      printf("\n");
	      
	      if((rc_[Instb2-1])>0 && pRank!=0)
		{

		  int ir2=(int) (cr->y2);
		  int ic2=(int) (cr->x2);
		  cr->Rank2=-1;
		  //		  printf("ir2=%d ic2=%d Instb2=%d \n",ir2,ic2,Instb2);
		  fiducial=false;

		  if(ir2>=0 && ir2<=NROWS[Instb2-1]-1 && ic2>=0 && ic2<NROWS[Instb-1]/2 &&N_phot[Instb2-1]>0 )
		    {
		      cr->Rank2=(Int_t) pRank->GetValue(Iew,Instb2,ir2,ic2);
		      TVector3 vxy2(cr->x2,cr->y2,0);
		      fiducial=CheckFiducial(2,Instb2,&vxy2,pRank,.5);
		    };
		  if(!fiducial)cr->Rank2=(cr->Rank2)|(0x0100000);
		};
	      //
	      //	      printf("X1=%f X2=%f Y1=%f Y2=%f \n",celldat.X1,celldat.X2,celldat.Y1,celldat.Y2);
	      //	      printf("x1=%f x2=%f y1=%f y2=%f \n",cr->x1,cr->x2,cr->y1,cr->y2);
	      //	      printf("j=%d Mcell=%f mr=%f \n",j,celldat.Mcell,cr->mr);


	      if(UseQcut)
		{
		  if( (cr->Rank2&0x0100000)==0 && (cr->Rank1&0x0100000)==0 )
		    {
		      //		      printf("Fill after RankCheck \n"); 
		      if(cntpass%NskipPrint ==0)
			{
			  printf("cnt=%d mr=%f er=%f e1=%f\n",
				 cntpass,cr->mr,cr->er,cr->e1);
			};
		      trr->Fill();

		    }
		  else
		    {
		      if(cntpass%NskipPrint ==0)
			{
			  printf("Fail to fill: cnt=%d mr=%f er=%f e1=%f\n",
				 cntpass,cr->mr,cr->er,cr->e1);
			};

		    }
		}
	      else
		{
		  if(cntpass%NskipPrint ==0)
		    {
		      printf("cnt=%d mr=%f er=%f e1=%f\n",
			     cntpass,cr->mr,cr->er,cr->e1);
		    };
		  
		  trr->Fill();
		};
	    };
	}
    };
  trr->Print();

  return trr;
};

Int_t Cell::ClNearX(Yiqun* prec,Float_t x1, Float_t y1)
{
  Int_t n=-1;
  Int_t Nr=prec->NRealClusts;
  Float_t dmin=100.;  
  for(int j=0;j<Nr;j++)
    {
      Float_t x0=prec->clust[j].x0*prec->widLG[0];
      Float_t y0=prec->clust[j].y0*prec->widLG[1];
      Float_t d=(x1-x0)*(x1-x0)+(y1-y0)*(y1-y0);
      if(d<dmin)
	{
	  n=j;
	  dmin=d;
	};
    };
  return n;
};
Int_t Cell::PhNearX(Yiqun* prec,Float_t x1,Float_t y1)
{
  Int_t n=-1;
  Float_t dmin=100000;
  for(int i=0;i<prec->NPh;i++)
    {
      Float_t d1x=prec->photons[i].xPos-x1;
      Float_t d1y=prec->photons[i].yPos-y1;
      Float_t d=d1x*d1x+d1y*d1y;
      if(d<dmin)
	{
	  dmin=d;
	  n=i;
	};      
    };
  return n;
};

Int_t Cell::RankCell(CalibStr* rank)
{
  /*
    0)  Dead Cell (ADC ~ 0)
    1)  Shifted gain
    100) Good Cell  
   */
  Int_t val=-1;
  Int_t nbins=0;

  if(p_adc)
    {
      nbins=p_adc->GetNbinsX();
      if(p_adc->Integral(6,100)<10)val=0;
    };
  if(p_adc && val!=0)
    {
      Float_t EvenOdd[2];
      EvenOdd[0]=EvenOdd[1]=0;
      Int_t nconseczero=0;
      Bool_t zerosend=false;
      for(int j=2;j<nbins;j=j+2)
	{
	  EvenOdd[0]+=p_adc->GetBinContent(j);
	  EvenOdd[1]+=p_adc->GetBinContent(j+1);
	  if(p_adc->GetBinContent(j)==0 &&!zerosend)nconseczero++;
	  if(p_adc->GetBinContent(j)>2)zerosend=true;
	  if(p_adc->GetBinContent(j+1)==0 && !zerosend)nconseczero++;
	  if(p_adc->GetBinContent(j+1)>2)zerosend=true;	  
	}
      if(EvenOdd[0]+EvenOdd[1]>10)
	{
	  Float_t asym=fabs((EvenOdd[0]-EvenOdd[1])/(EvenOdd[0]+EvenOdd[1]));
	  Float_t dasm=1/sqrt(EvenOdd[0]+EvenOdd[1]);
	  if(asym>.9 && dasm<.2)
	    {
	      val=1;
	      if(nconseczero>32)nconseczero=32;
	      if(nconseczero>1)val=nconseczero;
	    };
	}
    };
  if(val<0)val=0x80;
  if(rank)
    {
      rank->SetValue(Iew,Instb,Row1-1,Col1-1,val);
    };
  return val;
};
Bool_t Cell::CheckFiducial(int iew0,int instb0,TVector3* vlocal, CalibStr* rank,Float_t cellfract)
{
  Int_t NROWS[4]={34,34,24,24};
  int EdgeCode=0x00000100;
  int Dead5Code=0x0001000;
  int Dead9Code=0x0010000;
  Int_t r1=vlocal->Y()+1;
  Int_t c1=vlocal->X()+1;
  //  printf("vlocal->X()=%f vlocal->Y()=%f r1=%d c1=%d \n",vlocal->X(),vlocal->Y(),r1,c1);
  if(r1<1 || r1>NROWS[instb0] || c1<1 || c1>NROWS[instb0]/2)return false;
  int thisrank=rank->GetValue(iew0,instb0,r1-1,c1-1);
  if(thisrank==-1)thisrank=0;
  if(((0x08F)&thisrank)==0)return false;
  if(!((EdgeCode | Dead5Code | Dead9Code)&thisrank))
    {
      return true;
    }
  else
    {
      for(int r=r1-1;r<r1+2;r++)
	{
	  for(int c=c1-1;c<c1+2;c++)
	    {
	      Bool_t checkthis=false;
	      if(c==0 || c==NROWS[instb0]/2+1 ||r==0 ||r==NROWS[instb0]+1)
		{
		  checkthis=true;
		  //		  printf("edge\n");
		}
	      else if(c<0 || c>NROWS[instb0]/2+1 ||r<0 ||r>NROWS[instb0]+1)
		{
		  int nrank=rank->GetValue(iew0,instb0,r-1,c-1);
		  if(nrank<0)nrank=0;
		  if((nrank&0x08F)==0)
		    {
		      checkthis=true;
		      //		      printf("r=%d ,c=%d  rankcheck rank=%x \n",r,c,nrank);
		    };
		};
	      if(checkthis)
		{
		  Float_t dist=.5+cellfract;
		  Float_t distx=fabs(vlocal->X()-(c-.5));
		  Float_t disty=fabs(vlocal->Y()-(r-.5));
		  //		  printf(" r=%d c=%d distx=%f disty=%f \n",r,c,distx,disty);
		  if((distx<dist*1.01) && (disty<dist*1.01))
		    {
		      //		      printf("Fail Vx=%f Vy=%f \n",vlocal->X(),vlocal->Y());
		      return false;
		    };
		};
	    };
	};
    };
  return true;
}
Int_t Cell::Quality(CalibStr* rank)
{
  //  Large Cell Edges need work
  Int_t NROWS[4]={34,34,24,24};
  Int_t HOLE[4]={8,8,5,5};
  Int_t q=0;
  if(!rank)return q;
  if(Iew!=2 || Instb<1  || Instb>4)return q;

  Bool_t Edge=false;
  Int_t nRows=NROWS[Instb-1];
  Int_t nCols=nRows/2;
  if(Row1<2 || Col1<2 ||Row1>(nRows-1) ||Col1>(nCols-1))Edge=true;
  if(fabs(-.5+Row1-NROWS[Instb-1]/2)<HOLE[Instb-1]+1 &&
     fabs(-.5+Col1)<HOLE[Instb-1]+1 )Edge=true;
  Bool_t Dead5=false;
  Bool_t Dead9=false;
  if(!Edge)
    {
      for(int ir=Row1-1;ir<Row1+2;ir++)
	{
	  for(int ic=Col1-1;ic<Col1+2;ic++)
	    {
	      Int_t v1=rank->GetValue(Iew,Instb,ir-1,ic-1);
	      Int_t rankmask= 0x000000FF & v1;
	      if(rankmask<1)
		{
		  Dead9=true;
		  if((ir-Row1)*(ic-Col1)==0)Dead5==true;
		}
	    }      
	}
    }
  else {Dead5=Dead9=true;};
  q= 0x000000FF & ((int) (rank->GetValue(Iew,Instb,Row1-1,Col1-1)));
  if(Edge)q=q+0x00000100;
  if(Dead5)q=q+0x0001000;
  if(Dead9)q=q+0x0010000;
  rank->SetValue(Iew,Instb,Row1-1,Col1-1,q);
  return q;
};
