#include "AnalTools.h"
#include <string>
using namespace std;

static TH2F* histFit_EPt;
static TH2F* histFit_EPtF;
static Int_t rowLim1,rowLim2,colLim1,colLim2,colMask[49];
static Geom* PGEOM;
static Int_t FitOption;

extern Double_t FunEPt_(Double_t *x,Double_t *par)
{ 
  Double_t pt,ptnom;
  Double_t xpos;
  Double_t xf;
  Double_t powx= par[1];
  Double_t powpt= par[2];
  Double_t Dphi=0;
  Double_t Theta=0;
  Double_t ebin=x[1];
  Double_t ch=x[0];
  xf=(ebin+.05)/100.;
  Int_t bin=(Int_t) ch;
  Int_t col=((Int_t) (bin))%7;
  Int_t row=((Int_t) (bin)/7);
  TVector3 vv((col+.5)*3.8,(row+.5)*3.8,0.);
  pt=sin(Theta=(PGEOM->GlobalXYZ(1,1,vv).Theta()))*xf*100.;
  if(Theta!=0)Dphi=1/Theta;
  ptnom=50*(20.4+3.8)/809;  

  powpt=powpt+par[3]*pt/ptnom+
    par[4]*TMath::Power(pt/ptnom,2)+
    par[5]*TMath::Power(pt/ptnom,3);

  Double_t fitval=par[0]
    *TMath::Power(1-xf,powx)*  TMath::Power(pt/ptnom,-powpt )*Dphi;
  return fitval; 
};

extern void FcnFitEPt_(Int_t & npar, Double_t *grad,  Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t pt,ptnom;
  Double_t xpos;
  Double_t ypos;
  Double_t xf;
  Double_t powx= par[1];
  Double_t powpt= par[2];
  Double_t Theta;
  Double_t De;
  Double_t Dpt;
  Double_t Dphi;
  //  printf("npar=%d par[0]=%f par[2]=%f\n",npar,par[0],par[2]);
  f=0;
  for(Int_t ch=colLim1;ch<colLim2+1;ch++)
  {
    for(Int_t ebin=rowLim1;ebin<rowLim2;ebin++)
    {

      //	  xf=(ebin-.5)/100.;
      xf=(ebin+.5)/100.;
      Int_t bin=(Int_t) ch;
      Int_t col=((Int_t) (bin))%7;
      Int_t row=((Int_t) (bin)/7);
      //	  if(col==0 || col==6)continue;
      //	  TVector3 vv((col+.5)*3.8,(row+.5)*3.8,0.);
      TVector3 vv((col+.1)*3.8,(row+.5)*3.8,0.);
      TVector3 vlow((col)*3.8,(row+.5)*3.8,0.);
      TVector3 vhi((col+1)*3.8,(row+.5)*3.8,0.);
      Double_t dY=fabs(
          PGEOM->GlobalXYZ(1,1,vlow).PseudoRapidity()-
          PGEOM->GlobalXYZ(1,1,vhi).PseudoRapidity());

      TVector3 gv=PGEOM->GlobalXYZ(1,1,vv);
      xpos=gv.X();
      ypos=gv.Y();
      pt=sin(Theta=(gv.Theta()))*xf*100.;
      ptnom=50*(20.4+2*3.8)/809;  
      Dphi=0;
      Dpt=sqrt(xpos*xpos+ypos*ypos)/gv.Z();
      if(Theta!=0)Dphi=fabs((xpos)/(xpos*xpos+ypos*ypos));
      De=1;
      Double_t Jacobian=xf;
      Double_t fitval=0;
      //fabs(xpos/(xpos*xpos+ypos*ypos)*pt);
      if(FitOption==1)
      {
        fitval=par[0]*exp(-par[1]*(ebin+.5-50));
      }
      else
      { 

        powpt=powpt+par[3]*pt/ptnom+
          par[4]*TMath::Power(pt/ptnom,2)+
          par[5]*TMath::Power(pt/ptnom,3);


        fitval=par[0]
          *TMath::Power(1-xf,powx)*  TMath::Power(pt/ptnom,-powpt )*Jacobian;
      };
      Int_t binn= histFit_EPtF->FindBin(ch+.5,ebin+.5);
      Double_t hval=histFit_EPt->GetBinContent(binn);
      Double_t err=histFit_EPt->GetBinError(binn);
      if(err<.1*hval)err=.1*hval;
      if(err<1.)err=1.;

      if(histFit_EPtF){
        histFit_EPtF->SetBinContent(binn,fitval);
      }
      //only increment f for channels with Mask set >0

      if(err>0 && colMask[ch]>0)f=f+(fitval-hval)*(fitval-hval)/err/err;
      if(iflag==5) printf("col=%d row=%d hval=%f fitval=%f f=%f\n",col,row,hval,fitval,f);

    };
  };

};

ClassImp(AnalTools)

AnalTools::AnalTools()
{
  pbbc=new PullBBC();
  prp=new PullRP();
  for(int j=0;j<4;j++)Eminpout[4]=6.;
  Enable_ReRec=false;
  Mgr=0;
  MaxSavedHits=0;
  RunDepBaseOverride=0;
  UseRunDepCor=false;
  Rdep=0;
  currentdset=0;
  TrigBitsSelect=0;
  TrigBitsVeto=0;
  OnlyLED=false;
  allLED=false;
  File_adcTr="adcTr.root";
  IdentifyLED=true;
  ERecMin[0]=12;
  ERecMin[1]=12;
  ERecMin[2]=20;
  ERecMin[3]=20;
  p_hist1=0;
  p_hist2=0;
  OutFileName="Output.root";
  NumberEventsToProcess=0;
  OutputToSingle=false;
  OutputFileName="tmpdata/OFile.root";
  NumberLCell4Vectors=0;
  p_Geom=0;
  Float_t OutOf10=1.;
  FcnFitEPt=FcnFitEPt_;
  FunEPt=FunEPt_;
  p_colMask=&(colMask[0]);
  SetFitOption();
  ReScaledHist=0;
  WriteADC=false;
};
AnalTools::~AnalTools()
{ 
  for(Int_t j=0;j<NumberLCell4Vectors;j++)
  {
    delete p_LCell4Vec[j];
  }; 
};
Bool_t AnalTools::SetFitOption(Int_t iopt)
{
  FitOption=iopt;
};
Int_t AnalTools::GetFitOption()
{
  return FitOption;
};
TH1F AnalTools::GetAsm(TH1F* h1,TH1F* h2,char* name)
{
  h1->Sumw2();
  h2->Sumw2();
  TH1F newasm( (*h1-*h2)/(*h1+*h2));
  newasm.SetName(name);
  char str[100];
  sprintf(str,"(%s -%s)/(%s+%s)",h1->GetName(),h2->GetName(),h1->GetName(),h2->GetName());
  newasm.SetTitle(str);
  return newasm;
};
TH1D AnalTools::GetRatio(Int_t ih,TH2F* h20, char* hname)
{
  TH1D* h[4];
  Int_t n0=(ih-1)/4;
  Int_t nsel=ih-4*n0-1;
  Int_t j;
  char* nms[4]={"tmp1","tmp2","tmp3","tmp4"};
  for(j=0;j<4;j++)
  {
    Int_t k=4*n0+1+j;
    h[j]=h20->ProjectionX(nms[j],k,k);
    h[j]->Sumw2();
    std::cout<<h[j]->GetName()<<j<<":"<<h[j]->Integral()<<"\n";
  };
  p_histogram=new TH1D( (*h[nsel])/(*h[0]+*h[1]+*h[2]+*h[3]));
  for(j=0;j<4;j++)delete h[j];
  p_histogram->SetName(hname);
  return *p_histogram;
};
TH1D AnalTools::GetAsm(Int_t jup, Int_t jdn, TH2F* h20,char* name)
{

  TH1D h1=(TH1D) *(h20->ProjectionX("tmp1",jup,jup));
  TH1D h2=(TH1D) *(h20->ProjectionX("tmp2",jdn,jdn));
  h1.Sumw2();
  h2.Sumw2();
  TH1D newdif( (TH1D) (h1-h2));
  TH1D newsum( (TH1D) (h1+h2));
  TH1D newasm(newdif/newsum);
  newasm.SetName(name);
  char s1[100];
  char s2[100];
  char str[100];
  sprintf(s1,"%s%d",h20->GetName(),jup);
  sprintf(s2,"%s%d",h20->GetName(),jdn);
  sprintf(str,"(%s -%s)/(%s+%s)",s1,s2,s1,s2);
  newasm.SetTitle(str);
  return newasm;
};
TH2F AnalTools::GetAsm(TH2F* h1,TH2F* h2,char* name)
{
  h1->Sumw2();
  h2->Sumw2();
  TH2F newdif( (TH2F) (*h1-*h2));
  TH2F newsum( (TH2F) (*h1+*h2));
  TH2F newasm(newdif/newsum);
  newasm.SetName(name);
  char str[100];
  sprintf(str,"(%s -%s)/(%s+%s)",h1->GetName(),h2->GetName(),h1->GetName(),h2->GetName());
  newasm.SetTitle(str);
  return newasm;
};


TH2F AnalTools::GetCross(TH2F* h1u,TH2F* h1d,TH2F* h2u, TH2F* h2d,char* name,Float_t ermin)
{
  Int_t nx=h1u->GetNbinsX();
  Int_t ny=h1u->GetNbinsY();
  Bool_t bad=false;
  if(nx != h1d->GetNbinsX())bad=true;
  if(ny != h1d->GetNbinsY())bad=true;
  if(nx != h2u->GetNbinsX())bad=true;
  if(ny != h2u->GetNbinsY())bad=true;
  if(nx != h2d->GetNbinsX())bad=true;
  if(ny != h2d->GetNbinsY())bad=true;
  Float_t lowx=h1u->GetXaxis()->GetXmin();
  Float_t highx=h1u->GetXaxis()->GetXmax();
  Float_t lowy=h1u->GetYaxis()->GetXmin();
  Float_t highy=h1u->GetYaxis()->GetXmax();

  if(bad){
    printf("Error: unlike histograms\n");
    return *h1u;
  };

  Float_t nuS=0;
  Float_t nuN=0;
  Float_t ndS=0;
  Float_t ndN=0;

  TH2F nasm(name,name,nx,lowx,highx,ny,lowy,highy);

  for(Int_t i=1;i<nx+1;i++){
    for(Int_t j=1;j<ny+1;j++){
      nuS=TMath::Max(h1u->GetBinContent(i,j),0.);
      ndS=TMath::Max(h1d->GetBinContent(i,j),0.);
      nuN=TMath::Max(h2u->GetBinContent(i,j),0.);
      ndN=TMath::Max(h2d->GetBinContent(i,j),0.);
      Float_t inverr2=nuS+ndS+nuN+ndN;
      Float_t value=0;
      Float_t err=0;
      Float_t denom=sqrt(nuS*ndN)+sqrt(ndS*nuN);
      if(denom>0 && inverr2>0){
        value=(sqrt(nuS*ndN)-sqrt(ndS*nuN))/denom;
        err=1/sqrt(inverr2);
        if(err>ermin)value=0;
      };
      nasm.SetBinContent(i,j,value);
      nasm.SetBinError(i,j,err);
    };
  };
  return nasm;
}
TH1F AnalTools::GetCross(Int_t ul, Int_t dl,Int_t ur,Int_t dr,TH2F* h20,char* name,Int_t secondoff)
{ 
  Int_t dc;
  TH1D hul=(TH1D) *(h20->ProjectionX("tmp1",ul,ul));
  TH1D hdl=(TH1D) *(h20->ProjectionX("tmp2",dl,dl));
  TH1D hur=(TH1D) *(h20->ProjectionX("tmp3",ur,ur));
  TH1D hdr=(TH1D) *(h20->ProjectionX("tmp4",dr,dr));
  if((dc=secondoff)>0)
  {
    hul=hul+(TH1D) *(h20->ProjectionX("tmp1a",ul+dc,ul+dc));
    hdl=hdl+(TH1D) *(h20->ProjectionX("tmp2a",dl+dc,dl+dc));
    hur=hur+(TH1D) *(h20->ProjectionX("tmp3a",ur+dc,ur+dc));
    hdr=hdr+(TH1D) *(h20->ProjectionX("tmp4a",dr+dc,dr+dc));
  }; 
  hul.Sumw2();
  hdl.Sumw2();
  hur.Sumw2();
  hdr.Sumw2();
  printf("hul=%f hdl=%f hur=%f hdr=%f \n",hul.Integral(),hdl.Integral(),hur.Integral(),hdr.Integral());
  Double_t xmin=  hul.GetXaxis()->GetXmin();
  Double_t xmax=  hul.GetXaxis()->GetXmax();
  Int_t nbins=hul.GetNbinsX();
  TH1F newasm(name,name,nbins,xmin,xmax);
  for(Int_t j=1;j<nbins+1;j++)
  {
    Float_t nul=hul.GetBinContent(j);
    Float_t ndl=hdl.GetBinContent(j);
    Float_t nur=hur.GetBinContent(j);
    Float_t ndr=hdr.GetBinContent(j);
    Float_t err=0;
    if(nul+ndl+nur+ndr>0)err=1./sqrt(nul+ndl+nur+ndr);
    Float_t v1=sqrt(nul*ndr);
    Float_t v2=sqrt(ndl*nur);
    Float_t val=0.;
    if(v1>0 && (v2>0))val=(v1-v2)/(v1+v2);
    newasm.SetBinContent(j,val);
    newasm.SetBinError(j,err);
  };
  return newasm;
}
TH1F AnalTools::GetCross(TH1F* h1u,TH1F* h1d,TH1F* h2u, TH1F* h2d,char* name,Float_t ermin)
{
  Int_t nx=h1u->GetNbinsX();
  Bool_t bad=false;
  if(nx != h1d->GetNbinsX())bad=true;
  if(nx != h2u->GetNbinsX())bad=true;
  if(nx != h2d->GetNbinsX())bad=true;
  Float_t lowx=h1u->GetXaxis()->GetXmin();
  Float_t highx=h1u->GetXaxis()->GetXmax();

  if(bad){
    printf("Error: unlike histograms\n");
    return *h1u;
  };

  Float_t nuS=0;
  Float_t nuN=0;
  Float_t ndS=0;
  Float_t ndN=0;

  TH1F nasm(name,name,nx,lowx,highx);

  for(Int_t i=1;i<nx+1;i++){
    nuS=TMath::Max(h1u->GetBinContent(i),0.);
    ndS=TMath::Max(h1d->GetBinContent(i),0.);
    nuN=TMath::Max(h2u->GetBinContent(i),0.);
    ndN=TMath::Max(h2d->GetBinContent(i),0.);
    Float_t inverr2=nuS+ndS+nuN+ndN;
    Float_t value=0;
    Float_t denom=sqrt(nuS*ndN)+sqrt(ndS*nuN);
    Float_t err=0;

    if(inverr2>0 && denom>0){
      value=(sqrt(nuS*ndN)-sqrt(ndS*nuN))/denom;
      err=1/sqrt(inverr2);
    };
    nasm.SetBinContent(i,value);
    nasm.SetBinError(i,err);
  };
  return nasm;
}


TH1F AnalTools::HistShift(TH1F h0,Float_t dbin)
{
  return h0;
};

Int_t AnalTools::NTowersAbove(TMatrixT<float>* tm, Float_t v)
{
  Int_t ne =tm->GetNoElements();
  Float_t* fe=tm->GetMatrixArray();
  Int_t cnt=0;
  for(Int_t j=0;j<ne;j++){if(fe[j]>v)cnt++;};
  return cnt;
};
TVector3 AnalTools::mom(TMatrixT<float>* tm,Geom* p_geom,Int_t EW, Int_t NSTB)
{
  return v3=mom( tm,-1000,1000,-1000,1000, p_geom,EW,NSTB);
}; 

/* TLorentzVector AnalTools::FourMom(TMatrixT<float>* tm,Geom* p_geom,Int_t EW, Int_t NSTB)
   {
   return v4=FourMom( tm,-1000,1000,-1000,1000, p_geom,EW,NSTB);
   };
   */
TLorentzVector AnalTools::FourMom(TMatrixT<float>* tm,Int_t rowlo ,Int_t rowhi,Int_t collo,Int_t colhi,Geom* p_geom,Int_t EW,Int_t NSTB)
{
  TVector3 mtemp= mom(tm,rowlo,rowhi,collo,colhi,p_geom,EW,NSTB);
  TLorentzVector lv(mtemp,EnergySum);
  v4=lv;
  return v4;
};

Int_t AnalTools::SegmentSearch(TMatrixT<float>* tm,Geom* p_geom,Int_t EW,Int_t NSTB)
{
  Float_t emin=.25;
  Int_t lowcol=0;
  Int_t highcol=-1000;
  Int_t nrows=tm->GetNrows();
  Int_t ncols=tm->GetNcols();
  TLorentzVector vcol;
  for(Int_t j=0;j<NumberLCell4Vectors;j++)
  {
    delete p_LCell4Vec[j];
  };
  Bool_t energyfound=false;
  for(Int_t jc=0;jc<ncols;jc++)
  {
    vcol= FourMom(tm,0,nrows-1,jc,jc,p_geom,EW,NSTB);
    if(vcol.E()<emin)
    {
      if(energyfound==true)
      {
        v4= FourMom(tm,0,nrows-1,lowcol,jc,p_geom,EW,NSTB);
        if(v4.E()>2*emin)
        {
          p_LCell4Vec[NumberLCell4Vectors]=new TLorentzVector(v4);
          if(NumberLCell4Vectors<9)NumberLCell4Vectors++;
          lowcol=jc+1;
          energyfound=false;
        };
      };	
      lowcol=jc;  
    }
    else
    {
      energyfound=true;
    };
  };
  return NumberLCell4Vectors;
};

TVector3 AnalTools::mom(TMatrixT<float>* tm,Int_t rowlo ,Int_t rowhi,Int_t collo,Int_t colhi,Geom* p_geom,Int_t EW,Int_t NSTB)
{
  TVector3 off(0,0,0);
  Float_t width=*(p_geom->FpdTowWid(EW,NSTB));
  off[2]=*(p_geom->ZFPD(EW,NSTB));
  off[1]=*(p_geom->yOffset(EW,NSTB));
  off[0]=(*(p_geom->xOffset(EW,NSTB)));
  Int_t signx=-1;
  if(NSTB==1 || (EW==2 && NSTB==3))
  {
    signx=1;
  };
  TVector3 momsum(0.,0.,0.);
  if(rowlo<=tm->GetRowLwb())rowlo=tm->GetRowLwb();
  if(rowhi>=tm->GetRowUpb())rowhi=tm->GetRowUpb();

  if(colhi>=tm->GetColUpb())colhi=tm->GetColUpb();
  TVector3 xyz(0.,0.,0.);
  EnergySum=0;
  Float_t Energy;
  for(Int_t ir=rowlo;ir<=rowhi;ir++)
  {
    for(Int_t ic=collo;ic<=colhi;ic++)
    {
      xyz[0]=off[0]-signx*ic*width;
      xyz[1]=off[1]-ir*width;
      xyz[2]=off[2];
      xyz=xyz*(1/xyz.Mag());
      Energy=TMath::Max((*tm)(ir,ic),(Float_t) 0.);
      xyz=xyz*Energy;
      EnergySum+=Energy;
      momsum=momsum+xyz;
    };
  };
  v3=momsum;
  return v3;   
};

Int_t AnalTools::PeakZeroDisplace(TH1D* h)
{
  Int_t zerobin=h->FindBin(0.);
  Float_t max=-1000.;
  Int_t maxloc=1;
  Int_t nbins=h->GetNbinsX();
  for(Int_t j=1;j<nbins+1;j++)
  {
    Float_t bin=h->GetBinContent(j);
    if(bin>max)
    {
      max=bin;
      maxloc=j;
    };
  };
  return maxloc-zerobin;
};
TH2F* AnalTools::Mat2NewHist2(TMatrixT<float>* tm,char* hname)
{
  Int_t nhx,nhy;
  nhx=tm->GetNcols();
  nhy=tm->GetNrows();
  TH2F* hist=new TH2F(hname,hname,nhx,0,nhx,nhy,0,nhy);
  hist->SetStats(0);
  Int_t ix,iy;
  for(iy=0;iy<nhy;iy++)
  {
    for(ix=0;ix<nhx;ix++)
    {
      Float_t value=(*tm)(iy,ix);
      hist->SetBinContent(ix+1,iy+1,value);
    };
  };
  return hist;
};
TMatrixT<float> AnalTools::MultMat(TMatrixT<float>* tm1,TMatrixT<float>* tm2)
{
  Int_t nhx1,nhy1,nhx2,nhy2;
  nhx1=tm1->GetNcols();
  nhy1=tm1->GetNrows();
  nhx2=tm2->GetNcols();
  nhy2=tm2->GetNrows();
  TMatrixT<float> nmat(nhy1,nhx1);
  nmat=0.;
  if(nhx1!=nhx2 || nhy1!=nhy2)
  {
    return nmat;
  };

  Int_t ix,iy;
  for(iy=0;iy<nhy1;iy++)
  {
    for(ix=0;ix<nhx1;ix++)
    {
      Float_t value=(*tm1)(iy,ix)*(*tm2)(iy,ix);
      nmat[iy][ix]=value;

    };
  };
  return nmat;

};
void AnalTools::CopyMatrix(TMatrixT<float>* ptm,TH2F* ph)
{
  Int_t nhx,nhy;
  Double_t value;
  nhx=ph->GetNbinsX();
  nhy=ph->GetNbinsY();
  for(Int_t ix=0;ix<nhx;ix++)
  {
    for(Int_t iy=0;iy<nhy;iy++)
    {
      value=(*ptm)(iy,ix);
      ph->SetBinContent(ix+1,iy+1,value);
    };
  };
};

Int_t AnalTools::readO(char* filelist,Int_t uu, Int_t ud, Int_t du, Int_t dd)
{
  Int_t DetMap[9]={0,1,2,0,0,0,0,3,4};
  Int_t spinindex=0;
  TChain* p_OF=new TChain("p_out");
  FILE* pf;
  Int_t knt=0;
  Int_t FileCnt=0;

  if(pf=fopen(filelist,"r"))
  {
    char fname[200];
    while(!feof(pf))
    {
      if(fscanf(pf,"%s\n",fname)>0)
      {
        knt++;
        FileCnt++;
        p_OF->Add(fname);
      };
      if(knt>39)
      {
        knt=0;
        printf("%d files added to TChain\n",FileCnt);
      };
      // if(FileCnt>20)break;// max number of segs for testing
    };
  };
  if(OutputToSingle) {
    p_OF->Merge( OutputFileName);
    return 0;
  };
  Float_t lowlim=0;
  Float_t highlim=120.;
  TFile* out=new TFile(OutFileName,"recreate");
  TH2F Ex("Ex","Ex",24,lowlim,highlim,20,0,20);
  TH2F Em1("Em1","Em1",10,0.,1.,20,0,20);
  TH2F Em2("Em2","Em2",40,0.,1.,20,0,20);
  TH2F Emm1("Emm1","Emm1",10,0,1.,20,0,20);
  TH2F Emm2("Emm2","Emm2",40,0.,1.,20,0,20);

  TH1F WSup("WSup","WSup",24,lowlim,highlim);
  TH1F WSdn("WSdn","WSdn",24,lowlim,highlim);
  TH1F WNup("WNup","WNup",24,lowlim,highlim);
  TH1F WNdn("WNdn","WNdn",24,lowlim,highlim);
  TH1F ESup("ESup","ESup",24,lowlim,highlim);
  TH1F ESdn("ESdn","ESdn",24,lowlim,highlim);
  TH1F ENup("ENup","ENup",24,lowlim,highlim);
  TH1F ENdn("ENdn","ENdn",24,lowlim,highlim);

  TH2F WSup2("WSup2","WSup2",12,lowlim,highlim,40,0,10.);
  TH2F WSdn2("WSdn2","WSdn2",12,lowlim,highlim,40,0,10);
  TH2F WNup2("WNup2","WNup2",12,lowlim,highlim,40,0,10);
  TH2F WNdn2("WNdn2","WNdn2",12,lowlim,highlim,40,0,10);
  TH2F ESup2("ESup2","ESup2",12,lowlim,highlim,40,0,10);
  TH2F ESdn2("ESdn2","ESdn2",12,lowlim,highlim,40,0,10);
  TH2F ENup2("ENup2","ENup2",12,lowlim,highlim,40,0,10);
  TH2F ENdn2("ENdn2","ENdn2",12,lowlim,highlim,40,0,10);


  TH2F WSup3("WSup3","WSup3",6,lowlim,highlim,20,0,10);
  TH2F WSdn3("WSdn3","WSdn3",6,lowlim,highlim,20,0,10);
  TH2F WNup3("WNup3","WNup3",6,lowlim,highlim,20,0,10);
  TH2F WNdn3("WNdn3","WNdn3",6,lowlim,highlim,20,0,10);
  TH2F ESup3("ESup3","ESup3",6,lowlim,highlim,20,0,10);
  TH2F ESdn3("ESdn3","ESdn3",6,lowlim,highlim,20,0,10);
  TH2F ENup3("ENup3","ENup3",6,lowlim,highlim,20,0,10);
  TH2F ENdn3("ENdn3","ENdn3",6,lowlim,highlim,20,0,10);


  TH1F WSup4("WSup4","WSup4",24,lowlim,highlim);
  TH1F WSdn4("WSdn4","WSdn4",24,lowlim,highlim);
  TH1F WNup4("WNup4","WNup4",24,lowlim,highlim);
  TH1F WNdn4("WNdn4","WNdn4",24,lowlim,highlim);
  TH1F ESup4("ESup4","ESup4",24,lowlim,highlim);
  TH1F ESdn4("ESdn4","ESdn4",24,lowlim,highlim);
  TH1F ENup4("ENup4","ENup4",24,lowlim,highlim);
  TH1F ENdn4("ENdn4","ENdn4",24,lowlim,highlim);
  lowlim=0.;
  highlim=2.;
  TH1F WSupm("WSupm","WSupm",20,lowlim,highlim);
  TH1F WSdnm("WSdnm","WSdnm",20,lowlim,highlim);
  TH1F WNupm("WNupm","WNupm",20,lowlim,highlim);
  TH1F WNdnm("WNdnm","WNdnm",20,lowlim,highlim);
  TH1F ESupm("ESupm","ESupm",20,lowlim,highlim);
  TH1F ESdnm("ESdnm","ESdnm",20,lowlim,highlim);
  TH1F ENupm("ENupm","ENupm",20,lowlim,highlim);
  TH1F ENdnm("ENdnm","ENdnm",20,lowlim,highlim);

  TH1F WSupmm("WSupmm","WSupmm",20,lowlim,highlim);
  TH1F WSdnmm("WSdnmm","WSdnmm",20,lowlim,highlim);
  TH1F WNupmm("WNupmm","WNupmm",20,lowlim,highlim);
  TH1F WNdnmm("WNdnmm","WNdnmm",20,lowlim,highlim);
  TH1F ESupmm("ESupmm","ESupmm",20,lowlim,highlim);
  TH1F ESdnmm("ESdnmm","ESdnmm",20,lowlim,highlim);
  TH1F ENupmm("ENupmm","ENupmm",20,lowlim,highlim);
  TH1F ENdnmm("ENdnmm","ENdnmm",20,lowlim,highlim); 
  lowlim=0;
  highlim=120.;


  TVector3 v3;
  TLorentzVector v4;
  TH1F  Nphots("Nphots","Nphots",20,0,20);
  TH1F mass("mass","mass",100,0,1.);
  TH1F massU("massU","massU",100,0,1.);

  Int_t tpes[40];
  Int_t spin;
  Float_t pxyzt[40];
  Int_t nwrds;
  Int_t nphotons;

  p_OF->SetBranchAddress("spin",&spin);
  p_OF->SetBranchAddress("nphotons",&nphotons);
  p_OF->SetBranchAddress("br_nwrds",&nwrds);
  p_OF->SetBranchAddress("br_types",tpes);
  p_OF->SetBranchAddress("br_pxyzt",pxyzt);

  TLorentzVector lv[10];
  TLorentzVector ph[10];

  Int_t nentries=p_OF->GetEntries();
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;
  printf("GetEntries=%d \n",nentries);
  knt=0;
  Bool_t starting=true;
  Float_t m0=.135;
  Float_t dm=.05;
  Float_t Enom=40.;
  Float_t dEnom=40;
  for(Int_t i=0;i<nentries;i++)
  {

    Int_t nbytes=p_OF->GetEntry(i);
    knt++;

    if(knt>100000 || (starting && knt>10000)){
      printf("Evt: %d \n",i);
      knt=1;
      if(i>100000)starting=false;
    };

    for(Int_t k=0;k<10;k++)
    {
      lv[k]=TLorentzVector(0.,0.,0.,0.);
      ph[k]=TLorentzVector(0.,0.,0.,0.);
    };
    Int_t kk=0;
    Int_t nwords=nwrds;
    for(Int_t k=0;k<nwords;k+=4){
      if(tpes[k]<9){
        kk=tpes[k];
        lv[kk].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
        //       printf(" non-photon : tpes[%d]=%d mass=%f\n",k,tpes[k],lv[kk].Mag()); 
        massU.Fill(lv[kk].Mag());
        Emm1.Fill(lv[kk].Mag(),spinindex);
        Emm2.Fill(lv[kk].Mag(),spinindex);
      }
      else if((tpes[k]>100) && (tpes[k]<109)){
        Int_t kk=tpes[k]-100;
        Int_t spn=0;
        if(spin>=0 && spin<4)spn=spin;
        spinindex=DetMap[kk]*4+spn;
        ph[kk].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
        //       printf(" photon : tpes[%d]=%d mass=%f\n",k,tpes[k],ph[kk].Mag()); 
        if(ph[kk].E()>65.)
        {
          mass.Fill(ph[kk].Mag());
          Em1.Fill(ph[kk].Mag(),spinindex);
          Em2.Fill(ph[kk].Mag(),spinindex);
        };
        if( (ph[kk].Mag()>.085)&&(ph[kk].Mag()<.185)){
          Ex.Fill(ph[kk].E(),spinindex);
        };
      };
    };
    Nphots.Fill(nphotons);
    if(nphotons==1)
    {
      if(spin==uu || spin==ud){
        if(fabs(lv[8].E())>0)WSdn2.Fill(lv[8].E(),fabs(lv[8].Px()));
        if(fabs(lv[7].E())>0)WNdn2.Fill(lv[7].E(),fabs(lv[7].Px()));
        if(fabs(lv[2].E())>0)ESdn2.Fill(lv[2].E(),fabs(lv[2].Px()));
        if(fabs(lv[1].E())>0)ENdn2.Fill(lv[1].E(),fabs(lv[1].Px()));
      };

      if(spin==du || spin==dd){
        if(fabs(lv[8].E())>0)WSup2.Fill(lv[8].E(),fabs(lv[8].Px()));
        if(fabs(lv[7].E())>0)WNup2.Fill(lv[7].E(),fabs(lv[7].Px()));
        if(fabs(lv[2].E())>0)ESup2.Fill(lv[2].E(),fabs(lv[2].Px()));
        if(fabs(lv[1].E())>0)ENup2.Fill(lv[1].E(),fabs(lv[1].Px()));
      };
      Float_t ptm=0.;
      if(spin==uu || spin==ud){
        if(fabs(lv[8].Px())>ptm)WSdn.Fill(lv[8].E());
        if(fabs(lv[7].Px())>ptm)WNdn.Fill(lv[7].E());
        if(fabs(lv[2].Px())>ptm)ESdn.Fill(lv[2].E());
        if(fabs(lv[1].Px())>ptm)ENdn.Fill(lv[1].E());
      };

      if(spin==du || spin==dd){
        if(fabs(lv[8].Px())>ptm)WSup.Fill(lv[8].E());
        if(fabs(lv[7].Px())>ptm)WNup.Fill(lv[7].E());
        if(fabs(lv[2].Px())>ptm)ESup.Fill(lv[2].E());
        if(fabs(lv[1].Px())>ptm)ENup.Fill(lv[1].E());
      };
    };
    if(nphotons==2)
    {

      if(spin==uu || spin==ud){
        if( fabs(ph[8].Mag()-m0)<dm )WSdn3.Fill(ph[8].E(),fabs(ph[8].Px()));
        if( fabs(ph[7].Mag()-m0)<dm )WNdn3.Fill(ph[7].E(),fabs(ph[7].Px()));
        if( fabs(ph[2].Mag()-m0)<dm )ESdn3.Fill(ph[2].E(),fabs(ph[2].Px()));
        if( fabs(ph[1].Mag()-m0)<dm )ENdn3.Fill(ph[1].E(),fabs(ph[1].Px()));
      };

      if(spin==du || spin==dd){
        if( fabs(ph[8].Mag()-m0)<dm )WSup3.Fill(ph[8].E(),fabs(ph[8].Px()));
        if( fabs(ph[7].Mag()-m0)<dm )WNup3.Fill(ph[7].E(),fabs(ph[7].Px()));
        if( fabs(ph[2].Mag()-m0)<dm )ESup3.Fill(ph[2].E(),fabs(ph[2].Px()));
        if( fabs(ph[1].Mag()-m0)<dm )ENup3.Fill(ph[1].E(),fabs(ph[1].Px()));
      };


      if(spin==uu || spin==ud){

        if(fabs(ph[8].E())>10 && (fabs(ph[8].Mag()-m0)<dm))WSdn4.Fill(ph[8].E());
        if(fabs(ph[7].E())>10 && (fabs(ph[7].Mag()-m0)<dm))WNdn4.Fill(ph[7].E());
        if(fabs(ph[2].E())>10 && (fabs(ph[2].Mag()-m0)<dm))ESdn4.Fill(ph[2].E());
        if(fabs(ph[1].E())>10 && (fabs(ph[1].Mag()-m0)<dm))ENdn4.Fill(ph[1].E());

        if(fabs((ph[8].E())-Enom)<dEnom)WSdnm.Fill(ph[8].Mag());
        if(fabs((ph[7].E())-Enom)<dEnom)WNdnm.Fill(ph[7].Mag());
        if(fabs((ph[2].E())-Enom)<dEnom)ESdnm.Fill(ph[2].Mag());
        if(fabs((ph[1].E())-Enom)<dEnom)ENdnm.Fill(ph[1].Mag());

      };

      if(spin==du || spin==dd){
        if(fabs(ph[8].E())>10 && (fabs(ph[8].Mag()-m0)<dm))WSup4.Fill(ph[8].E());
        if(fabs(ph[7].E())>10 && (fabs(ph[7].Mag()-m0)<dm))WNup4.Fill(ph[7].E());
        if(fabs(ph[2].E())>10 && (fabs(ph[2].Mag()-m0)<dm))ESup4.Fill(ph[2].E());
        if(fabs(ph[1].E())>10 && (fabs(ph[1].Mag()-m0)<dm))ENup4.Fill(ph[1].E());

        if(fabs((ph[8].E())-Enom)<dEnom)WSupm.Fill(ph[8].Mag());
        if(fabs((ph[7].E())-Enom)<dEnom)WNupm.Fill(ph[7].Mag());
        if(fabs((ph[2].E())-Enom)<dEnom)ESupm.Fill(ph[2].Mag());
        if(fabs((ph[1].E())-Enom)<dEnom)ENupm.Fill(ph[1].Mag());
      };
    };
    if(nphotons==1 && (spin==uu || spin==ud) ){
      if(fabs(lv[8].E())>10)WSdnmm.Fill(lv[8].Mag());
      if(fabs(lv[7].E())>10)WNdnmm.Fill(lv[7].Mag());
      if(fabs(lv[2].E())>10)ESdnmm.Fill(lv[2].Mag());
      if(fabs(lv[1].E())>10)ENdnmm.Fill(lv[1].Mag());
    }
    if(nphotons==1 && (spin==du || spin==dd) ){
      if(fabs(lv[8].E())>10)WSupmm.Fill(lv[8].Mag());
      if(fabs(lv[7].E())>10)WNupmm.Fill(lv[7].Mag());
      if(fabs(lv[2].E())>10)ESupmm.Fill(lv[2].Mag());
      if(fabs(lv[1].E())>10)ENupmm.Fill(lv[1].Mag());
    }      
  };
  Ex.Write();
  Em1.Write();
  Em2.Write();
  Emm1.Write();
  Emm2.Write();
  mass.Write();
  massU.Write();
  Nphots.Write();
  WSup.Sumw2();
  WSdn.Sumw2();
  WNup.Sumw2();
  WNdn.Sumw2();

  WSupmm.Write();
  WSdnmm.Write();
  WNupmm.Write();
  WNdnmm.Write();

  ESupmm.Write();
  ESdnmm.Write();
  ENupmm.Write();
  ENdnmm.Write();

  TH1F asmCrossW("asmCrossW","asmCrossW",24,lowlim,highlim);
  TH1F asmCrossWp("asmCrossWp","asmCrossWp",24,lowlim,highlim);
  for(int j=1;j<=20;j++){
    Int_t nWSup=(Int_t) WSup.GetBinContent(j);
    Int_t nWSdn=(Int_t) WSdn.GetBinContent(j);
    Int_t nWNup=(Int_t) WNup.GetBinContent(j);
    Int_t nWNdn=(Int_t) WNdn.GetBinContent(j);

    Int_t npWSup=(Int_t) WSup4.GetBinContent(j);
    Int_t npWSdn=(Int_t) WSdn4.GetBinContent(j);
    Int_t npWNup=(Int_t) WNup4.GetBinContent(j);
    Int_t npWNdn=(Int_t) WNdn4.GetBinContent(j);

    Float_t rootn=sqrt(nWSup+nWSdn+nWNup+nWNdn);
    Float_t rootnp=sqrt(npWSup+npWSdn+npWNup+npWNdn);
    Float_t err=0;
    Float_t v=0;
    Float_t vv;
    if(rootn>0){
      err=1/rootn;
      vv=(sqrt(1.*nWSup*nWNdn)+sqrt(1.*nWNup*nWSdn));
      v=0;
      if(vv!=0)
      {
        v=(sqrt(1.*nWSup*nWNdn)-sqrt(1.*nWNup*nWSdn))/vv;
      };
    };

    asmCrossW.SetBinContent(j,v);
    asmCrossW.SetBinError(j,err);

    if(rootnp>0){
      err=1/rootnp;
      vv=(sqrt(1.*npWSup*npWNdn)+sqrt(1.*npWNup*npWSdn));
      v=0;
      if(vv!=0)
      {
        v=(sqrt(1.*npWSup*npWNdn)-sqrt(1.*npWNup*npWSdn))/vv;
      };
    };

    asmCrossWp.SetBinContent(j,v);
    asmCrossWp.SetBinError(j,err);

  };
  WSup.Write();
  WNup.Write();
  ESup.Write();
  ENup.Write();
  WSdn.Write();
  WNdn.Write();
  ESdn.Write();
  ENdn.Write();

  WSupm.Write();
  WNupm.Write();
  ESupm.Write();
  ENupm.Write();
  WSdnm.Write();
  WNdnm.Write();
  ESdnm.Write();
  ENdnm.Write();



  WSup2.Write();
  WNup2.Write();
  ESup2.Write();
  ENup2.Write();
  WSdn2.Write();
  WNdn2.Write();
  ESdn2.Write();
  ENdn2.Write();

  WSup3.Write();
  WNup3.Write();
  ESup3.Write();
  ENup3.Write();
  WSdn3.Write();
  WNdn3.Write();
  ESdn3.Write();
  ENdn3.Write();

  WSup4.Write();
  WNup4.Write();
  ESup4.Write();
  ENup4.Write();
  WSdn4.Write();
  WNdn4.Write();
  ESdn4.Write();
  ENdn4.Write();

  asmCrossW.Write();
  asmCrossWp.Write();
  ESup.Sumw2();
  ESdn.Sumw2();
  ENup.Sumw2();
  ENdn.Sumw2();

  TH1F asmCrossE("asmCrossE","asmCrossE",24,lowlim,highlim);
  TH1F asmCrossEp("asmCrossEp","asmCrossEp",24,lowlim,highlim);
  for(int j=1;j<=40;j++){
    Int_t nESup=(Int_t) ESup.GetBinContent(j);
    Int_t nESdn=(Int_t) ESdn.GetBinContent(j);
    Int_t nENup=(Int_t) ENup.GetBinContent(j);
    Int_t nENdn=(Int_t) ENdn.GetBinContent(j);

    Int_t npESup=(Int_t) ESup4.GetBinContent(j);
    Int_t npESdn=(Int_t) ESdn4.GetBinContent(j);
    Int_t npENup=(Int_t) ENup4.GetBinContent(j);
    Int_t npENdn=(Int_t) ENdn4.GetBinContent(j);

    Float_t rootn=sqrt(nESup+nESdn+nENup+nENdn);
    Float_t rootnp=sqrt(nESup+nESdn+nENup+nENdn);
    Float_t err=0;
    Float_t v=0;
    Float_t vv;
    if(rootn>0){
      err=1/rootn;
      vv=(sqrt(1.*nESup*nENdn)+sqrt(1.*nENup*nESdn));
      if(vv!=0){
        v=(sqrt(1.*nESup*nENdn)-sqrt(1.*nENup*nESdn))/vv;
      }
      else {
        v=0;
      };
    };
    asmCrossE.SetBinContent(j,v);
    asmCrossE.SetBinError(j,err);

    if(rootnp>0){
      err=1/rootnp;
      vv=(sqrt(1.*npESup*npENdn)+sqrt(1.*npENup*npESdn));
      if(vv!=0){
        v=(sqrt(1.*npESup*npENdn)-sqrt(1.*npENup*npESdn))/vv;
      }
      else {
        v=0;
      };
    };

    asmCrossEp.SetBinContent(j,v);
    asmCrossEp.SetBinError(j,err);

  };
  asmCrossE.Write();
  asmCrossEp.Write();
  return 0;
};

Int_t AnalTools::read1(char* filelist,Int_t set)
{
  Int_t DetMap[9]={0,1,2,0,0,5,6,3,4};
  Int_t spinindex=0;
  FilesSet* p_files=0;
  if(set>=61 && set <=63)
  {
    //transverse runs
    if(set==62)
    {
      //			   "root06T/fpdcorr.txt",
      p_files=new FilesSet("../",
          "run6_export/fpdped.txt",
          "run6_export/fpdgain.txt", 
          "run6_export/fpdcorr_itr10.2day.txt",
          "run6_export/fill.txt",
          "Fake",
          "run6_export/spinpat",
          "run6_export/geom_fy06_trans_survey.txt");
    }
    else
    {
      p_files=new FilesSet("../",
          "run6_export/fpdped.txt",
          "run6_export/fpdgain.txt", 
          "run6_export/gain_set1/fpdcorr_itr11.txt",
          "run6_export/fill.txt",
          "Fake",
          "run6_export/spinpat",
          "run6_export/geom_fy06_trans_survey.txt");
    };
    p_Geom=new Geom(p_files);
  }
  if(set>=65 && set<=68)
  {
    // longitudinal runs
    p_files=new FilesSet("/star/u/leun/fpd06/",
        "fpd++ped_20060425_corr.txt",
        "fpdgain_2006pp.txt", 
        "gain_set601/fpdcorr_itr5.txt",
        "fill.txt",
        "Fake",
        "spinpat",
        "geom_fy06_trans_survey.txt");
  };

  if(set==60)
  {
    p_files=new FilesSet("../../","fpd++ped_20060425.txt",
        "nogach06/fpdgain_2006pp.txt", 
        "nogach06/fpdcorr_2006pp_set1.txt",
        "fill_7.txt","Fake","spinpat6","geom.txt");
  };

  TString CalName=p_files->p_fpdgain()->Path();
  char calname[50];
  strcpy(calname,(const char*) CalName);
  CalibStr*  Rgain=new CalibStr(7116050,calname);
  CalName=p_files->p_fpdgaincorr()->Path();
  strcpy(calname,(const char*) CalName);
  CalibStr*  Rgainorigcorr=new CalibStr(7116050,calname);
  CalibStr*  Rgaincorr=new CalibStr(7116050,calname);
  CalibStr*  ptOnePct=new CalibStr(7116050,"ptonepct.txt");

  Rgainorigcorr->Print();
  Rgaincorr->Print();
  Rgainorigcorr->tm(2,2)->Print();
  for(Int_t iz=0;iz<14;iz++){ for(Int_t izz=0;izz<14;izz++){ printf("iz,izz,val=%d %d %f \n",iz,izz,
      Rgaincorr->GetValue(2,2,iz,izz));}; }
  Fill* RFill=0;
  TChain* p_OF=new TChain("p_out");
  FILE* pf;
  Int_t knt=0;
  Int_t FileCnt=0;
  Double_t Pi=TMath::Pi();
  std::cout<<" filelist= "<<filelist<<"\n";

  if(pf=fopen(filelist,"r"))
  {
    char fname[200];
    while(!feof(pf))
    {
      if(fscanf(pf,"%s\n",fname)>0)
      {
        if(knt==0)printf("%d files added to TChain\n",FileCnt);

        knt++;
        FileCnt++;
        p_OF->Add(fname);
      };
      if(knt>20)
      {
        knt=0;
      };
      // if(FileCnt>20)break;// max number of segs for testing
    };
  };

  if(OutputToSingle) {
    std::cout<<"Output to "<<OutputFileName<<"\n";
    p_OF->Merge( OutputFileName);
    return 0;
  };

  std::cout<<"OuputToSingle not selected\n";
  Float_t lowlim=0;
  Float_t highlim=120.;

  TH2F BBCsumEvsW("BBCsumEvsW","BBCsumEvsW",256,0,256,256,0,256);
  TH2F BBCsumEvsWetac("BBCsumEvsWetac","BBCsumEvsWetac",256,0,256,256,0,256);
  TH2F BBCtacEvsW("BBCtacEvsW","BBCtacEvsW",256,0,256,256,0,256);
  TH2F BBCtacEvsWetac("BBCtacEvsWetac","BBCtacEvsWetac",256,0,256,256,0,256);

  TH2F BBCsumEvseEast("BBCsumEvseEast","BBCsumEvseEast",256,0,256,10,0.,100);
  TH2F BBCsumWvseEast("BBCsumWvseEast","BBCsumWvseEast",256,0,256,10,0.,100);
  TH2F BBCsumEvseEast0("BBCsumEvseEast0","BBCsumEvseEast0",256,0,256,10,0.,100);
  TH2F BBCsumWvseEast0("BBCsumWvseEast0","BBCsumWvseEast0",256,0,256,10,0.,100);
  TH2F BBCsumEvsWeEn("BBCsumEvsWeEn","BBCsumEvsWeEn",256,0,256,256,0,256);
  TH2F BBCsumEvsWwEn("BBCsumEvsWwEn","BBCsumEvsWwEn",256,0,256,256,0,256);

  TH2F BBCsumEEta("BBCsumEEta","BBCsumEEta",256,0,256,30,0,30);
  TH2F BBCsumWEta("BBCsumWEta","BBCsumWEta",256,0,256,30,0,30);
  TH2F BBCsumEg150("BBCsumEg150","BBCsumEg150",256,0,256,30,0,30);
  TH2F BBCsumWg150("BBCsumWg150","BBCsumWg150",256,0,256,30,0,30);

  TH2F BBCsumEg1("BBCsumEg1","BBCsumEg1",256,0,256,30,0,30);
  TH2F BBCsumWg1("BBCsumWg1","BBCsumWg1",256,0,256,30,0,30);
  TH2F BBCsumEg0("BBCsumEg0","BBCsumEg0",256,0,256,30,0,30);
  TH2F BBCsumWg0("BBCsumWg0","BBCsumWg0",256,0,256,30,0,30);

  TH2F BBCsumEEta55("BBCsumEEta55","BBCsumEEta55",256,0,256,30,0,30);
  TH2F BBCsumWEta55("BBCsumWEta55","BBCsumWEta55",256,0,256,30,0,30);


  TH2F BBCsumEPi0("BBCsumEPi0","BBCsumEPi0",256,0,256,30,0,30);
  TH2F BBCsumWPi0("BBCsumWPi0","BBCsumWPi0",256,0,256,30,0,30);
  TH2F MvsTrig("MvsTrig","MvsTrig",12,0.,1.2,2,0,2);
  TH2F EvsTrig("EvsTrig","EvsTrig",12,0.,120,2,0,2);
  printf("Opening %s \n",(const char*) OutFileName);
  TFile* out=new TFile(OutFileName,"recreate");

  TH2F Ex("Ex","Ex",24,lowlim,highlim,30,0,30);

  TH2F Exw("Exw","Exw",24,lowlim,highlim,30,0,30);
  TH2F ExwPic("ExwPic","ExwPic",24,lowlim,highlim,30,0,30);
  TH2F ExwPicZ5("ExwPicZ5","ExwPicZ5",24,lowlim,highlim,30,0,30);
  TH2F ExwEta("ExwEta","ExEta",24,lowlim,highlim,30,0,30);
  TH2F ExwEtacZ5("ExwEtacZ5","ExEtacZ5",24,lowlim,highlim,30,0,30);
  TH2F ExwPicf("ExwPicf","ExwPicf",24,lowlim,highlim,30,0,30);
  TH2F ExwEtacf("ExwEtacf","ExwEtacf",24,lowlim,highlim,30,0,30);

  TH2F ExwEtac("ExwEtac","ExEtac",24,lowlim,highlim,30,0,30);


  TH2F ExIso("ExIso","ExIso",24,lowlim,highlim,30,0,30);
  TH2F ExNIso("ExNIso","ExNIso",24,lowlim,highlim,30,0,30);
  TH2F Ex1("Ex1","Ex1",24,lowlim,highlim,30,0,30);
  TH2F ExW("ExW","ExW",24,lowlim,highlim,30,0,30);
  TH2F En1("En1","En1",24,lowlim,highlim,30,0,30);
  TH2F En2("En2","En2",24,lowlim,highlim,30,0,30);
  TH2F Ex2En2("Ex2En2","Ex2En2",100,0,100,100,0,100);
  TH2F ExEta("ExEta","ExEta",24,lowlim,highlim,30,0,30);
  TH2F ExEtac("ExEtac","ExEtac",24,lowlim,highlim,30,0,30);
  TH2F ExPic("ExPic","ExPic",24,lowlim,highlim,30,0,30);

  TH2F ExEtacf("ExEtacf","ExEtacf",24,lowlim,highlim,30,0,30);
  TH2F ExPicf("ExPicf","ExPicf",24,lowlim,highlim,30,0,30);


  TH2F Em1("Em1","Em1",12,0.,1.2,30,0,30);
  TH2F Em2("Em2","Em2",48,0.,1.2,30,0,30);
  TH2F Em2EN("Em2EN","Em2EN",48,0.,1.2,81,0.,81);
  TH2F Em2ES("Em2ES","Em2ES",48,0.,1.2,81,0.,81);
  TH2F Em2ENwf("Em2ENwf","Em2ENwf",48,0.,1.2,81,0.,81);
  TH2F Em2ESwf("Em2ESwf","Em2ESwf",48,0.,1.2,81,0.,81);
  TH2F T7x7("T7x7","T7x7",7,0,7,7,0,7);
  TH3F T14x14("T14x14","T14x14",14,0,14,14,0,14,40,0,1.);
  TH3F TP14x14("TP14x14","TP14x14",14,0,14,14,0,14,40,0,1.);
  TH3F SP14x14("SP14x14","SP14x14",14,0,14,14,0,14,120,0,1.2);
  TH3F EP14x14("EP14x14","EP14x14",14,0,14,14,0,14,40,0,20.);
  TH3F TQ14x14("TQ14x14","TQ14x14",14,0,14,14,0,14,40,0,1.);
  TH3F SQ14x14("SQ14x14","SQ14x14",14,0,14,14,0,14,120,0,1.2);
  TH3F EQ14x14("EQ14x14","EQ14x14",14,0,14,14,0,14,40,0,20.);
  TH2F T014x14("T014x14","T014x14",14,0,14,14,0,14);
  TH2F Em2w("Em2w","Em2w",48,0.,1.2,30,0,30);
  TH2F Em230("Em230","Em230",48,0.,1.2,30,0,30);
  TH2F Em240("Em240","Em240",48,0.,1.2,30,0,30);
  TH2F Em250("Em250","Em250",48,0.,1.2,30,0,30);
  TH2F Em2w40("Em2w40","Em2w40",48,0.,1.2,30,0,30);
  TH2F Em2w45("Em2w45","Em2w45",48,0.,1.2,30,0,30);
  TH2F Em2w50("Em2w50","Em2w50",48,0.,1.2,30,0,30);

  TH2F E2phvsEm2w("E2phvsEm2w","E2phvsEm2w",48,0,1.2,20,0,100);
  TH2F E2phvsEm2wf("E2phvsEm2wf","E2phvsEm2wf",48,0,1.2,20,0,100);
  TH2F E2phvsEm2wfc("E2phvsEm2wfc","E2phvsEm2wfc",48,0,1.2,20,0,100);
  TH2F E2phvsEm2wc("E2phvsEm2wc","E2phvsEm2wc",48,0,1.2,20,0,100);
  TH2F E2phvsEm2wcZ5("E2phvsEm2wcZ5","E2phvsEm2wcZ5",48,0,1.2,20,0,100);
  TH2F E2phvsEm2wcZg5("E2phvsEm2wcZg5","E2phvsEm2wcZg5",48,0,1.2,20,0,100);

  TH2F E2phvsEm2("E2phvsEm2","E2phvsEm2",48,0,1.2,20,0,100);
  TH2F E2phvsEm2c("E2phvsEm2c","E2phvsEm2c",48,0,1.2,20,0,100);
  TH2F E2phvsEm2f("E2phvsEm2f","E2phvsEm2f",48,0,1.2,20,0,100);
  TH2F E2phvsEm2fc("E2phvsEm2fc","E2phvsEm2fc",48,0,1.2,20,0,100);

  TH2F Em2c("Em2c","Em2c",48,0.,1.2,30,0,30);
  TH2F Em2c40("Em2c40","Em2c40",48,0.,1.2,30,0,30);
  TH2F Em2c50("Em2c50","Em2c50",48,0.,1.2,30,0,30);
  TH2F Em2c55("Em2c55","Em2c55",48,0.,1.2,30,0,30);
  TH2F Em2c40w("Em2c40w","Em2c40w",48,0.,1.2,30,0,30);
  TH2F Em2c50w("Em2c50w","Em2c50w",48,0.,1.2,30,0,30);
  TH2F Em2c55w("Em2c55w","Em2c55w",48,0.,1.2,30,0,30);

  TH2F Em1c("Em1c","Em1c",24,0.,1.2,30,0,30);
  TH2F Phi2G("Phi2G","Phi2G",40,0.,Pi,30,0,30);
  TH2F Z("Z","Z",40,0.,1.,30,0,30);
  TH2F YphiPi("YphiPi","YphiPi",100,-5,5,100,-Pi,Pi);
  TH2F YphiEta("YphiEta","YphiEta",100,-5,5,100,-Pi,Pi);
  TH2F Emm1("Emm1","Emm1",120,0,12.,30,0,30);
  TH2F Emm2("Emm2","Emm2",480,0.,12.,30,0,30);
  TH2F Phi2GCEta("Phi2GCEta","Phi2GCEta",40,0.,Pi,30,0,30);
  TH2F ZCEta("ZCEta","ZCEta",40,0.,1.,30,0,30);
  TH2F ZCEta50("ZCEta50","ZCEta50",40,0.,1.,30,0,30);
  TH2F Phi2GCPi("Phi2GCPi","Phi2GCPi",40,0.,Pi,30,0,30);
  TH2F ZCPi("ZCPi","ZCPi",40,0.,1.,30,0,30);
  TH2F ZCPi50("ZCPi50","ZCPi50",40,0.,1.,30,0,30);

  TH2F ZMid("ZMid","ZMid",40,0.,1.,30,0,30);
  TH2F ZMid40("ZMid40","ZMid40",40,0.,1.,30,0,30);
  TH2F ZMid40c("ZMid40c","ZMid40c",40,0.,1.,30,0,30);
  TH2F ZMid50("ZMid50","ZMid50",40,0.,1.,30,0,30);

  TH2F ZMida("ZMida","ZMida",40,0.,1.,30,0,30);
  TH2F ZMid40a("ZMid40a","ZMid40a",40,0.,1.,30,0,30);
  TH2F ZMid40ac("ZMid40ac","ZMid40ac",40,0.,1.,30,0,30);
  TH2F ZMid50a("ZMid50a","ZMid50a",40,0.,1.,30,0,30);

  TH2F ZMidb("ZMidb","ZMidb",40,0.,1.,30,0,30);
  TH2F ZMid40b("ZMid40b","ZMid40b",40,0.,1.,30,0,30);
  TH2F ZMid40bc("ZMid40bc","ZMid40bc",40,0.,1.,30,0,30);
  TH2F ZMid50b("ZMid50b","ZMid50b",40,0.,1.,30,0,30);

  TH2F YTphiEtac("YTphiEtac","YTphiEtac",1000,-5,5,157,-Pi/4,Pi/4);
  TH2F YTphiPic("YTphiPic","YTphiPic",1000,-5,5,157,-Pi/4,Pi/4);
  TH2F YPhi1p20("YPhi1p20","YPhi1p20",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhi1p20f("YPhi1p20f","YPhi1p20f",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhi1p4050("YPhi1p4050","YPhi1p4050",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhi1p5060("YPhi1p5060","YPhi1p5060",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi4050f("YPhiPi4050f","YPhiPi4050f",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi5060f("YPhiPi5060f","YPhiPi5060f",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi6070f("YPhiPi6070f","YPhiPi6070f",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi4050("YPhiPi4050","YPhiPi4050",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi5060("YPhiPi5060","YPhiPi5060",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F YPhiPi6070("YPhiPi6070","YPhiPi6070",1000,-5,5,471,-1*Pi/4,5*Pi/4);
  TH2F Pt1Pt2("Pt1Pt2","Pt1Pt2",100,-5,5,100,-5,5);
  TH2F AbsPt1("AbsPt1","AbsPt1",50,0,7,28,0,28);
  TH2F m2p("m2p","m2p",50,0,5,30,0,30);
  TH2F RNum("RNUM","RNUM",200,7110000,7130000,30,0,30);
  TH2F Bnchid("Bnchid","Bnchid",120,0,120,30,0,30);
  TH2F BnchidvsEW("BnchidvsEW","BnchidvsEW",120,0,120,10,0,100);
  TH2F BnchidvsEE("BnchidvsEE","BnchidvsEE",120,0,120,10,0,100);
  TH2F BnchidvsEc("BnchidvsEc","BnchidvsEc",120,0,120,10,0,100);
  TH2F BnchidvsEEtac("BnchidvsEEtac","BnchidvsEEtac",120,0,120,10,0,100);
  TH2F hFill("hFill","hFill",100,7750,7850,30,0,30);
  TH2F YFillN50("YFillN50","YFillN50",100,7750,7850,200,-5,5);
  TH2F YFillS50("YFillS50","YFillS50",100,7750,7850,200,-5,5);
  TH2F ExFill("Ex2Fill","Ex2Fill",100,7750,7850,100,0,100.);
  TH2F ESumFill("ESumFill","ESumFill",100,7750,7850,100,0,100.);
  TH2F Em2Fill("Em2Fill","Em2Fill",100,7750,7850,100,0.,1.);
  TH2F ZvsE1("ZvsE1","ZvsE1",40,0.,1.,24,0.,120.);
  TH2F ZvsEPic("ZvsEPic","ZvsEPic",40,0.,1.,24,0.,120.);
  TH2F ZvsEEtac("ZvsEEtac","ZvsEEtac",40,0.,1.,24,0.,120.);
  TH2F EP1_XY("EP1_XY","EP1_XY",28,-14,14,28,-14,14);
  TH2F EP2_XY("EP2_XY","EP2_XY",28,-14,14,28,-14,14);
  TH2F EP1Pi0_XY("EP1Pi0_XY","EP1Pi0_XY",28,-14,14,28,-14,14);
  TH2F EP2Pi0_XY("EP2Pi0_XY","EP2Pi0_XY",28,-14,14,28,-14,14);
  TH2F EP1Eta_XY("EP1Eta_XY","EP1Eta_XY",28,-14,14,28,-14,14);
  TH2F EP2Eta_XY("EP2Eta_XY","EP2Eta_XY",28,-14,14,28,-14,14);
  TH1F Mio("Mio","Mio",120,0,1.2);
  TH2F Miosp("Miosp","Miosp",120,0,1.2,30,0,30);
  TH2F MioL1("MioL1","MioL1",50,0,50,120,0,1.2);
  TH2F MioL2("MioL2","MioL2",50,0,50,120,0,1.2);
  TH2F MioL3("MioL3","MioL3",50,0,50,120,0,1.2);
  TH2F WNinvsWSout("WNinvsWSout","WNinvsWSout",100,0,100,100,0,100);
  TH2F WSinvsWNout("WSinvsWNout","WSinvsWNout",100,0,100,100,0,100);
  TH2F WNinvsWNout("WNinvsWNout","WNinvsWNout",100,0,100,100,0,100);
  TH2F WSinvsWSout("WSinvsWSout","WSinvsWSout",100,0,100,100,0,100);
  TH2F WNinvsWSin("WNinvsWSin","WNinvsWSin",100,0,100,100,0,100);
  TH2F Y5_7("Y5_7","Y5_7",14,0,14,6,0,6);
  TH1F mv7("mv7","mv7",200,0.,4.);
  TH1F mv5("mv5","mv5",200,0.,4.);
  TH1F mv57("mv57","mv57",200,0.,4.);
  TH2F sepvse("sepvse","sepvse",100,0,14,50,0,100);
  TH2F m5_42_44("m5_42_44","m5_42_44",50,0,50,50,0.,1.);
  TH3F d5_44("d5_44","d5_44",14,0,14,14,0,14,60,0,1.2);
  TH3F e5_44("e5_44","e5_44",14,0,14,14,0,14,120,0,60);
  TH2F m2comp("m2comp","m2comp",100,0,1.,100,0,1.);
  TH2F EcompA("EcompA","EcompA",100,0,100.,100,0,100);
  TH2F EcompB("EcompB","EcompB",100,0,100.,100,0,100);
  //1 tr
  TH2F Em1c401Tr("Em1c401Tr","Em1c401Tr",48,0.,1.2,30,0,30);
  TH2F Em1c501Tr("Em1c501Tr","Em1c501Tr",48,0.,1.2,30,0,30);
  TH2F Em1c601Tr("Em1c601Tr","Em1c601Tr",48,0.,1.2,30,0,30);
  TH2F Exc1Trw("Exc1Trw","Exc1Trw",24,0.,120,30,0,30);
  TH2F Em1cvsEn1_1Tr("Em1cvsEn1_1Tr","Em1cvsEn1_1Tr",48,0.,1.2,40,1,100);
  TH2F Em1cvsEn1_2Tr("Em1cvsEn1_2Tr","Em1cvsEn1_2Tr",48,0.,1.2,40,1,100);
  TH2F Em2cvsEmm2_2Tr("Em2cvsEmm2_2Tr","Em2cvsEmm2_2Tr",48,0.,1.2,48,0,1.2);
  TH2F Em2cvsEmm2_2Tr40("Em2cvsEmm2_2Tr40","Em2cvsEmm2_2Tr40",48,0.,1.2,48,0,1.2);
  TH2F Em2cvsEmm2_2Tr50("Em2cvsEmm2_2Tr50","Em2cvsEmm2_2Tr50",48,0.,1.2,48,0,1.2);
  TH2F Em2cvsEmm2_2Tr60("Em2cvsEmm2_2Tr60","Em2cvsEmm2_2Tr60",48,0.,1.2,48,0,1.2);
  lowlim=0.;
  highlim=2.;

  Int_t pcnt=0;
  TVector3 v3;
  TLorentzVector v4;
  TH2F  Nphots("Nphots","Nphots",30,0,30,30,0,30);
  TH2F  NphvsEsum("NphvsEsum","NphvsEsum",6,0,6,24,1,120);
  TH2F  NphvsMsum("NphvsMsum","NphvsMsum",6,0,6,48,0,1.2);

  Int_t tpes[320];
  Int_t spin;
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t nphotons;
  Int_t Rnum;
  Int_t FillNumber=0;
  Int_t LastRun=0;
  Int_t Bunchid7bit;
  Int_t EventN=0;
  UChar_t BBcSums[5];
  p_OF->SetBranchAddress("spin",&spin);
  p_OF->SetBranchAddress("nphotons",&nphotons);
  p_OF->SetBranchAddress("br_nwrds",&nwrds);
  p_OF->SetBranchAddress("br_types",tpes);
  p_OF->SetBranchAddress("br_pxyzt",pxyzt);
  p_OF->SetBranchAddress("br_Rnum",&(Rnum));
  p_OF->SetBranchAddress("br_Bunchid7bit",&(Bunchid7bit));
  p_OF->SetBranchAddress("br_BBcSums",BBcSums);
  if(p_OF->FindBranch("br_EventN"))
  {
    p_OF->SetBranchAddress("br_EventN",&EventN);
  };

  Double_t NearUnique;
  //  TFile* new_P_File=new TFile("OF2.root","recreate");
  //  TTree* new_p_OF=p_OF->CloneTree(0);
  //  new_p_OF->Branch("br_NearUnique",&NearUnique,"NearUnique/D");
  out->cd();

  TLorentzVector vec[80];
  Float_t DepFactor[10];
  Int_t nentries=p_OF->GetEntries();
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;
  printf("GetEntries=%d \n",nentries);
  knt=0;
  Bool_t starting=true;
  Float_t m0=.135;
  Float_t dm=.05;
  Float_t Enom=40.;
  Float_t dEnom=40;
  TLorentzVector LastV6(0,0,0,0);
  TLorentzVector LastV5(0,0,0,0);
  for(Int_t i=0;i<nentries;i++)
  {

    Int_t nbytes=p_OF->GetEntry(i);
    knt++;
    NearUnique=Rnum*1000000000.;
    NearUnique+=Bunchid7bit*1000000+BBcSums[2]*1000+BBcSums[4];

    if(knt>100000 || (starting && knt>10000)){
      printf("Evt: %d \n",i);
      knt=1;
      if(i>100000)starting=false;

    };

    Int_t kk=0;
    TLorentzVector ph[10];
    Int_t nlv[10];
    Int_t plv[10][10];
    Int_t nph2[10];
    Int_t pph2[10][10];
    Int_t nph3[10];
    Int_t pph3[10][10];
    Float_t hiEsum[10],loEsum[10];

    for(Int_t i1=0;i1<10;i1++)
    {
      nlv[i1]=0;
      nph2[i1]=0;
      nph3[i1]=0;
      DepFactor[i1]=1.;
      hiEsum[i1]=0;
      loEsum[i1]=0;
    };
    Int_t vtyp;
    Int_t nwords=nwrds;

    if(Rnum!=LastRun)
    {
      if(RFill!=0)delete RFill;
      RFill=new Fill(p_files,6900,8000,Rnum,Rnum);
      RFill->SetFillNumberforRun(Rnum);
      Int_t fNum=RFill->GetFillNumber();
      if(fNum!=FillNumber)RFill->Print();
      FillNumber=fNum;
      pcnt=0;
    };

    LastRun=Rnum;
    Int_t bSpin= RFill->BlueSpin(Bunchid7bit);
    Int_t ySpin= RFill->YellowSpin(Bunchid7bit);
    Int_t newspin= bSpin+1+(ySpin+1)/2;
    if(RFill->pfd->kicked(Bunchid7bit) || bSpin*ySpin==0)newspin+=29;
    if(newspin<21 && newspin!=spin && (pcnt++)<10)std::cout<<Rnum<<" : "<<i<<" "<<"spin="<<spin<<" newspin="<<newspin<<"bunchid="<<Bunchid7bit<<"\n";
    spin=newspin;

    //Trigger Selection ******************
    Bool_t trigsel=false;
    if(BBcSums[0]==1)trigsel=true;
    if(BBcSums[2]>10){
      if(BBcSums[4]>20 && BBcSums[4]<175)trigsel=true;
    };
    //*beg trigsel
    if(trigsel||true){
      //                ******************
      /*	if(((UChar_t) BBcSums[1])<10 ||((UChar_t)  BBcSums[2])<10){
                printf("BBcSums:%x (%x ) (%x )  (%x ) (%x ) \n",(UChar_t) BBcSums[0],
                (UChar_t) BBcSums[1],(UChar_t) BBcSums[2],
                (UChar_t) BBcSums[3],(UChar_t) BBcSums[4]);
                };
                */
      for(Int_t k=0;k<nwords;k+=4){
        Bool_t cut56=true;
        if(tpes[k]<9)
        {
          vtyp=tpes[k];
          plv[vtyp][nlv[vtyp]]=k;
          vec[k].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);

          if(nlv[vtyp]==0)
          {
            DepFactor[vtyp]=SingleEdepFactor(vtyp,vec[k].E());
          };
          nlv[vtyp]++;
        }
        else if((tpes[k]>300) && (tpes[k]<309))
        {
          vtyp=tpes[k]-300;
          cut56=true;
          pph2[vtyp][nph2[vtyp]]=k;
          pph3[vtyp][nph3[vtyp]]=k;
          vec[k].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
          nph2[vtyp]++;
        };
      };


      if(nph2[5]==4 && nph2[6]==4)continue;
      if(nph2[5]==4 && nph2[6]==4)printf("this should not be printing \n");
      if(nph3[6]<3)WNinvsWSout.Fill(hiEsum[7],hiEsum[6]);
      if(nph3[5]<3)WSinvsWNout.Fill(hiEsum[8],hiEsum[5]);
      if(nph3[5]<3)WNinvsWNout.Fill(hiEsum[7],hiEsum[5]);
      if(nph3[6]<3)WSinvsWSout.Fill(hiEsum[8],hiEsum[6]);
      WNinvsWSin.Fill(hiEsum[7],hiEsum[8]);
      /*
         if((nph2[5]>1 && nph2[7]>1) || (nph2[6]>1 && nph2[8]>1) )
         {
         TLorentzVector v7;
         for(Int_t ixx=0;ixx<nph2[7];ixx++)v7+=vec[pph2[7][ixx]];
         TLorentzVector v5;
         for(Int_t ixx=0;ixx<nph2[5];ixx++)v5+=vec[pph2[5][ixx]];
         if(v7.E()>20)mv7.Fill(v7.Mag());
         if(v5.E()>10)mv5.Fill(v5.Mag());
         if((v5+v7).E()>40)mv57.Fill((v5+v7).Mag());
         Int_t tspn;
         if(spin>=0 && spin<4)tspn=spin;
         spinindex=DetMap[5]*4+tspn;
         ExW.Fill(hiEsum[5]+hiEsum[7]+loEsum[5]+loEsum[7],spinindex);
         spinindex=DetMap[6]*4+tspn;
         ExW.Fill(hiEsum[6]+hiEsum[8]+loEsum[6]+loEsum[8],spinindex);
         };
         if(nph2[5]==1 && nph2[7]==1)
         {
         Float_t Etot5=vec[plv[5][0]].E();
         Float_t Rt5=0;
         if(Etot5>0)Rt5=vec[pph2[5][0]].E()/Etot5;
         Bool_t ok5=fabs(Rt5-1.1)<.15;
         for(Int_t iz=0;iz<nph2[5]&&ok5;iz++)
         {
         Int_t tspn;
         if(spin>=0 && spin<4)tspn=spin;
         spinindex=DetMap[5]*4+tspn;
         TLorentzVector vsm=vec[pph2[5][iz]];
      //TLorentzVector vsm=LastV5;
      //LastV5=vec[pph2[5][iz]];
      TVector3 Sa=vsm.Vect();
      TVector3 Sb=vec[pph2[7][0]].Vect();
      Sa=Sa*(1./Sa.Z()*(*p_Geom->ZFPD(2,1)));
      Sb=Sb*(1./Sb.Z()*(*p_Geom->ZFPD(2,3)));
      TVector3 v51=p_Geom->LocalXYZ(2,1,Sa,true);
      TVector3 v71=p_Geom->LocalXYZ(2,3,Sb,true);
      Int_t ixa1=(Int_t) v51.X();
      Int_t iyb1=(Int_t) v51.Y();
      Float_t fudge=1.;
      Float_t Elimptone=2.;
      if(ixa1>=0 && iyb1>=0 && ixa1<14 && iyb1<14)
      {

      Float_t rfa=Rgainorigcorr->GetValue(2,1,iyb1,ixa1);
      Float_t rfb=Rgaincorr->GetValue(2,1,iyb1,ixa1);
      if(rfa>0)fudge=sqrt(rfb/rfa);
      Float_t ptone=ptOnePct->GetValue(2,1,iyb1,ixa1);
      Elimptone=ptone*Rgain->GetValue(2,1,iyb1,ixa1)*rfb;
      //		    printf("row=%d col=%d Elimptone=%f\n",iyb1,ixa1,Elimptone);
      };
      TLorentzVector vsm7=vec[pph2[7][0]];
      TLorentzVector v0sm=vec[plv[5][iz]];
      TLorentzVector  vmio=(vsm+vec[pph2[7][0]]);
      TVector3 vcent5(7.,7.,v51.Z());
      if(vsm.E()*fudge*fudge>Elimptone*.25 && vmio.E()*fudge>10. && v71.X()>.5 && v71.X()<5.5 && (v51-vcent5).Mag()>2.5   )
      {
      Mio.Fill(vmio.Mag()*fudge);
      SP14x14.Fill(v51.X(),v51.Y(),vmio.Mag()*fudge);
      Miosp.Fill(vmio.Mag()*fudge,spinindex);
      TVector3 v3d=vec[pph2[5][iz]].Vect();
      TVector3 vv3d=v3d*(1./v3d.Z()*(*(p_Geom->ZFPD(2,1))));
      TVector3 v3d7=vec[pph2[7][0]].Vect();
      TVector3 vv3d7=v3d7*(1./v3d7.Z()*(*(p_Geom->ZFPD(2,3))));

      if(fabs(vsm.E()-2*fudge*fudge)<1)MioL1.Fill((vv3d7-vv3d).Mag(),vmio.Mag()*fudge);
      if(fabs(vsm.E()-5*fudge*fudge)<3)MioL2.Fill((vv3d7-vv3d).Mag(),vmio.Mag()*fudge);
      if(fabs(vsm.E()-11)<3*fudge*fudge)MioL3.Fill((vv3d7-vv3d).Mag(),vmio.Mag()*fudge);
      TVector3 v5=p_Geom->LocalXYZ(2,1,vv3d,true);
      T014x14.Fill(v5.X(),v5.Y());
      T14x14.Fill(v5.X(),v5.Y(),vmio.Mag());
      TVector3 v7=p_Geom->LocalXYZ(2,3,vv3d7,true);
      if( ((Int_t) v7.X()) ==0 && ((Int_t) v5.X())==4)Y5_7.Fill(v5.Y(),v7.Y());

    };
    };

    };
    if(nph2[6]==1 && nph2[8]==1)
    {  
      Float_t Etot6=vec[plv[6][0]].E();
      Float_t Rt6=0;
      if(Etot6>0)Rt6=vec[pph2[6][0]].E()/Etot6;
      Bool_t ok6=fabs(Rt6-1.1)<.15;

      for(Int_t iz=0;iz<nph2[6]&&ok6;iz++)
      { 
        Int_t tspn;
        if(spin>=0 && spin<4)tspn=spin;
        spinindex=DetMap[6]*4+tspn;
        TLorentzVector vsm=vec[pph2[6][iz]];
        //	    TLorentzVector vsm=LastV6;
        //	    LastV6=vec[pph2[6][0]];
        TLorentzVector v0sm=vec[plv[6][iz]];
        TLorentzVector  vmio=(vsm+vec[pph2[8][0]]);
        TVector3 Sa=vsm.Vect();
        TVector3 Sb=vec[pph2[8][0]].Vect();
        Sa=Sa*(1./Sa.Z()*(*p_Geom->ZFPD(2,2)));
        Sb=Sb*(1./Sb.Z()*(*p_Geom->ZFPD(2,4)));
        TVector3 v51=p_Geom->LocalXYZ(2,2,Sa,true);
        TVector3 v71=p_Geom->LocalXYZ(2,4,Sb,true);

        Int_t ixa1=(Int_t) v51.X();
        Int_t iyb1=(Int_t) v51.Y();
        Float_t fudge=1.;
        Float_t Elimptone=2.;
        if(ixa1>=0 && iyb1>=0 && ixa1<14 && iyb1<14)
        {

          Float_t rfa=Rgainorigcorr->GetValue(2,2,iyb1,ixa1);
          Float_t rfb=Rgaincorr->GetValue(2,2,iyb1,ixa1);
          if(rfa>0)fudge=sqrt(rfb/rfa);
          Float_t ptone=ptOnePct->GetValue(2,2,iyb1,ixa1);
          Elimptone=ptone*Rgain->GetValue(2,2,iyb1,ixa1)*rfb;
          //		    printf("row=%d col=%d Elimptone=%f\n",iyb1,ixa1,Elimptone);
        };
        TVector3 vcent5(7.,7.,v51.Z());
        if(vsm.E()*fudge*fudge>Elimptone*.25 && vmio.E()*fudge>10 && v71.X()>.5 && v71.X()<5.5 && (v51-vcent5).Mag()>2.5  )	    
        {
          Mio.Fill((vmio.Mag())*fudge);
          SQ14x14.Fill(v51.X(),v51.Y(),vmio.Mag()*fudge);
          Miosp.Fill(vmio.Mag()*fudge,spinindex);
        };
      };

    };
    */

      //*beg kloop
      for(Int_t k=0;k<9;k++)
      {
        for(Int_t zk=0;zk<nlv[k];zk++)
        {
          Float_t en00=vec[plv[k][zk]].E();
          vec[plv[k][zk]]*=DepFactor[k];
          Float_t en01=vec[plv[k][zk]].E();
          EcompA.Fill(en00,en01);
        };

        for(Int_t zk=0;zk<nph2[k];zk++)
        {
          Float_t en00=vec[pph2[k][zk]].E();
          vec[pph2[k][zk]]*=DepFactor[k];
          Float_t en01=vec[pph2[k][zk]].E();
          EcompB.Fill(en00,en01);
        };

        Int_t cellN=0;
        Bool_t fcut1=true;
        Bool_t fcut2=true;
        Int_t EW=1;
        Float_t zDet=0;
        Float_t yDet=0;
        Float_t xDet=0;
        Float_t cellwidth=1.0;


        Int_t NSTB=k;
        if(k>4)
        {
          EW=2;
          NSTB=k-4;
        };
        if(k>0){
          zDet=*(p_Geom->ZFPD(EW,NSTB));
          xDet=*(p_Geom->xOffset(EW,NSTB));
          yDet=*(p_Geom->yOffset(EW,NSTB));
          cellwidth=*(p_Geom->FpdTowWid(EW,NSTB));
        };
        Bool_t Center1=false;
        Bool_t Centercut=false;
        ph[k]=TLorentzVector(0.,0.,0.,0.);
        Int_t spn=8;
        if(spin>=0 && spin<4)spn=spin;
        if(newspin>4)spn=100;
        spinindex=DetMap[k]*4+spn;
        TLorentzVector vsum(0.,0.,0.,0.);
        Nphots.Fill(nph2[k],spinindex);

        //* beg nlv[k]>0

        if(nlv[k]>0)
        {

          for(Int_t jk=0;jk<nlv[k];jk++)
          {
            vsum+=vec[plv[k][jk]];
          };
          if(k>3){
            BnchidvsEW.Fill(Bunchid7bit,vsum.E());
            BBCsumEvsWwEn.Fill(BBcSums[2],BBcSums[1]);
          };
          if(k<=3){
            BnchidvsEE.Fill(Bunchid7bit,vsum.E());
            BBCsumEvsWeEn.Fill(BBcSums[2],BBcSums[1]);
            BBCsumEvseEast0.Fill(BBcSums[1],vsum.E());
            BBCsumWvseEast0.Fill(BBcSums[2],vsum.E());

          };

          if(vsum.E()>50){
            BBCsumEg150.Fill(BBcSums[1],spinindex);
            BBCsumWg150.Fill(BBcSums[2],spinindex);
            if(fabs(vsum.Phi())<.5){
              YFillN50.Fill(FillNumber,vsum.Eta());
            } else if( (vsum.Phi()-Pi)<.5 ) {
              YFillS50.Fill(FillNumber,vsum.Eta());
            };
          }
          else if (vsum.E()>40){
            BBCsumEg1.Fill(BBcSums[1],spinindex);
            BBCsumWg1.Fill(BBcSums[2],spinindex);
          }
          else {
            BBCsumEg0.Fill(BBcSums[1],spinindex);
            BBCsumWg0.Fill(BBcSums[2],spinindex);
          };
          NphvsEsum.Fill(nph2[k],vsum.E());
          NphvsMsum.Fill(nph2[k],vsum.Mag());
          Emm1.Fill(vsum.Mag(),spinindex);
          Emm2.Fill(vsum.Mag(),spinindex);
          if(true)
          {
            Bnchid.Fill(Bunchid7bit,spinindex);
            RNum.Fill(Rnum,spinindex);
            hFill.Fill(FillNumber,spinindex);
            ESumFill.Fill(FillNumber,vsum.E());
          };
        };
        //* end nlv[k]>0
        //* beg nph2[k]>0
        if(nph2[k]>0)
        {	      
          TLorentzVector v1;
          for(Int_t jph=0;jph<nph2[k];jph++)
          {
            v1=vec[pph2[k][jph]];

            if(k==5 && hiEsum[8]>20)
            {

              TVector3 vn1=v1.Vect();
              TVector3 vnorm1=vn1*(1./vn1.Z())*(*(p_Geom->ZFPD(2,1)));
              TVector3 v5=p_Geom->LocalXYZ(2,1,vnorm1,true);
              Int_t ixa1=(Int_t) v5.X();
              Int_t iyb1=(Int_t) v5.Y();
              if(ixa1>=0 && ixa1<14 && iyb1>=0 && iyb1<14)
              {
                Float_t rfa=Rgainorigcorr->GetValue(2,1,iyb1,ixa1);
                Float_t rfb=Rgaincorr->GetValue(2,1,iyb1,ixa1);
                Float_t gain=1;
                if(rfa>0)gain=rfb/rfa;
                EP14x14.Fill(v5.X(),v5.Y(),v1.E()*gain);
              };
            };
            if(k==6 && hiEsum[7]>20)
            {

              TVector3 vn1=v1.Vect();

              TVector3 vnorm1=vn1*(1./vn1.Z())*(*(p_Geom->ZFPD(2,2)));
              TVector3 v5=p_Geom->LocalXYZ(2,2,vnorm1,true);
              Int_t ixa1=(Int_t) v5.X();
              Int_t iyb1=(Int_t) v5.Y();
              if(ixa1>=0 && ixa1<14 && iyb1>=0 && iyb1<14)
              {
                Float_t rfa=Rgainorigcorr->GetValue(2,2,iyb1,ixa1);
                Float_t rfb=Rgaincorr->GetValue(2,2,iyb1,ixa1);
                //			    printf("iyb1=%d,ixa1=%d,rfb=%f rfa=%f \n",iyb1,ixa1,rfb,rfa);
                Float_t gain=1;
                if(rfa>0)gain=rfb/rfa;
                EQ14x14.Fill(v5.X(),v5.Y(),v1.E()*gain);
                Float_t tmean=EQ14x14.ProjectionZ("tmpEQ",4,4,8,8)->GetMean();
                //			    if(iyb1==7 && ixa1==3)printf("EQ14x14(7,3): gain=%f tmean=%f \n",gain,tmean);
              };
            };
          };
          TLorentzVector v2(0.,0.,0.,0.);
          TVector3 V3_1(0.,0.,0.);
          TVector3 V3_2(0.,0.,0.);
          if(v1.Z()!=0)
          {
            V3_1=TVector3(v1.X(),v1.Y(),v1.Z());
            V3_1.SetMag(V3_1.Mag()*zDet/v1.Z());
            Float_t fposx,fposy;
            fposx=((V3_1.X()-xDet)/cellwidth);
            fposy=((V3_1.Y()-yDet)/cellwidth);
            if(k==1)EP1_XY.Fill(fposx,fposy);
            if(k==2)EP2_XY.Fill(fposx,fposy);
            Float_t dxpos=-3.5;
            Float_t dxrange=3;
            Float_t dypos=-3.5;
            if(zDet>0){dxpos=-3;dypos=-3;dxrange=2.5;};
            if(xDet>0 ){dxpos*=-1;};
            if(k==1)cellN=T7x7.FindBin(-fposx,-fposy);
            if(k==2)cellN=T7x7.FindBin(fposx,-fposy);
            //	printf("k %d fposx=%f dxpos=%f fposy=%f cellN=%d \n",  k,fposx,dxpos,fposy,cellN);
            if(fabs(fposx-dxpos)>dxrange)fcut1=false;
            if(fabs(fposy-dypos)>dxrange)fcut1=false;
            if(nph2[k]>1){
              v2=vec[pph2[k][1]];
              V3_2=TVector3(v2.X(),v2.Y(),v2.Z());
              if(v2.Z()!=0){
                V3_2.SetMag(V3_2.Mag()*zDet/v2.Z());
                fposx=((V3_2.X()-xDet)/cellwidth);
                fposy=((V3_2.Y()-yDet)/cellwidth);
                if(fabs(fposx-dxpos)>dxrange)fcut2=false;
                if(fabs(fposy-dypos)>dxrange)fcut2=false;
              };
            };
          };
          Ex1.Fill(v1.E(),spinindex);
          En1.Fill(vsum.E(),spinindex);
          if(v1.E()>20)
          {
            Float_t phiph=v1.Phi();
            phiph=fmod(phiph+2.5*Pi,2*Pi)-Pi/2.;
            YPhi1p20.Fill(v1.Eta(),phiph);
            if(fcut1)YPhi1p20f.Fill(v1.Eta(),phiph);
            if(v1.E()>40&&v1.E()<50)YPhi1p4050.Fill(v1.Eta(),phiph);
            if(v1.E()>50&&v1.E()<60)YPhi1p5060.Fill(v1.Eta(),phiph);
          }; 

          TLorentzVector va,vz;
          TLorentzVector vb(0.,0.,0.,0.);
          Float_t z;
          va=vec[pph2[k][0]];
          vz=va;
          Center1= (pow(vz.Eta()+3.65,2)+pow(tan(vz.Phi()),2)<pow(.15,2));

          //* beg nph2[k]>1
          if(nph2[k]>1)
          {
            vb=vec[pph2[k][1]];
            vz=va+vb;

            Centercut=
              (pow(vz.Eta()+3.65,2)+pow(tan(vz.Phi()),2)<pow(.15,2));

            if(nph2[k]==2)
            {
              Ex2En2.Fill(vz.E(),vsum.E());
              En2.Fill(vsum.E(),spinindex);
            };
          };
          //* end nph2[k]>1

          if(Centercut)
          {
            if( vz.E()>0){
              z=(va.E()-vb.E())/vz.E();
              ZvsE1.Fill(z,vz.E());
              BnchidvsEc.Fill(Bunchid7bit,vz.E());
              if(fabs( vz.Mag()-.135)<.05)
              {
                ZvsEPic.Fill(z,vz.E());
              };
              if(fabs( vz.Mag()-.55)<.07)
              {
                BnchidvsEEtac.Fill(Bunchid7bit,vz.E());
                ZvsEEtac.Fill(z,vz.E());
              };
            };
          };


          //*beg nph2[k]==2

          if(nph2[k]==2)
          {
            Float_t WidthEta=.07;
            TLorentzVector va,vb;
            va=vec[pph2[k][0]];
            vb=vec[pph2[k][1]];
            TVector3 Va(va.X(),va.Y(),va.Z());
            TVector3 Vb(vb.X(),vb.Y(),vb.Z());
            Double_t fphi=(Va.Cross(Vb)).Phi();
            Double_t zu=vec[pph2[k][0]].E()-vec[pph2[k][1]].E();
            Double_t zb=vec[pph2[k][0]].E()+vec[pph2[k][1]].E();
            Bool_t zcut=fabs(zu/zb)<.85;
            ph[k]=vec[pph2[k][0]]+vec[pph2[k][1]];
            if(!zcut)continue;

            if(k==5 )
            {

              Float_t m5=ph[k].Mag();
              for(Int_t k5=0;k5<2;k5++)
              {
                TVector3 v3d=vec[pph2[k][k5]].Vect();
                TVector3 vv3d=v3d*(1./v3d.Z()*(*(p_Geom->ZFPD(2,1))));
                TVector3 v5=p_Geom->LocalXYZ(2,1,vv3d,true);
                T14x14.Fill(v5.X(),v5.Y(),m5);
              };

            };
            if(k<3){
              MvsTrig.Fill(ph[k].Mag(),BBcSums[0]);
              EvsTrig.Fill(ph[k].E(),BBcSums[0]);
              BBCsumEvsWeEn.Fill(BBcSums[2],BBcSums[1]);
              BBCsumEvseEast.Fill(BBcSums[1],ph[k].E());
              BBCsumWvseEast.Fill(BBcSums[2],ph[k].E());
              BBCtacEvsW.Fill(BBcSums[4],BBcSums[3]);
            }
            else {
            };
            if(Centercut)
            {
              Em2cvsEmm2_2Tr.Fill(ph[k].Mag(),vsum.Mag());
              if(vsum.E()>40)Em2cvsEmm2_2Tr40.Fill(ph[k].Mag(),vsum.Mag());
              if(vsum.E()>50)Em2cvsEmm2_2Tr50.Fill(ph[k].Mag(),vsum.Mag());
              if(vsum.E()>60)Em2cvsEmm2_2Tr60.Fill(ph[k].Mag(),vsum.Mag());
            };
            if( (ph[k].Mag()>.085)&&(ph[k].Mag()<.185))
            {

              YphiPi.Fill(ph[k].Eta(),tan(ph[k].Phi()));
              Float_t phiph=v1.Phi();
              phiph=fmod(phiph+2.5*Pi,2*Pi)-Pi/2.;
              if(fabs(ph[k].E()-45)<5.)YPhiPi4050.Fill(ph[k].Eta(),phiph);
              if(fabs(ph[k].E()-55)<5.)YPhiPi5060.Fill(ph[k].Eta(),phiph);
              if(fabs(ph[k].E()-65)<5.)YPhiPi6070.Fill(ph[k].Eta(),phiph);
              if(fcut1&&fcut2){
                if(fabs(vsum.E()-45)<5.)YPhiPi4050f.Fill(ph[k].Eta(),phiph);
                if(fabs(vsum.E()-55)<5.)YPhiPi5060f.Fill(ph[k].Eta(),phiph);
                if(fabs(vsum.E()-65)<5.)YPhiPi6070f.Fill(ph[k].Eta(),phiph);
              };
              ExFill.Fill(FillNumber,ph[k].E());
              Ex.Fill(ph[k].E(),spinindex);
              Exw.Fill(vsum.E(),spinindex);

              if(Centercut){
                ExPic.Fill(ph[k].E(),spinindex);
                ExwPic.Fill(vsum.E(),spinindex);
                if(fabs(zu/zb)<.4)ExwPicZ5.Fill(vsum.E(),spinindex);
                if(fcut1&&fcut2)
                {
                  ExPicf.Fill(ph[k].E(),spinindex);
                  ExwPicf.Fill(vsum.E(),spinindex);
                };
              };
            };
            if( (fabs(ph[k].Mag()-.55)<WidthEta))
              //			  &&fabs(vsum.E()-ph[k].E())<.2*ph[k].E())
            {
              ExEta.Fill(ph[k].E(),spinindex);
              ExwEta.Fill(vsum.E(),spinindex);

              if(Centercut)
              {
                ExEtac.Fill(ph[k].E(),spinindex);
                ExwEtac.Fill(vsum.E(),spinindex);
                if(fabs(zu/zb)<.4)ExwEtacZ5.Fill(vsum.E(),spinindex);

                if(fcut1&&fcut2)
                {
                  ExEtacf.Fill(ph[k].E(),spinindex);
                  ExwEtacf.Fill(vsum.E(),spinindex);
                };

              };

              BBCsumEvsWetac.Fill(BBcSums[2],BBcSums[1]);
              BBCtacEvsWetac.Fill(BBcSums[2],BBcSums[1]);
            };
            if(vec[plv[k][0]].E()>1 )
            {
              if(ph[k].E()>0 && Centercut){
                E2phvsEm2wc.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                    ph[k].E());	
                if(fabs(zu/zb)<.4)
                {
                  E2phvsEm2wcZ5.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                      ph[k].E());	
                }
                else
                {
                  E2phvsEm2wcZg5.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                      ph[k].E());	
                };
                E2phvsEm2c.Fill(ph[k].Mag(),
                    ph[k].E());	
                if(fcut1&&fcut2)
                {
                  E2phvsEm2fc.Fill(ph[k].Mag(),
                      ph[k].E());
                  E2phvsEm2wfc.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                      ph[k].E());	


                };

              };

              E2phvsEm2w.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                  ph[k].E());
              E2phvsEm2.Fill(ph[k].Mag(),
                  ph[k].E());
              if(fcut1&&fcut2)
              {
                E2phvsEm2f.Fill(ph[k].Mag(),
                    ph[k].E());
                E2phvsEm2wf.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),
                    ph[k].E());	

              };
              Float_t zz=fabs( ( vec[pph2[k][0]].E()-vec[pph2[k][1]].E() )/(ph[k].E()) );
              if(ph[k].E()>2. && zz<.85)
              {
                if(k==1){
                  Em2EN.Fill(ph[k].Mag(),cellN);
                  Em2ENwf.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),cellN);
                };
                if(k==2){
                  Em2ES.Fill(ph[k].Mag(),cellN);
                  Em2ESwf.Fill(ph[k].Mag()/ph[k].E()*vsum.E(),cellN);
                };
                Em1.Fill(ph[k].Mag(),spinindex);
                Float_t fudge=1.;
                Bool_t Ein8=false;
                Bool_t Ein7=false;
                if(nlv[7]>0)
                {
                  if(vec[plv[7][0]].E()>20)Ein7=true;
                };

                if(nlv[8]>0)
                {
                  if(vec[plv[8][0]].E()>20)Ein8=true;
                };
                Bool_t veto=false;
                if(k==5 )
                {
                  veto=true;
                  TVector3 Sa=Va *(1./Va.Z()*(*p_Geom->ZFPD(2,1)));
                  TVector3 Sb=Vb *(1./Vb.Z()*(*p_Geom->ZFPD(2,1)));
                  TVector3 v51=p_Geom->LocalXYZ(2,1,Sa,true);
                  Int_t ixa1=(Int_t) v51.X();
                  Int_t iyb1=(Int_t) v51.Y();
                  if(ixa1>=0 && iyb1>=0 && ixa1<14 && iyb1<14)
                  {
                    Float_t rfa=Rgainorigcorr->GetValue(2,1,iyb1,ixa1);
                    Float_t rfb=Rgaincorr->GetValue(2,1,iyb1,ixa1);
                    if(rfa>0)fudge=fudge*sqrt(rfb/rfa);
                  };
                  TVector3 v52=p_Geom->LocalXYZ(2,1,Sb,true);
                  Int_t ixa2=(Int_t) v52.X();
                  Int_t iyb2=(Int_t) v52.Y();
                  if(ixa1==4 && iyb1==4){d5_44.Fill(ixa2,iyb2,fudge*ph[k].Mag());e5_44.Fill(ixa2,iyb2,ph[k].E()*fudge);};
                  if(ixa2==4 && iyb2==4){d5_44.Fill(ixa1,iyb1,ph[k].Mag()*fudge);e5_44.Fill(ixa1,iyb1,fudge*ph[k].E());};

                  if(ixa1==4 && iyb1==4 && ixa2==2 && iyb2==4)
                  {
                    m5_42_44.Fill(ph[k].Mag(),ph[k].E());
                  };
                  if(ixa1==2 && iyb1==4 && ixa2==4 && iyb2==4)m5_42_44.Fill(ph[k].Mag(),ph[k].E());
                  if(ixa2>=0 && iyb2>=0 && ixa2<14 && iyb2<14)
                  {
                    Float_t rfa=Rgainorigcorr->GetValue(2,1,iyb2,ixa2);
                    Float_t rfb=Rgaincorr->GetValue(2,1,iyb2,ixa2);
                    if(rfa>0)fudge=fudge*sqrt(rfb/rfa);
                  };
                  sepvse.Fill((v52-v51).Mag()*fudge,ph[k].E());
                  if( (v52-v51).Mag() >.5 && ((v52-v51).Mag()<16.) &&ph[k].E()*fudge>3)
                  {
                    veto=false;
                    TP14x14.Fill(v51.X(),v51.Y(),ph[k].Mag()*fudge);
                    TP14x14.Fill(v52.X(),v52.Y(),ph[k].Mag()*fudge);
                  };
                  if(nph2[7]==2)
                  {
                    m2comp.Fill( 
                        (vec[pph2[7][1]]+vec[pph2[7][0]]).Mag(),
                        ph[k].Mag()*fudge);
                  }
                  else if(nph2[7]<2)
                  { 
                    m2comp.Fill(0., ph[k].Mag()*fudge);
                  }
                  else
                  { 
                    m2comp.Fill(2.,ph[k].Mag()*fudge);
                  };
                };
                if(k==6 )
                {
                  veto=true;
                  TVector3 Sa=Va *(1./Va.Z()*(*p_Geom->ZFPD(2,2)));
                  TVector3 Sb=Vb *(1./Vb.Z()*(*p_Geom->ZFPD(2,2)));
                  TVector3 v51=p_Geom->LocalXYZ(2,2,Sa,true);
                  Int_t ixa1=(Int_t) v51.X();
                  Int_t iyb1=(Int_t) v51.Y();
                  if(ixa1>=0 && iyb1>=0 && ixa1<14 && iyb1<14)
                  {
                    Float_t rfa=Rgainorigcorr->GetValue(2,2,iyb1,ixa1);
                    Float_t rfb=Rgaincorr->GetValue(2,2,iyb1,ixa1);
                    if(rfa>0)fudge=fudge*sqrt(rfb/rfa);
                  };
                  TVector3 v52=p_Geom->LocalXYZ(2,2,Sb,true);
                  Int_t ixa2=(Int_t) v52.X();
                  Int_t iyb2=(Int_t) v52.Y();
                  if(ixa2>=0 && iyb2>=0 && ixa2<14 && iyb2<14)
                  {
                    Float_t rfa=Rgainorigcorr->GetValue(2,2,iyb2,ixa2);
                    Float_t rfb=Rgaincorr->GetValue(2,2,iyb2,ixa2);
                    if(rfa>0)fudge=fudge*sqrt(rfb/rfa);
                  };
                  sepvse.Fill((v52-v51).Mag(),ph[k].E());
                  if( (v52-v51).Mag() >.5 && ((v52-v51).Mag()<16.) &&ph[k].E()*fudge>3 )
                  {
                    veto=false;
                    TQ14x14.Fill(v51.X(),v51.Y(),ph[k].Mag()*fudge);
                    TQ14x14.Fill(v52.X(),v52.Y(),ph[k].Mag()*fudge);
                  };
                };
                if(!veto)
                {
                  Em2Fill.Fill(FillNumber,ph[k].Mag());				  
                  Em2.Fill(ph[k].Mag()*fudge,spinindex);
                  Em2w.Fill(ph[k].Mag()/ph[k].E()*vsum.E()*fudge
                      ,spinindex);

                  if(vsum.E()>40)Em2w40.Fill(ph[k].Mag()/ph[k].E()*vsum.E()
                      ,spinindex);
                  if(vsum.E()>45) Em2w45.Fill(ph[k].Mag()/ph[k].E()*vsum.E()
                      ,spinindex);
                  if(vsum.E()>50)Em2w50.Fill(ph[k].Mag()/ph[k].E()*vsum.E()
                      ,spinindex);
                  if(ph[k].E()>30.)Em230.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>40.)Em240.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>50.)Em250.Fill(ph[k].Mag(),spinindex);
                };

                if(fabs(ph[k].Mag()-.3)<.1)
                {
                  ZMid.Fill(fabs(zu/zb),spinindex);

                  if(ph[k].E()>40.)ZMid40.Fill(fabs(zu/zb),spinindex);
                  if(ph[k].E()>40.&& Centercut)ZMid40c.Fill(fabs(zu/zb),spinindex);
                  if(ph[k].E()>50.)ZMid50.Fill(fabs(zu/zb),spinindex);
                  if(ph[k].Mag()<.3)
                  {
                    ZMida.Fill(fabs(zu/zb),spinindex);

                    if(ph[k].E()>40.)ZMid40a.Fill(fabs(zu/zb),spinindex);
                    if(ph[k].E()>50.)ZMid50a.Fill(fabs(zu/zb),spinindex);				      
                    if(ph[k].E()>40.&&Centercut)ZMid40ac.Fill(fabs(zu/zb),spinindex);

                  }
                  else
                  {
                    ZMidb.Fill(fabs(zu/zb),spinindex);

                    if(ph[k].E()>40.)ZMid40b.Fill(fabs(zu/zb),spinindex);
                    if(ph[k].E()>50.)ZMid50b.Fill(fabs(zu/zb),spinindex);
                    if(ph[k].E()>40.&& Centercut)ZMid40bc.Fill(fabs(zu/zb),spinindex);
                  };

                };
                if(Centercut){

                  Em1c.Fill(ph[k].Mag(),spinindex);
                  Em2c.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>40.)Em2c40.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>50.)Em2c50.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>55.)Em2c55.Fill(ph[k].Mag(),spinindex);
                  if(ph[k].E()>40.)Em2c40w.Fill(ph[k].Mag()*vsum.E()/ph[k].E(),spinindex);
                  if(ph[k].E()>50.)Em2c50w.Fill(ph[k].Mag()*vsum.E()/ph[k].E(),spinindex);
                  if(ph[k].E()>55.)Em2c55w.Fill(ph[k].Mag()*vsum.E()/ph[k].E(),spinindex);

                  if( (ph[k].Mag()>.085)&&(ph[k].Mag()<.185))
                  {
                    YTphiPic.Fill(ph[k].Eta(),tan(ph[k].Phi()));
                    if(zb!=0)
                    {
                      ZCPi.Fill(fabs(zu/zb),spinindex);
                      if(ph[k].E()>50)ZCPi50.Fill(fabs(zu/zb),spinindex);
                    };
                    Phi2GCPi.Fill(fmod(fphi,Pi),spinindex);
                    if(k==1){

                      EP1Pi0_XY.Fill((V3_1.X()-xDet)/cellwidth,
                          (V3_1.Y()-yDet)/cellwidth);
                      EP1Pi0_XY.Fill((V3_2.X()-xDet)/cellwidth,
                          (V3_2.Y()-yDet)/cellwidth);
                    };
                    if(k==2){
                      EP2Pi0_XY.Fill((V3_1.X()-xDet)/cellwidth,
                          (V3_1.Y()-yDet)/cellwidth);
                      EP2Pi0_XY.Fill((V3_2.X()-xDet)/cellwidth,
                          (V3_2.Y()-yDet)/cellwidth);
                    };
                    BBCsumEPi0.Fill(BBcSums[1],spinindex);
                    BBCsumWPi0.Fill(BBcSums[2],spinindex);

                  };
                  if( (fabs(ph[k].Mag()-.55)<WidthEta)){
                    YTphiEtac.Fill(ph[k].Eta(),tan(ph[k].Phi()));
                    if(zb!=0)
                    {
                      ZCEta.Fill(fabs(zu/zb),spinindex);
                      if(ph[k].E()>50)ZCEta50.Fill(fabs(zu/zb),spinindex);
                    };
                    Phi2GCEta.Fill(fmod(fphi,Pi),spinindex);

                    if(k==1){
                      EP1Eta_XY.Fill((V3_1.X()-xDet)/cellwidth,
                          (V3_1.Y()-yDet)/cellwidth);
                      EP1Eta_XY.Fill((V3_2.X()-xDet)/cellwidth,
                          (V3_2.Y()-yDet)/cellwidth);};
                      if(k==2){
                        EP2Eta_XY.Fill((V3_1.X()-xDet)/cellwidth,
                            (V3_1.Y()-yDet)/cellwidth);
                        EP2Eta_XY.Fill((V3_2.X()-xDet)/cellwidth,
                            (V3_2.Y()-yDet)/cellwidth);
                      };
                      BBCsumEEta.Fill(BBcSums[1],spinindex);
                      BBCsumWEta.Fill(BBcSums[2],spinindex);
                      if(ph[k].E()>55){
                        BBCsumEEta55.Fill(BBcSums[1],spinindex);
                        BBCsumWEta55.Fill(BBcSums[2],spinindex);
                      };
                  };
                };
              };
              if( (ph[k].Mag()>.085)&&(ph[k].Mag()<.185)){

                Phi2G.Fill(fmod(fphi,Pi),spinindex);
                if(zb!=0)Z.Fill(fabs(zu/zb),spinindex);
              } else 
                if( fabs(ph[k].Mag()-.55)<WidthEta){
                  YphiEta.Fill(ph[k].Eta(),tan(ph[k].Phi()));
                };
            };
          };
          if(Centercut)Em1cvsEn1_2Tr.Fill(vsum.Mag(),vsum.E());
          //*end  nph2[k]==2
          //*beg nph2[k]==1
          if(nph2[k]==1)
          {
            AbsPt1.Fill(fabs(vec[pph2[k][0]].Px() ),spinindex);
            if(k==1 && nph2[2]>0){
              Pt1Pt2.Fill(vec[pph2[1][0]].Px(),vec[pph2[2][0]].Px());
              TLorentzVector vm2p=vec[pph2[1][0]]+vec[pph2[2][0]];
              m2p.Fill(vm2p.Mag(),spinindex);

            };
            if(k==7 && nph2[8]>0 && nph2[5]==0 && nph2[6]==0){
              Pt1Pt2.Fill(vec[pph2[7][0]].Px(),vec[pph2[8][0]].Px());
              TLorentzVector vm2p=vec[pph2[7][0]]+vec[pph2[8][0]];
              m2p.Fill(vm2p.Mag(),spinindex);
            };
            if((k==7 && loEsum[5]+hiEsum[5]<2.) || 
                (k==8 && loEsum[6]+hiEsum[6]<2.) || (k!=7 && k!=8)
              )
            {
              ExIso.Fill(vec[pph2[k][0]].E(),spinindex);
            }
            else
            {
              ExNIso.Fill(vec[pph2[k][0]].E(),spinindex);
            };
            if(Center1)
            {
              if(50>vsum.E() && vsum.E()>40)
              {
                Em1c401Tr.Fill(vsum.Mag(),spinindex);
              };
              if(60>vsum.E() && vsum.E()>50)
              {
                Em1c501Tr.Fill(vsum.Mag(),spinindex);
              };
              if(70>vsum.E() && vsum.E()>60)
              {
                Em1c601Tr.Fill(vsum.Mag(),spinindex);
              };
              Exc1Trw.Fill(vsum.E(),spinindex);
              Em1cvsEn1_1Tr.Fill(vsum.Mag(),vsum.E());
            };
          };
          //*end nph2[k]==1
        };
        //*end nph2[k]>0

      };
    //*end kloop
    };
    //*end trigsel
  };
  //*end Entries
  NphvsEsum.Write();
  NphvsMsum.Write();
  Em1c401Tr.Write();
  Em1c501Tr.Write();
  Em1c601Tr.Write();
  Exc1Trw.Write();
  Em1cvsEn1_1Tr.Write();
  Em1cvsEn1_2Tr.Write();
  Em2cvsEmm2_2Tr.Write();
  Em2cvsEmm2_2Tr40.Write();
  Em2cvsEmm2_2Tr50.Write();
  Em2cvsEmm2_2Tr60.Write();


  ZMid40a.Write();
  ZMid40ac.Write();
  ZMid50a.Write();

  ZMid40b.Write();
  ZMid40bc.Write();
  ZMid50b.Write();


  ZMid.Write();
  ZMid40.Write();
  ZMid40c.Write();
  ZMid50.Write();

  EcompA.Write();
  EcompB.Write();
  ExPicf.Write();
  ExEtacf.Write();
  ExwPicf.Write();
  ExwEtacf.Write();
  Em240.Write();
  Em250.Write();
  Em230.Write();
  Em2w40.Write();
  Em2w45.Write();
  Em2w50.Write();
  Em2c40.Write();
  Em2c50.Write();
  Em2c55.Write();
  Em2c40w.Write();
  Em2c50w.Write();
  Em2c55w.Write();

  ESumFill.Write();
  Em2Fill.Write();
  ExFill.Write();
  m2comp.Write();
  e5_44.Write();
  d5_44.Write();
  m5_42_44.Write();
  sepvse.Write();
  mv7.Write();
  mv5.Write();
  mv57.Write();
  Y5_7.Write();
  ExW.Write();
  ExIso.Write();
  ExNIso.Write();
  WNinvsWSin.Write();
  WNinvsWSout.Write();
  WNinvsWNout.Write();
  WSinvsWSout.Write();
  WSinvsWNout.Write();
  Miosp.Write();
  Mio.Write();
  MioL1.Write();
  MioL2.Write();
  MioL3.Write();
  YFillN50.Write();
  YFillS50.Write();
  Em2EN.Write();
  Em2ES.Write();
  Em2ENwf.Write();
  Em2ESwf.Write();
  YPhiPi4050.Write();
  YPhiPi5060.Write();
  YPhiPi6070.Write();
  YPhiPi4050f.Write();
  YPhiPi5060f.Write();
  YPhiPi6070f.Write();
  YPhi1p4050.Write();
  YPhi1p5060.Write();

  YPhi1p20f.Write();
  E2phvsEm2fc.Write();
  E2phvsEm2f.Write();
  E2phvsEm2wf.Write();
  E2phvsEm2wfc.Write();

  E2phvsEm2c.Write();
  E2phvsEm2.Write();
  E2phvsEm2wc.Write();
  E2phvsEm2wcZ5.Write();
  E2phvsEm2wcZg5.Write();
  E2phvsEm2w.Write();

  Em2w.Write();
  BBCsumEvseEast0.Write();
  BBCsumWvseEast0.Write();
  BBCsumEg150.Write();
  BBCsumWg150.Write();
  BBCsumEg1.Write();
  BBCsumWg1.Write();
  BBCsumEg0.Write();
  BBCsumWg0.Write();
  BBCsumEEta55.Write();
  BBCsumWEta55.Write();  
  BBCsumEvseEast.Write();
  BBCsumWvseEast.Write();
  BBCsumEvsWeEn.Write();
  BBCsumEvsWwEn.Write();
  BBCsumEEta.Write();
  BBCsumWEta.Write();
  BBCsumEPi0.Write();
  BBCsumWPi0.Write();
  BnchidvsEW.Write();
  BnchidvsEE.Write();
  BnchidvsEEtac.Write();
  BnchidvsEc.Write();
  hFill.Write();
  RNum.Write();
  Bnchid.Write();
  Ex.Write();
  Exw.Write();
  ExwPic.Write();
  ExwPicZ5.Write();
  ExwEtac.Write();
  ExwEtacZ5.Write();
  ExwEta.Write();
  Ex1.Write();
  ExEta.Write();
  ExEtac.Write();
  Ex2En2.Write();
  En1.Write();
  En2.Write();
  ExPic.Write();
  Em1.Write();
  Em2.Write();
  Em2c.Write();
  Em1c.Write();
  Emm1.Write();
  Emm2.Write();
  Z.Write();
  Phi2G.Write();
  Nphots.Write();
  YphiPi.Write();
  YphiEta.Write();
  Phi2GCEta.Write();
  ZCEta.Write();
  ZCEta50.Write();
  Phi2GCPi.Write();
  ZCPi.Write();
  ZCPi50.Write();
  YTphiEtac.Write();
  YTphiPic.Write();
  YPhi1p20.Write();
  Pt1Pt2.Write();
  AbsPt1.Write();
  m2p.Write();
  ZvsE1.Write();
  ZvsEPic.Write();
  ZvsEEtac.Write();
  EP1_XY.Write();
  EP2_XY.Write();
  EP1Pi0_XY.Write();
  EP2Pi0_XY.Write();
  EP1Eta_XY.Write();
  EP2Eta_XY.Write();
  MvsTrig.Write();
  EvsTrig.Write();
  BBCsumEvsW.Write();
  BBCtacEvsW.Write();
  BBCsumEvsWetac.Write();
  BBCtacEvsWetac.Write();
  T14x14.Write();

  SP14x14.Write();  
  SQ14x14.Write();
  TP14x14.Write();  
  TQ14x14.Write();
  EP14x14.Write();  
  EQ14x14.Write();
  T014x14.Write();
  T7x7.Write();
  //  new_P_File->cd();
  //  new_p_OF->Write();
  return 0;

};

Int_t AnalTools::read1a(char* filelist,Int_t set,FilesSet* p_files)
{

  Int_t DetMap[9]={0,1,2,0,0,5,6,3,4};
  Geom* p_geom=new Geom(p_files);
  pout=new poutTree(filelist);  
  pout->EnableEdepCorr=0;//not Enable Edep corrections

  for(Int_t kk=1;kk<9;kk++)
  {
    pout->MinEnergy[kk]=3.;
  };
  Float_t lowlim=0;
  Float_t highlim=120;
  TH2F Em2("Em2","Em2",40,0.,1.,20,0,20);
  TH2F ExEtac("ExEtac","ExEtac",24,lowlim,highlim,30,0,30);
  TH2F ExPic("ExPic","ExPic",24,lowlim,highlim,30,0,30);

  Int_t nentries=0;
  TH2F E2phvsEm2("E2phvsEm2","E2phvsEm2",48,0,1.2,20,0,100);
  TH2F EtaPhi("EtaPhi","EtaPhi",   100,-4.2,-3.2,100,-.5,.5);
  TH2F EtaPhic("EtaPhic","EtaPhic",100,-4.2,-3.2,100,-.5,.5);
  TH2F Em2c("Em2c","Em2c",48,0.,1.2,30,0,30);  
  TH1F  Nphotons("Nphotons","Nphotons",20,0,20);
  typedef struct {
    Float_t M12;
    Float_t Phi12;
    Float_t Eta12;
    Float_t E12;
    Float_t Eta1;
    Float_t Phi1;
    Float_t E1;
    Float_t Eta2;
    Float_t Phi2;
    Float_t E2;
    Int_t Run;
    Int_t Event;
    Int_t SpinIndex;
    Int_t Bunchid7bit;
  } tp_eta_anal;
  tp_eta_anal tet;
  TFile Outputz("Outputz.root","recreate");
  TTree EtaAnal("EtaAnal","EtaAnal events");
  EtaAnal.Branch("M12",&tet.M12,"M12/F");
  EtaAnal.Branch("Phi12",&tet.Phi12,"Phi12/F");
  EtaAnal.Branch("Eta12",&tet.Eta12,"Eta12/F");
  EtaAnal.Branch("E12",&tet.E12,"E12/F");
  EtaAnal.Branch("Eta1",&tet.Eta1,"Eta1/F");
  EtaAnal.Branch("Phi1",&tet.Phi1,"Phi1/F");
  EtaAnal.Branch("E1",&tet.E1,"E1/F");
  EtaAnal.Branch("Eta2",&tet.Eta2,"Eta2/F");
  EtaAnal.Branch("Phi2",&tet.Phi2,"Phi2/F");
  EtaAnal.Branch("E2",&tet.E2,"E2/F");
  EtaAnal.Branch("Run",&tet.Run,"Run/I");
  EtaAnal.Branch("Event",&tet.Event,"Event/I");
  EtaAnal.Branch("SpinIndex",&tet.SpinIndex,"SpinIndex/I");
  EtaAnal.Branch("Bunchid7bit",&tet.Bunchid7bit,"Bunchid7bit/I");
  nentries=pout->nentries;

  Int_t TwoCnt=0;

  for(Int_t iev=0;iev<nentries;iev++)
  {
    if( (iev%100000) ==0)printf( "%d Events Read (%d two photon)\n ",iev,TwoCnt);
    pout->GetEntry(iev);  
    pout->AllToScratch(false); //make a scratch list of all hard 
    int chkcnt=0;
    TObjArray* pd[2];
    Int_t spin=pout->spin;
    for(int k=0;k<2;k++)
    {
      int kk=0;
      tet.Bunchid7bit=pout->Bunchid7bit;
      if(k==0)kk=1;
      if(k==1)kk=2;
      Int_t spn=100;
      if(spin>=0 && spin<4)spn=spin;
      Int_t spinindex=DetMap[kk]*4+spn;
      tet.SpinIndex=spinindex;
      pout->ClearScratch();
      pd[k]=pout->AddToScratch(kk);
      Int_t nph=pout->scratchlist->GetEntries();
      Nphotons.Fill(nph);

      if(nph==2)
      {

        LVec* p1=(LVec*) pd[k]->First();
        LVec* p2=(LVec*) pd[k]->After(p1);
        TLorentzVector ph=(TLorentzVector) (*p1+*p2);
        tet.E12=ph.E();
        tet.M12=ph.Mag();
        tet.Eta12=ph.PseudoRapidity();
        tet.Phi12=ph.Phi();
        tet.E1=p1->E();
        tet.E2=p2->E();
        tet.Phi1=p1->Phi();
        tet.Phi2=p2->Phi();
        tet.Eta1=p1->PseudoRapidity();
        tet.Eta2=p2->PseudoRapidity();
        tet.Event=pout->EventN;
        tet.Run=pout->Rnum;

        EtaPhi.Fill(ph.Eta(),tan(ph.Phi()));
        Bool_t Centercut=
          (pow(ph.Eta()+3.65,2)+pow(tan(ph.Phi()),2)<pow(.15,2));

        Em2.Fill(ph.Mag(),spinindex);

        if(Centercut)
        {
          if(tet.E12>30) EtaAnal.Fill();

          EtaPhic.Fill(ph.Eta(),tan(ph.Phi()));
          E2phvsEm2.Fill(ph.Mag(),ph.E()); 
          TwoCnt++;
          Em2c.Fill(ph.Mag(),spinindex);

          if(fabs(ph.Mag()-.55)<.07)
          {
            ExEtac.Fill(ph.E(),spinindex);
          };
          if(fabs(ph.Mag()-.135)<.05)
          {
            ExPic.Fill(ph.E(),spinindex);
          };
        };
      };
    };
  };
  EtaAnal.Write();
  ExPic.Write();
  Em2.Write();
  ExEtac.Write();
  Em2c.Write();
  E2phvsEm2.Write();
  EtaPhi.Write();
  EtaPhic.Write();
  Nphotons.Write();
};


Int_t AnalTools::Examp_script4()
{

  FilesSet* p_files=new FilesSet("../../",
      "run6_export/fpdped.txt",
      "run6_export/fpdgain.txt", 
      "run6_export/fpdcorr.txt",
      "run6_export/fill.txt",
      "Fake",
      "run6_export/spinpat",
      "run6_export/geom.txt");

  p_files->p_fpdrunped()->Directory="../../run6_export/runped";
  //p_files->p_fpdrunped->Directory="./";
  p_files->p_fpdlumi()->Directory="./";
  p_files->Print();


  dataSet d("./run*.root",p_files,"h111");

  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  //AnalTools* at = new AnalTools();
  // FpdMap is a class that figures out maps into data blocks
  FpdMap* pmap=new FpdMap(Rnum);

  TFile* p_OFile=new TFile("OFile.root","recreate");

  // Now return the 7x7  map matrix for EW=2 (west) and NSTP=2 (south)

  TMatrix map_pwn=pmap->GetMatrix(2,3);
  TMatrix map_pws=pmap->GetMatrix(2,4);
  TMatrix map_pen=pmap->GetMatrix(1,1);
  TMatrix map_pes=pmap->GetMatrix(1,2);
  TMatrix n199(14,14);
  n199=199;
  TMatrix map_pbwn=pmap->GetMatrix(2,1);
  TMatrix map_pbws=pmap->GetMatrix(2,2);
  map_pbwn=map_pbwn-n199;
  map_pbws=map_pbws-n199;


  Float_t  EsumWS,EsumWN,EsumES,EsumEN,EsumbWS,EsumbWN;
  EsumbWS=EsumbWN=0;
  Short_t EsumWS_s,EsumWN_s,EsumES_s,EsumEN_s;
  Int_t tpes[40];
  Int_t spin;
  Int_t nphotons;
  Float_t px[10];
  Float_t py[10];
  Float_t pz[10];
  Float_t pE[10];
  Float_t pxyzt[40];
  Int_t nwrds;
  //Define a Tree
  TTree* p_out = new TTree("p_out","Output Tree");

  p_out->Branch("spin",&(spin),"spin/I");
  p_out->Branch("nphotons",&(nphotons),"nphotons/I");
  p_out->Branch("br_nwrds",&(nwrds),"nwrds/I");
  p_out->Branch("br_types",&(tpes),"tpes[nwrds]/I");
  p_out->Branch("br_pxyzt",&(pxyzt),"pxyzt[nwrds]/F");

  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  // event loop
  for(Int_t i=0;i<nentries;i++)
  {
    Int_t nbytes=d.GetEntry(i); 
    pcnt++;
    if(pcnt>5000){printf("cnt=%d \n",i);pcnt=1;}; //print every 5000 events

    spin=(d.BlueSpin+1)+(d.YellowSpin+1)/2;
    if(d.kicked)spin+=29;
    if(d.BlueSpin*d.YellowSpin==0)spin=40;
    TMatrix tmp0wn= d.dMatrix(&(map_pwn),d.Fpdwns);
    TMatrix tmp0ws= d.dMatrix(&(map_pws),d.Fpdwns);
    TMatrix tmp0en= d.dMatrix(&(map_pen),d.Fpdens);
    TMatrix tmp0es= d.dMatrix(&(map_pes),d.Fpdens);
    TMatrix tmp0bws= d.dMatrix(&(map_pbws),d.fpdwadc);
    TMatrix tmp0bwn= d.dMatrix(&(map_pbwn),d.fpdwadc);

    TMatrixF energydata_2_4=d.Em(&tmp0ws,2,4);  // WS energy data for event 
    TMatrixF energydata_2_3=d.Em(&tmp0wn,2,3);  // WN energy data for event 
    TMatrixF energydata_2_2=d.Em(&tmp0bws,2,2);  // bWS energy data for event 
    TMatrixF energydata_2_1=d.Em(&tmp0bwn,2,1);  // bWN energy data for event 
    TMatrixF energydata_1_1=d.Em(&tmp0en,1,1);  // EN energy data for event 
    TMatrixF energydata_1_2=d.Em(&tmp0es,1,2);  // ES energy data for event 

    EsumWS=energydata_2_4.Sum(); // sum over 36 WS channels
    EsumWN=energydata_2_3.Sum(); // sum over 36 WN channels
    EsumES=energydata_1_2.Sum(); // sum over 49 ES channels
    EsumEN=energydata_1_1.Sum(); // sum over 49 EN channels

    //count the W channels with more than 0.25 GeV Energy

    Int_t highcnt4= NTowersAbove(&energydata_2_4,.25);
    Int_t highcnt3= NTowersAbove(&energydata_2_3,.25);
    EnergySum=0;
    nwrds=0;
    nphotons=0;

    TLorentzVector vecs[10];
    for(Int_t kj=0;kj<40;kj++)tpes[kj]=0; 
    if(highcnt3+highcnt4<40)
    {
      if(EsumWS>25){
        TLorentzVector lv(FourMom(&energydata_2_4,0,100,0,100,d.pGeom,2,4));
        tpes[nwrds*4]=8;
        vecs[nwrds]=lv;
        nwrds++;
        Yiqun recon(&energydata_2_4,d.pGeom,2,4);
        nphotons+=recon.NPh;
        if(recon.NPh==2)
        {
          TLorentzVector v4;
          v4=recon.mom(0)+recon.mom(1);
          tpes[nwrds*4]=108;
          vecs[nwrds]=v4;
          nwrds++;
        };
        //	TLorentzVector lv2(FourMom(&energydata_2_2,0,100,0,100,d.pGeom,2,2));
      };
      if(EsumWN>25){
        TLorentzVector lv(FourMom(&energydata_2_3,0,100,0,100,d.pGeom,2,3));
        tpes[nwrds*4]=7;
        vecs[nwrds]=lv;
        nwrds++;
        Yiqun recon(&energydata_2_3,d.pGeom,2,3);
        nphotons+=recon.NPh;
        if(recon.NPh==2)
        {
          TLorentzVector v4;
          v4=recon.mom(0)+recon.mom(1);
          tpes[nwrds*4]=107;
          vecs[nwrds]=v4;
          nwrds++;
        };

        //	TLorentzVector lv(FourMom(&energydata_2_1,0,100,0,100,d.pGeom,2,1));
      };

      if(EsumES>35){
        TLorentzVector lv(FourMom(&energydata_1_2,0,100,0,100,d.pGeom,1,2));
        tpes[nwrds*4]=2;
        vecs[nwrds]=lv;
        nwrds++;
        Yiqun recon(&energydata_1_2,d.pGeom,1,2);
        nphotons+=recon.NPh;
        if(recon.NPh==2)
        {
          TLorentzVector v4;
          v4=recon.mom(0)+recon.mom(1);
          tpes[nwrds*4]=102;
          vecs[nwrds]=v4;
          nwrds++;
        };
      };
      if(EsumEN>35){
        TLorentzVector lv(FourMom(&energydata_1_1,0,100,0,100,d.pGeom,1,1));
        tpes[nwrds*4]=1;
        vecs[nwrds]=lv;
        nwrds++;
        Yiqun recon(&energydata_1_1,d.pGeom,1,1);
        nphotons+=recon.NPh;	
        if(recon.NPh==2)
        {
          TLorentzVector v4;
          v4=recon.mom(0)+recon.mom(1);
          tpes[nwrds*4]=101;
          vecs[nwrds]=v4;
          nwrds++;
        };
      };
      Int_t wcnt=0;
      for(Int_t k=0;k<nwrds;k++)
      {
        pxyzt[wcnt]=vecs[k].Px();
        wcnt++;
        pxyzt[wcnt]=vecs[k].Py();
        wcnt++;
        pxyzt[wcnt]=vecs[k].Pz();
        wcnt++;
        pxyzt[wcnt]=vecs[k].E();
        wcnt++;
      };
      nwrds=nwrds*4;
      if(nwrds>0)p_out->Fill();
    };
  };

  //end event loop
  p_out->Print();
  p_out->Write();
  return 0;
};
Int_t AnalTools:: Select_ADC_QT(Int_t set,FilesSet* P_Files,char* hntName)
{
  printf("first here \n");
  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=false;
  AnalDet[0][1]=false;
  AnalDet[1][0]=true;
  AnalDet[1][1]=true;
  AnalDet[1][2]=true;
  AnalDet[1][3]=true;
  if(P_Files==0)
  {
    if(set==80)
    {
      //transverse run 8
      p_files=new FilesSet("/star/u/heppel/BatchDir/rootfms",
          "fpdped.txt",
          "fpdgain.txt", 
          "fpdcorr.txt",
          "fill.txt",
          "Fake",
          "fpd08/spinpat",
          "geom.txt",
          "qtmap.txt",
          "qtmap2pp.txt");


      p_files->Print();
      hnt_name="h111_08";

    }
    else  if(set==90)
    {
      hnt_name=hntName;
    }
    else
    {
      printf("no valid set selected for Select_scriptQT1\n");
      return 0;
    };
  }
  else
  {
    p_files=P_Files;
    hnt_name=hntName;
    std::cout<<"Defining "<<hnt_name<<"\n";
  };
  printf("got here \n");
  char hntp[20];
  strcpy(hntp,(const char*) hnt_name);
  printf("dataSet for %s \n",hntp);
  dataSet d("./run*.root",p_files,hntp);
  TFile Tr_File(File_adcTr,"recreate");
  TTree Tr_adc("Tr_adc","ADCRates");
  TFile adcout("adcout.root","recreate");

  Int_t ADC;
  Int_t rowADC;
  Int_t colADC;
  Int_t nstbADC;
  Int_t evtnum;
  Int_t runnum=d.CurrentRunNumber;
  Int_t segnum=d.CurrentSegNumber;
  Int_t sumADC=0;
  Int_t led=0;
  Int_t Bunchid7bit=0;
  Tr_adc.Branch("br_ADC",&ADC,"ADC/I");
  Tr_adc.Branch("br_rowADC",&rowADC,"rowADC/I");
  Tr_adc.Branch("br_colADC",&colADC,"colADC/I");
  Tr_adc.Branch("br_nstbADC",&nstbADC,"nstbADC/I");
  Tr_adc.Branch("br_runnum",&runnum,"runnum/I");
  Tr_adc.Branch("br_evtnum",&evtnum,"evtnum/I");
  Tr_adc.Branch("br_segnum",&segnum,"segnum/I");
  Tr_adc.Branch("br_led",&led,"led/I");
  Tr_adc.Branch("br_Bunchid7bit",&Bunchid7bit,"Bunchid7bit/I");
  Tr_adc.Branch("br_sumADC",&sumADC,"SumADC/I");


  //  d.pQt->EnableqtHist(); //make ADC histograms
  d.RFill->ListFills();
  TH2F ndatvssum("ndatvssum","ndatvssum",1000,0,1000,1000,0,50000);
  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  printf("d.CurrentRunNumber=%d \n",Rnum);
  //  Rnum=12010001; //force run number for now SH
  Bunchid7bit=d.bunchid7bit;

  // FpdMap is a class that figures out maps into data blocks
  //  FpdMap* pmap=new FpdMap(Rnum);


  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  // event loop
  Int_t skpcnt=0;
  UInt_t lastLED=0;
  for(Int_t i=0;i<nentries;i++)
  {
    Int_t nbytes=d.GetEntry(i);
    Bunchid7bit=d.bunchid7bit;  

    Bool_t LEDevent=false;
    Bool_t LEDecho=false;
    UInt_t LEDbit=d.lastdsm[7];
    Int_t bc=d.Bclo+d.Bchi*(Int_t)pow(2,32);
    if(Rnum>13000000)bc=d.Bc[1];

    if(Rnum>10000000)
    {
      if((LEDbit)%16>0)
      {
        LEDevent=true; 
      };
    }

    else
    {
      if((LEDbit/2048)%2==1)
      {
        LEDevent=true; 
        lastLED=bc;
        if(lastLED!=0 && ((bc-lastLED)%2000001<10) || (bc-lastLED)%2000001>1999997)
        {
          LEDevent=true;
        };  
      };
    };

    if(Rnum>12000000)//LKE 20110211
    {
      LEDevent=false;
      //	  printf("setting TrigBits for Rnum=%d \n",Rnum);
      TrigBits=d.lastdsm[5];	  
      LEDbit=(d.lastdsm[4] & 0x0001);
      if(LEDbit==1)
      {
        LEDevent=true; 
      };
      //TrigBits=(fpdebit<<14)+(dijetbit<<11)+(L2_jp_bit[2]<<10)+(L2_jp_bit[1]<<9)+(L2_jp_bit[0]<<8)
      //+(lrgclbit2<<7)+(lrgclbit1<<6)+(lrgclbit0<<5)+(L1_clbit[2][0]<<4)+(L1_clbit[1][0]<<3)+(L1_clbit[0][0]<<2)
      //+(fmshtbit1<<1)+fmshtbit0;
    };



    if(IdentifyLED==false)LEDevent=false;
    if(AllLed())LEDevent=true;

    //      printf("TrigBitsSelect=%x TrigBits=%x \n",TrigBitsSelect,TrigBits);
    if((TrigBitsSelect!=0) && (TrigBitsSelect&TrigBits)==0)continue;
    if( (TrigBitsVeto&TrigBits)!=0)continue;
    if(!d.decodeQT()){printf("qt Error \n"); continue;};
    TMatrix* madc;
    led=0;
    Float_t AdcSum=0;
    for(int iadc=0;iadc<4;iadc++)
    {
      madc=(TMatrix*) d.pQt->qtMat->At(iadc);
      AdcSum+=madc->Sum();
    }
    //      if(AdcSum>1000)printf("ev=%d AdcSum=%f  \n",i,AdcSum);
    if(LEDevent)led=1;
    if(AdcSum>8000)
    {
      LEDevent=true;
      led=1;
    }
    if (AdcSum>100 && abs((Int_t)(bc-lastLED))<2)
    {
      LEDevent=true;
      LEDecho=true;
      led=2;
    }
    lastLED=bc;
    sumADC=AdcSum;
    if(!LEDevent && OnlyLED )continue;
    evtnum=i;
    for(int iadc=0;iadc<4;iadc++)
    {
      nstbADC=iadc+1;
      madc=(TMatrix*) d.pQt->qtMat->At(iadc);
      int nrow=madc->GetNrows();
      int ncol=madc->GetNcols();
      for(int row=0;row<nrow;row++)
      {
        rowADC=row;
        for(int col=0;col<ncol;col++)
        {
          colADC=col;
          ADC=(*madc)(row,col);
          if(ADC>0)
          {
            Tr_adc.Fill();
          };
        };
      };
    };
    /*
       if(d.pQt->FmsSumAdc>1300)
       {
       printf("LEDbit=%X sumadc=%d \n",d.lastdsm[7],d.pQt->FmsSumAdc); 
       };
       if(!LEDevent)ndatvssum.Fill(d.pQt->ndata,d.pQt->sumadc);
       */
    pcnt++;
    if(pcnt>100){printf("cnt=%d \n",i);pcnt=1;}; //print every 5000 events
  };
  printf("Done with %d  event loop \n",nentries);
  //  ndatvssum.Write();
  printf("ndatvsssu.Write() done \n");
  Tr_File.cd();
  Tr_adc.Write();
  adcout.cd();
  printf("Tr_adc.Write() done \n");
  /*
     d.pQt->qtHist[0]->Write();
     d.pQt->qtHist[1]->Write();
     d.pQt->qtHist[2]->Write();
     d.pQt->qtHist[3]->Write();
     */
  TH2F* qtHist0[4];
  TH2F* qtHist1[4];
  qtHist0[0]=new TH2F("qtHistLN","qtHistLN",578,0,578,4096,0,4096);
  qtHist0[1]=new TH2F("qtHistLS","qtHistLS",578,0,578,4096,0,4096);
  qtHist0[2]=new TH2F("qtHistSN","qtHistSN",288,0,288,4096,0,4096);
  qtHist0[3]=new TH2F("qtHistSS","qtHistSS",288,0,288,4096,0,4096);
  qtHist1[0]=new TH2F("qtHistLedLN","qtHistLedLN",578,0,578,4096,0,4096);
  qtHist1[1]=new TH2F("qtHistLedLS","qtHistLedLS",578,0,578,4096,0,4096);
  qtHist1[2]=new TH2F("qtHistLedSN","qtHistLedSN",288,0,288,4096,0,4096);
  qtHist1[3]=new TH2F("qtHistLedSS","qtHistLedSS",288,0,288,4096,0,4096);
  printf("8 qthists defined \n");
  Int_t NCOL[4]={17,17,12,12};
  for(int j=0;j<Tr_adc.GetEntries();j++)
  {
    Int_t Nb=Tr_adc.GetEntry(j);
    if( (j%1000 == 0 ))
      printf("read Tr_adc j=%d  nbytes = %d \n",j,Nb);

    Int_t nstb=nstbADC;
    if(nstb>0 && nstb<5)
    {
      Int_t chan0=(rowADC)*NCOL[nstb-1]+colADC;
      if(led==0)
      {
        qtHist0[nstb-1]->Fill(chan0,ADC);
      }
      else
      {
        qtHist1[nstb-1]->Fill(chan0,ADC);
      };
    };
  };
  for(int nh=0;nh<4;nh++)
  {
    qtHist0[nh]->Write();
    qtHist1[nh]->Write();
  };
  printf("d.pQt->qtHist[3]->.Write() done \n");
  printf("ready to return \n");
  //end event loop
  return 0;
};

Int_t AnalTools:: Select_ADC_QT2(Int_t set,FilesSet* P_Files,char* hntName)
{
  printf("first here \n");
  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=false;
  AnalDet[0][1]=false;
  AnalDet[1][0]=true;
  AnalDet[1][1]=true;
  AnalDet[1][2]=true;
  AnalDet[1][3]=true;
  if(P_Files==0)
  {
    if(set==80)
    {
      //transverse run 8
      p_files=new FilesSet("/star/u/heppel/BatchDir/rootfms",
          "fpdped.txt",
          "fpdgain.txt", 
          "fpdcorr.txt",
          "fill.txt",
          "Fake",
          "fpd08/spinpat",
          "geom.txt",
          "qtmap.txt",
          "qtmap2pp.txt");


      p_files->Print();
      hnt_name="h111_08";

    }
    else  if(set==90)
    {
      hnt_name=hntName;
    }
    else
    {
      printf("no valid set selected for Select_scriptQT1\n");
      return 0;
    };
  }
  else
  {
    p_files=P_Files;
    hnt_name=hntName;
    std::cout<<"Defining "<<hnt_name<<"\n";
  };
  printf("got here \n");
  char hntp[20];
  strcpy(hntp,(const char*) hnt_name);
  printf("dataSet for %s \n",hntp);
  dataSet d("./run*.root",p_files,hntp);
  TFile Tr_File(File_adcTr,"recreate");
  TTree Tr_adc("Tr_adc","ADCRates");
  TFile adcout("adcout.root","recreate");

  Int_t ADC;
  Int_t rowADC;
  Int_t colADC;
  Int_t nstbADC;
  Int_t evtnum;
  Int_t runnum=d.CurrentRunNumber;
  Int_t segnum=d.CurrentSegNumber;

  Int_t led=0;
  Int_t Bunchid7bit=0;
  Int_t ClusterCount=0;
  Int_t ClusterOrder=0;
  Float_t CellE=0.;
  Float_t ClusterE=0;
  Float_t minang=0;
  Float_t minE=0;
  Float_t minM=0;
  Float_t HiCluE=0;
  Float_t NextHiCluE=0;
  Float_t dToHiTower=0;
  Float_t PseudoR=0;
  Float_t Edet=0;
  Float_t ADCDet=0;
  cout<<"set up branches\n";
  Tr_adc.Branch("br_ADC",&ADC,"ADC/I");
  Tr_adc.Branch("br_rowADC",&rowADC,"rowADC/I");
  Tr_adc.Branch("br_colADC",&colADC,"colADC/I");
  Tr_adc.Branch("br_nstbADC",&nstbADC,"nstbADC/I");
  Tr_adc.Branch("br_runnum",&runnum,"runnum/I");
  Tr_adc.Branch("br_evtnum",&evtnum,"evtnum/I");
  Tr_adc.Branch("br_segnum",&segnum,"segnum/I");
  Tr_adc.Branch("br_led",&led,"led/I");
  Tr_adc.Branch("br_Bunchid7bit",&Bunchid7bit,"Bunchid7bit/I");
  Tr_adc.Branch("br_CellE",&CellE,"CellE/F");
  Tr_adc.Branch("br_ClusterCount",&ClusterCount,"ClusterCount/I");
  Tr_adc.Branch("br_ClusterE",&ClusterE,"ClusterE/F");
  Tr_adc.Branch("br_ClusterOrder",&ClusterOrder,"ClusterOrder/I");
  Tr_adc.Branch("br_minang",&minang,"minang/F");
  Tr_adc.Branch("br_minE",&minE,"minE/F");
  Tr_adc.Branch("br_minM",&minM,"minM/F");
  Tr_adc.Branch("br_HiCluE",&HiCluE,"HiCluE/F");
  Tr_adc.Branch("br_NextHiCluE",&NextHiCluE,"NextHiCluE/F");
  Tr_adc.Branch("br_dToHiTower",&dToHiTower,"dToHiTower/F");
  Tr_adc.Branch("br_PseudoR",&PseudoR,"PseudoR/F");
  Tr_adc.Branch("br_EDet",&Edet,"EDet/F");
  Tr_adc.Branch("br_ADCDet",&ADCDet,"ADCDet/F");
  Tr_adc.Branch("br_TrigBits",&TrigBits,"TrigBits/I");
  cout<<"setup done\n";
  Yiqun* recon[4]={0,0,0,0};
  //  d.pQt->EnableqtHist(); //make ADC histograms
  d.RFill->ListFills();
  TH2F ndatvssum("ndatvssum","ndatvssum",1000,0,1000,1000,0,50000);
  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  printf("d.CurrentRunNumber=%d \n",Rnum);
  //  Rnum=12010001; //force run number for now SH
  Bunchid7bit=d.bunchid7bit;

  // FpdMap is a class that figures out maps into data blocks
  //  FpdMap* pmap=new FpdMap(Rnum);


  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  // event loop
  Int_t skpcnt=0;
  UInt_t lastLED=0;

  TMatrix* p_adc[4]={0,0,0,0};
  TMatrix* p_Emat[4]={0,0,0,0};
  TObjArray* clusts=0;
  TObjArray* threeV=0;
  for(Int_t i=0;i<nentries;i++)
  {
    Edet=0;
    ADCDet=0;
    Int_t nbytes=d.GetEntry(i);
    Bunchid7bit=d.bunchid7bit;  
    Bool_t LEDevent=false;
    Bool_t LEDecho=false;
    UInt_t LEDbit=d.lastdsm[7];
    Int_t bc=d.Bclo+d.Bchi*(Int_t)pow(2,32);
    if(Rnum>13000000)bc=d.Bc[1];
    if(Rnum>10000000)
    {
      if((LEDbit)%16>0)
      {
        LEDevent=true; 
      };
    }

    else
    {
      if((LEDbit/2048)%2==1)
      {
        LEDevent=true; 
        lastLED=bc;
        if(lastLED!=0 && ((bc-lastLED)%2000001<10) || (bc-lastLED)%2000001>1999997)
        {
          LEDevent=true;
        };  
      };
    };

    if(Rnum>12000000)//LKE 20110211
    {
      LEDevent=false;
      //	  printf("setting TrigBits for Rnum=%d \n",Rnum);
      TrigBits=d.lastdsm[5];	  
      LEDbit=(d.lastdsm[4] & 0x0001);
      if(LEDbit==1)
      {
        LEDevent=true; 
      };
      //TrigBits=(fpdebit<<14)+(dijetbit<<11)+(L2_jp_bit[2]<<10)+(L2_jp_bit[1]<<9)+(L2_jp_bit[0]<<8)
      //+(lrgclbit2<<7)+(lrgclbit1<<6)+(lrgclbit0<<5)+(L1_clbit[2][0]<<4)+(L1_clbit[1][0]<<3)+(L1_clbit[0][0]<<2)
      //+(fmshtbit1<<1)+fmshtbit0;
    };




    if(IdentifyLED==false)LEDevent=false;
    if(AllLed())LEDevent=true;
    if(!LEDevent && OnlyLED )continue;

    //      printf("TrigBitsSelect=%x TrigBits=%x \n",TrigBitsSelect,TrigBits);
    if((TrigBitsSelect!=0) && (TrigBitsSelect&TrigBits)==0)continue;
    if( (TrigBitsVeto&TrigBits)!=0)continue;
    if(!d.decodeQT()){printf("qt Error \n"); continue;};
    TMatrix* madc;
    led=0;
    if(LEDevent)led=1;
    evtnum=i;
    if(clusts)delete clusts;
    clusts=new TObjArray(10,0);
    clusts->SetOwner();
    if(threeV)delete threeV;
    threeV=new TObjArray(10,0);
    threeV->SetOwner();
    int allclcnt=0;
    HitCluster* HiClu=0;
    HitCluster* NextHiClu=0;
    HiCluE=NextHiCluE=0.;
    if(Rnum>13000000)
    { 
      if(d.lastdsm[4]&1>0)LEDevent=true;
      if(abs((Int_t) (bc-lastLED))<2)
      {
        LEDevent=true;
        LEDecho=true;
      }; 
      for(int iadc=0;iadc<4;iadc++)
      {
        nstbADC=iadc+1;

        madc=(TMatrix*) d.pQt->qtMat->At(iadc);
        int nrow=madc->GetNrows();
        int ncol=madc->GetNcols();
        int k=iadc;
        p_adc[k]=new TMatrix(d.dMatrixQt(k));
        ADCDet+=p_adc[k]->Sum();
        p_Emat[k]=new TMatrix(d.Em(p_adc[k],2,k+1));
        Edet+=p_Emat[k]->Sum();;
      }
      if(Edet>1000)
      {
        LEDevent=true;
      };
      if(LEDevent)led=1;
      if(LEDecho)led=2;
    };
    if(LEDevent)lastLED=bc;

    Edet=0;
    for(int iadc=0;iadc<4;iadc++)
    {
      nstbADC=iadc+1;

      madc=(TMatrix*) d.pQt->qtMat->At(iadc);
      int nrow=madc->GetNrows();
      int ncol=madc->GetNcols();
      int k=iadc;
      d.dMatrixQt(k);
      p_adc[k]=new TMatrix(d.dMatrixQt(k));
      p_Emat[k]=new TMatrix(d.Em(p_adc[k],2,k+1));
      Edet+=p_Emat[k]->Sum();;

      if(UseRunDepCor)RCorrect(p_adc[k],p_Emat[k],k);
      if(recon[k])delete recon[k];
      recon[k]=new Yiqun(p_Emat[k],d.pGeom,d.Rgain,d.Rgaincorr,2,k+1);
      //	  printf("recon: NSTB=%d Nclust=%d \n",k+1,recon[k]->Clu2->GetEntries());
      TObjArray* clu2=recon[k]->Clu2; 
      TIter nex0(clu2);
      HitCluster* hcl=0;

      while(hcl=(HitCluster*) nex0())
      {

        hcl->CalClusterMoment(.001);
        hcl->EDepUpdate(recon[k]->pwe);
        HitCluster* nhcl=new HitCluster();
        *nhcl=*hcl;
        nhcl->IEW=2;
        nhcl->NSTB=k+1;
        if(!HiClu)
        {
          HiClu=nhcl;
        } 
        else
        {
          if(nhcl->energy>HiClu->energy)
          {
            NextHiClu=HiClu;
            HiClu=nhcl;
          }
          else if(NextHiClu)
          { 
            if(nhcl->energy > NextHiClu->energy)
            { NextHiClu=nhcl;};
          }
          else 	{ NextHiClu=nhcl;};
        };

        if(HiClu)HiCluE=HiClu->energy;
        if(NextHiClu)NextHiCluE=NextHiClu->energy;
        //   printf("HiCluE=%f NextHiCluE=%f \n",HiCluE,NextHiCluE);
        clusts->AddAtAndExpand(nhcl,allclcnt);
        Float_t* wid=d.pGeom->FpdTowWid(2,k+1);
        TVector3* vec=new TVector3(hcl->x0*wid[0],hcl->y0*wid[1],0);
        *vec=d.pGeom->GlobalXYZ(2,k+1,*vec);
        //	      printf("x0=%f,y0=%f lab: x0=%f y0=%f z0=%f pseudo=%f\n",hcl->x0,hcl->y0,vec->X(),vec->Y(),vec->Z(),vec->PseudoRapidity());

        //	      printf("count=%d NSTB=%d\n",allclcnt,k+1);
        threeV->AddAtAndExpand(vec,allclcnt++);
        //((HitCluster*) (*clusts)[allclcnt-1])->Print();
      };
    };

    minang=1000.;
    minE=-1.;
    minM=-1.;
    for(int jj=0;jj<threeV->GetEntries();jj++)
    {

      TVector3* vec= (TVector3*) (*threeV)[jj];
      PseudoR=vec->PseudoRapidity();
      TIter nex1(threeV);
      HitCluster* hcl=(HitCluster*) (*clusts)[jj];
      TVector3 mom=vec->Unit()*hcl->energy;
      TLorentzVector lv1(mom,mom.Mag());

      for(int jk=0;jk<threeV->GetEntries();jk++)
      {
        TVector3* vec2=(TVector3*) (*threeV)[jk];
        HitCluster* hcl2=(HitCluster*) (*clusts)[jk];
        float ang=vec->Angle(*vec2);
        if(vec!=vec2 && ang<minang && hcl2->energy>5.)
        {
          TVector3 mom2=vec2->Unit()*hcl2->energy;
          minang=ang;
          TLorentzVector lv1(mom,mom.Mag());
          TLorentzVector lv2(mom2,mom2.Mag());
          minE=(lv1+lv2).E();
          minM=(lv1+lv2).Mag();
        }
      };
      nstbADC=hcl->NSTB;
      TObjArray* tow=hcl->tow;
      ClusterCount=tow->GetEntries();
      ClusterE=hcl->energy;
      TIter next(tow);
      TowerFPD* to=0;
      int cntr=0;
      dToHiTower=0.;
      TowerFPD* HiTower=0;
      //	  printf("Towercnt=%d\n",tow->GetEntries());
      while(to=(TowerFPD*) next())
      {
        if(cntr==0)HiTower=to;;
        ClusterOrder=cntr++;
        rowADC=to->row-1;
        colADC=to->col-1;
        //	      cout<<"dTowHi..\n";
        dToHiTower=sqrt(pow(1.*(rowADC+1-HiTower->row),2)+
            pow(1.*(colADC+1-HiTower->col),2));
        //	      cout<<"dToHiTower="<<dToHiTower<<"\n";
        ADC=to->adc_over_ped;
        CellE=to->energy;
        Tr_adc.Fill();
      };
    };

    pcnt++;

    if(pcnt>100){printf("cnt=%d \n",i);pcnt=1;}; //print every 5000 events
  };
  printf("Done with %d  event loop \n",nentries);
  //  ndatvssum.Write();
  printf("ndatvsssu.Write() done \n");
  Tr_File.cd();
  Tr_adc.Write();
  adcout.cd();
  printf("Tr_adc.Write() done \n");
  /*
     d.pQt->qtHist[0]->Write();
     d.pQt->qtHist[1]->Write();
     d.pQt->qtHist[2]->Write();
     d.pQt->qtHist[3]->Write();
     */
  TH2F* qtHist0[4];
  TH2F* qtHist1[4];
  qtHist0[0]=new TH2F("qtHistLN","qtHistLN",578,0,578,4096,0,4096);
  qtHist0[1]=new TH2F("qtHistLS","qtHistLS",578,0,578,4096,0,4096);
  qtHist0[2]=new TH2F("qtHistSN","qtHistSN",288,0,288,4096,0,4096);
  qtHist0[3]=new TH2F("qtHistSS","qtHistSS",288,0,288,4096,0,4096);
  qtHist1[0]=new TH2F("qtHistLedLN","qtHistLedLN",578,0,578,4096,0,4096);
  qtHist1[1]=new TH2F("qtHistLedLS","qtHistLedLS",578,0,578,4096,0,4096);
  qtHist1[2]=new TH2F("qtHistLedSN","qtHistLedSN",288,0,288,4096,0,4096);
  qtHist1[3]=new TH2F("qtHistLedSS","qtHistLedSS",288,0,288,4096,0,4096);
  printf("8 qthists defined \n");
  Int_t NCOL[4]={17,17,12,12};
  for(int j=0;j<Tr_adc.GetEntries();j++)
  {
    Int_t Nb=Tr_adc.GetEntry(j);
    if( (j%5000 == 0 ))
      printf("read Tr_adc j=%d  nbytes = %d \n",j,Nb);

    Int_t nstb=nstbADC;
    if(nstb>0 && nstb<5)
    {
      Int_t chan0=(rowADC)*NCOL[nstb-1]+colADC;
      if(led==0)
      {
        qtHist0[nstb-1]->Fill(chan0,ADC);
      }
      else
      {
        qtHist1[nstb-1]->Fill(chan0,ADC);
      };
    };
  };
  for(int nh=0;nh<4;nh++)
  {
    qtHist0[nh]->Write();
    qtHist1[nh]->Write();
  };
  printf("d.pQt->qtHist[3]->.Write() done \n");
  printf("ready to return \n");
  //end event loop
  return 0;
};

Int_t AnalTools::Select_scriptQT1(Int_t set,FilesSet* P_Files,char* hntName)
{

  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=false;
  AnalDet[0][1]=false;
  AnalDet[1][0]=true;
  AnalDet[1][1]=true;
  AnalDet[1][2]=true;
  AnalDet[1][3]=true;
  if(P_Files==0)
  {
    if(set==80)
    {
      //transverse run 8
      p_files=new FilesSet("/star/u/heppel/BatchDir/rootfms",
          "fpdped.txt",
          "fpdgain.txt", 
          "fpdcorr.txt",
          "fill.txt",
          "Fake",
          "fpd08/spinpat",
          "geom.txt",
          "qtmap.txt",
          "qtmap2pp.txt");


      p_files->Print();
      hnt_name="h111_08";

    }
    else
    {
      printf("no valid set selected for Select_scriptQT1\n");
      return 0;
    };
  }
  else    {
    p_files=P_Files;
    hnt_name=hntName;
  };



  char hntp[20];
  strcpy(hntp,hnt_name);
  printf("dataSet d(./run*.root,p_files,%s)\n",hntp);
  dataSet d("./run*.root",p_files,hntp);
  currentdset=&d;
  std::cout<<"ok\n";
  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  Int_t Bunchid7bit=d.bunchid7bit;

  // FpdMap is a class that figures out maps into data blocks
  //  FpdMap* pmap=new FpdMap(Rnum);

  TFile* p_OFile=new TFile("OFile.root","recreate");

  //
  Int_t tpes[320];
  Int_t spin;
  Int_t nphotons;
  Float_t px[10];
  Float_t py[10];
  Float_t pz[10];
  Float_t pE[10];
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t strig;
  Char_t BBcSums[5];
  Float_t BBcVertex[7];
  Int_t EventN;
  Int_t ievt;

  QTBBCInfo QBI;
  QTRPInfo QRPI;
  //Define a Tree
  TTree* p_out = new TTree("p_out","Output Tree");
  TH1F Emat[5];
  Emat[0]=TH1F("Emat0","Emat0",2000,0,2000);
  Emat[1]=TH1F("Emat1","Emat1",2000,0,2000);
  Emat[2]=TH1F("Emat2","Emat2",2000,0,2000);
  Emat[3]=TH1F("Emat3","Emat3",2000,0,2000);
  Emat[4]=TH1F("EmatAll","EmatAll",2000,0,2000);
  TH1F NHi[5];
  NHi[0]=TH1F("NHi0","NHi0",300,0,300);
  NHi[1]=TH1F("NHi1","NHi1",300,0,300);
  NHi[2]=TH1F("NHi2","NHi2",300,0,300);
  NHi[3]=TH1F("NHi3","NHi3",300,0,300);
  NHi[4]=TH1F("NHiAll","NHiAll",300,0,300);

  p_out->Branch("spin",&(spin),"spin/I"); // spin of protons in bXing (0=--, 1=-+, 2=+1, 3=++ (notation (B)(Y)), 40=otherwise)
  p_out->Branch("nphotons",&(nphotons),"nphotons/I"); //  number of photons found by yiqun
  p_out->Branch("br_nwrds",&(nwrds),"nwrds/I"); // number of "words"
                                                // incremented by 1 for each NSTB which passes "esum-cut"
                                                // incremented also by 1 for each photon found 
                                                // finally, multiplied by 4 after all increments above have occurred for all 4 NSTBs
  p_out->Branch("br_types",(tpes),"tpes[nwrds]/I"); // "types array" -- the type of int[4]-length "word" (either new NSTB or new photon... see below)
                                                    // the way the code loops below after the "event-cut" is that it loops through the 4 NSTBs, and then
                                                    // through the photons found after the reconstruction
                                                    // tpes is filled with "5+NSTB" at entry nwrds*4, then nwrds is incremented by 1
                                                    // tpes is then filled with "305+NSTB" at entry nwrds+4 for each new photon, then nwrds is incremented by 1
                                                    // all other entries are 0
                                                    // example -- say 2 photons in nstb0 and 1 in nstb1, zero in the others --> nwrds = [ (4 nstbs) + (3 total photons) ] * 4 = 7*4 = 28
                                                    // index   entry
                                                    // 0       5 + 0 = 5
                                                    // 1       0
                                                    // 2       0
                                                    // 3       0
                                                    // 4       305 + 0 = 305
                                                    // 5       0
                                                    // 6       0
                                                    // 7       0
                                                    // 8       305 + 0 = 305 
                                                    // 9       0
                                                    // 10      0
                                                    // 11      0
                                                    // 12      5 + 1 = 6
                                                    // 13      0
                                                    // 14      0
                                                    // 15      0
                                                    // 16      305 + 1 = 306
                                                    // 17      0
                                                    // 18      0
                                                    // 19      0
                                                    // 20      5 + 2 = 7
                                                    // 21      0
                                                    // 22      0
                                                    // 23      0
                                                    // 24      5 + 3 = 8
                                                    // 25      0
                                                    // 26      0
                                                    // 27      0
  p_out->Branch("br_pxyzt",(pxyzt),"pxyzt[nwrds]/F"); // momentum 4-vectors for each word
                                                      // if type<300, then the 4-vector is total 4-momentum of all of the hits in the NSTB  (see AnalTools::FourMom)
                                                      // if type>300, then the 4-vector is the 4-momentum  of the photon
                                                      // the reason the types are filled every 4 entries above is because we store 4-vector info in this array: 
                                                      // example -- say 2 photons in nstb0 and 1 in nstb1, zero in the others --> nwrds = [ (4 nstbs) + (3 total photons) ] * 4 = 7*4 = 28
                                                      // index   entry
                                                      // 0       p_x of NSTB 0
                                                      // 1       p_y of NSTB 0
                                                      // 2       p_z of NSTB 0
                                                      // 3       Esum of NSTB 0
                                                      // 4       p_x of 1st photon in NSTB 0
                                                      // 5       p_y of 1st photon in NSTB 0
                                                      // 6       p_z of 1st photon in NSTB 0
                                                      // 7       E of 1st photon in NSTB 0
                                                      // 8       p_x of 2nd photon in NSTB 0
                                                      // 9       p_y of 2nd photon in NSTB 0
                                                      // 10      p_z of 2nd photon in NSTB 0
                                                      // 11      E of 2nd photon in NSTB 0
                                                      // 12      p_x of NSTB 1
                                                      // 13      p_y of NSTB 1
                                                      // 14      p_z of NSTB 1
                                                      // 15      E of NSTB 1
                                                      // 16      p_x of 1st photon in NSTB 1
                                                      // 17      p_y of 1st photon in NSTB 1
                                                      // 18      p_z of 1st photon in NSTB 1
                                                      // 19      E of 1st photon in NSTB 1
                                                      // 20      p_x of NSTB 2
                                                      // 21      p_y of NSTB 2
                                                      // 22      p_z of NSTB 2
                                                      // 23      E of NSTB 2
                                                      // 24      p_x of NSTB 3
                                                      // 25      p_y of NSTB 3
                                                      // 26      p_z of NSTB 3
                                                      // 27      E of NSTB 3
  p_out->Branch("br_Rnum",&(Rnum),"Rnum/I"); // run number 
  p_out->Branch("br_Bunchid7bit",&(Bunchid7bit),"Bunchid7bit/I"); // bXing number
  p_out->Branch("br_BBcSums",BBcSums,"BBcSums[5]/b"); // not needed for TwoTr creation; just passed along
  p_out->Branch("br_BBcVertex",BBcVertex,"BBcVertex[7]/F"); // not needed for TwoTr creation; just passed along


  p_out->Branch("br_EventN",&(EventN),"EventN/I"); // deprecated ?; comes from dataSet::event, which is from branch "event" from much older version h111 trees
  p_out->Branch("br_ievt",&(ievt),"ievt/I"); // event number (assuming each segment has 10,000 events
  p_out->Branch("br_nSavedHits",&(nSavedHits),"nSavedHits/I"); // number of saved hits, see AnalTools::storeCluster, which calls AnalTools::SaveClHit, which increments this number
                                                               // for call to AnalTools::storeCluser, search for "STORE-CLUSTER" --> minimum ADC = 1; min cluster energy = 0.5 
  p_out->Branch("br_SavedHits",(SavedHits),"SavedHits[nSavedHits]/I"); // encoding of a single hit: encodes east/west, nstb, row, column, and ADC
                                                                       // see AnalTools::storeCluster which calls AnalTools::SaveClHit, which does this encoding for each hit tower 
                                                                       //
                                                                       // the 4-byte encoding is layed out in the following way:
                                                                       //  - E = EW
                                                                       //  - N = NSTB
                                                                       //  - R = Row
                                                                       //  - C = Col
                                                                       //  - A = ADC
                                                                       //
                                                                       //   |  E N N N  |  R R R R  |  R R C C  |  C C C C  |  A A A A  |  A A A A  |  A A A A  |  A A A A  | 
                                                                       //  
                                                                       // decode it by using the following (decimal places put in hex numbers for readability!)
                                                                       // --> EW =   ( SavedHits & 0x.8000.0000 ) / 0x.8000.0000
                                                                       // --> NSTB = ( SavedHits & 0x.7000.0000 ) / 0x.1000.0000
                                                                       // --> Row  = ( SavedHits & 0x.0F90.0000 ) / 0x.0040.0000
                                                                       // --> Col =  ( SavedHits & 0x.003F.0000 ) / 0x.0001.0000
                                                                       // --> ADC =    SavedHits & 0x.0000.FFFF 
                                                                       //
                                                                       // ARE THESE HITS ORDERED?? --> orderd by highest-energy cluster
  p_out->Branch("br_TrigBits",&TrigBits,"TrigBits/I"); // FMS layer 2 (DSMs) output bits (i.e., input bits to TCU, known as lasatdsm[5] in the h111 tree)
                                                       // to filter, use (TrigBits >> (L2_output_bit)) & 0x1
                                                       // 2015 pp output bits:
                                                       // 0 - smbs3
                                                       // 1 - smbs2
                                                       // 2 - smbs1
                                                       // 3 - lgbs3
                                                       // 4 - lgbs2
                                                       // 5 - lgbs1
                                                       // 6 - dibs
                                                       // 7 - jp2
                                                       // 8 - jp1
                                                       // 9 - jp0
                                                       // 10 - dijet
                                                       // 11 - unused
                                                       // 12 - unused
                                                       // 13 - unused
                                                       // 14 - unused
                                                       // 15 - unused

  p_out->Branch("br_nCluster",&(nCluster),"nCluster/I"); // number of clusters; incremented every time AnalTools::storeCluster is called
  p_out->Branch("br_nPhotonClu",&(nPhotonClu),"nPhotonClu/I"); // total number of photons within clusters (equivalent to StMuFmsPoints?)
                                                               // n.b. nPhotonClu-nCluster >= 0 for all events
  p_out->Branch("br_SCIndex",(SavedCluHitIndex),"SCIndex[nCluster]/I"); // number of saved hits (nSavedHits) for current cluster, but... 
                                                                        // nSavedHits is incremented as we loop through each cluster, so basically:
                                                                        // to find out how many hits are in cluster number 3, we do: SCIndex[3]-SCIndex[2]
                                                                        // to find out how many hits are in the first cluster, we do: SCIndex[1]-SCindex[0] (n.b. SCindex[0]==0 for all events)
                                                                        // to find out how many hits are in the last cluster, we do: nSavedHits-SCindex[nCluster-1]
  p_out->Branch("br_SPCIndex",(SavedPhotonCluIndex),"SPCIndex[nPhotonClu]/I"); // cluster index (starting at 0) of nth photon (ordered by highest-energy clusters)
  p_out->Branch("br_SPCEnergy",(SavedPhotonCluEnergy),"SPCEnergy[nPhotonClu]/F"); // photon energy
  p_out->Branch("br_L2sum",L2sum,"L2sum[2]/i"); // essentially the trigger ID
  p_out->Branch("br_lastdsm",lastdsm,"lastdsm[8]/i"); // full TCU input bitstream
  p_out->Branch("br_Fpde",Fpde,"Fpde[8]/i"); // I think this is just passed from h111; it's not in heppel's OFiles, so this may just be something I added

  //QTB branch

  p_out->Branch("br_QTNE",&(QBI.QTNE),"BQTNE/B");
  p_out->Branch("br_QTNW",&(QBI.QTNW),"BQTNW/B");   
  p_out->Branch("br_QTEBBCInd",(QBI.QTEBBCInd),"BQTEInd[BQTNE]/B");
  p_out->Branch("br_QTWBBCInd",(QBI.QTWBBCInd),"BQTWInd[BQTNW]/B");

  p_out->Branch("br_QTEBBCTAC",(QBI.QTEBBCTAC),"QTEBTAC[BQTNE]/S");
  p_out->Branch("br_QTWBBCTAC",(QBI.QTWBBCTAC),"QTWBTAC[BQTNW]/S");

  p_out->Branch("br_QTEBBCADC",(QBI.QTEBBCADC),"QTEBADC[BQTNE]/S");
  p_out->Branch("br_QTWBBCADC",(QBI.QTWBBCADC),"QTWBADC[BQTNW]/S");

  p_out->Branch("br_QTBVertex",&(QBI.vertex),"QTBVertex/F");

  // QTRP branch

  p_out->Branch("br_RPE_QTN",&(QRPI.NE),"RPE_QTN/B");
  p_out->Branch("br_RPW_QTN",&(QRPI.NW),"RPW_QTN/B");
  p_out->Branch("br_RPE_Idx",(QRPI.RPE_Idx),"RPE_Idx[RPE_QTN]/B");
  p_out->Branch("br_RPW_Idx",(QRPI.RPW_Idx),"RPW_Idx[RPW_QTN]/B");
  p_out->Branch("br_RPE_TAC",(QRPI.RPE_TAC),"RPE_TAC[RPE_QTN]/S");
  p_out->Branch("br_RPW_TAC",(QRPI.RPW_TAC),"RPW_TAC[RPW_QTN]/S");
  p_out->Branch("br_RPE_ADC",(QRPI.RPE_ADC),"RPE_ADC[RPE_QTN]/S");
  p_out->Branch("br_RPW_ADC",(QRPI.RPW_ADC),"RPW_ADC[RPW_QTN]/S");
  p_out->Branch("br_RPvertex",&(QRPI.vertex),"RPvertex/F");


  //
  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  printf("NumberEventsToProcess found = %d ",NumberEventsToProcess);
  if((NumberEventsToProcess>0) && (NumberEventsToProcess<nentries))nentries=NumberEventsToProcess;
  printf(" nentries will be = %d \n",nentries);
  // event loop
  Int_t skpcnt=0;
  UInt_t lastLED=0;

  TMatrix* p_adc[4]={0,0,0,0};
  TMatrix* p_Emat[4]={0,0,0,0};
  Vertex* vtx=new Vertex();
  //vtx->Bbc_threshold=15;

  for(Int_t i=0;i<nentries;i++)
  {
    nSavedHits=0;
    nPhotonClu=0;
    nCluster=0;
    TrigBits=0;

    Int_t nbytes=d.GetEntry(i);
    Bool_t trgXing=true;
    Rnum=d.CurrentRunNumber;


    if(Rnum<10000000)
    {
      if(Rnum<10000000 &&(d.ipre!=0 || d.ipost!=0)){trgXing=false; continue;};
    };
    EventN=d.event;
    ievt=i+(d.CurrentSegNumber-1)*10000;//assume 10k per segmant
    if(!d.decodeQT()){printf("qt Error \n"); continue;};
    pbbc->ADCmin=70;
    pbbc->scale=1.;
    pbbc->VtxRange[0]=1;
    pbbc->VtxRange[1]=16;
    pbbc-> UpdateQTBBCinfo(&QBI,d.pQt);
    prp->scale=1.;
    QRPI.SetZero(); //+++
    prp-> UpdateQTRPinfo(&QRPI,d.pQt);
    Bunchid7bit=d.bunchid7bit;
    pcnt++;
    if(pcnt%100==1){printf("cnt=%d \n",i);pcnt=1;}; //print every 100 events
    Double_t NearUnique;

    spin=40;
    if(abs(d.BlueSpin*d.YellowSpin)==1)
    {
      spin=(d.BlueSpin+1)+(d.YellowSpin+1)/2;
    };
    //{0,1,2,3}={bd_yd,bd_yu,bu_yd,bu,yu};
    if(d.kicked)spin+=29;
    if(d.BlueSpin*d.YellowSpin==0)spin=40;
    //Spin repair
    spin=0;
    int bspn=d.BlueSpin;
    int yspn=d.YellowSpin;
    if(abs(bspn)==1)spin+=(bspn+1);
    if(abs(yspn)==1)spin+=(yspn+1)/2;
    if(bspn*yspn==0)spin+=32;
    if(bspn==0&&yspn==0)spin=64;
    //end spin repair
    Trigger trig(d.Bbcl1);
    for(Int_t ki=0;ki<5;ki++)BBcSums[ki]=trig.BBcSums[ki];
    BBcVertex[0]=vtx->GetBbcVertex(d.Bbc);
    BBcVertex[1]=vtx->maxtace;
    BBcVertex[2]=vtx->maxtacw;
    BBcVertex[3]=vtx->iemax;
    BBcVertex[4]=vtx->iwmax;
    BBcVertex[5]=vtx->qe;
    BBcVertex[6]=vtx->qw;


    Float_t Esum[4]={0.,0.,0.,0.};
    Float_t Esum4=0.;
    Int_t Nhigh[4]={0,0,0,0};
    Int_t Nhigh4=0;
    Float_t ECorRat[4]={0.,0.,0.,0.};
    for(Int_t k=0;k<4;k++)
    {
      if(AnalDet[1][k]){
        p_adc[k]=new TMatrix(d.dMatrixQt(k));
        p_Emat[k]=new TMatrix(d.Em(p_adc[k],2,k+1));
        if(UseRunDepCor&&Mgr)
        {
          ECorRat[k]=p_Emat[k]->Sum();
          Mgr->RunDepBaseOverride=RunDepBaseOverride;
          Mgr->Rcorrect(p_adc[k],p_Emat[k],k,2,&d);
          if(ECorRat[k]!=0)ECorRat[k]=p_Emat[k]->Sum()/ECorRat[k];
        }
        Esum[k]=p_Emat[k]->Sum();
        Emat[k].Fill(Esum[k]);

        Esum4+=Esum[k];
        //	    printf("accum Esum4[%d]=%f\n",k,Esum4); 
        Nhigh[k]= NTowersAbove(p_Emat[k],.25);
        Nhigh4+=Nhigh[k];
        NHi[k].Fill(Nhigh[k]);
      };
    };
    Emat[4].Fill(Esum4);
    NHi[4].Fill(Nhigh4);
    hardtrig=1;
    nwrds=0;
    nphotons=0;

    Bool_t LEDevent=false;
    Bool_t LEDecho=false;
    UInt_t LEDbit=d.lastdsm[7];
    TrigBits=d.lastdsm[7]; // overridden later!
    Int_t bc=d.Bclo+d.Bchi*(Int_t)pow(2,32);
    if(Rnum>13000000)bc=d.Bc[1];

    //printf("%d %d %d \n",bc,(UInt_t)d.ipre,(UInt_t)d.ipost);
    Bool_t FpdTrig=false;
    Bool_t HiTow=false;
    Bool_t Jet=false;

    if(Rnum>10000000)
    {
      if((LEDbit)%16>0)
      {
        LEDevent=true; 
      };
    }

    else
    {
      if((LEDbit/2048)%2==1)
      {
        LEDevent=true; 
        lastLED=bc;
        if(lastLED!=0 && ((bc-lastLED)%2000001<10) || (bc-lastLED)%2000001>1999997)
        {
          LEDevent=true;
        };  
      };
    };

    if(Rnum>12000000)//LKE 20110211
    {
      for(Int_t qq=0; qq<2; qq++) L2sum[qq]=d.L2sum[qq];
      for(Int_t qq=0; qq<8; qq++) lastdsm[qq]=d.lastdsm[qq];
      for(Int_t qq=0; qq<8; qq++) Fpde[qq]=d.Fpde[qq];
      LEDevent=false;
      TrigBits=d.lastdsm[5]; // TCU input 5 corresponds to FMS L2 output
      LEDbit=(d.lastdsm[4] & 0x0001);
      if(LEDbit==1)
      {
        LEDevent=true; 
      };
      //TrigBits=(fpdebit<<14)+(dijetbit<<11)+(L2_jp_bit[2]<<10)+(L2_jp_bit[1]<<9)+(L2_jp_bit[0]<<8)
      //+(lrgclbit2<<7)+(lrgclbit1<<6)+(lrgclbit0<<5)+(L1_clbit[2][0]<<4)+(L1_clbit[1][0]<<3)+(L1_clbit[0][0]<<2)
      //+(fmshtbit1<<1)+fmshtbit0;
    };

    // until we fix the bit

    if((Rnum>13000000 && Esum4>1000) || (d.lastdsm[4]&1) >0)
    {
      LEDevent==true;
      lastLED=bc;
    }
    else if(Rnum>13000000 && abs((Int_t) (lastLED-bc))<2)
    {
      lastLED=bc;
      LEDevent=true;
      LEDecho=true;
    }

    Bool_t Debug=false;
    if(IdentifyLED==false)LEDevent=false;
    if(Debug)
    {
      if(LEDevent)printf("lastLED=%d bc=%d mod_bc=%d Nhigh4=%d \n",lastLED,bc,(bc-lastLED)%2000001,Nhigh4);      
    };

    TLorentzVector vecs[80];
    for(Int_t kj=0;kj<320;kj++)tpes[kj]=0; 
    Bool_t pass=Esum4>40 || (Esum4>30 && (fmod((skpcnt++),10)<OutOf10-.1));
    pass=Esum[0]>ERecMin[0] || Esum[1]>ERecMin[1] || Esum[2]>ERecMin[2] || Esum[3]>ERecMin[3];

    hardtrig=1;
    if((TrigBitsSelect!=0) && (TrigBits &TrigBitsSelect)==0)hardtrig=0;
    if( (TrigBitsVeto&TrigBits)!=0)hardtrig=0;

    //      printf("led=%d Nhigh4=%d Esum4=%f\n",LEDevent,Nhigh4,Esum4);
    //      if((!LEDevent) &&hardtrig==1 && Nhigh4<70 && Esum4<400  && pass )

    // event-cut
    if((!LEDevent) &&hardtrig==1 && Nhigh4<300 && Esum4<3000  && pass )
    {
      // nstb loop
      for(Int_t k=0;k<4;k++)
      {
        // esum-cut
        //	      if(Esum[k]>1. && Esum[k]<350.)
        if(Esum[k]>1. && Esum[k]<700.)
        {
          TLorentzVector lv(FourMom(p_Emat[k],0,100,0,100,d.pGeom,2,k+1));
          tpes[nwrds*4]=k+5; // new NSTB
          vecs[nwrds]=lv; // 4 vector based on energy in NSTB
          nwrds++; // itereated every NSTB with energy sum >1 && <700 for every event which passes event cut above
          if(i<5)
          {
            p_Emat[k]->Print();
          };
          FitTower::SetSubclu2(false);
          if(i==204)FitTower::SetYPrintLevel(2);      
         
          //////////////////////////////////////////////////////////
          Yiqun recon(p_Emat[k],d.pGeom,d.Rgain,d.Rgaincorr,2,k+1);
          //////////////////////////////////////////////////////////
          
          FitTower::SetYPrintLevel(0);      
          if(i==204)
          {

            p_Emat[k]->Print();
            for(int jpho=0;jpho<recon.NPh;jpho++)
            {
              printf("Event = %d photon for k:jpho=%d:%d",i,k,jpho);
              printf(" Energy=%f \n",recon.mom(jpho).E());
              HitCluster* aClust=(HitCluster*) recon.photons[jpho].pHitCluster;
              if(aClust)aClust->Print();
            }
            for(int jcc=0;jcc<recon.NRealClusts;jcc++)
            {
              printf("Event = %d Cluster for k=%d\n",i,k);
              recon.clust[jcc].Print();
            } 
          }
          Int_t nz=recon.NPh;
          if(nz>10)nz=10;

          for(Int_t n=0;n<nz;n++)
          {
            vecs[nwrds]=recon.mom(n); // photon 4 momentum
            float myphi=recon.mom(n).Phi();
            tpes[nwrds*4]=305+k; // new photon
            nwrds++; // now iterates for every photon
            nphotons++;
          };
          //		  if(Esum[k]>10 && k>1)CheckCluster(&recon);
          //		  if(Esum[k]>10 && k>-1)storeCluster(&recon,1,1.5);
          //		  if(Esum[k]>2 && k>-1)storeCluster(&recon,1,.5);
          //		  if(Esum[k]>1 && k>-1)storerec(&recon,1);
          if(Esum[k]>1 && k>-1)storeCluster(&recon,1,.5); // STORE-CLUSTER -----------------
        };
      }; // eo NSTB loop
    }; // eo event cut
    //      if(nSavedHits>0)printf("number stored =%d \n",nSavedHits);
    Int_t wcnt=0;
    for(Int_t k=0;k<nwrds;k++)
    {
      pxyzt[wcnt]=vecs[k].Px();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Py();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Pz();
      wcnt++;
      pxyzt[wcnt]=vecs[k].E();
      wcnt++;
    };
    nwrds=nwrds*4; // multiply by 4

    if(nwrds>0)p_out->Fill();
    /*      printf(" i=%d EventN=%d Correction Ratios %f %f %f %f \n adc= %f %f %f %f Energy= %f %f %f %f \n",i,EventN,
            ECorRat[0],ECorRat[1],ECorRat[2],ECorRat[3],
            p_adc[0]->Sum(),p_adc[1]->Sum(),p_adc[2]->Sum(),p_adc[3]->Sum(),
            p_Emat[0]->Sum(),p_Emat[1]->Sum(),p_Emat[2]->Sum(),p_Emat[3]->Sum());
            */
    for(Int_t k=0;k<4;k++)
    {
      if(p_adc[k]){delete p_adc[k];p_adc[k]=0;};
      if(p_Emat[k]){delete p_Emat[k];p_Emat[k]=0;};
    };
  };

  //end event loop
  p_out->Print();
  p_out->Write();
  /*
     for(int jj=0;jj<5;jj++)
     {
     Emat[jj].Write();
     NHi[jj].Write();
     };
     */
  printf("about to return\n");
  return 0;
};


// for testing RP branches
Int_t AnalTools::Select_script_RPtest(Int_t set,FilesSet* P_Files,char* hntName)
{

  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=false;
  AnalDet[0][1]=false;
  AnalDet[1][0]=true;
  AnalDet[1][1]=true;
  AnalDet[1][2]=true;
  AnalDet[1][3]=true;
  if(P_Files==0)
  {
    if(set==80)
    {
      //transverse run 8
      p_files=new FilesSet("/star/u/heppel/BatchDir/rootfms",
          "fpdped.txt",
          "fpdgain.txt", 
          "fpdcorr.txt",
          "fill.txt",
          "Fake",
          "fpd08/spinpat",
          "geom.txt",
          "qtmap.txt",
          "qtmap2pp.txt");


      p_files->Print();
      hnt_name="h111_08";

    }
    else
    {
      printf("no valid set selected for Select_scriptQT1\n");
      return 0;
    };
  }
  else    {
    p_files=P_Files;
    hnt_name=hntName;
  };



  char hntp[20];
  strcpy(hntp,hnt_name);
  printf("dataSet d(./run*.root,p_files,%s)\n",hntp);
  dataSet d("./run*.root",p_files,hntp);
  currentdset=&d;
  std::cout<<"ok\n";
  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  Int_t Bunchid7bit=d.bunchid7bit;

  // FpdMap is a class that figures out maps into data blocks
  //  FpdMap* pmap=new FpdMap(Rnum);

  TFile* p_OFile=new TFile("OFile.root","recreate");

  //
  Int_t tpes[320];
  Int_t spin;
  Int_t nphotons;
  Float_t px[10];
  Float_t py[10];
  Float_t pz[10];
  Float_t pE[10];
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t strig;
  Char_t BBcSums[5];
  Float_t BBcVertex[7];
  Int_t EventN;
  Int_t ievt;

  QTBBCInfo QBI;
  QTRPInfo QRPI;
  //Define a Tree
  TTree* p_out = new TTree("p_out","Output Tree");
  TH1F Emat[5];
  Emat[0]=TH1F("Emat0","Emat0",2000,0,2000);
  Emat[1]=TH1F("Emat1","Emat1",2000,0,2000);
  Emat[2]=TH1F("Emat2","Emat2",2000,0,2000);
  Emat[3]=TH1F("Emat3","Emat3",2000,0,2000);
  Emat[4]=TH1F("EmatAll","EmatAll",2000,0,2000);
  TH1F NHi[5];
  NHi[0]=TH1F("NHi0","NHi0",300,0,300);
  NHi[1]=TH1F("NHi1","NHi1",300,0,300);
  NHi[2]=TH1F("NHi2","NHi2",300,0,300);
  NHi[3]=TH1F("NHi3","NHi3",300,0,300);
  NHi[4]=TH1F("NHiAll","NHiAll",300,0,300);

  //p_out->Branch("spin",&(spin),"spin/I");
  //p_out->Branch("nphotons",&(nphotons),"nphotons/I");
  //p_out->Branch("br_nwrds",&(nwrds),"nwrds/I");
  //p_out->Branch("br_types",(tpes),"tpes[nwrds]/I");
  //p_out->Branch("br_pxyzt",(pxyzt),"pxyzt[nwrds]/F");
  p_out->Branch("br_Rnum",&(Rnum),"Rnum/I");
  p_out->Branch("br_Bunchid7bit",&(Bunchid7bit),"Bunchid7bit/I");
  //p_out->Branch("br_BBcSums",BBcSums,"BBcSums[5]/b");
  //p_out->Branch("br_BBcVertex",BBcVertex,"BBcVertex[7]/F");


  p_out->Branch("br_EventN",&(EventN),"EventN/I");
  p_out->Branch("br_ievt",&(ievt),"ievt/I");
  //p_out->Branch("br_nSavedHits",&(nSavedHits),"nSavedHits/I");
  //p_out->Branch("br_SavedHits",(SavedHits),"SavedHits[nSavedHits]/I");
  //p_out->Branch("br_TrigBits",&TrigBits,"TrigBits/I");
  p_out->Branch("br_L2sum",L2sum,"L2sum[2]/i");
  p_out->Branch("br_lastdsm",lastdsm,"lastdsm[8]/i");
  p_out->Branch("br_Fpde",Fpde,"Fpde[8]/i");

  /*
     p_out->Branch("br_nCluster",&(nCluster),"nCluster/I");
     p_out->Branch("br_nPhotonClu",&(nPhotonClu),"nPhotonClu/I");
     p_out->Branch("br_SCIndex",(SavedCluHitIndex),"SCIndex[nCluster]/I");
     p_out->Branch("br_SPCIndex",(SavedPhotonCluIndex),"SPCIndex[nPhotonClu]/I");
     p_out->Branch("br_SPCEnergy",(SavedPhotonCluEnergy),"SPCEnergy[nPhotonClu]/F");

  //QTB branch

  p_out->Branch("br_QTNE",&(QBI.QTNE),"BQTNE/B");
  p_out->Branch("br_QTNW",&(QBI.QTNW),"BQTNW/B");   
  p_out->Branch("br_QTEBBCInd",(QBI.QTEBBCInd),"BQTEInd[BQTNE]/B");
  p_out->Branch("br_QTWBBCInd",(QBI.QTWBBCInd),"BQTWInd[BQTNW]/B");

  p_out->Branch("br_QTEBBCTAC",(QBI.QTEBBCTAC),"QTEBTAC[BQTNE]/S");
  p_out->Branch("br_QTWBBCTAC",(QBI.QTWBBCTAC),"QTWBTAC[BQTNW]/S");

  p_out->Branch("br_QTEBBCADC",(QBI.QTEBBCADC),"QTEBADC[BQTNE]/S");
  p_out->Branch("br_QTWBBCADC",(QBI.QTWBBCADC),"QTWBADC[BQTNW]/S");

  p_out->Branch("br_QTBVertex",&(QBI.vertex),"QTBVertex/F");
  */

  // QTRP branch

  p_out->Branch("br_RPE_QTN",&(QRPI.NE),"RPE_QTN/B");
  p_out->Branch("br_RPW_QTN",&(QRPI.NW),"RPW_QTN/B");
  p_out->Branch("br_RPE_Idx",(QRPI.RPE_Idx),"RPE_Idx[RPE_QTN]/B");
  p_out->Branch("br_RPW_Idx",(QRPI.RPW_Idx),"RPW_Idx[RPW_QTN]/B");
  p_out->Branch("br_RPE_TAC",(QRPI.RPE_TAC),"RPE_TAC[RPE_QTN]/S");
  p_out->Branch("br_RPW_TAC",(QRPI.RPW_TAC),"RPW_TAC[RPW_QTN]/S");
  p_out->Branch("br_RPE_ADC",(QRPI.RPE_ADC),"RPE_ADC[RPE_QTN]/S");
  p_out->Branch("br_RPW_ADC",(QRPI.RPW_ADC),"RPW_ADC[RPW_QTN]/S");
  p_out->Branch("br_RPvertex",&(QRPI.vertex),"RPvertex/F");


  //
  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  printf("NumberEventsToProcess found = %d ",NumberEventsToProcess);
  if((NumberEventsToProcess>0) && (NumberEventsToProcess<nentries))nentries=NumberEventsToProcess;
  printf(" nentries will be = %d \n",nentries);
  // event loop
  Int_t skpcnt=0;
  UInt_t lastLED=0;

  TMatrix* p_adc[4]={0,0,0,0};
  TMatrix* p_Emat[4]={0,0,0,0};
  Vertex* vtx=new Vertex();
  //vtx->Bbc_threshold=15;

  for(Int_t i=0;i<nentries;i++)
  {
    nSavedHits=0;
    nPhotonClu=0;
    nCluster=0;
    TrigBits=0;

    Int_t nbytes=d.GetEntry(i);
    Bool_t trgXing=true;
    Rnum=d.CurrentRunNumber;


    if(Rnum<10000000)
    {
      if(Rnum<10000000 &&(d.ipre!=0 || d.ipost!=0)){trgXing=false; continue;};
    };
    EventN=d.event;
    ievt=i+(d.CurrentSegNumber-1)*10000;//assume 10k per segmane
    if(!d.decodeQT()){printf("qt Error \n"); continue;};
    /*
       pbbc->ADCmin=70;
       pbbc->scale=1.;
       pbbc->VtxRange[0]=1;
       pbbc->VtxRange[1]=16;
       pbbc-> UpdateQTBBCinfo(&QBI,d.pQt);
       */
    prp->scale=1.;
    QRPI.SetZero(); //+++
    prp-> UpdateQTRPinfo(&QRPI,d.pQt);
    Bunchid7bit=d.bunchid7bit;
    pcnt++;
    if(pcnt%100==1){printf("cnt=%d \n",i);pcnt=1;}; //print every 100 events


    Double_t NearUnique;

    spin=40;
    if(abs(d.BlueSpin*d.YellowSpin)==1)
    {
      spin=(d.BlueSpin+1)+(d.YellowSpin+1)/2;
    };

    //{0,1,2,3}={bd_yd,bd_yu,bu_yd,bu,yu};
    if(d.kicked)spin+=29;
    if(d.BlueSpin*d.YellowSpin==0)spin=40;
    //Spin repair
    spin=0;
    int bspn=d.BlueSpin;
    int yspn=d.YellowSpin;
    if(abs(bspn)==1)spin+=(bspn+1);
    if(abs(yspn)==1)spin+=(yspn+1)/2;
    if(bspn*yspn==0)spin+=32;
    if(bspn==0&&yspn==0)spin=64;
    //end spin repair

    /*

       Trigger trig(d.Bbcl1);
       for(Int_t ki=0;ki<5;ki++)BBcSums[ki]=trig.BBcSums[ki];
       BBcVertex[0]=vtx->GetBbcVertex(d.Bbc);
       BBcVertex[1]=vtx->maxtace;
       BBcVertex[2]=vtx->maxtacw;
       BBcVertex[3]=vtx->iemax;
       BBcVertex[4]=vtx->iwmax;
       BBcVertex[5]=vtx->qe;
       BBcVertex[6]=vtx->qw;
       Float_t Esum[4]={0.,0.,0.,0.};
       Float_t Esum4=0.;
       Int_t Nhigh[4]={0,0,0,0};
       Int_t Nhigh4=0;
       Float_t ECorRat[4]={0.,0.,0.,0.};
       for(Int_t k=0;k<4;k++)
       {
       if(AnalDet[1][k]){
       p_adc[k]=new TMatrix(d.dMatrixQt(k));
       p_Emat[k]=new TMatrix(d.Em(p_adc[k],2,k+1));
       if(UseRunDepCor&&Mgr)
       {
       ECorRat[k]=p_Emat[k]->Sum();
       Mgr->RunDepBaseOverride=RunDepBaseOverride;
       Mgr->Rcorrect(p_adc[k],p_Emat[k],k,2,&d);
       if(ECorRat[k]!=0)ECorRat[k]=p_Emat[k]->Sum()/ECorRat[k];
       }
       Esum[k]=p_Emat[k]->Sum();
       Emat[k].Fill(Esum[k]);

       Esum4+=Esum[k];
    //	    printf("accum Esum4[%d]=%f\n",k,Esum4); 
    Nhigh[k]= NTowersAbove(p_Emat[k],.25);
    Nhigh4+=Nhigh[k];
    NHi[k].Fill(Nhigh[k]);
    };
    };
    Emat[4].Fill(Esum4);
    NHi[4].Fill(Nhigh4);
    hardtrig=1;
    nwrds=0;
    nphotons=0;

    Bool_t LEDevent=false;
    Bool_t LEDecho=false;
    UInt_t LEDbit=d.lastdsm[7];
    TrigBits=d.lastdsm[7];
    Int_t bc=d.Bclo+d.Bchi*(Int_t)pow(2,32);
    if(Rnum>13000000)bc=d.Bc[1];

    //printf("%d %d %d \n",bc,(UInt_t)d.ipre,(UInt_t)d.ipost);
    Bool_t FpdTrig=false;
    Bool_t HiTow=false;
    Bool_t Jet=false;

    if(Rnum>10000000)
    {
    if((LEDbit)%16>0)
    {
    LEDevent=true; 
    };
    }

    else
    {
    if((LEDbit/2048)%2==1)
    {
    LEDevent=true; 
    lastLED=bc;
    if(lastLED!=0 && ((bc-lastLED)%2000001<10) || (bc-lastLED)%2000001>1999997)
    {
      LEDevent=true;
    };  
  };
  };
  */

    if(Rnum>12000000)//LKE 20110211
    {
      for(Int_t qq=0; qq<2; qq++) L2sum[qq]=d.L2sum[qq];
      for(Int_t qq=0; qq<8; qq++) lastdsm[qq]=d.lastdsm[qq];
      for(Int_t qq=0; qq<8; qq++) Fpde[qq]=d.Fpde[qq];
      //LEDevent=false;
      TrigBits=d.lastdsm[5];	  
      //LEDbit=(d.lastdsm[4] & 0x0001);
      /*
         if(LEDbit==1)
         {
         LEDevent=true; 
         };
         */
      //TrigBits=(fpdebit<<14)+(dijetbit<<11)+(L2_jp_bit[2]<<10)+(L2_jp_bit[1]<<9)+(L2_jp_bit[0]<<8)
      //+(lrgclbit2<<7)+(lrgclbit1<<6)+(lrgclbit0<<5)+(L1_clbit[2][0]<<4)+(L1_clbit[1][0]<<3)+(L1_clbit[0][0]<<2)
      //+(fmshtbit1<<1)+fmshtbit0;
    };

  // until we fix the bit

  /*
     if((Rnum>13000000 && Esum4>1000) || (d.lastdsm[4]&1) >0)
     {
     LEDevent==true;
     lastLED=bc;
     }
     else if(Rnum>13000000 && abs((Int_t) (lastLED-bc))<2)
     {
     lastLED=bc;
     LEDevent=true;
     LEDecho=true;
     }

     Bool_t Debug=false;
     if(IdentifyLED==false)LEDevent=false;
     if(Debug)
     {
     if(LEDevent)printf("lastLED=%d bc=%d mod_bc=%d Nhigh4=%d \n",lastLED,bc,(bc-lastLED)%2000001,Nhigh4);      
     };
     TLorentzVector vecs[80];
     for(Int_t kj=0;kj<320;kj++)tpes[kj]=0; 
     Bool_t pass=Esum4>40 || (Esum4>30 && (fmod((skpcnt++),10)<OutOf10-.1));
     pass=Esum[0]>ERecMin[0] || Esum[1]>ERecMin[1] || Esum[2]>ERecMin[2] || Esum[3]>ERecMin[3];
     hardtrig=1;
     if((TrigBitsSelect!=0) && (TrigBits &TrigBitsSelect)==0)hardtrig=0;
     if( (TrigBitsVeto&TrigBits)!=0)hardtrig=0;
  //      printf("led=%d Nhigh4=%d Esum4=%f\n",LEDevent,Nhigh4,Esum4);
  //      if((!LEDevent) &&hardtrig==1 && Nhigh4<70 && Esum4<400  && pass )
  if((!LEDevent) &&hardtrig==1 && Nhigh4<300 && Esum4<3000  && pass )
  {
  for(Int_t k=0;k<4;k++)
  {
  //	      if(Esum[k]>1. && Esum[k]<350.)
  if(Esum[k]>1. && Esum[k]<700.)
  {
  TLorentzVector lv(FourMom(p_Emat[k],0,100,0,100,d.pGeom,2,k+1));
  tpes[nwrds*4]=k+5;
  vecs[nwrds]=lv;
  nwrds++;
  if(i<5)
  {
  p_Emat[k]->Print();
  };
  FitTower::SetSubclu2(false);
  if(i==204)FitTower::SetYPrintLevel(2);      
  Yiqun recon(p_Emat[k],d.pGeom,d.Rgain,d.Rgaincorr,2,k+1);
  FitTower::SetYPrintLevel(0);      
  if(i==204)
  {

  p_Emat[k]->Print();
  for(int jpho=0;jpho<recon.NPh;jpho++)
  {
  printf("Event = %d photon for k:jpho=%d:%d",i,k,jpho);
  printf(" Energy=%f \n",recon.mom(jpho).E());
  HitCluster* aClust=(HitCluster*) recon.photons[jpho].pHitCluster;
  if(aClust)aClust->Print();
  }
  for(int jcc=0;jcc<recon.NRealClusts;jcc++)
  {
  printf("Event = %d Cluster for k=%d\n",i,k);
  recon.clust[jcc].Print();
  } 
  }
  Int_t nz=recon.NPh;
  if(nz>10)nz=10;

  for(Int_t n=0;n<nz;n++)
  {
  vecs[nwrds]=recon.mom(n);
  float myphi=recon.mom(n).Phi();
  tpes[nwrds*4]=305+k;
  nwrds++;
  nphotons++;
  };
  //		  if(Esum[k]>10 && k>1)CheckCluster(&recon);
  //		  if(Esum[k]>10 && k>-1)storeCluster(&recon,1,1.5);
  //		  if(Esum[k]>2 && k>-1)storeCluster(&recon,1,.5);
  //		  if(Esum[k]>1 && k>-1)storerec(&recon,1);
  if(Esum[k]>1 && k>-1)storeCluster(&recon,1,.5);
  };
  };
  };
  //      if(nSavedHits>0)printf("number stored =%d \n",nSavedHits);
  Int_t wcnt=0;
  for(Int_t k=0;k<nwrds;k++)
  {
    pxyzt[wcnt]=vecs[k].Px();
    wcnt++;
    pxyzt[wcnt]=vecs[k].Py();
    wcnt++;
    pxyzt[wcnt]=vecs[k].Pz();
    wcnt++;
    pxyzt[wcnt]=vecs[k].E();
    wcnt++;
  };
  nwrds=nwrds*4;
  */

    /*if(nwrds>0)*/p_out->Fill();
  /*      printf(" i=%d EventN=%d Correction Ratios %f %f %f %f \n adc= %f %f %f %f Energy= %f %f %f %f \n",i,EventN,
          ECorRat[0],ECorRat[1],ECorRat[2],ECorRat[3],
          p_adc[0]->Sum(),p_adc[1]->Sum(),p_adc[2]->Sum(),p_adc[3]->Sum(),
          p_Emat[0]->Sum(),p_Emat[1]->Sum(),p_Emat[2]->Sum(),p_Emat[3]->Sum());
          */
  /*
     for(Int_t k=0;k<4;k++)
     {
     if(p_adc[k]){delete p_adc[k];p_adc[k]=0;};
     if(p_Emat[k]){delete p_Emat[k];p_Emat[k]=0;};
     };
     */
  };

  //end event loop
  p_out->Print();
  p_out->Write();
  /*
     for(int jj=0;jj<5;jj++)
     {
     Emat[jj].Write();
     NHi[jj].Write();
     };
     */
  printf("about to return\n");
  return 0;
};

Int_t AnalTools::Select_script1(Int_t set)
{
  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=true;
  AnalDet[0][1]=true;
  AnalDet[1][0]=false;
  AnalDet[1][1]=false;
  AnalDet[1][2]=false;
  AnalDet[1][3]=false;

  if(set==60)
  {
    p_files=new FilesSet("../../","fpd++ped_20060425.txt",
        "nogach06/fpdgain_2006pp.txt", 
        "root112/calwork/fpdcorr_2006pp_set1.txt",
        "fill_7.txt","Fake","spinpat6","geom.txt");


    p_files->p_fpdrunped()->Directory="./";
    p_files->p_fpdlumi()->Directory="./";
    p_files->Print();
    hnt_name="h112";
    clip=true;
  }


  if(set>=61 && set <=63)
  {
    // Transverse Runs
    if(set==62)
    {
      //			   "root06T/fpdcorr.txt",

      p_files=new FilesSet("../../",
          "run6_export/fpdped.txt",
          "run6_export/fpdgain.txt", 
          "run6_export/fpdcorr_itr10.2day.txt",
          "run6_export/fill.txt",
          "Fake",
          "run6_export/spinpat",
          "run6_export/geom_fy06_trans_survey.txt");
    }
    else
    {
      p_files=new FilesSet("../../",
          "run6_export/fpdped.txt",
          "run6_export/fpdgain.txt", 
          "run6_export/gain_set1/fpdcorr_itr11.txt",
          "run6_export/fill.txt",
          "Fake",
          "run6_export/spinpat",
          "run6_export/geom_fy06_trans_survey.txt");
    };
    p_files->p_fpdrunped()->Directory="../../run6_export/runped";
    //p_files->p_fpdrunped()->Directory="./";
    p_files->p_fpdlumi()->Directory="./";
    p_files->Print();
    hnt_name="h111";
  };
  // longitudinal runs

  if(set>=65 && set<=68)
  {
    p_files=new FilesSet("/star/u/leun/fpd06/",
        "fpd++ped_20060425_corr.txt",
        "fpdgain_2006pp.txt", 
        "gain_set601/fpdcorr_itr5.txt",
        "fill.txt",
        "Fake",
        "spinpat",
        "geom_fy06_trans_survey.txt");

    p_files->p_fpdrunped()->Directory="/star/u/leun/fpd06/runped";
    //p_files->p_fpdrunped()->Directory="./";
    p_files->p_fpdlumi()->Directory="./";
    p_files->Print();
    hnt_name="h111";
  };

  if(set==80)
  {
    //transverse run 8
    p_files=new FilesSet("/star/u/heppel/BatchDir/root08",
        "fpd++ped_20060425_corr.txt",
        "fpd08/fpdgain_2006pp.txt", 
        "fpdcorr.txt_VSH",
        "fpd08/fill.txt",
        "Fake",
        "fpd08/spinpat",
        "fpd08/geom_fy08.txt");
    p_files->p_fpdrunped()->Directory="/star/u/leun/fpd08/runped"; 
    p_files->p_fpdlumi()->Directory="./";  

    p_files->Print();
    hnt_name="h111_08";
    //Look only East for now
    for(Int_t j=0;j<6;j++)AnalDet[1][j]=false;
  }

  char hntp[20];
  strcpy(hntp,hnt_name);
  printf("Calling d()\n");
  dataSet d("./run*.root",p_files,hntp);
  printf("exit d()\n");

  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };
  printf("Call GetEntry(1)\n");
  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };
  Int_t Rnum=d.CurrentRunNumber;
  printf("Run Number=%d \n",Rnum);
  Int_t Bunchid7bit=d.bunchid7bit;
  //AnalTools* at = new AnalTools();
  // FpdMap is a class that figures out maps into data blocks
  FpdMap* pmap=new FpdMap(Rnum);

  TFile* p_OFile=new TFile("OFile.root","recreate");

  // Now return the 7x7  map matrix for EW=2 (west) and NSTP=2 (south)

  TMatrix map_pwn=pmap->GetMatrix(2,3);
  TMatrix map_pws=pmap->GetMatrix(2,4);
  TMatrix map_pen=pmap->GetMatrix(1,1);
  TMatrix map_pes=pmap->GetMatrix(1,2);
  TMatrix n199(14,14);
  n199=199;

  TMatrix map_pbwn=pmap->GetMatrix(2,1);
  TMatrix map_pbws=pmap->GetMatrix(2,2);
  if(AnalDet[1][0])map_pbwn=map_pbwn-n199;
  if(AnalDet[1][1])map_pbws=map_pbws-n199;



  Float_t  EsumWS,EsumWN,EsumES,EsumEN,EsumbWS,EsumbWN;
  EsumbWS=EsumbWN=0;
  Short_t EsumWS_s,EsumWN_s,EsumES_s,EsumEN_s;
  //
  Int_t tpes[320];
  Int_t spin;
  Int_t nphotons;
  Float_t px[10];
  Float_t py[10];
  Float_t pz[10];
  Float_t pE[10];
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t strig;
  Char_t BBcSums[5];
  Float_t adc[98];
  Int_t EventN;
  //Define a Tree
  TTree* p_out = new TTree("p_out","Output Tree");

  p_out->Branch("spin",&(spin),"spin/I");
  p_out->Branch("nphotons",&(nphotons),"nphotons/I");
  p_out->Branch("br_nwrds",&(nwrds),"nwrds/I");
  p_out->Branch("br_types",&(tpes),"tpes[nwrds]/I");
  p_out->Branch("br_pxyzt",&(pxyzt),"pxyzt[nwrds]/F");
  p_out->Branch("br_Rnum",&(Rnum),"Rnum/I");
  p_out->Branch("br_Bunchid7bit",&(Bunchid7bit),"Bunchid7bit/I");
  p_out->Branch("br_BBcSums",BBcSums,"BBcSums[5]/b");

  if(WriteADC)p_out->Branch("br_adc",&(adc),"adc[98]/F");
  p_out->Branch("br_EventN",&(EventN),"EventN/I");

  //
  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  if(NumberEventsToProcess>0 && NumberEventsToProcess>nentries)nentries=NumberEventsToProcess;
  // event loop
  Int_t skpcntEN=0;
  Int_t skpcntES=0;
  Int_t skpcntWN=0;
  Int_t skpcntWS=0;
  map_pen.Print();
  map_pes.Print();
  printf("nentriies=%d\n",nentries);
  for(Int_t i=0;i<nentries;i++)
  {
    Int_t nbytes=d.GetEntry(i);  

    Rnum=d.CurrentRunNumber;

    EventN=d.event+(d.CurrentSegNumber)*100000;

    Bunchid7bit=d.bunchid7bit;
    pcnt++;
    if(pcnt>100){printf("cnt=%d \n",i);pcnt=1;}; //print every 100 events
    Double_t NearUnique;

    spin=(d.BlueSpin+1)+(d.YellowSpin+1)/2;
    if(d.kicked)spin+=29;
    if(d.BlueSpin*d.YellowSpin==0)spin=40;
    Trigger trig(d.Bbcl1);
    for(Int_t ki=0;ki<5;ki++)BBcSums[ki]=trig.BBcSums[ki];

    TMatrix tmp0wn= d.dMatrix(&(map_pwn),d.Fpdwns);
    TMatrix tmp0ws= d.dMatrix(&(map_pws),d.Fpdwns);
    TMatrix tmp0en= d.dMatrix(&(map_pen),d.Fpdens);
    TMatrix tmp0es= d.dMatrix(&(map_pes),d.Fpdens);

    TMatrix tmp0bwn= d.dMatrix(&(map_pbwn),d.fpdwadc);
    TMatrix tmp0bws= d.dMatrix(&(map_pbws),d.fpdwadc);
    for(Int_t ji=0;ji<49;ji++)
    {
      Int_t col=ji%7;
      Int_t row=ji/7;
      adc[ji]=tmp0en(row,col);
      adc[ji+49]=tmp0es(row,col);
    };
    //      if(tmp0en(4,6)>10)tmp0en.Print();
    if(AnalDet[0][3]){
      tmp0bws(10,10)=0.;}
    hardtrig=1;
    TMatrix energydata_2_4(tmp0ws.GetNrows(),tmp0ws.GetNcols());
    TMatrix energydata_2_3(tmp0wn.GetNrows(),tmp0wn.GetNcols());
    TMatrix energydata_2_2(tmp0bws.GetNrows(),tmp0bws.GetNcols());
    TMatrix energydata_2_1(tmp0bwn.GetNrows(),tmp0bwn.GetNcols());
    TMatrix energydata_1_2(tmp0es.GetNrows(),tmp0es.GetNcols());
    TMatrix energydata_1_1(tmp0en.GetNrows(),tmp0en.GetNcols());

    if(AnalDet[1][3])energydata_2_4=d.Em(&tmp0ws,2,4);  // WS energy data for event 
    if(AnalDet[1][2])energydata_2_3=d.Em(&tmp0wn,2,3);  // WN energy data for event 
    if(AnalDet[1][1])energydata_2_2=d.Em(&tmp0bws,2,2);  // bWS energy data for event 
    if(AnalDet[1][0])energydata_2_1=d.Em(&tmp0bwn,2,1);  // bWN energy data for event 
    if(AnalDet[0][0])energydata_1_1=d.Em(&tmp0en,1,1);  // EN energy data for event 
    if(AnalDet[0][1])energydata_1_2=d.Em(&tmp0es,1,2);  // ES energy data for event 

    EsumWS=energydata_2_4.Sum(); // sum over 36 WS channels
    EsumWN=energydata_2_3.Sum(); // sum over 36 WN channels
    EsumES=energydata_1_2.Sum(); // sum over 49 ES channels
    EsumEN=energydata_1_1.Sum(); // sum over 49 EN channels
    EsumbWN=energydata_2_1.Sum();//sum over 196 WS large channels;
    EsumbWS=energydata_2_2.Sum();//sum over 196 WS large channels;
    hardtrig=0;
    if(EsumES>50 || EsumEN>50)hardtrig=1;

    //count the W channels with more than 0.25 GeV Energy
    Int_t highcnt4=0;
    Int_t highcnt3=0;
    if(AnalDet[1][3])Int_t highcnt4= NTowersAbove(&energydata_2_4,.25);
    if(AnalDet[1][2])Int_t highcnt3= NTowersAbove(&energydata_2_3,.25);
    EnergySum=0;
    nwrds=0;
    nphotons=0;
    Bool_t WB=EsumbWN>2 || EsumbWS>2;
    TLorentzVector vecs[80];
    for(Int_t kj=0;kj<40;kj++)tpes[kj]=0; 
    if(hardtrig==1)
    {
      if(highcnt3+highcnt4<40)
      {
        if(EsumWS>30 || (EsumWS>20 && (fmod((skpcntWS++),10)<OutOf10-.1))|| WB){
          TLorentzVector lv(FourMom(&energydata_2_4,0,100,0,100,d.pGeom,2,4));
          tpes[nwrds*4]=8;
          vecs[nwrds]=lv;
          nwrds++;
          Yiqun recon(&energydata_2_4,d.pGeom,d.Rgain,d.Rgaincorr,2,4);
          nphotons+=recon.NPh;
          for(Int_t n=0;n<recon.NPh && n<4;n++)
          {
            vecs[nwrds]=recon.mom(n);
            tpes[nwrds*4]=308;
            nwrds++;
          };

          if(tmp0bws.Sum()!=0 && EsumbWS>2.)
          {
            TLorentzVector lvws(FourMom(&energydata_2_2,0,100,0,100,d.pGeom,2,2));
            tpes[nwrds*4]=6;
            vecs[nwrds]=lvws;
            nwrds++;
            Yiqun reconWS(&energydata_2_2,d.pGeom,d.Rgain,d.Rgaincorr,2,2);
            nphotons+=reconWS.NPh;
            for(Int_t n=0;n<reconWS.NPh && n<4;n++)
            {
              vecs[nwrds]=reconWS.mom(n);
              tpes[nwrds*4]=306;
              nwrds++;
            };
          };
        };
        if(EsumWN>30 || (EsumWN>20 && (fmod((skpcntWN++),10)<OutOf10-.1))|| WB){
          TLorentzVector lv(FourMom(&energydata_2_3,0,100,0,100,d.pGeom,2,3));
          tpes[nwrds*4]=7;
          vecs[nwrds]=lv;
          nwrds++;
          Yiqun recon(&energydata_2_3,d.pGeom,d.Rgain,d.Rgaincorr,2,3);
          nphotons+=recon.NPh;
          for(Int_t n=0;n<recon.NPh && n<4;n++)
          {
            vecs[nwrds]=recon.mom(n);
            tpes[nwrds*4]=307;
            nwrds++;
          };
          if(tmp0bwn.Sum()!=0 && EsumbWN>2.)
          {

            TLorentzVector lvwn(FourMom(&energydata_2_1,0,100,0,100,d.pGeom,2,1));
            tpes[nwrds*4]=5;
            vecs[nwrds]=lvwn;
            nwrds++;
            Yiqun reconWN(&energydata_2_1,d.pGeom,d.Rgain,d.Rgaincorr,2,1);
            nphotons+=reconWN.NPh;
            for(Int_t n=0;n<reconWN.NPh && n<4;n++)
            {
              vecs[nwrds]=reconWN.mom(n);
              tpes[nwrds*4]=305;
              nwrds++;
            };
          };
        };
        if(EsumES>35 || (EsumES>20 && (fmod((skpcntES++),20)<OutOf10-.1))){
          TLorentzVector lv(FourMom(&energydata_1_2,0,100,0,100,d.pGeom,1,2));
          tpes[nwrds*4]=2;
          vecs[nwrds]=lv;
          nwrds++;
          Yiqun recon(&energydata_1_2,d.pGeom,d.Rgain,d.Rgaincorr,1,2);
          nphotons+=recon.NPh;
          for(Int_t n=0;n<recon.NPh && n<4;n++)
          {
            vecs[nwrds]=recon.mom(n);
            tpes[nwrds*4]=302;
            nwrds++;
          };
        };
        if(EsumEN>35 || (EsumEN>20 && (fmod((skpcntEN++),20)<OutOf10-.1))){
          TLorentzVector lv(FourMom(&energydata_1_1,0,100,0,100,d.pGeom,1,1));
          tpes[nwrds*4]=1;
          vecs[nwrds]=lv;
          nwrds++;
          Yiqun recon(&energydata_1_1,d.pGeom,d.Rgain,d.Rgaincorr,1,1);
          nphotons+=recon.NPh;
          for(Int_t n=0;n<recon.NPh && n<4;n++)
          {
            vecs[nwrds]=recon.mom(n);
            tpes[nwrds*4]=301;
            nwrds++;
          }; 
        };
      };
    };
    Int_t wcnt=0;
    for(Int_t k=0;k<nwrds;k++)
    {
      pxyzt[wcnt]=vecs[k].Px();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Py();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Pz();
      wcnt++;
      pxyzt[wcnt]=vecs[k].E();
      wcnt++;
    };
    nwrds=nwrds*4;
    if(nwrds>0)p_out->Fill();
    //p_out->Fill();
  };

  //end event loop
  p_out->Print();
  p_out->Write();
  return 0;
};

Bool_t AnalTools::SimScript2(char* path)
{
  FilesSet* p_files=new FilesSet("../",
      "run6_export/fpdped.txt",
      "run6_export/fpdgain.txt", 
      "run6_export/fpdcorr.txt",
      "run6_export/fill.txt",
      "Fake",
      "run6_export/spinpat",
      "geomMC.txt");


  p_files->p_fpdrunped()->Directory="../run6_export/runped";
  p_files->Print();

  Geom* p_Geom=new Geom(p_files);

  Sim* sh=new Sim(path,p_Geom);
  TH1F* Mass=new TH1F("Mass","Mass",140,.0,1.2);
  TH2F* MPivsMall=new TH2F("MPivsMall","MPivsMall",140,0.,1.2,140,0.,1.2);
  TH2F* EWSbvsEWNb=new TH2F("EWSbvsEWNb","EWSbvEWNb",200,0,50.,200,0.,50.);
  TH2F* ZvsM=new TH2F("ZvsM","ZvsM",140,0,1.4,100,0,1);
  TH2F* EN = new TH2F("EN","EN",14,0,14,14,0,14);
  TH2F* ENdet=0;
  TH2F* ENdetPi0=0;
  pSim=sh;
  Int_t nevents=sh->Input->GetEntries();
  //  TCanvas* c1=new TCanvas("c1","c1",600,600);
  //  TCanvas* c2=new TCanvas("c2","c2",600,600);

  Int_t skip=0; 
  printf("nevents=%d \n",nevents);

  for(Int_t i=0;i<nevents;i++)
  {
    if(skip==0)printf("count=%d\n",i);
    skip++;
    if(skip>100)skip=0;
    sh->Input->GetEntry(i);
    sh->SelectTest=1;
    sh->GenADC(-.00001,1000.);
    sh->Fillshape();
    TMatrix pi0mat=(*sh->newene[0][1]);
    sh->SelectTest=0;
    sh->GenADC(-.00001,1000.);

    Float_t En;
    Float_t Ewest=sh->newene[0][1]->Sum()+sh->newene[1][1]->Sum()+
      sh->newene[2][1]->Sum()+sh->newene[3][1]->Sum();
    Float_t Eeast=sh->newene[0][0]->Sum()+sh->newene[1][0]->Sum();
    TMatrix tm00(*sh->newene[0][0]);
    TMatrix tm01(*sh->newene[0][1]);
    TH2F dEN(tm01);
    *EN=*EN+dEN;
    TMatrix tm11(*sh->newene[1][1]);
    TMatrix tm21(*sh->newene[2][0]);
    TMatrix tm31(*sh->newene[3][0]);
    EWSbvsEWNb->Fill(tm01.Sum(),tm11.Sum());
    if(tm11.Sum()>1. && tm01.Sum()>5.)
    {
      //	  tm01.Draw("box");


      Yiqun recon(&tm01,p_Geom,2,1);
      if(recon.NPh>=2)
      { 
        /*	      printf("Epi=%f \n",pi0mat.Sum());
                      printf("EAll=%f \n",sh->newene[0][1]->Sum());
                      */
        for(Int_t i1=0;i1<recon.NPh;i1++)
        {
          for(Int_t i2=0;i2<i1;i2++)
          {
            TLorentzVector v1a=recon.mom(i1);
            TLorentzVector v1b=recon.mom(i2);
            char cstr[20];
            //		      printf(" m%d_%d=%f",i1,i2,(v1a+v1b).Mag());
          };
        };
        //	      printf("\n");
        Int_t ncl=recon.NRealClusts;
        for(Int_t ii=0;ii<ncl;ii++)
        { 
          Int_t catag=recon.clust[ii].catag;
          Int_t Ntow=recon.clust[ii].tow->GetEntries();
          //		  printf("Ntow=%d\n",Ntow);
          TObjArray* P_Oarray=recon.clust[ii].tow;
          TIterator* iter=P_Oarray->MakeIterator();
          TowerFPD* ftow;
          /*
             while( ftow=(TowerFPD*)iter->Next())
             {
             Int_t col=ftow->col;
             Int_t row=ftow->row;
             Float_t adcval=tm01(row-1,col-1);
             printf(
             "cluster:%d catag=%d col=%d row=%d energy=%f \n", ii,catag,col,row,ftow->energy);
             };
             */
        };
        //	      printf("tm01.E()=%f tm11.E()=%f \n",tm01.Sum(),tm11.Sum());

        /*
           c1->cd();pi0mat.Draw("box");
           if(ENdet!=0)delete ENdet;
           if(ENdetPi0!=0)delete ENdetPi0;

           ENdetPi0=new TH2F("ENdetPi0","ENDetPi0",14,0,14,14,0,14);
           ENdet=new TH2F("ENdet","ENdet",14,0,14,14,0,14);
           TH2F tmph(pi0mat);
         *ENdetPi0= *ENdetPi0+tmph;
         TH2F tmph2(tm01);
         *ENdet= *ENdet+tmph2;
         ENdetPi0->SetLineColor(2);
         ENdetPi0->GetZaxis()->SetRangeUser(0.,15.);
         ENdet->Draw("box");
         ENdetPi0->Draw("boxsame");
         c1->Print("SimEvt.eps");
         c1->Update();
         */
        //	      c2->cd();tm01.Draw("box");
        //	      c2->Update();
        /*	      
                      for(Int_t j=0;j<recon.NPh;j++)
                      {
                      printf("x=%f y=%f e=%f \n",recon.photons[j].xPos/5.8,recon.photons[j].yPos/5.8,recon.photons[j].energy);
                      };
                      */
        if(recon.NPh==2)
        {
          Float_t z,e1,e2;
          e1=recon.mom(0).E();
          e2=recon.mom(1).E();
          Float_t Eratio=(e1+e2)/tm01.Sum();
          if(Eratio>.7)
          {
            //	      tm01.Draw("lego");
            Float_t mpi=(recon.mom(0)+recon.mom(1)).Mag();
            Mass->Fill(mpi);
            //	      printf("mass=%f \n",mpi);
            z=-1;
            if(e1+e2>0)z=fabs((e1-e2)/(e1+e2));
            ZvsM->Fill(mpi,z);
            Yiqun reconPi(&pi0mat,p_Geom,2,1);
            TLorentzVector vpisum(0.,0.,0.,0.);
            for(Int_t picnt=0;picnt<reconPi.NPh;picnt++)vpisum+=reconPi.mom(picnt);
            MPivsMall->Fill(vpisum.Mag(),mpi);
          };
        };
        //	  c1->Update();
        //	  	  system("sleep 2");
        //  printf("_\n");
      };
    };
  };
  Mass->Draw();
  //  c1->Update();
  EWSbvsEWNb->Draw();
  TFile tfile("simout.root","update");
  Mass->Write("MassAll8");
  EWSbvsEWNb->Write("EWSbvsEWNbAll8");
};

Bool_t AnalTools::SimScript(char* path)
{
  FilesSet* p_files=new FilesSet("../",
      "run6_export/fpdped.txt",
      "run6_export/fpdgain.txt", 
      "run6_export/fpdcorr.txt",
      "run6_export/fill.txt",
      "Fake",
      "run6_export/spinpat",
      "geomMC.txt");

  p_files->p_fpdrunped()->Directory="../run6_export/runped";
  p_files->Print();
  Geom* p_Geom=new Geom(p_files);
  Sim* sh=new Sim(path,p_Geom);
  pSim=sh;
  Int_t nevents=sh->Input->GetEntries();
  TH1F* hadc11=new TH1F("hadc11","hadc11",200,0,200);
  TH1F* hene11=new TH1F("hene11","hene11",200,0,100);
  TH1F* hene12=new TH1F("hene12","hene12",200,0,100);
  TH1F* hene15=new TH1F("hene15","hene13",200,0,100);
  TH1F* hene21=new TH1F("hene21","hene21",200,0,100);
  TH1F* hene22=new TH1F("hene22","hene22",200,0,100);
  TH1F* hene23=new TH1F("hene23","hene23",200,0,100);
  TH1F* hene24=new TH1F("hene24","hene24",200,0,100);
  TH2F* heneEW=new TH2F("heneEW","heneEW",200,0,100,200,0,100);
  TH1F* hmass11=new TH1F("hmass11","hmass11",60,0,1.2);
  TH1F* hmass21=new TH1F("hmass21","hmass21",60,0,1.2);
  TH1F* hmass23=new TH1F("hmass23","hmass23",60,0,1.2);
  TH1F* hemass11=new TH1F("hemass11","hemass11",60,0,1.2);
  TH1F* hnphot11=new TH1F("hnphot11","hnphot11",10,0,10);
  TH1F* hnphot21=new TH1F("hnphot21","hnphot21",10,0,10);
  TH1F* hnphot23=new TH1F("hnphot23","hnphot23",10,0,10);
  TH2F* GenMassPzph=new TH2F("GenMassPzph","GenMassPzph",60,0,1.2,100,-100,100);  Double_t Pi=TMath::Pi();
  TH2F* MothPi0EtaPhi=new TH2F("MothPi0EtaPhi","MothPi0EtaPhi",
      120,-5,7,628,-Pi/2,3*Pi/2);
  TH2F* PhotEtaPhi=new TH2F("PhotEtaPhi","PhotEtaPhi",
      120,-5,7,628,-Pi/2,3*Pi/2);
  TH2F* XYrec11=new TH2F("XYrec11","XYrec11",100,-1,9,100,-1,9);
  TH2F* XYgen11e=new TH2F("XYgen11e","XYgen11e",100,-1,9,100,-1,9);
  TH2F* XYgen11=new TH2F("XYgen11","XYgen11",100,-1,9,100,-1,9);
  TH2F* XYgenb11=new TH2F("XYgenb11","XYgenb11",100,-1,9,100,-1,9);
  TH2F* XYdif11=new TH2F("XYdif11","XYdif11",100,-5,5,100,-5,5);
  TH1F* HZrec=new TH1F("HZrec","HZrec",20,0,1);
  TH1F* HZgen=new TH1F("HZgen","HZgen",20,0,1);
  TCanvas* c1=new TCanvas("c1","c1",600,600);
  c1->Divide(1,2);
  Int_t skip=0; 

  printf("nevents=%d \n",nevents);
  for(Int_t i=0;i<nevents;i++)
  {
    if(skip==0)printf("count=%d\n",i);
    skip++;
    if(skip>100)skip=0;
    sh->Input->GetEntry(i);
    sh->FillHighFour();
    for( Int_t ii=0;ii<4;ii++){
      if(ii>=0){
        Int_t hid=sh->HighFour[ii];	
        if(hid>=0){
          //printf("trk%d id =%d mom=(%f,%f,%f)  mother: trk=%d  mass=%f \n",ii, hid,(sh->p)[hid][0],(sh->p)[hid][1],(sh->p)[hid][2],sh->mo[hid][0]-1,sh->m[sh->mo[hid][0]-1]);
          GenMassPzph->Fill(sh->m[sh->mo[hid][0]-1],(sh->p)[sh->mo[hid][0]-1][2]);
          Bool_t newmother=true;
          for(Int_t im=0;im<ii;im++){
            Int_t motherid=sh->mo[sh->HighFour[im]][0]-1;
            if(sh->mo[hid][0]-1==motherid)newmother=false;
          };
          if(newmother)
          {
            Int_t motherid=sh->mo[hid][0]-1;
            if(sh->e[motherid]>10)GenMassPzph->Fill(sh->m[motherid],(sh->p)[motherid][2]);
            //look at pi0;
            if(fabs(sh->m[motherid]-.135)<.02 && sh->e[motherid]>13)
            {
              Float_t phi=fmod(2.5*Pi+sh->phi[motherid],2*Pi)-.5*Pi;
              MothPi0EtaPhi->Fill(sh->eta[motherid],phi);
            };
          };
        };
      };
    };

    sh->GenADC();
    sh->Fillshape();

    //      sh->FillHits(-.00001,.00001);// photons only
    Float_t En;
    Float_t Ewest=sh->newene[0][1]->Sum()+sh->newene[1][1]->Sum()+
      sh->newene[2][1]->Sum()+sh->newene[3][1]->Sum();
    Float_t Eeast=sh->newene[0][0]->Sum()+sh->newene[1][0]->Sum();
    heneEW->Fill(Eeast,Ewest);
    TMatrix tm00(*sh->newene[0][0]);
    TMatrix tm01(*sh->newene[0][0]);
    TMatrix tm21(*sh->newene[2][0]);
    En=tm00.Sum();
    hadc11->Fill(sh->adc[0][0]->Sum());
    //      sh->FillHits();
    hene11->Fill(En);
    hene12->Fill(sh->newene[1][0]->Sum());
    hene22->Fill(sh->newene[1][1]->Sum());
    hene24->Fill(sh->newene[3][1]->Sum());
    Float_t Eps;
    TMatrix tm40=(*sh->newene[4][0]);
    Eps=tm40.Sum();
    hene15->Fill(Eps);
    En+=Eps;
    if(Ewest>10.)
    {
      printf("ok\n");
      Yiqun recon(&tm21,p_Geom,3,2);
      Bool_t cut=true;
      for(Int_t k=0;k<recon.NPh;k++){
        if(recon.mom(k).E()>.001){
          Float_t phi=fmod(2.5*Pi+recon.mom(k).Phi(),2*Pi)-.5*Pi;
          PhotEtaPhi->Fill(recon.mom(k).Eta(),phi);
        };
      };
      Float_t x=-(recon.ph_coord_lab(0)[0]-*recon.p_geom->xOffset(3,2))
        /(*recon.p_geom->FpdTowWid(3,2));
      Float_t y=-(recon.ph_coord_lab(0)[1]-*recon.p_geom->yOffset(3,2))
        /(*recon.p_geom->FpdTowWid(3,2));  
      XYrec11->Fill(x,y);    

      TVector3 v3= pSim->ProjectTrack(pSim->HighFour[0],
          *recon.p_geom->ZFPD(3,2));
      Float_t x0,y0,e0;
      e0=pSim->e[pSim->HighFour[0]];
      x0=-(v3.X()-*recon.p_geom->xOffset(3,2))
        /(*recon.p_geom->FpdTowWid(3,2));
      y0=-(v3.Y()-*recon.p_geom->yOffset(3,2))
        /(*recon.p_geom->FpdTowWid(3,2));
      // printf("x =%f x0=%f y=%f  y0=%f \n",x,x0,y,y0);
      XYgen11->Fill(x0,y0);
      Int_t io=pSim->da[ 
        pSim->mo[pSim->HighFour[0]][0]-1][0]-1;
      if(io==pSim->HighFour[0])
        io=pSim->da[ pSim->mo[pSim->HighFour[0]][0]-1][1]-1;

      TVector3 v3b= pSim->ProjectTrack(io,
          *recon.p_geom->ZFPD(1,1));
      Float_t x1,y1,e1;
      e1=pSim->e[io];

      x1=-(v3b.X()-*recon.p_geom->xOffset(1,1))
        /(*recon.p_geom->FpdTowWid(1,1));
      y1=-(v3b.Y()-*recon.p_geom->yOffset(1,1))
        /(*recon.p_geom->FpdTowWid(1,1));
      //	     printf("x =%f x0=%f y=%f  y0=%f \n",x,x0,y,y0);
      if(x0 < 0 || x0>7  || y0<0 || y0 > 7)cut=false;
      if(x1 < 0 || x1>7  || y1<0 || y1 > 7)cut=false;

      XYgenb11->Fill(x1,y1);
      if(e0+e1>0 && cut)HZgen->Fill(fabs((e0-e1)/(e0+e1)));
      XYdif11->Fill(x-x0,y-y0);

      if(recon.NPh==2)
      {
        printf("x0=%f y0=%f x1=%f y1=%f \n",x0,y0,x1,y1);
        c1->cd(1);
        XYgen11e->Reset();
        XYgen11e->Fill(x0,y0,e0);
        XYgen11e->Fill(x1,y1,e1);
        XYgen11e->Draw("box");
      };


      if(recon.NPh==2)
      {
        c1->cd(2);
        pSim->newene[0][0]->Draw("box");
        c1->Update();
        //	    recon.FittedMat().Draw("box");
        TH2F fmat(recon.FittedMat());
        fmat.SetLineColor(2);
        fmat.Draw("samebox");
        c1->cd(1);
        fmat.Draw("boxtext");
        c1->Update();
        Float_t WD=*p_Geom->FpdTowWid(1,1);
        printf("xpos0=%f ypos0=%f xpos1=%f ypos1=%f \n",
            recon.photons[0].xPos/WD,recon.photons[0].yPos/WD,recon.photons[1].xPos/WD,recon.photons[1].yPos/WD);
        //	    Int_t sysret=system("sleep 2");
        printf("_\n");
        //	    sysret=system("sleep 2");
      };

      hnphot11->Fill(recon.NPh);
      if(recon.NPh==2)
      {
        Float_t mass=(recon.mom(0)+recon.mom(1)).Mag();
        Float_t efit=(recon.mom(0)+recon.mom(1)).E();
        Float_t er0,er1;
        er0=recon.mom(0).E();
        er1=recon.mom(1).E();
        if(er0+er1>0 && cut)HZrec->Fill(fabs((er0-er1)/(er0+er1)));
        hmass11->Fill(mass);
        hemass11->Fill(mass*En/efit);
      };
    };
    En=tm01.Sum(); 
    hene21->Fill(En);
    if(En>3)
    {	  
      Yiqun recon(&tm01,p_Geom,2,1);
      hnphot21->Fill(recon.NPh);
      for(Int_t k=0;k<recon.NPh;k++){
        if(recon.mom(k).E()>1.){
          Float_t phi=fmod(2.5*Pi+recon.mom(k).Phi(),2*Pi)-.5*Pi;
          PhotEtaPhi->Fill(recon.mom(k).Eta(),phi);
        };
      };

      if(recon.NPh==2)
      {
        Float_t mass=(recon.mom(0)+recon.mom(1)).Mag();
        Float_t efit=(recon.mom(0)+recon.mom(1)).E();
        hmass21->Fill(mass);
        //	      hemass11->Fill(mass*En/efit);
        if(recon.NPh==2 && false)
        {
          c1->cd(2);
          pSim->newene[0][1]->Draw("box");
          c1->Update();
          //		  system("sleep 2");
          c1->cd(1);
          recon.FittedMat().Draw("box");
          c1->Update();
          //		  system("sleep 2");
          printf("-\n");
        };
      };

    };
    En=tm21.Sum(); 
    hene23->Fill(En);
    if(En>5)
    {	  
      Yiqun recon(&tm21,p_Geom,2,3);
      hnphot23->Fill(recon.NPh);
      if(recon.NPh==2)
      {
        Float_t mass=(recon.mom(0)+recon.mom(1)).Mag();
        Float_t efit=(recon.mom(0)+recon.mom(1)).E();
        hmass23->Fill(mass);
        //	      hemass11->Fill(mass*En/efit);
      };

    };
  };

  hmass11->Draw();

};

Int_t AnalTools::readqM(char* filelist,Int_t set,FilesSet* p_files)
{
  Qt qt(p_files);
  Geom* p_geom=new Geom(p_files);
  TFile* out=new TFile(OutFileName,"recreate");
  pout=new poutTree(filelist);
  pout->GetEntry(0);
  Int_t RunNum=pout->Rnum;
  p_files->Print();
  std::cout<<"path="<<p_files->p_fpdgaincorr()->path<<"\n";
  CalibStr gcor(RunNum,(const char*) p_files->p_fpdgaincorr()->path);
  CalibStr gain(RunNum,(const char*) p_files->p_fpdgain()->path);
  gcor.Print();
  Cell* cell[34][17][4];
  TObjArray mcells(2400,0);
  for(Int_t i=0;i<4;i++)
  {
    Int_t up=34;
    if(i>1)up=24;
    for(Int_t j=0;j<34;j++)
    {
      for(Int_t k=0;k<17;k++)
      {
        //-----------------------------

        cell[j][k][i]=0;
        if(j>=up || k>=up/2)continue;
        if((i<2)&&(k<8) && (j>8) && (j<25))continue;
        if((i>1)&&(k<5)&&(j>6)&&(j<17))continue;

        char nam[40];

        sprintf(nam,"Cellr%d_c%d_%d",j,k,i);
        std::cout<<"about to make: "<<nam<<"\n";
        out->cd();
        Cell* pcell=new Cell(2,i+1,j+1,k+1,RunNum,p_geom,&gain,&gcor);
        pcell->SetName(nam);
        cell[j][k][i]=pcell;
        mcells.Add(pcell);	
        sprintf(nam,"Mcellr%d_c%d_%d",j,k,i);
        TH1F* pmcell=new TH1F(nam,nam,240,0.,2.4);

        pcell->p_Mcell=pmcell;

        sprintf(nam,"dxdyMr%d_c%d_%d",j,k,i);
        TH3F* p_dxdyM=new TH3F(nam,nam,13,k-6.5,k+6.5,13,j-6.5,j+6.5,60,0.,.8);
        pcell->p_dxdyM=p_dxdyM;

        sprintf(nam,"dxdyEr%d_c%d_%d",j,k,i);
        TH3F* p_dxdyE=new TH3F(nam,nam,13,k-6.5,k+6.5,13,j-6.5,j+6.5,20,0.,100);
        pcell->p_dxdyE=p_dxdyE;

        sprintf(nam,"dxdyERr%d_c%d_%d",j,k,i);
        TH3F* p_dxdyER=new TH3F(nam,nam,13,k-6.5,k+6.5,13,j-6.5,j+6.5,40,0.,2.);
        pcell->p_dxdyER=p_dxdyER;

        sprintf(nam,"Exyr%d_c%d_%d",j,k,i);
        TH1F* pexy=new TH1F(nam,nam,120,0.,120.);
        pcell->p_EphXY=pexy;
        sprintf(nam,"MvsZr%d_c%d_%d",j,k,i);
        TH2F* pmvz=new TH2F(nam,nam,120,0,1.2,10,0,1.);
        pcell->p_MvsZ=pmvz;
        sprintf(nam,"MvsDr%d_c%d_%d",j,k,i);
        TH2F* pmvd=new TH2F(nam,nam,120,0,1.2,20,0,20*3.82);
        pcell->p_MvsD=pmvd;

        sprintf(nam,"MvsEr%d_c%d_%d",j,k,i);
        TH2F* pmve=new TH2F(nam,nam,120,0,1.2,20,0,100);
        pcell->p_MvsE=pmve;

        sprintf(nam,"MvsY%d_c%d_%d",j,k,i);
        TH2F* pmvy=new TH2F(nam,nam,120,0,1.2,30,2.,5.);
        pcell->p_MvsY=pmvy;

        sprintf(nam,"Nphthis%d_c%d_%d",j,k,i);
        TH1F* pnphthis=new TH1F(nam,nam,20,0,20);
        pcell->p_Nphthis=pnphthis;

        sprintf(nam,"NphAll%d_c%d_%d",j,k,i);
        TH1F* pnphall=new TH1F(nam,nam,20,0,20);
        pcell->p_NphAll=pnphall;

        sprintf(nam,"DetE%d_c%d_%d",j,k,i);
        TH1F* DetE=new TH1F(nam,nam,100,0,100);
        pcell->p_DetE=DetE;

        sprintf(nam,"FmsE%d_c%d_%d",j,k,i);
        TH1F* p_FmsE=new TH1F(nam,nam,100,0,100);
        pcell->p_FmsE=p_FmsE;

        sprintf(nam,"CellD%d_c%d_%d",j,k,i);
        pcell->InitTree(nam);

        //-------------------------
      };
    };
  };
  TH2F* qtHist[4];
  TH2F* qtHistLed[4];
  TFile adcroot("adc.root");
  if(adcroot.IsOpen())
  {

    qtHist[0]=(TH2F*) adcroot.GetKey("qtHistLN")->ReadObj();
    qtHist[1]=(TH2F*) adcroot.GetKey("qtHistLS")->ReadObj();
    qtHist[2]=(TH2F*) adcroot.GetKey("qtHistSN")->ReadObj();
    qtHist[3]=(TH2F*) adcroot.GetKey("qtHistSS")->ReadObj();

    qtHistLed[0]=(TH2F*) adcroot.GetKey("qtHistLedLN")->ReadObj();
    qtHistLed[1]=(TH2F*) adcroot.GetKey("qtHistLedLS")->ReadObj();
    qtHistLed[2]=(TH2F*) adcroot.GetKey("qtHistLedSN")->ReadObj();
    qtHistLed[3]=(TH2F*) adcroot.GetKey("qtHistLedSS")->ReadObj();

  };
  TFile nfl("CellSet.root");
  TObjArray dnoa(1,0);
  Bool_t ProcessAll;
  TObjArray* noa;
  if(nfl.IsOpen())
  {
    printf("nfl open\n");
    ProcessAll=false;
    noa=(TObjArray*) nfl.GetKey("NewCells")->ReadObj();
  }
  else
  {
    printf("nfl not open\n");
    ProcessAll=true;
    noa=&dnoa;
  };
  out->cd();
  Cell* p_cl;
  TIter next(noa);
  while(p_cl=(Cell*) next())
  {
    if(p_cl->Iew==2 && p_cl->Instb<5 && p_cl->Instb>0)
    {
      Int_t k1=p_cl->Row1-1;
      Int_t k2=p_cl->Col1-1;
      Int_t k3=p_cl->Instb-1;
      cell[k1][k2][k3]=p_cl;
      p_cl->p_Mcell->Reset();
      p_cl->p_EphXY->Reset();
      p_cl->p_MvsZ->Reset();
      p_cl->p_MvsD->Reset();
      p_cl->p_MvsE->Reset();
      p_cl->p_dxdyM->Reset();
      p_cl->p_dxdyE->Reset();
      p_cl->p_dxdyER->Reset();
      p_cl->p_Nphthis->Reset();
      p_cl->p_NphAll->Reset();
      p_cl->p_DetE->Reset();
      p_cl->p_FmsE->Reset();
      mcells.Add(p_cl);
    };
  };

  for(Int_t kk=1;kk<9;kk++)pout->MinEnergy[kk]=4;
  pout->MinEnergy[5]=2.;
  pout->MinEnergy[6]=2.;

  Int_t nentries=pout->nentries;
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;  
  if(NumberEventsToProcess>pout->nentries)nentries=pout->nentries;
  printf("start loop of %d entries \n",nentries);
  p_geom->FMSGeom=true;
  Float_t EinNSTB[4];
  TCanvas* c1=new TCanvas("c1","c1",600,600);
  gStyle->SetPalette(1);
  for(Int_t iev=0;iev<nentries;iev++)
  {
    int evprnt=1000;

    if( (iev%evprnt) ==0)printf( "%d Events Read \n ",iev);
    Int_t npo=pout->GetEntry(iev);
    /* comment out below      
       c1->Clear();
       c1->Divide(2,2);
       for(int dj=1;dj<5;dj++)
       {
       if(pout->FillFMSADC(dj).Sum()>5)
       {
       c1->cd(dj);
       pout->DrawFMSADC(2,dj,p_geom);
       c1->GetPad(dj)->SetLogz();
       c1->Update();
       };
       };
       int iret= system("sleep 2");
       */
    //comment out above ^^^^^^^^^  
    RunNum=pout->Rnum; 
    Float_t mass;

    for(int zi=0;zi<4;zi++)EinNSTB[zi]=0;
    TIter nex(pout->dlist);
    while(LVec* d_=(LVec*) nex())
    {
      Int_t nstb_=pout->WhatNSTB(d_);
      if(nstb_>=0 && nstb_<4)
      {
        EinNSTB[nstb_]=d_->E(); //save total NSTB det energy
      };
    };
    TIter next(pout->vlist);


    while(LVec* v_=(LVec*) next())
    {
      int _ew,_nstb,_row,_col;
      _ew=pout->WhatEW(v_);
      _nstb= pout->WhatNSTB(v_);
      _row=pout->WhatRow0(v_,p_geom,0.);
      _col=pout->WhatCol0(v_,p_geom,0.);  


      if(_ew==2 &&(_nstb)>0)
      {
        if( _row>=0 &&  _col>=0)
        {
          if(cell[_row][_col][_nstb-1]){

            cell[_row][_col][_nstb-1]->p_Nphthis->Fill(pout->GetNPhotVec(2,_nstb));
            cell[_row][_col][_nstb-1]->p_NphAll->Fill(pout->vlist->GetEntries());
            cell[_row][_col][_nstb-1]->p_DetE->Fill(v_->E());
            cell[_row][_col][_nstb-1]->p_FmsE->Fill(pout->TotalphotE);
          }; 
        };

      };
    };
    pout->ClearScratch();
    pout->AllToScratch(false); // make a scratch list of all hard phots
    pout->ClusterwithYPhi(false);// cluster for angular size
    TObjArray* Clust=pout->ClusterScratch(.045);
    TLorentzVector Vscr=pout->SumScratch();
    TIter nxt(Clust);
    TObjArray* ctmp;
    while(ctmp=(TObjArray*) nxt())
    {
      Int_t ClusterIndexOf=Clust->IndexOf(ctmp);


      if(ctmp->GetEntries()==2)
      {

        LVec* vab[2];
        vab[0]=(LVec*) ctmp->First();//first photon
        vab[1]=(LVec*) ctmp->After(vab[0]);//other photon
        TLorentzVector twophot=*(vab[0])+*(vab[1]);
        mass=twophot.Mag();	      

        if(vab[0]->Nstb==0)continue;
        if(vab[1]->Nstb==0)continue;
        if(vab[0]->Nstb!=vab[1]->Nstb)continue; // remove events crossing boundary


        TVector3 va3=pout->PosInDet(vab[0],p_geom,0.);
        TVector3 xyab[2];
        xyab[0]=p_geom->LocalXYZ(2,vab[0]->Nstb,va3,true);
        TVector3 vb3=pout->PosInDet(vab[1],p_geom,0.);
        xyab[1]=p_geom->LocalXYZ(2,vab[1]->Nstb,vb3,true);
        Int_t rowab[2],colab[2];
        rowab[0]=(Int_t) xyab[0].Y();
        rowab[1]=(Int_t) xyab[1].Y();
        colab[0]=(Int_t) xyab[0].X();
        colab[1]=(Int_t) xyab[1].X();
        for(Int_t iz=0;iz<2;iz++)
        {


          Int_t iz2=1-iz;
          Bool_t good=true;
          if(rowab[iz]<0 ||colab[iz]<0)good=false;
          if(rowab[iz]>=qt.ROW_NUM[vab[iz]->Nstb-1] ||
              colab[iz]>=qt.COL_NUM[vab[iz]->Nstb-1])good=false;
          if(rowab[iz2]<0 ||colab[iz2]<0)good=false;
          if(rowab[iz2]>=qt.ROW_NUM[vab[iz2]->Nstb-1] ||
              colab[iz2]>=qt.COL_NUM[vab[iz2]->Nstb-1])good=false;


          if(vab[iz]->Iew==2)
          {
            if(vab[iz]->Nstb<3)
            {
              if(colab[iz]<8 && rowab[iz]>8 && rowab[iz]<25)good=false;
            }
            else
            {
              if(colab[iz]<5 && rowab[iz]>6 && rowab[iz]<17)good=false;
            };
          };


          if(good)
          {

            if(cell[rowab[iz]][colab[iz]][vab[iz]->Nstb-1]==0)
            {


              std::cout<<"------------------????? should not be here \n";
              char nam[40];

              out->cd();
              Cell* pcell=new Cell(2,vab[iz]->Nstb,rowab[iz]+1,
                  colab[iz]+1,RunNum,p_geom,&gain,&gcor);
              sprintf(nam,"Cellr%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);
              pcell->SetName(nam);

              cell[rowab[iz]][colab[iz]][vab[iz]->Nstb-1]=pcell;

              mcells.Add(pcell);
              sprintf(nam,"Mcellr%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);
              TH1F* pmcell=new TH1F(nam,nam,240,0.,2.4);
              pcell->p_Mcell=pmcell;
              sprintf(nam,"Exyr%d_c%d_%d",rowab[iz],colab[iz],
                  vab[iz]->Nstb-1);
              TH1F* pexy=new TH1F(nam,nam,120,0.,120.);
              pcell->p_EphXY=pexy;
              sprintf(nam,"MvsZr%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);		      
              TH2F* pmvz=new TH2F(nam,nam,120,0,1.2,10,0,1.);
              pcell->p_MvsZ=pmvz;
              sprintf(nam,"MvsDr%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);		      
              TH2F* pmvd=new TH2F(nam,nam,120,0,1.2,20,0,20*3.82);
              pcell->p_MvsD=pmvd;

              sprintf(nam,"MvsEr%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);		      
              TH2F* pmve=new TH2F(nam,nam,120,0,1.2,20,0,100);
              pcell->p_MvsE=pmve;

              sprintf(nam,"MvsY%d_c%d_%d",rowab[iz],
                  colab[iz],vab[iz]->Nstb-1);		      		      
              TH2F* pmvy=new TH2F(nam,nam,120,0,1.2,30,2.,5.);
              pcell->p_MvsY=pmvy;

            };
            Bool_t FillIt=false;

            Cell* ocl=cell[rowab[iz2]][colab[iz2]][vab[iz2]->Nstb-1];
            Cell* tcl=cell[rowab[iz]][colab[iz]][vab[iz]->Nstb-1];

            if(ocl!=0 && tcl!=0)
            {

              if((ocl->Masspeakfraction>.15 &&
                    ocl->Massmaxcontents>10) || ProcessAll)
              {
                FillIt=true;
              };
            };

            if(tcl!=0 && ProcessAll)FillIt=true;
            if(FillIt)
            {

              Int_t oiz=1-iz;
              Float_t eto=vab[iz]->E()+vab[oiz]->E();
              Float_t y2=(*vab[iz]+ *vab[oiz]).PseudoRapidity();
              Float_t myz=0;
              if(eto>0)myz=fabs(vab[iz]->E()-vab[oiz]->E())/eto;
              if(ocl){
                if(tcl->Geom_ok && ocl->Geom_ok)
                {
                  Float_t dx=ocl->xLocal-tcl->xLocal;
                  Float_t dy=ocl->yLocal-tcl->yLocal;
                  Float_t dxy=sqrt(dx*dx+dy*dy);

                  if(tcl->Instb>2)
                  {
                    if(dxy<8 && dxy>2)tcl->p_Mcell->Fill(mass);
                  }
                  else
                  {
                    if(dxy<28 && dxy>15)tcl->p_Mcell->Fill(mass);			    
                  };
                  tpCellDat cd;
                  cd.Mcell=mass;
                  cd.Epair=eto;
                  cd.Ypair=y2;
                  cd.E1=vab[iz]->E();
                  cd.E2=vab[oiz]->E();
                  cd.X1=xyab[iz].X();
                  cd.X2=xyab[oiz].X();
                  cd.Y1=xyab[iz].Y();
                  cd.Y2=xyab[oiz].Y();
                  cd.NSTB=vab[iz]->Nstb;
                  cd.ievt=pout->ievt;
                  cd.Rnum=pout->Rnum;

                  cd.nSavedHits=pout->nSavedHits;

                  if(cd.nSavedHits>100)cd.nSavedHits=MaxSavedHits;
                  for(int in=0;in<cd.nSavedHits;in++)
                  {
                    cd.SavedHits[in]=pout->SavedHits[in];
                  };
                  tcl->FillTree(cd);

                  tcl->p_EphXY->Fill(vab[iz]->E());
                  tcl->p_MvsZ->Fill(mass,myz);
                  tcl->p_MvsE->Fill(mass,vab[iz]->E());
                  tcl->p_MvsY->Fill(mass,y2);
                  tcl->p_MvsD->Fill(mass,sqrt(dx*dx+dy*dy));

                  Float_t twidth=*(tcl->p_Geom->FpdTowWid(tcl->Iew,tcl->Instb));

                  if(eto> 10 && twidth>0 && mass<.8)tcl->p_dxdyM->Fill(ocl->xLocal/twidth-.5,ocl->yLocal/twidth-.5,mass);
                  if(eto> 5 && twidth>0)tcl->p_dxdyE->Fill(ocl->xLocal/twidth-.5,ocl->yLocal/twidth-.5,vab[oiz]->E());
                  Float_t eall=EinNSTB[ocl->Instb];
                  if(eall>5 && twidth>0 && mass<.8)tcl->p_dxdyER->Fill(ocl->xLocal/twidth-.5,ocl->yLocal/twidth-.5,eto/eall);
                };
              };
            };

          };
        };
      };
    };
  };
  std::cout<< " Out of loop\n";
  TIter next2(&mcells);
  Cell* pcl;
  int docnt=0;
  while(pcl=(Cell*) next2())
  {

    Int_t idet=1;
    if(pcl->Iew==2&& ((idet=pcl->Instb)<5))
    {
      Int_t chan=(pcl->Row1-1)*qt.COL_NUM[idet-1] + pcl->Col1;
      if(pcl->p_adc)delete pcl->p_adc;
      if(pcl->p_adcLed)delete pcl->p_adcLed;
      char nam[40];
      sprintf(nam,"adc_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
      printf("set adc %s \n",nam);
      pcl->p_adc=qtHist[idet-1]->ProjectionY(nam,chan,chan);
      sprintf(nam,"adcLed_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
      printf("set adcLed %s \n",nam);
      pcl->p_adcLed=qtHistLed[idet-1]->ProjectionY(nam,chan,chan);
      if(TestCellTree(pcl))
      {
        printf("%s ok entries =%d \n",(const char*) pcl->GetName(),(int) pcl->CellD->GetEntries());
        TString fncell="CodeSmall/data2/_";
        fncell=fncell+pcl->GetName()+".root";
        cout<<"opening "<<fncell<<"\n";
        TFile* savecd=gROOT->GetFile();
        TFile* Fc=new TFile(fncell,"recreate");
        cout<<"open\n";
        TTree* trn=pcl->CellD->CopyTree("E1>0");
        cout<<"trn created with entries="<<trn->GetEntries()<<"linked to dir="<<trn->GetDirectory()->GetName()<<"\n";
        trn->Write("trn");
        TString cellname=pcl->GetName();
        TString treen=pcl->CellD->GetName();
        const char* ctrenm=(const char*) treen;
        int jmax=pcl->CellD->GetEntries();
        Cell* ncl=new Cell();
        *ncl=*pcl;
        ncl->CellD->SetDirectory(gDirectory);
        ncl->CellD=trn->CopyTree("E1>0");
        cout<<"ncl created with entries="<<ncl->CellD->GetEntries()<<"linked to dir="<<ncl->CellD->GetDirectory()->GetName()<<"\n";
        pcl->CellD=0;
        pcl->Write(cellname);	      
        Fc->Close();
        delete Fc;
        cout<<fncell<<" wrote \n";
        if(savecd)
        {
          savecd->cd();
        };
      }
      else
      {
        if(pcl->CellD)
        {
          printf("%s bad entries=%d\n",(const char*) pcl->GetName(),(int)pcl->CellD->GetEntries()); 
        } else { printf("%s bad no CellD\n",(const char*) pcl->GetName());}; 
      };
    };

  };
  std::cout<<"Finish adc setup\n";
  out->cd();
  for(Int_t j=0;j<4;j++)
  {
    qtHist[j]->Write();
    printf("qtHist[%d] Dimensions %d by %d \n",j,
        qtHist[j]->GetNbinsX(),qtHist[j]->GetNbinsY());
  };
  //  PrintCell(cell[3][3][1],20,2,"testCellr3_c3_1.ps");

  mcells.Print();

  mcells.Write("Mycells",kSingleKey);

  Cell* cm=(Cell*) mcells.FindObject("Cellr3_c3_1");
  out->Close();
  out=new TFile(OutFileName,"read");
  TObjArray* mcells2=(TObjArray*) out->Get("Mycells");
  Cell* cc=0;
  TIter cnext(mcells2);
  while(Cell* cc=(Cell*) cnext())
  {
    TString cnm=cc->GetName();
    if(TestCellTree(cc)){printf("check %s ok\n",(const char*)cnm);}
    else {printf("check %s bad\n",(const char*)cnm);};
  }
  // Cell* cm2=(Cell*) mcells2->FindObject("Cellr3_c3_1");
  //PrintCell(cm2,20,2,"testCellr3_c3_1c.ps");

};


Int_t AnalTools::readq2(char* filelist,Int_t set,FilesSet* p_files)
{

  FitTower::SetSubclu2(true);//Yiqun will substitute clu2 for clusters
  printf("readq2 called\n");
  Geom* p_geom=new Geom(p_files);
  CalibStr* pgain=new CalibStr(10170000,p_files->p_fpdgain()->Path());
  CalibStr* pgaincorr=new CalibStr(10170000,p_files->p_fpdgaincorr()->Path());

  pout=new poutTree(filelist);
  typedef struct {
    Int_t NRealCluster1;
    Int_t NRealCluster2;
    HitCluster cl1[2];
    HitCluster cl2[2];
    Int_t NPh;
    Float_t M12;
    Float_t SigmaMax0[5];
    Float_t EDet;
    Float_t m2_3;
    Float_t mm2_3;
    Int_t HiTowNstb;
    Int_t HiTowrow;
    Int_t HiTowcol;

    Int_t Hi2ndTowrow;
    Int_t Hi2ndTowcol;

    Float_t HiPairFrac;
    Float_t HiTowFrac;
    Float_t Near4Frac;
    Float_t Diag4Frac;
    Float_t predictHiFrac;
    Float_t predictNearFrac;
    Float_t predictDiagFrac;
    Float_t predictHiPairFrac;
    Float_t LeftTowSum;
    Float_t RightTowSum;
    Float_t UpTowSum;
    Float_t DownTowSum;
    Float_t m2gl;
    Float_t e1gl;
    Float_t e2gl;
    Float_t pephot;
    Float_t Hi2Sep;
    Float_t LHi2Sep;
    Float_t L2ndHiFrac;
    Float_t L2ndHirow;
    Float_t L2ndHicol;
  } tp_0;

  tp_0 tr0;
  TTree rab("rab","rab");

  rab.Branch("p_LeftTowSum",&tr0.LeftTowSum,"LeftTowSum/F");
  rab.Branch("p_RightTowSum",&tr0.RightTowSum,"RightTowSum/F");
  rab.Branch("p_UpTowSum",&tr0.UpTowSum,"UpTowSum/F");
  rab.Branch("p_DownTowSum",&tr0.DownTowSum,"DownTowSum/F");

  rab.Branch("p_HiTowNstb",&tr0.HiTowNstb,"HiTowNstb/I");
  rab.Branch("p_HiTowrow",&tr0.HiTowrow,"HiTowrow/I");
  rab.Branch("p_HiTowcol",&tr0.HiTowcol,"HiTowcol/I");
  rab.Branch("p_Hi2ndTowrow",&tr0.Hi2ndTowrow,"Hi2ndTowrow/I");
  rab.Branch("p_Hi2ndTowcol",&tr0.Hi2ndTowcol,"Hi2ndTowcol/I");

  rab.Branch("p_HiTowFrac",&tr0.HiTowFrac,"HiTowFrac/F");
  rab.Branch("p_Near4Frac",&tr0.Near4Frac,"Near4Frac/F");
  rab.Branch("p_Diag4Frac",&tr0.Diag4Frac,"Diag4Frac/F");
  rab.Branch("p_HiPairFrac",&tr0.HiPairFrac,"HiPairFrac/F");

  rab.Branch("p_pHiTowFrac",&tr0.predictHiFrac,"pHiTowFrac/F");
  rab.Branch("p_pHiPairFrac",&tr0.predictHiPairFrac,"pHiPairFrac/F");
  rab.Branch("p_pNearFrac",&tr0.predictNearFrac,"pNearFrac/F");
  rab.Branch("p_pDiagFrac",&tr0.predictDiagFrac,"pDiagFrac/F");

  rab.Branch("p_NReal1",&tr0.NRealCluster1,"NR1/I");
  rab.Branch("p_NReal2",&tr0.NRealCluster2,"NR2/I");
  rab.Branch("p_N1phc1",&tr0.cl1[0].nPhoton,"N1phc1/I");
  rab.Branch("p_N1phc2",&tr0.cl1[1].nPhoton,"N1phc2/I");
  rab.Branch("p_N2phc1",&tr0.cl2[0].nPhoton,"N2phc1/I");
  rab.Branch("p_N2phc2",&tr0.cl2[1].nPhoton,"N2phc2/I");

  rab.Branch("p_N1towc1",&tr0.cl1[0].numbTower,"N1towc1/I");
  rab.Branch("p_N1towc2",&tr0.cl1[1].numbTower,"N1towc2/I");
  rab.Branch("p_N2towc1",&tr0.cl2[0].numbTower,"N2towc1/I");
  rab.Branch("p_N2towc2",&tr0.cl2[1].numbTower,"N2towc2/I");

  rab.Branch("p_1clu1E",&tr0.cl1[0].energy,"e1c1/F");
  rab.Branch("p_1clu2E",&tr0.cl1[1].energy,"e1c2/F");
  rab.Branch("p_2clu1E",&tr0.cl2[0].energy,"e2c1/F");
  rab.Branch("p_2clu2E",&tr0.cl2[1].energy,"e2c2/F");

  rab.Branch("p_1chic1",&tr0.cl1[0].chiSquare,"chi1c1/F");
  rab.Branch("p_1chic2",&tr0.cl1[1].chiSquare,"chi1c2/F");
  rab.Branch("p_2chic1",&tr0.cl2[0].chiSquare,"chi2c1/F");
  rab.Branch("p_2chic2",&tr0.cl2[1].chiSquare,"chi2c2/F");

  rab.Branch("p_1x1_1",&tr0.cl1[0].photon[0].xPos,"x11_1/F");
  rab.Branch("p_1x1_2",&tr0.cl1[0].photon[1].xPos,"x11_2/F");
  rab.Branch("p_1x2_1",&tr0.cl1[1].photon[0].xPos,"x12_1/F");
  rab.Branch("p_1x2_2",&tr0.cl1[1].photon[1].xPos,"x12_2/F");

  rab.Branch("p_1y1_1",&tr0.cl1[0].photon[0].yPos,"y11_1/F");
  rab.Branch("p_1y1_2",&tr0.cl1[0].photon[1].yPos,"y11_2/F");
  rab.Branch("p_1y2_1",&tr0.cl1[1].photon[0].yPos,"y12_1/F");
  rab.Branch("p_1y2_2",&tr0.cl1[1].photon[1].yPos,"y12_2/F");

  rab.Branch("p_1e1_1",&tr0.cl1[0].photon[0].energy,"e11_1/F");
  rab.Branch("p_1e1_2",&tr0.cl1[0].photon[1].energy,"e11_2/F");
  rab.Branch("p_1e2_1",&tr0.cl1[1].photon[0].energy,"e12_1/F");
  rab.Branch("p_1e2_2",&tr0.cl1[1].photon[1].energy,"e12_2/F");

  rab.Branch("p_2x1_1",&tr0.cl2[0].photon[0].xPos,"x21_1/F");
  rab.Branch("p_2x1_2",&tr0.cl2[0].photon[1].xPos,"x21_2/F");
  rab.Branch("p_2x2_1",&tr0.cl2[1].photon[0].xPos,"x22_1/F");
  rab.Branch("p_2x2_2",&tr0.cl2[1].photon[1].xPos,"x22_2/F");

  rab.Branch("p_2y1_1",&tr0.cl2[0].photon[0].yPos,"y21_1/F");
  rab.Branch("p_2y1_2",&tr0.cl2[0].photon[1].yPos,"y21_2/F");
  rab.Branch("p_2y2_1",&tr0.cl2[1].photon[0].yPos,"y22_1/F");
  rab.Branch("p_2y2_2",&tr0.cl2[1].photon[1].yPos,"y22_2/F");

  rab.Branch("p_2e1_1",&tr0.cl2[0].photon[0].energy,"e21_1/F");
  rab.Branch("p_2e1_2",&tr0.cl2[0].photon[1].energy,"e21_2/F");
  rab.Branch("p_2e2_1",&tr0.cl2[1].photon[0].energy,"e22_1/F");
  rab.Branch("p_2e2_2",&tr0.cl2[1].photon[1].energy,"e22_2/F");
  rab.Branch("p_SigmaMax0",tr0.SigmaMax0,"SigmaMax0[5]/F");
  rab.Branch("p_NPh",&tr0.NPh,"NPh/I");
  rab.Branch("p_M12",&tr0.M12,"M12/F");
  rab.Branch("p_EDet",&tr0.EDet,"EDet/F");
  rab.Branch("p_m2_3",&tr0.m2_3,"m2_3/F");
  rab.Branch("p_mm2_3",&tr0.mm2_3,"mm2_3/F");
  rab.Branch("p_m2gl",&tr0.m2gl,"m2gl/F");
  rab.Branch("p_e1gl",&tr0.e1gl,"e1gl/F");
  rab.Branch("p_e2gl",&tr0.e2gl,"e2gl/F");
  rab.Branch("p_pephot",&tr0.pephot,"pephot/F");
  rab.Branch("p_Hi2Sep",&tr0.Hi2Sep,"Hi2Sep/F");
  rab.Branch("p_LHi2Sep",&tr0.LHi2Sep,"LHi2Sep/F");
  rab.Branch("p_L2ndHiFrac",&tr0.L2ndHiFrac,"L2ndHiFrac/F");
  rab.Branch("p_L2ndHirow",&tr0.L2ndHirow,"L2ndHirow");
  rab.Branch("p_L2ndHicol",&tr0.L2ndHicol,"L2ndHicol");

  for(Int_t kk=1;kk<9;kk++)
  {
    pout->MinEnergy[kk]=6.;
    //Large Cell threshold
    if(kk==5 || kk==6)pout->MinEnergy[kk]=4.;
  };  

  TFile Out("Outputq2.root","recreate");
  int cnt=0;
  Int_t nentries=pout->nentries;
  printf("Events to Process=%d \n",nentries);
  if(NumberEventsToProcess>pout->nentries)nentries=pout->nentries;
  else if (NumberEventsToProcess>0)nentries=NumberEventsToProcess;
  TCanvas* cfms=new TCanvas("cfms","cfms",600,600);
  cfms->Divide(2,2);
  TH2F* matfms[4];
  for(int jj=0;jj<4;jj++)matfms[jj]=0;
  cfms->Print("cfms.ps(");
  cfms->GetPad(4)->SetLogz();
  TH1F* chi2 = new TH1F("chi2","chi2",500,0,50);
  TH2F* Chi2D=new TH2F("Chi2D","Chi2D",1000,0,400,1000,0,400);
  TH2F* devVSr=new TH2F("devVSr","devVSr",50,0.,10.,200,-.1,.1);
  Chi2D->GetXaxis()->SetRangeUser(0,50);
  TH3F* Rdev=new TH3F("Rdev","Rdev",40,0.,10,40,0,10.,1000,0,1.);
  TH2F* hdev2=0;
  TH2F* hhdev2=0;
  TH2F* hfit=0;
  TH2F* hhfit=0;
  TGraph* hbr=0;
  TGraph* hbr2=0;
  Int_t hbrcnt=0;
  Int_t hbr2cnt=0.;
  TH2F* mhist=0;
  TH2D* rproj=0;
  TH2F* posxy=new TH2F("posxy","posxy",120,0,12.,240,0,24);
  TH2F* unitxy=new TH2F("unitxy","unitxy",30,-1,2.,30,-1,2);
  for(Int_t iev=0;iev<nentries;iev++)
  {
    //      printf("------------iev=%d\n",iev);
    Int_t selectk=1;
    tr0.e1gl=0;
    tr0.e2gl=0;
    tr0.m2gl=-1.;
    tr0.Hi2Sep=0;
    tr0.LHi2Sep=0;
    tr0.L2ndHiFrac=0;
    tr0.L2ndHirow-1;
    tr0.L2ndHicol-1;

    tr0.LeftTowSum=0;
    tr0.RightTowSum=0;
    tr0.UpTowSum=0;
    tr0.DownTowSum=0;

    tr0.HiTowNstb=-1;
    tr0.HiTowrow=-1;
    tr0.HiTowcol=-1;
    tr0.Hi2ndTowrow=-1;
    tr0.Hi2ndTowcol=-1;
    tr0.HiTowFrac=0.;
    tr0.HiPairFrac=0.;
    tr0.Near4Frac=0.;
    tr0.Diag4Frac=0.;

    tr0.predictHiFrac=0.;
    tr0.predictHiPairFrac=0.;
    tr0.predictNearFrac=0.;
    tr0.predictDiagFrac=0.;
    for(int j=0;j<5;j++)tr0.SigmaMax0[j]=0;
    tr0.m2_3=0.;
    tr0.mm2_3=0.;

    if( (iev%25000) ==0)printf( "%d Events Read \n ",iev);
    pout->GetEntry(iev);
    Float_t Ed1[4],En1[4],En2[4],En[4];
    Int_t Nphh[4];
    Float_t m2[4]; 
    for(Int_t j0=0;j0<4;j0++)
    {
      Ed1[j0]=En1[j0]=En2[j0]=En[j0]=0.;
      Nphh[j0]=0;
      m2[j0]=0;
    };
    TIter diter(pout->dlist);
    while(LVec* v=(LVec*) diter())
    {
      if(v->Nstb>=0 && v->Nstb<5 )
      {
        Ed1[v->Nstb-1]=v->E();
      };
    };

    bool veto=false;
    TIter viter(pout->vlist);

    pout->ClearScratch();
    while(LVec* v=(LVec*) viter())
    {
      Int_t j0;
      if((j0=v->Nstb-1)>=0 && v->Nstb<5)
      { 
        TVector3 posindet=pout->PosInDet(v,p_geom);
        En[j0]+=v->E();
        Nphh[j0]++;
        if(pout->NearEdge(v,p_geom,1.25)!=0)veto=true;
        if(j0==selectk)pout->scratchlist->Add(v);
        char* str="true";
        if(!veto)str="false";
        //	      printf("Accum En[%d],Nphh[%d]=%f %d veto=%s  \n",j0,j0,En[j0],Nphh[j0],str);
        if( j0==selectk && Nphh[j0]==2 && !veto)
        {
          m2[selectk]=pout->SumScratch().Mag();
        };
      };
    }; 
    pout->ClearScratch();
    tr0.m2_3=m2[selectk];
    //     printf("En[%d]=%f\n",selectk,En[selectk]);
    bool cuta1=(En[selectk]-Energy_study)>-25 && En[selectk]<130 && !veto && (Nphh[selectk]>0);
    bool cuta1b=En[selectk]>2 && En[selectk]<50 && !veto && (Nphh[selectk]>0);
    bool cuta3=fabs(En[selectk]-Energy_study) <15 && !veto && (Nphh[selectk]==2)&&fabs(m2[selectk]-.135)<.12;
    //look for etas
    bool cuta2=fabs(En[selectk]-Energy_study-10)<15. && !veto&& (Nphh[selectk]==2 &&fabs(m2[selectk]-.55)<.3);
    //if(fabs(En[3]-Energy_study-15)<10. && !veto && (Nphh[3]==2))
    //if(fabs(En[3]-Energy_study-10)<15. && !veto&& (Nphh[3]==2 &&fabs(m2[3]-.55)<.3) )

    if(cuta1b )
    {
      bool seeall=false;
      if(seeall)
      {
        for(int k=0;k<4;k++)
        {
          if(matfms[k])delete matfms[k];
          char snam[100];
          sprintf(snam,"matfms%d",k);
          TMatrix m=pout->FillFMSADC(k+1);
          for(int row=0;row<m.GetNrows();row++)
          {
            for(int col=0;col<m.GetNcols();col++)
            {
              (m[row][col]) *=pgain->GetValue(2,k+1,row,col);
              (m[row][col]) *=pgaincorr->GetValue(2,k+1,row,col);
            };
          };
          matfms[k]=new TH2F(m);
          matfms[k]->SetName(snam);
          cfms->cd(k+1);
          printf("cd(%d), Esum=%f\n",k+1,matfms[k]->Integral());
          matfms[k]->Draw("zcoltext");
        };
        cfms->Update();
        int sretur=system("sleep 2");
        printf("update matfms\n");
      };
      for(int k=selectk;k<selectk+1;k++)
      {
        cfms->cd(1);
        cfms->GetPad(1)->SetLogz();
        TMatrix m=pout->FillFMSADC(k+1);
        for(int row=0;row<m.GetNrows();row++)
        {
          for(int col=0;col<m.GetNcols();col++)
          {
            (m[row][col]) *=pgain->GetValue(2,k+1,row,col);
            (m[row][col]) *=pgaincorr->GetValue(2,k+1,row,col);
          };
        };
        tr0.EDet=m.Sum();
        FitTower::SetNoCatag(true,1);
        Yiqun rec(&m,p_geom,pgain,pgaincorr,2,k+1);
        //	      rec.PrintClu();
        FitTower::SetNoCatag(false);	      
        int nc=rec.NRealClusts;
        cfms->cd(3);
        if(nc>2)nc=2;
        tr0.NRealCluster1=rec.NRealClusts;
        //	      printf("tr0.NRealCluster1=%d NPh=%d \n",nc,rec.NPh);
        for(int jc=0;jc<nc;jc++)
        {
          tr0.cl1[jc]=rec.clust[jc];
        };
        tr0.NPh=rec.NPh;
        tr0.M12=0;
        if(rec.NPh>1)
        {
          tr0.M12=(rec.mom(0)+rec.mom(1)).Mag();
          //printf("NPh=%d tr0.M12=%f\n",rec.NPh,tr0.M12);
        };
        Int_t clust50;
        if(hbr)delete hbr;
        hbr=new TGraph();
        hbrcnt=0;
        Float_t E1,E2,x1,x2,y1,y2;
        if(rec.NRealClusts>0)
        {
          E1=rec.clust[0].energy;
          x1=rec.clust[0].x0;
          y1=rec.clust[0].y0;
        }

        if(rec.NRealClusts==2)
        {
          E2=rec.clust[1].energy;
          x2=rec.clust[1].x0;
          y2=rec.clust[1].y0;
        }

        bool cutb1=rec.NPh==2 &&rec.NRealClusts==2 &&fabs(m2[selectk]-.135)<.2 ;
        bool cutb3=fabs(rec.NRealClusts-2)<1.5 &&E1>5 ;
        bool cutb2=rec.NPh==2 &&rec.NRealClusts==2 &&fabs(m2[selectk]-.55)<.3;
        if(cutb3 && m.Sum()<130)
        {
          //		  printf(" . ");
          float ephot=rec.mom(0).E();
          if(hbr!=0)delete hbr;
          hbr=new TGraph();
          hbrcnt=0;
          if(rec.NPh==1)
          {
            hbr->SetPoint(hbrcnt++,
                rec.photons[0].xPos/rec.widLG[0],
                rec.photons[0].yPos/rec.widLG[1]);
            clust50=0;
          };
          if(rec.NPh==2)
          {
            if(rec.NRealClusts==2)hbr->SetPoint(hbrcnt++,
                rec.clust[0].photon[0].xPos/rec.widLG[0],
                rec.clust[0].photon[0].yPos/rec.widLG[1]);
            clust50=0;
            if(rec.clust[1].energy>rec.clust[0].energy)clust50=1;
          };

          //	      chi2->Fill(rec.ChiSqG);
          if(clust50>=0)
          {
            if(rec.clust[clust50].ShadowPointers)
            {
              TIter shnext(rec.clust[clust50].ShadowPointers);
              HitCluster* sclu;
              int shcnt=0;
              while(sclu=(HitCluster*) shnext())
              {
                if(shcnt<5)tr0.SigmaMax0[shcnt]=sclu->sigmaMax;
                shcnt++;
              };
            };
            if(rec.clust[clust50].photon[0].energy>40)
            {
              if(rec.clust[clust50].chiSquare<=.0000001)
              {
                printf("En[%d]=%f\n",selectk,En[selectk]);
                rec.PrintClu();
              };
              chi2->Fill(rec.clust[clust50].chiSquare);
              //			  printf("rec.clust[clust50.chiSquare]=%f\n",rec.clust[clust50].chiSquare);
            };
          };
          ephot=0;
          int Ntow50=0;
          if(clust50>=0)
          {
            ephot=rec.clust[clust50].energy;
            Ntow50=rec.clust[clust50].numbTower;
          };
          if(Ntow50<2)ephot=0;
          float xx=rec.clust[clust50].photon[0].xPos/rec.widLG[0];
          float yy=rec.clust[clust50].photon[0].yPos/rec.widLG[1];
          //		  printf("clust50=%d, ephot=%f ntow=%d\n",clust50,ephot,Ntow50);
          if(ephot>5 
              //		  &&   abs(fmod(xx,1)-.5)<.3 
              //		  &&   abs(fmod(yy,1)-.5)<.3 
            )
          {
            cnt++;
            //		      printf("ephot(cnt=%d)=%f\n",cnt,ephot);
            posxy->Fill(xx,yy);
            unitxy->Fill(fmod(xx,1.),fmod(yy,1.));
            TMatrix nm=(1/ephot)*m;

            TowerFPD* p_tow;
            TowerFPD* tw1=0;
            TowerFPD* tw2=0;
            TowerFPD* ltw2=0;
            if(rec.clust[clust50].numbTower>1)
            {
              TIter towiter(rec.clust[clust50].tow);
              TIter ltowiter(rec.clust[clust50].tow);
              tw1=(TowerFPD*) ltowiter();
              tw1=(TowerFPD*) towiter();
              tw2=(TowerFPD*) towiter();
              tr0.LHi2Sep=0.;
              tr0.Hi2ndTowrow=tw2->row-1;
              tr0.Hi2ndTowcol=tw2->col-1;
              tr0.HiPairFrac=(tw1->energy+tw2->energy)/ephot;
              tr0.HiTowrow=tw1->row-1;
              tr0.HiTowcol=tw1->col-1;
              tr0.HiTowFrac=tw1->energy/ephot;
              bool FindLHi2=false;
              while(ltw2=(TowerFPD*) ltowiter())
              {
                if(ltw2->LocalHighTower() && !FindLHi2)
                {
                  tr0.L2ndHirow=ltw2->row-1;
                  tr0.L2ndHicol=ltw2->col-1;
                  tr0.L2ndHiFrac=ltw2->energy/ephot;
                  FindLHi2=true;
                  tr0.LHi2Sep=sqrt(pow(tr0.HiTowrow-tr0.L2ndHirow,2)
                      +pow(tr0.HiTowcol-tr0.L2ndHicol,2)
                      );
                };
              };
              tr0.Hi2Sep=sqrt(pow(tr0.HiTowrow-tr0.Hi2ndTowrow,2)
                  +pow(tr0.HiTowcol-tr0.Hi2ndTowcol,2)
                  );
            };
            TIter next0(rec.clust[clust50].tow);
            Float_t pephot=0;
            while(p_tow=(TowerFPD*) next0())
            {
              int row=p_tow->row-1;
              int col=p_tow->col-1;
              float ecell=m[row][col];
              float mydev=(*(rec.pwe->dev))[row][col];
              float mypredict=ecell-mydev;
              pephot+=mypredict;
              //			  printf("r,c,ec,myp=%d,%d,%f,%f\n",row,col,ecell,mypredict);
            };
            tr0.pephot=pephot;

            //		      printf("pephot=%f \n",pephot);
            TIter next(rec.clust[clust50].tow);

            while(p_tow=(TowerFPD*) next())
            {
              int row=p_tow->row-1;
              int col=p_tow->col-1;
              if(col< tw1->col-1)tr0.LeftTowSum+=p_tow->energy/ephot;
              if(col> tw1->col-1)tr0.RightTowSum+=p_tow->energy/ephot;
              if(row< tw1->row-1)tr0.DownTowSum+=p_tow->energy/ephot;
              if(row> tw1->row-1)tr0.UpTowSum+=p_tow->energy/ephot;
              float xcell=(col+.5)*rec.widLG[0]-rec.clust[clust50].photon[0].xPos;
              float ycell=(row+.5)*rec.widLG[1]-rec.clust[clust50].photon[0].yPos;
              float r=sqrt(xcell*xcell+ycell*ycell);
              float ecell=m[row][col];
              float rcell=sqrt(xcell*xcell+ycell*ycell);
              float mydev=(*(rec.pwe->dev))[row][col];
              if(rec.clust[clust50].chiSquare<500)
              {
                devVSr->Fill(r,mydev/ephot);
                float mypredict=ecell-mydev;
                if(mypredict!=0)
                {
                  tr0.HiTowNstb=4;
                  int rowh=tw1->row-1;
                  int colh=tw1->col-1;
                  if(colh==col && rowh==row)
                  {
                    tr0.predictHiFrac=mypredict/pephot;
                  };
                  if(p_tow==tw2 ||p_tow==tw1)
                  {
                    tr0.predictHiPairFrac+=mypredict/pephot;
                  };
                  if(abs(colh-col)==1 && abs(rowh-row)==1)
                  {
                    tr0.Diag4Frac+=ecell/ephot; 
                    tr0.predictDiagFrac+=mypredict/pephot;
                  }
                  else if (abs(colh-col)==1 || abs(rowh-row)==1)
                  {
                    tr0.Near4Frac+=ecell/ephot;
                    tr0.predictNearFrac+=mypredict/pephot;
                  };
                  Float_t e11_1=tr0.cl1[0].photon[0].energy;
                  //Sigmamax cut .55 small cells
                  //Sigmamax cut .75 large cells
                  float_t Sigmacut=.55;
                  if(k<2)Sigmacut=.8;
                  //single photon cut for large cells
                  Float_t e1c1=rec.clust[clust50].energy;
                  if(fabs(e1c1-10)<5 && tr0.LHi2Sep<1.2 &&tr0.mm2_3<.1 ) 
                  {
                    Rdev->Fill(xcell,ycell,
                        ecell/ephot*pephot/e11_1);
                  };
                };
              };

            };

            bool do_update=false;
            if(cnt%1000==0 )
              //if(tr0.HiTowFrac<.8)
            {		
              do_update=true;
              cfms->cd(1);
              pout->DrawFMSADC(2,k+1,p_geom);
              if(mhist)delete mhist;
              mhist=new TH2F(m);
              mhist->SetName("m");
              mhist->SetStats(0);
              char Tstr[100];
              sprintf(Tstr,"E=%f chi1=%f chi2=%f SigmaMax[2]=%f HiTowFrac=%f \n",
                  tr0.cl1[0].energy, tr0.cl1[0].chiSquare,
                  tr0.cl2[0].chiSquare,tr0.SigmaMax0[2],
                  tr0.HiTowFrac);
              mhist->SetTitle(Tstr);
              mhist->Draw("zcol");
              mhist->Draw("sametext");
              cfms->cd(2);
              cfms->GetPad(2)->SetLogz();

              char titl[100];
              sprintf(titl,"chi2=%f4.2",rec.ChiSqG);
              if(hdev2!=0)delete hdev2;
              if(hfit!=0)delete hfit;
              hfit=new TH2F(rec.FittedMat());
              hdev2=new TH2F(*(rec.pwe->dchi2));
              hdev2->SetTitle(titl);
              hdev2->SetMinimum(-5);
              hdev2->SetMaximum(5);
              //			  hdev2->Draw("box");
              hfit->SetStats(0);
              hfit->SetTitle(titl);
              hfit->Draw("zcol");
              hfit->Draw("textsame");
              //			  rec.pwe->dev->Draw("zcol");
              //			  rec.pwe->dev->Draw("boxsame");
              hbr->Draw("*");
              cfms->cd(3);
              chi2->Draw();
              cfms->cd(4);
              cfms->GetPad(4)->SetLogz();
              if(rproj!=0)delete rproj;
              rproj=(TH2D*) Rdev->Project3DProfile("yx");
              rproj->Draw("zcol");

            };
            FitTower::SetNoCatag(true,2);
            FitTower::SetForceMass(-1);
            FitTower::SetDoGlobal(false);
            Yiqun rec2(&m,p_geom,pgain,pgaincorr,2,k+1);
            FitTower::SetForceMass(-1.);
            FitTower::SetNoCatag(false);
            FitTower::SetDoGlobal(true);

            if(rec2.NPh<2)rec2.PrintClu();
            if(rec2.NPh>1 && do_update)
            {
              if(hhdev2!=0)delete hhdev2;
              if(hhfit!=0)delete hhfit;
              hhdev2=new TH2F(*(rec2.pwe->dchi2));
              hhfit=new TH2F(rec2.FittedMat());
              cfms->cd(3);
              hhdev2->SetLineColor(4);
              hhdev2->SetMaximum(5);
              hhdev2->SetMinimum(-5);
              hhdev2->Draw("boxtext");
              cfms->cd(4);
              hhfit->Draw("zcol");
              hhfit->Draw("sametext");
              //			  printf("NPh=%d tr0.mm2=%f \n",rec2.NPh,tr0.mm2_3);
              //			  printf("rec2.clust[0/1].chiSquare=%f %f\n",rec2.clust[0].chiSquare,
              //rec2.clust[1].chiSquare);

            };
            if(rec2.NPh>1)tr0.mm2_3=(rec2.mom(0)+rec2.mom(1)).Mag();
            int nc2=rec.NRealClusts;
            if(nc2>2)nc2=2;
            tr0.NRealCluster2=rec2.NRealClusts;
            for(int jc=0;jc<nc2;jc++)
            {
              tr0.cl2[jc]=rec2.clust[jc];
            };

            if(hbr2)delete hbr2;
            hbr2cnt=0;
            hbr2=new TGraph();
            hbr2->SetMarkerColor(4);
            //		      if( rec2.NPh>=2  && rec2.ChiSqG>.1 && rec2.NRealClusts>=1 )
            if( rec2.NPh>=2  && rec2.NRealClusts>=1 )
            {
              //			  printf("NPh=%d m12=%f ntowers=%d\n", rec2.NPh,(rec2.mom(0)+rec2.mom(1)).Mag(),rec2.clust[0].numbTower);
              for(int ncl=0;ncl<rec2.clust[0].nPhoton;ncl++)
              {
                float x,y;
                x=rec2.clust[0].photon[ncl].xPos;
                y=rec2.clust[0].photon[ncl].yPos;
                hbr2->SetPoint(hbr2cnt++,x/rec2.widLG[0],y/rec2.widLG[1]);
                //			      printf("cl_ph#=%d en=%f x0=%f y0=%f  \n",ncl,
                //	     rec2.clust[0].photon[ncl].energy,x,y);

              };
              //			  printf("chi_1=%f chi_2 %f \n",rec.clust[clust50].chiSquare,rec2.clust[0].chiSquare);
              cfms->cd(3);
              hbr2->Draw("*");
              float ea,eb;
              ea=rec2.clust[0].photon[0].energy;
              eb=rec2.clust[0].photon[1].energy;
              /*
                 cfms->cd(3);
                 Chi2D->GetXaxis()->SetRangeUser(0,50);
                 Chi2D->GetYaxis()->SetRangeUser(0,50);
                 Chi2D->Draw("box");
                 */
              FitTower::SetForceMass(-1.);
              FitTower::SetNoCatag(false);
              FitTower::SetDoGlobal(true);
              Yiqun rec3(&m,p_geom,pgain,pgaincorr,2,k+1);
              if(rec3.NPh>1)
              {
                tr0.e1gl=rec3.mom(0).E();
                tr0.e2gl=rec3.mom(1).E();
                tr0.m2gl=(rec3.mom(0)+rec3.mom(1)).Mag();
              };
              if(rec2.clust[0].numbTower>6 )
              {
                Chi2D->Fill(rec.clust[clust50].chiSquare,rec2.clust[0].chiSquare);
              };
              rab.Fill();

            };
            if(do_update)
            {
              cfms->Print("cfms.ps");
              cfms->Update();
              //			  system("sleep 3");
              //			  printf("end sleep\n");
            };

          };
        };
      };

    };
    if(iev%9000==0)
    {
      Rdev->Write("Rdev");
      devVSr->Write("devVSr");
      TGraph2DErrors gr;
      int grcnt=0;
      TH2D* gyx=(TH2D*) Rdev->Project3DProfile("yx");
      int nx=gyx->GetNbinsX();
      int ny=gyx->GetNbinsY();
      TGraph sh;
      int shcnt=0;
      for(int j=1;j<nx;j++)
      {
        for(int k=1;k<ny;k++)
        {
          float w=.25;
          float x=w*(1.*j+.5);
          float y=w*(1.*k+.5);
          float r=sqrt(x*x+y*y);
          int bin=gyx->FindBin(x,y);
          float val=gyx->GetBinContent(bin);
          if( r<10)
          {
            sh.SetPoint(shcnt++,r,val);
            gr.SetPoint(grcnt++,x,y,val);
            //	      gr.SetPoint(grcnt++,a-x,y,val);
            //	      gr.SetPoint(grcnt++,-x,-y,val);
            //	      gr.SetPoint(grcnt++,x,-y,val);

          };
        };
      };
      gyx->Write("gyx");
      gr.Write("gr");
      sh.Write("sh");
      posxy->Write("posxy");
      unitxy->Write("unitxy");
      chi2->Write("chi2");
      Chi2D->Write("chi2D");
      rab.Write("rab");
      Out.Purge();
    };
  }
};
Int_t AnalTools::readq1(char* filelist,Int_t set,FilesSet* p_files)
{
  Fill fill(p_files,18000,20000,16050001,16170001);//get fill patterns for run15
  Geom* p_geom=new Geom(p_files);
  CalibStr* pgain=new CalibStr(10170000,p_files->p_fpdgain()->Path());
  CalibStr* pgaincorr=new CalibStr(10170000,p_files->p_fpdgaincorr()->Path());

  pout=new poutTree(filelist);
  pout->EnableEdepCorr=0;//not Enable Edep corrections
  if(OutputToSingle) { 
    std::cout<<"Output to "<<OutputFileName<<"\n";
    TFile* NewOut=new TFile("OutputFileName","recreate");
    TTree* newp_out= pout->p_out->CloneTree();
    newp_out->Write();
    return 0;
  };
  std::cout<<"OuputToSingle not selected\n";

  //set energy threshold to include photons
  /*  for(Int_t kk=1;kk<9;kk++) 
      {
      pout->MinEnergy[kk]=6.;
  //Large Cell threshold
  if(kk==5 || kk==6)pout->MinEnergy[kk]=4.;
  }
  */
  for(Int_t kk=1;kk<9;kk++)
  {
    pout->MinEnergy[kk]=6.;
    //Large Cell threshold
    if(kk==5 || kk==6)pout->MinEnergy[kk]=6.;
    if(kk>4&&kk<9)pout->MinEnergy[kk]=Eminpout[kk-5];
  }
  Geom* p_Geom=new Geom(p_files);

  TH1F MassH("MassH","MassH",500,0,5.);
  TH1F MassH0("MassH0","MassH0",500,0,5.);
  TH2F MassHSvsZ("MassHSvsZ","MassHSvsZ",500,0,5.,50,0,1.);
  TH1F MassHS("MassHS","MassHS",500,0,5.);
  TH1F EDetH("EDetH","EDetH",120,0,120.);
  TH1F E2phH("E2phH","E2ph",120,0,120.);
  TH1F NphH("NphH","NphH",20,0,20);
  TH1F NphHAll("NphHAll","NphHAll",20,0,20);
  TH2F Mass2vsE2("Mass2vsE2","Mass2vsE2",500,0,5.,120,0,120.);
  TH2F MvEtwo_13("MvEtwo_13","MvEtwo_13",500,0,5.,120,0,120.);
  TH2F EtavEtwo_13("EtavEtwo_13","EtavEtwo_13",100,-5.,5.,120,0,120.);
  TH2F EtavMtwo_13("EtavMtwo_13","EtavMtwo_13",100,-5.,5.,500,0,5.);
  TH2F MvE_Eta_33_38("MvE_Eta_33_38","MvE_Eta_33_38",100,0,2.,100,0.,100.);
  TH2F MvE_Eta_38_43("MvE_Eta_38_43","MvE_Eta_38_43",100,0,2.,100,0.,100.);

  TH2F MvEtwo_24("MvEtwo_24","MvEtwo_24",500,0,5.,120,0,120.);
  TH2F EtavEtwo_24("EtavEtwo_24","EtavEtwo_24",100,-5.,5.,120,0,120.);
  TH2F EtavMtwo_24("EtavMtwo_24","EtavMtwo_24",100,-5.,5.,500,0,5.);
  TH1F NearPi("NearPi","NearPi",120,0,1.2);
  TH2F Near2Pi("Near2Pi","Near2Pi",120,0,1.2,120,0,1.2);
  TH2F M2PivE("M2PivE","M2PivE",60,0.,3.,100,0,100.);
  TH2F MPiEtavE("MPiEtavE","MPiEtavE",20,2.,5.,100,0,100.);
  TH2F MPiEtavZ_E30("MPiEtavZ_E30","MPiEtavE_E30",30,2.,5.,20,0,1);
  TH2F MPiEtavZ_E40("MPiEtavZ_E40","MPiEtavE_E40",30,2.,5.,20,0,1);
  TH2F MPiEtavZ_E50("MPiEtavZ_E50","MPiEtavE_E50",30,2.,5.,20,0,1);
  TH2F MPiEtavZ_E60("MPiEtavZ_E60","MPiEtavE_E60",30,2.,5.,20,0,1);

  TH2F MEtaEtavZ_E30("MEtaEtavZ_E30","MEtaEtavE_E30",30,2.,5.,20,0,1);
  TH2F MEtaEtavZ_E40("MEtaEtavZ_E40","MEtaEtavE_E40",30,2.,5.,20,0,1);
  TH2F MEtaEtavZ_E50("MEtaEtavZ_E50","MEtaEtavE_E50",30,2.,5.,20,0,1);
  TH2F MEtaEtavZ_E60("MEtaEtavZ_E60","MEtaEtavE_E60",30,2.,5.,20,0,1);

  TH2F XYsoft("XYsoft","XYsoft",100,-100.,100.,100,-100.,100);
  TH2F XY("XY","XY",100,-100.,100.,100,-100.,100);
  TH2F XYpi("XYpi","XYpi",100,-100.,100.,100,-100.,100);
  TH2F XYpiHS("XYpiHS","XYpiHS",100,-100.,100.,100,-100.,100);
  TH2F XYEta("XYEta","XYEta",100,-100.,100.,100,-100.,100);
  TH2F XY4Eta("XY4Eta","XY4Eta",100,-100.,100.,100,-100.,100);
  TH2F MvEta("MvEta","MvEta",100,0.,2.,30,2.,5.);
  TH2F MvEta_E30("MvEta_E30","MvEta_E30",100,0.,2.,30,2.,5.);
  TH2F MvEta_E40("MvEta_E40","MvEta_E40",100,0.,2.,30,2.,5.);
  TH2F MvEta_E50("MvEta_E50","MvEta_E50",100,0.,2.,30,2.,5.);
  TH2F MvEta_E60("MvEta_E60","MvEta_E60",100,0.,2.,30,2.,5.);
  TH2F MvEta_E65("MvEta_E65","MvEta_E65",100,0.,2.,30,2.,5.);
  TH2F HSpairMvsE("HSpairMvsE","HSpairMvsE",120,0.,1.2,100,0,100.);
  TH2F HSpairMvsm0("HSpairMvsm0","HSpairMvsm0",120,0.,1.2,41,.0,.82);
  TH2F Mhi3("Mhi3","Mhi3",200,0.,5.,20,0.,100);
  TH2F Mho3("Mho3","Mho3",200,0.,5.,20,0.,100);
  TH2F Mhi3Pt("Mhi3Pt","Mhi3Pt",200,0.,5.,40,0.,4.);
  TH3F Momega("Momega","Momega",12,0,120,12,0,6.,200,0.,5.);
  TH3F MomegaL("MomegaL","MomegaL",12,0,120,12,0,6.,200,0.,5.);
  TH3F MomegaN("MomegaN","MomegaN",12,0,120,12,0,6.,200,0.,5.);
  TH3F Mdal1("Mdal1","Mdal1",120,0.,1.2,120,0.,1.2,40,0.,1.);
  TH3F Mdal2("Mdal2","Mdal2",120,0.,1.2,120,0.,1.2,50,0.,100.);
  TH1F M3BPairM("M3BPairM","M3BPairM",120,0.,1.2);
  TH2F M3BPairEM("M3BPairEM","M3BPairEM",120,0.,1.2,12,0,120.);
  TH2F M3BEtaEM("M3BEtaEM","M3BEtaEM",120,0.,1.2,12,0,120.);
  TH2F M3BEtaEMNPi("M3BEtaEMNPi","M3BEtaEMNPi",120,0.,1.2,12,0,120.);
  TH3F M3Eta("M3Eta","M3Eta",12,0,120,12,0,6.,200,0.,5.);
  TH3F M3EtaL("M3EtaL","M3EtaL",12,0,120,12,0,6.,200,0.,5.);
  TH3F M2vsPt("M2vsPt","M2vsPt",10,0,100,30,0.,6.,200,0.,4.);
  TH3F MomegEPh("MomegEPh","MomegEPh",200,0.,4.,12,0,120,160,-3.2,3.2);
  TH3F MomegEPhN("MomegEPhN","MomegEPhN",200,0.,4.,12,0,120,160,-3.2,3.2);
  TH2F cl2MvE("cl2MvE","cl2MvE",120,0,2.4,120,0,120.);
  TH2F Ecomp1[4];
  TH2F Ecomp2[4];
  TH2F Ecomp[4];
  for(Int_t j0=0;j0<4;j0++)
  {
    char nm0[100];
    sprintf(nm0,"Ecomp1_%d",j0);
    Ecomp1[j0]=TH2F(nm0,nm0,100,1,100,100,1,100);
    sprintf(nm0,"Ecomp2_%d",j0);
    Ecomp2[j0]=TH2F(nm0,nm0,100,1,100,100,1,100);
    sprintf(nm0,"Ecomp_%d",j0);
    Ecomp[j0]=TH2F(nm0,nm0,100,1,100,100,1,100);      
  };
  TFile* out=new TFile(OutFileName,"recreate");

  TTree ThreeTr("ThreeTr","omega event");
  TTree TwoTr("TwoTr","pair event");
  TTree TwoC("TwoC"," Two photon event");
  typedef struct {
    Float_t mpi; 
    Float_t meta; 
    Float_t mlo;
    Float_t mmid;
    Float_t mhi;

    Float_t elo;
    Float_t emid;
    Float_t ehi;

    Float_t Ylo;
    Float_t Ymid;
    Float_t Yhi;

    Float_t Philo;
    Float_t Phimid;
    Float_t Phihi;

    Float_t Ptlo;
    Float_t Ptmid;
    Float_t Pthi;

    Float_t m; 
    Float_t Y; 
    Float_t Phi;
    Float_t Pt;
    Float_t e;

    Float_t dthetalo; 
    Float_t dphilo;
    Float_t dthetamid; 
    Float_t dphimid;
    Int_t spin;
    Char_t BBcSums[5];
    Float_t BBcVertex;
    Float_t Eloa;
    Float_t Phloa;
    Float_t Yloa;
    Float_t Ptloa;
  } tp3_t;
  tp3_t t3o;
  Int_t ClIndex;
  Int_t NExtra;
  Int_t IndCl[20];
  Float_t PtClu[20];
  Float_t PzClu[20];
  Float_t EClu[20];
  Float_t PhiClu[20];
  Float_t MClu[20];
  Float_t M12Clu[20];
  Float_t EtaClu[20];
  Int_t N12Clu[20];
  typedef struct {
    Int_t Nclust;
    Float_t N12;
    Float_t M12;
    Float_t MC12;
    Float_t E12;
    Float_t EC12;
    Float_t C12;
    Float_t Y1;
    Float_t Y2;
    Float_t det12;
    Int_t edge;
    Float_t Phi1;
    Float_t Phi2;
    Float_t Z;
    Float_t Phi;
    Float_t Eta;
    Float_t Pt;
    Float_t Phicm;
    Float_t Ntracks;
    //total variables
    Float_t NMass;
    Float_t NE;
    Float_t NY;
    Float_t NPhicm;
    Float_t NPhi;
    Float_t NPt;
    Float_t Esoft;
    //away 
    Float_t Phiaway;
    Float_t Maway;
    Float_t Yaway;
    Float_t Ptaway;
    Int_t spin;
    Char_t BBcSums[5];
    Float_t BBcVertex;
    Float_t Edet[4];
    Int_t TrigBits;
    Int_t Bunchid7bit;
    Int_t Rnum;
    Int_t Evnum;
    Float_t dthx;
    Float_t dthy;
    Float_t sigmaMaxF1;
    Float_t chiSquareF1;
    Float_t chiSquareF2;
    UInt_t L2sum[2];
    UInt_t lastdsm[8];
    UInt_t Fpde[8];
  } tp2_t;

  tp2_t t2o;
  Int_t NAB;
  Float_t MAB[12];
  Float_t EAB[12];
  Float_t ZAB[12];
  Float_t PtAB[12];
  typedef struct {
    Int_t R1N12;
    Int_t R2N12;
    Int_t RN12C;
    Float_t R1M12;
    Float_t R2M12;
    Float_t R1E12;
    Float_t R2E12;
    Float_t R1Pt;
    Float_t R2Pt;
    Float_t RChiSq1;
    Float_t RChiSq2;
    Float_t RCorRat;
    Int_t RCatag;
    Int_t RNDF1;
    Int_t RNDF2;
    Float_t R1SigMin;
    Float_t R1SigMax;
  } tp2R_t;
  tp2R_t t2oR;
  tp2R_t t2oReset={0,0,0,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1,-1,-1,-1.,-1.};
  typedef struct {
    Float_t M;
    Float_t z;
    Float_t E;
    Float_t Y;
    Float_t Pt;
    Float_t Phi;
    Float_t Cth1;
    Float_t Cth2;
    Float_t Phi1;
    Float_t Phi2;
    Float_t iso;
  } tp2c_t;

  tp2c_t t2c;
  TwoC.Branch("M",&t2c.M,"M/F");
  TwoC.Branch("E",&t2c.E,"E/F");
  TwoC.Branch("z",&t2c.z,"z/F");
  TwoC.Branch("Y",&t2c.Y,"Y/F");
  TwoC.Branch("Pt",&t2c.Pt,"Pt/F");
  TwoC.Branch("Phi",&t2c.Phi,"Phi/F");
  TwoC.Branch("Cth1",&t2c.Cth1,"Cth1/F");
  TwoC.Branch("Cth2",&t2c.Cth2,"Cth2/F");
  TwoC.Branch("Phi1",&t2c.Phi1,"Phi1/F");
  TwoC.Branch("Phi2",&t2c.Phi2,"Phi2/F");
  TwoC.Branch("iso",&t2c.iso,"iso/F");

  TwoTr.Branch("M12",&t2o.M12,"M12/F");
  TwoTr.Branch("MC12",&t2o.MC12,"MC12/F");
  TwoTr.Branch("N12",&t2o.N12,"N12/F");
  TwoTr.Branch("Nclust",&t2o.Nclust,"Nclust/I");
  TwoTr.Branch("NExtra",&NExtra,"NExtra/I");
  TwoTr.Branch("ClIndex",&ClIndex,"ClIndex/I");
  TwoTr.Branch("PhiClu",PhiClu,"PhiClu[NExtra]/F");
  TwoTr.Branch("PtClu",PtClu,"PtClu[NExtra]/F");
  TwoTr.Branch("EClu",EClu,"EClu[NExtra]/F");
  TwoTr.Branch("MClu",MClu,"MClu[NExtra]/F");
  TwoTr.Branch("M12Clu",M12Clu,"M12Clu[NExtra]/F");
  TwoTr.Branch("EtaClu",EtaClu,"EtaClu[NExtra]/F");
  TwoTr.Branch("N12Clu",N12Clu,"N12Clu[NExtra]/I");

  //QTB branch

  TwoTr.Branch("QTNE",&(qtbinfo.QTNE),"BQTNE/B");
  TwoTr.Branch("QTNW",&(qtbinfo.QTNW),"BQTNW/B");   
  TwoTr.Branch("QTEBBCInd",(qtbinfo.QTEBBCInd),"BQTEInd[BQTNE]/B");
  TwoTr.Branch("QTWBBCInd",(qtbinfo.QTWBBCInd),"BQTWInd[BQTNW]/B");

  TwoTr.Branch("QTEBBCTAC",(qtbinfo.QTEBBCTAC),"QTEBTAC[BQTNE]/S");
  TwoTr.Branch("QTWBBCTAC",(qtbinfo.QTWBBCTAC),"QTWBTAC[BQTNW]/S");

  TwoTr.Branch("QTEBBCADC",(qtbinfo.QTEBBCADC),"QTEBADC[BQTNE]/S");
  TwoTr.Branch("QTWBBCADC",(qtbinfo.QTWBBCADC),"QTWBADC[BQTNW]/S");

  TwoTr.Branch("QTBVertex",&(qtbinfo.vertex),"QTBVertex/F");


  //QTRP branch
  TwoTr.Branch("RPE_QTN",&(qtrpinfo.NE),"RPE_QTN/B");
  TwoTr.Branch("RPW_QTN",&(qtrpinfo.NW),"RPW_QTN/B");
  TwoTr.Branch("RPE_Idx",(qtrpinfo.RPE_Idx),"RPE_Idx[RPE_QTN]/B");
  TwoTr.Branch("RPW_Idx",(qtrpinfo.RPW_Idx),"RPW_Idx[RPW_QTN]/B");
  TwoTr.Branch("RPE_TAC",(qtrpinfo.RPE_TAC),"RPE_TAC[RPE_QTN]/S");
  TwoTr.Branch("RPW_TAC",(qtrpinfo.RPW_TAC),"RPW_TAC[RPW_QTN]/S");
  TwoTr.Branch("RPE_ADC",(qtrpinfo.RPE_ADC),"RPE_ADC[RPE_QTN]/S");
  TwoTr.Branch("RPW_ADC",(qtrpinfo.RPW_ADC),"RPW_ADC[RPW_QTN]/S");
  TwoTr.Branch("RPvertex",&(qtrpinfo.vertex),"RPvertex/F");



  TwoTr.Branch("EC12",&t2o.EC12,"EC12/F");
  TwoTr.Branch("E12",&t2o.E12,"E12/F");
  TwoTr.Branch("C12",&t2o.C12,"C12/F");
  TwoTr.Branch("Y1",&t2o.Y1,"Y1/F");
  TwoTr.Branch("Y2",&t2o.Y2,"Y2/F");
  TwoTr.Branch("Phi1",&t2o.Phi1,"Phi1/F");
  TwoTr.Branch("Phi2",&t2o.Phi2,"Phi2/F");
  TwoTr.Branch("Phi",&t2o.Phi,"Phi/F");
  TwoTr.Branch("Z",&t2o.Z,"Z/F");
  TwoTr.Branch("Eta",&t2o.Eta,"Eta/F");
  TwoTr.Branch("Pt",&t2o.Pt,"Pt/F");
  TwoTr.Branch("Phicm",&t2o.Phicm,"Phicm/F");
  TwoTr.Branch("Ntracks",&t2o.Ntracks,"Ntracks/F");
  TwoTr.Branch("det12",&t2o.det12,"det12/F");
  TwoTr.Branch("edge",&t2o.edge,"edge/I");
  //
  TwoTr.Branch("Esoft",&t2o.Esoft,"Esoft/F");
  TwoTr.Branch("NMass",&t2o.NMass,"NMass/F");
  TwoTr.Branch("NE",&t2o.NE,"NE/F");
  TwoTr.Branch("NY",&t2o.NY,"NY/F");
  TwoTr.Branch("NPhicm",&t2o.NPhicm,"NPhicm/F");
  TwoTr.Branch("NPhi",&t2o.NPhi,"NPhi/F");
  TwoTr.Branch("NPt",&t2o.NPt,"NPt/F");

  TwoTr.Branch("Maway",&t2o.Maway,"Maway/F");
  TwoTr.Branch("Phiaway",&t2o.Phiaway,"Phiaway/F");
  TwoTr.Branch("Yaway",&t2o.Yaway,"Yaway/F");
  TwoTr.Branch("Ptaway",&t2o.Ptaway,"NPtaway/F");
  TwoTr.Branch("spin",&t2o.spin,"spin/I");
  TwoTr.Branch("BBcSums",t2o.BBcSums,"BBcSums[5]/b");
  TwoTr.Branch("BBcVertex",&t2o.BBcVertex,"BBcVertex/F");
  TwoTr.Branch("EDet",t2o.Edet,"Edet[4]/F");
  TwoTr.Branch("TrigBits",&t2o.TrigBits,"TrigBits/I");
  TwoTr.Branch("Bunchid7bit",&t2o.Bunchid7bit,"Bunchid7bit/I");
  TwoTr.Branch("Rnum",&t2o.Rnum,"Rnum/I");
  TwoTr.Branch("Evnum",&t2o.Evnum,"Evnum/I");
  TwoTr.Branch("dthx",&t2o.dthx,"dthx/F");
  TwoTr.Branch("dthy",&t2o.dthy,"dthy/F");
  TwoTr.Branch("sigmaMaxF1",&t2o.sigmaMaxF1,"sigmaMaxF1/F");
  TwoTr.Branch("chiSquareF1",&t2o.chiSquareF1,"chiSquareF1/F");
  TwoTr.Branch("chiSquareF2",&t2o.chiSquareF2,"chiSquareF2/F");
  TwoTr.Branch("L2sum",t2o.L2sum,"L2sum[2]/I");
  TwoTr.Branch("lastdsm",t2o.lastdsm,"lastdsm[8]/i");
  TwoTr.Branch("Fpde",t2o.Fpde,"Fpde[8]/i");

  TwoTr.Branch("NAB",&NAB,"NAB/I");
  TwoTr.Branch("MAB",MAB,"MAB[NAB]/F");
  TwoTr.Branch("EAB",EAB,"EAB[NAB]/F");
  TwoTr.Branch("ZAB",ZAB,"ZAB[NAB]/F");
  TwoTr.Branch("PtAB",PtAB,"PtAB[NAB]/F");


  /*
     ReRec
     */  
  TwoTr.Branch("R1N12",&t2oR.R1N12,"R1N12/I");
  TwoTr.Branch("R2N12",&t2oR.R2N12,"R2N12/I");
  TwoTr.Branch("RN12C",&t2oR.RN12C,"RN12C/I");
  TwoTr.Branch("R1M12",&t2oR.R1M12,"R1M12/F");
  TwoTr.Branch("R2M12",&t2oR.R2M12,"R2M12/F");
  TwoTr.Branch("R1E12",&t2oR.R1E12,"R1E12/F");
  TwoTr.Branch("R2E12",&t2oR.R2E12,"R2E12/F");
  TwoTr.Branch("R1Pt",&t2oR.R1Pt,"R1Pt/F");
  TwoTr.Branch("R2Pt",&t2oR.R2Pt,"R2Pt/F");
  TwoTr.Branch("RChiSq1",&t2oR.RChiSq1,"RChiSq1/F");
  TwoTr.Branch("RChiSq2",&t2oR.RChiSq2,"RChiSq2/F");
  TwoTr.Branch("RCorRat",&t2oR.RCorRat,"RCorRat/F");
  TwoTr.Branch("RCatag",&t2oR.RCatag,"RCatag/I");
  TwoTr.Branch("RNDF1",&t2oR.RNDF1,"RNDF1/I");
  TwoTr.Branch("RNDF2",&t2oR.RNDF2,"RNDF2/I");
  TwoTr.Branch("R1SigMax",&t2oR.R1SigMax,"R1SigMax/F");
  TwoTr.Branch("R1SigMin",&t2oR.R1SigMin,"R1SigMin/F");


  // RP data from MuDsts
  /*
   * BEFORE MFILES UPGRADE OCTOBER 2016
   *
  TwoTr.Branch("RP_n_tracks",&(pout->RP_n_tracks),"RP_n_tracks/I");
  TwoTr.Branch("RP_n_trackpoints",&(pout->RP_n_trackpoints),"RP_n_trackpoints/I");
  TwoTr.Branch("RP_t_index",pout->RP_t_index,"RP_t_index[RP_n_tracks]/I"); // track index number
  TwoTr.Branch("RP_t_branch",pout->RP_t_branch,"RP_t_branch[RP_n_tracks]/I"); // RP branch (0=EU 1=ED 2=WU 3=WD)
  TwoTr.Branch("RP_t_type",pout->RP_t_type,"RP_t_type[RP_n_tracks]/I"); // track type (see next line)
  TwoTr.Branch("RP_t_planesUsed",pout->RP_t_planesUsed,"RP_t_planesUsed[RP_n_tracks]/I"); // number of SSD planes hit by track points in track
  TwoTr.Branch("RP_t_p",pout->RP_t_p,"RP_t_p[RP_n_tracks]/D"); // momentum
  TwoTr.Branch("RP_t_pt",pout->RP_t_pt,"RP_t_pt[RP_n_tracks]/D"); // transverse momentum
  TwoTr.Branch("RP_t_eta",pout->RP_t_eta,"RP_t_eta[RP_n_tracks]/D"); // pseudorapidity
  TwoTr.Branch("RP_t_time",pout->RP_t_time,"RP_t_time[RP_n_tracks]/D"); // time of track detection
  TwoTr.Branch("RP_t_theta",pout->RP_t_theta,"RP_t_theta[RP_n_tracks]/D"); // polar angle at RP according to STAR coord sys
  TwoTr.Branch("RP_t_thetaRP",pout->RP_t_thetaRP,"RP_t_thetaRP[RP_n_tracks]/D"); // polar angle at RP according to STAR survey
  TwoTr.Branch("RP_t_phi",pout->RP_t_phi,"RP_t_phi[RP_n_tracks]/D"); // azimuth at RP according to STAR coord sys
  TwoTr.Branch("RP_t_phiRP",pout->RP_t_phiRP,"RP_t_phiRP[RP_n_tracks]/D"); // azimuth at RP according to STAR survey
  TwoTr.Branch("RP_t_t",pout->RP_t_t,"RP_t_t[RP_n_tracks]/D"); // squared 4-momentum transfer
  TwoTr.Branch("RP_t_xi",pout->RP_t_xi,"RP_t_xi[RP_n_tracks]/D"); // fractional momentum loss (pbeam-p)/pbeam
  TwoTr.Branch("RP_t_gold",pout->RP_t_gold,"RP_t_gold[RP_n_tracks]/O"); // my track quality variable (2 track points in all 2x4=8 Si planes)
  TwoTr.Branch("RP_p0_tpExists",pout->RP_p_tpExists[0],"RP_p0_tpExists[RP_n_tracks]/O"); // true if track point 0 exists
  TwoTr.Branch("RP_p0_RPid",pout->RP_p_RPid[0],"RP_p0_RPid[RP_n_tracks]/I"); // RP id (see next line)
  TwoTr.Branch("RP_p0_quality",pout->RP_p_quality[0],"RP_p0_quality[RP_n_tracks]/I"); // track point quality (see next line)
  TwoTr.Branch("RP_p0_x",pout->RP_p_x[0],"RP_p0_x[RP_n_tracks]/D"); // STAR survey coords x-position
  TwoTr.Branch("RP_p0_y",pout->RP_p_y[0],"RP_p0_y[RP_n_tracks]/D"); // STAR survey coords y-position
  TwoTr.Branch("RP_p0_z",pout->RP_p_z[0],"RP_p0_z[RP_n_tracks]/D"); // STAR survey coords z-position
  TwoTr.Branch("RP_p1_tpExists",pout->RP_p_tpExists[1],"RP_p1_tpExists[RP_n_tracks]/O"); // true if track point 1 exists
  TwoTr.Branch("RP_p1_RPid",pout->RP_p_RPid[1],"RP_p1_RPid[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_quality",pout->RP_p_quality[1],"RP_p1_quality[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_x",pout->RP_p_x[1],"RP_p1_x[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_y",pout->RP_p_y[1],"RP_p1_y[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_z",pout->RP_p_z[1],"RP_p1_z[RP_n_tracks]/D");
  */

  TwoTr.Branch("RP_n_tracks",&(pout->RP_n_tracks),"RP_n_tracks/I");
  TwoTr.Branch("RP_n_trackpoints",&(pout->RP_n_trackpoints),"RP_n_trackpoints/I");

  // tracks
  TwoTr.Branch("RP_t_index",pout->RP_t_index,"RP_t_index[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_t_branch",pout->RP_t_branch,"RP_t_branch[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_t_type",pout->RP_t_type,"RP_t_type[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_t_planesUsed",pout->RP_t_planesUsed,"RP_t_planesUsed[RP_n_tracks]/i"); 
  TwoTr.Branch("RP_t_p",pout->RP_t_p,"RP_t_p[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_pt",pout->RP_t_pt,"RP_t_pt[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_eta",pout->RP_t_eta,"RP_t_eta[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_time",pout->RP_t_time,"RP_t_time[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_theta",pout->RP_t_theta,"RP_t_theta[RP_n_tracks][3]/D"); 
  TwoTr.Branch("RP_t_thetaRP",pout->RP_t_thetaRP,"RP_t_thetaRP[RP_n_tracks][3]/D"); 
  TwoTr.Branch("RP_t_phi",pout->RP_t_phi,"RP_t_phi[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_phiRP",pout->RP_t_phiRP,"RP_t_phiRP[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_t",pout->RP_t_t,"RP_t_t[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_xi",pout->RP_t_xi,"RP_t_xi[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_t_gold",pout->RP_t_gold,"RP_t_gold[RP_n_tracks]/O"); 
  TwoTr.Branch("RP_t_isBad",pout->RP_t_isBad,"RP_t_isBad[RP_n_tracks]/O"); 
  TwoTr.Branch("RP_t_qualHash",pout->RP_t_qualHash,"RP_t_qualHash[RP_n_tracks]/D"); 

  // track point 0
  TwoTr.Branch("RP_p0_tpExists",pout->RP_p_tpExists[0],"RP_p0_tpExists[RP_n_tracks]/O"); 
  TwoTr.Branch("RP_p0_RPid",pout->RP_p_RPid[0],"RP_p0_RPid[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_p0_clustid_s1",pout->RP_p_clustid_s1[0],"RP_p0_clustid_s1[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_p0_clustid_s2",pout->RP_p_clustid_s2[0],"RP_p0_clustid_s2[RP_n_tracks]/I");
  TwoTr.Branch("RP_p0_clustid_s3",pout->RP_p_clustid_s3[0],"RP_p0_clustid_s3[RP_n_tracks]/I");
  TwoTr.Branch("RP_p0_clustid_s4",pout->RP_p_clustid_s4[0],"RP_p0_clustid_s4[RP_n_tracks]/I");
  TwoTr.Branch("RP_p0_quality",pout->RP_p_quality[0],"RP_p0_quality[RP_n_tracks]/I"); 
  TwoTr.Branch("RP_p0_planesUsed",pout->RP_p_planesUsed[0],"RP_p0_planesUsed[RP_n_tracks]/i"); 
  TwoTr.Branch("RP_p0_x",pout->RP_p_x[0],"RP_p0_x[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_p0_y",pout->RP_p_y[0],"RP_p0_y[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_p0_z",pout->RP_p_z[0],"RP_p0_z[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_p0_time_pmt1",pout->RP_p_time_pmt1[0],"RP_p0_time_pmt1[RP_n_tracks]/D"); 
  TwoTr.Branch("RP_p0_time_pmt2",pout->RP_p_time_pmt2[0],"RP_p0_time_pmt2[RP_n_tracks]/D"); 
  // track point 1
  TwoTr.Branch("RP_p1_tpExists",pout->RP_p_tpExists[1],"RP_p1_tpExists[RP_n_tracks]/O");
  TwoTr.Branch("RP_p1_RPid",pout->RP_p_RPid[1],"RP_p1_RPid[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_clustid_s1",pout->RP_p_clustid_s1[1],"RP_p1_clustid_s1[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_clustid_s2",pout->RP_p_clustid_s2[1],"RP_p1_clustid_s2[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_clustid_s3",pout->RP_p_clustid_s3[1],"RP_p1_clustid_s3[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_clustid_s4",pout->RP_p_clustid_s4[1],"RP_p1_clustid_s4[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_quality",pout->RP_p_quality[1],"RP_p1_quality[RP_n_tracks]/I");
  TwoTr.Branch("RP_p1_planesUsed",pout->RP_p_planesUsed[1],"RP_p1_planesUsed[RP_n_tracks]/i");
  TwoTr.Branch("RP_p1_x",pout->RP_p_x[1],"RP_p1_x[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_y",pout->RP_p_y[1],"RP_p1_y[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_z",pout->RP_p_z[1],"RP_p1_z[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_time_pmt1",pout->RP_p_time_pmt1[1],"RP_p1_time_pmt1[RP_n_tracks]/D");
  TwoTr.Branch("RP_p1_time_pmt2",pout->RP_p_time_pmt2[1],"RP_p1_time_pmt2[RP_n_tracks]/D");


  //
  ThreeTr.Branch("mpi",&t3o.mpi,"mpi/F");
  ThreeTr.Branch("meta",&t3o.meta,"meta/F");
  ThreeTr.Branch("mlo",&t3o.mlo,"mlo/F");
  ThreeTr.Branch("mmid",&t3o.mmid,"mmid/F");
  ThreeTr.Branch("mhi",&t3o.mhi,"mhi/F");

  ThreeTr.Branch("m",&t3o.m,"m/F");
  ThreeTr.Branch("Y",&t3o.Y,"Y/F");
  ThreeTr.Branch("Pt",&t3o.Pt,"Pt/F");
  ThreeTr.Branch("Phi",&t3o.Phi,"Phi/F");
  ThreeTr.Branch("e",&t3o.e,"e/F");

  ThreeTr.Branch("dthetalo",&t3o.dthetalo,"dthetalo/F");
  ThreeTr.Branch("dphilo",&t3o.dphilo,"dphilo/F");
  ThreeTr.Branch("dthetamid",&t3o.dthetamid,"dthetamid/F");
  ThreeTr.Branch("dphimid",&t3o.dphimid,"dphimid/F");

  ThreeTr.Branch("Ylo",&t3o.Ylo,"Ylo/F");
  ThreeTr.Branch("Ymid",&t3o.Ymid,"Ymid/F");
  ThreeTr.Branch("Yhi",&t3o.Yhi,"Yhi/F");

  ThreeTr.Branch("elo",&t3o.elo,"elo/F");
  ThreeTr.Branch("emid",&t3o.emid,"emid/F");
  ThreeTr.Branch("ehi",&t3o.ehi,"ehi/F");

  ThreeTr.Branch("Ptlo",&t3o.Ptlo,"Ptlo/F");
  ThreeTr.Branch("Ptmid",&t3o.Ptmid,"Ptmid/F");
  ThreeTr.Branch("Pthi",&t3o.Pthi,"Pthi/F");

  ThreeTr.Branch("Philo",&t3o.Philo,"Philo/F");
  ThreeTr.Branch("Phimid",&t3o.Phimid,"Phimid/F");
  ThreeTr.Branch("Phihi",&t3o.Phihi,"Phihi/F");
  ThreeTr.Branch("spin",&t3o.spin,"spin/I");
  ThreeTr.Branch("BBcSums",t3o.BBcSums,"BBcSums[5]/b");
  ThreeTr.Branch("BBcVertex",&t3o.BBcVertex,"BBcVertex/F");
  ThreeTr.Branch("Eloa",&t3o.Eloa,"Eloa/F");
  ThreeTr.Branch("Phloa",&t3o.Phloa,"Phloa/F");
  ThreeTr.Branch("Yloa",&t3o.Yloa,"Yloa/F");
  ThreeTr.Branch("Ptloa",&t3o.Ptloa,"Ptloa/F");
  Int_t nentries=pout->nentries;
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;  
  printf("Events to Process=%d \n",nentries);
  if(NumberEventsToProcess>pout->nentries)nentries=pout->nentries;
  for(Int_t iev=0;iev<nentries;iev++)
  {
    if( (iev%10000) ==0)printf( "%d Events Read \n ",iev);
    //      printf( "%d Events Read \n ",iev);
    if( (iev%100000)==0&&iev>50000)
    {
      TwoTr.Write();
      cout<<"Root File Name ="<<gROOT->GetFile()->GetName()<<"\n";
      gROOT->GetFile()->Purge();
    };
    pout->GetEntry(iev);


    t2o.Bunchid7bit=pout->Bunchid7bit;
    t2o.Rnum=pout->Rnum;
    t2o.Evnum=pout->ievt;
    t2o.TrigBits=pout->TrigBits;	  
    Bool_t TrigPass=true;
    if(TrigBitsSelect!=0)
    {
      if((TrigBitsSelect&(pout->TrigBits))==0)TrigPass=false;
      if((TrigBitsVeto&(pout->TrigBits))!=0)TrigPass=false;

    };
    if(TrigPass)
    {
      NphH.Fill(pout->Nphotons); // histogram hard photons
      Float_t Ed1[4],En1[4],En2[4],En[4];
      Int_t Nphh[4];
      for(Int_t j0=0;j0<4;j0++)
      {
        Ed1[j0]=En1[j0]=En2[j0]=En[j0]=0.;
        Nphh[j0]=0;
      };
      TIter diter(pout->dlist);
      while(LVec* v=(LVec*) diter())
      {
        if(v->Nstb>=0 && v->Nstb<4 &&Ed1[v->Nstb]==0.)
        {
          Ed1[v->Nstb]=v->E();
        };
      };
      TIter viter(pout->vlist);
      while(LVec* v=(LVec*) viter())
      {
        Int_t j0;
        if((j0=v->Nstb)>=0 && v->Nstb<4)
        {
          if(En1[j0]==0){En1[j0]+=v->E();}
          else if (En2[j0]==0)En2[j0]=En1[j0]+v->E();
          En[j0]+=v->E();
          Nphh[j0]++;
        };
      };
      for(Int_t j0=0;j0<4;j0++)
      {
        if(Nphh[j0]==1)Ecomp1[j0].Fill(Ed1[j0],En1[j0]);
        if(Nphh[j0]==2)Ecomp2[j0].Fill(Ed1[j0],En2[j0]);
        Ecomp[j0].Fill(Ed1[j0],En[j0]);
      };
      TIter next(pout->vlist);
      while(LVec* v=(LVec*) next())
      {
        TVector3 vhit=pout->PosInDet(v,p_geom);
        XY.Fill(vhit.X(),vhit.Y());
      };
      TIter nexts(pout->softlist);
      while(LVec* v=(LVec*) nexts())
      {
        TVector3 vhit=pout->PosInDet(v,p_geom);
        XYsoft.Fill(vhit.X(),vhit.Y());
      };
      pout->AllToScratch(false); //make a scratch list of all hard 
      Int_t nphall;

      if((nphall=pout->scratchlist->GetEntries())==2)
      {
        Float_t mhardsoft=pout->SumScratch().Mag();
        MassH0.Fill(mhardsoft);
        LVec* v2pa;
        LVec* v2pb;
        v2pa=(LVec*) pout->scratchlist->First();//first photon
        v2pb=(LVec*) pout->scratchlist->After(v2pa);//other photon
        if(v2pa->E()<v2pb->E())
        {
          LVec* ltemp=v2pa;
          v2pa=v2pb;
          v2pb=ltemp;
        };

        t2c.M=((*v2pa)+(*v2pb)).Mag();
        t2c.E=((*v2pa)+(*v2pb)).E();
        t2c.Y=((*v2pa)+(*v2pb)).PseudoRapidity();
        t2c.Pt=((*v2pa)+(*v2pb)).Pt();
        t2c.Phi=((*v2pa)+(*v2pb)).Phi();
        t2c.Cth1=cos(v2pa->Theta());
        t2c.Cth2=cos(v2pb->Theta());
        t2c.Phi1=v2pa->Phi();
        t2c.Phi2=v2pb->Phi();

        t2c.iso=1000.;
        TIter nsft(pout->softlist);
        LVec* sphot;
        while(sphot=(LVec*) nsft())
        {
          t2c.iso=TMath::Min((Float_t) t2c.iso,
              (Float_t) sqrt(pow(sphot->Phi()-v2pa->Phi(),2)
                +pow(sphot->PseudoRapidity()-
                  v2pa->PseudoRapidity(),2)));
          t2c.iso=TMath::Min((Float_t) t2c.iso,
              (Float_t) sqrt(pow(sphot->Phi()-v2pb->Phi(),2)
                +pow(sphot->PseudoRapidity()-
                  v2pb->PseudoRapidity(),2)));
        };
        //TwoC.Fill();
      };
      NphHAll.Fill(nphall);
      pout->ClearScratch();
      TString CLYPHI=gSystem->Getenv("CLYPHI");
      bool clyphi=false;
      if(CLYPHI!="" && CLYPHI!="false")
      {
        clyphi=true;
      }

      TString ClusterSizeRad=gSystem->Getenv("ClusterSizeRad");
      float clsize=0.1;
      if(ClusterSizeRad != "")
      {
        sscanf((const char*) ClusterSizeRad,"%f",&clsize);
        if(clsize<.001 || clsize>1.)return -1;
        if(iev<2){
          printf("AnalTools::clsize=%f \n",clsize);
          if(clyphi)
          {
            printf("AnalTools::clyphi=true\n");
          } else printf("AnalTools::clyphi=false\n");
        };
      };


      pout->AllToScratch(false); // make a scratch list of all hard phots
      pout->scratchlist->Sort();

      pout->ClusterwithYPhi(clyphi);// cluster for angular size
      //	  TObjArray* Clust=pout->ClusterScratch(clsize);
      TObjArray* Clust=pout->ClusterScratch(clsize);
      if(iev<2 && clyphi)printf("pout->ClusterwithYPhi(true)\n");

      //   TObjArray* Clust=pout->ClusterScratch(.2);
      //	     pout->ClusterwithYPhi(false);// cluster for angular size
      // pout->ClusterwithYPhi(true);// cluster for angular size
      //      TObjArray* Clust=pout->ClusterScratch(10.);
      //TObjArray* Clust=pout->ClusterScratch(.07); //pre Jun 2011
      //	  TObjArray* Clust=pout->ClusterScratch(.035);
      // Oct 6	  TObjArray* Clust=pout->ClusterScratch(.07);
      // Nov 1 TObjArray* Clust=pout->ClusterScratch(.1);
      //	  TObjArray* Clust=pout->ClusterScratch(.035);

      TLorentzVector Vscr=pout->SumScratch();
      t2o.Nclust=Clust->GetEntries();
      TObjArray* ctmp;
      LVec* ctmpL;
      TObjArray* ClustL=new TObjArray(Clust->GetEntries());
      ClustL->SetOwner(1);
      TIter nxt0(Clust);
      int vlorind=100;
      while((ctmp=(TObjArray*) nxt0()))
      {
        ctmp->Sort();
        LVec* lvor=new LVec(vlorind++,TLorentzVector(pout->SumList(ctmp)));
        ClustL->Add(lvor);	      
      };
      ClustL->Sort();
      TIter nxt(ClustL);
      Bool_t PiFound=false;
      TLorentzVector vz(0,0,1,1);
      ClIndex=-1;
      while((ctmpL=(LVec*) nxt())) //loop over coneclusters (Energy sorted)
      {
        ClIndex++;
        NExtra=0;
        ctmp=(TObjArray*) Clust->At(ctmpL->Vtype);//next conecluster (E s)
        ctmp->Sort();
        //  QT BBC info in first cluster
        if(ClIndex==0)
        {
          qtbinfo=pout->qtbbc;
          qtrpinfo=pout->qtrp;
          //printf("Event =%d :\n",iev);
          //		  qtbinfo.Print();

        }
        else
        { 
          qtbinfo.SetZero();
          qtrpinfo.SetZero();
        };
        //	      printf("ctmp is elem %d of Clust\n",ctmpL->Vtype);
        if(ClIndex==0&&t2o.Nclust>1)
        {//only for high energy conecluster include other cluster 4 mom
          NExtra=t2o.Nclust;
          if(NExtra>19)NExtra=19;
          LVec* ctmpL1;
          TIter nxt1(ClustL);
          Int_t clindex=0;
          while((ctmpL1=(LVec*) nxt1()))
          {//other conecluster 4mom list is energy sorted
            if(clindex<NExtra)
            {
              PtClu[clindex]=ctmpL1->Pt();
              EClu[clindex]=ctmpL1->E();
              PhiClu[clindex]=ctmpL1->Phi();
              MClu[clindex]=ctmpL1->Mag();
              EtaClu[clindex]=ctmpL1->PseudoRapidity();
              clindex++;
              TObjArray* ctmp1=(TObjArray*) Clust->At(ctmpL1->Vtype);//next conecluster (E s)
              //printf("  num=%d    ctmp1 is elem %d of Clust\n",clindex-1,ctmpL1->Vtype);

              N12Clu[clindex-1]=ctmp1->GetEntries();
              // printf("N12Clu[%d]=%d \n",clindex,N12Clu[clindex-1]);
              M12Clu[clindex]=0;
              if(N12Clu[clindex-1]>1)
              {
                M12Clu[clindex-1]=(*((LVec*) ctmp1->At(0))+*((LVec*) ctmp1->At(1))).Mag();
                //printf("    M12Clu=%f\n",M12Clu[clindex]);
              }
            }
          };


        };
        //	      printf("Getting cluster (Vtype) %d ctmpL E=%f ctmp E=%F \n",ctmpL->Vtype,ctmpL->E(),pout->SumList(ctmp).E());
        NAB=0;
        if(int ctmpN=ctmp->GetEntries())
        {
          for(int jctmp=0;(jctmp<ctmpN)&&(jctmp<5);jctmp++)
          {
            for(int kctmp=0;kctmp<jctmp;kctmp++)
            {
              TLorentzVector* cva=(TLorentzVector*) ctmp->At(jctmp);
              TLorentzVector* cvb=(TLorentzVector*) ctmp->At(kctmp);
              MAB[NAB]=(*cva +*cvb).Mag();
              EAB[NAB]=(*cva+*cvb).E();
              PtAB[NAB]=(*cva+*cvb).Pt();
              ZAB[NAB]=0;
              if(EAB[NAB]>0)ZAB[NAB]=((*cvb-*cva).E())/EAB[NAB];
              NAB++;
            }


          }
        };
        Int_t ClusterIndexOf=Clust->IndexOf(ctmp);
        t2o.Yaway=0.;
        t2o.Maway=0.;
        t2o.Phiaway=0.;
        t2o.Ptaway=0.;
        t2o.Ntracks=pout->scratchlist->GetEntries();
        t2o.N12=0;
        t2o.M12=0;
        t2o.E12=0;
        t2o.MC12=0;
        t2o.EC12=0;
        t2o.Pt=0; 
        //	      t2o.Nclust=0;
        t2o.det12=0;
        t2o.edge=0.;
        t2o.Phi1=0.;
        t2o.Phi2=0.;
        t2o.Z=0.;
        t2o.Phi=0.;
        t2o.Eta=0.;
        t2o.Phicm=0.;
        t2o.NMass=0.;
        t2o.NE=0.;
        t2o.NY=0.;
        t2o.NPhicm=0.;
        t2o.NPhi=0.;
        t2o.NPt=0.;
        t2o.Esoft=0.;
        t2o.spin=-100.;
        t2o.dthx=0.;
        t2o.dthy=0.;
        for(Int_t qq=0; qq<2; qq++) t2o.L2sum[qq]=0;
        for(Int_t qq=0; qq<8; qq++) t2o.lastdsm[qq]=0;
        for(Int_t qq=0; qq<8; qq++) t2o.Fpde[qq]=0;
        t2oR=t2oReset;
        for(int jjj=0;jjj<5;jjj++)t2o.BBcSums[jjj]=0.;
        t2o.BBcVertex=0;
        for(int jjj=0;jjj<4;jjj++)t2o.Edet[jjj]=0.;
        t2o.Bunchid7bit=pout->Bunchid7bit;
        t2o.Rnum=pout->Rnum;
        if(t2o.Rnum!=fill.CurrentRunNumber)
        {fill.SetFillNumberforRun(t2o.Rnum);};

        t2o.TrigBits=pout->TrigBits;	  
        t2o.NMass=Vscr.Mag();
        t2o.NE=Vscr.E();
        t2o.NPt=Vscr.Pt();
        t2o.NY=0;
        if(t2o.NPt!=0)t2o.NY=Vscr.PseudoRapidity();
        t2o.NPhi=Vscr.Phi();
        t2o.sigmaMaxF1=-1.;
        t2o.chiSquareF1=-1.;
        t2o.chiSquareF2=-1.;
        for(Int_t qq=0; qq<2; qq++) t2o.L2sum[qq]=pout->L2sum[qq];
        for(Int_t qq=0; qq<8; qq++) t2o.lastdsm[qq]=pout->lastdsm[qq];
        for(Int_t qq=0; qq<8; qq++) t2o.Fpde[qq]=pout->Fpde[qq];

        if(ctmp->GetEntries()>0  )
        {
          TLorentzVector tv=pout->SumList(ctmp);
          t2o.NPhicm=pout->reframe(tv,Vscr,vz).Phi();
          t2o.Ntracks=pout->scratchlist->GetEntries();
          t2o.N12=ctmp->GetEntries();
          t2o.spin=pout->spin;
          //Spin repair
          t2o.spin=0;
          int bspn=  fill.BlueSpin(t2o.Bunchid7bit);
          int yspn=fill.YellowSpin(t2o.Bunchid7bit);
          //		  printf("run=%d bnch=%d bspn=%d yspn=%d \n",t2o.Rnum,t2o.Bunchid7bit,bspn,yspn);
          if(abs(bspn)==1)t2o.spin+=(bspn+1);
          if(abs(yspn)==1)t2o.spin+=(yspn+1)/2;
          if(bspn*yspn==0)t2o.spin+=32;
          if(bspn==0&&yspn==0)t2o.spin=64;
          //end Spin repair
          t2o.BBcVertex=pout->BBcVertex;
          TIter dinxt(pout->dlist);
          LVec* detE;
          while(detE=(LVec*) dinxt())
          {
            if(pout->WhatEW(detE)==2)
            {
              Int_t nstb0=pout->WhatNSTB(detE);
              if(nstb0>0 && nstb0<5)
              {
                t2o.Edet[nstb0-1]=detE->E();
              };
            };
          };	      
          if(tv.Pt()==0)continue;

          t2o.M12=tv.Mag(); //M12 set to cluster mass
          t2o.E12=tv.E(); //E12 set to cluster Energy
          t2o.MC12=t2o.M12;
          t2o.EC12=t2o.E12;
          t2o.Pt=tv.Pt();
          LVec* tv1=(LVec*) ctmp->First();

          t2o.Esoft=pout->ClusterSoftE(ClusterIndexOf,clsize).E();
          TVector3 uv1=tv1->Vect();
          uv1.SetMag(1.);
          t2o.Y1=0;
          if(tv1->Pt()==0)continue;
          LVec* EOrder[2];
          EOrder[0]=0;
          EOrder[1]=0;
          Float_t EV[2];
          EV[0]=0.;
          EV[1]=0;
          TIter sset(ctmp);
          LVec* tvv1;
          t2o.Z=1.;
          while(tvv1=(LVec*) sset() )
          {
            if(EOrder[0]==0 ||tvv1->E()>EV[0])
            {
              if(EV[0]>0){EV[1]=EV[0];EOrder[1]=EOrder[0];};
              EV[0]=tvv1->E();
              EOrder[0]=tvv1;
            }
            else if(EOrder[1]==0 || tvv1->E()>EV[1])
            {
              EV[1]=tvv1->E();
              EOrder[1]=tvv1;
            }

          }; 
          if(EOrder[0])
          {
            t2o.Y1=EOrder[0]->PseudoRapidity();
            t2o.Phi1=EOrder[0]->Phi();
          };
          if(EOrder[1])
          {
            t2o.E12=EOrder[0]->E()+EOrder[1]->E();
            t2o.Y2=EOrder[1]->PseudoRapidity();
            if(t2o.E12>0)t2o.Z=fabs(EOrder[0]->E()-EOrder[1]->E())/t2o.E12;
            t2o.Phi2=EOrder[1]->Phi();

            t2o.M12=(*EOrder[0]+*EOrder[1]).Mag();
            TVector3 uv=tv.Vect();
            TVector3 uv1=EOrder[0]->Vect();
            TVector3 uv2=EOrder[1]->Vect();
            uv.SetMag(1.);
            uv1.SetMag(1.);
            uv2.SetMag(1.);
            t2o.dthx=uv1.X()-uv.X();
            t2o.dthy=uv1.Y()-uv.Y();

          };
          t2o.Eta=tv.PseudoRapidity();	      
          t2o.Phi=tv.Phi();

          LVec* tv2;
          if(ctmp->GetEntries()==2)
          {
            tv1=EOrder[0];
            tv2=EOrder[1];
            if(tv1->E()<tv2->E())
            {
              printf("Order of 1 and 2 messed up\n");
              LVec* ltemp=tv1;
              tv1=tv2;
              tv2=ltemp;
            };
            // check for pion;
            /*
               tv1 - high energy photon
               tv2 - low energy photon
               */

            TIter nsoft(pout->softlist);
            LVec* lso;
            // lso is the LVec for the sum of all soft photons
            Bool_t bloop=true;
            TLorentzVector xtv1=*tv1;
            TLorentzVector xtv2=*tv2;

            while((lso=(LVec*) nsoft())&&bloop)
            {

              if(lso->E()>2)
              {
                if(fabs((*lso+*tv1).Mag()-.135)<.02)
                {
                  xtv1=(*tv1+*lso);
                  PiFound=true;
                  bloop=false;
                }
                else
                {
                  if(fabs((*lso+*tv2).Mag()-.135)<.02)
                  {
                    xtv2=(*tv2+*lso);
                    PiFound=true;
                    bloop=false;
                  };
                };
              };
            };
            //		  tv=xtv1+xtv2;
            uv1=tv1->Vect();
            uv1.SetMag(1.);
            TVector3 uv2=tv2->Vect();
            uv2.SetMag(1.);
            t2o.C12=uv1.Dot(uv2);
            cl2MvE.Fill(tv.Mag(),tv.E());
            t2o.M12=tv.Mag();
            //		  PiFound=false;
            if(PiFound)t2o.M12=-t2o.M12;
            t2o.E12=tv.E();
            if(tv1->Pt()==0)continue;
            t2o.Y1=tv1->PseudoRapidity();
            if(tv2->Pt()==0)continue;
            t2o.Y2=tv2->PseudoRapidity();
            if(tv.Pt()==0)continue;		  
            t2o.Eta=tv.PseudoRapidity();	      

            t2o.Phi=tv.Phi();
            t2o.Z=1.;
            if(t2o.E12>0)t2o.Z=fabs(tv1->E()-tv2->E())/t2o.E12;
            t2o.det12=(tv1->Iew-1)*6+tv1->Nstb-1+
              ((tv2->Iew-1)*6+tv2->Nstb-1)*100;

            t2o.edge=pout->NearEdge(tv1 ,p_geom,.52)+
              10*pout->NearEdge(tv2,p_geom,.52);
            t2o.Phi1=tv1->Phi();
            t2o.Phi2=tv2->Phi();
            t2o.Pt=tv.Pt();

            TVector3 ttv=pout->reframe(*tv1,tv,vz);
            t2o.Phicm=ttv.Phi();
            t2o.Ntracks=pout->scratchlist->GetEntries();
            for(Int_t jj=0;jj<5;jj++){t2o.BBcSums[jj]=pout->BBcSums[jj];};
          };
          if(t2o.Ntracks>t2o.N12)
          {		  
            t2o.Ptaway=(Vscr-tv).Pt();
            t2o.Phiaway=(Vscr-tv).Phi();
            if(!(t2o.Ptaway==0))
            {
              t2o.Yaway=(Vscr-tv).PseudoRapidity();
              t2o.Maway=(Vscr-tv).Mag();
            };
          };
          //  printf("ctmp->GetEntries()=%d\n",ctmp->GetEntries());

          if(true&&(fabs(ctmp->GetEntries()-1.5))<.75 &&Enable_ReRec && fabs(t2o.M12)<1.2 )
          {
            bool N12eq2=false;
            if(ctmp->GetEntries()==2)N12eq2=true;
            t2oR=t2oReset;		      
            pout->SavedCluHitIndex[pout->nCluster]=pout->nSavedHits;
            /*
               printf("numb clust=%d\n",pout->nCluster);
               for(int jcl=0;jcl<=pout->nCluster;jcl++)
               {
               printf("Cl[%d] hit index=%d  ",jcl,pout->SavedCluHitIndex[jcl]);
               };
               printf("\n");
               for(int jph=0;jph<pout->nPhotonClu;jph++)
               {
               printf("CLI[%d]=%d ",jph,pout->SavedPhotonCluIndex[jph]);
               };
               printf("\n");
               */
            LVec* tv1=(LVec*) ctmp->First();
            LVec* tv2=tv1;
            if(N12eq2)tv2=(LVec*) ctmp->After(tv1);
            Int_t npor1,npor2;

            //printf("tv1->E()=%f ,tv2->E()=%f \n ",tv1->E(),tv2->E());
            // printf("tv1->PhotOrder=%d  tv2->PhotOrder=%d \n",npor1=tv1->PhotOrder,npor2=tv2->PhotOrder);
            //		      printf("pout->nCluster =%d , pout->nPhotonClu=%d,pout->nSavedHits=%d\n",pout->nCluster, pout->nPhotonClu,pout->nSavedHits);
            //		      printf("SavedPhotonCluIndex[%d]=%d E=%f SavedPhotonCluIndex[%d]=%d E=%f \n",npor1,pout->SavedPhotonCluIndex[npor1],pout->SavedPhotonCluEnergy[npor1],npor2,pout->SavedPhotonCluIndex[npor2],pout->SavedPhotonCluEnergy[npor2]);  




            if(t2o.N12<=2 &&t2o.E12>35  )
            {
              LVec* tv1=(LVec*) ctmp->First();
              if(t2o.N12==2) LVec* tv2=(LVec*) ctmp->After(tv1);
              Int_t npor1,npor2;
              //			  pout->SavedCluHitIndex[nCluster]=pout->nSavedHits;
              npor1=npor2==tv1->PhotOrder;
              if(t2o.N12==2) npor2=tv2->PhotOrder;
              Int_t llim1=pout->SavedCluHitIndex[pout->SavedPhotonCluIndex[npor1]];
              Int_t ulim1=pout->SavedCluHitIndex[pout->SavedPhotonCluIndex[npor1]+1];
              Int_t llim2=pout->SavedCluHitIndex[pout->SavedPhotonCluIndex[npor2]];
              Int_t ulim2=pout->SavedCluHitIndex[pout->SavedPhotonCluIndex[npor2]+1];
              //			  printf("lim1 = (%d, %d) lim2=(%d,%d)\n",llim1,ulim1,llim2,ulim2);  

              Int_t  nstb=pout->WhatNSTB(tv1);
              TMatrix m=pout->FillFMSADC(nstb,llim1,ulim1);
              TMatrix mm=pout->FillFMSADC(nstb,llim2,ulim2);
              FitTower::SetNoCatag(true,1);
              t2oR.RCatag=1;
              if(llim2==llim1)
              {
                mm=0;
                FitTower::SetNoCatag(true,1);
                t2oR.RCatag=0;
              }; 
              m=m+mm;
              /*
                 Int_t  nstb=pout->WhatNSTB(tv1);
                 TMatrix m=pout->FillFMSADC(nstb);
                 */
              TMatrix me=m;
              Int_t row0=pout->WhatRow0(tv1,p_geom);
              Int_t col0=pout->WhatCol0(tv1,p_geom);
              for(int row=0;row<m.GetNrows();row++)
              {
                for(int col=0;col<m.GetNcols();col++)
                {
                  if(abs(row-row0)<40 &&abs(col-col0)<40)
                  {
                    (me[row][col]) =(m[row][col])*pgain->GetValue(2,nstb,row,col);
                    (me[row][col]) *=pgaincorr->GetValue(2,nstb,row,col);
                  }
                  else 
                  {
                    m[row][col]=0;
                    me[row][col]=0.;
                  }
                };
              };

              int pouttoievt=pout->ievt%10000;
              FitTower::SetYPrintLevel(0);      

              //	    printf("Evnum=%d \n",t2o.Evnum);
              //	    me.Print();

              if(t2o.Evnum==204)
              {
                printf("Evnum=%d \n",t2o.Evnum);
                me.Print();
                pgaincorr->tm(2,nstb)->Print();
                pgain->tm(2,nstb)->Print();
              };

              Float_t ECorRat[4]={0.,0.,0.,0.};
              ECorRat[nstb-1]=me.Sum();
              DAtaEv da((int) pout->Rnum,(int) pout->ievt/10000,(int) pout->ievt%10000);
              if(UseRunDepCor&&Mgr)
              {
                Mgr->RunDepBaseOverride=RunDepBaseOverride;
                Mgr->Rcorrect(&m,&me,nstb-1,2,&da);
                if(!(Mgr->Rdep))printf("Rcorrect called Mgr->Rdep undefined\n");
              }
              else
              {
                // printf(" No call to Rcorrect\n");
                if(UseRunDepCor&&(!Mgr)) printf("Mgr = 0\n");

              }
              if(ECorRat[nstb-1]!=0)ECorRat[nstb-1]=me.Sum()/ECorRat[nstb-1];
              t2oR.RCorRat=ECorRat[nstb-1];
              /*
                 printf("Evnum=%d ECorRat=%f %f %f %f ",
                 t2o.Evnum,ECorRat[0],ECorRat[1],ECorRat[2],
                 ECorRat[3]);
                 printf("AdcSum[%d]=%f Esum[%d]=%f \n",nstb-1,m.Sum(),nstb-1,me.Sum());
                 */

              if(t2o.Evnum==204)
              {
                printf("corrected\n");
                if(!Mgr)printf("Mgr Undefined\n");
                if(Mgr&&Mgr->Rdep)
                {printf("run=%d base=%d\n",(int) Mgr->Rdep->RunNumber,(int) Mgr->Rdep->BaseRunNumber);}
                else printf("Mgr->Rdep not defined\n");

                printf("Ph1=%d Ph2=%d \n",tv1->PhotOrder,tv2->PhotOrder);
                me.Print();
              };
              Yiqun* recon[2]={0,0};
              FitTower::SetSubclu2(false);
              if(t2o.Evnum==204) FitTower::SetYPrintLevel(2);
              FitTower::SetNoCatag(true,1);
              recon[0]=new Yiqun(&me,p_geom,pgain,pgaincorr,2,nstb);
              FitTower::SetNoCatag(true,2);

              recon[1]=new Yiqun(&me,p_geom,pgain,pgaincorr,2,nstb);
              FitTower::SetForceMass(-1);
              FitTower::SetNoCatag(false,2);
              FitTower::SetYPrintLevel(0);

              if(t2o.Evnum==204)
              {
                printf("found recon[0]->NPh=%d \n",recon[0]->NPh);
                for(int jpho=0;jpho<recon[0]->NPh;jpho++)
                {
                  printf("Event = %d photon for k:jpho=%d:%d",t2o.Evnum,nstb-1,jpho);
                  printf(" Energy0=%f \n",recon[0]->mom(jpho).E());
                };
                for(int jcc=0;jcc<recon[0]->NRealClusts;jcc++)
                {			  
                  printf("Event =%d Cluster for DEV=%d CLUST=%d \n",t2o.Evnum,nstb-1,jcc);
                  recon[0]->clust[jcc].Print();
                }
                printf("found recon[1]->NPh=%d \n",recon[1]->NPh);
                for(int jpho=0;jpho<recon[1]->NPh;jpho++)
                {
                  printf("Event = %d photon for k:jpho=%d:%d",t2o.Evnum,nstb-1,jpho);
                  printf(" Energy0=%f \n",recon[1]->mom(jpho).E());
                };
                for(int jcc=0;jcc<recon[1]->NRealClusts;jcc++)
                {			  
                  printf("Event =%d Cluster for DEV=%d CLUST=%d \n",t2o.Evnum,nstb-1,jcc);
                  recon[1]->clust[jcc].Print();
                }
              };

              t2oR.R1N12=recon[0]->NPh;
              t2oR.R2N12=recon[1]->NPh;
              // printf("NRealCluster=%d \n",rec.NRealClusts);
              if(recon[0]->NRealClusts>0 &&recon[0]->NPh<3&&recon[0]->clust[0].energy>35 )
              {
                TLorentzVector vLtmp(0,0,0,0);
                for(int jph=0;jph<recon[0]->clust[0].nPhoton;jph++)
                {
                  vLtmp+=recon[0]->mom(&(recon[0]->clust[0].photon[jph]));
                };
                t2oR.RN12C=recon[0]->clust[0].nPhoton;
                t2oR.R1M12=vLtmp.Mag();
                t2oR.R1E12=vLtmp.E();
                t2oR.R1Pt=vLtmp.Pt();
                t2oR.RChiSq1=recon[0]->clust[0].chiSquare;
                t2oR.RNDF1=recon[0]->clust[0].NDFchiSq;
                t2oR.R1SigMax=recon[0]->clust[0].sigmaMax;
                t2oR.R1SigMin=recon[0]->clust[0].sigmaMin;
              };
              if(recon[1]->NRealClusts>0 &&recon[1]->NPh<3&&recon[1]->clust[0].energy>35 )
              {

                TLorentzVector vLtmp(0,0,0,0);
                for(int jph=0;jph<recon[1]->clust[0].nPhoton;jph++)
                {
                  vLtmp+=recon[1]->mom(&(recon[1]->clust[0].photon[jph]));
                };

                t2oR.R2M12=vLtmp.Mag();
                t2oR.R2E12=vLtmp.E();
                t2oR.R2Pt=vLtmp.Pt();
                t2oR.RChiSq2=recon[1]->clust[0].chiSquare;
                t2oR.RNDF2=recon[1]->clust[0].NDFchiSq;
              };
              /*
                 if(rec.NRealClusts==2 &&rec.NPh==2 )
                 {
                 t2o.sigmaMaxF1=rec.clust[0].sigmaMax;
                 t2o.chiSquareF1=rec.clust[0].chiSquare;
                 t2o.M12=(rec.mom(0)+rec.mom(1)).Mag();
                 t2o.E12=(rec.mom(0)+rec.mom(1)).E();
                 t2o.Eta=(rec.mom(0)+rec.mom(1)).PseudoRapidity();
                 t2o.Phi=(rec.mom(0)+rec.mom(1)).Phi();
                 };
                 */
              if(recon[0])delete recon[0];
              if(recon[1])delete recon[1];

            }; 
          };


          if(ctmp->GetEntries()==1 &&Enable_ReRec &&false)
          {
            LVec* tv1=(LVec*) ctmp->First();
            if(tv1->E()>50) 
            {
              Int_t  nstb=pout->WhatNSTB(tv1);
              TMatrix m=pout->FillFMSADC(nstb);
              TMatrix me=m;

              Int_t row0=pout->WhatRow0(tv1,p_geom);
              Int_t col0=pout->WhatCol0(tv1,p_geom);
              for(int row=0;row<m.GetNrows();row++)
              {
                for(int col=0;col<m.GetNcols();col++)
                {
                  if(abs(row-row0)<300 &&abs(col-col0)<300)
                  {
                    (me[row][col]) *=pgain->GetValue(2,nstb,row,col);
                    (me[row][col]) *=pgaincorr->GetValue(2,nstb,row,col);
                  }
                  else me[row][col]=0.;
                };
              };

              FitTower::SetSubclu2(true);//Yiqun will substitute clu2 for cluster
              FitTower::SetNoCatag(true,1); 
              Yiqun rec(&me,p_geom,pgain,pgaincorr,2,nstb);
              FitTower::SetNoCatag(false);

              if(rec.NRealClusts==1)
              {
                t2o.sigmaMaxF1=rec.clust[0].sigmaMax;
                t2o.chiSquareF1=rec.clust[0].chiSquare;
              }
              FitTower::SetNoCatag(true,2);
              FitTower::SetForceMass(.135);
              FitTower::SetDoGlobal(false);
              Yiqun rec2(&me,p_geom,pgain,pgaincorr,2,nstb);
              FitTower::SetForceMass(-1);
              FitTower::SetNoCatag(false);
              FitTower::SetDoGlobal(true);
              if(rec2.NRealClusts==1)
              {
                t2o.chiSquareF2=rec2.clust[0].chiSquare;
              }
            };
          };
          //  printf("NAB=%d MAB=%f EAB=%f\n",NAB,MAB[NAB-1],EAB[NAB-1]);		  
          TwoTr.Fill();

        }; 
      };
      Int_t nphhard;
      if((nphhard=pout->scratchlist->GetEntries())==2 )
      {
        TLorentzVector twophot=pout->SumScratch();
        Float_t m2,e2;
        MassH.Fill(m2=twophot.Mag());
        M2vsPt.Fill(twophot.E(),twophot.Pt(),twophot.Mag());
        MvEta.Fill(m2,twophot.PseudoRapidity());
        if(e2>30)MvEta_E30.Fill(m2,twophot.PseudoRapidity());
        if(e2>40)MvEta_E40.Fill(m2,twophot.PseudoRapidity());
        if(e2>50)MvEta_E50.Fill(m2,twophot.PseudoRapidity());
        if(e2>60)MvEta_E60.Fill(m2,twophot.PseudoRapidity());
        if(e2>65)MvEta_E65.Fill(m2,twophot.PseudoRapidity());
        LVec* va=(LVec*) pout->vlist->First();//first photon
        LVec* vb=(LVec*) pout->vlist->After(va);//other photon
        va->SetPartner(vb);
        vb->Partner=va;

        Float_t Zz=(*va+*vb).E();
        if(Zz>0){Zz=fabs((va->E()-vb->E())/Zz);} else {Zz=1;};


        if(fabs(m2-.135)<.08)
        {
          LVec ltwo(5,twophot);//fake an LVec to get some z position (ie NL
          TVector3 vpi=pout->PosInDet(&ltwo,p_geom);
          XYpi.Fill(vpi.X(),vpi.Y()); //histogram position for 2 phot pi

          pout->histpair(va,p_geom,0.,0);//fill histogram for 2 photon locations
          MPiEtavE.Fill((*va+*vb).PseudoRapidity(),(*va+*vb).E());

          if((*va+*vb).E()>30)MPiEtavZ_E30.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>40)MPiEtavZ_E40.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>50)MPiEtavZ_E50.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>60)MPiEtavZ_E60.Fill((*va+*vb).PseudoRapidity(),Zz);

        };

        if(fabs(m2-.55)<.25)
        {
          LVec ltwo(5,twophot);
          TVector3 vpi=pout->PosInDet(&ltwo,p_geom);
          XYEta.Fill(vpi.X(),vpi.Y());//histogram for eta (high) mass photons
          if((*va+*vb).E()>30)MEtaEtavZ_E30.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>40)MEtaEtavZ_E40.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>50)MEtaEtavZ_E50.Fill((*va+*vb).PseudoRapidity(),Zz);
          if((*va+*vb).E()>60)MEtaEtavZ_E60.Fill((*va+*vb).PseudoRapidity(),Zz);

        };
        E2phH.Fill(e2=twophot.E());
        EDetH.Fill(pout->TotaldetE);
        Mass2vsE2.Fill(m2,e2);
        if(fabs(twophot.PseudoRapidity()-3.55)<.25) MvE_Eta_33_38.Fill(m2,e2);		
        if(fabs(twophot.PseudoRapidity()-4.05)<.25) MvE_Eta_38_43.Fill(m2,e2);		
      };
      Int_t n1=pout->GetNPhotVec(2,1);
      Int_t n2=pout->GetNPhotVec(2,2);
      Int_t n3=pout->GetNPhotVec(2,3);
      Int_t n4=pout->GetNPhotVec(2,4);
      if((nphhard=pout->scratchlist->GetEntries())==3 )
      {
        t3o.spin=pout->spin;
        t3o.BBcVertex=pout->BBcVertex;

        LVec* lphi=pout->GetPairNearMass(pout->scratchlist,4.);
        Mhi3.Fill(lphi->PairMass(),lphi->E());
        Mhi3Pt.Fill(lphi->PairMass(),lphi->Pt());

        LVec* lpho=pout->GetPairNearMass(pout->scratchlist,3.);
        Mho3.Fill(lpho->PairMass(),lpho->E());	
        LVec lp=*pout->GetPairNearMass(pout->scratchlist,.135);
        LVec* lpartner=lp.Partner;
        LVec lpEta=*pout->GetPairNearMass(pout->scratchlist,.55);
        TLorentzVector v3=pout->SumScratch();

        M3BPairM.Fill(lp.PairMass());
        M3BPairEM.Fill(lp.PairMass(),v3.E());
        M3BEtaEM.Fill(lpEta.PairMass(),v3.E());
        if(fabs(lp.PairMass()-.135)>.1)
        {
          M3BEtaEMNPi.Fill(lpEta.PairMass(),v3.E());
        };
        TLorentzVector veta1=lpEta;
        TLorentzVector veta2=*lpEta.Partner;
        TLorentzVector v00=v3-veta1-veta2;

        LVec* v3_1=pout->GetPairNearMass(pout->scratchlist,.0);
        TLorentzVector vlo=pout->Pair4V(v3_1);
        TLorentzVector* v3_1a=v3_1;
        TLorentzVector* v3_1b=v3_1->Partner;
        if((v3_1a->E()) < (v3_1b->E()))
        {
          v3_1b=v3_1;
          v3_1a=v3_1->Partner;
          if(fabs(vlo.E()-v3_1a->E()-v3_1b->E())>.001)
          {
            printf("!!!!!!!!!!!!!! inconsistency of 3 vectors\n");
          };
        };


        LVec* v3_3=pout->GetPairNearMass(pout->scratchlist,1000.);
        TLorentzVector vhi=pout->Pair4V(v3_3);
        LVec* v3_2=pout->GetPairNearMass(pout->scratchlist,
            (vlo.Mag()+vhi.Mag())/2.);	
        TLorentzVector vmid=pout->Pair4V(v3_2);

        Mdal1.Fill(pow((vlo).Mag(),2),pow((vmid).Mag(),2),
            (vhi).Mag());
        Mdal2.Fill(pow((vlo).Mag(),2),pow((vmid).Mag(),2),
            (v3).E());
        TLorentzVector vz(0,0,1,1);
        TVector3 tvv=pout->reframe(vlo,v3,vz);
        t3o.dphilo=tvv.Phi();
        t3o.dthetalo=tvv.Theta();
        t3o.Ptlo=vlo.Pt();
        t3o.elo=vlo.E();
        t3o.mlo=vlo.Mag();
        t3o.Ylo=vlo.PseudoRapidity();
        t3o.Philo=vlo.Phi();

        t3o.Eloa=v3_1a->E();
        t3o.Phloa=v3_1a->Phi();
        t3o.Yloa=v3_1a->PseudoRapidity();
        t3o.Ptloa=v3_1a->Pt();

        TVector3 tvv0=pout->reframe(vmid,v3,vz);
        t3o.dphimid=tvv0.Phi();
        t3o.dthetamid=tvv0.Theta();
        t3o.Ptmid=vmid.Pt();
        t3o.emid=vmid.E();
        t3o.mmid=vmid.Mag();
        t3o.Ymid=vmid.PseudoRapidity();
        t3o.Phimid=vmid.Phi();

        t3o.mpi=lp.PairMass();
        t3o.meta=lpEta.PairMass();
        t3o.mhi=vhi.Mag();
        t3o.mmid=vmid.Mag();

        t3o.e=v3.E();
        t3o.m=v3.Mag();
        t3o.Phi=v3.Phi();
        t3o.Y=v3.PseudoRapidity();
        t3o.Pt=v3.Pt();

        t3o.Pthi=vhi.Pt();
        t3o.ehi=vhi.E();
        t3o.mhi=vhi.Mag();
        t3o.Yhi=vhi.PseudoRapidity();
        t3o.Phihi=vhi.Phi();
        t3o.spin=pout->spin;

        if((fabs(lpEta.PairMass()-.55)<.12) &&(fabs(lp.PairMass()-.135)>.1) )
        {
          if(v00.E()>(veta1+veta2).E())
          {
            M3Eta.Fill(v3.E(),v3.Pt(),v3.Mag());
          }
          else
          {
            M3EtaL.Fill(v3.E(),v3.Pt(),v3.Mag());
          };
        };
        if((fabs(t3o.mlo -.135) < .08))
        {

          if(t3o.e-t3o.elo>t3o.elo)
          {
            Momega.Fill(v3.E(),v3.Pt(),v3.Mag());
          }
          else
          {
            MomegaL.Fill(v3.E(),v3.Pt(),v3.Mag());		
          };
          MomegEPh.Fill(v3.Mag(),v3.E(),t3o.dphilo);
        }
        else
        {
          MomegEPhN.Fill(v3.Mag(),v3.E(),t3o.dphilo);
          MomegaN.Fill(v3.E(),v3.Pt(),v3.Mag());
        };

        //ThreeTr.Fill();    

      };
      pout->AllToScratch(true);//full hard and soft to scratch
      Int_t nhardsoft=pout->scratchlist->GetEntries();

      if(nhardsoft==3){
        LVec* lp=pout->GetPairNearMass(pout->scratchlist,.135);
        NearPi.Fill( lp->PairMass());
        Float_t EbPair=lp->E();
        EbPair+=(lp->Partner)->E();
        HSpairMvsE.Fill(lp->PairMass(),EbPair); 

        if(fabs(lp->PairMass() -.135) < .08)
        {
          pout->histpair(lp,p_geom,0.,1);
        };
        for(Float_t m0=.025;m0<.75;m0+=.02)
        {
          LVec* lp0=pout->GetPairNearMass(pout->scratchlist,m0);
          HSpairMvsm0.Fill(lp0->PairMass(),m0);
        };
      };
      if(nhardsoft==2 && (n3+n4)==1)
      {
        Float_t mtmp;
        MassHS.Fill(mtmp=(pout->SumScratch().Mag()));
        LVec* va=(LVec*) pout->scratchlist->First();//first photon
        LVec* vb=(LVec*) pout->scratchlist->After(va);//other photon
        Float_t E2=va->E()+vb->E();
        Float_t z=0;
        z=fabs((va->E()- vb->E())/E2 );
        MassHSvsZ.Fill(mtmp,z);
        if(fabs(mtmp-.135)>.08 || z>.8)continue;
        TVector3 vhit=pout->PosInDet(va,p_geom);
        XYpiHS.Fill(vhit.X(),vhit.Y());
        vhit=pout->PosInDet(vb,p_geom);
        XYpiHS.Fill(vhit.X(),vhit.Y());
      };
      if(pout->Nphotons==4)
      {
        pout->AllToScratch();
        if(pout->scratchlist->GetEntries()!=4)printf("scratchlist error\n");
        Float_t Mass4=pout->SumScratch().Mag();
        Float_t E4=pout->SumScratch().E();
        LVec* pair=pout->GetPairNearMass(pout->scratchlist,.135);
        Float_t M4_1=pair->PairMass();
        pout->RemoveFromScratch(pair);
        pout->RemoveFromScratch(pair->Partner);
        if(pout->scratchlist->GetEntries()!=2)printf("scratchlist error\n");
        Float_t M4_2=pout->SumScratch().Mag();
        Near2Pi.Fill(M4_1,M4_2);
        if(fabs(M4_2-.55)<.15 && fabs(M4_1-.135)<.04)
        {
          pout->AllToScratch();
          M2PivE.Fill(Mass4,E4);
          TIter next(pout->scratchlist);

          while(LVec* vv=(LVec*) next())
          {
            TVector3 vhit=pout->PosInDet(vv,p_geom);
            XY4Eta.Fill(vhit.X(),vhit.Y(),vv->E());
          };

        };
      };
      if(n3==2)
      {
        TLorentzVector v13(0,0,0,0);
        //	for(Int_t i=0;i<n1;i++)v13+=pout->GetPhotVec(2,1,i);
        for(Int_t i=0;i<n3;i++)v13+=pout->GetPhotVec(2,3,i);
        MvEtwo_13.Fill(v13.Mag(),v13.E());
        EtavEtwo_13.Fill(v13.PseudoRapidity(),v13.E());
        if(v13.E()>40)EtavMtwo_13.Fill(v13.PseudoRapidity(),v13.Mag());
      };
      if(n4==2)
      {
        TLorentzVector v24(0,0,0,0);
        //for(Int_t i=0;i<n2;i++)v24+=pout->GetPhotVec(2,2,i);
        for(Int_t i=0;i<n4;i++)v24+=pout->GetPhotVec(2,4,i);
        MvEtwo_24.Fill(v24.Mag(),v24.E());
        EtavEtwo_24.Fill(v24.PseudoRapidity(),v24.E());
        if(v24.E()>40)EtavMtwo_24.Fill(v24.PseudoRapidity(),v24.Mag());

      };
      if(ClustL)delete ClustL;
    };
  };
  out->Print();
  for(Int_t j0=0;j0<4;j0++)
  {
    Ecomp1[j0].Write();
    Ecomp2[j0].Write();
    Ecomp[j0].Write();
  };
  cl2MvE.Write();
  TwoC.Write();
  TwoTr.Write();
  ThreeTr.Write();
  M3Eta.Write();
  M3EtaL.Write();
  M3BEtaEM.Write();
  M3BEtaEMNPi.Write();
  MomegEPh.Write();
  MomegEPhN.Write();
  M2vsPt.Write();
  M3BPairM.Write();
  M3BPairEM.Write();
  Mdal1.Write();
  Mdal2.Write();
  Momega.Write();
  MomegaL.Write();
  MomegaN.Write();
  Mhi3.Write();
  Mhi3Pt.Write();
  Mho3.Write();
  MassH.Write();
  MassH0.Write();
  MassHS.Write();
  MassHSvsZ.Write();
  E2phH.Write();
  EDetH.Write();
  Mass2vsE2.Write();
  MvEtwo_13.Write();
  NphH.Write();
  NphHAll.Write();
  EtavEtwo_13.Write();
  EtavMtwo_13.Write();
  EtavEtwo_24.Write();
  EtavMtwo_24.Write();
  NearPi.Write();
  Near2Pi.Write();
  M2PivE.Write();
  MvE_Eta_38_43.Write();
  MvE_Eta_33_38.Write();
  XYpi.Write();
  XYpiHS.Write();
  XYsoft.Write();
  XYEta.Write();
  XY4Eta.Write();
  XY.Write();
  MvEta.Write();
  MvEta_E30.Write();
  MvEta_E40.Write();
  MvEta_E50.Write();
  MvEta_E60.Write();
  MvEta_E65.Write();
  HSpairMvsE.Write();
  MPiEtavE.Write();
  MPiEtavZ_E30.Write();
  MPiEtavZ_E40.Write();
  MPiEtavZ_E50.Write();
  MPiEtavZ_E60.Write();

  MEtaEtavZ_E30.Write();
  MEtaEtavZ_E40.Write();
  MEtaEtavZ_E50.Write();
  MEtaEtavZ_E60.Write();

  pout->p_fms[0][0]->Write();
  pout->p_fms[1][0]->Write();
  pout->p_fms[2][0]->Write();
  pout->p_fms[3][0]->Write();
  pout->p_fms[0][1]->Write();
  pout->p_fms[1][1]->Write();
  pout->p_fms[2][1]->Write();
  pout->p_fms[3][1]->Write();
  HSpairMvsm0.Write();
  delete p_geom;
  return 0;
};
Int_t AnalTools::readSimple(char* filelist,Int_t set)
{
  Int_t DetMap[9]={0,1,2,0,0,5,6,3,4};

  //**********Get Access Calib files
  FilesSet* p_files=0;
  Int_t spinindex=0;
  if(set>=61 && set <=63)
  {
    //transverse runs
    p_files=new FilesSet("../",
        "run6_export/fpdped.txt",
        "run6_export/fpdgain.txt", 
        "root06T/fpdcorr.txt",
        "run6_export/fill.txt",
        "Fake",
        "run6_export/spinpat",
        "root06T/geom.txt");
  };

  if(set==80)
  {
    //transverse run 8
    p_files=new FilesSet("/star/u/heppel/BatchDir/root08",
        "fpdped.txt",
        "fpdgain.txt", 
        "fpdcorr.txt",
        "fill.txt",
        "Fake",
        "spinpat",
        "geom.txt");
  }
  //****** Reading in the chain and perhaps joining into 
  //****** one file (OutputFileName)
  p_Geom=new Geom(p_files);

  TChain* p_OF=new TChain("p_out");
  FILE* pf;
  Int_t knt=0;
  Int_t FileCnt=0;

  if(pf=fopen(filelist,"r"))
  {
    char fname[200];
    while(!feof(pf))
    {
      if(fscanf(pf,"%s\n",fname)>0)
      {
        knt++;
        FileCnt++;
        p_OF->Add(fname);
      };
      if(knt>20)
      {
        knt=0;
        printf("%d files added to TChain\n",FileCnt);
      };
      // if(FileCnt>20)break;// max number of segs for testing
    };
  };
  if(OutputToSingle) {
    std::cout<<"Output to "<<OutputFileName<<"\n";
    p_OF->Merge( OutputFileName);
    return 0;
  };

  std::cout<<"OuputToSingle not selected\n";

  // **********end of  Reading in the chain and perhaps ...


  // ***********Now set up the loop for reading tree "p_OF"
  GetRunIndex(0,true);// reset run index
  Int_t pcnt=0;
  TVector3 v3;
  TLorentzVector v4;
  TH2F  Nphots("Nphots","Nphots",30,0,30,30,0,30);

  Int_t tpes[320];
  Int_t spin;
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t nphotons;
  Int_t Rnum;
  Int_t FillNumber=0;
  Int_t LastRun=0;
  Int_t Bunchid7bit;
  Float_t raw_adc[98];

  UChar_t BBcSums[5];
  p_OF->SetBranchAddress("spin",&spin);
  p_OF->SetBranchAddress("nphotons",&nphotons);
  p_OF->SetBranchAddress("br_nwrds",&nwrds);
  p_OF->SetBranchAddress("br_types",tpes);
  p_OF->SetBranchAddress("br_pxyzt",pxyzt);
  p_OF->SetBranchAddress("br_Rnum",&(Rnum));
  p_OF->SetBranchAddress("br_Bunchid7bit",&(Bunchid7bit));
  p_OF->SetBranchAddress("br_BBcSums",BBcSums);
  if(p_OF->FindBranch("br_adc"))
  {  
    p_OF->SetBranchAddress("br_adc",raw_adc);
  };

  TLorentzVector vec[80];
  TH1F MassEN("MassEN","MassEN",100,0,1.);  
  TH2F MassENvsRun("MassENvsRun","MassENvsRun",1000,0,1000,100,0,1.);  
  TH1F MassES("MassES","MassES",100,0,1.);
  TH2F MassESvsRun("MassESvsRun","MassESvsRun",1000,0,1000,100,0,1.);  

  TH1F* MhArray[9];  MhArray[1]=&MassEN;  MhArray[2]=&MassES;

  TH1F PTYup("PTYup","PTYup",200,-10,10.);
  TH1F PTYdn("PTYdn","PTYdn",200,-10,10.);

  TH2F PhotonENCells("PhotonENCells","PhotonENCells",7,0.,7.,7,0.,7.);
  TH2F PhotonESCells("PhotonESCells","PhotonESCells",7,0.,7.,7,0.,7.);
  TH2F P2ENCells("P2ENCells","P2ENCells",7,0.,7.,7,0.,7.);
  TH2F P2ESCells("P2ESCells","P2ESCells",7,0.,7.,7,0.,7.);
  TH2F PCENCells("PCENCells","PCENCells",7,0.,7.,7,0.,7.);
  TH2F PCESCells("PCESCells","PCESCells",7,0.,7.,7,0.,7.);
  TH2F PUENCells("PUENCells","PUENCells",7,0.,7.,7,0.,7.);
  TH2F PUESCells("PUESCells","PUESCells",7,0.,7.,7,0.,7.);

  TH2F PCENcCells("PCENcCells","PCENcCells",7,0.,7.,7,0.,7.);
  TH2F PCEScCells("PCEScCells","PCEScCells",7,0.,7.,7,0.,7.);
  TH2F PUENcCells("PUENcCells","PUENcCells",7,0.,7.,7,0.,7.);
  TH2F PUEScCells("PUEScCells","PUEScCells",7,0.,7.,7,0.,7.);


  TH2F* PhArray[9]; PhArray[1]=&PhotonENCells;PhArray[2]=&PhotonESCells;

  TH1F EdistN("EdistN","EdistN",100,0.,100.);
  TH1F EdistS("EdistS","EdistS",100,0.,100.);
  TH1F* EhArray[9];EhArray[1]=&EdistN;EhArray[2]=&EdistS; 

  TH3F PtEnSp("PtEnSp","PtEnSp",50,-5,5,20,0,100.,5,0,5);
  Float_t lowlim=0;
  Float_t highlim=120.;
  TH2F absPtPi("absPtPi","absPtPi",24,1.0,4.0,30,0,30);
  TH2F Ex("Ex","Ex",24,lowlim,highlim,30,0,30);
  TH2F ExPi("ExPi","ExPi",24,lowlim,highlim,30,0,30);
  TH2F ExPiPt2("ExPiPt2","ExPiPt2",24,lowlim,highlim,30,0,30);
  TH2F ExPiPt15("ExPiPt15","ExPiPt15",24,lowlim,highlim,30,0,30);
  TH2F ExEta("ExEta","ExEta",24,lowlim,highlim,30,0,30);
  TH2F ExBetw("ExBetw","ExBetw",24,lowlim,highlim,30,0,30);
  TH2F ENEs("ENEs","ENEs",100,0,100,100,0,100);
  TH1F ADCnor("ADCnor","ADCnor",256,0,256);
  TH1F ADCsou("ADCsou","ADCsou",256,0,256);
  TH2F ADCNS("ADCNS","ADCNS",256,0,256,256,0,256);
  TH2F ADCsumNS("ADCsumNS","ADCsumNS",600,0,600,600,0,600);
  TH2F ADCsumNS2("ADCsumNS2","ADCsumNS2",600,0,600,600,0,600);
  TH2F* sphxy[2];
  sphxy[0]=new TH2F("sphxyN","sphxyN",70,0,7.,70,0,7.);
  sphxy[1]=new TH2F("sphxyS","sphxyS",70,0,7.,70,0,7.);
  TH2F* ptvsx[2];
  ptvsx[0]=new TH2F("ptvsxN","ptvsxN",50,0,1.,20,0.,4.);
  ptvsx[1]=new TH2F("ptvsxS","ptvsxS",50,0,1.,20,0.,4.);
  TH1F* EdistXY[7][7][4];
  TH1F* RawAdc[7][7][4];
  TH1F* Mcell[7][7][2];
  TH1F* Mcell2[7][7][2];
  TH2F* McellZ[7][7][2];
  TH2F* EphN=new TH2F("EphN","EphN",49,0,49,120,0,120);
  TH2F* EphS=new TH2F("EphS","EphS",49,0,49,120,0,120);
  TH2F* EpiN=new TH2F("EpiN","EpiN",49,0,49,120,0,120);
  TH2F* EpiS=new TH2F("EpiS","EpiS",49,0,49,120,0,120);
  TH2F* EpiN1=new TH2F("EpiN1","EpiN1",49,0,49,120,0,120);
  TH2F* EpiS1=new TH2F("EpiS1","EpiS1",49,0,49,120,0,120);
  TH2F* EvsZ[7][7][2];
  for(Int_t ix=0;ix<7;ix++)
  {
    for(Int_t iy=0;iy<7;iy++)
    {
      char nam[20];
      sprintf(nam,"EphXY%d_%d_N",ix+1,iy+1);
      EdistXY[ix][iy][0]=new TH1F(nam,nam,100,0,100);
      sprintf(nam,"EphXY%d_%d_S",ix+1,iy+1);
      EdistXY[ix][iy][1]=new TH1F(nam,nam,100,0,100);
      sprintf(nam,"RawAdc%d_%d_N",ix+1,iy+1);
      RawAdc[ix][iy][0]=new TH1F(nam,nam,256,0,256);
      sprintf(nam,"RawAdc%d_%d_S",ix+1,iy+1);
      RawAdc[ix][iy][1]=new TH1F(nam,nam,256,0,256);
      sprintf(nam,"RawAdcL%d_%d_N",ix+1,iy+1);
      RawAdc[ix][iy][2]=new TH1F(nam,nam,256,0,256);
      sprintf(nam,"RawAdcL%d_%d_S",ix+1,iy+1);
      RawAdc[ix][iy][3]=new TH1F(nam,nam,256,0,256);
      sprintf(nam,"Mcell%d_%d_N",ix+1,iy+1);
      Mcell[ix][iy][0]=new TH1F(nam,nam,100,0,1.);
      sprintf(nam,"Mcell%d_%d_S",ix+1,iy+1);
      Mcell[ix][iy][1]=new TH1F(nam,nam,100,0,1.);

      sprintf(nam,"Mcell2%d_%d_N",ix+1,iy+1);
      Mcell2[ix][iy][0]=new TH1F(nam,nam,100,0,1.);
      sprintf(nam,"Mcell2%d_%d_S",ix+1,iy+1);
      Mcell2[ix][iy][1]=new TH1F(nam,nam,100,0,1.);

      sprintf(nam,"McellZ%d_%d_N",ix+1,iy+1);
      McellZ[ix][iy][0]=new TH2F(nam,nam,60,0,.3,6,30.,60.);
      sprintf(nam,"McellZ%d_%d_S",ix+1,iy+1);
      McellZ[ix][iy][1]=new TH2F(nam,nam,60,0,.3,6,30.,60.);

      sprintf(nam,"EvsZ%d_%d_N",ix+1,iy+1);
      EvsZ[ix][iy][0]=new TH2F(nam,nam,10,0.,1.,10,0.,60.);
      sprintf(nam,"EvsZ%d_%d_S",ix+1,iy+1);
      EvsZ[ix][iy][1]=new TH2F(nam,nam,10,0.,1.,10,0.,60.);

    };
  };

  Int_t nentries=p_OF->GetEntries();
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;
  printf("GetEntries=%d \n",nentries);
  knt=0;
  Bool_t YellowUp,YellowDn,BlueUp,BlueDn;
  Bool_t starting=true;

  //********************Now Start the loop

  for(Int_t i=0;i<nentries;i++)
  {
    YellowUp=YellowDn=BlueUp=BlueDn=false;
    Int_t nbytes=p_OF->GetEntry(i);
    if(spin==1 || spin==3)YellowUp=true;
    if(spin==0 || spin==2)YellowDn=true;
    if(spin==2 || spin==3)BlueUp=true;
    if(spin==0 || spin==1)BlueDn=true;

    knt++;

    if(knt>20000 || (starting && knt>10000)){
      printf("Evt: %d \n",i);
      knt=1;
      if(i>100000)starting=false;
    };

    Int_t kk=0;
    TLorentzVector ph[10];
    Int_t nlv[10];
    Int_t plv[10][10];
    Int_t nph2[10];
    Int_t pph2[10][10];
    Int_t nph3[10];
    Int_t pph3[10][10];
    Float_t hiEsum[10],loEsum[10];

    for(Int_t i1=0;i1<10;i1++)
    {
      nlv[i1]=0;
      nph2[i1]=0;
      nph3[i1]=0;
      hiEsum[i1]=0;
      loEsum[i1]=0;
    };
    Int_t vtyp;
    Int_t nwords=nwrds;

    //*****************Now fill the vec array with 4 vectors

    for(Int_t k=0;k<nwords;k+=4){
      Bool_t cut56=true;
      if(tpes[k]<9)
      {
        vtyp=tpes[k];
        plv[vtyp][nlv[vtyp]]=k;
        vec[k].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
        nlv[vtyp]++;
      }
      else if((tpes[k]>300) && (tpes[k]<309))
      {
        vtyp=tpes[k]-300;
        cut56=true;
        pph2[vtyp][nph2[vtyp]]=k;
        pph3[vtyp][nph3[vtyp]]=k;
        vec[k].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);  
        nph2[vtyp]++;  
      };
    };

    // now loop over all possible detectors even though we will only look
    // at  k=1 which is EN
    Float_t Enor=0;
    Float_t Esou=0;
    if(nlv[1]>0)Enor=vec[plv[1][0]].E();
    if(nlv[2]>0)Esou=vec[plv[2][0]].E();
    ENEs.Fill(Enor*1.,Esou*1.);
    Float_t NADCsum=0.;
    Float_t SADCsum=0.;
    for(Int_t j=0;j<49;j++)
    {
      NADCsum+=raw_adc[j];
      SADCsum+=raw_adc[j+49];
    };

    for(Int_t j=0;j<49;j++)
    {
      Int_t irow=j/7;
      Int_t icol=j-7*irow;
      Int_t rowl=irow-1;
      if(rowl<0)rowl=0;
      Int_t rowh=irow+1;
      if(rowh>6)rowh=6;
      Int_t coll=icol-1;
      if(coll<0)coll=0;
      Int_t colh=icol+1;
      if(colh>6)colh=6;
      Float_t corv;
      for(Int_t ir=rowl;ir<=rowh;ir++)
      {
        for(Int_t ic=coll;ic<=colh;ic++)
        {
          Int_t cntr=7*ir+ic;
          if(ic !=icol ||ir != irow)
          {
            if(raw_adc[j]>20){
              PCENCells.Fill(ic+.5,ir+.5,raw_adc[cntr]*raw_adc[j]);
              if(raw_adc[cntr]>10)PCENcCells.Fill(ic+.5,ir+.5);
            };
            if(raw_adc[j+49]>20){
              PCESCells.Fill(ic+.5,ir+.5,raw_adc[cntr+49]*raw_adc[j+49]);
              if(raw_adc[cntr+49]>10)PCEScCells.Fill(ic+.5,ir+.5);
            };
          }
          else
          {
            if(raw_adc[j]>20){
              PUENCells.Fill(ic+.5,ir+.5,raw_adc[cntr]*raw_adc[j]);
              PUENcCells.Fill(ic+.5,ir+.5);
            };
            if(raw_adc[j+49]>20){
              PUESCells.Fill(ic+.5,ir+.5,raw_adc[cntr+49]*raw_adc[j+49]);
              PUEScCells.Fill(ic+.5,ir+.5);

            };

          };
        };
      };
      ADCnor.Fill(raw_adc[j]);
      ADCsou.Fill(raw_adc[j+49]);
      ADCNS.Fill(raw_adc[j],raw_adc[j+49]);
      if(NADCsum>124)RawAdc[icol][irow][0]->Fill(raw_adc[j]);
      if(SADCsum>124)RawAdc[icol][irow][1]->Fill(raw_adc[j+49]);
      if(SADCsum>124.)RawAdc[icol][irow][2]->Fill(raw_adc[j]);
      if(NADCsum>124.)RawAdc[icol][irow][3]->Fill(raw_adc[j+49]);
    };
    ADCsumNS.Fill(NADCsum,SADCsum);
    ADCsumNS2.Fill(NADCsum-raw_adc[34],SADCsum);
    Bool_t goodtrig=true;
    if(NADCsum<120 && SADCsum<120)goodtrig=false;
    if(!goodtrig)
    {
      for(Int_t j=0;j<49;j++)
      {
        Int_t irow=j/7;
        Int_t icol=j-7*irow;
        P2ENCells.Fill(icol+.5,irow+.5,raw_adc[j]);
        P2ESCells.Fill(icol+.5,irow+.5,raw_adc[j+49]);
      };
    };
    for(Int_t k=0;k<9;k++)
    {
      Int_t tspn,spinindex;
      if(spin>=0 && spin<4)tspn=spin;
      spinindex=DetMap[k]*4+tspn;

      if(k==1 || k==2) // EN or ES FPD
      {

        Float_t zz=0;
        Int_t kbar=1;
        if(k==1)kbar=2;
        if(nph2[k]==1)
        { 
          TLorentzVector phot=vec[pph2[k][0]];
          TVector3 Sa=phot.Vect();
          Sa=Sa*(1./Sa.Z()*(*p_Geom->ZFPD(1,k)));
          TVector3 va=p_Geom->LocalXYZ(1,k,Sa,true);
          if(phot.E()>45.)sphxy[k-1]->Fill(va.X(),va.Y());
          ptvsx[k-1]->Fill(phot.E()/100.,phot.Pt());
          Float_t e1=phot.E();
          Int_t cnum=(Int_t) va.X()+7*((Int_t) va.Y());

          if(k==1)EpiN1->Fill(cnum,(e1)*EdepFactor(e1));
          if(k==2)EpiS1->Fill(cnum,(e1)*EdepFactor(e1));



        };
        Float_t ESUM=0;
        if(nph2[k]>0)
        {
          ESUM=vec[plv[k][0]].E();

          Int_t ix=-1;Int_t iy=-1;Float_t photMaxE=0.;
          Bool_t skip=true;
          //	  for(Int_t jp=0;jp<nph2[k];jp++)
          if(nph2[k]==2)
          {

            if(fabs((vec[pph2[k][0]]+vec[pph2[k][1]]).Mag()-.135)>.5)continue;
            zz=fabs((vec[pph2[k][0]].E()-vec[pph2[k][1]].E())/(vec[pph2[k][0]].E()+vec[pph2[k][1]].E()));
          };
          for(Int_t jp=0;jp<nph2[k];jp++)
          {

            TLorentzVector phot=vec[pph2[k][jp]];
            if(ix==-1)skip=true;
            if( phot.E()>photMaxE)
            {
              photMaxE=phot.E();
              TVector3 Sa=phot.Vect();
              Sa=Sa*(1./Sa.Z()*(*p_Geom->ZFPD(1,k)));
              TVector3 va=p_Geom->LocalXYZ(1,k,Sa,true);
              ix=(Int_t) va.X();
              iy=(Int_t) va.Y();
              skip=false;
              if(ix<0)skip=true;
              if(ix>6)skip=true;
              if(iy<0)skip=true;
              if(iy>6)skip=true;
            };
          };

          if((!skip) && photMaxE>0 &&ESUM>40)
          {
            Int_t cnum=ix+7*iy;
            if(k==1)EphN->Fill(cnum,photMaxE);
            if(k==2)EphS->Fill(cnum,photMaxE);
            EdistXY[ix][iy][k-1]->Fill(photMaxE);
          };
        };

        if(nlv[k]>0)
        {
          TLorentzVector det=vec[plv[k][0]];
          Float_t EneN=det.E();
          EhArray[k]->Fill(EneN);
        };
        if(nph2[k]==2) // 2 Photon Event
        {
          TLorentzVector vphot1=vec[pph2[k][0]];
          TLorentzVector vphot2=vec[pph2[k][1]];
          TLorentzVector vpi=vphot1+vphot2;
          MhArray[k]->Fill(vpi.Mag());
          if(k==1)MassENvsRun.Fill(GetRunIndex(Rnum),vpi.Mag());
          if(k==2)MassESvsRun.Fill(GetRunIndex(Rnum),vpi.Mag());


          // for fun, find 3 vectors locations 
          // for photons in frame of detector with units of 
          // detector cell size (vaEN & vbEN)

          TVector3 Sa=vphot1.Vect();
          TVector3 Sb=vphot2.Vect();
          TVector3 Sc=vpi.Vect();
          Sa=Sa*(1./Sa.Z()*(*p_Geom->ZFPD(1,k)));
          Sb=Sb*(1./Sb.Z()*(*p_Geom->ZFPD(1,k)));
          Sc=Sc*(1./Sc.Z()*(*p_Geom->ZFPD(1,k)));
          TVector3 va=p_Geom->LocalXYZ(1,k,Sa,true);
          TVector3 vb=p_Geom->LocalXYZ(1,k,Sb,true);
          TVector3 vc=p_Geom->LocalXYZ(1,k,Sc,true);
          // now accum photon energy in cells

          PhArray[k]->Fill(va.X(),va.Y(),vphot1.E());
          PhArray[k]->Fill(vb.X(),vb.Y(),vphot2.E());
          Bool_t Inside=true;
          Int_t iax=(Int_t) va.X();
          Int_t iay=(Int_t) va.Y();
          Int_t ibx=(Int_t) vb.X();
          Int_t iby=(Int_t) vb.Y();
          Int_t icx=(Int_t) vc.X();
          Int_t icy=(Int_t) vc.Y();
          if(va.X()<.5 || va.X()>6.5)Inside=false;
          if(va.Y()<.5 || va.Y()>6.5)Inside=false;
          if(vb.X()<.5 || vb.X()>6.5)Inside=false;
          if(vb.Y()<.5 || vb.Y()>6.5)Inside=false;
          if(vc.X()<.5 || vc.X()>6.5)Inside=false;
          if(vc.Y()<.5 || vc.Y()>6.5)Inside=false;
          Float_t e1=vphot1.E();
          Float_t e2=vphot2.E();
          Float_t Etot=e1+e2;
          Float_t SumRat=ESUM/Etot;
          //		  SumRat=1.;

          Float_t z=0;
          if(e1+e2>0)z=fabs((e1-e2)/(e1+e2));
          //		  if(z<.85 && Inside && (e1+e2)>30)
          //		  if(z<.85 && Inside )
          if(z<.95 )
          {

            if(icx>=0 && icx<7 && icy>=0 && icy<7&& (e1+e2)>40)
            {
              if(iax>=0 && iax<=7 && iay>=0 && iay<=7 &&
                  ibx>=0 && ibx<=7 && iby>=0 && iby<=7)
              {
                EvsZ[icx][icy][k-1]->Fill(z,Etot*EdepFactor(Etot));
                Int_t cnum=icx+7*icy;
                if(fabs(vpi.Mag()*EdepFactor(e1+e2)-.135)<.1)
                {
                  if(k==1)
                  {
                    EpiN->Fill(cnum,(e1+e2)*EdepFactor(e1+e2));
                    EpiN1->Fill(cnum,(e1+e2)*EdepFactor(e1+e2)*SumRat);
                  }
                  if(k==2)
                  {
                    EpiS->Fill(cnum,(e1+e2)*EdepFactor(e1+e2));
                    EpiS1->Fill(cnum,(e1+e2)*EdepFactor(e1+e2)*SumRat);
                  };
                };

              };
            };

            if(iax>=0 && iax<7 && iay>=0 && iay<7&& (e1+e2)>40)
            {
              if(k==1 || k==2)
              {
                Float_t mm=vpi.Mag();
                Mcell2[iax][iay][k-1]->Fill(mm*EdepFactor(e1+e2)*SumRat);
                Mcell[iax][iay][k-1]->Fill(mm*EdepFactor(e1+e2));
                McellZ[iax][iay][k-1]->Fill(mm*EdepFactor(e1+e2),
                    (e1+e2)*EdepFactor(e1+e2));
              };
            };

            if(ibx>=0 && ibx<7 && iby>=0 && iby<7 && (e1+e2)>40)
            {
              if(k==1 || k==2)
              {
                Float_t mm=vpi.Mag();
                Mcell2[ibx][iby][k-1]->Fill(mm*EdepFactor(e1+e2)*SumRat);
                Mcell[ibx][iby][k-1]->Fill(mm*EdepFactor(e1+e2));
                McellZ[ibx][iby][k-1]->Fill(mm*EdepFactor(e1+e2),
                    (e1+e2)*EdepFactor(e1+e2));

              }
            };
          };
          if(YellowUp && vpi.E()>50.)PTYup.Fill(vpi.X());
          if(YellowDn && vpi.E()>50.)PTYdn.Fill(vpi.X()); 
          PtEnSp.Fill(vpi.Px(),vpi.E(),spin+.5);
          Ex.Fill(vpi.E(),spinindex);
          if(fabs(vpi.Mag()-.135)<.08){
            ExPi.Fill(vpi.E(),spinindex);
            if(vpi.E()>35)absPtPi.Fill(fabs(vpi.Px()),spinindex);
            if(fabs(vpi.Px())>2.)ExPiPt2.Fill(vpi.E(),spinindex);
            if(fabs(vpi.Px())>1.5)ExPiPt15.Fill(vpi.E(),spinindex);
          }
          if(fabs(vpi.Mag()-.55)<.15)ExEta.Fill(vpi.E(),spinindex);    
          if(fabs(vpi.Mag()-.30)<.08)ExBetw.Fill(vpi.E(),spinindex);

        };
      };

    };
  };

  //Now Save Histograms to Output.root

  TFile* out=new TFile(OutFileName,"recreate");
  EpiN->Write();
  EpiS->Write();
  EpiN1->Write();
  EpiS1->Write();
  EphN->Write();
  EphS->Write();
  ptvsx[0]->Write();
  ptvsx[1]->Write();
  sphxy[0]->Write();
  sphxy[1]->Write();
  PUENCells.Write();
  PUESCells.Write();
  PCENCells.Write();
  PCESCells.Write();

  PUENcCells.Write();
  PUEScCells.Write();
  PCENcCells.Write();
  PCEScCells.Write();
  MassENvsRun.Write();
  MassESvsRun.Write();




  P2ENCells.Write();
  P2ESCells.Write();
  ADCsumNS.Write();
  ADCsumNS2.Write();
  ENEs.Write();
  ADCNS.Write();
  ADCnor.Write();
  ADCsou.Write();
  MassEN.Write(); 
  MassES.Write(); 
  PhotonENCells.Write();
  PhotonESCells.Write();
  PTYup.Write();
  PTYdn.Write();
  EdistN.Write();
  EdistS.Write();
  PtEnSp.Write();
  Ex.Write();
  ExPi.Write();
  ExPiPt2.Write();
  ExPiPt15.Write();
  ExEta.Write();
  ExBetw.Write();
  absPtPi.Write();
  for(Int_t ix=0;ix<7;ix++)
  {
    for(Int_t iy=0;iy<7;iy++)
    {
      char nam[20];
      EdistXY[ix][iy][0]->Write();
      EdistXY[ix][iy][1]->Write();
      RawAdc[ix][iy][0]->Write();
      RawAdc[ix][iy][1]->Write();
      RawAdc[ix][iy][2]->Write();
      RawAdc[ix][iy][3]->Write();
      Mcell[ix][iy][0]->Write();
      Mcell[ix][iy][1]->Write();
      Mcell2[ix][iy][0]->Write();
      Mcell2[ix][iy][1]->Write();
      McellZ[ix][iy][0]->Write();
      McellZ[ix][iy][1]->Write();
      EvsZ[ix][iy][0]->Write();
      EvsZ[ix][iy][1]->Write();
    };
  };

  return 0;

};

TH2D* AnalTools::CutOnZ(TH3F* ph3,Int_t lowz1,Int_t hiz1,Int_t lowz2,Int_t hiz2)
{

  Int_t nbinx=ph3->GetNbinsX();
  Int_t nbiny=ph3->GetNbinsY();
  Int_t nbinz=ph3->GetNbinsZ();

  Float_t minx,miny,minz,maxx,maxy,maxz;
  minx=ph3->GetXaxis()->GetXmin();
  miny=ph3->GetYaxis()->GetXmin();
  minz=ph3->GetZaxis()->GetXmin();
  maxx=ph3->GetXaxis()->GetXmax();
  maxy=ph3->GetYaxis()->GetXmax();
  maxz=ph3->GetZaxis()->GetXmax();
  TString nnam=ph3->GetName();
  nnam=nnam+"_";nnam+=lowz1;nnam+="_";nnam+=hiz1;


  TH2D* nhist=new TH2D((const char*) nnam,(const char*) nnam,nbinx,minx,maxx,nbiny,miny,maxy);
  for(Int_t ix=0;ix<nbinx;ix++)
  {
    for(Int_t iy=0;iy<nbiny;iy++)
    {
      for(Int_t iz=0;iz<nbinz;iz++)
      {
        if((iz>=lowz1 && iz<=hiz1) || (iz>=lowz2 && iz<=hiz2) )
        {

          Float_t val=ph3->GetBinContent(ix+1,iy+1,iz+1);
          Float_t nval=nhist->GetBinContent(ix+1,iy+1)+val;
          nhist->SetBinContent(ix+1,iy+1,nval);
        };
      };
    };
  };

  return nhist;
};

Bool_t  AnalTools::SetHistFitEpt(TH2F* hist,Int_t rLim1,Int_t rLim2,Int_t cLim1,Int_t cLim2,TH2F* fhi, Geom* pGlobalGeom )
{
  PGEOM=pGlobalGeom;
  histFit_EPt=hist;
  ps_Hist=&histFit_EPt;
  histFit_EPtF=fhi;
  rowLim1=rLim1;
  rowLim2=rLim2;
  colLim1=cLim1;
  colLim2=cLim2;
  return true;
};

TH1F*  AnalTools::ReScaleHist(TH1F* p_hist,Float_t xfactor,Int_t minbin,Bool_t Normalize)
{
  /*
     based on the bins of index>=minbin rescale so bin width is multiplied by xfactor
     */
  printf("ReScaleHist called with xfactor=%f minbin=%d \n",xfactor,minbin);
  p_hist->Print();
  if(ReScaledHist)delete ReScaledHist;

  Int_t nb1x=p_hist->GetNbinsX();
  Float_t lowx=p_hist->GetBinLowEdge(1);
  Float_t highx=p_hist->GetBinLowEdge(nb1x+1);
  ReScaledHist=new TH1F("ReScaledHist","ReScaledHist",nb1x,lowx,lowx+(highx-lowx)*xfactor);
  Float_t binw0=p_hist->GetBinWidth(1);
  Float_t binw=ReScaledHist->GetBinWidth(1);
  Float_t ReScaledBins[10000];
  Float_t Sumsq=0;
  for(Int_t j=1;j<nb1x+1;j++)
  {
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
  ReScaledHist=new TH1F("ReScaledHist","ReScaledHist",nb1x,lowx,highx);

  Float_t sqrsumsq=1;
  if(Normalize)sqrsumsq=sqrt(Sumsq);

  for(Int_t j=1;j<nb1x+1;j++)ReScaledHist->SetBinContent(j,ReScaledBins[j]/sqrsumsq);
  return ReScaledHist;  
};

Float_t AnalTools::GetHistPairRescale(TH1F* p0_hist1,TH1F* p0_hist2,Float_t lowfac,Float_t highfac,Int_t firstbin)
{
  /*
     Determine the scale factor so that rescaling the p0_hist2 by the factor will match p0_hist1
     */

  Int_t nb1x=p0_hist1->GetNbinsX();
  Int_t nb2x=p0_hist2->GetNbinsX();
  printf("nb1x=%d nb2x=%d \n",nb1x,nb2x);
  Float_t lbin=p0_hist1->GetBinLowEdge(1);
  Float_t hbin=p0_hist1->GetBinLowEdge(nb2x+1);
  if(p_hist1)delete p_hist1;
  if(p_hist2)delete p_hist2;
  p_hist1=new TH1F("tmp_hist1","tmp_hist1",nb1x,lbin,hbin);
  p_hist2=new TH1F("tmp_hist2","tmp_hist2",nb1x,lbin,hbin);
  if(nb1x!=nb2x)return 0;
  Float_t sumsq1=0;
  Float_t sumsq2=0;
  for(Int_t j=0;j<nb1x+1;j++)
  {
    Float_t val=p0_hist1->GetBinContent(j);
    if(val<1)val=1;
    val=log(val);
    sumsq1+=val*val;
    p_hist1->SetBinContent(j,val);
    val=p0_hist2->GetBinContent(j);
    if(val<1)val=1;
    val=log(val);
    sumsq2+=val*val;
    p_hist2->SetBinContent(j,val);
  };
  if(sumsq1>0)p_hist1->Scale(1/sqrt(sumsq1));
  if(sumsq2>0)p_hist2->Scale(1/sqrt(sumsq2));
  printf("logs done\n");

  Double_t scalf=lowfac;
  Double_t dot=0;
  Double_t oldscalf=0.;
  Double_t olddot=0;
  Double_t del=(highfac-lowfac)/10.;
  for(Int_t j=0;j<10;j++)
  {
    scalf=scalf+del;
    dot=((*p_hist1)*(*(ReScaleHist(p_hist2,scalf,firstbin)))).Integral();
    printf( "A(%f,%f) \n",scalf,dot);
    if(olddot>dot)break;
    oldscalf=scalf;
    olddot=dot;
  };
  scalf=oldscalf-del;
  olddot=0;
  oldscalf=0;
  del=del/20.;
  for(Int_t j=0;j<40;j++)
  {
    scalf=scalf+del;
    dot=((*p_hist1)*(*(ReScaleHist(p_hist2,scalf,firstbin)))).Integral();
    printf( "B(%f,%f)\n",scalf,dot);
    if(olddot>dot)break;
    oldscalf=scalf;
    olddot=dot;
  };
  scalf=oldscalf;

  printf("scalf=%f\n",scalf);
  return scalf;
};

Int_t AnalTools::GetRunIndex(Long_t runnum,Bool_t init)
{
  if(init)
  {
    RunListLength=0;
    return 0;
  };
  for(Int_t j=0;j<RunListLength;j++)
  {
    if(runnum==RunList[RunListLength-j-1])return RunListLength-j;
  };
  RunListLength++;
  return RunListLength;
};
Float_t AnalTools::SingleEdepFactor(Int_t vtyp,Float_t energy)
{
  //  return 1.;
  if(vtyp==1)
  {
    return 1.004-.0045*(energy-60.);
  };
  if(vtyp==2)
  {
    return .984-.0045*(energy-60.);
  };
  return 1.;
};
void AnalTools::chkEta(char* filelist, float eta0,float dcenter,
    float e0, float de,float m0, float dm)
{
  printf("eta0=%f dcenter=%f e0=%f de=%f m0=%f dm=%f\n",eta0,dcenter,
      e0,de,m0,dm);
  TString home= "/star/u/heppel/ABatch";

  FilesSet* p_files61=new FilesSet(home,
      "run6_export/fpdped.txt",
      "run6_export/fpdgain.txt", 
      "run6_export/gain_set1/fpdcorr_itr11.txt",
      "run6_export/fill.txt",
      "Fake",
      "run6_export/spinpat",
      "run6_export/geom_fy06_trans_survey.txt");
  TString RUNPED=home+"/run6_export/runped/";
  p_files61->p_fpdrunped()->Directory=RUNPED;

  FilesSet* p_files62=new FilesSet("../../",
      "run6_export/fpdped.txt",
      "run6_export/fpdgain.txt", 
      "run6_export/fpdcorr_itr10.2day.txt",
      "run6_export/fill.txt",
      "Fake",
      "run6_export/spinpat",
      "run6_export/geom_fy06_trans_survey.txt");
  p_files62->p_fpdrunped()->Directory=RUNPED;
  FilesSet* p_files=p_files61;
  Fill* RFill=new Fill(p_files,6900,11000,7038002,7171011);
  poutTree* pout=new poutTree(filelist);
  TString line[10];
  TCanvas* c1=new TCanvas("c1","c1",600,600);
  c1->cd();
  c1->Divide(2,2);
  Geom* p_geom=0;
  CalibStr* p_gain=0;
  CalibStr* p_gaincorr=0;
  CalibStr* p_ped=0;
  CalibStr* p_runped=0;
  TObjArray RunsArray(100,0);
  RunData* thisRunData=0;
  Int_t runcnt=0;
  TH1F spcntH("spcntH","spcntH",5,-2.5,2.5);
  spcntH.SetTitle("Running Yellow Spin Accumulation");
  spcntH.GetXaxis()->SetTitle("[Yellow S] dot [L direction] (+1/-1)");
  spcntH.GetYaxis()->SetTitle("Number of Events");

  Long64_t nentries =pout->nentries;
  Float_t fnentries=nentries;
  printf("nentries=%f\n",fnentries);
  Int_t CurrentRun=0;
  Float_t* gn;
  Float_t* gncorr;
  Float_t* ped;
  Float_t* rped=0;
  Float_t En[98];
  TH2F* HN=0;
  TH2F* HS=0;
  Int_t buI=0;
  Int_t nI=0;
  Int_t NPass=0;
  Long64_t nbytes = 0;
  TString PSFileName="";  
  TPaveText* txt=0;
  TMatrixT<float>* tmN=new TMatrixT<float>(7,7);
  TMatrixT<float>* tmS=new TMatrixT<float>(7,7);
  Yiqun* recN=0;
  Yiqun* recS=0;
  for (Long64_t i=0; i<nentries;i++) 
  {
    pout->GetEntry(i);
    TString sptxt="";
    if(pout->YellowUp){sptxt="Yu";buI=1;};
    if(pout->YellowDn){sptxt="Yd";buI=-1;};
    pout->AllToScratch(false); //make a scratch list of all hard 
    if(i%1000==0)printf("Event %f \n",(float) i);
    Int_t nfound=0;
    Bool_t centercut=false;
    TIter next(pout->scratchlist);
    float M12,E12,Y,phi,Z12;
    M12=E12=Y=phi=Z12=0;
    if(pout->scratchlist->GetEntries()==2)
    {
      TLorentzVector* v1=(LVec*) pout->scratchlist->First();
      TLorentzVector* v2=(LVec*) pout->scratchlist->After(pout->scratchlist->First());
      Z12=1;
      if((v1->E()+v2->E())>0)Z12=fabs((v1->E()-v2->E())/(v1->E()+v2->E()));
      TLorentzVector vs=pout->SumScratch();
      phi=vs.Phi();
      if(cos(phi)<0)sptxt=sptxt+"-S";
      nI=-1;
      if(cos(phi)>0)nI=1;
      if(cos(phi)>0)sptxt=sptxt+"-N";
      Y=vs.PseudoRapidity();
      M12=vs.Mag();
      E12=vs.E();
      if(sqrt(pow(Y-eta0,2)+pow(tan(phi),2))<dcenter)centercut=true;
    };

    if(fabs(M12-m0)>dm)continue;
    if(fabs(E12-e0)>de)continue;
    if(fabs(Z12)>.8)continue;
    if(pout->spin<0 || pout->spin>4)continue;
    if(!centercut)continue;
    if(nI*buI>0)sptxt=sptxt+"[+]";
    if(nI*buI<0)sptxt=sptxt+"[-]";
    spcntH.Fill(nI*buI);
    NPass++;
    Float_t Nminus=spcntH.GetBinContent(2);
    char str[100];
    LVec* v;
    TLorentzVector vsum(0,0,0,0);
    const char* spt=(const char*) sptxt;
    sprintf(str,"%d) Rnum=%d EventN=%d (%d %s)\n",NPass,
        pout->Rnum,pout->EventN,pout->spin,spt);
    Int_t nlines=0;
    line[nlines]=str;
    nlines++;
    while(v=(LVec*) next())
    {
      vsum=vsum+*v;
      nfound++;
      sprintf(str,"(%d) Energy=%f  Y=%f Phi=%f \n",
          nfound,v->E(),v->PseudoRapidity(),v->Phi());
      line[nlines]=str;
      nlines++;
      if(nfound==2)
      {
        sprintf(str,"M12=%5.3f  Y12=%5.3f Phi12=%5.3f E12=%5.3f\n",
            vsum.Mag(),vsum.PseudoRapidity(),vsum.Phi(),E12);
        line[nlines]=str;
        nlines++;
      };

    }; 
    sprintf(str,"SumN(-) = %5.0f SumN(+)=%5.0f  AN/.6=%2.3f \n",
        1.* Nminus,NPass-Nminus,(NPass-2.*Nminus)/NPass/0.6);
    line[nlines]=str;
    nlines++;

    if(pout->Rnum!=CurrentRun)
    {

      CurrentRun=pout->Rnum;
      p_files=p_files61;
      if(CurrentRun>=7128001)p_files=p_files62;
      thisRunData=0;
      TIter next(&RunsArray);
      RunData* tmpRD;
      while(tmpRD=(RunData*) next())
      {
        if(tmpRD->RunNumber==CurrentRun)
        {
          thisRunData=tmpRD;
        };
      };
      if(thisRunData==0)
      {
        if( (CurrentRun>=RFill->FirstRun) && (CurrentRun<=RFill->LastRun) ) 
        {
          if(RFill->SetFillNumberforRun(CurrentRun)>0)
          {
            RunsArray.AddAtAndExpand(thisRunData=
                new RunData(CurrentRun,RFill,p_files)
                ,runcnt);
            runcnt++;
          };
        };

      };

      if(p_geom)delete p_geom; 
      if(p_gain)delete p_gain;
      if(p_gaincorr)delete p_gaincorr;
      if(p_ped)delete p_ped;
      if(p_runped)delete p_runped;
      p_geom=new Geom(p_files);
      p_gain=new CalibStr(CurrentRun,p_files->p_fpdgain()->Path());
      p_gaincorr=new CalibStr(CurrentRun,p_files->p_fpdgaincorr()->Path());
      p_ped=new CalibStr(CurrentRun,p_files->p_fpdped()->Path());
      if(thisRunData)p_runped=thisRunData->RunPed;
      printf("Reset Run Number to:%d \n",CurrentRun);
    };

    for(int k=0;k<2;k++)
    {

      gn=        p_gain->tm(1,k+1)->GetMatrixArray();
      gncorr=p_gaincorr->tm(1,k+1)->GetMatrixArray();
      ped=        p_ped->tm(1,k+1)->GetMatrixArray();
      rped=0;
      if(thisRunData)
      {
        rped=thisRunData->RunPed->tm(1,k+1)->GetMatrixArray();
      }
      for(int j=0;j<49;j++)
      {
        if(rped)ped[j]=ped[j]+rped[j];
        En[j+k*49]=(pout->adc[j+k*49]-ped[j])*gn[j]*gncorr[j];
      };
    };

    if(HN)delete HN;
    if(HS)delete HS;
    if(tmN)delete tmN;
    if(tmS)delete tmS;
    tmN=new TMatrixT<float>(7,7);
    tmS=new TMatrixT<float>(7,7);
    tmN->SetMatrixArray(En);
    tmS->SetMatrixArray(&(En[49]));
    if(recN)delete recN;
    if(recS)delete recS;
    printf("calling recN\n");
    recN=new Yiqun(tmN,p_geom,p_gain,p_gaincorr,1,1);
    printf("calling recS\n");
    recS=new Yiqun(tmS,p_geom,p_gain,p_gaincorr,1,2);
    printf("nphrecN=%d nphrecS=%d \n",recN->NPh,recS->NPh);
    float phE[4]={0,0,0,0};
    for(int iz=0;iz<recN->NPh;iz++)
    {
      phE[iz]=recN->mom(iz).E();
    };
    sprintf(str,"NEn=%8.3f NEne=%8.3f %8.3f %8.3f %8.3f \n",tmN->Sum(),phE[0],phE[1],phE[2],phE[3]);
    line[nlines]=str;
    nlines++;
    for(int iz=0;iz<4;iz++)
    {
      phE[iz]=0;
      if(iz<recS->NPh){
        phE[iz]=recS->mom(iz).E();
      };
    };
    sprintf(str,"SEne=%8.3f NSne=%8.3f %8.3f %8.3f %8.3f \n",tmS->Sum(),phE[0],phE[1],phE[2],phE[3]);
    line[nlines]=str;
    nlines++;
    HN=new TH2F(*tmN);
    HS=new TH2F(*tmS);
    HN->SetMaximum(60);
    HS->SetMaximum(60);
    HS->SetTitle("South FPD");
    HS->GetXaxis()->SetTitle("col");
    HS->GetYaxis()->SetTitle("row");
    HN->SetTitle("North FPD");
    HN->GetXaxis()->SetTitle("col");
    HN->GetYaxis()->SetTitle("row");

    HN->SetStats(0);
    HS->SetStats(0);
    c1->cd(3);
    HN->Draw("zcol");
    HN->Draw("sametext");
    c1->cd(4);
    c1->GetPad(3)->SetLogz();
    c1->GetPad(4)->SetLogz();
    HS->Draw("zcol");
    HS->Draw("sametext");
    c1->cd(1);
    if(txt)delete txt;
    txt=new TPaveText(.0,.3,.99,1.);  
    txt->SetTextSize(.04);
    if(nlines>8)nlines=8;

    for(Int_t k=0;k<nlines;k++)
    {
      txt->AddText(.5,.9-k*.12,line[k]);
    };
    txt->Draw();
    c1->cd(2);
    spcntH.SetFillColor(3);
    spcntH.Draw("e1");
    spcntH.Draw("same");
    c1->Update();
    if(PSFileName=="")
    {
      PSFileName="EtaEvents.ps(";
    }
    else
    {
      PSFileName="EtaEvents.ps";
    }
    if(NPass<600){ c1->Print(PSFileName);}
    else if(NPass%50==0){ c1->Print(PSFileName);};
    printf("event=%d\n Nsum=%f Ssum=%f \n",pout->EventN,tmN->Sum(),tmS->Sum());
  };
  c1->Clear();
  c1->Divide(2,1);
  c1->cd(1);
  txt->Draw();
  c1->cd(2);
  spcntH.SetFillColor(3);
  spcntH.Draw("e1");
  spcntH.Draw("same");

  c1->Print((PSFileName+")"));
};

Bool_t AnalTools::IsPhysicalFMS(Cell* cl)
{
  if(!cl)return false;
  Bool_t ret=true;

  if((cl->Instb<3) && (cl->Col1<9) && (cl->Row1>9) && (cl->Row1<26))
  {
    ret=false;
  };
  if((cl->Instb>2) && (cl->Col1<6) && (cl->Row1>7) && (cl->Row1<18))
  {
    ret=false;
  };	      
  if((cl->Instb<3) && (cl->Col1 +abs((int)(cl->Row1-17.5))>27))
  {
    ret=false;
  };
  return ret;
};

Bool_t AnalTools::IsPopulatedFMS(Cell* cl)
{
  Bool_t ret=false;

  if(AnalTools::IsPhysicalFMS(cl))
  {
    ret=true;
    Int_t adccnt=cl->p_adc->Integral( 4,1024);//Number of counts with adc > 3
    if(adccnt<5)ret=false;
  };
  return ret;
};

void AnalTools::TrimFMS(CalibStr* FF)
{
  CalibStr* F=(CalibStr*) FF;
  printf("hi\n");
  F->tm(2,1)->Print();
  // zero gcorr for unphysical cells
  for(int id=0;id<4;id++)
  {
    int nrow=34;
    if(id>1)nrow=24;

    for(int ir=0;ir<nrow;ir++)
    {
      for(int ic=0;ic<nrow/2;ic++)
      {

        printf("id=%d ic=%d ir=%d \n",id,ic,ir);
        if(id<2 && ic<8&& ir>8 && ir<25)
        {
          F->SetValue(2,id+1,ir,ic,0.);
        };
        if(id>1 && ic<5&& ir>6 && ir<17)
        {
          F->SetValue(2,id+1,ir,ic,0.);
        };
        if(ic+abs((int)(ir-16.5))>26)
          F->SetValue(2,id+1,ir,ic,0.);
      };
    };
  };

};

void AnalTools::SetCorrLimitFMS(CalibStr* FpdCorr,Float_t max,Float_t min)
{
  Int_t nr[4]={34,34,24,24};
  for(int k=0;k<4;k++)
  {
    for(int r=0;r<nr[k];r++)
    {
      for(int c=0;c<nr[k]/2;c++)
      {
        Float_t val= FpdCorr->GetValue(2,k+1,r,c);
        if(val>max)val=max;
        FpdCorr->SetValue(2,k+1,r,c,val); 
      };
    };

  };


};
void AnalTools::CheckCluster(Yiqun* rec)
{
  Int_t NSTB=rec->NSTB;
  Int_t NPh=rec->NPh;
  for(int j=0;j<NPh;j++)
  {
    printf("E%d = %f x,y=(%f,%f) ",j,rec->mom(j).E(),rec->photons[j].xPos,
        rec->photons[j].yPos);
  };
  printf("\n");
  if(rec->Clu2)
  {
    printf("NSTB=%d\n",NSTB);
    TIter next(rec->Clu2);
    HitCluster* hcl=0;
    Int_t cnt=0;

    while( hcl=(HitCluster*) next())
    {
      rec->fitter->pTowerUtil->CalClusterMoment(hcl);
      hcl->FindClusterAxis();
      printf("Cluster cnt=%d \n",cnt++);
      TowerFPD* last=(TowerFPD*)  hcl->tow->Last();
      Int_t scol=last->col;
      Int_t srow=last->row;
      Int_t sew=rec->EW;
      Int_t snstb=rec->NSTB;
      Float_t xtow=(scol-.5)*rec->widLG[0];
      Float_t ytow=(srow-.5)*rec->widLG[1];
      printf("x_clu=%f y_clu=%f \n",xtow,ytow);

      hcl->Print();
      PrClHit(SaveClHit(last,sew,snstb,srow,scol));
    };
  };
};


unsigned int AnalTools::SaveClHit(TowerFPD* h,Int_t IEW,Int_t INSTB,Int_t row,Int_t col)
{
  unsigned int p1,p2,p3,p4,adc;
  p1=((IEW-1)&1)*     0x80000000;
  p2=((INSTB-1)&7)*   0x10000000;
  p3=((row-1)&0x3F)*  0x00400000;
  p4=((col-1)&0x3F)*  0x00010000;
  adc=(h->adc_over_ped&0xFFFF)+p1+p2+p3+p4;
  if(nSavedHits<MaxSavedHits)
  {
    SavedHits[nSavedHits]=adc;
    nSavedHits++;
  };
  return adc;
};
void AnalTools::PrClHit(unsigned int s0)
{
  unsigned int s=s0;
  Int_t sew,snstb,srow,scol,sadc;
  sew=1;
  if(s&&0x80000000)sew=2;
  s=s&0x7FFFFFFF;
  snstb=((s/0x10000000)&7)+1;
  srow= ((s/0x00400000)&0x3F)+1;
  scol= ((s/0x00010000)&0x3F)+1;
  sadc= s&0xFFF;
  printf("SavedHit= %x: ew=%d nstb=%d row=%d col=%d adc=%d \n",s0,sew,snstb,srow,scol,sadc);
}
void AnalTools::storerec(Yiqun* rec,unsigned int minADC)
{
  if(rec)
  {
    if(rec->tow_Arr)
    {
      TowerFPD* tow=0;
      TIter nexthit(rec->tow_Arr);
      while(tow=(TowerFPD*) nexthit())
      {
        if(tow->adc_over_ped>=minADC)
        {
          Int_t scol=tow->col;
          Int_t srow=tow->row;
          Int_t sew=rec->EW;
          Int_t snstb=rec->NSTB;
          SaveClHit(tow,sew,snstb,srow,scol);
        }; 
      };
    };
  }
};
void AnalTools::storeCluster(Yiqun* rec,unsigned int minADC,Float_t minclusterE)
{
  // loop through clusters in NSTB of event (since Yiqun reconstruction called for each NSTB of each event separately)
  for(int ncl=0;ncl<rec->NRealClusts;ncl++)
  {

    // loop through photons within the cluster
    for(int np=0;np<rec->clust[ncl].nPhoton;np++)
    {
      SavedCluHitIndex[nCluster]=nSavedHits;
      SavedPhotonCluIndex[nPhotonClu]=nCluster;
      Float_t ene=rec->clust[ncl].photon[np].energy;
      Float_t xp=rec->clust[ncl].photon[np].xPos;
      Float_t yp=rec->clust[ncl].photon[np].yPos;
      SavedPhotonCluEnergy[nPhotonClu]=ene;
      nPhotonClu++;
      //	  printf("nstb=%d cluster #=%d ph#=%d E=%f x=%f y=%f \n",rec->NSTB,ncl,np,ene,xp,yp);

    };
    if(rec->clust[ncl].energy>minclusterE)
    {
      // loop through hits which comprise the cluster
      TIter next(rec->clust[ncl].tow);
      while(TowerFPD* tow=(TowerFPD*) next())
      {
        if(tow->adc_over_ped>=minADC)
        {
          Int_t scol=tow->col;
          Int_t srow=tow->row;
          Int_t sew=rec->EW;
          Int_t snstb=rec->NSTB;
          SaveClHit(tow,sew,snstb,srow,scol);
        };
      };
    };
    nCluster++;
    SavedCluHitIndex[nCluster]=nSavedHits;
  };
}

Int_t AnalTools::ClNearX(Yiqun* prec,Float_t x1, Float_t y1)
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

Int_t AnalTools::PhNearX(Yiqun* prec,Float_t x1,Float_t y1)
{
  Int_t n=-1;
  Float_t dmin=100;
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


Int_t AnalTools::Select_Trig(char* infile,Int_t set,FilesSet* P_Files,char* hntName)
{
  FilesSet* p_files=0;
  TString hnt_name="";
  Bool_t clip=false;
  Bool_t hardtrig=false;
  Bool_t AnalDet[2][6];
  for(Int_t i=0;i<2;i++){for(Int_t j=0;j<6;j++)AnalDet[i][j]=false;};
  AnalDet[0][0]=false;
  AnalDet[0][1]=false;
  AnalDet[1][0]=true;
  AnalDet[1][1]=true;
  AnalDet[1][2]=true;
  AnalDet[1][3]=true;
  if(P_Files==0)
  {
    if(set==80)
    {
      //transverse run 8
      p_files=new FilesSet("/star/u/heppel/BatchDir/rootfms",
          "fpdped.txt",
          "fpdgain.txt", 
          "fpdcorr.txt",
          "fill.txt",
          "Fake",
          "fpd08/spinpat",
          "geom.txt",
          "qtmap.txt",
          "qtmap2pp.txt");


      p_files->Print();
      hnt_name="h111_08";

    }
    else
    {
      printf("no valid set selected for Select_scriptQT1\n");
      return 0;
    };
  }
  else
  {
    p_files=P_Files;
    hnt_name=hntName;
  };

  char hntp[256];
  strcpy(hntp,hnt_name);
  printf("dataSet d(%s,p_files,%s)\n",infile,hntp);
  dataSet d(infile,p_files,hntp);
  std::cout<<"ok\n";
  if(p_files==0){printf( " Exit unknown data set\n");return -2;};
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -1;
  };

  d.GetEntry(1);
  if(d.error!=0)
  {
    printf("error %d making data set\n",d.error);
    return -2;
  };

  Int_t Rnum=d.CurrentRunNumber;
  Int_t Bunchid7bit=d.bunchid7bit;

  // FpdMap is a class that figures out maps into data blocks
  //  FpdMap* pmap=new FpdMap(Rnum);

  string s(infile);
  string s1 = s.substr(s.find("run10"));
  sprintf(hntp, "output/OFile.%s", &s1[0]);
  TFile* p_OFile=new TFile(hntp,"recreate");

  //
  Int_t tpes[320];
  Int_t spin;
  Int_t nphotons;
  Float_t px[10];
  Float_t py[10];
  Float_t pz[10];
  Float_t pE[10];
  Float_t pxyzt[320];
  Int_t nwrds;
  Int_t strig;
  Char_t BBcSums[5];
  Float_t BBcVertex[7];
  Int_t EventN;
  Int_t ievt;
  //Define a Tree
  TTree* p_out = new TTree("p_out","Output Tree");
  TH1F Emat[5];
  Emat[0]=TH1F("Emat0","Emat0",2000,0,2000);
  Emat[1]=TH1F("Emat1","Emat1",2000,0,2000);
  Emat[2]=TH1F("Emat2","Emat2",2000,0,2000);
  Emat[3]=TH1F("Emat3","Emat3",2000,0,2000);
  Emat[4]=TH1F("EmatAll","EmatAll",2000,0,2000);
  TH1F NHi[5];
  NHi[0]=TH1F("NHi0","NHi0",300,0,300);
  NHi[1]=TH1F("NHi1","NHi1",300,0,300);
  NHi[2]=TH1F("NHi2","NHi2",300,0,300);
  NHi[3]=TH1F("NHi3","NHi3",300,0,300);
  NHi[4]=TH1F("NHiAll","NHiAll",300,0,300);


  TH2F hSml0("hSml0","hSml0",24,-12.,12.,24,0.,24.);
  TH2F hSml1("hSml1","hSml1",24,-12.,12.,24,0.,24.);
  TH2F hSml2("hSml2","hSml2",24,-12.,12.,24,0.,24.);
  TH2F* hSml[3];
  hSml[0]=&hSml0;hSml[1]=&hSml1;hSml[2]=&hSml2;

  TH2F hLrg0("hLrg0","hLrg0",34,-17.,17.,34,0.,34.);
  TH2F hLrg1("hLrg1","hLrg1",34,-17.,17.,34,0.,34.);
  TH2F hLrg2("hLrg2","hLrg2",34,-17.,17.,34,0.,34.);
  TH2F* hLrg[3];
  hLrg[0]=&hLrg0;hLrg[1]=&hLrg1;hLrg[2]=&hLrg2;

  TH1F hBits("hBits","hBits",16,0,16);

  p_out->Branch("spin",&(spin),"spin/I");
  p_out->Branch("nphotons",&(nphotons),"nphotons/I");
  p_out->Branch("br_nwrds",&(nwrds),"nwrds/I");
  p_out->Branch("br_types",(tpes),"tpes[nwrds]/I");
  p_out->Branch("br_pxyzt",(pxyzt),"pxyzt[nwrds]/F");
  p_out->Branch("br_Rnum",&(Rnum),"Rnum/I");
  p_out->Branch("br_Bunchid7bit",&(Bunchid7bit),"Bunchid7bit/I");
  p_out->Branch("br_BBcSums",BBcSums,"BBcSums[5]/b");
  p_out->Branch("br_BBcVertex",BBcVertex,"BBcVertex[7]/F");


  p_out->Branch("br_EventN",&(EventN),"EventN/I");
  p_out->Branch("br_ievt",&(ievt),"ievt/I");
  p_out->Branch("br_nSavedHits",&(nSavedHits),"nSavedHits/I");
  p_out->Branch("br_SavedHits",(SavedHits),"SavedHits[nSavedHits]/I");


  //
  Int_t pcnt=1;
  Int_t nentries= d.Input->GetEntries();
  printf("NumberEventsToProcess found = %d ",NumberEventsToProcess);
  if((NumberEventsToProcess>0) && (NumberEventsToProcess<nentries))nentries=NumberEventsToProcess;
  printf(" nentries will be = %d \n",nentries);
  // event loop
  Int_t skpcnt=0;
  UInt_t lastLED=0;

  TMatrix* p_adc[4]={0,0,0,0};
  TMatrix* p_Emat[4]={0,0,0,0};
  Vertex* vtx=new Vertex();
  //vtx->Bbc_threshold=15;

  for(Int_t i=0;i<nentries;i++)
  {
    //if(i%1000==0)
    printf("Processing event %d **********\n",i);
    nSavedHits=0;
    Int_t nbytes=d.GetEntry(i);
    Bool_t trgXing=true;
    Rnum=d.CurrentRunNumber;
    if(Rnum<10000000)
    {
      if(Rnum<10000000 &&(d.ipre!=0 || d.ipost!=0)){trgXing=false; continue;};
    };
    EventN=d.event;
    ievt=i;
    if(!d.decodeQT()){printf("qt Error \n"); continue;};
    if(!d.decodeDSM()){printf("qt Error \n"); continue;};


    Bunchid7bit=d.bunchid7bit;
    pcnt++;
    if(pcnt>1000){printf("event=%d \n",i);pcnt=1;}; //print every 5000 events
    Double_t NearUnique;

    spin=40;
    if(abs(d.BlueSpin*d.YellowSpin)==1)
    {
      spin=(d.BlueSpin+1)+(d.YellowSpin+1)/2;
    };
    if(d.kicked)spin+=29;
    if(d.BlueSpin*d.YellowSpin==0)spin=40;
    Trigger trig(d.Bbcl1);
    for(Int_t ki=0;ki<5;ki++)BBcSums[ki]=trig.BBcSums[ki];
    BBcVertex[0]=vtx->GetBbcVertex(d.Bbc);
    BBcVertex[1]=vtx->maxtace;
    BBcVertex[2]=vtx->maxtacw;
    BBcVertex[3]=vtx->iemax;
    BBcVertex[4]=vtx->iwmax;
    BBcVertex[5]=vtx->qe;
    BBcVertex[6]=vtx->qw;
    Float_t Esum[4]={0.,0.,0.,0.};
    Float_t Esum4=0.;
    Int_t Nhigh[4]={0,0,0,0};
    Int_t Nhigh4=0;
    for(Int_t k=0;k<4;k++)
    {
      if(AnalDet[1][k]){
        p_adc[k]=new TMatrix(d.dMatrixQt(k));
        p_Emat[k]=new TMatrix(d.Em(p_adc[k],2,k+1));
        Esum[k]=p_Emat[k]->Sum();
        Emat[k].Fill(Esum[k]);
        Esum4+=Esum[k];
        Nhigh[k]= NTowersAbove(p_Emat[k],.25);
        Nhigh4+=Nhigh[k];
        NHi[k].Fill(Nhigh[k]);
      };
    };
    Emat[4].Fill(Esum4);
    NHi[4].Fill(Nhigh4);
    hardtrig=1;
    nwrds=0;
    nphotons=0;

    Bool_t LEDevent=false;
    UInt_t LEDbit=d.lastdsm[7];
    Int_t bc=d.Bclo+d.Bchi*(Int_t)pow(2,32);
    //printf("%d %d %d \n",bc,(UInt_t)d.ipre,(UInt_t)d.ipost);
    if(Rnum<10000000)
    {
      if((LEDbit/2048)%2==1)
      {
        LEDevent=true; 
        lastLED=bc;
      };
      if(lastLED!=0 && ((bc-lastLED)%2000001<10) || (bc-lastLED)%2000001>1999997)
      {
        LEDevent=true;
      };
    }     
    else
    {
      if((LEDbit)%16>0)
      {
        LEDevent=true; 
      };	  
    };
    Bool_t Debug=false;
    if(Debug)
    {
      if(LEDevent)printf("lastLED=%d bc=%d mod_bc=%d Nhigh4=%d \n",lastLED,bc,(bc-lastLED)%2000001,Nhigh4);      
    };

    TLorentzVector vecs[80];
    for(Int_t kj=0;kj<320;kj++)tpes[kj]=0; 
    Bool_t pass=Esum4>40 || (Esum4>30 && (fmod((skpcnt++),10)<OutOf10-.1));
    pass=Esum[0]>8 || Esum[1]>8 || Esum[2]>15 || Esum[3]>15;

    UInt_t fastbit = (d.Dsm & 0x1000)>>12;
    UInt_t slowbit = (d.Dsm & 0x2000)>>13;
    UInt_t ledbit =  (d.Dsm & 0x4000)>>14;
    //cout<<"fastbit="<<fastbit<<" slowbit="<<slowbit<<" ledbit="<<ledbit<<endl;
    if(ledbit==0)
      //if(!LEDevent)
    {
      //cout<<hex<<*(d.Trgwd)<<" "<<*(d.Trgwd+1)<<" "<<*(d.Trgwd+2)<<dec<<endl;
      //(d.pQt)->printADC();
      //(d.pQt)->printDSM();
      //(d.pQt)->printQt();
      TrigQt* trg = new TrigQt(d.pQt,true);
      trg->L0dsm();
      //trg->printDSM();
      trg->AddRealEdgeClusters();
      trg->L12dsm();

      //cout<<"data L1 trigger bits: "<<hex<<(d.pQt)->qtdsm2_trgbits[0]<<" "<<(d.pQt)->qtdsm2_trgbits[1]
      //<<" "<<(d.pQt)->qtdsm2_trgbits[2]<<dec<<endl;
      //cout<<"data L2 trigger bits: "<<hex<<(d.pQt)->qtdsm_trgbits<<dec<<endl;
      //cout<<"simu L1 trigger bits: "<<hex<<trg->L1_TrigBits[0]<<" "<<trg->L1_TrigBits[1]<<" "<<trg->L1_TrigBits[2]<<dec<<endl;
      //cout<<"simu L2 trigger bits: "<<hex<<trg->L2_TrigBits<<dec<<endl;
      if((d.pQt)->qtdsm_trgbits != trg->L2_TrigBits) cout<<"Warning: unmatched trigger bits on event "<<i<<endl;
      UInt_t prescale = 1; //tmp
      UInt_t selectbits;
      selectbits = ((d.pQt)->qtdsm_trgbits>>1&0x01) | ((d.pQt)->qtdsm_trgbits>>9&0x01);
      selectbits += (((d.pQt)->qtdsm_trgbits>>6&0x01) | ((d.pQt)->qtdsm_trgbits>>13&0x01))<<1;
      selectbits += ((((d.pQt)->qtdsm_trgbits>>3&0x01)) | (((d.pQt)->qtdsm_trgbits>>10&0x01)) |
          (((d.pQt)->qtdsm_trgbits&0x01) & ((d.pQt)->qtdsm_trgbits>>7&0x01)))<<2;
      selectbits += ((((d.pQt)->qtdsm_trgbits&0x01) | ((d.pQt)->qtdsm_trgbits>>7&0x01)) & (prescale))<<3;
      UInt_t selectb = selectbits&0x0d;

      for(int ii=0; ii<14; ii++)
      {
        int bit = ((d.pQt)->qtdsm_trgbits>>ii) & 0x1;
        if(bit==fastbit) hBits.Fill(ii);
      }

      TIter next(trg->TrigClusters);
      for(Int_t j=0;j<40;j++)//loop over Qt32 boards
      {
        ClusterQt* Cl = (ClusterQt*)next();
        if(Cl->Defined)
        {
          //Cl->Print();
          Int_t crt=Cl->crate;
          Int_t slt=Cl->slot;
          Int_t dchan=Cl->detchan-1;
          Int_t dnstb=Cl->detnstb-1;
          Int_t clsum=Cl->sum;
          Int_t cladr=Cl->qthtadr;//0~31
          Int_t row=-1;
          Int_t col=-1;
          Int_t dsm1 = slt/4;
          if(dsm1==0)
          {
            row=dchan/12;
            if(dnstb==2)col=-1-dchan%12;
            if(dnstb==3)col=dchan%12;
            if(clsum>trg->thS[0] || clsum>trg->thS[1] || clsum>trg->thS[2])
            {
              hSml[0]->Fill(col,row);
            }
          }
          else
          {
            row=dchan/17;
            if(dnstb==0)col=-1-dchan%17;
            if(dnstb==1)col=dchan%17;
            if(clsum>trg->thL[0] || clsum>trg->thL[1] || clsum>trg->thL[2])
            {
              hLrg[0]->Fill(col,row);
            }
          }
        }
      }

      delete trg;
    }

    if(!LEDevent &&hardtrig==1 && Nhigh4<70 && Esum4<400  && pass )
    {
      for(Int_t k=0;k<4;k++)
      {
        if(Esum[k]>1. && Esum[k]<110.)
        {
          TLorentzVector lv(FourMom(p_Emat[k],0,100,0,100,d.pGeom,2,k+1));
          tpes[nwrds*4]=k+5;
          vecs[nwrds]=lv;
          nwrds++;
          Yiqun recon(p_Emat[k],d.pGeom,d.Rgain,d.Rgaincorr,2,k+1);
          Int_t nz=recon.NPh;
          if(nz>10)nz=10;
          nphotons+=nz;
          for(Int_t n=0;n<nz;n++)
          {
            vecs[nwrds]=recon.mom(n);
            tpes[nwrds*4]=305+k;
            nwrds++;
          };
          //		  if(Esum[k]>10 && k>1)CheckCluster(&recon);
          if(Esum[k]>10 && k>-1)storeCluster(&recon,1,1.5);
        };
      };
    };
    //      if(nSavedHits>0)printf("number stored =%d \n",nSavedHits);
    Int_t wcnt=0;
    for(Int_t k=0;k<nwrds;k++)
    {
      pxyzt[wcnt]=vecs[k].Px();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Py();
      wcnt++;
      pxyzt[wcnt]=vecs[k].Pz();
      wcnt++;
      pxyzt[wcnt]=vecs[k].E();
      wcnt++;
    };
    nwrds=nwrds*4;

    if(nwrds>0)p_out->Fill();
    for(Int_t k=0;k<4;k++)
    {
      if(p_adc[k]){delete p_adc[k];p_adc[k]=0;};
      if(p_Emat[k]){delete p_Emat[k];p_Emat[k]=0;};
    };
  };

  //end event loop
  p_out->Print();
  p_out->Write();
  for(int jj=0;jj<3;jj++)
  {
    hSml[jj]->Write();
    hLrg[jj]->Write();
  }
  hBits.Write();
  /*
     for(int jj=0;jj<5;jj++)
     {
     Emat[jj].Write();
     NHi[jj].Write();
     };
     */
  return 0;
};
void AnalTools::fcn_gra(Int_t &npar,Double_t* grad,Double_t& fval,Double_t* par,Int_t& iflag)
{
  //fcn for fit of p_gra graph to parameters a1,a2,a3,b1,b2,b3
  Double_t val=0;

  Double_t w=5.81;
  Double_t FCN=0;
  Double_t zval=AnalTools::ishape(par,0.,.0,w);
  Double_t gzval;
  //drive origin value to FreeGlobals[0]
  if( (gzval=FitTower::we.FreeGlobals[0])>0)FCN=pow((gzval-zval)/.001,4);
  // vnorm is the integral of the fshpe function over 7x7 array
  printf("par=%f,%f,%f,%f,%f,%f\n",par[0],par[1],1.-par[0]-par[1],par[2],par[3],par[4]);

  Double_t vnorm=AnalTools::ishape(par,0.,0.,w*3.5);
  printf("vnorm=%f zval=%f \n",vnorm,zval);
  TH2D* t2d=p_gra2;
  vnorm=1.;// not now
  int cntc=0;
  for(float r2=0;r2<3.5;r2+=.25) 
  {
    float r1=r2;
    for(float th=0.001;th<=3.14159/2+.01;th+=3.14159/12.1)
    {
      float x=r1*cos(th);
      float y=r1*sin(th);
      Int_t binx,biny,binz;
      Int_t bin=t2d->FindBin(x,y);
      t2d->GetBinXYZ(bin,binx,biny,binz);
      //	    if(bin>0 && binx>0 && binx<40 && biny>0 && biny<40)
      if(bin>0 && binx>0 && binx<13 && biny>0 && biny<13)
      {
        x=(1.*binx-.5)*.25;
        y=(1.*biny-.5)*.25;
        cntc++;
        val=0;
        if(p_gra==0) printf("p_gra zero\n");
        //	    Double_t fgr=p_gra->Interpolate(x,y);

        Double_t fgr=p_gra->Interpolate(x,y);
        if(fgr==0)fgr= t2d->GetBinContent(bin);
        //	    printf("fgr=%f \n",fgr);
        val=AnalTools::ishape(par,x,y,w);
        val=val/vnorm;
        float enom=(float) FitTower::we.Energy_study;
        float p1=FitTower::we.Power1;
        float p2=FitTower::we.Power2;
        float efact=FitTower::we.errFactor;
        float errQ=FitTower::we.errQ;
        float err2=efact*pow(val,p1)*pow(1-val,p2)*enom+errQ;
        err2=err2/enom/enom;
        if(sqrt(x*x+y*y)<1.2)err2=err2/20;
        if(fabs(x-3*w/4)<1 || fabs(y-3*w/4)<1)err2=err2/3.;
        if(p_doverr)
        {
          int nbins=p_doverr->GetNbinsX()*p_doverr->GetNbinsY();
          int bnum=p_doverr->FindBin(x,y);
          if(bnum>0 && bnum<=nbins)
          {
            p_doverr->SetBinContent(bnum,(val-fgr)/sqrt(err2));
          };
        };
        FCN+=((val-fgr)*(val-fgr)/err2);

        //		printf("x=%f y=%f val=%f err=%f fgr=%f  FCN=%f \n",x,y,val,sqrt(err2),fgr,FCN);

        //	    printf("x =%f y=%f val=%f fgr=%f FCN=%f \n",x,y,val,fgr,FCN);
      };
    };
  };
  fval= FCN/cntc;
};
void AnalTools::Set_gra2(TH2D* Gra2,TH2D* gdoe)
{p_gra2=Gra2;p_doverr=gdoe;};
void AnalTools::Set_gra(TGraph2DErrors* Gra)
{p_gra=Gra;lp_gra=Gra;};

Double_t AnalTools::ishape(Double_t* par,Double_t x, Double_t y,Double_t w0)
{
  Double_t w=w0/2.;
  Double_t val=0;
  val+=AnalTools::fshape(par,x+w,y+w);
  val+=AnalTools::fshape(par,x-w,y-w);
  val+=-AnalTools::fshape(par,x+w,y-w);
  val+=-AnalTools::fshape(par,x-w,y+w);
  return val;
};

Double_t AnalTools::fshape(Double_t* par,Double_t x, Double_t y)
{
  Double_t a[3];
  Double_t b[3];
  a[0]=par[0];
  a[1]=par[1];
  a[2]=1.-par[0]-par[1];
  b[0]=par[2];
  b[1]=par[3];
  b[2]=par[4];

  Double_t ff[3];
  ff[0]=atan2((x*y)/b[0],sqrt(x*x+y*y+b[0]*b[0]))/2./3.1415926;
  ff[1]=atan2((x*y)/b[1],sqrt(x*x+y*y+b[1]*b[1]))/2./3.1415926;
  ff[2]=atan2((x*y)/b[2],sqrt(x*x+y*y+b[2]*b[2]))/2./3.1415926;
  Double_t rval=a[0]*ff[0]+a[1]*ff[1]+a[2]*ff[2];
  return rval;

};
void AnalTools::RCorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW)
{
  if(currentdset==0)return;
  Int_t RN=currentdset->CurrentRunNumber;
  Int_t SN=currentdset->CurrentSegNumber;
  Int_t EvN= currentdset->CurrentEventNumber;
  RCorrect(padc,pEmat,det,iEW,RN,EvN,SN);
}
void AnalTools::RCorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW,Int_t runnum,Int_t EvNum,Int_t SegNM)
{
  if(runnum<12000000)return;
  Int_t NROWS[4]={34,34,24,24};
  if(Rdep)
  {
    if(runnum!=Rdep->RunNumber)
    {
      delete Rdep;
      Rdep=0;
    };
  };
  if(Rdep==0)
  {
    printf("Run changed to %d\n",runnum);
    TFile* SaveDefaultFile=gROOT->GetFile();
    TString RDPath=gSystem->Getenv("RunDepPath");
    TFile* RdepFile=new TFile(RDPath);//set in SetFMSEnv
    if(RdepFile)
    {
      printf("opened %s \n",(const char*) RDPath);
      if(Mgr==0)
      {
        if(! RdepFile->IsOpen())
        {
          delete RdepFile;
        };
        printf("Getting new Mgr from RdepFile\n");
        Mgr=(RunDepMgr*) RdepFile->Get("Mgr");
        Mgr->RunDepPath=RDPath;
        if(RunDepBaseOverride>12000000)
        {
          Mgr->RunDepBaseOverride=RunDepBaseOverride;
        };
      }
      delete RdepFile;
      char sfl[200];	  
      sprintf(sfl,"RDC%d",runnum);
      printf("%s\n",sfl);
      Rdep=Mgr->SetRdep(runnum);
      printf("Return from Mgr->SetRdep()\n");
      if(!Mgr->Rdep)return;//failed to find file
      if(!Rdep)printf("file %s not found \n",sfl);

      if(RunDepBaseOverride && Rdep)
      {
        if(RunDepBaseOverride<12000000)
        {
          BaseRdep=Rdep;
          Mgr->SetBase(runnum);
        }
        else
        {
          char sflbase[200];
          // sprintf(sflbase,"RDC%d",RunDepBaseOverride);
          // printf("opening %s\n",sflbase);
          // BaseRdep=(RunDepCor*) RdepFile->Get(sflbase);
          Mgr->SetBase(RunDepBaseOverride);
          if(!Mgr->Base)
          {
            printf("Base not found\n");
            return;
          }
          BaseRdep=Mgr->Base;
          printf("Return from SetBase \n");
          if(BaseRdep==0)
          {
            printf("Requested base run %d not available, using self\n",RunDepBaseOverride);
            BaseRdep=Rdep;
          };
        };
      };	  
    };
    if(SaveDefaultFile)SaveDefaultFile->cd();
  };
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
  CellTDep* ctm;


  for(int ro=0;ro<nr;ro++)
  {
    for(int co=0;co<nc;co++)
    {
      if(((*padc)(ro,co))>1)
      {
        if(Rdep->Legal(2,det+1,ro,co))
        {
          ct=Rdep->Cdep(det+1,ro,co);
          ctm= Mgr->Rdep->Cdep(det+1,ro,co);
          Float_t ffac=ct->FixFactor(SegNM,EvNum,2);
          Float_t mfac=ctm->FixFactor(SegNM,EvNum,2);
          Float_t ledfac=Mgr->ledFactor[det][ro][co];
          printf("RCorrect:EvNum=%d\n",EvNum);
          if(EvNum==31)
          {
            printf("RunN= %d Base->RunN= %d Mgr->Base \n",
                (int)  Rdep->RunNumber,(int) Rdep->BaseRunNumber);

            printf("BaseRdep->RunNumber=%d ct->GetnominalLedADC=%d ct->GetthisLedADC()=%d \n",
                (int) BaseRdep->RunNumber,(int) ct->GetNominalLedADC(),(int) ct->GetthisLedADC());

            (*pEmat)(ro,co)=((*pEmat)(ro,co))*ffac;
            printf("nstb=%d ro=%d co=%d ffac=%f\n",det+1,ro,co,ffac);

            printf("                    (mfac=%f) x  (ledfac=%f) =%f  \n",mfac,ledfac,mfac*ledfac);
          };
          (*pEmat)(ro,co)=((*pEmat)(ro,co))*mfac*ledfac;
        };
      }
    };
  };
  //  printf("\nCorrected Esum=%f \n",pEmat->Sum());

};
Bool_t AnalTools::AddRunToRunDepCorFile(Int_t runNumber,char* adcTrLpath,Int_t baseRunNumber,char* scratchDir)
{
  const char* RunDepPath=gSystem->Getenv("RunDepPath");

  if( RunDepPath=="" )
  {
    printf("Error making RunDependent file: environment RunDepPath variable must point to a root file");
    return false;
  };
  TString RDpath=RunDepPath;
  if(! gSystem->FindFile("./",RDpath))
  {
    printf("Error making RunDependent file: $RunDepPath file not found\n");
    return false;
  };
  TString ADCLpath=adcTrLpath;
  if(! gSystem->FindFile("./",ADCLpath))
  {
    printf("Error making RunDependent file: %s does not exist\n",adcTrLpath);
    return false;
  };

  TString FMSTXT =gSystem->ExpandPathName("${FMSTXT}");
  FMSTXT=FMSTXT+"/";
  TString FMSCORR =gSystem->ExpandPathName("${FMSCORR}");
  TString Tgc=FMSTXT+"/"+FMSCORR;
  CalibStr* gcor=new CalibStr(9010000,Tgc);
  gSystem->Setenv("RDCtmpDir",scratchDir);

  TFile* RdepFile=new TFile("$RunDepPath","update");//set in SetFMSEnv;
  RunDepCor* newRDC=0;
  RunDepCor* baseRDC=0;
  char nRDCname[200];
  char baseRDCname[200];
  sprintf(nRDCname,"RDC%d",runNumber);
  sprintf(baseRDCname,"RDC%d",baseRunNumber);
  if(baseRunNumber==runNumber)
  {
    newRDC=new RunDepCor(runNumber,adcTrLpath,gcor);
  }
  else
  {
    baseRDC=(RunDepCor*) RdepFile->Get(baseRDCname);
    newRDC=new RunDepCor(runNumber,adcTrLpath,gcor,baseRDC);      
  };
  RdepFile->cd();
  newRDC->Write();
  RdepFile->Close();
  delete RdepFile;
  delete gcor;
};



Int_t AnalTools::OutputTwoTr(char* filelist, FilesSet* p_files) {
  // creates a small Output.root file (minimum) from a OFile -- similar to readq1 but with much fewer branches in the output file - by Saroj Adhikari - for learning purposes. Only one tree TwoTr (similar to the TwoTr of readq1) is created with only some properties thought of to be useful in preliminary spin analysis//
  TFile* out = new TFile(OutFileName,"recreate");  
  TTree TwoTr("TwoTr", "consists of pair events and their various properties");

  Int_t keyIndex=kSingleKey;
  Float_t maxsep=0.07; 
  Float_t minsep=0.0;
  Float_t ECutMin=60.;
  Float_t ECutMax=90.;
  Float_t ZMax=0.75;


  typedef struct{
    Int_t Bunchid7bit;
    Int_t spin;
    Int_t N12;
    Float_t M12;
    Float_t E12;
    Float_t Pt;
    Float_t Eta;
    Float_t Ntracks;
    Float_t Phiaway;
    Float_t Ptaway;
    Float_t Maway;
    Float_t Yaway;
    Float_t Y1;
    Float_t Phi;   
    Float_t Z;
  } Properties_t;
  Properties_t p2;

  TwoTr.Branch("spin", &p2.spin, "spin/I");
  TwoTr.Branch("Bunchid7bit", &p2.Bunchid7bit, "Bunchid7bit/I");
  TwoTr.Branch("N12", &p2.N12, "N12/I");
  TwoTr.Branch("M12", &p2.M12, "M12/F");
  TwoTr.Branch("E12", &p2.E12, "E12/F");
  TwoTr.Branch("Pt", &p2.Pt, "Pt/F");
  //TwoTr.Branch("Pttest", &p2.Pttest, "Pt/F");
  TwoTr.Branch("Eta", &p2.Eta, "Eta/F");
  TwoTr.Branch("Ntracks", &p2.Ntracks, "Ntracks/F");
  TwoTr.Branch("Phiaway", &p2.Phiaway, "Phiaway/F");
  TwoTr.Branch("Ptaway", &p2.Ptaway, "Ptaway/F");
  TwoTr.Branch("Maway", &p2.Maway, "Maway/F");
  TwoTr.Branch("Yaway", &p2.Yaway, "Yaway/F");
  TwoTr.Branch("Y1", &p2.Y1, "Y1/F");
  TwoTr.Branch("Phi", &p2.Phi, "Phi/F");
  TwoTr.Branch("Z", &p2.Z, "Z/F"); //Z is a measure of the energy differences between two photons in case of N12==2; if one takes most of the energy then, Z is large and it is more likely that the event is not a real one as it is very easy to get a low energy photon in the background. a value of 1 is assigned for N12!=2.

  pout = new poutTree(filelist);
  pout->EnableEdepCorr=0; // do not enable Edep corrections

  //set energy threshold to include photons
  for(Int_t kk=1;kk<9;kk++) {
    pout->MinEnergy[kk]=6.;
    //Large Cell threshold
    if(kk==5 || kk==6)pout->MinEnergy[kk]=4.;
  };

  Int_t nentries = pout->nentries;
  if (NumberEventsToProcess>0) nentries = NumberEventsToProcess;  

  for (Int_t ievent=0; ievent<nentries; ievent++) {
    pout->GetEntry(ievent);
    p2.Bunchid7bit = pout->Bunchid7bit;

    pout->ClearScratch();
    pout->AllToScratch(false); // make a scratch list of all hard photons -- if the argument were true then the scratch list would also include soft photons.
    pout->ClusterwithYPhi(false);// cluster for angular size

    TObjArray* Clust=pout->ClusterScratch(maxsep);
    TLorentzVector Vscr=pout->SumScratch();
    TIter nxt(Clust);
    TObjArray* ctmp;

    if((ievent%100000)==0) printf("%d Events Read \n", ievent);
    while(ctmp=(TObjArray*) nxt())
    {
      Int_t ClusterIndexOf=Clust->IndexOf(ctmp);

      if (ctmp->GetEntries()>0)
      {
        p2.Z = 1.0;
        p2.N12=ctmp->GetEntries();
        p2.spin=pout->spin;

        TLorentzVector TotalVector = pout->SumList(ctmp);
        if (TotalVector.Pt()==0) continue;

        p2.E12 = TotalVector.E();
        p2.M12 = TotalVector.M();
        //p2.Pttest = TotalVector.Pt();
        p2.Eta = TotalVector.PseudoRapidity();
        p2.Ntracks=pout->scratchlist->GetEntries();

        LVec* tv1=(LVec*) ctmp->First();

        //p2.Phi1=tv1->Phi();
        //p2.Esoft=pout->ClusterSoftE(ClusterIndexOf,.9).E();

        TVector3 uv1=tv1->Vect();
        uv1.SetMag(1.);
        //p2.Y1=0;

        if(tv1->Pt()==0)continue;

        p2.Y1=tv1->PseudoRapidity();
        p2.Phi=TotalVector.Phi();

        //if (ctmp->GetEntries()==2) p2.Pt = TotalVector.Pt();
        p2.Pt = TotalVector.Pt();

        if(p2.Ntracks>p2.N12) // otherwise (i.e. if Ntracks=N12) these should be zero anyway. 
          /* added later: but it turns out that what gets written is not zero but the last value set to the variable; so, its necessary to set the values to 0 explicitly. */
        {
          p2.Ptaway=(Vscr-TotalVector).Pt();
          p2.Phiaway=(Vscr-TotalVector).Phi();
          if(p2.Ptaway==0)continue;
          p2.Yaway=(Vscr-TotalVector).PseudoRapidity();
          p2.Maway=(Vscr-TotalVector).Mag();
        }
        else
        {
          // if these are not set to zero explicitly here then they get the values of what was last assigned to them.
          p2.Ptaway = 0;
          p2.Phiaway = 0;
          p2.Yaway=0;
          p2.Maway=0;
        }
        if (p2.N12==2) {
          LVec* tv2 = (LVec*) ctmp->Last();
          p2.Z = fabs(tv1->E() - tv2->E())/p2.E12;
        }

        TwoTr.Fill(); //Fill the tree
      }
    }   
  };
  TwoTr.Write();	/* write to the file OutFileName that is  open */
  return 0;
};

Int_t AnalTools::OutputMycells(char* filelist, FilesSet* p_files, Bool_t doSingle) {
  // creates a small output file (minimum) from a OFile -- similar to readqM but with much fewer branches in the output file - by Saroj Adhikari - for learning purposes; the output file contains Cell instances for each of the cells which can then be used to create _Cellr%_c%_%.root files needed for calibration analysis; this code will only create an output file with minimum information necessary for calibration and will be expanded later as needed //
  Int_t keyIndex=kSingleKey;
  Float_t maxsep=0.07; 
  Float_t minsep=0.0;
  Float_t ECutMin=60.;
  Float_t ECutMax=90.;
  Float_t ZMax=0.75;

  TFile* out = new TFile(OutFileName, "recreate");
  Qt qt(p_files);
  Geom* p_geom=new Geom(p_files);
  pout = new poutTree(filelist);
  pout->GetEntry(0);
  Int_t RunNum = pout->Rnum;
  p_files->Print(); // necessary to initialize ->path below or else use Path()
  CalibStr gcor(RunNum, (const char*) p_files->p_fpdgaincorr()->path);
  CalibStr gain(RunNum, (const char*) p_files->p_fpdgain()->path);

  Cell* cell[34][17][4];
  TObjArray mcells(2400,0);  
  Bool_t ProcessAll = false;

  for (Int_t i=0; i<4; i++) {
    Int_t up=34;
    if (i>1) up=24;
    for (Int_t j=0; j<34; j++) {
      for (Int_t k=0; k<17; k++) {
        cell[j][k][i] = 0;
        if (j>=up || k>=up/2) continue;
        if ((i<2) && (k<8) && (j>8) && (j<25)) continue;
        if ((i>1) && (k<5) && (j>6) && (j<17)) continue;

        char nam[40];
        sprintf(nam,"Cellr%d_c%d_%d",j,k,i);
        std::cout<<"about to make: "<<nam<<"\n";        
        out->cd();
        Cell* pcell = new Cell(2, i+1, j+1, k+1, RunNum,p_geom, &gain, &gcor);
        pcell->SetName(nam);
        cell[j][k][i] = pcell;
        mcells.Add(pcell);

        sprintf(nam,"CellD%d_c%d_%d",j,k,i);
        pcell->InitTree(nam);
      };
    };
  };

  TH2F* qtHist[4];
  TH2F* qtHistLed[4];
  TFile adcroot("adc.root");
  if(adcroot.IsOpen()) {
    qtHist[0]=(TH2F*) adcroot.GetKey("qtHistLN")->ReadObj();
    qtHist[1]=(TH2F*) adcroot.GetKey("qtHistLS")->ReadObj();
    qtHist[2]=(TH2F*) adcroot.GetKey("qtHistSN")->ReadObj();
    qtHist[3]=(TH2F*) adcroot.GetKey("qtHistSS")->ReadObj();

    qtHistLed[0]=(TH2F*) adcroot.GetKey("qtHistLedLN")->ReadObj();
    qtHistLed[1]=(TH2F*) adcroot.GetKey("qtHistLedLS")->ReadObj();
    qtHistLed[2]=(TH2F*) adcroot.GetKey("qtHistLedSN")->ReadObj();
    qtHistLed[3]=(TH2F*) adcroot.GetKey("qtHistLedSS")->ReadObj();
  }

  remove ("data.txt");
  FILE* fp=fopen("data.txt","a");  
  ProcessAll=true;
  out->cd();

  for(Int_t kk=1;kk<9;kk++)pout->MinEnergy[kk]=4;
  pout->MinEnergy[5]=2.;
  pout->MinEnergy[6]=2.;

  Int_t nentries = pout->nentries;
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;  
  if(NumberEventsToProcess>nentries)nentries=pout->nentries;
  p_geom->FMSGeom=true;
  Float_t EnergyNSTB[4];
  printf("Number of Events to Process: %d \n", nentries);
  for (Int_t ievent=0; ievent<nentries; ievent++) {
    if((ievent%100000)==0) printf("%d Events Read \n", ievent);    
    pout->GetEntry(ievent);
    RunNum = pout->Rnum;
    Float_t mass;
    Float_t energy;
    Float_t Z;
    //printf("here0 \n");
    pout->ClearScratch();
    pout->AllToScratch(false);
    pout->ClusterwithYPhi(false);
    TObjArray* Clust = pout->ClusterScratch(maxsep); //maxseparation = .045 for pi0 reconstruction.
    //TLorentzVector Vscr = pout->SumScratch(); //not used
    TIter nxt(Clust);
    TObjArray* ctmp;
    //printf("here-1 \n");
    while(ctmp=(TObjArray*) nxt()) {
      Int_t ClusterIndexOf = Clust->IndexOf(ctmp);
      /* the claibration works by looking at pi0 masses, which is detected as two photons close together on the detector (since pi0 decays into two photons separated by a small angle from its path of motion); so, we just look for two photons. The small angular separation is characterized by the .045 anglular separation above in pout->ClusterScratch */
      //printf("here \n");
      if (!doSingle && ctmp->GetEntries()==2) {
        LVec* photon[2];
        photon[0] = (LVec*) ctmp->First(); //first photon
        photon[1] = (LVec*) ctmp->After(photon[0]); //other photon
        TLorentzVector twophotons = *(photon[0]) + *(photon[1]);
        mass = twophotons.Mag();
        energy = twophotons.E();
        Float_t y2 = (*photon[0] + *photon[1]).PseudoRapidity();

        // eta selection based on energy and pseudorapidity
        if (energy<ECutMin || energy>ECutMax) continue;
        if (fabs(y2-3.7)>0.2) continue;
        /*
           if (fabs(energy-85.0)<5.0) {
           if (fabs(y2-3.75)>0.15) continue;
           }
           if (fabs(energy-95.0)<5.0) {
           if (fabs(y2-3.75)>0.25) continue;
           }
           if (fabs(energy-110.0)<10.0) {
           if (fabs(y2-3.85)>0.35) continue;
           }
           if (fabs(energy-135.0)<15.0) {
           if (fabs(y2-3.9)>0.3) continue;
           }
           */

        if (photon[0]->Nstb==0) continue;
        if (photon[1]->Nstb==0) continue;
        if (photon[0]->Nstb != photon[1]->Nstb) continue; // remove events crossing boundary
        TVector3 phpos0 = pout->PosInDet(photon[0], p_geom, 0.);
        TVector3 phpos1 = pout->PosInDet(photon[1], p_geom, 0.);
        TVector3 xyab[2];
        xyab[0] = p_geom->LocalXYZ(2, photon[0]->Nstb, phpos0, true);
        xyab[1] = p_geom->LocalXYZ(2, photon[1]->Nstb, phpos1, true);
        Int_t phRow[2], phCol[2];
        phRow[0] = (Int_t) xyab[0].Y();
        phRow[1] = (Int_t) xyab[1].Y();
        phCol[0] = (Int_t) xyab[0].X();
        phCol[1] = (Int_t) xyab[1].X();

        for (Int_t iz=0; iz<2; iz++) {
          /* check if the Row and Column numbers are within the detector range */
          Bool_t good = true;

          Int_t iz2 = 1 - iz;
          if (phRow[iz]<0 || phCol[iz]<0 || phRow[iz2]<0 || phCol[iz2]<0) good = false;
          if (phRow[iz]>=qt.ROW_NUM[photon[iz]->Nstb-1] || phCol[iz]>=qt.COL_NUM[photon[iz]->Nstb-1]) good = false;
          if (phRow[iz2]>=qt.ROW_NUM[photon[iz2]->Nstb-1] || phCol[iz2]>=qt.COL_NUM[photon[iz2]->Nstb-1]) good = false;

          if (photon[iz]->Iew==2) {
            if (photon[iz]->Nstb<3) {
              if (phCol[iz]<8 && phRow[iz]>8 && phRow[iz]<25) good = false;
            }
            else {
              if (phCol[iz]<5 && phRow[iz]>6 && phRow[iz]<17) good = false;
            }
          }   

          if (good) {

            Bool_t FillIt = false; 
            Cell* ocl = cell[phRow[iz2]][phCol[iz2]][photon[iz2]->Nstb-1];
            Cell* tcl = cell[phRow[iz]][phCol[iz]][photon[iz]->Nstb-1];

            if (ocl!=0 && tcl!=0) {

              if ((ocl->Masspeakfraction>.15 && ocl->Massmaxcontents>10) || ProcessAll) {
                FillIt = true;
              }
            }

            if (tcl!=0 && ProcessAll) FillIt = true;

            Z = (fabs(photon[0]->E()-photon[1]->E()))/(photon[0]->E()+photon[1]->E());
            if (Z>ZMax) FillIt = false;

            if (FillIt) {  
              //printf("fillit is true \n");          

              tpCellDat cd;
              cd.Mcell = mass;
              //printf("mass: %f \n", mass);
              cd.Epair = photon[iz]->E() + photon[iz2]->E();
              cd.Ypair = y2;
              cd.E1 = photon[iz]->E();
              cd.E2 = photon[iz2]->E();
              cd.X1 = xyab[iz].X();
              cd.X2 = xyab[iz2].X();
              cd.Y1 = xyab[iz].Y();
              cd.Y2 = xyab[iz2].Y();
              cd.NSTB = photon[iz]->Nstb; // detector number for the first photon
              cd.ievt = pout->ievt;
              cd.Rnum = pout->Rnum;
              cd.nSavedHits = pout->nSavedHits;

              if (cd.nSavedHits > 100) cd.nSavedHits = 100;
              for (int in = 0; in < cd.nSavedHits; in++) {
                cd.SavedHits[in] = pout->SavedHits[in];
              };

              tcl->FillTree(cd);
              char setvalchar[100];
              sprintf(setvalchar,"%d,%d,%d\n",phRow[iz],phCol[iz],cd.NSTB);
              fputs(setvalchar,fp);

              //printf("here filling tree %d %d \n", phRow[iz], phCol[iz]);
            };
          };
        };      
      }
      else {
        if (doSingle && ctmp->GetEntries()==1) {
          // create Output.root for single photon events Ntracks==1 (and then obviously N12==1;
          LVec* photon; 
          photon = (LVec*) ctmp->First(); //first and only photon 	    
          mass = photon->Mag();
          energy = photon->E();
          if (photon->Nstb==0) continue;
          TVector3 phpos = pout->PosInDet(photon, p_geom, 0.);
          TVector3 xyab;
          xyab = p_geom->LocalXYZ(2, photon->Nstb, phpos, true);
          Int_t phRow, phCol;
          phRow = (Int_t) xyab.Y();
          phCol = (Int_t) xyab.X();

          TLorentzVector TotalVector = pout->SumList(ctmp);
          if (TotalVector.Pt()==0) continue;  	      

          Bool_t good = true;
          //check if the row and cols and within the detecor range
          if (phRow<0 || phCol<0 || phRow>=qt.ROW_NUM[photon->Nstb-1] || phCol >= qt.COL_NUM[photon->Nstb-1] ) good = false;

          if (photon->Iew==2) {
            if (photon->Nstb<3) {
              if (phCol<8 && phRow>8 && phRow<25) good = false;
            }
            else {
              if (phCol<5 && phRow>6 && phRow<17) good = false;
            }
          }

          if (good) {
            Bool_t FillIt = false; 
            Cell* ocl = cell[phRow][phCol][photon->Nstb-1];
            Cell* tcl = cell[phRow][phCol][photon->Nstb-1];

            if (ocl!=0 && tcl!=0) {

              if ((ocl->Masspeakfraction>.15 && ocl->Massmaxcontents>10) || ProcessAll) {
                FillIt = true;
              }
            }

            if (tcl!=0 && ProcessAll) FillIt = true;

            if (FillIt) {  
              //printf("fillit is true \n");          

              tpCellDat cd;
              cd.Mcell = mass;
              //printf("mass: %f \n", mass);
              cd.Epair = energy;
              cd.Ypair = photon->PseudoRapidity();
              cd.E1 = photon->E();
              cd.E2 = pout->scratchlist->GetEntries(); //the Ntracks variable set to E2 which is of no use here.
              cd.X1 = xyab.X();
              cd.X2 = TotalVector.Pt();
              cd.Y1 = xyab.Y();
              cd.Y2 = pout->spin;
              cd.NSTB = photon->Nstb; // detector number for the first photon
              cd.ievt = pout->ievt;
              cd.Rnum = pout->Rnum;
              cd.nSavedHits = pout->nSavedHits;

              if (cd.nSavedHits > 100) cd.nSavedHits = 100;
              for (int in = 0; in < cd.nSavedHits; in++) {
                cd.SavedHits[in] = pout->SavedHits[in];
              };

              tcl->FillTree(cd);
              char setvalchar[100];
              sprintf(setvalchar,"%d,%d,%d\n",phRow,phCol,cd.NSTB);
              fputs(setvalchar,fp);

              //printf("here filling tree %d %d \n", phRow[iz], phCol[iz]);
            };
          };	          	      
        }
      }
    };    
  };
  fclose(fp);
  // out of nentries loop
  printf("out of nentries loop");

  TIter next2(&mcells);
  Cell* pcl;
  while(pcl=(Cell*) next2()) {
    Int_t idet=1;
    if(pcl->Iew==2&& ((idet=pcl->Instb)<5)) {
      Int_t chan=(pcl->Row1-1)*qt.COL_NUM[idet-1] + pcl->Col1;
      if(pcl->p_adc)delete pcl->p_adc;
      if(pcl->p_adcLed)delete pcl->p_adcLed;
      char nam[40];
      sprintf(nam,"adc_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
      printf("set adc %s \n",nam);
      pcl->p_adc=qtHist[idet-1]->ProjectionY(nam,chan,chan);
      sprintf(nam,"adcLed_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
      printf("set adcLed %s \n",nam);
      pcl->p_adcLed=qtHistLed[idet-1]->ProjectionY(nam,chan,chan);
    };
  };
  std::cout<<"Finish adc setup\n";
  out->cd();
  for (Int_t j=0; j<4; j++) {
    qtHist[j]->Write();
    qtHist[j]->GetNbinsX(), qtHist[j]->GetNbinsY();
  }

  mcells.Print();
  printf("%d", keyIndex);
  mcells.Write("Mycells", keyIndex);
};

RunDepMgr* AnalTools::OpenRunDepMgr(TString RDPath)
{
  printf("enter OpenRunDepMgr\n");
  if(RDPath=="")
  {
    RDPath=gSystem->Getenv("RunDepPath");//set in SetFMSEnv
  };
  TFile* SaveDefaultFile=gROOT->GetFile();
  TFile* RdepFile=new TFile(RDPath);
  if(Mgr)delete Mgr;
  Mgr=0;
  printf("enter OpenRunDepMgr 1\n");
  if(RdepFile)
  {
    printf("opened %s \n",(const char*) RDPath);
    if(! RdepFile->IsOpen())
    {
      printf("Getting new Mgr from RdepFile\n");
      delete RdepFile;
      return 0;
    };
    Mgr=(RunDepMgr*) RdepFile->Get("Mgr");
    Mgr->RunDepPath=RDPath;
    if(RdepFile)delete RdepFile;
    if(SaveDefaultFile)SaveDefaultFile->cd();
    return Mgr;
  };
  if(SaveDefaultFile)SaveDefaultFile->cd();
  return 0;
};
void AnalTools::PrintCell(Cell* c ,Int_t np,Int_t nstb1,TString psname)
{
  c->InitBranches();
  int nentries=c->CellD->GetEntries();
  TCanvas* cc1=new TCanvas("cc1","cc1",500,800);
  cc1->cd();
  cc1->Print(psname+"(");
  printf("Printcell: nentries=%d\n",nentries);
  if(np<nentries)nentries=np;
  for(int j=0;j<nentries;j++)
  {
    c->CellD->GetEntry(j);
    printf("PrintCell %d hits \n",j);
    ((TMatrixT<float>) c->FillFMSADC(nstb1)).Draw("box");
    Float_t sumadc=((TMatrixT<float>) c->FillFMSADC(nstb1)).Sum();;
    char tit[200];
    sprintf(tit,"%s Index=%d Sum=%d \n",(const char*) psname,j,(int) sumadc);
    cc1->SetTitle(tit);
    cc1->Print(psname);
  };
  cc1->Clear();
  cc1->Print(psname+")");
  delete cc1;
};
void AnalTools::Extract1(TObjArray* mycells , TString cellname,TString dirname="data2")
{
  TFile* SaveDefaultFile=gROOT->GetFile();
  TString fname=dirname;
  fname=fname+"/_"+cellname+".root";
  cout<<"to filename="<<fname<<"\n";
  TFile* flout=new TFile(fname,"recreate");
  Cell* c=(Cell*) mycells->FindObject(cellname);
  c->Write(cellname);
  flout->Close();
  if(SaveDefaultFile)SaveDefaultFile->cd();
  delete flout;
}
bool AnalTools::TestCellTree(Cell* c)
{
  if(!c)return false;
  if(!c->CellD)return false;
  c->InitBranches();
  int entries=c->CellD->GetEntries();
  int cntrep=0;
  int adclast=0;
  for(int j=0;j<entries;j++)
  {
    c->CellD->GetEntry(j);
    int adcv=c->FillFMSADC(c->Instb).Sum();
    if(adcv>10&&adcv==adclast)cntrep++;
    adclast=adcv;
  };
  if(cntrep<entries/30 || cntrep<20)return true;
  return false;
}

