#include "Asymmetry.h"
using namespace std;
ClassImp(Asymmetry);
Asymmetry::Asymmetry(TH2F* h,TString name, Float_t norm,TString titleString)
{
  TitleString=titleString;
  TitleString.Replace(0,13,"");
  printf("start\n");
  
  h2dAsm=new TH2F(*h);
  printf("step 1\n");
  
  h2dAsm->SetName(name);
  printf("step 1\n");
  
  TString sn[8]={"N_yu_bu","N_yd_bu","N_yu_bd","N_yd_bd","S_yu_bu","S_yd_bu","S_yu_bd","S_yd_bd"};
  TString sn2[4]={"yu_bu","yd_bu","yu_bd","yd_bd"};
  printf("got here\n");
  
  for(int j=0;j<8;j++)
    {
      printf("j=%d \n",j);
      TCanvas* ctmp=(TCanvas*) gROOT->GetListOfCanvases()->First();
      if(ctmp){ctmp->cd();}
      else
	ctmp=new TCanvas("ctmp","ctmp",600,600);
      
      SpinNames[j]=sn[j];
      SpinLoc[j]=j+1;
      cout<<SpinNames[j]<<" "<<SpinLoc[j]<<"\n ";
      
	SpinH[j]=new TH1D(*( h2dAsm->ProjectionX(SpinNames[j],SpinLoc[j],SpinLoc[j])));

	
	SpinH[j]->SetName(name+"_"+SpinNames[j]);
	SpinH[j]->SetTitle(name+"_"+SpinNames[j]);
	SpinH[j]->Sumw2();
      };
    for(int j=0;j<4;j++)
      {
	NSsumNames[j]=sn2[j];
	NSSpinH[j]=new TH1D((*SpinH[j])+(*SpinH[j+4]));
	NSSpinH[j]->SetName(NSsumNames[j]);
	NSSpinH[j]->SetTitle(NSsumNames[j]);
      };

    Norm=1.;
    BN[0]=new TH1D( (*SpinH[0])+(*SpinH[1]));
    BN[0]->SetName("BNUp");
    BN[1]=new TH1D( (*SpinH[2])+(*SpinH[3]));
    BN[1]->SetName("BNDn");
    
    BS[0]=new TH1D( (*SpinH[4])+(*SpinH[5]));
    BS[0]->SetName("BNUp");
    BS[1]=new TH1D( (*SpinH[6])+(*SpinH[7]));
    BS[1]->SetName("BNDn");
    
    YN[0]=new TH1D( (*SpinH[0])+(*SpinH[2]));
    YN[0]->SetName("BNUp");
    YN[1]=new TH1D( (*SpinH[1])+(*SpinH[3]));
    YN[1]->SetName("BNDn");
    
    YS[0]=new TH1D( (*SpinH[4])+(*SpinH[6]));
    YS[0]->SetName("BNUp");
    YS[1]=new TH1D( (*SpinH[5])+(*SpinH[7]));
    YS[1]->SetName("BNDn");
    BNorth=new TH1D((*BN[0])+(*BN[1]));
    BNorth->SetName("BNorth");

    BSouth=new TH1D((*BS[0])+(*BS[1]));
    BSouth->SetName("BSouth");
    YNorth=new TH1D((*YN[0])+(*YN[1]));
    YNorth->SetName("YNorth");
    YSouth=new TH1D((*YS[0])+(*YS[1]));
    YSouth->SetName("YSouth");
    North=new TH1D((*BNorth)+(*YNorth));
    South=new TH1D((*BSouth)+(*YSouth));
    All=new TH1D((*North)+(*South));
    South->SetName("South");
    North->SetName("North");
    All->SetName("All");
    AYN=new TH1D(((*YN[0])-(*YN[1]))/((*YN[0])+(*YN[1])));
    AYS=new TH1D(((*YS[0])-(*YS[1]))/((*YS[0])+(*YS[1])));
    ABN=new TH1D(((*BN[0])-(*BN[1]))/((*BN[0])+(*BN[1])));
    ABS=new TH1D(((*BS[0])-(*BS[1]))/((*BS[0])+(*BS[1])));

    AS_xx_yx=new TH1D(((*SpinH[4])+(*SpinH[7])-((*SpinH[5])+(*SpinH[6])))/(*South));
    AN_xx_yx=new TH1D(((*SpinH[0])+(*SpinH[3])-((*SpinH[1])+(*SpinH[2])))/(*North));
    AS_xx_yx->SetName("AS_xx_yx");
    AS_xx_yx->SetTitle("ASouth uu +dd-ud-du");
    AN_xx_yx->SetName("AN_xx_yx");
    AN_xx_yx->SetTitle("ANorth uu +dd-ud-du");

    AS_xx_yxp=new TH1D(((*SpinH[4])-(*SpinH[7])+((*SpinH[5])-(*SpinH[6])))/(*South));
    AN_xx_yxp=new TH1D(((*SpinH[0])-(*SpinH[3])+((*SpinH[1])-(*SpinH[2])))/(*North));


    AS_xx_yxp->SetName("AS_xx_yx");
    AS_xx_yxp->SetTitle("ASouth uu -dd-ud+du");
    AN_xx_yxp->SetName("AN_xx_yx");
    AN_xx_yxp->SetTitle("ANorth uu -dd-ud+du");
    AS_pp_mm=new TH1D(((*SpinH[4])-(*SpinH[7]))/((*SpinH[4])+(*SpinH[7])));
    AN_pp_mm=new TH1D(((*SpinH[0])-(*SpinH[3]))/((*SpinH[0])+(*SpinH[3])));

    CrossY=0;
    CrossB=0;
    TH1D* tmp=CrossAN("yellow");
    tmp=CrossAN("blue");
    for(int m=0;m<4;m++)
      {
	ANS[m]=new TH1D( (*(SpinH[m])-(*(SpinH[m+4])))/ (*(SpinH[m])+(*(SpinH[m+4]))) );
      }
    ReNorm(norm);
};


void Asymmetry::ReNorm(Float_t norm)
{
  AYN->Scale(norm);
  AYS->Scale(norm);
  ABN->Scale(norm);
  ABS->Scale(norm);
  CrossY->Scale(norm);
  CrossB->Scale(norm);
  AN_xx_yx->Scale(norm);
  AS_xx_yx->Scale(norm);
  AN_xx_yxp->Scale(norm);
  AS_xx_yxp->Scale(norm);
  AN_pp_mm->Scale(norm);
  AS_pp_mm->Scale(norm);


  Norm=Norm*norm;
}
Asymmetry::~Asymmetry()
{
  for(int j=0;j<8;j++)
    {
      if(SpinH[j])delete SpinH[j];
      if(j<4)
	{
	  if(NSSpinH[j])delete NSSpinH[j];
	};
      if(j<2)
	{
	  if(BN[j])delete BN[j];
	  if(BS[j])delete BS[j];
	  if(YN[j])delete YN[j];
	  if(YS[j])delete YS[j];
	};
      if(AYN)delete AYN;
      if(AYS)delete AYS;
      if(ABN)delete ABN;
      if(ABS)delete ABS;
      if(CrossY)delete CrossY;
      if(CrossB)delete CrossB;
      if(AS_xx_yx)delete AS_xx_yx;      
      if(AN_xx_yx)delete AN_xx_yx;
      if(AS_xx_yxp)delete AS_xx_yxp;      
      if(AN_xx_yxp)delete AN_xx_yxp;
      if(AS_pp_mm)delete AS_pp_mm;      
      if(AN_pp_mm)delete AN_pp_mm;
      for(int j=0;j,4;j++){if(ANS[j]){delete ANS[j];};};
    };
};
TGraphErrors Asymmetry::ValueDist(TH1D h1,TString name)
{
  Double_t mx=h1.GetMaximum();
  Double_t mn=h1.GetMinimum();
  mx=1.1 * mx;
  mn=.9*mn;
  if(mx<0)mx=.8*mx;
  if(mn<0)mn=1.2*mn;

  TH1D vh(name,name,100,mn,mx);
  int nb=h1.GetNbinsX();
  for(int j=0;j<nb;j++)
    {
      Double_t bc=h1.GetBinContent(j);
      if(bc!=0)vh.Fill(bc);
    };
  TGraphErrors gr;
  Int_t cnt=0;
  for(int j=0;j<vh.GetNbinsX();j++)
    {
      Double_t bx,by,be;
      by=vh.GetBinContent(j);
      bx=vh.GetBinCenter(j);
      be=sqrt(by);
      if(by!=0)
	{
	  gr.SetPoint(cnt,bx,by);
	  gr.SetPointError(cnt++,0,be);
	};
    }
  return gr;
}
TH1D Asymmetry::ErrorDist(TH1D h1,TString name)
{
  TH1D vh(name,name,100,h1.GetMinimum(),h1.GetMaximum());
  int nb=h1.GetNbinsX();
  for(int j=0;j<nb;j++)
    {
      Double_t bc=h1.GetBinContent(j);
      if(bc!=0)vh.Fill(h1.GetBinError(j));
    };
  return vh;
}
TH1D* Asymmetry::CrossAN(TString beam)
{
  TH1D* btmp[4];
  TH1D* cross=0;
  if(beam=="yellow")
    {
      //nu,nd,su,sd
      btmp[0]=YS[0];
      btmp[1]=YS[1];
      btmp[2]=YN[0];
      btmp[3]=YN[1];
      if(CrossY)delete CrossY;
      CrossY=new TH1D("CrossY","CrossY",NbinX(),LowX(),HiX());
      cross=CrossY;
    }
  else if(beam=="blue")
    {
      btmp[0]=BS[0];
      btmp[1]=BS[1];
      btmp[2]=BN[0];
      btmp[3]=BN[1];
      if(CrossB)delete CrossB;
      CrossB=new TH1D("CrossB","CrossB",NbinX(),LowX(),HiX());
      cross=CrossB;
    }
  else {return 0;};
  Float_t nuS=0;
  Float_t nuN=0;
  Float_t ndS=0;
  Float_t ndN=0;

  for(Int_t i=1;i<NbinX()+1;i++){
    nuS=TMath::Max(btmp[0]->GetBinContent(i),0.);
    ndS=TMath::Max(btmp[1]->GetBinContent(i),0.);
    nuN=TMath::Max(btmp[2]->GetBinContent(i),0.);
    ndN=TMath::Max(btmp[3]->GetBinContent(i),0.);

    Float_t inverr2=nuS+ndS+nuN+ndN;
    Float_t value=0;
    Float_t denom=sqrt(nuS*ndN)+sqrt(ndS*nuN);
    Float_t err=0;

    if(inverr2>0 && denom>0){
      value=(sqrt(nuS*ndN)-sqrt(ndS*nuN))/denom;
      err=1/sqrt(inverr2);
    };
    cross->SetBinContent(i,value);
    cross->SetBinError(i,err);
  };
  cross->Scale(Norm);
 return cross;
};
void Asymmetry::PlotYBNS(TCanvas* cc,TString XAxisName,Float_t xmin,Float_t xmax)
{
  if(cc==0)return;
  cc->cd();
  cc->Clear();
  North->SetLineColor(2);
  North->GetXaxis()->SetTitle(XAxisName);			      
  North->GetXaxis()->SetRangeUser(xmin,xmax);
  North->SetTitle(TitleString);
  North->Draw();
  South->Draw();
  cc->Print("cnorthsouth.gif");
  cc->Clear();
  CrossB->SetLineColor(4);
  CrossB->Draw();
  cc->GetListOfPrimitives()->ls();
  CrossB->GetYaxis()->SetTitle("AN corrected");
  CrossB->GetYaxis()->SetRangeUser(-.2,.2);
  CrossB->GetYaxis()->SetTitleOffset(1.3);
  CrossB->GetXaxis()->SetRangeUser(xmin,xmax);
  
  CrossB->GetXaxis()->SetTitle(XAxisName);
  CrossY->SetLineColor(5);
  CrossY->Draw("same");
  cc->cd(2);
  ABN->SetLineColor(2);
  ABN->Draw("same");
  ABS->Draw("same");
  cc->Update();
  TFrame* fr=(TFrame*) cc->GetListOfPrimitives()->FindObject("TFrame");
  if(fr)
    {
      fr->SetFillColor(27);
      printf("fill color set to 27\n");
    };
  TPaveText* titl=(TPaveText*) cc->GetListOfPrimitives()->FindObject("title");
  if(titl)
    {
      titl->Clear();

      titl->AddText(TitleString);
      titl->AddText(TitleString);
      titl ->AddText("Blue: CrossRatio blue;  Yellow: CrossRatio yellow");
      titl ->AddText("Red: AN (blue-north) ;  Black: AN (blue-south)");

    };
  cc->Update();
  cc->Print("cYBNS.gif");
}
