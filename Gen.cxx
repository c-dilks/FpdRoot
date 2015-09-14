#include "Gen.h"
ClassImp(Gen);
Gen::Gen(Int_t IEW, Int_t NSTB, Int_t nev,
	 Geom* p_Geom,CalibStr* RGain,CalibStr* RGcor,Int_t seed=12345)
{
  CallFillLimitE=10.;
  g_out=new TTree("g_out","g_out");
  g_out->Branch("Iew",&Iew,"Iew/I");
  g_out->Branch("Instb",&Instb,"Instb/I");
  g_out->Branch("pimass",&pimass,"pimass/F");
  g_out->Branch("Eph1",&Eph1,"Eph1/F");
  g_out->Branch("Eph2",&Eph2,"Eph2/F");
  g_out->Branch("piPt",&piPt,"piPt/F");
  g_out->Branch("piEta",&piEta,"piEta/F");
  g_out->Branch("EnergySum",&EnergySum,"EnergySum/F");
  g_out->Branch("inFiducial",&inFiducial,"inFiducial/B");
  g_out->Branch("piPhi",&piPhi,"piPhi/B");
  g_out->Branch("genEnergy",&Energy,"genEnergy/D");
  g_out->Branch("genpt",&pt,"genpt/D");
  g_out->Branch("genZ",&Z,"genZ/F");
  g_out->Branch("genphi",&phi,"genphi/F");
  g_out->Branch("NPh",&NPh,"NPh/I");
  g_out->Branch("genY",&Y,"genY/D");
  g_out->Branch("genmass",&mass,"genmass/F");
  g_out->Branch("PreDigitizedEnergy",&PreDigitizedEnergy,"PreDigitizedEnergy/F");

  

  GenCnt=0;
  NominalEta=4.;
  NominalPhi=0.;
  p_pi=0;
  pGeom=p_Geom;
  Rgain=RGain;
  Rgaincorr=RGcor;
  Nevents=nev;
  Iew=IEW;
  Instb=NSTB;
  Zlimit=.85;
  DoRecon=false;
  p_gam1=0;
  p_gam2=0;
  Error=0;
  Eratch1=new TH2F("Eratch1","Eratch1",49,0,49.,140,0.,1.4);
  Eratch1F=new TH2F("Eratch1F","Eratch1F",49,0,49.,140,0.,1.4);
  Eratch2=new TH2F("Eratch2","Eratch2",49,0,49.,140,0.,1.4);
  EratchP=new TH2F("EratchP","EratchP",49,0,49.,140,0.,1.4);
  Eratch40=new TH2F("Eratch40","Eratch40",49,0,49.,140,0.,1.4);
  
  Episum=new TH2F("Episum","EpiSum",49,0,49.,120,0.,120.);
  Epigsum=new TH2F("Epigsum","EpigSum",49,0,49.,120,0.,120.);
  Epipsum=new TH2F("Epipsum","EpipSum",49,0,49.,120,0.,120.);
  EpipRsum=new TH2F("EpipRsum","EpipRSum",49,0,49.,120,0.,120.);
  EpipbRsum=new TH2F("EpipbRsum","EpipbRSum",49,0,49.,120,0.,120.);

  GenEvsY=new TH2F("GenEvsY","GenEvsY",100,3.5,4.8,100,0.,100);
  GenpEvsY=new TH2F("GenpEvsy","GenpEvsY",100,3.5,4.8,100,0.,100);


  GenEvsPt=new TH2F("GenEvsPt","GenEvsPt",180,0.,6.,100,0.,100);
  GenPEvsPt=new TH2F("GenPEvsPt","GenPEvsPt",180,0.,6.,100,0.,100);
  GenThvsPh=new TH2F("GenThvsPh","GenThvsPh",64,-3.14,3.14,100,0.,.1);
  GenPThvsPh=new TH2F("GenPThvsPh","GenPThvsPh",64,-3.14,3.14,100,0.,.1);
  GenYvsPtWPt=new TH2F("GenYvsPtWPt","GenYvsPtWPt",180,0,6,150,2.5,7.);
  GenPYvsPtWPt=new TH2F("GenPYvsPtWPt","GenPYvsPtWPt",180,0,6,150,2.5,7.);

  GenPPEvsPt=new TH2F("GenPPEvsPt","GenPPEvsPt",180,0.,6.,100,0.,100);
  GenPPThvsPh=new TH2F("GenPPThvsPh","GenPPThvsPh",64,-3.14,3.14,100,0.,.1);
  GenPPYvsPtWPt=new TH2F("GenPPYvsPtWPt","GenPPYvsPtWPt",180,0,6,150,2.5,7);

  GenZbeam=new TH1F("GenZbeam","GenZbeam",200,-200,200);
  GenPZbeam=new TH1F("GenPZbeam","GenPZbeam",200,-200,200);
  phi1=new TH1F("phi1","phi1",64,-6.28,6.28);
  phi2=new TH1F("phi2","phi2",64,-6.28,6.28);
  phi3=new TH1F("phi3","phi3",64,-6.28,6.28);
  NextEvt=new TMatrix(*(RGain->tm(Iew,Instb)));
  (*NextEvt)=0.;
  Digitized=new TMatrix(*(NextEvt));
  UnDigitized=new TMatrix(*(NextEvt));
  rand=new TRandom(seed);
  
  Nrows=NextEvt->GetNrows();
  Ncols=NextEvt->GetNcols();
  if(Iew==1)
    {
      EPGen=new TF2("EPGen",EPdist,3.2,4.8,.5,5.,2);
      PtGen=new TF1("PtGen",Ptdist,.5,5.,1);
      XGen=new TF1("XGen",Xdist,.3,.95,1);
      NominalEta=4.0;
      NominalPhi=0.;
      if(Instb==2)NominalPhi=TMath::Pi();
    }
  else
    {
      printf("failed to initialize");
      Error=0;
      return;
    };
  EPGen->SetParameter(0,5);
  PtGen->SetParameter(0,6);
  EPGen->SetParameter(1,6);
  hEgen =new TH1F("hEgen","hEgen",100,0,100);
  hEpass =new TH1F("hEpass","hEpass",100,0,100);
  Mass=new TH1F("Mass","Mass",120,0.,1.2);
  Nphotons=new TH1F("Nphotons","Nphotons",20,0.,20.);
  Double_t Pi_d=TMath::Pi();
  for(Int_t ir=0;ir<Nrows;ir++)
    {
      for(Int_t ic=0;ic<Ncols;ic++)
	{
	  char nam[30];
	  sprintf(nam,"DEc%d_r%d",ic+1,ir+1);
	  Exy[ir][ic]=new TH1F(nam,nam,256,0,256.);
	  sprintf(nam,"EphXY%d_%d",ic+1,ir+1);
	  EphXY[ir][ic]=new TH1F(nam,nam,100,0,100.);
	  sprintf(nam,"Z_ADCc%d_r%d",ic+1,ir+1);
	  ZvsADC[ir][ic]=new TH2F(nam,nam,256,0,256.,10,0,1.);
	  sprintf(nam,"PhiADCc%d_r%d",ic+1,ir+1);
	  PhivsADC[ir][ic]=new TH2F(nam,nam,256,0,256.,32,-Pi_d,Pi_d);

	  sprintf(nam,"LocalPosc%d_r%d",ic+1,ir+1);
	  LocalPos[ir][ic]=new TH2F(nam,nam,26,-3.,10.,26,-3.,10.);

	  sprintf(nam,"LocalPos125c%d_r%d",ic+1,ir+1);
	  LocalPos125[ir][ic]=new TH2F(nam,nam,26,-3.,10.,26,-3.,10.);

	  sprintf(nam,"Mcellc%d_r%d",ic+1,ir+1);
	  Mcell[ir][ic]=new TH1F(nam,nam,80,0,.8);

	  sprintf(nam,"NphotvsEc%d_r%d",ic+1,ir+1);
	  NphotvsE[ir][ic]=new TH2F(nam,nam,20,0.,100.,5,0,5);

	  sprintf(nam,"phiZ40c%d_r%d",ic+1,ir+1);
	  Z_Phi_40[ir][ic]=new TH2F(nam,nam,10,-TMath::Pi(),TMath::Pi(),10,0.,1.);
	};
    };
  fitter=new FitTower(NextEvt,pGeom,Iew,Instb);
  pGGams=&(fitter->GGams);
  p_fcnSS=fitter->we.fcnSS;
  TowerWidth=*(pGeom->FpdTowWid(Iew,Instb));
 
};
Float_t Gen::Fill()
{
  Float_t return_val;
  Float_t EnergyDeposited=0.;
  GenCnt++;

  if(rec){delete rec;rec=0;};  
  (*NextEvt)=0.;
  (*Digitized)=0.;
  (*UnDigitized)=0.;
  Bool_t GoodGen=false;
  Double_t px=0;
  Double_t py=0;
  Double_t pz=0;
  mass=.135;

  while(!GoodGen)
    {
      phi=.8*(rand->Rndm()*2-1.)+NominalPhi;
      BeamZ=rand->Gaus(0.,60.);
      GenZbeam->Fill(BeamZ);
      if(rand->Rndm()<.38)mass=.55;
      // EPGen->GetRandom2(Y,pt);
      Double_t rndv=1.;
      while(rndv>.1)
	{
	  pt=PtGen->GetRandom();
	  rndv=rand->Rndm();
	};
      rndv=1.;
      while(rndv>.1)
	{
	  xF=XGen->GetRandom();
	  rndv=rand->Rndm();
	};
      
      Energy=xF*100.;
      pz=sqrt(Energy*Energy-mass*mass-pt*pt);
      if(Iew==2)pz=-pz;
      px=pt*cos(phi);
      py=pt*sin(phi);
      TVector3 pgen(px,py,pz);

      Y=pgen.PseudoRapidity();
      GoodGen=true;
    };
  /*
  Double_t tanh=TMath::TanH(Y);
  Energy=sqrt((mass*mass+pt*pt)/(1-tanh*tanh));
  */
  GenEvsY->Fill(Y,Energy);
  GenEvsPt->Fill(pt,Energy);
  if(p_pi)delete p_pi;
  p_pi=new TLorentzVector(px,py,pz,Energy);
  EnergyDeposited=Energy;
  GenThvsPh->Fill(p_pi->Phi(),p_pi->Theta());
  GenYvsPtWPt->Fill(p_pi->Pt(),p_pi->Eta());
  TVector3 Sab=p_pi->Vect();
  Sab=Sab*(1./Sab.Z()*(*pGeom->ZFPD(Iew,Instb) - BeamZ));
  Sab.SetZ(Sab.Z()+BeamZ); // This is now the vector from the origin to the hit location
  TVector3 vab=pGeom->LocalXYZ(Iew,Instb,Sab,true); 
  //vab is now that vector in the Chamber Local Frame
  
  if(vab.X()<-5. ||vab.X()>12)return 0.;
  if(vab.Y()<-5. ||vab.Y()>12)return 0.;
  //  if(fabs(vab.X()-3.5)>.5)return 0.;
  //  if(fabs(vab.Y()-.5)>.5)return 0.;
  
  piInside=true;
  if(vab.X()<0. ||vab.X()>7)piInside=false;
  if(vab.Y()<0. ||vab.Y()>7)piInside=false;
  
  Float_t xc=-1;
  Float_t xa=-1;
  Float_t xb=-1;
  if(piInside)xc=((Int_t) vab.Y())*7.0+vab.X();
  Float_t cs=2*rand->Rndm()-1;
  theta_s=acos(cs);
  phi_s=2*(rand->Rndm()-.5)*TMath::Pi();
  TVector3 pstar;
  pstar.SetMagThetaPhi(mass/2,theta_s,phi_s);
  TVector3 boost=p_pi->BoostVector();
  if(p_gam1)delete p_gam1;
  if(p_gam2)delete p_gam2;
  p_gam1=new TLorentzVector(pstar,mass/2.);
  p_gam2=new TLorentzVector(-1*pstar,mass/2.);
  
  p_gam1->Boost(boost);
  p_gam2->Boost(boost);
  Epigsum->Fill(xc,p_pi->E());


  TVector3 Sa=p_gam1->Vect();
  TVector3 Sb=p_gam2->Vect();
  Sa=Sa*(1./Sa.Z()*(*pGeom->ZFPD(Iew,Instb)- BeamZ));
  Sa.SetZ(Sa.Z()+BeamZ);
  Sb=Sb*(1./Sb.Z()*(*pGeom->ZFPD(Iew,Instb)- BeamZ));
  Sb.SetZ(Sb.Z()+BeamZ);
  TVector3 va=pGeom->LocalXYZ(Iew,Instb,Sa,true);  
  TVector3 vb=pGeom->LocalXYZ(Iew,Instb,Sb,true);
  LocalGamPos[0]=va;
  LocalGamPos[1]=vb;
  inFiducial=piInside;
  phi1->Fill(phi_s);
  if(va.X()<0.5 || va.X()>6.5)inFiducial=false;
  if(va.Y()<0.5 || va.Y()>6.5)inFiducial=false;
  if(vb.X()<0.5 || vb.X()>6.5)inFiducial=false;
  if(vb.Y()<0.5 || vb.Y()>6.5)inFiducial=false;

  Int_t nphots_in=0;
  if(va.X()>0 && va.X()<7 && va.Y()>0 && va.Y()<7)
    {  
      nphots_in++;
      EnergyFromPhoton(&va,p_gam1->E());
      xa=((Int_t) va.Y())*7.0+va.X();
    };
  
  if(vb.X()>0 && vb.X()<7 && vb.Y()>0 && vb.Y()<7)  
    {
      nphots_in++;
      EnergyFromPhoton(&vb,p_gam2->E());
      xb=((Int_t) vb.Y())*7.0+vb.X();
    };

  if(nphots_in<1)return -1.;
  phi2->Fill(phi_s);

  Z=fabs((p_gam2->E()-p_gam1->E())/(p_gam2->E()+p_gam1->E()));
  if(Z>Zlimit)return -1;
  if(nphots_in==2)Epipsum->Fill(xc,p_pi->E());  
  GenPThvsPh->Fill(p_pi->Phi(),p_pi->Theta());
  GenPEvsPt->Fill(p_pi->Pt(),p_pi->E());
  GenPYvsPtWPt->Fill(piPt=p_pi->Pt(),piEta=p_pi->Eta());
  piPhi=p_pi->Phi();
  
  Float_t ADCsum=DigitizeEnergy()->Sum();
  EnergySum=NextEvt->Sum();
  if(EnergySum>1)hEgen->Fill(EnergySum);
  if(ADCsum<125)return -1.;

  GenPPThvsPh->Fill(p_pi->Phi(),p_pi->Theta());
  GenPPEvsPt->Fill(p_pi->Pt(),p_pi->E());
  GenPPYvsPtWPt->Fill(p_pi->Pt(),p_pi->Eta());
  GenPZbeam->Fill(BeamZ);

  Epipsum->Fill(xc,EnergyDeposited=UnDigitized->Sum());
  GenpEvsY->Fill(Y,EnergyDeposited=UnDigitized->Sum());
  hEpass->Fill(EnergySum);
  for(Int_t ir=0;ir<Nrows;ir++)
    {
      for(Int_t ic=0;ic<Ncols;ic++)
	{
	  if((*NextEvt)(ir,ic)>.01) 
	    {
	      Double_t AdC=(*Digitized)(ir,ic);
	      Exy[ir][ic]->Fill(AdC);
	      ZvsADC[ir][ic]->Fill(AdC,Z);
	      PhivsADC[ir][ic]->Fill(AdC,phi_s);
	      if(AdC>10){
		//		printf("ir=%d ic=%d AdC=%f \n",ir,ic,AdC);
		LocalPos[ir][ic]->Fill(LocalGamPos[0].X(),LocalGamPos[0].Y());
		LocalPos[ir][ic]->Fill(LocalGamPos[1].X(),LocalGamPos[1].Y());
		if(AdC>125){
		LocalPos125[ir][ic]->Fill(LocalGamPos[0].X(),LocalGamPos[0].Y());
		LocalPos125[ir][ic]->Fill(LocalGamPos[1].X(),LocalGamPos[1].Y());
		};
	      };
	    };
	};
      
    };
  EpipbRsum->Fill(xc,UnDigitized->Sum());
  if(piInside)Eratch1->Fill( xc,(UnDigitized->Sum())/p_pi->E());

  if(DoRecon&& EnergySum>40.&& piInside)
    {
      rec=new Yiqun(UnDigitized,pGeom,Rgain,Rgaincorr,Iew,Instb);
      NPh=rec->NPh;
      Nphotons->Fill(NPh);
      Double_t lgwidth;
      rec->fitter->GetTWidthCM(lgwidth);
      Int_t iHiE=-1;
      Int_t HiRow=0;
      Int_t HiCol=0;
      Float_t HiE=0.;
      for(Int_t j=0;j<NPh;j++)
	{
	  Float_t xph=rec->photons[j].xPos/lgwidth;
	  Float_t yph=rec->photons[j].yPos/lgwidth;	  

	  Int_t rph=(Int_t) yph;
	  Int_t cph=(Int_t) xph;
	  if(rph>=0 && rph<7 && cph>=0 && cph<7)
	    {
	      NphotvsE[rph][cph]->Fill(EnergySum,NPh);
	      if(rec->mom(j).E()>HiE)
		{
		  HiE=rec->mom(j).E();
		  iHiE=j;
		  HiRow=rph;
		  HiCol=cph;
		};
	      
	    };
	  
	};
      if(iHiE>=0)
	{
	  EphXY[HiRow][HiCol]->Fill(HiE);
	};
      
      phi3->Fill(phi_s);
      Int_t ia=((Int_t) xa);
      Int_t ib=((Int_t) xb);
      Int_t ic=((Int_t) xc);
      Int_t rw1=ia/7;
      Int_t cl1=ia%7;
      Int_t rw2=ib/7;
      Int_t cl2=ib%7;
      Int_t rw3=ic/7;
      Int_t cl3=ic%7;
      pimass=0;
      Eph1=0;
      Eph2=0;
      if(rec->NPh==1)Eph1=rec->mom(0).E();
      if(rec->NPh==2) 
	
	{
	  Zrec=1.;
          Eph1=rec->mom(0).E();
          Eph2=rec->mom(1).E();
	  Float_t Esum=(rec->mom(0)+rec->mom(1)).E();
	  Float_t Ediff=(rec->mom(0)-rec->mom(1)).E();
	  if(Esum>0)Zrec=fabs(Ediff/Esum);
	  Mass->Fill( pimass=(rec->mom(0)+rec->mom(1)).Mag()  );
	  if(fabs(pimass-.135)<.1&& Zrec<.95)EpipRsum->Fill(xc,UnDigitized->Sum());
	  
	  if(inFiducial)
	    {
	      Eratch2->Fill(xc,(EnergyDeposited=(rec->mom(0)+rec->mom(1)).E())
			    /p_pi->E());
	      Float_t ratio=(p_gam1->E()/(p_pi->E()));
	      EratchP->Fill(xa,(EnergyDeposited)
			    /p_pi->E(),ratio);
	      EratchP->Fill(xb,(EnergyDeposited)
			    /p_pi->E(),1-ratio);
	      Eratch1F->Fill( xc,(UnDigitized->Sum())/p_pi->E());

	    };
	  if(Esum>40 && Zrec<.95)
	    {
	      if(rw1>=0 &&rw1<7 &&cl1>=0 && cl1<7)Mcell[rw1][cl1]->Fill(pimass);
	      
	      if(rw2>=0 &&rw2<7 &&cl2>=0 && cl2<7)Mcell[rw2][cl2]->Fill(pimass);
	    };
	  if(p_pi->E()>38)
	    {
	      Z_Phi_40[rw3][cl3]->Fill(phi_s,Z);
	      Eratch2->Fill(xc,(EnergyDeposited=(rec->mom(0)+rec->mom(1)).E())
			    /p_pi->E());
	      
	    };
	      
	  
	  
	};
      g_out->Fill();

    };
  return EnergyDeposited;
};
Double_t Gen::Ptdist(Double_t *x, Double_t *par)
{
  Double_t pt=x[0];
  return TMath::Power(pt,-par[0])*pt; // extra pt for invarient xsec
};
Double_t Gen::Xdist(Double_t *x, Double_t *par)
{
  Double_t Xf=x[0];
  return TMath::Power((1-Xf),par[0])/Xf;//extra Xf for Jacobian factor
};


Double_t Gen::EPdist(Double_t *x,Double_t *par)
{
  // funciton of y and pt
  // returns f(y,pt) where dN=(dpt)(dy)*f(y,pt)
  // dN/(pt dpt)/dy = (1-xf)^par[0]*(pt)^(-par[1])
  //

  Double_t y=x[0];
  Double_t pt=x[1];
  Double_t thy=TMath::TanH(y);
  Double_t yfudge=1.+2*(fabs(y)-3.95)*(fabs(y)-3.95)/.25/.25;
  Double_t Supress=1.;
  Double_t E=sqrt((pt*pt+.135*.135)/(1-thy*thy));
  if(E<20)
    {
      if(E<19)return 0.;
      Float_t dE=E-19.; 
      E=20.;
      Supress=dE;      
    };
  if(E>99.9)return 0.;

  Float_t xf=E/100.;

  // Generate with powers par[0] and par[1] multiplied by invariant phase space
  return TMath::Power(1-xf,par[0])/TMath::Power(pt,par[1])*(pt)*Supress*yfudge;

};
TMatrix* Gen::DigitizeEnergy()
{
  (*Digitized)=0.;
  for(Int_t nr=0;nr<Nrows;nr++)
    {
      for(Int_t nc=0;nc<Ncols;nc++)

	{
	  Float_t factor=Rgain->GetValue(Iew,Instb,nr,nc);
	  factor*=Rgaincorr->GetValue(Iew,Instb,nr,nc);
	  if(factor>0) (*Digitized)(nr,nc) =((Float_t)  ((Int_t) 
							 ((*NextEvt)(nr,nc)/factor) 
							 ));//+rand->Poisson(4./(nc+4));
	  (*UnDigitized)(nr,nc)=(*Digitized)(nr,nc)*factor;
	  
	};
	
    };
  return Digitized;
};
void Gen::CallFill(Int_t ncalls)
{
  for(Int_t j=0;j<ncalls;j++)
    {
      
      
      Float_t Ereturn=Fill();
      if(GenCnt>100000)
	{
	  printf("Count j=%d\n",j);
	  GenCnt=1;
	};
      if(Ereturn>=CallFillLimitE)break;
 
    };
  
};
void Gen::EnergyFromPhoton(TVector3* p_v3,Double_t photE)
{
  for(Int_t col=0;col<Ncols;col++)
    {
      xx = (( (Double_t) col + 0.5 )  - p_v3->X());

      if(fabs(xx)>5)continue;
      xx=xx*TowerWidth ;
      for(Int_t row=0;row<Nrows;row++)
	{

	  yy = (( (Double_t) row + 0.5 )  - p_v3->Y());
	  if(fabs(yy)>5)continue;
	  yy=yy*TowerWidth ;
	  eSS=photE*(FitTower::we.fcnSS)->Eval(xx,yy);
	  (*NextEvt)(row,col)+=eSS;	  
	};
    };
  
};
void Gen::WriteAll()
{
  for(Int_t ir=0;ir<Nrows;ir++)
    {
      for(Int_t ic=0;ic<Ncols;ic++)
	{
	  Exy[ir][ic]->Write();
	  ZvsADC[ir][ic]->Write();
	  PhivsADC[ir][ic]->Write();
	  Z_Phi_40[ir][ic]->Write();
	  LocalPos[ir][ic]->Write();
	  LocalPos125[ir][ic]->Write();
	  Mcell[ir][ic]->Write(); 
	  NphotvsE[ir][ic]->Write();
	  EphXY[ir][ic]->Write();
	};
      
    };
  g_out->Write();
  phi1->Write();
  phi2->Write();
  phi3->Write();
  Eratch1F->Write();
  EratchP->Write();
  Eratch1->Write();
  Eratch2->Write();
  Episum->Write();
  Epipsum->Write();
  EpipRsum->Write();
  EpipbRsum->Write();
  Epigsum->Write();
  GenEvsY->Write();
  GenpEvsY->Write();
  hEgen->Write();
  hEpass->Write();
  GenZbeam->Write();
  GenPZbeam->Write();
  GenEvsPt->Write();
  GenPEvsPt->Write();
  GenPPEvsPt->Write();
  GenThvsPh->Write();
  GenPThvsPh->Write(); 
  GenPPThvsPh->Write(); 
  GenYvsPtWPt->Write();
  GenPYvsPtWPt->Write();
  GenPPYvsPtWPt->Write();
  Mass->Write();
  Nphotons->Write();
};
