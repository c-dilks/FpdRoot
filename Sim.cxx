#include "Sim.h"

ClassImp(Sim)

Sim::Sim(char* files,Geom* pgeom):simh112(files)
{
  SelectTest=0;
  SimType="fmssim";
  isSim=3;
  pGeom=pgeom;
  for(Int_t j=0;j<6;j++)
    {
      if(j<2)
	{
	  ene[j][0]=new TMatrix(1,1);
	  ene[j][1]=new TMatrix(34,17);
	}
      else if(j<4)
	{
	  ene[j][0]=new TMatrix(1,1);	
	  ene[j][1]=new TMatrix(24,12);	
	}
      else
	{
	  ene[j][0]=new TMatrix(1,1);	
	  ene[j][1]=new TMatrix(1,1);	
	}
      (*ene[j][0])=0;
      (*ene[j][1])=0;
      noise[j][0]=new TMatrix(*ene[j][0]);
      noise[j][1]=new TMatrix(*ene[j][1]);
      adc[j][0]=new TMatrix(*ene[j][0]);
      adc[j][1]=new TMatrix(*ene[j][1]);
      sumene[j][0]=new TMatrix(*ene[j][0]);
      sumene[j][1]=new TMatrix(*ene[j][1]);
      newene[j][0]=new TMatrix(*ene[j][0]);
      newene[j][1]=new TMatrix(*ene[j][1]);

      sigped[1][j]=Esigped[j];
      sigped[1][j]=Wsigped[j];
      //      printf("sigped[0][%d]=%f sigped[1][%d]=%f \n",j,sigped[0][j],j,sigped[1][j]);
    }
  GAIN=1.045;
  p_ran=new TRandom1(0);
  for(Int_t k=0;k<4;k++)HighFour[k]=-1;
  shape=new TH3F("shape","shape",40,-2,2,40,-2,2,50,0,1.);
  skipID=-999;
  
}

Sim::~Sim()
{
};

Float_t Sim::digitBits(Int_t nstb,Int_t ew)
{
  return .05;
  return .2;
  if(pGeom->datatype[6]=="digitBits")
    {
     
    }
  return 0;
};

Bool_t Sim::FillNoise()
{
  for(Int_t iew=0;iew<2;iew++)
    {
      for(Int_t instb=0;instb<6;instb++)
	{
	  TMatrix* pn=noise[instb][iew];
	  for(Int_t ir=0;ir<pn->GetNrows();ir++)
	    {
	      for(Int_t ic=0;ic<pn->GetNcols();ic++)
		{
		  (*pn)(ir,ic)=p_ran->Gaus(0,sigped[iew][instb])*
		    digitBits(instb,iew);
		};
	    };
	};
    };
};

Bool_t Sim::AddNoise()
{
  FillNoise();
  for(Int_t i=0;i<6;i++)
    {
      (*ene[i][0])= (*(ene[i][0]))+(*(noise[i][0]));
      (*ene[i][1])= (*(ene[i][1]))+(*(noise[i][1]));
    };
  return true;
};

Bool_t Sim::GenADC(Float_t Mmin,Float_t Mmax)
{
  if(!FillHits(Mmin,Mmax))return false;
  if(!AddNoise())return false;
  for(Int_t i=0;i<6;i++)
    {
      for(Int_t j=0;j<2;j++)
	{
	  TMatrix* pm=ene[i][j];
	  TMatrix* pa=adc[i][j];
	  TMatrix* pb=newene[i][j];
	  TMatrix* pc=sumene[i][j];
	  Int_t nr=pm->GetNrows();
	  Int_t nc=pm->GetNcols();
	  for(Int_t ir=0;ir<nr;ir++)
	    {
	      for(Int_t ic=0;ic<nc;ic++)
		{
		  Double_t d0=1./digitBits(i,j);
		  d0= (*pm)(ir,ic)*d0;
		  (*pa)(ir,ic)=(Int_t) d0;
		  (*pb)(ir,ic)=((Int_t) d0)*digitBits(i,j);
		  (*pc)(ir,ic)+=(*pb)(ir,ic);
		};
	    };
	};
    };
  return true;
};

Bool_t Sim::Fillshape()
{
   for(Int_t i=0;i<Nhit;i++)
    {
      if(HitTest(i) && getvolid(id[i]) )
	{
	  Int_t tr=trk[i]-1;
	  
	  if(ew==1 && nstb==1&& tr==HighFour[0] && m[tr]<.001)
	    {
	      TVector3 v3=ProjectTrack(tr,*pGeom->ZFPD(1,1));	  
	      Float_t x0=-(v3.X()-*pGeom->xOffset(1,1))
		/(*pGeom->FpdTowWid(1,1));
	      Float_t y0=-(v3.Y()-*pGeom->yOffset(1,1))
		/(*pGeom->FpdTowWid(1,1));
	      Float_t xch=(fmod(ch-1,7))+.5;
	      Float_t ych=((Int_t) (ch-1)/7)+.5;
	      Float_t en=e[tr];
	      if(en>20 && en<40 && fabs(y0-3.5)<2.5 && fabs(x0-3.5)<2.5){
		shape->Fill(xch-x0,ych-y0,hit[i]/en);
		//	printf("xch-x0=%f ych-y0=%f \n",xch-x0,ych-y0);
	      };
	    };
	};
    };  
};

Bool_t Sim::FillHits(Float_t Mmin,Float_t Mmax)
{
  Bool_t retval=true;
  for(Int_t i=0;i<6;i++)
    {
      *(ene[i][0])=0.;
      *(ene[i][1])=0.;
    };
  //  printf("Nhit=%d \n",Nhit);
  for(Int_t i=0;i<Nhit;i++)
    {
      //printf("test here %d %d %d\n", id[i],HitTest(i),getvolid(id[i]));

      if(HitTest(i) && getvolid(id[i]))
	{
	  //	  printf("test here %d %d %d\n", id[i],trk[i],mo[trk[i]-1][0]-1);
	  // printf("ew=%d nstb=%d m=%f \n",ew,nstb,m[trk[i]-1]);

	  Float_t TrackMass=-200.;
	  Int_t mother=-1;
	  if(trk[i]>0 && trk[i]<=500){
	    TrackMass=m[trk[i]-1];
	    mother=mo[trk[i]-1][0]-1;
	  };
	  Bool_t skip=false;
	  if(mother>=0)
	    {
	      if(pdg[mother]==skipID)skip=true;
	    };
	  if(skip) printf("Skip hits for mother id=%d trk=%d daughter id=%d trk=%d ch=%d \n",pdg[mother],mother,pdg[trk[i]-1],trk[i]-1,ch);	    
	  if(towsmd==1 && TrackMass>Mmin && TrackMass<Mmax && !skip)
	    {
	      TMatrix* mm=ene[nstb-1][ew-1];
	      Int_t chan=ch-1;
	      Int_t nrows=mm->GetNrows();
	      Int_t ncols=mm->GetNcols();
	      Int_t row=chan/ncols;
	      //Int_t col=chan-row*nrows;
	      Int_t col=chan%ncols;
	      //	      printf("volid=%d ew=%d nstb=%d chan=%d nrows=%d ncols=%d row=%d col=%d\n",id[i],ew,nstb,chan,nrows,ncols,row,col);
	      if(row>=0 && col>=0 && row<nrows && col<ncols)
		{
		  Float_t hval=(*mm)(row,col);
		  hval+=hit[i]*GAIN;
		  (*mm)(row,col)=hval;
		  //		  printf("Hitno=%d Ew=%d NSTB=%d ch=%d row=%d col=%d hit=%f hval=%f \n",i,ew,nstb,ch,row,col,hit[i],hval);
		}
	      else
		{
		  //printf("sim hit error: ew=%d nstb=%d ch=%d row=%d col=%d nrows=%d ncols=%d \n",ew,nstb,ch,row,col,nrows,ncols);
		  retval=false;
		};
	    };
	};
    };

  return retval;
}

Bool_t Sim::PrintHits(Float_t hitmin)
{
  for(Int_t j=0;j<Nhit;j++)
    {
      if(getvolid(id[j])){
	Int_t trknum=trk[j];
	Int_t ptrknum=ptrk[trknum];
	if(hit[j]>hitmin){
	  
	  printf("j=%d vid=%d towsmd=%d  ew=%d nstb=%d ch=%d hit=%f \n",j,vid,towsmd,ew,nstb,ch,hit[j]);
	  printf("     trknum-1=%d m=%f e=%f p0=%f p1=%f p2=%f  v1=%f v2=%f v3=%f v4=%f \n",trknum-1, m[trknum-1],e[trknum-1],p[trknum-1][0],p[trknum-1][1],p[trknum-1][2],v[trknum-1][0],v[trknum-1][1],v[trknum-1][2],v[trknum-1][3]);
	  printf("     ptrk-1=%d da0-1=%d da1-1=%d \n",ptrk[trknum-1]-1, da[trknum-1][0]-1,da[trknum-1][1]-1);
	};
      };
    };
};
TVector3 Sim::ProjectTrack(Int_t tracknum,Float_t zdet)
{
  Int_t tn=tracknum;
  TVector3 orig(v[tn][0],v[tn][1],v[tn][2]);
  TVector3 dis(p[tn][0],p[tn][1],p[tn][2]);
  dis.SetMag(1.);
  TVector3 ret(0.,0.,0.);
  if(dis.Z()*zdet<=0)return ret;
  dis.SetMag((zdet-orig.Z())/dis.Z());
  TVector3 vv=orig+dis;
  /*  printf("dis=(%f, %f, %f) orig+dis= (%f, %f, %f) \n",
	 dis.X(),dis.Y(),dis.Z(),vv.X(),vv.Y(),vv.Z());
  */
  return vv;
};

Bool_t Sim::getvolid(Int_t vid_in)
{
  vid=vid_in;
  nstb=0;
  towsmd=0;
  ew=0;
  vh=0;
  if(isSim!=3)return false;
  if(vid<100000)
    {
      towsmd=1;
      ew=(int) fmod(vid/10000,10);
      nstb=(int) fmod(vid/1000,10);
      vh=0;
      ch=(int) fmod(vid,1000);
    }
  else if(vid<1000000)
    {     
      towsmd=2;
      ew=(int) fmod(vid/10000,10);
      nstb=(int) fmod(vid/1000,10);
      vh=(int) fmod(vid/100,10);
      ch=(int) fmod(vid,100);
    }
  else 
    {
      towsmd=3;
      ew=(int) fmod((vid-1000000)/1000,10);
      nstb=0;      
      vh=(int) fmod((vid-1000000)/100,10);
      ch=(int) fmod(vid,100);
    };
  return true;
};

Bool_t Sim::FillHighFour(Float_t tstmass,Float_t dmass)
{
  for(Int_t i=0;i<4;i++)HighFour[i]=-1;
  for(Int_t i=0;i<Nhit;i++)
    {
      CheckForHigh(trk[i]-1);
    };
};

Bool_t Sim::CheckForHigh(Int_t id_chk,Float_t tstmass,Float_t dmass)
{
  // id_chk starts at zero
  if(m[id_chk]>tstmass-dmass && m[id_chk]<tstmass+dmass)
    {
      for(Int_t k=0;k<4;k++)
	{
	  if(HighFour[k]<0){
	    HighFour[k]=id_chk;
	    return true;
	  };
	  if(id_chk==HighFour[k])return true; // Ignore repeat track
	  if(e[id_chk]>e[HighFour[k]])
	    {
	      for(Int_t j=0;j<3-k;j++)
		{
		  HighFour[3-j]=HighFour[2-j];
		};
	      HighFour[k]=id_chk;
	      return true;
	    };
	};
      return false;
    }
  else return false;
}

Bool_t Sim::HitTest(Int_t hitid)
{
  if(SelectTest==0)return true;
  /*
  return true;
  
  Int_t parent = mo[trk[hitid]-1][0]-1;
  Int_t Hipar = mo[HighFour[0]][0]-1;
  if(parent>=0 && parent< Ntrack)
    {
      if(parent==Hipar && e[Hipar]>10)return true;
    };
  return false;
  if(trk[hitid]-1==HighFour[0]&&e[HighFour[0]]>10)return true;// check for highest en photon
  if(hitid>0 && hitid<=Nhit)
    {
      Float_t Pi=3.1415926;
      Int_t trkN=trk[hitid]-1;
      TLorentzVector vec(p[trkN][0],p[trkN][1],p[trkN][2],e[trkN]);
      if(fabs(vec.Eta()-(-3.8))<.05 && fabs(vec.Phi()-Pi)<.05)return true; 
    };
  */
  //  printf("hitid=%d\n",hitid);
  if(hitid>0 && hitid<=Nhit)
    {
      // only if from pi0...

      //      printf("test  ****  %d", id[hitid]);
      Int_t itr=trk[hitid]-1;
      if(itr<0 || itr>Ntrack-1)return false;
      Int_t io=mo[itr][0]-1;
      
      if(io>=0 && io< Ntrack){
      if(!getvolid(id[hitid]))return false;

      if(e[io]>0. && ew==2 &&nstb>0 && nstb<5 && fabs(m[io]-.135)<.01)
	  {
	    return true;
	    Int_t ioa=da[ mo[HighFour[0]][0]-1][0]-1;
	    Int_t iob=da[ mo[HighFour[0]][0]-1][1]-1;
	    TVector3 v3=ProjectTrack(io,*pGeom->ZFPD(ew,nstb));
	    Float_t x0=-(v3.X()-*pGeom->xOffset(ew,nstb))
	      /(*pGeom->FpdTowWid(ew,nstb));
	    Float_t y0=-(v3.Y()-*pGeom->yOffset(ew,nstb))
	      /(*pGeom->FpdTowWid(ew,nstb));
	    //	    printf("ew=%d nstb %d m=%f x0=%f y0=%f \n",ew,nstb,m[io],x0,y0);

	    if(x0 > 0 && x0<14 )return true;
	  };
      }
      return false;
    };
};
