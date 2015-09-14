#include "poutTree.h"
ClassImp(poutTree)
poutTree::poutTree(char* filelist)
{
  for(int k=0;k<4;k++)
    {
      gr[k]=0;
      gradc[k]=0;
    };
  ClusterwithYPhi(true);
  spin=0;
  EnableEdepCorr=0;// initially turn off Edep corrections
  nphotons=0;
  nwrds=0;
  strig=0;
  for(int ib=0;ib<5;ib++){BBcSums[ib]=0;};
  BBcVertex=0;
  EventN=-1;
  Rnum=-1;  
  Bunchid7bit=0;
  vlist=new TObjArray(80,0);
  dlist=new TObjArray(80,0);
  softlist=new TObjArray(80,0);
  ClusterList=new TObjArray(80,0);
  Float_t eDepInt[8]= {1.0,1.0,1.0,1.0,  1.18,1.183,1.1346,1.1549};
  //Float_t eDepInt[8]= {1.0,1.0,1.0,1.0,  1.119, 1.166,  1.159,1.159};
  //Float_t eDepSlope[8]={0,0,0,0,-.00338,-.00542,-.00291,-.00279};
   Float_t eDepSlope[8]={0,0,0,0,  -.00886,-.00918,-.0079,-.0078};
  //Float_t eDepSlope[8]={0,0,0,0,-.012,-.015,-.010,-.010};

  for(Int_t ii=0;ii<8;ii++){EDepInt[ii]=eDepInt[ii];EDepSlope[ii]=eDepSlope[ii];};
  for(Int_t ii=0;ii<8;ii++){EDepInt[ii]=1.;EDepSlope[ii]=0;};
  p_fms[2][0]=new TH2F("SNfms0","SNfms0",12,0,12,24,0,24);
  p_fms[3][0]=new TH2F("SSfms0","SSfms0",12,0,12,24,0,24);
  p_fms[0][0]=new TH2F("LNfms0","LNfms0",17,0,17,34,0,34);
  p_fms[1][0]=new TH2F("LSfms0","LSfms0",17,0,17,34,0,34);
  p_fms[2][1]=new TH2F("SNfms1","SNfms1",12,0,12,24,0,24);
  p_fms[3][1]=new TH2F("SSfms1","SSfms1",12,0,12,24,0,24);
  p_fms[0][1]=new TH2F("LNfms1","LNfms1",17,0,17,34,0,34);
  p_fms[1][1]=new TH2F("LSfms1","LSfms1",17,0,17,34,0,34);
  scratchlist=new TObjArray(80,0);
  for(Int_t ix=1;ix<3;ix++){MinEnergy[ix]=0.;};
  p_out= new TChain("p_out");
  FILE* pf; 
  Int_t FileCnt=0;
  if(pf=fopen(filelist,"r"))
    {
      char fname[200];
      while(!feof(pf))
	{
	  if(fscanf(pf,"%s\n",fname)>0)
	    {
	      FileCnt++;
	      p_out->Add(fname);
	    };
	  if(FileCnt%20==0)printf("%d files added to TChain\n",FileCnt);
	};
      nentries=p_out->GetEntries();
      printf("Number files read=%d, Number of Events=%d \n",FileCnt,nentries);
      fclose(pf);
    };

  p_out->SetBranchAddress("spin",&spin);
  p_out->SetBranchAddress("nphotons",&nphotons);
  p_out->SetBranchAddress("br_nwrds",&nwrds);
  p_out->SetBranchAddress("br_types",tpes);
  p_out->SetBranchAddress("br_pxyzt",pxyzt);
  p_out->SetBranchAddress("br_Rnum",&(Rnum));
  p_out->SetBranchAddress("br_Bunchid7bit",&(Bunchid7bit));
  p_out->SetBranchAddress("br_BBcSums",BBcSums);
  p_out->SetBranchAddress("br_BBcVertex",&BBcVertex);
  p_out->SetBranchAddress("br_EventN",&EventN);
  p_out->SetBranchAddress("br_ievt",&ievt);
  p_out->SetBranchAddress("br_L2sum",L2sum);
  p_out->SetBranchAddress("br_lastdsm",lastdsm);
  

  TBranch        *b_BQTNE;   //!
  TBranch        *b_BQTNW;   //!
  TBranch        *b_br_QTEBBCInd;   //!
  TBranch        *b_br_QTWBBCInd;   //!
  TBranch        *b_br_QTEBBCTAC;   //!
  TBranch        *b_br_QTWBBCTAC;   //!
  TBranch        *b_br_QTEBBCADC;   //!
  TBranch        *b_br_QTWBBCADC;   //!
  TBranch        *b_QTBVertex;   //!

  TBranch        *b_RPE_QTN;   //!
  TBranch        *b_RPW_QTN;   //!
  TBranch        *b_br_RPE_Idx;   //!
  TBranch        *b_br_RPW_Idx;   //!
  TBranch        *b_br_RPE_TAC;   //!
  TBranch        *b_br_RPW_TAC;   //!
  TBranch        *b_br_RPE_ADC;   //!
  TBranch        *b_br_RPW_ADC;   //!
  TBranch        *b_RPvertex;   //!
  

  p_out->SetBranchAddress("br_QTNE",&(qtbbc.QTNE),&b_BQTNE);
  p_out->SetBranchAddress("br_QTNW",&(qtbbc.QTNW),&b_BQTNW);
  if(p_out->GetBranch("br_QTEBBCInd"))p_out->SetBranchAddress("br_QTEBBCInd",qtbbc.QTEBBCInd);
  if(p_out->GetBranch("br_QTWBBCInd"))p_out->SetBranchAddress("br_QTWBBCInd",qtbbc.QTWBBCInd);
  if(p_out->GetBranch("br_QTEBBCADC"))p_out->SetBranchAddress("br_QTEBBCADC",qtbbc.QTEBBCADC);
  if(p_out->GetBranch("br_QTWBBCADC"))p_out->SetBranchAddress("br_QTWBBCADC",qtbbc.QTWBBCADC);
  if(p_out->GetBranch("br_QTEBBCTAC"))p_out->SetBranchAddress("br_QTEBBCTAC",qtbbc.QTEBBCTAC);
  if(p_out->GetBranch("br_QTWBBCTAC"))p_out->SetBranchAddress("br_QTWBBCTAC",qtbbc.QTWBBCTAC);
  if(p_out->GetBranch("br_QTBVertex"))p_out->SetBranchAddress("br_QTBVertex",&(qtbbc.vertex));


  p_out->SetBranchAddress("br_RPE_QTN",&(qtrp.NE),&b_RPE_QTN);
  p_out->SetBranchAddress("br_RPW_QTN",&(qtrp.NW),&b_RPW_QTN);
  if(p_out->GetBranch("br_RPE_Idx"))p_out->SetBranchAddress("br_RPE_Idx",qtrp.RPE_Idx);
  if(p_out->GetBranch("br_RPW_Idx"))p_out->SetBranchAddress("br_RPW_Idx",qtrp.RPW_Idx);
  if(p_out->GetBranch("br_RPE_ADC"))p_out->SetBranchAddress("br_RPE_ADC",qtrp.RPE_ADC);
  if(p_out->GetBranch("br_RPW_ADC"))p_out->SetBranchAddress("br_RPW_ADC",qtrp.RPW_ADC);
  if(p_out->GetBranch("br_RPE_TAC"))p_out->SetBranchAddress("br_RPE_TAC",qtrp.RPE_TAC);
  if(p_out->GetBranch("br_RPW_TAC"))p_out->SetBranchAddress("br_RPW_TAC",qtrp.RPW_TAC);
  if(p_out->GetBranch("br_RPvertex"))p_out->SetBranchAddress("br_RPvertex",&(qtrp.vertex));
  

  if(p_out->GetBranch("br_adc"))p_out->SetBranchAddress("br_adc",&adc);
  if(p_out->GetBranch("br_nSavedHits"))
    {
      p_out->SetBranchAddress("br_nSavedHits",&nSavedHits);
      p_out->SetBranchAddress("br_SavedHits",&SavedHits);
    };
  if(p_out->GetBranch("br_nCluster"))
    {
      p_out->SetBranchAddress("br_nCluster",&(nCluster));
      p_out->SetBranchAddress("br_SCIndex",(SavedCluHitIndex));
    };

  if(p_out->GetBranch("br_nPhotonClu"))
    {
      p_out->SetBranchAddress("br_nPhotonClu",&(nPhotonClu));
      p_out->SetBranchAddress("br_SPCIndex",(SavedPhotonCluIndex));
      p_out->SetBranchAddress("br_SPCEnergy",(SavedPhotonCluEnergy));
    };


  if(p_out->GetBranch("br_TrigBits"))p_out->SetBranchAddress("br_TrigBits",&TrigBits);

  

};
poutTree::~poutTree()
{
  delete p_out;
  if(vlist!=0){
    vlist->Delete();
  };
  for(int i=0;i<4;i++)
    {
      if(gr[i])delete gr[i];
      if(gradc[i])delete gradc[i];
    };
  if(dlist!=0){
    dlist->Delete();
  };
  if(softlist!=0){
    softlist->Delete();
  };
  delete softlist;
  delete vlist;
  delete scratchlist;
  delete dlist;
  if(ClusterList!=0)ClusterList->Delete();
  delete ClusterList;
  for(int jj=0;jj<4;jj++)
    {
      if(p_fms[jj][0])delete p_fms[jj][0];
      if(p_fms[jj][1])delete p_fms[jj][1];
    }
};
void poutTree::histpair(LVec* vv,Geom* p_geom,Float_t zinteraction,Int_t hnum)
{
  if(hnum>=0 && hnum<2)
    {
      TVector3 v3=PosInDet(vv,p_geom,zinteraction);
      
      TVector3 v13=p_geom->LocalXYZ(2,vv->Nstb,v3,true);
      p_fms[vv->Nstb-1][hnum]->Fill(v13.X(),v13.Y());
      LVec* vq=vv->Partner;
      TVector3 v4=PosInDet(vq,p_geom,zinteraction);
      TVector3 v14=p_geom->LocalXYZ(2,vq->Nstb,v4,true);
      p_fms[vq->Nstb-1][hnum]->Fill(v14.X(),v14.Y());
    };
};
Int_t poutTree::NearEdge(LVec* lv,Geom* pgeom,Float_t d,Float_t zinteraction)
{
  //return 1 if within d of inner edge  (interface if small)
  //return 2 if within d of outter edge (interface if large)
  //
  Int_t iedge[2][2][2][2]={{{{0,0},{0,0}},{{0,0},{9,7}}}, 
			   {{{17,12},{34,24}},{{8,5},{25,17}}}};
  //[min=0;max=1][out=0/in=1][h=0,v=1][l=0/s=1]
  
  TVector3 v3=PosInDet(lv,pgeom,zinteraction);
  Int_t min,out,h;
  Int_t max,in,v;
  min=out=h=0;
  max=in=v=1;
  TVector3 v13=pgeom->LocalXYZ(2,lv->Nstb,v3,true);
  Int_t result=0;
  if(lv->Iew==2 )
    {
      Int_t ns=lv->Nstb;
      Int_t ls=(ns-1)/2;
      
      //Near inner hole
      if( 
	 v13.X()<(iedge[max][in][h][ls]+d) && 
	 v13.Y()>(iedge[min][in][v][ls]-d) &&
	 v13.Y()<(iedge[max][in][v][ls]+d)
	 )result=result|1;
      
      //near horizontal boundary
      if(v13.X()<(iedge[min][out][h][ls]+d))result=result|2;
      
      //outside
      
      if(v13.X()>(iedge[max][out][h][ls]-d))result=result|4;
      if(v13.Y()<(iedge[min][out][v][ls]+d))result=result|4;
      if(v13.Y()>(iedge[max][out][v][ls]-d))result=result|4;
      
    };
  return result;
};
Int_t poutTree::GetEntry(Int_t evt)
{
  /*
  int j0,j1,j2,j3,j4,j5;
  j0=(int) (qtbbc.QTNE);
  std::cout<<"evt="<<evt<<" QTNE="<<j0<<" \n";
  for(int j1=0;j1<j0;j1++) 
    {
      j2=(int) (qtbbc.QTEBBCInd[j1]);
    std::cout<<"j="<<j1<<"_"<<evt<<" QTEBBCInd="<<j2<<" ";
      j3=(int) (qtbbc.QTEBBCTAC[j1]);
    std::cout<<"j="<<j1<<"_"<<evt<<" QTEBBCTAC="<<j3<<" ";
    };
  std::cout<<" \n";
  */
  ClearScratch();
  vlist->Delete();
  dlist->Delete();
  softlist->Delete();
  Int_t nbytes=p_out->GetEntry(evt); 
  entry=evt;


  YellowUp=YellowDn=BlueUp=BlueDn=false;
  
  if(spin==1 || spin==3)YellowUp=true;
  if(spin==0 || spin==2)YellowDn=true;
  if(spin==2 || spin==3)BlueUp=true;
  if(spin==0 || spin==1)BlueDn=true;
  for(Int_t i1=0;i1<10;i1++)
    {
      nlv[i1]=0;
      nph2[i1]=0;
      nph3[i1]=0;
    };
  Nphotons=0;
  Ndet=0;
  TotaldetE=0.;
  TotalphotE=0.;
  Int_t photcnt=0;
  Int_t vtyp;
  Int_t nwords=nwrds;
  LVec* tmpLV=0;

  for(Int_t k=0;k<nwords;k+=4){
    if(tpes[k]<9)
      {
	Int_t kk=k/4;
	vtyp=tpes[k];
	plv[vtyp][nlv[vtyp]]=kk;
	vec[kk].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
	//vec[k]=EdepCorr(vec[k],vtyp); // Energy dependent Correction
	if(nlv[vtyp]<10)
	  {
	    nlv[vtyp]++;
	    dlist->AddAt(tmpLV=new LVec(vtyp,vec[kk]),Ndet);	    	
	    p_detvec[Ndet]=&(vec[kk]);
	    Ndet++;
	    TotaldetE+=vec[kk].E();
	  };
      }
    else if((tpes[k]>300) && (tpes[k]<309))
      {
	Int_t kk=k/4;
	vtyp=tpes[k]-300;
	if(nph2[vtyp]<10 && nph2[vtyp]<10)
	  {
	    pph2[vtyp][nph2[vtyp]]=kk;
	    pph3[vtyp][nph3[vtyp]]=kk;
	    vec[kk].SetXYZT(pxyzt[k],pxyzt[k+1],pxyzt[k+2],pxyzt[k+3]);
	    //	vec[kk]=EdepCorr(vec[k],vtyp+300);// Energy dependent Correction
	    if(vec[kk].E()>MinEnergy[vtyp])
	      {
		nph2[vtyp]++;
		p_photvec[Nphotons]=&(vec[kk]);
		vlist->AddAt(tmpLV=new LVec(vtyp,vec[kk]),Nphotons);
		tmpLV->PhotOrder=photcnt++;
		Nphotons++;
		TotalphotE+=vec[kk].E();
	      }
	    else
	      {
		Int_t kk=k/4;
		softlist->Add(tmpLV=new LVec(vtyp,vec[kk]));
		tmpLV->PhotOrder=photcnt++;
	      };
	  };
      };
  };
  return nbytes;
};
TVector3 poutTree::PosInDet(LVec* lv,Geom* pgeom,Float_t zinteraction)
{
  TVector3 unit=lv->Vect();
  unit.SetMag(1.);
  Float_t znew=*pgeom->ZFPD(lv->Iew,lv->Nstb)-zinteraction;
  unit=(znew/unit.Z())*unit;
  unit.SetZ(unit.Z()+zinteraction);
  return unit;
};
Int_t poutTree::GetNPhotVec(Int_t iew,Int_t nstb)
{
  if(iew<1 || iew>2)return 0;
  if(nstb<1|| nstb>4)return 0;
  Int_t k=(iew-1)*4+nstb;
  return nph2[k];
  
};
Int_t poutTree::GetNDetVec(Int_t iew,Int_t nstb)
{
  if(iew<1 || iew>2)return 0;
  if(nstb<1 || nstb>4)return 0;
  Int_t k=(iew-1)*4+nstb;
  return nlv[k];
};

Int_t poutTree:: WhatEW(LVec* lv)
{
  if(!lv)return -1;
  return lv->Iew;
};

Int_t poutTree:: WhatNSTB(LVec* lv)
{
  if(!lv)return -1;
  return lv->Nstb;
};
  
Int_t poutTree::WhatRow0(LVec* lv,Geom* pgeom,Float_t zinteraction)
{
  if(!pgeom->FMSGeom)return -1;
  if(!lv)return -1;
  TVector3 v3=PosInDet(lv,pgeom,zinteraction);
  TVector3 v_L=pgeom->LocalXYZ(2,lv->Nstb,v3,true);
  Int_t nrows=24;
  if(lv->Nstb>4)return -1;
  if(lv->Nstb<3)nrows=34;
  if(v_L.Y()<0 || v_L.Y()>=nrows)return -1;
  return (int) floor(v_L.Y()); 
};
Int_t poutTree::WhatCol0(LVec* lv,Geom* pgeom,Float_t zinteraction)
{
  if(!pgeom->FMSGeom)return -1;
  if(!lv)return -1;
  TVector3 v3=PosInDet(lv,pgeom,zinteraction);
  TVector3 v_L=pgeom->LocalXYZ(2,lv->Nstb,v3,true);
  Int_t ncols=12;
  if(lv->Nstb>4)return -1;
  if(lv->Nstb<3)ncols=17;
  if(v_L.X()<0 || v_L.X()>=ncols)return -1;
  return (int) floor(v_L.X()); 
};

TLorentzVector poutTree::GetPhotVec(Int_t iew,Int_t nstb,Int_t index)
{
  if(iew<1 || iew>2)return nullvec;
  if(nstb<1 || nstb>4)return nullvec;
  Int_t k=(iew-1)*4+nstb;
  if(index<0 || ( index >= nph2[k]))return nullvec;
  /*
  if( false)
  {
    int jcnt=0;
    while(jcnt<80){printf("vec[jcnt]=%f\n",vec[jcnt].E());jcnt++;};
    printf("k=%d index=%d nph2[k]=%d \n",k,index,nph2[k]);
    printf("pph2=%d",pph2[k][index]);
    printf("vec=%f\n",vec[pph2[k][index]].E());
  }
  */
  return vec[  pph2[k][index ] ];  
};
void poutTree::Print()
{
  printf("entry Number = %d EventN=%d (of %d) Run Number=%d  Nphotons=%d\n",
	 entry,EventN,nentries,Rnum,Nphotons);
  for(Int_t iew=1;iew<3;iew++)
    {
      for(Int_t nstb=1;nstb<5;nstb++)
	{
	  if(GetNDetVec(iew,nstb)>0)
	    {
	      printf("iew=%d nstb=%d Energy=%f Nphot=%d \n",iew,nstb,
		     GetDetVec(iew,nstb,0).E(),
		     GetNPhotVec(iew,nstb));
	    };
	};
    }
};
void poutTree::PrintList(TObjArray* plist)
{
  TIter next(plist);
  Int_t cnt=0;
    while(TObject* tmpob= next())
    {
      if(tmpob->ClassName()=="LVec" || true)
	{
	  LVec* v=(LVec*) tmpob;
	  printf("%d) Iew=%d Nstb=%d \n",cnt,v->Iew,v->Nstb);
	  printf("  %d) Px=%f Py=%f Pz=%f E=%f Y= %f phi=%f\n",cnt++,v->Px(),v->Py(),v->Pz(),v->E(),v->PseudoRapidity(),v->Phi()); 
	  
	};
    };
};
TLorentzVector poutTree::GetDetVec(Int_t iew,Int_t nstb,Int_t index)
{

  if(iew<1 || iew>2)return nullvec;
  if(nstb<1 || nstb>4)return nullvec;
  Int_t k=(iew-1)*4+nstb;
  if(index<0 || (index>= nlv[k]))return nullvec;
  return vec[  plv[k][index ] ];
};

 
LVec* poutTree::GetPairNearMass(TObjArray* plist,Float_t massgoal)
{

  Int_t len=plist->GetEntries();
  LVec* pbst1;
  if(len<2)return 0;
  LVec* pbst2;
  Float_t diff=1000000.;
  Float_t d;
  LVec* pv2;
  LVec* pv1;
  for(Int_t i1=0;i1<len;i1++)
    {
      pv1=(LVec*) plist->At(i1);
      for(Int_t j1=0;j1<i1;j1++)
	{
	  pv2=(LVec*) plist->At(j1);
	  if((d=fabs(pv1->PairMass(pv2)-massgoal))<diff)
	    {
	      pbst1=pv1;
	      pbst2=pv2;
	      diff=d;
	    };
	};
    };

  pbst1->SetPartner(pbst2);
  pbst2->Partner=pbst1;

  return pbst1;
};
TObjArray* poutTree::AddToScratch(Int_t vtyp)
{
  Int_t n=vlist->GetEntries();
  TIter next(vlist);
  while(LVec* v=(LVec*) next())
    {
      if(v->Vtype==vtyp)scratchlist->Add(v);
    };

  return scratchlist;
};

TObjArray* poutTree::RemoveFromScratch(LVec* v)
{
  scratchlist->Remove(v);
  return scratchlist;
};

TObjArray* poutTree::AllToScratch(Bool_t includesoft)
{
  ClearScratch();
  TIter next(vlist);
  while(LVec* v=(LVec*) next())
    {
      scratchlist->Add(v);
    };

  if(includesoft)
    {
      TIter next1(softlist);
      while(LVec* v=(LVec*) next1())
	{
	  scratchlist->Add(v);
	};
    };
  
  return scratchlist;
};

TLorentzVector poutTree::SumScratch()
{
  TLorentzVector sum(0,0,0,0);
  TIter next(scratchlist);
  while(LVec* v=(LVec*) next())
    {
      sum+=*((TLorentzVector*) v);
    };
  return sum;
}; 
TLorentzVector poutTree::EdepCorr(TLorentzVector vold, Int_t vtype)
{
  if(EnableEdepCorr==0) return vold;// disable EdepCorr
  return vold;//never used
  Float_t oldslope=0;
  Float_t oldint=1.;
  Float_t nlgain=1.;
  if(vtype>1 && vtype<9)
    {
      oldslope=EDepSlope[vtype-1];
      oldint=EDepInt[vtype-1];
      nlgain=oldint+oldslope*vold.E();
      return vold*nlgain;
    }
  else if(vtype>301 && vtype<309)
    {
      oldslope=EDepSlope[vtype-301];
      oldint=EDepInt[vtype-301];
      oldslope=-.003;
      oldint=1.08;
      nlgain=oldint+oldslope*vold.E();
      Float_t e0=vold.E();
      if(vtype>306)
	{
	  nlgain=1.303-.0277*pow(e0,1)+.00094*pow(e0,2)-.0000247*pow(e0,3)+.000000333*pow(e0,4)-.00000000163*pow(e0,5);
	}
      else
	{
	  if(vtype==305)nlgain=.792+.018*pow(e0,1)-.000427*pow(e0,2)+.000003053*pow(e0,3);
	  if(vtype==306)nlgain=1.34-.0401*pow(e0,1)+.000814*pow(e0,2);
	};

      TLorentzVector vret=vold*nlgain;
      return vret;
    }
  else
    { 
      // no change
      return vold;
    };
};
TVector3 poutTree::reframe(TLorentzVector vec,TLorentzVector zvec,TLorentzVector bvec)
{
  // bvec represents beam, 
  // zvec represents scattered object, 
  //  vec represents a decay fragment of zvec
  // This routing transforms vec to zvec rest frame 
  // new z axis will be in the zvec direction
  // bvec and zvec will be in the new x-z plane (zvec scatters to the x direction)
  //  printf("\nnew Event\n uz,vec(E=%f),zvec(E=%f)\n",vec.E(),zvec.E());
  TVector3 uz=bvec.Vect();
  //uz.Print();
  //  vec.Vect().Print();
  //  zvec.Vect().Print();
  uz.SetMag(1.);
  Float_t phi,theta;
  phi=uz.Phi();
  uz.RotateZ(-phi);
  theta=uz.Theta();
  uz.RotateY(-theta);
  //  printf("step2\n");
  zvec.RotateZ(-phi);
  zvec.RotateY(-theta);
  vec.RotateZ(-phi);
  vec.RotateY(-theta);
  //  vec.Vect().Print();
  //  zvec.Vect().Print();
  //  printf("step3\n");
  Float_t phi2=zvec.Phi();
  vec.RotateZ(-phi2);
  zvec.RotateZ(-phi2);
  //  vec.Vect().Print();
  //zvec.Vect().Print();
  vec.Boost(-zvec.BoostVector());
  //  printf("result=");
  //  vec.Vect().Print();
  return vec.Vect();
  /*
  TVector3 vy=bvec.Vect().Cross(zvec.Vect());
  TVector3 vz=zvec.Vect();
  if(vy.Mag()<=0 || vz.Mag()<=0)return TVector3(0,0,1);
  vy.SetMag(1.);
  vz.SetMag(1.);
  TVector3 vx=vy.Cross(vz);
  TVector3 v=vec.Vect();
  TVector3 nv(v.Dot(vx),v.Dot(vy),v.Dot(vz));
  return nv;

  */

};
TLorentzVector poutTree::ClusterHardE(Int_t ClusterListNumber )
{
  TLorentzVector vsum(0,0,0,0);
  if(ClusterList->GetEntries()>ClusterListNumber)
    {
      TObjArray* cl= (TObjArray*) ClusterList->At(ClusterListNumber);
      TIter next(cl);
      while(LVec* v=(LVec*) next())
	{
	  vsum=*v+vsum;
	};
    };
  return vsum;  
};
TLorentzVector poutTree::ClusterSoftE(Int_t ClusterListNumber,Float_t maxsep)
{
  TLorentzVector vh=ClusterHardE(ClusterListNumber);
  TLorentzVector vs(0,0,0,0);
  TIter next(softlist);
  while(LVec* v=(LVec*) next())
    {
      TLorentzVector vtmp=*v;
      Float_t sep=pow(vtmp.PseudoRapidity()-vh.PseudoRapidity(),2);
      vtmp.RotateZ(-vh.Phi());
      sep+=vtmp.Phi()*vtmp.Phi();
      sep=sqrt(sep);
      if(sep<maxsep)vs=vs+*v;
    };
  return vs;
};
TObjArray* poutTree::ClusterScratch(Float_t maxsep)
{
  TIter next(scratchlist);
  ClusterList->Delete();
  while(LVec* v=(LVec*) next())
    {
      Bool_t used=false;
      TIter nxt(ClusterList);
      
      while(TObjArray* pcl=(TObjArray*) nxt())
	{
	  if(used)continue;
	  TLorentzVector vec=SumList(pcl);
	  TLorentzVector v0;
	  Float_t sep=1000.;
	  v0=*v;
	  if((v0.Vect().Mag()>0) &&( vec.Vect().Mag()>0))
	    {
	      if(ClusterYPhi)
		{
		  sep=pow(vec.PseudoRapidity()-v0.PseudoRapidity(),2);
		  v0.RotateZ(-vec.Phi());
		  sep=sep+v0.Phi()*v0.Phi();
		}
	      else
		{
		  Float_t cth=v0.Vect().Dot(vec.Vect())/v0.Vect().Mag()/vec.Vect().Mag();
		  Float_t th=acos(cth);
		  sep=th*th;
		};
	    };
	  if(sep<maxsep*maxsep)
	    {
	      pcl->Add(v);
	      used=true;
	    };
	};
      if(!used)
	{
	  TObjArray* ntob=new TObjArray(10,0);
	  ntob->Add(v);
	  ClusterList->Add(ntob);
	};
    };
  return ClusterList;
};
TLorentzVector poutTree::SumList(TObjArray* list)
{
  TLorentzVector vv(0,0,0,0);
  TIter next(list);
  while(LVec* v=(LVec*) next())
    {
      vv=vv+*v;
    }
  return vv;
};
TLorentzVector poutTree::Pair4V(LVec* lv)
{
  if(lv->Partner)
    {
      return (*lv->Partner+*lv);
    }
  else
    {
      return *lv;
    };
};
TMatrix poutTree::FillFMSADC(Int_t NSTB)
{
  return FillFMSADC(NSTB,0,0);
};
TMatrix poutTree::FillFMSADC(Int_t NSTB,Int_t startADC,Int_t nextADC)
 {
   int nrows=34;
   int ncols=17;
   if(NSTB>2)
     {
       nrows=24;
       ncols=12;
     };
   TMatrix tm(nrows,ncols);
   Int_t Startadc=startADC;
   Int_t Nextadc=nextADC;
   if(startADC==0 && nextADC==0)
   {
     Startadc=0;
     Nextadc=nSavedHits;
   };

   for(int nadc=Startadc;nadc<Nextadc;nadc++)
     {
       unsigned int s=(unsigned int) SavedHits[nadc];
       //       printf("nstb=%d s=%x \n",NSTB,s);

       Int_t sew,snstb,srow,scol,sadc;
       sew=1;
       if(s&&0x80000000)sew=2;
       //       printf("eq=%d \n",sew);
       if(sew!=2)continue;
       s=s&0x7FFFFFFF;
       snstb=((s/0x10000000)&7)+1;
       //       printf("snstb=%d \n", snstb);
       if(snstb!=NSTB)continue;
       srow= ((s/0x00400000)&0x3F)+1;
       scol= ((s/0x00010000)&0x3F)+1;
       sadc= s&0xFFF;
       tm(srow-1,scol-1)=sadc;
       //       printf("nstb=%d srow=%d scol=%d sadc=%d\n",snstb,srow,scol,sadc);

     };
   return tm;
 };
void poutTree::DrawFMSADC(Int_t iew, Int_t instb,Geom* p_geom)
{
  if(iew!=2 || instb<1 || instb>4)return;

  if(gr[instb-1])delete gr[instb-1];
  if(gradc[instb-1])delete gradc[instb-1];
  gr[instb-1]=0;
  gradc[instb-1]=0;
  TMatrix tm=FillFMSADC(instb);
  gradc[instb-1]=new TH2F(tm);
  
  gradc[instb-1]->Draw("zcol");
  gradc[instb-1]->Draw("samebox");
  ClearScratch();
  AddToScratch((iew-1)*4+instb);	      
  Int_t n3=scratchlist->GetEntries();
  Float_t xhit[10],yhit[10];
  for(int j=0;j<n3;j++)
    {
      LVec* vv= (LVec*) scratchlist->At(j);
      TVector3 v3=PosInDet(vv,p_geom,0);
      TVector3 v13=p_geom->LocalXYZ(iew,instb,v3,true);
      xhit[j]=v13.X();
      yhit[j]=v13.Y();
      printf("(E(%d)=%f x=%f y=%f) ",j,vv->E(),xhit[j],yhit[j]);
    };
  if(n3>0)
    {
      gr[instb-1]=new TGraph(n3,xhit,yhit);
      if(n3==2)printf(" m2=%f",
		      SumScratch().Mag());
      printf("\n");
      gr[instb-1]->Draw("*");
    };
};
