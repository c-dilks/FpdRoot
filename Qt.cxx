#include "Qt.h"
using namespace std;

ClassImp(Qt)


Qt::Qt(FilesSet* p_files)
{
  qtsumadcmcnt=0;
  qtmaxch[0]=578;
  qtmaxch[1]=578;
  qtmaxch[2]=288;
  qtmaxch[3]=288;

  COL_NUM[0]=COL_NUM[1]=17;
  COL_NUM[2]=COL_NUM[3]=12;
  ROW_NUM[0]=ROW_NUM[1]=34;
  ROW_NUM[2]=ROW_NUM[3]=24;
  
  qtadd_off=15;
  // qtcrate_off[crate conf num - 1] = "root12fms crate num"
  qtcrate_off[0]=0;  // L1
  qtcrate_off[1]=0;  // BC1
  qtcrate_off[2]=5;  // MXQ
  qtcrate_off[3]=0;  // MIX
  qtcrate_off[4]=0;  // BCW
  qtcrate_off[5]=0;  // BCE
  qtcrate_off[6]=6;  // FEQ
  qtcrate_off[7]=0;  // BBC DSM
  qtcrate_off[8]=7;  // BBQ (BBC QT)
  qtcrate_off[9]=0;  // FMS DSM
  qtcrate_off[10]=1; // FMS QT1
  qtcrate_off[11]=2; // FMS QT2
  qtcrate_off[12]=3; // FMS QT3
  qtcrate_off[13]=4; // FMS QT4

  memset(pr,0,sizeof(pr));
  memset(rr,0,sizeof(rr));
  memset(cr,0,sizeof(cr));
  memset(rcrate,0,sizeof(rcrate));
  memset(rslot,0,sizeof(rslot));
  memset(rdcd,-1,sizeof(rdcd));
  memset(rdch,-1,sizeof(rdch));
  memset(dcrate,0,sizeof(dcrate));
  memset(dslot,0,sizeof(dslot));
  memset(ddcd,-1,sizeof(ddcd));
  memset(ddch,-1,sizeof(ddch));

  memset(detchan,0,sizeof(detchan));
  memset(detnstb,0,sizeof(detnstb));

  memset(NTrigMap,0,sizeof(NTrigMap));
  memset(TrigMapslot,0,sizeof(TrigMapslot));
  memset(TrigMapqt8,0,sizeof(TrigMapqt8));

  memset(qtdsm0sum,0,sizeof(qtdsm0sum));
  memset(qtdsm0sum2,0,sizeof(qtdsm0sum2));
  memset(qtdsm0dcsum,0,sizeof(qtdsm0dcsum));
  memset(qtdsm0ht,0,sizeof(qtdsm0ht));
  memset(qtdsm0htadr,0,sizeof(qtdsm0htadr));
  memset(qtdsm1sum,0,sizeof(qtdsm1sum));
  memset(qtdsm1clsum,0,sizeof(qtdsm1clsum));
  memset(qtdsm1clsum_apnd,0,sizeof(qtdsm1clsum_apnd));
  memset(qtdsm1cladr,0,sizeof(qtdsm1cladr));
  memset(qtdsm1dcsuma,0,sizeof(qtdsm1dcsuma));
  memset(qtdsm1dcsumd,0,sizeof(qtdsm1dcsumd));
  memset(qtdsm1ht,0,sizeof(qtdsm1ht));
  memset(qtdsm2_trgbits,0,sizeof(qtdsm2_trgbits));
  qtdsm_trgbits=0;

  trgadj=false; //only used for trigger analysis

  qtMat=new TObjArray(4);
  for(Int_t ii=0;ii<4;ii++)
    {
      qtHist[ii]=0;
      qtMat->AddAtAndExpand(new TMatrix(ROW_NUM[ii],COL_NUM[ii]),ii);
    };
  set_mm_order();
  QtError=0;
  p_qtmap=p_files->p_qtmap();
  p_qtmap2pp=p_files->p_qtmap2pp();
  if((p_qtmap) || (p_qtmap2pp) )
    {
      if(p_qtmap->check())
	{
	  if(p_qtmap->open())
	    {
	      FILE* fp=p_qtmap->fp;
	      char line[200];
	      for(Int_t ins=1;ins<3;ins++)
		{

		  for(Int_t ip=1;ip<3;ip++)
		    {
		      for(Int_t ir=1;ir<21;ir++)
			{
			  for(Int_t ic=1;ic<17;ic++)
			    {
			      if(!feof(fp))
				{
				  char* dummy=fgets(line,200,fp);
				  if(strlen(line)<2)QtError=101;
				  if(QtError){
				    std::cout<<"readfail: "<<line<<"\n";
				    return;
				  };
				  Int_t rns,rp,rr,rc,rcrt,rslt,rch;
				  sscanf(line,"%d %d %d %d %d %d %d",
					 &rns,&rp,&rr,&rc,&rcrt,&rslt,&rch);
				  if(rns!=-1)
				    {
				      //temporary fix
				      if(rslt==12)rslt=11;
				      //end temporary
				      if(rns!=ins)QtError=102;
				      if(rp!=ip)QtError=103;
				      if(rr!=ir)QtError=104;
				      if(rc!=ic)QtError=105;
				      if(QtError)return;
				      rcrate[rc-1][rr-1][rp-1][rns-1]=rcrt;//counts from 1
				      rslot[rc-1][rr-1][rp-1][rns-1]=rslt;//counts from 1
				      rdcd[rc-1][rr-1][rp-1][rns-1]=rch/8;//counts from 0
				      rdch[rc-1][rr-1][rp-1][rns-1]=rch%8;//counts from 0
				    };
				  
				}
			      else
				{
				  QtError=110;
				  std::cout<<p_qtmap->FileName<<" read Unexpected EOF\n";
				};
			    };
			};
		    };
		};
	    }
	  else
	    {
	      QtError=120;
	    };
	}
      else
	{
	  std::cout<<p_qtmap->FileName<<" file open failed\n";
	  QtError=130;
	  return;
	};

      if(p_qtmap2pp->check())
	{
	  if(p_qtmap2pp->open())
	    {
	      FILE*fp=p_qtmap2pp->fp;
	      char line[200];
	      for(Int_t i=1;i<5;i++)//only for fms
		{
		  Int_t m=mm[i-1];//order in which qtmap2pp is organized 
		  for(Int_t j=1;j<=qtmaxch[m-1];j++)
		    {
		      if(!feof(fp))
			{
			  if(!fgets(line,200,fp))QtError=201;
			  if(QtError)return;
			  Int_t modr,chr;
			  sscanf(line,"%d %d %d %d %d",&modr,&chr,
				&( pr[j-1][m-1]),&(rr[j-1][m-1]),
				 &(cr[j-1][m-1]));
			  if(m==modr && j==chr)
			    {
			      Int_t pr0,rr0,cr0,chan0,slot0,crt0,ns0;
			      pr0=pr[j-1][m-1];
			      rr0=rr[j-1][m-1];
			      cr0=cr[j-1][m-1];
			      ns0=(m-1)%2+1;
			      
			      dcrate[j-1][m-1]=0;
			      dslot[j-1][m-1]=0;
			      ddcd[j-1][m-1]=-1;
			      ddch[j-1][m-1]=-1;
			      if(cr0>0 && rr0>0 && pr0>0)
				{
				  dcrate[j-1][m-1]=rcrate[cr0-1][rr0-1][pr0-1][ns0-1];
				  dslot[j-1][m-1]=rslot[cr0-1][rr0-1][pr0-1][ns0-1];
				  ddcd[j-1][m-1]=rdcd[cr0-1][rr0-1][pr0-1][ns0-1];
				  ddch[j-1][m-1]=rdch[cr0-1][rr0-1][pr0-1][ns0-1];
				  
				  chan0=rdcd[cr0-1][rr0-1][pr0-1][ns0-1]*8+rdch[cr0-1][rr0-1][pr0-1][ns0-1];
				  slot0=rslot[cr0-1][rr0-1][pr0-1][ns0-1]-1;
				  crt0=rcrate[cr0-1][rr0-1][pr0-1][ns0-1]-1;
				  detchan[chan0][slot0][crt0]=chr;
				  detnstb[slot0][crt0]=m;
				  //cout<<"crt="<<crt0+1<<" slot="<<slot0+1<<" qtch="<<chan0<<" detchan="<<detchan[chan0][slot0][crt0]<<endl;
				};
			    }else
			      {
				QtError=202;
				//std::cout<<"m="<<m<<" j="<<j<<"\n";
				//std::cout<<line<<"\n";
			      };
			  if(QtError)return;
			}
		      else
			{
			  QtError=203;
			  std::cout<<p_qtmap2pp->FileName<<" read Unexpected EOF\n";
			  return;
			};		     
		    };
		};
	    };
	}
      else
	{
	  std::cout<<p_qtmap2pp->FileName<<" open failed\n";
	  QtError=204;
	  return;
	};
    }
  else
    {
      QtError=205;
      std::cout<<"FilesSet input to Qt does not include qtmap or qtmap2pp \n";
      return;
    };

  //DefineTrigMap();
  DefineRealTrigMap();

};

void Qt::initQT()
{
  for(Int_t i=1;i<5;i++)
    {
      for(Int_t j=1;j<17;j++)
	{
	  for(Int_t k=1;k<33;k++)
	    {
	      qtadc[k-1][j-1][i-1]=-1;
	      qttac[k-1][j-1][i-1]=-1;
	      if(k<5)
		{
		  qt8sum[k-1][j-1][i-1]=0;
		};
	    };
	  qtsum[j-1][i-1]=0;
	  qtht[j-1][i-1]=0;
	  qthtadr[j-1][i-1]=-1;
	};
      fmssum[i-1]=0;
    };
  fmssum[4]=0;
  errfff=0;
};

Bool_t Qt::decodeQT(Int_t nQTdata,UInt_t* QTdata, Int_t IDNEVT)
{
  initQT();
  if(nQTdata==0){return false;}
  else if(QTdata[nQTdata-1]!=0xAC10){return false;};
  Int_t nboard=0;
  ndata=0;
  sumadc=0;
  mcnt=0;
  Int_t crate,crt,add,addr,nentry,j;
  Bool_t header=true;

  for(Int_t i=0;i<nQTdata-1;i++)
    {
      Int_t dat1,dat2;
      if(header)
	{
	  dat1=QTdata[i]/0x10000;
	  dat2=QTdata[i]%0x10000;
	  crate=dat1/0x100;
	  addr=dat1%0x100;
	  crt=qtcrate_off[crate-1];//counts from 1
	  add=addr-qtadd_off;//counts from 1
	  nentry=dat2%0x100;
	  nboard=nboard+1;
	  j=0;
	  if(nentry>0)header=false;
	  //	  cout<<"head nqt="<<i<<" QTdata="<<hex<<QTdata[i]<<dec<<" crt(crate)="<<crt<<"("<<crate<<") slot="<<add<<" nline="<<nentry<<endl;
	}
      else
	{
	  j++;
	  ndata++;
	  dat1=QTdata[i]/0x10000;
	  dat2=QTdata[i]%0x10000;
	  Int_t adc=dat2%0x1000;
	  sumadc+=adc;
	  Int_t tdc=dat1%0x20;
	  Int_t dc=dat1/0x4000;
	  Int_t dcc=(dat1/0x800)%8;
	  Int_t ch=dcc+dc*8;    //counts from 0
	  Int_t unused=((dat1%0x800)/32)*8+dat2/0x1000;
	  Int_t qt8=ch/8;
	  //	  if(crate==9)cout<<"nqt="<<i<<hex<<" QTdata="<<QTdata[i]<<dec<<" dcard="<<dc<<" dcardch="<<dcc<<" tdc="<<tdc<<" adc="<<hex<<adc<<endl;
	  if((QTdata[i]==(UInt_t) (-1)) || (unused|=0))
	    {
	      errfff++;
	    };
	  Float_t trgR=1.;
	  if(trgadj)
	    {
	      trgR=GetTrgAdj(ch,add,crt);
	    };
	  adc=(Int_t)(trgR*adc);
	  qtadc[ch][add-1][crt-1]=adc;
	  if(trgR==0.)printf("ch=%d slt=%d crt=%d adc=%d adjadc=%d \n",ch,add,crt,adc,qtadc[ch][add-1][crt-1]);
	  if(adc>10 && crt<5){mcnt++;};
	  qttac[ch][add-1][crt-1]=tdc;
	  qtsum[add-1][crt-1]+=adc;	  
	  qt8sum[qt8][add-1][crt-1]+=adc;
	  //std::cout<< "crt="<<crt;
	  //std::cout<< " add="<<add;
	  //std::cout<< " ch="<<ch;
	  //std::cout<< " adc="<<adc<<"\n";;
	  if(adc/32>qtht[add-1][crt-1]/32)
	    {
	      qtnht[add-1][crt-1]=1;
	      qtht[add-1][crt-1]=adc;
	      qthtadr[add-1][crt-1]=ch;
	      //cout<<"crt="<<crt-1<<" slot="<<add-1<<" ch="<<ch<<" qtht="<<adc<<endl;
	    }
	  else if(adc/32>0 && adc/32==qtht[add-1][crt-1]/32)
	    {
	      qtnht[add-1][crt-1]++;
	      if(ch/8 <= qthtadr[add-1][crt-1]/8)
		{
		  qtht[add-1][crt-1]=adc;
		  qthtadr[add-1][crt-1]=ch;
		};
	    };
	  if(add!=12)
	    {
	      fmssum[crt-1]+=adc;
	      fmssum[4]+=adc;
	    };
	  if(j==nentry)header=true;	  
	};  
    };
  NFmsAdc=0;
  FmsSumAdc=0;
  for(Int_t m=0;m<4;m++)
    {      
      Bool_t flag=false;
      for(Int_t ch=0;ch<qtmaxch[m];ch++)
	{
	  Int_t p=pr[ch][m];
	  Int_t r=rr[ch][m];
	  Int_t c=cr[ch][m];
	  if(p>0 && p<3 && c>0 && r>0 && c<17 && r<21)
	    {
	      Int_t ins=m%2;	      
	      Int_t qt=rcrate[c-1][r-1][p-1][ins];
	      Int_t slot=rslot[c-1][r-1][p-1][ins];
	      Int_t dcd=rdcd[c-1][r-1][p-1][ins];
	      Int_t dch=rdch[c-1][r-1][p-1][ins];
	      Int_t row=ch/COL_NUM[m];
	      Int_t col=ch%COL_NUM[m];
	      Int_t adcval=qtadc[dcd*8+dch][slot-1][qt-1];
	      if(adcval<0)adcval=0;
	      (*((TMatrix*) (qtMat->At(m)) ))(row,col)=adcval;
	      if(adcval>10)NFmsAdc++;
	      FmsSumAdc+=adcval;
	      //cout<<"decode nstb="<<m+1<<" row="<<row<<" col="<<col<<" adc="<<adcval<<endl;
	      //if(qtHist[m] && sumadc<50000)
	      //if(qtHist[m])
	      if(qtHist[m] && FmsSumAdc<50000)
		{
		  if(sumadc>5000000)printf("sumadc large %d \n",sumadc);
		  qtHist[m]->Fill(ch,adcval);
		};

	      //std::cout<<"adc="<<qtadc[dcd*8+dch][slot-1][qt-1]<<"\n";
	    };
	};
      //if(flag)(*((TMatrix*) (qtMat->At(m)) )).Print();
    };
  if(qtsumadcmcnt)
    {
      Float_t xmcnt=NFmsAdc;
      Float_t xsumadc=FmsSumAdc;
      if(xmcnt>399)xmcnt=399;
      if(xsumadc>99999)xsumadc=99999;
      qtsumadcmcnt->Fill(xsumadc,xmcnt);
    };
  return true;
};

Bool_t Qt::encodeQT(TMatrix** m_adc,UInt_t* nQTdata,UInt_t* QTdata)
{

  Int_t nline[4][16];
  memset(nline, 0, sizeof(nline));

  for(Int_t iNSTB=0; iNSTB<4; iNSTB++)
    {
      Int_t nrows = m_adc[iNSTB]->GetNrows();
      Int_t ncols = m_adc[iNSTB]->GetNcols();
	
      for(Int_t i=0; i<nrows; i++)
	{
	  for(Int_t j=0; j<ncols; j++)
	    {
	      Int_t ch = i*ncols + j;
	      Int_t crt = dcrate[ch][iNSTB]-1;
	      Int_t slot = dslot[ch][iNSTB]-1;
	      if(crt<0 || slot<0) continue;
	      if((*(m_adc[iNSTB]))(i,j) > 0)
		nline[crt][slot]++;
	    }
	}
    }

  Int_t nqt = 0;
  for(Int_t iCRT=0; iCRT<4; iCRT++)
    {
      for(Int_t iSLT=0; iSLT<16; iSLT++)
	{
	  //cout<<"crt="<<iCRT+1<<" slot="<<iSLT+1<<" detnstb="<<detnstb[iSLT][iCRT]<<endl;
	  Int_t iiCRT = iCRT;
	  Int_t iiSLT = iSLT;
	  if(nline[iCRT][iSLT]==0) continue;
	  QTdata[nqt] = ((iiCRT+11) << 24) + ((iiSLT+16) << 16) + nline[iCRT][iSLT];
	  nqt++;

	  //cout<<"head nqt="<<nqt-1<<" QTdata="<<QTdata[nqt-1]<<" crt="<<iiCRT+1<<" slot="<<iiSLT+1<<" nline="<<nline[iCRT][iSLT]<<endl;
	  Int_t nstb = detnstb[iSLT][iCRT]-1;
	  if(nstb<0) continue;
	  Int_t nrows = m_adc[nstb]->GetNrows();
	  Int_t ncols = m_adc[nstb]->GetNcols();
	  for(Int_t iCh=0; iCh<32; iCh++)
	    {
	      //Int_t ch = 0;
	      Int_t ch = detchan[iCh][iSLT][iCRT]-1;
	      if(ch<0) continue;
	      //cout<<"crt="<<iCRT+1<<" slot="<<iSLT+1<<" qtch="<<iCh<<endl;
	      Int_t dcard = iCh/8;
	      Int_t dch = iCh%8;
	      Int_t row = ch/ncols;
	      Int_t col = ch%ncols;
	      Int_t tdc = 0;
	      //Int_t adc = 0;
	      Int_t adc = (*(m_adc[nstb]))(row,col);
	      //cout<<"encode nstb="<<nstb+1<<" row="<<row<<" col="<<col<<" adc="<<adc<<endl;
	      //cout<<"crt="<<iCRT+1<<" slot="<<iSLT+1<<" qtch="<<iCh<<" dcard="<<dcard<<" dch="<<dch<<" nstb="<<nstb+1<<" ch="<<ch+1<<" row="<<row<<" col="<<col<<" adc="<<adc<<endl;
	      if(adc<=0) continue;
	      QTdata[nqt] = (dcard << 30) + (dch << 27) + (tdc << 16) + adc;
	      nqt++;
	      //cout<<"nqt="<<nqt-1<<" QTdata="<<QTdata[nqt-1]<<" dcard="<<dcard<<" dcardch="<<dch<<" tdc="<<tdc<<" adc="<<adc<<endl;
	    }
	}
    }
  QTdata[nqt] = 0xAC10;
  *nQTdata = nqt+1;

  return true;
}

TMatrix* Qt::tm(Int_t instb)
{
  if(instb>=0 && instb<4)
    {
      return ((TMatrix*) qtMat->At(instb));
    };
  return 0;
};

void Qt::EnableqtHist()
{
  qtHist[0]=new TH2F("qtHistLN","qtHistLN",578,0,578,4096,0,4096);
  qtHist[1]=new TH2F("qtHistLS","qtHistLS",578,0,578,4096,0,4096);
  qtHist[2]=new TH2F("qtHistSN","qtHistSN",288,0,288,4096,0,4096);
  qtHist[3]=new TH2F("qtHistSS","qtHistSS",288,0,288,4096,0,4096);
  //qtsumadcmcnt=new TH2F("qtsumadcmcnt","qtsumadcmcnt",1000,0,100000,400,0,400);
};

Qt::~Qt()
{
  for(Int_t i=0;i<4;i++)
    {
      if(qtHist[i])delete qtHist[i];
    };
  if(qtMat!=0)
    {
      qtMat->Delete();
      delete qtMat;
    };
  
};

void Qt::DefineTrigMap()
{  
  AddToTrigMapSq("A",1,5,"A0",2);
  AddToTrigMapSq("A",9,13,"A0",3);
  AddToTrigMapSq("A",17,21,"A1",3);
  AddToTrigMapSq("A",25,29,"A2",3);

  AddToTrigMapSq("B",1,5,"A3",3);
  AddToTrigMapSq("B",10,14,"B0",3);
  AddToTrigMap("B",16,16,"B1","B2","B3","D0");
  AddToTrigMapSq("B",17,22,"B1",3);
  AddToTrigMap("B",24,24,"B2","B3","C0","D0");
  AddToTrigMapSq("B",25,30,"B2",3);

  AddToTrigMap("C",0,0,"B3","C0","C1","D0");
  AddToTrigMapSq("C",1,6,"B3",3);
  AddToTrigMap("C",8,8,"C0","C1","C2","D0");
  AddToTrigMapSq("C",9,14,"C0",3);
  AddToTrigMapSq("C",16,16,"C1",4);
  AddToTrigMapSq("C",17,22,"C1",3);

  AddToTrigMapSq("D",1,1,"D0",2);
  AddToTrigMapSq("D",2,2,"D0",2);
  AddToTrigMapSq("D",3,3,"D0",2);
  AddToTrigMapSq("D",4,4,"D0",2);
  AddToTrigMapSq("D",5,5,"D0",2);
  AddToTrigMapSq("D",9,13,"D0",3);
  AddToTrigMapSq("D",17,21,"D1",3);
  AddToTrigMapSq("D",25,29,"D2",2);

  AddToTrigMapSq("E",9,14,"E0",3);
  AddToTrigMapSq("E",17,22,"E1",3);
  AddToTrigMapSq("E",25,30,"E2",3);

  AddToTrigMapSq("F",1,6,"E3",3);
  AddToTrigMapSq("F",9,14,"F0",3);
  AddToTrigMapSq("F",17,22,"F1",3);
  AddToTrigMapSq("F",25,25,"F2",2);
  AddToTrigMapSq("F",26,26,"F2",2);
  AddToTrigMapSq("F",27,27,"F2",2);
  AddToTrigMapSq("F",28,28,"F2",2);
  AddToTrigMapSq("F",29,29,"F2",2);
  AddToTrigMapSq("F",30,30,"F2",2);

  AddToTrigMapSq("G",10,10,"G0",3);
  AddToTrigMap("G",11,11,"G0","G1","G2","F3");
  AddToTrigMapSq("G",17,18,"G1",3);
  AddToTrigMap("G",19,19,"G1","G2","G3","F3");
  AddToTrigMapSq("G",25,26,"G2",3);
  AddToTrigMap("G",27,27,"G2","G3","H0","F3");

  AddToTrigMap("H",6,6,"F3","G3","H0","H1");
  AddToTrigMapSq("H",4,5,"G3",3);
  AddToTrigMap("H",8,8,"F3","H0","H1","H2");
  AddToTrigMapSq("H",9,13,"H0",3);
  AddToTrigMap("H",16,16,"F3","H1","H2","H3");
  AddToTrigMapSq("H",17,22,"H1",3);
  AddToTrigMapSq("H",25,30,"H2",2);

  AddToTrigMapSq("I",1,6,"I0",2);
  AddToTrigMapSq("I",9,14,"I0",3);
  AddToTrigMapSq("I",17,22,"I1",3);
  AddToTrigMapSq("I",25,30,"I2",3);

  AddToTrigMapSq("J",1,6,"I3",3);
  AddToTrigMapSq("J",9,14,"J0",3);
  AddToTrigMapSq("J",17,22,"J1",3);
  AddToTrigMapSq("J",25,30,"J2",2);

}; 

void Qt::DefineRealTrigMap()
{  
  AddToTrigMapSq("A",0,6,"A0",2);
  AddToTrigMapSq("A",8,14,"A0",3);
  AddToTrigMapSq("A",16,22,"A1",3);
  AddToTrigMapSq("A",24,30,"A2",3);

  AddToTrigMapSq("B",0,6,"A3",3);
  AddToTrigMap("B",8,8,"B0","B1","B2","D0");
  AddToTrigMapSq("B",9,15,"B0",3);
  AddToTrigMap("B",16,16,"B1","B2","B3","D0");
  AddToTrigMapSq("B",17,23,"B1",3);
  AddToTrigMap("B",24,24,"B2","B3","C0","D0");
  AddToTrigMapSq("B",25,31,"B2",3);

  AddToTrigMap("C",0,0,"B3","C0","C1","D0");
  AddToTrigMapSq("C",1,7,"B3",3);
  AddToTrigMap("C",8,8,"C0","C1","C2","D0");
  AddToTrigMapSq("C",9,15,"C0",3);
  AddToTrigMapSq("C",16,16,"C1",4);
  AddToTrigMapSq("C",17,23,"C1",3);

  AddToTrigMapSq("D",0,0,"D0",2);
  AddToTrigMapSq("D",1,1,"D0",2);
  AddToTrigMapSq("D",2,2,"D0",2);
  AddToTrigMapSq("D",3,3,"D0",2);
  AddToTrigMapSq("D",4,4,"D0",2);
  AddToTrigMapSq("D",5,5,"D0",2);
  AddToTrigMapSq("D",6,6,"D0",2);
  AddToTrigMapSq("D",8,14,"D0",3);
  AddToTrigMapSq("D",16,22,"D1",3);
  AddToTrigMapSq("D",24,31,"D2",2);

  AddToTrigMapSq("E",8,15,"E0",3);
  AddToTrigMapSq("E",16,23,"E1",3);
  AddToTrigMapSq("E",24,31,"E2",3);

  AddToTrigMapSq("F",0,7,"E3",3);
  AddToTrigMapSq("F",8,15,"F0",3);
  AddToTrigMapSq("F",16,23,"F1",3);
  AddToTrigMapSq("F",24,24,"F2",2);
  AddToTrigMapSq("F",25,25,"F2",2);
  AddToTrigMapSq("F",26,26,"F2",2);
  AddToTrigMapSq("F",27,27,"F2",2);
  AddToTrigMapSq("F",28,28,"F2",2);
  AddToTrigMapSq("F",29,29,"F2",2);
  AddToTrigMapSq("F",30,30,"F2",2);
  AddToTrigMapSq("F",31,31,"F2",2);

  AddToTrigMapSq("G",8,10,"G0",3);
  AddToTrigMap("G",11,11,"G0","G1","G2","F3");
  AddToTrigMapSq("G",16,18,"G1",3);
  AddToTrigMap("G",19,19,"G1","G2","G3","F3");
  AddToTrigMapSq("G",24,26,"G2",3);
  AddToTrigMap("G",27,27,"G2","G3","H0","F3");

  AddToTrigMapSq("H",0,5,"G3",3);
  AddToTrigMap("H",6,6,"F3","G3","H0","H1");
  AddToTrigMap("H",8,8,"F3","H0","H1","H2");
  AddToTrigMapSq("H",9,15,"H0",3);
  AddToTrigMap("H",16,16,"F3","H1","H2","H3");
  AddToTrigMapSq("H",17,23,"H1",3);
  AddToTrigMapSq("H",24,31,"H2",2);

  AddToTrigMapSq("I",0,7,"I0",2);
  AddToTrigMapSq("I",8,15,"I0",3);
  AddToTrigMapSq("I",16,23,"I1",3);
  AddToTrigMapSq("I",24,31,"I2",3);

  AddToTrigMapSq("J",0,7,"I3",3);
  AddToTrigMapSq("J",8,15,"J0",3);
  AddToTrigMapSq("J",16,23,"J1",3);
  AddToTrigMapSq("J",24,31,"J2",2);

}; 

void Qt::AddToTrigMap(const char* slot,Int_t ch_min,Int_t ch_max, 
		      const char* s0,const char* s1,const char* s2,const char* s3)
{
  const char* s[4];
  s[0]=s0; s[1]=s1; s[2]=s2; s[3]=s3; 
  const char* sA="A";
  const char* szero="0";
  Int_t ht_slt=(Int_t)(*slot)-(Int_t)(*sA);

  for(Int_t ch=ch_min;ch<ch_max+1;ch++)
    {
      Int_t nqt8=0;
      for(Int_t j=0;j<4;j++)
	{
	  if(s[j]!="")
	    {
	      Int_t cl_slt=(Int_t)(*(s[j]))-(Int_t)(*sA);//slot for qt8 in the cluster sum
	      Int_t cl_qt8=(Int_t)(*(s[j]+1))-(Int_t)(*szero);
	      TrigMapslot[ch][ht_slt][j]=cl_slt;
	      TrigMapqt8[ch][ht_slt][j]=cl_qt8;
	      //cout<<"ch="<<ch<<" ht_slt="<<ht_slt<<" j="<<j<<" cl_slt="<<cl_slt<<" cl_qt8="<<cl_qt8<<endl;
	      nqt8++;
	    };
	};
      NTrigMap[ch][ht_slt]=nqt8;

      //cout<<"ch="<<ch<<" ht_slt="<<ht_slt<<"nqt8="<<nqt8<<endl;
    };

}; 

void Qt::AddToTrigMapSq(const char* slot,Int_t ch_min,Int_t ch_max, 
		      const char* s0,Int_t nqt8)
{ 
  const char* sA="A";
  const char* szero="0";
  Int_t ht_slt=(Int_t)(*slot)-(Int_t)(*sA);
  char ss[2];
  ss[0]=*s0;
  ss[1]=*(s0+1);
  Int_t icl_slt=(Int_t)(ss[0])-(Int_t)(*sA);
  Int_t icl_qt8=(Int_t)(ss[1])-(Int_t)(*szero);
    
  for(Int_t ch=ch_min;ch<ch_max+1;ch++)
    {
      NTrigMap[ch][ht_slt]=nqt8; 
      Int_t ii=icl_slt*4+icl_qt8;
      for(Int_t j=0;j<nqt8;j++)
	{	  		  
	  Int_t cl_slt=ii/4;
	  Int_t cl_qt8=ii%4; 
	  //printf("ch %d  j %d  cl_slt %d  cl_qt8 %d  ii%d  ht_slt %d \n",ch,j,cl_slt,cl_qt8,ii,ht_slt);
	  TrigMapslot[ch][ht_slt][j]=cl_slt;
	  TrigMapqt8[ch][ht_slt][j]=cl_qt8;
	  ii++;
	};
    };
}; 

Float_t Qt::GetTrgAdj(Int_t ch, Int_t slot, Int_t crate)
{
  Int_t chan=detchan[ch][slot-1][crate-1]-1;
  Int_t nstb=detnstb[slot-1][crate-1];
  Float_t trgR=1.;
  if(chan>-1 && nstb>0)
    {
      Int_t row=chan/COL_NUM[nstb-1];
      Int_t col=chan%COL_NUM[nstb-1];     
      //printf("row=%d col=%d \n",row,col);
      
      if(nstb==3 && row==5 && col==1)trgR=.5;
      if(nstb==3 && row==5 && col==0)trgR=.7;
      if(nstb==3 && row==4 && col==1)trgR=.7;
      if(nstb==3 && row==4 && col==0)trgR=.7;
      if(nstb==3 && row==10 && col==6)trgR=.7;
      if(nstb==3 && row==12 && col==5)trgR=.5;
      if(nstb==3 && row==14 && col==11)trgR=.001;
      if(nstb==3 && row==9 && col==9)trgR=.5;
      //if(nstb==3 && row==3 && col==2)trgR=.7;
      //if(nstb==3 && row==4 && col==1)trgR=.7;
      
      if(nstb==4 && row==5 && col==0)trgR=.9;
      if(nstb==4 && row==5 && col==1)trgR=.9;
      if(nstb==4 && row==5 && col==2)trgR=.9;
      if(nstb==4 && row==5 && col==3)trgR=.9;
      if(nstb==4 && row==4 && col==0)trgR=.7;
      if(nstb==4 && row==4 && col==1)trgR=.7;
      //if(nstb==4 && row==4 && col==1)trgR=.7;
      //if(nstb==4 && row==4 && col==2)trgR=.7;
      if(nstb==4 && row==3 && col==0)trgR=.5;
      //if(nstb==4 && row==0 && col==1)trgR=.5;
      if(nstb==4 && row==7 && col==6)trgR=.8;
      if(nstb==4 && row==6 && col==10)trgR=.5;
      if(nstb==4 && row==11 && col==6)trgR=.8;
      if(nstb==4 && row==12 && col==5)trgR=.8;
      if(nstb==4 && row==12 && col==6)trgR=.8;
      if(nstb==4 && row==15 && col==7)trgR=.5;
      if(nstb==4 && row==18 && col==0)trgR=.8;
      if(nstb==4 && row==21 && col==6)trgR=.5;
      if(nstb==4 && row==22 && col==1)trgR=.5;
      if(nstb==4 && row==1 && col==0)trgR=.7;
      if(nstb==4 && row==1 && col==1)trgR=.7;
      if(nstb==4 && row==2 && col==1)trgR=.8;
      if(nstb==4 && row==2 && col==2)trgR=.8;
      if(nstb==4 && row==2 && col==4)trgR=.9;
    }
  else
    {
      //printf("ch=%d slot=%d crate=%d chan=%d nstb=%d \n",ch,slot,crate,chan,nstb);
    };

  //if(trgR!=1.)printf("trgR = %f \n",trgR); 
 
  return trgR;
};

Bool_t Qt::decodeDSM(UChar_t* DSMdata, UInt_t* dsm2, UInt_t dsm)
{
  Int_t dsm0map_bd[40] = {1, 1, 1, 1, 6, 6, 6, 6, 7, 7, //10xQT1 qt32 boards map to L0dsm
			  2, 2, 2, 2, 8, 8, 8, 8, 9, 9, //10xQT2 qt32 boards map to L0dsm
			  3, 3, 3, 3,11,11,11,11,12,12, //10xQT3 qt32 boards map to L0dsm
			  4, 4, 4, 4,13,13,13,13,14,14};//10xQT4 qt32 boards map to L0dsm
  Int_t dsm0map_add[40]= {3, 4, 1, 2, 3, 4, 1, 2, 1, 2, //qt32 board index within L0dsm
			  3, 4, 1, 2, 3, 4, 1, 2, 1, 2,
			  3, 4, 1, 2, 3, 4, 1, 2, 1, 2,
			  3, 4, 1, 2, 3, 4, 1, 2, 1, 2};
  Int_t dsm1map_bd[12] = {5, 5, 5, 5,10,10,10,10,15,15,15,15};//12xL0dsm boards map to L1dsm
  Int_t dsm1map_add[12]= {2, 1, 4, 3, 2, 1, 4, 3, 2, 1, 4, 3};//L0dsm index within L1dsm
  
  //L0dsm inputs
  for(Int_t i=0;i<4;i++)//crate
    {
      for(Int_t j=0;j<10;j++)//slot (only 10/16 in the trigger)
	{
	  Int_t k = 16*(dsm0map_bd[j+i*10]-1)+4*(dsm0map_add[j+i*10]-1);
	  qtdsm0sum[j][i]  = DSMdata[k+3] + DSMdata[k+2]*256;
	  qtdsm0sum2[j][i] = DSMdata[k+1] + DSMdata[k]*256;

	  qtdsm0dcsum[0][j][i] = DSMdata[k+3]%32;
	  qtdsm0dcsum[1][j][i] = DSMdata[k+3]/32  + (DSMdata[k+2]%4)*8; 
	  qtdsm0dcsum[2][j][i] =(DSMdata[k+2]/4)%32; 
	  qtdsm0dcsum[3][j][i] = DSMdata[k+2]/128 + (DSMdata[k+1]%16)*2;
	  qtdsm0ht[j][i]    = DSMdata[k+1]/16  + (DSMdata[k]%8)*16;
	  qtdsm0htadr[j][i] = DSMdata[k]/8;
	};
    };

  //L1dsm inputs
  for(Int_t i=0;i<3;i++)//L1dsm: All_small, N_large, S_large
    {
      for(Int_t j=0;j<4;j++)//4 L0dsm per L1dsm 
	{
	  Int_t k = i*4 + j;            
	  Int_t l = 16*(dsm1map_bd[k]-1)+4*(dsm1map_add[k]-1);
	  qtdsm1sum[0][j][i] = DSMdata[l+3] + DSMdata[l+2]*256;
	  qtdsm1sum[1][j][i] = DSMdata[l+1] + DSMdata[l]*256;
	  qtdsm1clsum[k]  = qtdsm1sum[0][j][i]%256;
	  qtdsm1cladr[k]  = (qtdsm1sum[0][j][i]/256)%128;
	  qtdsm1dcsuma[k] = qtdsm1sum[1][j][i]%32;
	  qtdsm1dcsumd[k] = (qtdsm1sum[1][j][i]/32)%32;
	  qtdsm1ht[k]     = (qtdsm1sum[1][j][i]/32/32)%2;
	};
      ///L1 internal cluster sums  
      for(Int_t j=0;j<4;j++)//4 L0dsm per L1dsm 
	{
	  Int_t k = i*4 + j;    
	  qtdsm1clsum_apnd[k]=qtdsm1clsum[k];
	  Int_t ch32 =qtdsm1cladr[k]%32;
	  Int_t slot4=qtdsm1cladr[k]/32;
	  if(i==0 && slot4==0 && ch32>0   && ch32<6 ){qtdsm1clsum_apnd[k]+=qtdsm1dcsuma[(k+2)%4];};
	  if(i==0 && slot4==3 && ch32>24  && ch32<30){qtdsm1clsum_apnd[k]+=qtdsm1dcsumd[k+(Int_t)pow(-1,k)];};
	  if(i>0 && (j==0 || j==2) && slot4==3 && ch32>24  && ch32<31){qtdsm1clsum_apnd[k]+=qtdsm1dcsuma[k+1];};
	  if(i>0 && (j==1 || j==3) && slot4==0 && ch32>0   && ch32<7 ){qtdsm1clsum_apnd[k]+=qtdsm1dcsumd[k-1];};
	  if(i>0 && (j==1 || j==3) && slot4==1 && ch32>24  && ch32<31){qtdsm1clsum_apnd[k]+=qtdsm1dcsumd[i*4+(j+2)%4];};
	};
    };  
  
  ///L2dsm inputs
  //for(Int_t i=0; i<8; i++)
  //{
  //cout<<hex<<dsm2[i]<<" ";
  //}
  //cout<<dec<<endl;
  qtdsm2_trgbits[0] = dsm2[3];
  qtdsm2_trgbits[1] = dsm2[1];
  qtdsm2_trgbits[2] = dsm2[7];

  //L2dsm output bits
  //qtdsm_trgbits = dsm&0x03FFF;
  qtdsm_trgbits = dsm;

  return true;
};

void Qt::printDSM()
{

  printf("\ncrate   slot   qtdsm0ht   qtdsm0htadr   sum0   sum1   sum2   sum3\n");
  for(Int_t i=0; i<4; i++)
    {
      for(Int_t j=0; j<10; j++)
	{
	  printf("%5d   %4d   %8d   %11d   %4d   %4d   %4d   %4d\n",i,j,qtdsm0ht[j][i],
		 qtdsm0htadr[j][i],qtdsm0dcsum[0][j][i],qtdsm0dcsum[1][j][i],
		 qtdsm0dcsum[2][j][i],qtdsm0dcsum[3][j][i]);
	}
    }
  //printf("dsm1   dsm0   qt   qtdsm1clsum   qtdsm1htbit   qtdsm1htadr   qtdsm1clsum_apnd\n");
  printf("dsm1   dsm0   qt   qtdsm1clsum   qtdsm1htbit   qtdsm1htadr\n");
  for(Int_t i=0; i<12; i++)
    {
      Int_t dsm1 = i/4;
      Int_t dsm0 = i%4;
      //printf("%4d   %4d   %2d   %11d   %8d   %11d   %16d\n",dsm1,dsm0,qtdsm1cladr[i]/32,qtdsm1clsum[i],
      //qtdsm1ht[i],qtdsm1cladr[i]%32,qtdsm1clsum_apnd[i]);
      if(qtdsm1clsum[i])printf("%4d   %4d   %2d   %11d   %8d   %11d\n",dsm1,dsm0,qtdsm1cladr[i]/32,qtdsm1clsum[i],
	     qtdsm1ht[i],qtdsm1cladr[i]%32);
    }

  //cout<<"L1 trigger bits: "<<hex<<qtdsm2_trgbits[0]<<" "<<qtdsm2_trgbits[1]<<" "<<qtdsm2_trgbits[2]<<dec<<endl;
  //cout<<"L2 trigger bits: "<<hex<<qtdsm_trgbits<<dec<<endl;
}

void Qt::printQt()
{
  printf("\ncrate   slot   qtdsm0ht   qtdsm0htadr   sum0   sum1   sum2   sum3\n");
  for(int i=0; i<4; i++)
    {
      for(int j=0; j<10; j++)
	{
	  printf("%5d   %4d   %8d   %11d   %4d   %4d   %4d   %4d\n",i,j,qtht[j][i]/32,qthtadr[j][i],
		 qt8sum[0][j][i]/32<=31?qt8sum[0][j][i]/32:31,
		 qt8sum[1][j][i]/32<=31?qt8sum[1][j][i]/32:31,
		 qt8sum[2][j][i]/32<=31?qt8sum[2][j][i]/32:31,
		 qt8sum[3][j][i]/32<=31?qt8sum[3][j][i]/32:31);
	}
    }
}
