#include "Fill.h"
ClassImp(Fill)

Fill::~Fill()
{
  if(pFillArray)delete pFillArray;
};

Fill::Fill(FilesSet* fset,Int_t firstFill=7010,Int_t lastFill=7327,
	   Int_t firstRun=0,Int_t lastRun=100000000)
{
  p_Files=fset;
  FillFile=fset->p_FillFile();
  pfd=0;
  p_Files=fset;
  FillNumber=0;
  FirstRun=7000000;
  if(FirstRun<firstRun)FirstRun=firstRun;
  FirstFill=firstFill;
  LastFill=lastFill;
  FILE* fp;
  char line[200];
  int len=100;
  error=0;
  BlueCnt=0;
  YellowCnt=0;
  ErrorMessage="";
  CurrentRunNumber=-1;
  pFillArray=new TObjArray(FirstFill,LastFill-FirstFill+1);
  FillData* pFillData;
  Filentry* BluePatternFile=p_Files->p_BluePattern();
  Filentry* YellowPatternFile=p_Files->p_YellowPattern();
  Int_t tmp;
  if(FillFile->open())
    {
      fp=FillFile->fp;
      while(!feof(fp))
	{
	  if(fgets(line,len,fp)>0)
	    {
	      Int_t fillNum=0;
	      Int_t runnum=0;
	      sscanf(line,"%d %d", &runnum, &fillNum);
	      if((fillNum>=FirstFill) && (fillNum<=LastFill)
		 &&(runnum>=firstRun) && (runnum<=lastRun) )
		{
		  if(runnum<=FirstRun)FirstRun=runnum;
		  if(runnum>=LastRun)LastRun=runnum;
		  if((*pFillArray)[fillNum]==0)
		    {
		      pFillArray->AddAt(pFillData=new FillData(fillNum),fillNum);
		      printf("l\n");
		      
		    }
		  else
		    {
		      pFillData=(FillData*) ((*pFillArray)[fillNum]);
		    };
		  if((pFillData->FirstRun)>runnum)pFillData->FirstRun=runnum;
		  if((pFillData->LastRun)<runnum)pFillData->LastRun=runnum;
		  if(!(pFillData->BlueGood))
		    {
		      BluePatternFile->GetRunFillBlueName(runnum,fillNum);
		      
		      pFillData->BluePatternFile=BluePatternFile;
		      FILE* blue;
		      if(BluePatternFile->open())
			{
			  blue=BluePatternFile->fp;
			  //  printf("blue open\n");
			  BlueCnt++;
			  printf("B");
			  Int_t* sbit;
			  sbit=(pFillData->BluePattern);
			  Int_t blue_offset=pFillData->blue_offset;;
			  char bline[100];
			  Int_t bnch=0;
			  /* two different files formats: for newer runs with fill > 14797 use the unedited cdev dat files and only read the third column of every third row */
			  if (fillNum>14797)
			    {
			      while((!feof(blue)) && bnch<360)
				{
				  if(fgets(bline,50,blue))
				    {
				      if(fmod(bnch,3)==0.0)
					{
					  sscanf(bline, "%*d %*d %d", &(sbit[(Int_t) fmod(blue_offset, 120)]));
					  blue_offset=blue_offset+1;
					}
				      bnch=bnch+1;
				    }
				};
			    }
			  else
			    {
			      while((!feof(blue)) && bnch<120)
				{
				  if(fgets(bline,50,blue))
				    {
				      sscanf(bline,"%d 0 0 %d 0 0 %d 0 0 %d 0 0",
					     &(sbit[(Int_t) fmod(blue_offset,120)]),
					     &(sbit[(Int_t) fmod(blue_offset+1,120)]),
					     &(sbit[(Int_t) fmod(blue_offset+2,120)]),
					     &(sbit[(Int_t) fmod(blue_offset+3,120)])
					     );
				      blue_offset=blue_offset+4;
				      bnch=bnch+4;
				    };
				};
			    }
			  BluePatternFile->close();  
			  pFillData->BlueGood=true;
			}
		      else
			{
			  printf("b\n");
			};
		    };
		  if(!(pFillData->YellowGood))
		    {
		      YellowPatternFile->GetRunFillYellowName(runnum,fillNum);
		      
		      pFillData->YellowPatternFile=YellowPatternFile;
		      
		      FILE* yellow;
		      if(YellowPatternFile->open())
			{
			  yellow=YellowPatternFile->fp;
			  YellowCnt++;
			  printf("Y");
			  Int_t* sbit;
			  sbit=(pFillData->YellowPattern);
			  Int_t yellow_offset=pFillData->yellow_offset;
			  char yline[100];
			  Int_t bnch=0;
			  
			  if (fillNum>14797)
			    {
			      while((!feof(yellow)) && bnch<360)
				{
				  if(fgets(yline,50,yellow))
				    {
				      if(fmod(bnch,3)==0.0)
					{
					  sscanf(yline, "%*d %*d %d", &(sbit[(Int_t) fmod(yellow_offset, 120)]));
					  yellow_offset=yellow_offset+1;
					}
				      bnch=bnch+1;
				    }
				};
			    }
			  else
			    {
			      while((!feof(yellow)) && bnch<120)
				{
				  if(fgets(yline,50,yellow))
				    {
				      sscanf(yline,"%d 0 0 %d 0 0 %d 0 0 %d 0 0",
					     &(sbit[(Int_t) fmod(yellow_offset,120)]),
					     &(sbit[(Int_t) fmod(yellow_offset+1,120)]),
					     &(sbit[(Int_t) fmod(yellow_offset+2,120)]),
					     &(sbit[(Int_t) fmod(yellow_offset+3,120)])
					     );
				      yellow_offset=yellow_offset+4;
				      bnch=bnch+4;
				    };
				};
			    }
			  YellowPatternFile->close();  
			  pFillData->YellowGood=true;
			}
		      else
			{
			  
			  printf("y(%d)\n",YellowPatternFile->Error);
			};
		    };
		};
	    };
	};
      printf("\n \n");
      FillFile->close();
    }
  else
    {
      error=1;
      ErrorMessage=TString(FillFile->Path())+"  file not found";
    };
  if(BlueCnt==0 || YellowCnt==0)
    {
      error=error+2;
      ErrorMessage=ErrorMessage+"\n Missing spinpat files\n";
    };
};
Int_t Fill::ListFills()
{
  if(pFillArray!=0)
    {
      TIter next(pFillArray);
      FillData* fd;
      while(fd=(FillData*) next())
	{
	  printf("Fill Number=%d Run Range (%d , %d)",fd->FillNumber,fd->FirstRun,fd->LastRun);
	  if(fd->BlueGood)printf(" Blue Found  ");
	  if(fd->YellowGood)printf(" Yellow Found ");
	  printf("\n");
	};
    };
  return 0;
};
Float_t* Fill::pBlueEvenOddCor(Float_t dt)
{
  TComplex dot=0;
  TComplex max=0;
  bcorr[0]=-1000.;
  if(dt<=0)return bcorr;
  TComplex w( 2 * TMath::Pi() / dt,0.);
  if(pfd==0)bcorr[0]=-2000.;
  if(!pfd->BlueGood)return bcorr;
  for(Int_t i=0;i<120;i++)
    {
      if(!(pfd->kicked(i)))
	{
	  dot=dot+TComplex::Exp(TComplex::I() * w* (TComplex) i) * (TComplex) pfd->BluePattern[i];
	  max=max+((TComplex) (pfd->BluePattern[i])).Rho();
	};
    };
  bcorr[0]=0.;
  if(max.Rho()==0.)return bcorr;
  dot= dot/max;
  bcorr[0]=dot.Re();
  bcorr[1]=dot.Im();
  bcorr[2]=dot.Rho();
  bcorr[3]=dot.Theta();
  return bcorr;
};

Float_t Fill::BlueEvenOddCor()
{
  return (Float_t) (pBlueEvenOddCor(2.)[0]);
};

Float_t* Fill::pYellowEvenOddCor(Float_t dt)
{
  TComplex dot=0;
  TComplex max=0;
  ycorr[0]=-10000.;
  if(dt<=0)return ycorr;
  TComplex w( 2 * TMath::Pi() / dt,0.);
  if(pfd==0)return ycorr;
  ycorr[0]=-2000;
  if(!pfd->YellowGood)return ycorr;
  for(Int_t i=0;i<120;i++)
    {
      if(!pfd->kicked(i))
	{
	  dot=dot+TComplex::Exp(TComplex::I() * w* (TComplex) i) * (TComplex)( pfd->YellowPattern[i]);
	  max=max+((TComplex) (pfd->YellowPattern[i])).Rho();
	};
    };
  ycorr[0]=0.;
  if( dot.Rho()==0.)return ycorr;
  dot= dot/max;
  ycorr[0]=dot.Re();
  ycorr[1]=dot.Im();
  ycorr[2]=dot.Rho();
  ycorr[3]=dot.Theta();
  return ycorr;
};

Float_t Fill::YellowEvenOddCor()
{
  return (Float_t) (pYellowEvenOddCor(2.)[0]);
};

Int_t Fill::BlueSpin(Int_t bnch7)
{
  if(pfd==0)return 0;
  if(bnch7<0 || bnch7>120)return 0;
  if(pfd->kicked(bnch7))return 0;
  if(pfd->BlueGood)
    return pfd->BluePattern[bnch7];
  return 0;
};

Int_t Fill::YellowSpin(Int_t bnch7)
{
  if(pfd==0)return 0;
  if(bnch7<0 || bnch7>120)return 0;
  if(pfd->kicked(bnch7))return 0;
  if(pfd->YellowGood)return pfd->YellowPattern[bnch7];
  return 0;
};



Int_t Fill::SetFillNumber(Int_t fnum)
{
  FillNumber=0;
  if( (fnum<FirstFill) || (fnum>LastFill))return FillNumber=-1;
  pfd= (FillData*) ((*pFillArray)[fnum]);
  if(pfd==0)return FillNumber=-2;
  return FillNumber=fnum;
};

Float_t* Fill::GetYellow()
{
  if( (FillNumber>=FirstFill) && (FillNumber<=LastFill))
    {
      if((pfd==0 ))
	{
	  ErrorMessage="Attempt to acces undefined Yellow Map";
	  error=-2;
	  return 0;
	}
      else
	{
	  error=0;
	  for(int j=0;j<120;j++)yellowSpin[j]=(Float_t) pfd->YellowPattern[j];
	  return yellowSpin;
	};
    }
  else
    {
      error=-2;
      ErrorMessage="Attempt to acces undefined Yellow Map";
      return 0;
    };
};
Int_t Fill::Print()
{
  if( (FillNumber>=FirstFill) && (FillNumber<=LastFill))
    {
      TString mp[4]={"-"," ","+","k"};
      printf("FillNumber=%d \n",FillNumber);
      if((pfd==0 ))
	{
	  printf("FillData not defined\n");
	  return -1;
	}
      else
	{
	  
	  printf("Run Range (%d , %d ) \n",pfd->FirstRun,pfd->LastRun);
	  Bluestring[0]=  "Blue 0-59    :";
	  Bluestring[1]=  "Blue 60-120  :";
	  Yellowstring[0]="Yellow 0-59  :";
	  Yellowstring[1]="Yellow 60-120:";
	  if(pfd->BluePattern!=0)
	    {
	      for(Int_t j=0;j<60;j++)
		{
		  Int_t pat=pfd->BluePattern[j]+1;
		  if(pfd->kicked(j))pat=3;
		  Bluestring[0]=Bluestring[0]+mp[pat];
		  
		  pat=pfd->BluePattern[j+60]+1;
		  if(pfd->kicked(j+60))pat=3;
		  Bluestring[1]=Bluestring[1]+mp[pat];
		};
	      Bluestring[0]=Bluestring[0]+"\n";
	      Bluestring[1]=Bluestring[1]+"\n";
	    };
	  if(pfd->YellowPattern!=0)
	    {
	      for(Int_t j=0;j<60;j++)
		{
		  Int_t pat=pfd->YellowPattern[j]+1;
		  if(pfd->kicked(j))pat=3;
		  Yellowstring[0]=Yellowstring[0]+mp[pat];
		  
		  pat=pfd->YellowPattern[j+60]+1;
		  if(pfd->kicked(j+60))pat=3;
		  Yellowstring[1]=Yellowstring[1]+mp[pat];
		};
	      Yellowstring[0]=Yellowstring[0]+"\n";
	      Yellowstring[1]=Yellowstring[1]+"\n";
	    };
	  std::cout<<Bluestring[0]<<Bluestring[1]<<Yellowstring[0]<<Yellowstring[1];
	  return 1;
	};
    }
  else
    {
      printf("FillNumber not in defined range (%d to %d)\n",FirstFill,LastFill);
      return -1;
    };
};

Float_t* Fill::GetBlue()
{
  if( (FillNumber>=FirstFill) && (FillNumber<=LastFill))
    {
      if((pfd==0 ))
	{
	  error=-3;
	  ErrorMessage="Attempt to acces undefined Blue  Map";
	  return 0;
	}
      else
	{
	  error=0;
	  for(int j=0;j<120;j++)blueSpin[j]=(Float_t) pfd->BluePattern[j];
	  return blueSpin;
	};
    }
  else
    {
      error=-2;
      ErrorMessage="Attempt to acces undefined Blue Map";
      return 0;
    };
};
Int_t Fill::SetFillNumberforRun(Int_t runnum)
{
  if(CurrentRunNumber==runnum)return runnum;
  Int_t TmpFillNumber;
  FillData* fd;
  if(pfd!=0)
    {
      if( (pfd->FirstRun<=runnum)&&(pfd->LastRun>=runnum))return runnum;
    };
  TIter next(pFillArray);
  
  while( fd=(FillData*) next() )
    {
      if( (fd->FirstRun<=runnum) && (fd->LastRun>=runnum) )
	{
	  CurrentRunNumber=runnum;
	  Int_t filln=fd->FillNumber;
	  return SetFillNumber(filln);
	};
    };
  
  return SetFillNumber(0);
};

TVectorT<float> Fill::BlueVec()
{
  TVectorT<float>  myvec(120,GetBlue());
  return myvec;
};
TVectorT<float> Fill::YellowVec()
{
  TVectorT<float>  myvec(120,GetYellow());
  return myvec;
};

Float_t Fill::Dot(TVectorT<float> v1,TVectorT<float> v2)
{
  Int_t cnt1=v1.GetNoElements();
  Int_t cnt2=v2.GetNoElements();
  if(cnt1!=cnt2)return 0.;
  Float_t dot=0;
  for(Int_t j=0;j<cnt1;j++)
    {
      dot+=v1[j]*v2[j];
      //      printf("v1[j]=%f v2[j]=%f dot=%f \n",v1[j],v2[j],dot);
    };
  return dot;
};
TVectorT<float> Fill::Mask()
{
  TVectorT<float> tmp(120);
  for(int j=0;j<120;j++)tmp[j]=BlueSpin(j)*YellowSpin(j);
  return Mult(tmp,tmp);
}
TVectorT<float> Fill::Mult(TVectorT<float> v1, TVectorT<float> v2)
{
  Int_t cnt1=v1.GetNoElements();
  Int_t cnt2=v2.GetNoElements();
  TVectorT<float> v3(cnt1);
  if(cnt1!=cnt2)
    {
      return v3;
    };
  for(Int_t j=0;j<cnt1;j++)v3[j]=v1[j]*v2[j];
  return v3;
};
