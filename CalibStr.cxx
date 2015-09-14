#include "CalibStr.h"
ClassImp(CalibStr)
  
  CalibStr::CalibStr(const char* FileName)
{
  CalibStrBuild(0,FileName);
};
CalibStr::CalibStr(Int_t run,const char* FileName)
{
  CalibStrBuild(run,FileName);
};
CalibStr::CalibStr()
{
  M_array=new TObjArray(12);
  M_array->SetOwner();
};
Int_t CalibStr::CalibStrBuild(Int_t run,const  char* FileName)
{
  error=0;
  ZERO=0.;
  for(Int_t i=0;i<2;i++)
    {
      MatrixDimensionRows[i][0]=7;
      MatrixDimensionRows[i+2][0]=5;
      MatrixDimensionRows[i+4][0]=1;
      MatrixDimensionRows[i][1]=7;
      MatrixDimensionRows[i+2][1]=5;
      MatrixDimensionRows[i+4][1]=1;

      MatrixDimensionCols[i][0]=7;
      MatrixDimensionCols[i+2][0]=5;
      MatrixDimensionCols[i+4][0]=7;
      MatrixDimensionCols[i][1]=7;
      MatrixDimensionCols[i+2][1]=5;
      MatrixDimensionCols[i+4][1]=7;
    };
  if(run>=7094066)
    {
      printf(" Run>= 7094066\n");
      MatrixDimensionRows[0][1]=14;
      MatrixDimensionRows[1][1]=14;
      MatrixDimensionRows[2][1]=6;
      MatrixDimensionRows[3][1]=6;
      MatrixDimensionRows[4][1]=0;
      MatrixDimensionRows[5][1]=0;

      MatrixDimensionCols[0][1]=14;
      MatrixDimensionCols[1][1]=14;
      MatrixDimensionCols[2][1]=6;
      MatrixDimensionCols[3][1]=6;
      MatrixDimensionCols[4][1]=0;
      MatrixDimensionCols[5][1]=0;
    };
  if(run>=8000000)
    {
      printf(" Run>= 7094066\n");
      MatrixDimensionRows[0][1]=34;
      MatrixDimensionRows[1][1]=34;
      MatrixDimensionRows[2][1]=24;
      MatrixDimensionRows[3][1]=24;

      MatrixDimensionCols[0][1]=17;
      MatrixDimensionCols[1][1]=17;
      MatrixDimensionCols[2][1]=12;
      MatrixDimensionCols[3][1]=12;
    };

  filename=TString(FileName);
  OutName=filename+TString("_V0");
  std::cout<<"CalibStr with: run#="<<run<<" :"<<filename<<"\n";
  M_array=new TObjArray(12);
  M_array->SetOwner();
  M_array->ls();
  for(Int_t ew=0;ew<2;ew++)
    {
      for(Int_t nsup=0;nsup<6;nsup++)
	{
	  Int_t md2;
	  Int_t mdr=MatrixDimensionRows[nsup][ew];
	  Int_t mdc=MatrixDimensionCols[nsup][ew];
	  Int_t count=ew*6+nsup;
	  md2=mdr*mdc;
	  Float_t vv[md2];
	  for(Int_t kz=0;kz<md2;kz++){vv[kz]=.001*nsup*kz;};
	  TMatrix* tm0;
	  Int_t kcnt=0;
	  M_array->AddAtAndExpand(tm0=new TMatrix(mdr,mdc,vv),count);
	  *tm0=0;
	};
     
    };
  if(GetCalibration()!=1)error=1;
}
Int_t CalibStr::SetValue(Int_t IEW,Int_t INSTB, Int_t Irow0, Int_t Icol0,Float_t val)
{
  Int_t dim=MatrixDimensionCols[INSTB-1][IEW-1];
  if(dim<Icol0)return 0;
  Int_t ichan=Irow0*dim+Icol0+1;
  return SetValue(IEW,INSTB,ichan,val);
}

Int_t CalibStr::SetValue(Int_t EW,Int_t NSUD,Int_t chan,Float_t newvalue)
{
// remember that in callib file, indices start at iew=1 and instb=1
  if(EW==1 || EW==2)
  {
    if(NSUD>0&& NSUD<7)
      { 
	Int_t mdr=MatrixDimensionRows[NSUD-1][EW-1];
	Int_t mdc=MatrixDimensionCols[NSUD-1][EW-1];
	Int_t md=0;
	if(NSUD<5)md=mdr*mdc;
	if(chan>0 && chan<md+1)
	  {
	    *(value(chan-1,NSUD-1,EW-1))=newvalue;
	    return 1;
	  };
      };
  };
  return 0;
};
Float_t CalibStr::GetValue(Int_t IEW,Int_t INSTB, Int_t Irow0, Int_t Icol0)
{
  //first Icol0 =0; first Irow0=0;

  return (*((TMatrix*) ( M_array->At((IEW-1)*6+INSTB-1))))(Irow0,Icol0);
};
Int_t CalibStr::UpdateFile(Int_t version=0)
{
  char tmp[100];
  sprintf(tmp,"_V%d",version);
  OutName=filename+TString(tmp);
  FILE* fp;
  char line[110];
  int len=100;
  char fname[200];
  strcpy(fname,(const char*) OutName);
  if(fp=fopen(fname,"w"))
    {
      for(Int_t ew=0;ew<2;ew++)
	{
	  for(Int_t nstb=0;nstb<6;nstb++)
	    {
	      Int_t mdr=MatrixDimensionRows[nstb][ew];
	      Int_t mdc=MatrixDimensionCols[nstb][ew];
	      Int_t md=mdr*mdc;
	      for(Int_t ch=0;ch<md;ch++)
		{
		  fprintf(fp,"%d %d %d %5.3f \n",ew+1,nstb+1,ch+1,*(value(ch,nstb,ew)));
		};
	    };
	};
      fclose(fp);
      return (Int_t)  1;
    };
  return 0;
};

Int_t CalibStr::GetCalibration()
{
    error=0;
    FILE* fp;
    char line[110];
    int len=100;
    char fname[200];
    strcpy(fname,(const char*) filename);
    std::cout<<"from GetCalibration="<<fname<<"\n";
    if(fp=fopen(fname,"r"))
      {
	while(!feof(fp))
	  {
	    int i=0;
	    int j=0;
	    int k=0;
	    int i1=0;
	    int i2=0;
	    Int_t count=0;
	    float val;
	    char comment[10]="# comment";
	    char* thechars=fgets(line,len,fp);
	    if(line[0]!=comment[0])
	      {	    
		sscanf(line,"%d %d %d %f",&i,&j,&k,&val);
		if(i>0 && j>0 && k>0&& i<3&& j<7&& k<579)
		  {
		    Int_t mdr=MatrixDimensionRows[j-1][i-1];
		    Int_t mdc=MatrixDimensionCols[j-1][i-1];
		    Int_t  md2=mdr*mdc;
		    if(k<=md2&& mdc!=0)
		      {
			i1=(k-1)/mdc;
			i2=(k-1)%mdc;
			if(i1>=0 && i1<mdr && i2>=0 && i2<mdc)
			  { 
			    (*tm(i,j))(i1,i2)=val;
			  };
		      }
		    else
		      {
			//			printf("bad: %d %d %d %f %d %d count=%d\n",i,j,k,val,i1,i2,count);

			//			error=-10;//file format error
			//			return 0;
		      };
		  };
	      };
	  };
	fclose(fp); 
	error=0;
	return 1;
      }
    else
      {
	error=-1; // file not found
	printf("fail to open calibstr file\n");
	return 0;
      };
};

    Bool_t CalibStr::Compare(Int_t IEW,Int_t INSTB,CalibStr* other)
{
  Int_t nr= tm(IEW,INSTB)->GetNrows();
  Int_t nc= tm(IEW,INSTB)->GetNcols();
  if(other->tm(IEW,INSTB)->GetNrows()!=nr ||
     other->tm(IEW,INSTB)->GetNcols()!=nc)
    {
      printf("I cannot compare two matrices of different size \n");
      return false;
    };
  Bool_t same=true;
  Int_t cnt=0;
  for(Int_t ic=0;ic<nc;ic++)
    {
      for(Int_t ir=0;ir<nr;ir++)
	{
	  cnt=ic+nc*ir+1;
	  Float_t val1=GetValue(IEW,INSTB,ir,ic);
	  Float_t val2=other->GetValue(IEW,INSTB,ir,ic);
	  Float_t ratio=.0;
	  if(val2>=0)ratio=val1/val2;
	  printf("n=%d col=%d row=%d this=%f that = %f (ratio=%f)\n",cnt,ic+1,ir+1,val1,val2,ratio);
	  if(val1!=val2)same=false;
	};
    };
  return same;
}
Float_t* CalibStr::value(Int_t index, Int_t i_ns,Int_t i_ew)
{
  TMatrix* tm0=((TMatrix*) ( M_array->At((i_ew)*6+i_ns)));
  Int_t ncol=tm0->GetNcols();
  if(ncol<0)return &ZERO;
  Int_t nrows=tm0->GetNrows();
  Int_t Irow0=index/ncol;
  Int_t Icol0=index%ncol;
  if(Irow0>=nrows)return &ZERO;
  return &((*tm0)(Irow0,Icol0));
};
