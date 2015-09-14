#include "Geom.h"
using namespace std;
ClassImp(Geom)

Geom::Geom(FilesSet* p_files)
{
  InitGeom(p_files);
};
Geom::Geom(TString dir,TString FileName)
{
  FilesSet* p_LocalFilesSet=new FilesSet();
  p_LocalFilesSet->AddFile(8,FileName,"Geometry File",dir);
  InitGeom(p_LocalFilesSet);
  delete p_LocalFilesSet;
  cout<<"ok after InitGeom\n";
};
Geom::~Geom()
{

}
void Geom::InitGeom(FilesSet* p_files)
{
  
  char line[200];
  for(Int_t i=0;i<14;i++)
    {
      dataSize[i]=0;
      datatype[i]="none";
    };
  char star[10]="*";
  datatype[0]="ZFPD";
  dataSize[0]=8;
  datatype[1]="xOffset";
  dataSize[1]=8;
  datatype[2]="yOffset";
  dataSize[2]=8;
  datatype[3]="FpdTowWid";
  dataSize[3]=8;
  datatype[4]="SmdSep";
  dataSize[4]=8;
  datatype[5]="SmdOffset";
  dataSize[5]=8;
  datatype[6]="digitBits";
  dataSize[6]=8;
  datatype[7]="numebin";
  dataSize[7]=1;
  datatype[8]="ebin";
  dataSize[8]=9;
  datatype[9]="numetabin";
  dataSize[9]=1;
  datatype[10]="etabin";
  dataSize[10]=4;
  datatype[11]="absphibin";
  dataSize[11]=1;
  datatype[12]="end";
  dataSize[12]=1;

  Defined=false;
  FMSGeom=false;
  Int_t counter=0;
  Filentry* p_geom=p_files->p_Geom();
  if(!p_geom->check()) return;
  if(p_geom->open())
    {
      FILE* fp=p_geom->fp;
      while(!feof(fp))
	{
	  if(fgets(line,200,fp))
	    {
	      if(line[0]!=star[0] && counter<13)
		{
		  printf("%s :%s",(const char*) datatype[counter],line);
		  Int_t nread= sscanf(line,
				      "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ",
				      &(data[counter][0]),&(data[counter][1]),
				      &(data[counter][2]),&(data[counter][3]),
				      &(data[counter][4]),&(data[counter][5]),
				      &(data[counter][6]),&(data[counter][7]),
				      &(data[counter][8]),&(data[counter][9]),
				      &(data[counter][10]),&(data[counter][11]),
				      &(data[counter][12]),&(data[counter][13]),
				      &(data[counter][14]),&(data[counter][15]));
		  printf("%d words read \n",nread);
		  dataSize[counter]=nread;
		  if(counter==3 && nread==16)FMSGeom=true;
		  counter++;
		};
	      
	    };
	     
	};
      p_geom->close();
    };
};

TVector3 Geom::LocalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Globalxyz,Bool_t BlockUnit)
{
  TVector3 DirX(1.,0,0);
  TVector3 DirY(0,-1.,0);
  TVector3 DirZ(0,0,1);
  if(NSTB16==1 || (EW12==2 && NSTB16==3))DirX=-1*DirX;
  TVector3 voffset(*xOffset(EW12,NSTB16),*yOffset(EW12,NSTB16),
		  *ZFPD(EW12,NSTB16));
  //cout<<"nstb16="<<NSTB16<<" xoff="<<*xOffset(EW12,NSTB16)<<" yoff="<<*yOffset(EW12,NSTB16)
  //<<" Z="<<*ZFPD(EW12,NSTB16)<<" xwidth="<<*FpdTowWid(EW12,NSTB16)<<" ywidth="
  //<<FpdTowWid(EW12,NSTB16)[1]<<endl;
  TVector3 shiftxyz=Globalxyz-voffset;
  TVector3 locxyz=shiftxyz;
  locxyz(0)=shiftxyz.Dot(DirX);
  locxyz(1)=shiftxyz.Dot(DirY);
  locxyz(2)=shiftxyz.Dot(DirZ);
  Float_t scalefactorx=1./(*FpdTowWid(EW12,NSTB16));
  if(FMSGeom)
    {
      if(BlockUnit)
	{
	  Float_t scalefactory=1./(FpdTowWid(EW12,NSTB16)[1]);
	  locxyz(0)*=scalefactorx;
	  locxyz(1)*=scalefactory;
	};
    }
  else
    {
      if(BlockUnit)
	{
	  locxyz*=scalefactorx;
	};
    };
  return locxyz; 
};

TVector3 Geom::GlobalXYZ(Int_t EW12,Int_t NSTB16,TVector3 Localxyz)
{
  TVector3 DirX(1.,0,0);
  TVector3 DirY(0,-1.,0);
  TVector3 DirZ(0,0,1);
  if(NSTB16==1 || (EW12==2 && NSTB16==3))DirX=-1*DirX;
  TVector3 voffset(*xOffset(EW12,NSTB16),*yOffset(EW12,NSTB16),
		  *ZFPD(EW12,NSTB16));
  TVector3 locvec=Localxyz.X()*DirX+Localxyz.Y()*DirY+Localxyz.Z()*DirZ;
  return locvec+voffset;
};

Int_t Geom::getNSTB(Float_t gx, Float_t gy)
{
  Float_t xoff[4],yoff[4];
  for(int i=0; i<4; i++)
    {
      xoff[i] = *xOffset(2,i+1);
      yoff[i] = *yOffset(2,i+1);
    }
  Float_t xywidth12 = *FpdTowWid(2,1);
  Float_t xwidth34 = *FpdTowWid(2,3);
  Float_t ywidth34 = FpdTowWid(2,3)[1];
  if(gx>xoff[1]+xywidth12*17 || gx<xoff[0]-xywidth12*17 || gy>yoff[0] || gy<yoff[0]-xywidth12*34)
    return 0;
  else{
    if(gx<xoff[0] && (gx<xoff[0]-xywidth12*8 || gy>yoff[0]-xywidth12*9 || gy<yoff[0]-xywidth12*25))
      return 1;
    else{
      if(gx>xoff[1] && (gx>xoff[1]+xywidth12*8 || gy>yoff[1]-xywidth12*9 || gy<yoff[1]-xywidth12*25))
	return 2;
      else{
	if(gx>xoff[3]+xwidth34*12 || gx<xoff[2]-xwidth34*12 || gy>yoff[2] || gy<yoff[2]-ywidth34*24)
	  return 0;
	else{
	  if(gx<xoff[2] && (gx<xoff[2]-xwidth34*5 || gy>yoff[2]-ywidth34*5 || gy<yoff[2]-ywidth34*17))
	    return 3;
	  else{
	    if(gx>xoff[3] && (gx>xoff[3]+xwidth34*5 || gy>yoff[3]-ywidth34*5 || gy<yoff[3]-ywidth34*17))
	      return 4;
	    else return 0;
	  }
	}
      }
    }
  }
}

Float_t* Geom::ZFPD(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  return &(data[0][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};

Float_t* Geom::xOffset(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  return &(data[1][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
Float_t* Geom::yOffset(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  
	  return &(data[2][wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
void Geom::Print()
{
  for(Int_t i=0;i<14;i++)
    {
      if(dataSize[i]>0)
	{
	  printf("%s : ",(const char*) datatype[i]);
	  for(Int_t j=0;j<dataSize[i];j++)
	    {
	      printf("%f ",(data[i][j]));
	    };
	  printf("\n");
	};
    };
};
Float_t* Geom::FpdTowWid(Int_t ew, Int_t nstb)
{
  Int_t wcnt;
  Int_t Nxy=1;
  if(FMSGeom)Nxy=2;
  // for FMS geometry file this routine returns an array of Float_t's 
  // first: horizontal width
  // second: vertical width
  if(ew==1 || ew==2)
    {
      if(nstb>0 && nstb<5)
	{
	  wcnt=(ew-1)*4+nstb-1;
	  return &(data[3][Nxy*wcnt]);
	};
    };
  printf("Undefined reference in Geom\n");
  return 0;
};
