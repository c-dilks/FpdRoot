#include "Filentry.h"

ClassImp(Filentry)

  Filentry::Filentry(TString directory,TString filename,TString type)
{
  Directory = directory;
  FileName=filename;
  Type=type;
  Description="";
  Error=0;
  isOpen=false;
  exists=false;
};
TString Filentry::GetRunFillBlueName(Int_t rnum,Int_t fnum)
{
   char fname[200];
/* for Run11 and newer data, use the unmodified dat files from cdev, while for older data, use the old modified fill files.*/
   if (fnum>14797)
   {
	sprintf(fname,"blue-buckets-polarization-fillpatterns.fill%d.dat",
	fnum);
   }
   else
   {
   sprintf(fname,
	   "Blue.polarizationFillPattern.fill%d.run%d.txt",
	   fnum,rnum);
   }
   FileName=fname;
   Error=0;
   return fname;
};
TString Filentry::GetRunFillYellowName(Int_t rnum,Int_t fnum)
{
   char fname[200];
   if (fnum>14797)
   {
	sprintf(fname, "yell-buckets-polarization-fillpatterns.fill%d.dat",
	fnum);
   }
   else
   {
   sprintf(fname,
	   "Yell.polarizationFillPattern.fill%d.run%d.txt",
	   fnum,rnum);
   }
   FileName=fname;
   Error=0;
   return fname;
};
Filentry::~Filentry()
{
  close();
};
TString Filentry::GetRunpedFileName(Int_t rnum)
{
  char fname[200];
  Error=0;
  sprintf(fname,"%d.ped.txt",rnum);
  FileName=fname;
  return fname;
};
const char* Filentry::Path()
{
  return path=(const char*) ( Directory+"/"+FileName);
};
Int_t Filentry::Print()
{
  TString itexists=    "     File not Found:=";
  if(check())itexists="     File Found    :=";
  std::cout<<"Type:"<<Description<<Type<<"\n";
  std::cout<<itexists<<path<<"\n";
  return 0;
};
Bool_t Filentry::open()
{
  path= Path();
  if(isOpen)fclose(fp);
  if(fp=fopen((const char*) path,"r"))
    {
      isOpen=true;
      exists=true;
      Error=0;
      return true;
    }
  else
    {
      isOpen=false;
      exists=false;
      Error=errno;
      return false;
    };
  return false;
};

Bool_t Filentry::close()
{
     if(isOpen)
      {
	if(fclose(fp)==0)
	  {
	    Error=0;
	    exists=true;
	    isOpen=false;
	    return true;
	  }
	else
	  {
	    isOpen=false;
	    Error=errno;
	    return false;
	  };
      }
     
     else
       
       {
	 Error=0;
	 return true;
       };
};
Bool_t Filentry::check()
{
  if(isOpen)
    {
      exists=true;
      return true;
    }
  else
    {
      if(open())
	{
	  if(close())
	    {
	      isOpen=false;
	      exists=true;
	    }
	  else
	    {
	      Error=-100;
	    };
	}
      else
	{
	  exists=false;
	  isOpen=false;
	};
    };
  return exists;
};
