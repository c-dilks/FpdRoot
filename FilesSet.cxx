#include "FilesSet.h"

ClassImp(FilesSet)
FilesSet::FilesSet()
{
  FileList=new TObjArray(12,0); 
};
FilesSet::FilesSet(TString defaultDir,TString s_fpdped,TString s_fpdgain,
		   TString s_fpdgaincorr,TString s_Fill,TString s_lumi,TString s_spinpatdir,
		   TString s_Geom,TString s_qtmap,TString s_qtmap2pp)
{
  FileList=new TObjArray(12,0);
  Filentry* p_fpdped=new Filentry((const char*) defaultDir,s_fpdped,"Run Independent Ped File");
  FileList->AddAt(p_fpdped,0);
  Filentry* p_fpdgain=new Filentry((const char*) defaultDir,s_fpdgain,"Run Independent Gain File");
  FileList->AddAt(p_fpdgain,1);
  Filentry*  p_fpdgaincorr=new Filentry((const char*) defaultDir,s_fpdgaincorr,"Run Independent Gain File");
  FileList->AddAt(p_fpdgaincorr,2);
  TString peddir=defaultDir+"/"+"runped";
  Filentry* p_fpdrunped=new Filentry((const char*) peddir,"PedDependsOnRun","Run Dependent Ped File");
  FileList->AddAt(p_fpdrunped,4);
  Filentry* p_FillFile=new Filentry((const char*) defaultDir,s_Fill,"List of Fills and Runs");
  FileList->AddAt(p_FillFile,5);
  s_spinpatdir=defaultDir+"/"+s_spinpatdir;  
  Filentry* p_BluePattern=new Filentry((const char*) s_spinpatdir,"Blue_DependsOnFill","Blue Spin Pattern File");
  FileList->AddAt(p_BluePattern,6);
  Filentry* p_YellowPattern=new Filentry((const char*) s_spinpatdir,"Yellow_DependsOnFill","Yellow Spin Pattern File");
  FileList->AddAt(p_YellowPattern,7);
  TString lumidir=defaultDir;
  if(s_lumi=="Fake")
    {
      UseFakeLumi();
      lumidir="./";
      s_lumi="FakeLumi.txt";
    };
  
  Filentry* p_fpdlumi=new Filentry((const char*) lumidir,s_lumi,"Run Independent Lumi File");
  FileList->AddAt(p_fpdlumi,3);
  Filentry* p_Geom=new Filentry((const char*) defaultDir,s_Geom,"Geometry File");
  FileList->AddAt(p_Geom,8);   
  if(s_qtmap!="none")
    {
      Filentry* p_qtmap=new Filentry((const char*) defaultDir,s_qtmap,"QtMap File");
      FileList->AddAt(p_qtmap,9);   
    };
  if(s_qtmap2pp!="none")
    {
      Filentry* p_qtmap2pp=new Filentry((const char*) defaultDir,s_qtmap2pp,"QtMap File");
      FileList->AddAt(p_qtmap2pp,10);   
    };
};
Int_t FilesSet::Print()
{
  TIter next(FileList);
  while(Filentry* afile=(Filentry*) next())
    {
      afile->Print();
    };
  return 0;
};
Filentry* FilesSet::AddFile(Int_t slot,TString s_filename,TString Description,TString defaultDir)
{
  if(slot<0 || slot>11)return 0;
  if(FileList->At(slot))
    {
      std::cout<<"Removing from slot = "<<slot<<"\n";
      FileList->RemoveAt(slot);
    };

  Filentry*  p_tmpEntry=new Filentry((const char*) defaultDir,s_filename,Description);
  FileList->AddAt(p_tmpEntry,slot);

  return p_tmpEntry;
};
Int_t FilesSet::UseFakeLumi()
{
  if(p_FillFile()!=0)
    {
      printf("Creating FakeLumi.txt from entries in %s \n",p_FillFile()->Path());
      FILE* outfp;
      if(!(outfp=fopen("FakeLumi.txt","w")))
	{
	  printf("failed to create FakeLumi.txt\n");
	  return -10;
	};
      if(p_FillFile()->open())
	{
	  char line[300];
	  char outline[200];
	  FILE* fp=p_FillFile()->fp;
	  while(!feof(fp))
	    {
	        if(fgets(line,300,fp))
		  {
		    Int_t runnum;
		    sscanf(line,"%d",&runnum);
		    sprintf(outline,"%d  0  11  1.0000 0.0005  1.0000 0.0005\n",runnum);
		    fputs(outline,outfp);
		  };
	    };
	  p_FillFile()->close();
	  fclose(outfp);
	  printf("hello\n");
	  return 0;
	}
      else
	{
	  return -1;
	};
    };
};
FilesSet::~FilesSet()
{
  FileList->Delete();
  delete FileList;
};

