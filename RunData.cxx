#include "RunData.h"

ClassImp(RunData)

  RunData::RunData(Int_t runnum,Fill* pFill,FilesSet* p_files)
{
  p_Files=p_files;
  RunNumber=runnum;
  FillNumber=0;
  RunPed=0;
  if(pFill!=0)
    {
      if(pFill->SetFillNumberforRun(runnum)>0)
	{
	  FillNumber=pFill->GetFillNumber();
	};
    };
  char nam[20];
  sprintf(nam,"%d",runnum);
  charname=nam;
  GetRunPed();
  
};

const char* RunData::GetName()
{
  return charname;
};

Int_t RunData::GetRunPed()
{
  if(RunPed!=0)
    {
      printf("Not allowed to read pedistal file more than once\n");
      return 0;
    };
  if(p_Files!=0)
    {
      if(p_Files->p_fpdrunped()==0)
	{
	  printf("p_Files structure missing in RunData \n");
	  RunPed=0;
	  return -10000;
	};
      p_Files->p_fpdrunped()->GetRunpedFileName(RunNumber);
      
      if(p_Files->p_fpdrunped()->Error==0)
	{
	  CalibStr* runped=new CalibStr(RunNumber,(char*) p_Files->p_fpdrunped()->Path() );
	  error=runped->error;
	  if(error>=0)
	    {
	      RunPed= runped;
	      return 0;
	    }
	  else
	    {
	      RunPed=0;
	      return 0;
	    };
	}
      else
	{
	  printf("p_Files Name Building Error in RunData \n");
	  RunPed=0;
	  return 0;
	};
    };
  printf("p_Files structure missing in RunData \n");
  error=-20000;
  return 0;
};
