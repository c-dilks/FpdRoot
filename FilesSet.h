#ifndef FilesSet_
#define FilesSet_
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <iostream>
#include "TObject.h"
#include "TString.h"
#include "TObjArray.h"
#include "Filentry.h"

class FilesSet : public TObject
{
 public: 
  FilesSet(TString defaultDir,TString s_fpdped,TString s_fpdgain,TString s_fpdgaincorr,TString s_Fill,TString s_lumi,TString s_spinpatdir="spinpat",TString s_Geom="geom.txt",
	   TString s_qtmap="none",TString s_qtmap2pp="none");

  FilesSet();
  ~FilesSet();
  Filentry* AddFile(Int_t slot,TString s_filename,TString Description="none",
		TString defaultDir="");
  TObjArray* FileList;
  Filentry* p_fpdped(){return (Filentry*) FileList->At(0);};
  Filentry* p_fpdgain(){return (Filentry*) FileList->At(1);};
  Filentry* p_fpdgaincorr(){return (Filentry*) FileList->At(2);};
  Filentry* p_fpdlumi(){return (Filentry*) FileList->At(3);};
  Filentry* p_fpdrunped(){return (Filentry*) FileList->At(4);};
  Filentry* p_BluePattern(){return (Filentry*) FileList->At(6);};
  Filentry* p_YellowPattern(){return (Filentry*) FileList->At(7);};
  Filentry* p_FillFile(){return (Filentry*) FileList->At(5);};
  Filentry* p_Geom(){return (Filentry*) FileList->At(8);};
  Filentry* p_qtmap(){return (Filentry*) FileList->At(9);};
  Filentry* p_qtmap2pp(){return (Filentry*) FileList->At(10);};

  Int_t Print();
  Int_t UseFakeLumi();
 private:
  ClassDef(FilesSet,1);
};

#endif
