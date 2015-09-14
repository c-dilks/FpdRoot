#ifndef ProjectGCor_
#include "RunDepMgr.h"

class ProjectGCor : public TObject
{
 public:
  ProjectGCor(TString oldgcorpath,TString RunDepPath,
			    TString newgcorpath="");
  ~ProjectGCor();
  RunDepMgr* Mgr;
  CalibStr* Cal1;
  CalibStr* Cal2;
  CalibStr* ProjectTo(Int_t RunNumber,Int_t OldRunNumber);
  void SetPrint(bool tf=true){kPrint=tf;};
  TH2F* ratio1[4];
  TH2F* ratio2[4];
  TH1F* rat1[4];
  TH1F* rat2[4];
  
 private:
  CalibStr* Cal1Out;
  bool kPrint;
  ClassDef(ProjectGCor,1);
};
#endif

