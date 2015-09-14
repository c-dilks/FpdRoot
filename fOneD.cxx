#include "fOneD.h"
//ClassImp(fOneD)
fOneD::fOneD()
{
  NumberElements=0;
  order=0;
  TA=0;
};
fOneD::fOneD(TVector& tv)
{
  TA=0;
  NumberElements=tv.GetNoElements();
  TA=new TArrayI(NumberElements);
  TSortedList sort(0);
  order=TA->GetArray();
  TList orig;
  for(Int_t n=0;n<NumberElements;n++)
    {
      TObjFloat_t* flo=new  TObjFloat_t(tv[n]);
      flo->OriginalIndex=n;
      orig.Add(flo);
      sort.Add(flo);
    };
  for(int n=0;n<NumberElements;n++)
    {
      TObjFloat_t* pobj=(TObjFloat_t*) (sort.At(n));
      order[n]=pobj->OriginalIndex;
    };
  sort.RemoveAll();
  for(int n=0;n<NumberElements;n++)
    {
      TObjFloat_t* flo=(TObjFloat_t*) orig.At(n);
      if(flo!=0)delete flo;
    };
  orig.RemoveAll();
}
fOneD::~fOneD()
{
  if(TA!=0)delete TA;
};

