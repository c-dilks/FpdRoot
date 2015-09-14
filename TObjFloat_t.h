#ifndef TObjFloat
#define TObjFloat

#include "TObject.h"
#include "TObjString.h"
class TObjFloat_t : public TObject
{
private:
   Float_t  num;

public:
   TObjFloat_t(Float_t fl = 0.) : num(fl){};
   TObjFloat_t(TObjFloat_t & fl ) { num=fl.GetNum();};
 
   ~TObjFloat_t() {};
   Int_t OriginalIndex;
   void    SetNum(Float_t f) { num = f; }
   Float_t    GetNum() { return num; }
   void    Print(Option_t *) const { Printf("TObjFloat_t = %f", num); }
   Bool_t  IsEqual(const TObject *obj) const { return num == ((TObjFloat_t*)obj)->num; }
   Bool_t  IsSortable() const { return kTRUE; }
   Int_t   Compare(const TObject *obj) const { if (num > ((TObjFloat_t*) obj)->num)
                                      return 1;
                                   else if (num < ((TObjFloat_t*) obj)->num)
                                      return -1;
                                   else
                                      return 0; }
   //   ClassDef(TObjFloat_t,1);
};

#endif
