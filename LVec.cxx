#include "LVec.h"
ClassImp(LVec)

  LVec::LVec(Int_t iew,Int_t nstb,TLorentzVector v):TLorentzVector(v)
{
  Iew=iew;
  Nstb=nstb;
  Partner=0;
  PhotOrder=-1;
  NewPartner=TLorentzVector(0,0,0,0);
  Vtype=(iew-1)*4+nstb;
  if(Vtype>9)Vtype=0;
};

LVec::LVec(Int_t vtyp,TLorentzVector v):TLorentzVector(v)
{
  Vtype=vtyp;
  if(vtyp>99) 
    {
      Vtype=vtyp-100; // here Vtype simply records order for after sorting
      Iew=Nstb=0;
      Partner=0;
    }
  else
    {
      Iew=vtyp/5+1;
      Nstb=vtyp-(Iew-1)*4;
      Partner=0;
    }
  PhotOrder=-1;
  NewPartner=TLorentzVector(0,0,0,0);
};

Float_t LVec::PairMass()
{
  if(Partner)
    {
      return (*this + NewPartner).Mag();
    }
  else return this->Mag();
};
Float_t LVec::PairMass(LVec* pv2)
{
  if(pv2!=0){
    return (*this+*pv2).Mag();
  };
  return 0;
};
void LVec::SetPartner(LVec* newpartner)
{
  Partner=newpartner;
  NewPartner=*newpartner;
  newpartner->Partner=this;
};

