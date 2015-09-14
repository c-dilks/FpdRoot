#ifndef Sim_h
#define Sim_h
#include <TString.h>
#include <TMatrix.h>
#include <TRandom1.h>
#include "simh112.h"
#include "Geom.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TH3F.h>

const Float_t Esigped[6]={0.9, 0.9, 0.3, 0.3, 0.6, 0.6};
const Float_t Wsigped[6]={0.3, 0.3, 0.3, 0.3, 0.0, 0.0};

class Sim: public simh112
{
 public:
  Sim(char* files,Geom* pGeom);
  ~Sim();
  TString SimType;
  Int_t isSim;

  TMatrix* ene[6][2];
  TMatrix* newene[6][2];
  TMatrix* sumene[6][2];
  TMatrix* noise[6][2];
  TMatrix* adc[6][2];

  TRandom1* p_ran;
  Bool_t AddNoise();
  Bool_t FillNoise();
  Float_t sigped[2][6];
  Bool_t FillHits(Float_t Mmin=-100.,Float_t Mmax=100);
  Bool_t PrintHits(Float_t hitmin);
  TVector3 ProjectTrack(Int_t tracknum,Float_t zdet);
  Geom* pGeom;
  Bool_t GenADC(Float_t Mmin=-100.,Float_t Mmax=100);
  Float_t digitBits(Int_t nstb,Int_t ew);
  Int_t HighFour[4]; //pointers starting at zero
  Bool_t CheckForHigh(Int_t id_chk,Float_t tstmass=0,Float_t dmass=.00001);
  Bool_t FillHighFour(Float_t tstmass=0,Float_t dmass=.00001);
  Bool_t HitTest(Int_t hitid);
  Int_t SelectTest;
  TH3F* shape;
  Bool_t Fillshape();
  Int_t skipID;

 private:
  Float_t GAIN;
  Bool_t getvolid(Int_t vid_in);
  Int_t vid;
  Int_t towsmd;
  Int_t vh;
  //following variables with min values of 1
  Int_t ew;
  Int_t nstb;
  Int_t ch;
  ClassDef(Sim,1);
};

#endif

