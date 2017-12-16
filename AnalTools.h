#ifndef AnalTools_
#define AnalTools_
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "TMath.h"
#include "FillData.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TObjArray.h"
#include "TMatrix.h"
#include "TString.h"
#include "TComplex.h"
#include "FilesSet.h"
#include "TVector.h"
#include "TVectorT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Geom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TFile.h"
#include "dataSet.h"
#include "FpdMap.h"
#include "CalibStr.h"
#include "FitTower.h"
#include "Yiqun.h"
#include "TMath.h"
#include "Trigger.h"
#include "Sim.h"
#include "TKey.h"
#include "TIterator.h"
#include "TMinuit.h"
#include "poutTree.h"
#include "Cell.h"
#include "Vertex.h"
#include "TrigQt.h"
#include "TowerFPD.h"
#include "TGraph2DErrors.h"
#include "RunDepMgr.h"
#include "PullBBC.h"
#include "PullRP.h"
static TGraph2DErrors* p_gra;
static TH2D* p_gra2;
static TH2D* p_doverr;


class AnalTools : public TObject
{
 public:
    AnalTools(); 
    ~AnalTools();
    //    Double_t EdepFactor(Double_t Ephoton){return 1.0/(1+.00562*(fabs(Ephoton)-32.85));};
    Bool_t AllLed(){return allLED;};// set to true for led runs else false
    static void fcn_gra(Int_t &npar,Double_t* grad,Double_t& fval,Double_t*par,Int_t& iflag);
    static Double_t fshape(Double_t* par,Double_t x, Double_t y);
    static Double_t ishape(Double_t* par,Double_t x, Double_t y,Double_t w=1.91);
    void Set_gra(TGraph2DErrors* Gra);
    void Set_gra2(TH2D* Gra2,TH2D* Gdoe);
    TGraph2DErrors* lp_gra;
    Double_t Energy_study;
    Double_t EdepFactor(Double_t Ephoton){return 1.0;};
    Float_t SingleEdepFactor(Int_t vtyp,Float_t energy);
    TH1D* p_histogram;
    Bool_t  SetHistFitEpt(TH2F* hist,Int_t rLim1,Int_t rLim2,Int_t cLim1,Int_t cLim2,TH2F* fhi,Geom* pGlobalGeom );
    Long_t RunList[1000];
    Int_t RunListLength;
    Int_t GetRunIndex(Long_t runnum,Bool_t init=false);
    TH2F** ps_Hist;
    TH1F* p_hist1;
    TH1F* p_hist2;
    Float_t GetHistPairRescale(TH1F* p_hist1,TH1F* p_hist2,Float_t lowfac=.9,Float_t highfac=1.1,Int_t minbin=1);
    TH1F* ReScaledHist;
    TH1F* ReScaleHist(TH1F* p_hist,Float_t xfactor,Int_t minbin=1,Bool_t Normalize=true);
    Geom* p_Geom;
    Int_t NumberLCell4Vectors;
    TLorentzVector* p_LCell4Vec[10];
    Int_t SegmentSearch(TMatrixT<float>* tm,Geom* p_geom,Int_t EW,Int_t NSTB);
    TH1D GetAsm(Int_t, Int_t, TH2F*, char*);
    TH1F GetAsm(TH1F*,TH1F*,char*);
    TH2F GetAsm(TH2F*,TH2F*,char*);
    TH1D GetRatio(Int_t ih,TH2F* h20, char* hame);
    TH2F GetCross(TH2F*, TH2F*,TH2F*, TH2F*, char*,Float_t ermin=1.);
    TH1F GetCross(TH1F*, TH1F*,TH1F*, TH1F*, char*,Float_t ermin=1.);
    TH1F GetCross(Int_t ul,Int_t dl,Int_t ur,Int_t dr ,TH2F* h20,char* name,Int_t secondoff=0);
    Float_t EnergySum;
    TLorentzVector v4;
    TVector3 v3;
    TH1F HistShift(TH1F,Float_t);
    TH2D* CutOnZ(TH3F* ph3,Int_t lowz1,Int_t hiz1,Int_t lowz2=0,Int_t hiz2=0);
    void CopyMatrix(TMatrixT<float>*,TH2F*);
    Int_t PeakZeroDisplace(TH1D*);
    Int_t NTowersAbove(TMatrixT<float>*,Float_t);
    TVector3 mom(TMatrixT<float>*,Int_t,Int_t,Int_t,Int_t,Geom*,Int_t, Int_t);
    TVector3 mom(TMatrixT<float>* ,Geom* ,Int_t, Int_t);
    TLorentzVector FourMom(TMatrixT<float>*,Int_t,Int_t,Int_t,Int_t,Geom*,Int_t, Int_t);
    Int_t readO(char*,Int_t , Int_t , Int_t , Int_t );
    Int_t read1(char* filelist,Int_t set=61);
    Int_t read1a(char* filelist,Int_t set,FilesSet* p_files); 
    Int_t readSimple(char* filelist,Int_t set=61);
    Int_t readq1(char* filelist,Int_t set,FilesSet* p_files);    
    Int_t readqM(char* filelist,Int_t set,FilesSet* p_files);    

    Int_t Examp_script4();
    Int_t Select_script1(Int_t set=80);
    Int_t Select_scriptQT1(Int_t set=80,FilesSet* P_Files=0,char* hntName=(char*)"h111_08");
    Int_t Select_script_RPtest(Int_t set=80,FilesSet* P_Files=0,char* hntName=(char*)"h111_08");
    Int_t Select_Trig(char* infile=(char*)"./run*.root",Int_t set=80,FilesSet* P_Files=0,char* hntName=(char*)"h111_08");
    Int_t Select_ADC_QT(Int_t set=80,FilesSet* P_Files=0,char* hntName=(char*)"h111_08");
    Int_t Select_ADC_QT2(Int_t set=80,FilesSet* P_Files=0,char* hntName=(char*)"h111_08");
    Int_t NumberEventsToProcess;
    Int_t FirstEventToProcess;
    Bool_t OutputToSingle;
    TString OutputFileName;
    poutTree* pout;
    Bool_t SimScript(char* path);
    Bool_t SimScript2(char* path);
    TString OutFileName;
    Sim* pSim;
    void (*FcnFitEPt)(Int_t & npar, Double_t *grad,  Double_t &f, Double_t *par, Int_t iflag);
    Float_t OutOf10;
    Double_t (*FunEPt)(Double_t *x, Double_t *par);
    //   TLorentzVector FourMom(TMatrixT<float>*,Geom*,Int_t, Int_t);
    Int_t* p_colMask;
    Bool_t SetFitOption(Int_t iopt=0);
    Int_t GetFitOption();
    Bool_t WriteADC;
    TH2F* Mat2NewHist2(TMatrixT<float>* tm,char* hname=(char*)"temph");
    TMatrixT<float> MultMat(TMatrixT<float>* tm1,TMatrixT<float>* tm2);
    void chkEta(char* filelist=(char*)"tmp.txt",float eta0=-3.65,float dcenter=.15,
		float e0=85.,float de=5.,float m0=.55,float dm=.08);
    void TrimFMS(CalibStr* FF);
    void SetCorrLimitFMS(CalibStr* fpdc,Float_t max=35.,Float_t min=0.);
    Bool_t IsPopulatedFMS(Cell* cl);
    Bool_t IsPhysicalFMS(Cell* cl);
    void CheckCluster(Yiqun*);
    unsigned int SaveClHit(TowerFPD* h,Int_t IEW,Int_t INSTB,Int_t row,Int_t col);
    Int_t MaxSavedHits;
    void PrClHit(unsigned int);
    unsigned int SavedHits[1000];// List of all Saved Hits
    Int_t  SavedCluHitIndex[50];//  Index into SavedHits for each Cluster [nCluster]
    Int_t SavedPhotonCluIndex[50];// Cluster Index for each Photon [nPhotonClu]
    Float_t SavedPhotonCluEnergy[50];//Energy for each Photon
    Int_t nCluster;//Number of Clusters
    Int_t nPhotonClu;//Number of photons
    Int_t nSavedHits;
    Int_t TrigBits;
    Int_t TrigBitsSelect;
    Int_t TrigBitsVeto;
    void storeCluster(Yiqun* rec,unsigned int minADC=1,Float_t minclusterE=5.);
    void storerec(Yiqun* rec,unsigned int minADC=1);
    Int_t ClNearX(Yiqun*,Float_t, Float_t);
    Int_t PhNearX(Yiqun*,Float_t,Float_t);
    Int_t readq2(char* filelist,Int_t set,FilesSet* p_files);
    Bool_t IdentifyLED;
    Bool_t OnlyLED;
    Bool_t allLED;
    TString File_adcTr;
    Float_t ERecMin[4];//low limits for first Recon
    dataSet* currentdset;
    RunDepCor* Rdep;
    RunDepCor* BaseRdep;
    Int_t RunDepBaseOverride; //if non-zero then this run to be used as  base
    Bool_t UseRunDepCor;
    Bool_t AddRunToRunDepCorFile(Int_t runNumber,char* adcTrLpath,Int_t baseRunNumber,char* scratchDir);
    void RCorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW,Int_t runnum,Int_t EvNum,Int_t SegNM);
    void RCorrect(TMatrix* padc,TMatrix* pEmat,Int_t det,Int_t iEW=2);
    RunDepMgr* Mgr;
    Bool_t Enable_ReRec;
    Int_t OutputTwoTr(char* filelist, FilesSet* p_files);
    Int_t OutputMycells(char* filelist, FilesSet* p_files, Bool_t doSingle=false);
    RunDepMgr* OpenRunDepMgr(TString);
    void PrintCell(Cell* ,Int_t,Int_t ,TString);
    void Extract1(TObjArray*,TString cellname, TString dirname);
    bool TestCellTree(Cell*);
    UInt_t L2sum[2];
    UInt_t lastdsm[8];
    UInt_t Fpde[8];
    float Eminpout[4];
    PullBBC*  pbbc;
    PullRP*  prp;
    QTBBCInfo qtbinfo;
    QTRPInfo qtrpinfo;


 private:
    ClassDef(AnalTools,5); 
};

#endif
