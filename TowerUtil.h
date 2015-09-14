#ifndef TOWER_UTIL_H
#define TOWER_UTIL_H

#define NUM_SMD_V 48
#define NUM_SMD_H 48

// North and South
//
#define NS_ETA_NUM 7
#define NS_PHI_NUM 7

// Top and Bottom
//
#define TB_ETA_NUM 5
#define TB_PHI_NUM 5
// total # of tower (North-South, or Top-Bottom)
//

#define MAX_NUMB_PHOTONS 7

// maximum number of clusters that will can be handled
//
#define MAX_NUMER_CLUSTERS 6
#define PRINT_CLUSTER_PARA
#include "TH1.h"
#include "TObjArray.h"
#include "TowerFPD.h"
#include "HitCluster.h"
#include "TPaveText.h"
#include "TH2.h"
#include "TF2.h"
#include "WasExternal.h"



class TowerUtil {

 public:
  TowerUtil(WasExternal* Pwe){TowerUtilCommon(Pwe);};
  TowerUtil(){TowerUtilCommon(0);};
  void  TowerUtilCommon(WasExternal* Pwe);

  ~TowerUtil();
  Int_t nNSTow;
  Int_t nTBTow;
  
  Int_t FindTowerCluster(TObjArray *arrTow, HitCluster *clust);
  
  void CalClusterMoment(HitCluster *clust);
  
  Int_t CatagBySigmXY(HitCluster *clust);
  void PrintTowers(const TowerFPD *tows);
  void VisualizeFit(const Int_t iEv, const Int_t nRow, const Float_t widthLG, const TowerFPD *tows, const Int_t nPh, const Float_t *xPh, const Float_t *yPh, const Float_t *ePh, const Float_t *chi2, TF2 *f2, char path[64]);
  void SetMomentEcutoff(Float_t ecoff=0.5){Ecutoff=ecoff;};
  Float_t GetMomentEcutoff(){return Ecutoff;};
 private:
  ClassDef(TowerUtil,3);
  Float_t Ecutoff;  
  Float_t maxDistanceFromPeak;
  
  Int_t minTowerCatag02;
  Float_t cutEcSigma[2][2];

  TPaveText* pPaveText;
  TObjArray* neighbor;
  TObjArray* arrValley;
  TH2F* h2;

  Float_t minEcSigma2Ph;
  Float_t maxEcSigma1Ph;
  Float_t minTowerEnergy;
  Float_t minRatioPeakTower;
  void  mySetHistTitleLabelProperty(TH1* h);
  WasExternal* pwe;
};

#endif

