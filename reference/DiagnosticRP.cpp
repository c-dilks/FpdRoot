// MUST COMPILE: execute as root -b -q DiagnosticRP.C++
//
// draws various plots for RPs
// -- ewb012 is a control number for looking at:
//    east RPs only if ewb012==0
//    west RPs only if ewb012==1
//    both RPs if ewb012==2
//
// -- branchMask
//    simply mask for branches using 
//
//    (1 << t_branch) & branchMask
//
//     3  2  1  0
//    WD WU ED EU
// 1   0  0  0  1 --> EU
// 2   0  0  1  0 --> ED
// 3   0  0  1  1 --> E
// 4   0  1  0  0 --> WU
// 5   0  1  0  1
// 6   0  1  1  0
// 7   0  1  1  1
// 8   1  0  0  0 --> WD
// 9   1  0  0  1
// 10  1  0  1  0
// 11  1  0  1  1
// 12  1  1  0  0 --> W
// 13  1  1  0  1
// 14  1  1  1  0
// 15  1  1  1  1 --> E+W

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
//#include "./StRoot/StSpinPool/StFmsAnalysisMaker/StFmsAnalysisMaker.h"
//#include "./cvs/StMuRpsTrack.h"
//#include "./cvs/StMuRpsTrackPoint.h"


// hash a pair of cluster IDs into a number representing their difference
Int_t HashCluIds(Int_t num0A, Int_t num1A, Int_t num0B, Int_t num1B) {
  Int_t numA = 10000 * num0A + num1A;
  Int_t numB = 10000 * num0B + num1B;
  Int_t digitA,digitB,digit;
  Int_t count = 0;
  Int_t retval = 0;
  for(Int_t di=0; di<8; di++) {
    digitA = numA / (Int_t)(TMath::Power(10,di)) % 10;
    digitB = numB / (Int_t)(TMath::Power(10,di)) % 10;
    digit = TMath::Abs(digitA-digitB);
    if(digit>0) {
      retval += (Int_t)(TMath::Power(10,count)) * digit;
      count++;
    };
  };
  return retval;
};


//void DiagnosticRP(TString filename="st_fms_16088020_raw_5500016.single.fmsan.root",
//void DiagnosticRP(TString filename="hist_pptrans/OFset08705.root",
void DiagnosticRP(TString filename="hist_pptrans/OFset08802.root",
                  Int_t ewb012=1
                 ) {
  printf("analysing %s with ewb012 = %d\n",filename.Data(),ewb012);

  // track enums, from StMuRpsTrack.h
	enum StMuRpsTrackType { rpsLocal, rpsGlobal, rpsUndefined };
  enum StMuRpsAngles { rpsAngleThetaX, rpsAngleThetaY, rpsAngleTheta, mNumberOfAngleTypes };
  enum {mNumberOfStationsInBranch = 2};

  // control variables
  Bool_t use_diag_dir = true; // if true, output is directed to *_diag directory
  Bool_t print_evids = false; // if true, prints out evids of tree (useful for comparing trees)
  Bool_t printRPevents=false; // if true, prints some event-by-event info

  TFile * infile = new TFile(filename.Data(),"READ");
  TTree * p_out = (TTree*) infile->Get("p_out");

  //gROOT->Macro("load.C");  // Load all required libraries
  //gSystem->Load("StFmsAnalysisMaker");
  //Int_t n_tracks_max_tmp = StFmsAnalysisMaker::NTracksMax();
  //const Int_t n_tracks_max = n_tracks_max_tmp;
  const Int_t n_tracks_max = 1000;

  TString BranchName[5] = {"all","EU","ED","WU","WD"};

  Int_t evid,n_tracks,n_trackpoints;

  // track vars
  Int_t t_index[n_tracks_max];
  Int_t t_branch[n_tracks_max];
  Int_t t_type[n_tracks_max];
  UInt_t t_planesUsed[n_tracks_max];
  Double_t t_p[n_tracks_max];
  Double_t t_pt[n_tracks_max];
  Double_t t_eta[n_tracks_max];
  Double_t t_time[n_tracks_max];
  Double_t t_theta[n_tracks_max][mNumberOfAngleTypes];
  Double_t t_thetaRP[n_tracks_max][mNumberOfAngleTypes];
  Double_t t_phi[n_tracks_max];
  Double_t t_phiRP[n_tracks_max];
  Double_t t_t[n_tracks_max];
  Double_t t_xi[n_tracks_max];
  Bool_t t_gold[n_tracks_max];
  Bool_t t_isBad[n_tracks_max];
  Double_t t_qualHash[n_tracks_max];

  // trackpoint vars
  Bool_t p_tpExists[2][n_tracks_max];
  Int_t p_RPid[2][n_tracks_max];
  Int_t p_clustid_s1[2][n_tracks_max]; 
  Int_t p_clustid_s2[2][n_tracks_max]; 
  Int_t p_clustid_s3[2][n_tracks_max]; 
  Int_t p_clustid_s4[2][n_tracks_max]; 
  Int_t p_quality[2][n_tracks_max];
  UInt_t p_planesUsed[2][n_tracks_max];
  Double_t p_x[2][n_tracks_max];
  Double_t p_y[2][n_tracks_max];
  Double_t p_z[2][n_tracks_max];

  ////////
  p_out->SetBranchAddress("br_ievt",&evid);
  p_out->SetBranchAddress("RP_n_tracks",&n_tracks);
  p_out->SetBranchAddress("RP_n_trackpoints",&n_trackpoints);

  // tracks
  p_out->SetBranchAddress("RP_t_index",t_index);
  p_out->SetBranchAddress("RP_t_branch",t_branch);
  p_out->SetBranchAddress("RP_t_type",t_type);
  /* 0=rpsLocal -- 1 track point
   * 1=rpsGlobal -- 2 track points
   * 2=rpsUndefined -- track not defined
   */
  p_out->SetBranchAddress("RP_t_planesUsed",t_planesUsed);
  p_out->SetBranchAddress("RP_t_p",t_p);
  p_out->SetBranchAddress("RP_t_pt",t_pt);
  p_out->SetBranchAddress("RP_t_eta",t_eta);
  p_out->SetBranchAddress("RP_t_time",t_time);
  p_out->SetBranchAddress("RP_t_theta",t_theta);
  p_out->SetBranchAddress("RP_t_thetaRP",t_thetaRP);
  p_out->SetBranchAddress("RP_t_phi",t_phi);
  p_out->SetBranchAddress("RP_t_phiRP",t_phiRP);
  p_out->SetBranchAddress("RP_t_t",t_t);
  p_out->SetBranchAddress("RP_t_xi",t_xi);
  p_out->SetBranchAddress("RP_t_gold",t_gold);
  p_out->SetBranchAddress("RP_t_isBad",t_isBad); 
  p_out->SetBranchAddress("RP_t_qualHash",t_qualHash); 

  // track point 0
  p_out->SetBranchAddress("RP_p0_tpExists",p_tpExists[0]);
  p_out->SetBranchAddress("RP_p0_RPid",p_RPid[0]);
  /* 0=E1U  1=E1D  2=E2U  3=E2D
   * 4=W1U  5=W1D  6=W2U  7=W2D
   */
  p_out->SetBranchAddress("RP_p0_clustid_s1",p_clustid_s1[0]);
  p_out->SetBranchAddress("RP_p0_clustid_s2",p_clustid_s2[0]);
  p_out->SetBranchAddress("RP_p0_clustid_s3",p_clustid_s3[0]);
  p_out->SetBranchAddress("RP_p0_clustid_s4",p_clustid_s4[0]);
  p_out->SetBranchAddress("RP_p0_quality",p_quality[0]);
  /* 0=rpsNormal -- not golden and not undefined
   * 1=rpsGolden -- single cluster in all 4 SSD planes
   * 2=rpsNotSet -- undefined track point
   */
  p_out->SetBranchAddress("RP_p0_planesUsed",p_planesUsed[0]);
  p_out->SetBranchAddress("RP_p0_x",p_x[0]);
  p_out->SetBranchAddress("RP_p0_y",p_y[0]);
  p_out->SetBranchAddress("RP_p0_z",p_z[0]);

  // track point 1
  p_out->SetBranchAddress("RP_p1_tpExists",p_tpExists[1]);
  p_out->SetBranchAddress("RP_p1_RPid",p_RPid[1]);
  p_out->SetBranchAddress("RP_p1_clustid_s1",p_clustid_s1[1]);
  p_out->SetBranchAddress("RP_p1_clustid_s2",p_clustid_s2[1]);
  p_out->SetBranchAddress("RP_p1_clustid_s3",p_clustid_s3[1]);
  p_out->SetBranchAddress("RP_p1_clustid_s4",p_clustid_s4[1]);
  p_out->SetBranchAddress("RP_p1_quality",p_quality[1]);
  p_out->SetBranchAddress("RP_p1_planesUsed",p_planesUsed[1]);
  p_out->SetBranchAddress("RP_p1_x",p_x[1]);
  p_out->SetBranchAddress("RP_p1_y",p_y[1]);
  p_out->SetBranchAddress("RP_p1_z",p_z[1]);


  ///////////////////////////////////////////////////////
 
  // track hists
  TH1D * NTracksPerEvent = new TH1D("NTracksPerEvent",
    "total num. tracks per event",
    n_tracks_max,0,n_tracks_max);
  TH1D * NGlobalTracksPerEvent = new TH1D("NGlobalTracksPerEvent",
    "total num. global tracks per event",
    n_tracks_max,0,n_tracks_max);
  TH1D * NGoldTracksPerEvent = new TH1D("NGoldTracksPerEvent",
    "total num. gold tracks per event",
    n_tracks_max,0,n_tracks_max);
  TH1D * NGlobalGeoTracksPerEvent = new TH1D("NGlobalGeoTracksPerEvent",
    "total num. global + geometry-cut tracks per event",
    n_tracks_max,0,n_tracks_max);
  TH1D * NGlobalGeoNpTracksPerEvent = new TH1D("NGlobalGeoNpTracksPerEvent",
    "total num. global + geometry-cut + np>6 tracks per event",
    n_tracks_max,0,n_tracks_max);


  // geometry hists
  TString TrkThetaYvX_n[5];
  TString TrkThetaYvX_t[5];
  TString TrkThetaRpYvX_n[5];
  TString TrkThetaRpYvX_t[5];
  TH2D * TrkThetaYvX[5];
  TH2D * TrkThetaRpYvX[5];

  for(int bb=0; bb<5; bb++) {
    TrkThetaYvX_n[bb] = "TrkThetaYvX_"+BranchName[bb];
    TrkThetaYvX_t[bb] = "[glo] track #theta_{Y} vs. #theta_{X} :: "+BranchName[bb]+" branch(es);#theta_{X};#theta_{Y}";
    TrkThetaRpYvX_n[bb] = "TrkThetaRpYvX_"+BranchName[bb];
    TrkThetaRpYvX_t[bb] = "[glo] track #theta_{Y}^{RP} vs. #theta_{X}^{RP} :: "+BranchName[bb]+" branch(es);#theta_{X}^{RP};#theta_{Y}^{RP}";

    TrkThetaYvX[bb] = new TH2D(TrkThetaYvX_n[bb],
                               TrkThetaYvX_t[bb],
                               400,-4e-2,4e-2,
                               400,-4e-2,4e-2);
    TrkThetaRpYvX[bb] = new TH2D(TrkThetaRpYvX_n[bb],
                                 TrkThetaRpYvX_t[bb],
                                 400,-6e-2,6e-2,
                                 400,-6e-2,6e-2);
  };

  TH2D * TrkTheta01v2 = new TH2D("TrkTheta01v2",
    "[glo] track #sqrt{#theta_{X}^{2}+#theta_{Y}^{2}} vs. #theta",
    50,-4e-2,4e-2,
    50,-4e-2,4e-2);
  TH2D * TrkThetaRp01v2 = new TH2D("TrkThetaRp01v2",
    "[glo] track #sqrt{#theta_{X}^{RP 2}+#theta_{Y}^{RP 2}} vs. #theta^{RP}",
    50,-6e-2,6e-2,
    50,-6e-2,6e-2);


  // momentum hists
  TH1D * TrkP = new TH1D("TrkP",
    "[glo+geo+np] track p;p",
    300,0,200);
  TH1D * TrkPt = new TH1D("TrkPt",
    "[glo+geo+np] track p_{T};p_{T}",
    300,0,2);


  // branch bar chart
  TH1D * TrkInBranch = new TH1D("TrkInBranch",
    "[glo+geo+np] tracks in branches",
    4,0,4);
  for(int bb=1; bb<5; bb++) 
    TrkInBranch->GetXaxis()->SetBinLabel(bb,BranchName[bb].Data());
  TH2D * TrkBranchVsNTracks = new TH2D("TrkBranchVsNTracks",
    "[glo+geo+np] Branch vs. num. tracks;num tracks;branch",
    200,0,200,
    4,0,4);
  for(int bb=1; bb<5; bb++) 
    TrkBranchVsNTracks->GetYaxis()->SetBinLabel(bb,BranchName[bb].Data());



  // eta phi, theta 
  TH2D * TrkEtaPhi = new TH2D(
    "TrkEtaPhi",
    "[glo+geo+np] track #phi^{RP} vs. #eta",
    100,-10,10,
    100,-3.14,3.14
  );
  TH2D * TrkThetaPhi = new TH2D(
    "TrkThetaPhi",
    "[glo+geo+np] track #phi^{RP} vs. #theta^{RP}",
    200,-3e-2,3e-2,
    100,-3.14,3.14
  );



  // t, xi, time
  TH1D * TrkT = new TH1D("TrkT",
    "[glo+geo+np] track -t;-t",
    200,0,2);
  TH1D * TrkXi = new TH1D("TrkXi",
    "[glo+geo+np] track #xi;#xi",
    500,-3,3); 
  TH1D * TrkTime = new TH1D("TrkTime",
    "[glo+geo+np] track <time> [ns];<time> [ns]",
    200,-100,100);
  TH2D * TrkTvPt = new TH2D("TrkTvPt",
    "[glo+geo+np] track -t vs. p_{T}^{2};p_{T}^{2};-t",
    200,0,2,
    200,0,2);


  // trackpoint hists
  TH2D * TptY0vX0 = new TH2D("TptY0vX0",
    "[glo+geo+np] 1st track pt. Y vs. X;X_{1};Y_{1}",
    50,-0.1,0.1,50,-0.1,0.1);
  TH2D * TptY1vX1 = new TH2D("TptY1vX1",
    "[glo+geo+np] 2nd track pt. Y vs. X;X_{2};Y_{2}",
    50,-0.1,0.1,50,-0.1,0.1);
  TH2D * TptXvX = new TH2D("TptXvX",
    "[glo+geo+np] 1st track pt. X vs. 2nd track pt. X;X_{2};X_{1}",
    50,-0.05,0.1,50,-0.05,0.1);
  TH2D * TptYvY = new TH2D("TptYvY",
    "[glo+geo+np] 1st track pt. Y vs. 2nd track pt. Y;Y_{2};Y_{1}",
    50,-0.1,0.1,50,-0.1,0.1);
  TH2D * TptYDvsXD = new TH2D("TptYDvsXD",
    "[glo+geo+np] track pt. Y_{2}-Y_{1} vs. X_{2}-X_{1};X_{2}-X_{1};Y_{2}-Y_{1}",
    50,-0.05,0.05,50,-0.05,0.05);


  // num planes vs. num tracks
  TH2D * TptNPlanesVsNTracks = new TH2D("TptNPlanesVsNTracks ",
    "[glo+geo] 10 #times (num. hit planes in 1^{st} vessel) + (num. hit planes in 2^{nd} vessel) vs. num. tracks;num. tracks;10N_{planes(1)}+N_{planes(2)}",
    n_tracks_max,0,n_tracks_max,
    50,0,50);
  TH2D * TrkNPlanesVsNTracks = new TH2D("TrkNPlanesVsNTracks ",
    "[glo+geo] number of hit planes in track vs. num. tracks;num. tracks;num. planes",
    n_tracks_max,0,n_tracks_max,
    9,0,9);


  // num planes vs. kinematics
  TH2D * TrkNPlanesVsP = new TH2D("TrkNPlanesVsP",
    "[glo+geo] number of hit planes in track vs. p;p;num. planes",
    100,0,200,
    5,4,9);
  TH2D * TrkNPlanesVsPt = new TH2D("TrkNPlanesVsPt",
    "[glo+geo] number of hit planes in track vs. p_{T};p_{T};num. planes",
    100,0,2,
    5,4,9);
  TH2D * TrkNPlanesVsPhi = new TH2D("TrkNPlanesVsPhi",
    "[glo+geo] number of hit planes in track vs. #phi^{RP};#phi^{RP};num. planes",
    100,-3.14,3.14,
    5,4,9);
  TH2D * TrkNPlanesVsPhiU = new TH2D("TrkNPlanesVsPhiU",
    "[glo+geo] number of hit planes in track vs. #phi_{up}^{RP};#phi_{up}^{RP};num. planes",
    100,-3.14,3.14,
    5,4,9);
  TH2D * TrkNPlanesVsPhiD = new TH2D("TrkNPlanesVsPhiD",
    "[glo+geo] number of hit planes in track vs. #phi_{down}^{RP};#phi_{down}^{RP};num. planes",
    100,-3.14,3.14,
    5,4,9);
  TH2D * TrkNPlanesVsEta = new TH2D("TrkNPlanesVsEta",
    "[glo+geo] number of hit planes in track vs. #eta;#eta;num. planes",
    100,-10,10,
    5,4,9);
  TH2D * TrkNPlanesVsTheta = new TH2D("TrkNPlanesVsTheta",
    "[glo] number of hit planes in track vs. #theta^{RP};#theta^{RP};num. planes",
    200,-0.1,0.1,
    5,4,9);
  TH2D * TrkNPlanesVsThetaX = new TH2D("TrkNPlanesVsThetaX",
    "[glo] number of hit planes in track vs. #theta_{x}^{RP};#theta_{x}^{RP};num. planes",
    100,-6e-2,6e-2,
    5,4,9);
  TH2D * TrkNPlanesVsThetaY = new TH2D("TrkNPlanesVsThetaY",
    "[glo] number of hit planes in track vs. #theta_{y}^{RP};#theta_{y}^{RP};num. planes",
    100,-6e-2,6e-2,
    5,4,9);


  // track index vs. kinematics
  TH2D * TrkIdxVsP = new TH2D("TrkIdxVsP",
    "[glo+geo+np] track index vs. p;p;track index",
    100,0,200,
    101,0,101);
  TH2D * TrkIdxVsPt = new TH2D("TrkIdxVsPt",
    "[glo+geo+np] track index vs. p_{T};p_{T};track index",
    100,0,2,
    101,0,101);
  TH2D * TrkIdxVsPhi = new TH2D("TrkIdxVsPhi",
    "[glo+geo+np] track index vs. #phi^{RP};#phi^{RP};track index",
    100,-3.14,3.14,
    101,0,101);
  TH2D * TrkIdxVsPhiU = new TH2D("TrkIdxVsPhiU",
    "[glo+geo+np] track index vs. #phi_{up}^{RP};#phi_{up}^{RP};track index",
    100,-3.14,3.14,
    101,0,101);
  TH2D * TrkIdxVsPhiD = new TH2D("TrkIdxVsPhiD",
    "[glo+geo+np] track index vs. #phi_{down}^{RP};#phi_{down}^{RP};track index",
    100,-3.14,3.14,
    101,0,101);
  TH2D * TrkIdxVsEta = new TH2D("TrkIdxVsEta",
    "[glo+geo+np] track index vs. #eta;#eta;track index",
    100,-10,10,
    101,0,101);
  TH2D * TrkIdxVsTheta = new TH2D("TrkIdxVsTheta",
    "[glo] track index vs. #theta^{RP};#theta^{RP};track index",
    200,-0.1,0.1,
    101,0,101);
  TH2D * TrkIdxVsThetaX = new TH2D("TrkIdxVsThetaX",
    "[glo] track index vs. #theta_{x}^{RP};#theta_{x}^{RP};track index",
    100,-6e-2,6e-2,
    101,0,101);
  TH2D * TrkIdxVsThetaY = new TH2D("TrkIdxVsThetaY",
    "[glo] track index vs. #theta_{y}^{RP};#theta_{y}^{RP};track index",
    100,-6e-2,6e-2,
    101,0,101);



  // track momentum correlations
  const Int_t NCORR = 5; // max track index for pairwise correlation plots
  TH2D * TrkCorrP[NCORR][NCORR];
  TString TrkCorrP_n;
  TString TrkCorrP_t;
  for(int cj=0; cj<NCORR; cj++) {
    for(int ck=0; ck<NCORR; ck++) {
      TrkCorrP_n = Form("TrkCorrP_%d_%d",cj,ck);
      TrkCorrP_t = Form("[glo+geo+np] track p correlation (track %d p vs. track %d p);track %d p;track %d p",cj,ck,ck,cj);
      TrkCorrP[cj][ck] = new TH2D(TrkCorrP_n.Data(),TrkCorrP_t.Data(),100,0,200,100,0,200);
    };
  };

  // define "lists"
  const Int_t NIDX = 20;  // max track index for momentum storage list
  Double_t MomList[NIDX];
  Int_t CluIdList[2][NIDX];
  Int_t BranchList[NIDX];
  Int_t diff;

  // subsequent track momentum differences
  if(NIDX<NCORR) {
    fprintf(stderr,"ERROR: NIDX<NCORR\n");
    return;
  };
  TH2D * TrkCorrPdiffVsIdx = new TH2D("TrkCorrPdiffVsIdx",
    "[glo+geo+np] p_{0}-p_{idx}(>5 GeV) vs. idx;idx;p_{0}-p_{idx}",
    NIDX-1,1,NIDX,
    100,5,60);
  Double_t zoom_cutoff=1;
  TH2D * TrkCorrPdiffVsIdxZoom = new TH2D("TrkCorrPdiffVsIdxZoom",
    "[glo+geo+np] p_{0}-p_{idx}+(10 MeV) (zoomed-in) vs. idx;idx;p_{0}-p_{idx}+0.01",
    NIDX-1,1,NIDX,
    300,0,zoom_cutoff);
  

  // cluster index differences
  TH1D * TrkCluIdDiff[NCORR][NCORR];
  TString TrkCluIdDiff_n;
  TString TrkCluIdDiff_t;
  for(int cj=0; cj<NCORR-1; cj++) {
    for(int ck=cj+1; ck<NCORR; ck++) {
      TrkCluIdDiff_n = Form("TrkCluIdDiff_%d_%d",cj,ck);
      TrkCluIdDiff_t = Form("[glo+geo+np] track %d and %d clustID difference hash;(both tracks in same branch)",ck,cj);
      TrkCluIdDiff[cj][ck] = new TH1D(
        TrkCluIdDiff_n.Data(),TrkCluIdDiff_t.Data(),
        3000,0,3000
      );
    };
  };

  TH2D * TrkCluIdDiffVsP[NCORR][NCORR];
  TString TrkCluIdDiffVsP_n;
  TString TrkCluIdDiffVsP_t;
  for(int cj=0; cj<NCORR-1; cj++) {
    for(int ck=cj+1; ck<NCORR; ck++) {
      TrkCluIdDiffVsP_n = Form("TrkCluIdDiffVsP_%d_%d",cj,ck);
      TrkCluIdDiffVsP_t = Form("[glo+geo+np] track %d and %d clustID difference hash vs track p (both tracks in same branch);#Delta p [GeV];hash",ck,cj);
      TrkCluIdDiffVsP[cj][ck] = new TH2D(
        TrkCluIdDiffVsP_n.Data(),TrkCluIdDiffVsP_t.Data(),
        100,-40,60,
        3000,0,3000
      );
    };
  };


  // quality hash distribution
  TH1D * TrkQualHash;
  TString TrkQualHash_n = "TrkQualHash";
  TString TrkQualHash_t = "track quality hash distribution";
  TrkQualHash = new TH1D(
    TrkQualHash_n.Data(),
    TrkQualHash_t.Data(),
    2000,0,40000
  );

    


  // track points vessel correlations (sanity check)
  TH2D * TptInVessels = new TH2D("TptInVessels",
    "[glo+geo+np] track point 0 vessel vs. track point 1 vessel;track point 1 vessel;track point 0 vessel",
    8,0,8,
    8,0,8);
  TptInVessels->GetXaxis()->SetBinLabel(1,"E1U");
  TptInVessels->GetXaxis()->SetBinLabel(2,"E1D");
  TptInVessels->GetXaxis()->SetBinLabel(3,"E2U");
  TptInVessels->GetXaxis()->SetBinLabel(4,"E2D");
  TptInVessels->GetXaxis()->SetBinLabel(5,"W1U");
  TptInVessels->GetXaxis()->SetBinLabel(6,"W1D");
  TptInVessels->GetXaxis()->SetBinLabel(7,"W2U");
  TptInVessels->GetXaxis()->SetBinLabel(8,"W2D");

  TptInVessels->GetYaxis()->SetBinLabel(1,"E1U");
  TptInVessels->GetYaxis()->SetBinLabel(2,"E1D");
  TptInVessels->GetYaxis()->SetBinLabel(3,"E2U");
  TptInVessels->GetYaxis()->SetBinLabel(4,"E2D");
  TptInVessels->GetYaxis()->SetBinLabel(5,"W1U");
  TptInVessels->GetYaxis()->SetBinLabel(6,"W1D");
  TptInVessels->GetYaxis()->SetBinLabel(7,"W2U");
  TptInVessels->GetYaxis()->SetBinLabel(8,"W2D");





  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  
  Int_t n_global=0;
  Int_t n_gold=0;
  Int_t n_global_geo=0;
  Int_t n_global_geo_np=0;

  Int_t trk_idx;

  const Double_t thetaRpLimits[2][2] = { 
    {-1.5e-3, 5.0e-3},
    { 1.0e-3, 5.5e-3} 
  }; // [rad]; theta limits from rafal
  /*
  TLine * limLine[2][4]; 
  limLine[0][0] = new TLine(thetaRpLimits[0][0],thetaRpLimits[1][1],thetaRpLimits[0][1],thetaRpLimits[1][1]);
  limLine[0][1] = new TLine(thetaRpLimits[0][0],thetaRpLimits[1][1],thetaRpLimits[0][0],thetaRpLimits[1][0]);
  limLine[0][2] = new TLine(thetaRpLimits[0][0],thetaRpLimits[1][0],thetaRpLimits[0][1],thetaRpLimits[1][0]);
  limLine[0][3] = new TLine(thetaRpLimits[0][1],thetaRpLimits[1][1],thetaRpLimits[0][1],thetaRpLimits[1][0]);
  limLine[1][0] = new TLine(thetaRpLimits[0][0],-1*thetaRpLimits[1][1],thetaRpLimits[0][1],-1*thetaRpLimits[1][1]);
  limLine[1][1] = new TLine(thetaRpLimits[0][0],-1*thetaRpLimits[1][1],thetaRpLimits[0][0],-1*thetaRpLimits[1][0]);
  limLine[1][2] = new TLine(thetaRpLimits[0][0],-1*thetaRpLimits[1][0],thetaRpLimits[0][1],-1*thetaRpLimits[1][0]);
  limLine[1][3] = new TLine(thetaRpLimits[0][1],-1*thetaRpLimits[1][1],thetaRpLimits[0][1],-1*thetaRpLimits[1][0]);
  for(int ii=0; ii<2; ii++) {
    for(int jj=0; jj<4; jj++) {
      limLine[ii][jj]->SetLineWidth(4);
      limLine[ii][jj]->SetLineColor(kBlack);
    };
  };
  */

  Bool_t vesselFilter = false;

  Bool_t printRPheader = true;

  Int_t vg_cnt=0;

  Int_t ENT = p_out->GetEntries();
  //ENT = 100000;
  for(int x=0; x<ENT; x++) {
    if(x%1000==0) printf("%.2f%%\n",(float)x/ENT*100);
    p_out->GetEntry(x);

    if(print_evids) printf("evid=%d\n",evid);

    // reset some counters / booleans / internal data structures
    printRPheader = true;
    n_global=0;
    n_gold=0;
    n_global_geo=0;
    n_global_geo_np=0;
    for(int cc=0; cc<NIDX; cc++) { 
      MomList[cc]=-1000;
      for(int ll=0; ll<2; ll++) CluIdList[ll][cc]=-1000;
      BranchList[cc]=-1000;
    };
    

    if(n_tracks>0) NTracksPerEvent->Fill(n_tracks);


    // -------------------------------- -------------------------------------------------
    // track loop
    // -------------------------------- -------------------------------------------------
    for(int t=0; t<n_tracks; t++) {

      // filter east or west or both (et al)
      vesselFilter = false;
      switch(ewb012) {
        case 0: 
          if(t_branch[t]==0 || t_branch[t]==1) vesselFilter = true;
          break;
        case 1:
          if(t_branch[t]==2 || t_branch[t]==3) vesselFilter = true;
          break;
        case 2:
          if(t_branch[t]>=0 || t_branch[t]<=3) vesselFilter = true;
          break;
        default:
          fprintf(stderr,"ERROR: ewb012 is not correctly defined\n");
          return;
      };

      if(vesselFilter && !(t_isBad[t])) TrkQualHash->Fill(t_qualHash[t]);


      // global & ewb012 cut & !isBad - - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - 
      if(t_type[t]==rpsGlobal && vesselFilter && !(t_isBad[t])) {

        n_global++;
        if(t_gold[t]) n_gold++; 


        // fill geometry (theta_x, theta_y, theta) plots
        TrkThetaYvX[0]->Fill(t_theta[t][rpsAngleThetaX],
                             t_theta[t][rpsAngleThetaY]);
        TrkThetaRpYvX[0]->Fill(t_thetaRP[t][rpsAngleThetaX],
                               t_thetaRP[t][rpsAngleThetaY]);
        if(t_branch[t]>=0 && t_branch[t]<4) { 
          TrkThetaYvX[t_branch[t]+1]->Fill(t_theta[t][rpsAngleThetaX],
                                           t_theta[t][rpsAngleThetaY]);
          TrkThetaRpYvX[t_branch[t]+1]->Fill(t_thetaRP[t][rpsAngleThetaX],
                                             t_thetaRP[t][rpsAngleThetaY]);
        };

        TrkTheta01v2->Fill(TMath::Sqrt(
                   TMath::Power(t_theta[t][rpsAngleThetaX],2)+
                   TMath::Power(t_theta[t][rpsAngleThetaY],2)),
                           t_theta[t][rpsAngleTheta]);
        TrkThetaRp01v2->Fill(TMath::Sqrt(
                   TMath::Power(t_thetaRP[t][rpsAngleThetaX],2)+
                   TMath::Power(t_thetaRP[t][rpsAngleThetaY],2)),
                           t_thetaRP[t][rpsAngleTheta]);
        

        // fill nplanes vs. theta plots
        TrkNPlanesVsTheta->Fill(
          t_thetaRP[t][rpsAngleTheta],
          t_planesUsed[t]
        );
        TrkNPlanesVsThetaX->Fill(
          t_thetaRP[t][rpsAngleThetaX],
          t_planesUsed[t]
        );
        TrkNPlanesVsThetaY->Fill(
          t_thetaRP[t][rpsAngleThetaY],
          t_planesUsed[t]
        );

        
        // fill idx vs. theta plots (done before geometry cut)
        TrkIdxVsTheta->Fill(
          t_thetaRP[t][rpsAngleTheta],
          t_index[t]
        );
        TrkIdxVsThetaX->Fill(
          t_thetaRP[t][rpsAngleThetaX],
          t_index[t]
        );
        TrkIdxVsThetaY->Fill(
          t_thetaRP[t][rpsAngleThetaY],
          t_index[t]
        );

        


        // geometry cut - - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - 
        if( 
          t_thetaRP[t][rpsAngleThetaX] > thetaRpLimits[0][0] &&
          t_thetaRP[t][rpsAngleThetaX] < thetaRpLimits[0][1] &&
          TMath::Abs(t_thetaRP[t][rpsAngleThetaY]) > thetaRpLimits[1][0] &&
          TMath::Abs(t_thetaRP[t][rpsAngleThetaY]) < thetaRpLimits[1][1]
        ) {

          n_global_geo++;


          // fill nplanes vs ntracks hist
          if( (p_RPid[0][t] % 4) <= 1 &&
              (p_RPid[1][t] % 4) > 1) 
            TptNPlanesVsNTracks->Fill(
              n_tracks,
              10*p_planesUsed[0][t]+p_planesUsed[1][t]
            );
          else if(p_RPid[1][t] % 4 <= 1 &&
                  p_RPid[0][t] % 4 > 1) 
          {
            TptNPlanesVsNTracks->Fill(
              n_tracks,
              10*p_planesUsed[1][t]+p_planesUsed[0][t]
            );
            fprintf(stderr,"weird case\n");
          }
          else fprintf(stderr,"weirder case\n");


          // fill nplanes vs kinematics
          TrkNPlanesVsNTracks->Fill(n_tracks,t_planesUsed[t]);
          TrkNPlanesVsP->Fill(t_p[t],t_planesUsed[t]);
          TrkNPlanesVsPt->Fill(t_pt[t],t_planesUsed[t]);
          TrkNPlanesVsPhi->Fill(t_phiRP[t],t_planesUsed[t]);
          if(t_branch[t]==0 || t_branch[t]==2) 
            TrkNPlanesVsPhiU->Fill(t_phiRP[t],t_planesUsed[t]);
          else if(t_branch[t]==1 || t_branch[t]==3) 
            TrkNPlanesVsPhiD->Fill(t_phiRP[t],t_planesUsed[t]);
          TrkNPlanesVsEta->Fill(t_eta[t],t_planesUsed[t]);


          // nplanes cut - - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - - - -- - 
          if(t_planesUsed[t]>6) {

            n_global_geo_np++;

            trk_idx = n_global_geo_np-1;


            TrkP->Fill(t_p[t]); 
            TrkPt->Fill(t_pt[t]);

            TrkInBranch->Fill(t_branch[t]);
            TrkBranchVsNTracks->Fill(n_tracks,t_branch[t]);

            TrkEtaPhi->Fill(t_eta[t],t_phiRP[t]);
            TrkThetaPhi->Fill(
              t_thetaRP[t][rpsAngleTheta],
              t_phiRP[t]
            );

            TrkT->Fill(-1*t_t[t]);
            TrkXi->Fill(t_xi[t]);
            TrkTime->Fill(t_time[t]*1e9);
            TrkTvPt->Fill(t_pt[t]*t_pt[t],-1*t_t[t]);

            TptY0vX0->Fill(p_x[0][t],p_y[0][t]);
            TptY1vX1->Fill(p_x[1][t],p_y[1][t]);
            TptXvX->Fill(p_x[1][t],p_x[0][t]);
            TptYvY->Fill(p_y[1][t],p_y[0][t]);
            TptYDvsXD->Fill(p_x[1][t]-p_x[0][t],p_y[1][t]-p_y[0][t]);


            // fill vessels correlations sanity check
            TptInVessels->Fill(p_RPid[1][t],p_RPid[0][t]);


            // fill track index vs. kinematics
            TrkIdxVsP->Fill(t_p[t],t_index[t]);
            TrkIdxVsPt->Fill(t_pt[t],t_index[t]);
            TrkIdxVsPhi->Fill(t_phiRP[t],t_index[t]);
            if(t_branch[t]==0 || t_branch[t]==2) 
              TrkIdxVsPhiU->Fill(t_phiRP[t],t_index[t]);
            else if(t_branch[t]==1 || t_branch[t]==3) 
              TrkIdxVsPhiD->Fill(t_phiRP[t],t_index[t]);
            TrkIdxVsEta->Fill(t_eta[t],t_index[t]);
          

            // add track momentum to momentum list; 
            // add cluster id to cluster id list;
            // use current n_global_geo_np-1 as index
            if(trk_idx < NIDX  && trk_idx >= 0) {
              MomList[trk_idx] = t_p[t];
              BranchList[trk_idx] = t_branch[t];

              for(int ll=0; ll<2; ll++) {
                if(p_clustid_s1[ll][t]<10 &&
                   p_clustid_s2[ll][t]<10 &&
                   p_clustid_s3[ll][t]<10 &&
                   p_clustid_s4[ll][t]<10) {
                  CluIdList[ll][trk_idx] = 1000*(p_clustid_s1[ll][t]+1)+
                                            100*(p_clustid_s2[ll][t]+1)+
                                             10*(p_clustid_s3[ll][t]+1)+
                                              1*(p_clustid_s4[ll][t]+1);
                };
              };
            };




            // print out RP track information
            if(printRPevents) {
              if(printRPheader) {
                printf("\n evid=%d  n_tracks=%d\n-----------------------------------------\n",
                  evid,n_tracks);
                printRPheader = false;
              };


              // this print statement really slows things down, even if printRPevents==false
              // it was originally written to look at cluster indices
              printf("idx=%d  p=%f  pt=%f  np=%d   clustid[0]=(%d, %d, %d, %d)   clustid[1]=(%d, %d, %d, %d)\n",
                t_index[t],t_p[t],t_pt[t],t_planesUsed[t],
                p_clustid_s1[0][t],p_clustid_s2[0][t],p_clustid_s3[0][t],p_clustid_s4[0][t],
                p_clustid_s1[1][t],p_clustid_s2[1][t],p_clustid_s3[1][t],p_clustid_s4[1][t]
              );
            };

          }; // eo track nplanes cut
        }; // eo track geometry cut
      }; // eo track global cut (and in vessel cut)
    }; // eo track loop



    // after enumerating tracks in a single event 
    if(n_global>0) NGlobalTracksPerEvent->Fill(n_global);
    if(n_gold>0) NGoldTracksPerEvent->Fill(n_gold);
    if(n_global_geo>0) NGlobalGeoTracksPerEvent->Fill(n_global_geo);
    if(n_global_geo_np>0) NGlobalGeoNpTracksPerEvent->Fill(n_global_geo_np);


    // fill plots which look at correlations between subsequent plots; for this we need to have at least 2 
    // "good" tracks to look at, which is why this n_global_geo_np>1 cut is used 
    if(n_global_geo_np>1) {

      // fill track momentum correlations (correlations between first, second, third, etc. track momenta)
      for(int cj=0; cj<NCORR-1; cj++) {
        if(MomList[cj]>=0) {
          for(int ck=cj+1; ck<NCORR; ck++) {
            if(MomList[ck]>=0) {
              TrkCorrP[cj][ck]->Fill(MomList[ck],MomList[cj]);
              //printf("cj=%d ck=%d\n",cj,ck);
              //printf("ck=%d\n",ck);
            };
          };
        };
      };

      // fill subsequent track momentum differences
      if(MomList[0]>0) {
        for(int cs=1; cs<NIDX; cs++) {
          if(MomList[cs]>0) {
            TrkCorrPdiffVsIdx->Fill(cs,MomList[0]-MomList[cs]);
            TrkCorrPdiffVsIdxZoom->Fill(cs,MomList[0]-MomList[cs]+0.01);
            //printf("cs=%d\n",cs);
          };
        };
      };


      // fill cluster index plots
      for(int cj=0; cj<NCORR-1; cj++) {
        if(CluIdList[0][cj]>0 && CluIdList[1][cj]>0) {
          for(int ck=cj+1; ck<NCORR; ck++) {
            if(CluIdList[0][ck]>0 && CluIdList[1][ck]>0) {
              if(BranchList[cj] == BranchList[ck]) {
                diff = HashCluIds(CluIdList[0][cj],CluIdList[1][cj],
                                  CluIdList[0][ck],CluIdList[1][ck]);
                TrkCluIdDiff[cj][ck]->Fill(diff);
                if(MomList[cj]>0 && MomList[ck]>0)
                  TrkCluIdDiffVsP[cj][ck]->Fill(MomList[cj]-MomList[ck],diff);
              };
            };
          };
        };
      };

      // print lists (for debugging) 
      if(false) {
        printf("CluIdList:\n");
        for(int ccc=0; ccc<NIDX; ccc++) {
          diff = 0;
          if(MomList[ccc]>0) {
            if(ccc+1<NIDX) {
              if(CluIdList[0][ccc+1]>0 && CluIdList[1][ccc+1]>0) {
                diff = HashCluIds(CluIdList[0][ccc],   CluIdList[1][ccc],
                                  CluIdList[0][ccc+1], CluIdList[1][ccc+1]);
              };
            };
          };
          printf(" %d --- %d . %d -- diff=%d -- br=%d -- p=%f\n",ccc,
            CluIdList[0][ccc],
            CluIdList[1][ccc],
            diff,
            BranchList[ccc],
            MomList[ccc]);
        };
      };


      // loop limiter, for valgrind tests
      //vg_cnt++;
      //if(vg_cnt>5) break;
    };


  }; // eo p_out tree loop



  // - - --- - --- - --- - --- - --- - --- - --- - --- - --- - --
  // - - --- - --- - --- - --- - --- - --- - --- - --- - --- - --
  // - - --- - --- - --- - --- - --- - --- - --- - --- - --- - --


  // write output ======================================================
  TString outfile_n = filename;
  TString filesuffix;
  switch(ewb012) {
    case 0: 
      filesuffix = ".diag.east.root";
      break;
    case 1:
      filesuffix = ".diag.west.root";
      break;
    case 2:
      filesuffix = ".diag.bothsides.root";
      break;
  };
  if(use_diag_dir) {
    outfile_n = outfile_n.ReplaceAll(".root",filesuffix.Data());
    outfile_n = outfile_n.ReplaceAll("/OFset","_diag/OFset");
  } else {
    filename.ReplaceAll(".root",filesuffix.Data());
  };
  printf("%s created\n",outfile_n.Data());
  TFile * outfile = new TFile(outfile_n.Data(),"RECREATE");

  NTracksPerEvent->Write();
  NGlobalTracksPerEvent->Write();
  NGoldTracksPerEvent->Write();
  NGlobalGeoTracksPerEvent->Write();
  NGlobalGeoNpTracksPerEvent->Write();

  TrkInBranch->Write();
  TrkBranchVsNTracks->Write();

  for(int bb=0; bb<5; bb++) TrkThetaYvX[bb]->Write();
  for(int bb=0; bb<5; bb++) TrkThetaRpYvX[bb]->Write();

  /*
  TCanvas * c_TrkThetaRpYvX = new TCanvas("c_TrkThetaRpYvX","c_TrkThetaRpYvX",800,800);
  c_TrkThetaRpYvX->SetLogz();
  TrkThetaRpYvX[0]->Draw("colz");
  for(int ii=0; ii<2; ii++) {
    for(int jj=0; jj<4; jj++) {
      limLine[ii][jj]->Draw();
    };
  };
  */
  // c_TrkThetaRpYvX->Write(); // deprecated; use DrawGeometryBoxes.C instead!

  TrkP->Write();
  TrkPt->Write();
  TrkEtaPhi->Write();
  TrkThetaPhi->Write();

  TrkT->Write();
  TrkXi->Write();
  TrkTime->Write();
  TrkTvPt->Write();

  TptNPlanesVsNTracks->Write();
  TrkNPlanesVsNTracks->Write();
  TrkNPlanesVsP->Write();
  TrkNPlanesVsPt->Write();
  TrkNPlanesVsPhi->Write();
  TrkNPlanesVsPhiU->Write();
  TrkNPlanesVsPhiD->Write();
  TrkNPlanesVsEta->Write();
  TrkNPlanesVsTheta->Write();
  TrkNPlanesVsThetaX->Write();
  TrkNPlanesVsThetaY->Write();

  TrkIdxVsP->Write();
  TrkIdxVsPt->Write();
  TrkIdxVsPhi->Write();
  TrkIdxVsPhiU->Write();
  TrkIdxVsPhiD->Write();
  TrkIdxVsEta->Write();
  TrkIdxVsTheta->Write();
  TrkIdxVsThetaX->Write();
  TrkIdxVsThetaY->Write();

  outfile->mkdir("p_correlations");
  outfile->cd("p_correlations");
  for(int cj=0; cj<NCORR-1; cj++) {
    for(int ck=cj+1; ck<NCORR; ck++) {
      TrkCorrP[cj][ck]->Write();
    };
  };

  outfile->mkdir("clustid_hists");
  outfile->cd("clustid_hists");
  for(int cj=0; cj<NCORR-1; cj++) {
    for(int ck=cj+1; ck<NCORR; ck++) {
      TrkCluIdDiff[cj][ck]->Write();
    };
  };

  outfile->mkdir("clustid_v_p");
  outfile->cd("clustid_v_p");
  for(int cj=0; cj<NCORR-1; cj++) {
    for(int ck=cj+1; ck<NCORR; ck++) {
      TrkCluIdDiffVsP[cj][ck]->Write();
    };
  };

  outfile->cd();
  TrkCorrPdiffVsIdx->Write();
  TrkCorrPdiffVsIdxZoom->Write();

  TrkQualHash->Write();
  outfile->mkdir("extra_plots");
  outfile->cd("extra_plots");
  TptInVessels->Write();
  TrkTheta01v2->Write(); 
  TrkThetaRp01v2->Write();
  TptY0vX0->Write();
  TptY1vX1->Write();
  TptXvX->Write();
  TptYvY->Write();
  TptYDvsXD->Write();
  outfile->cd();


  // print statistics
  Double_t total = (Double_t)(p_out->GetEntries());
  Double_t subtotal = (Double_t)(NTracksPerEvent->GetEntries());
  printf("num. events with n_tracks>0: %.2f%% of total num. events\n",
    100.0*subtotal/total);
  printf("-----------------------------\ncut: (%% of n_tracks>0 events) (%% of total num. tracks)\n-----------------------------\n");
  printf("global tracks: %.2f%% (%.2f%% of total)\n",
    100.0*NGlobalTracksPerEvent->GetEntries()/subtotal,
    100.0*NGlobalTracksPerEvent->GetEntries()/total);
  /*
  printf("global + golden trackpoints tracks: %.2f%% (%.2f%% of total)\n",
    100.0*NGoldTracksPerEvent->GetEntries()/subtotal,
    100.0*NGoldTracksPerEvent->GetEntries()/total);
  */
  printf("global + geometry cut tracks: %.2f%% (%.2f%% of total)\n",
    100.0*NGlobalGeoTracksPerEvent->GetEntries()/subtotal,
    100.0*NGlobalGeoTracksPerEvent->GetEntries()/total);
};

#ifndef __CINT__
int main(int argc, char ** argv) {
  if(argc==1) DiagnosticRP();
  else if(argc==2) DiagnosticRP(TString(argv[1]));
  else DiagnosticRP(TString(argv[1]),(Int_t)strtof(argv[2],NULL));
};
#endif
