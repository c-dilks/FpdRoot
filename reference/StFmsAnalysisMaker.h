// \class StFmsAnalysisMaker
// \author Akio Ogawa, modified by dilks
//
//   This is analysis for FMS-FPS correlations.
// 
//  $Id: StFmsAnalysisMaker.h,v 1.5 2015/12/08 17:00:03 akio Exp $
//  $Log: StFmsAnalysisMaker.h,v $
//

#ifndef STAR_StFmsAnalysisMaker_HH
#define STAR_StFmsAnalysisMaker_HH

#include "StMaker.h"
#include "StEnumerations.h"
#include "StLorentzVectorF.hh"
#include "TMatrix.h"
#include "TVector3.h"


class StMuDstMaker;
class StFmsDbMaker;
class StFmsCollection;
class StMuFmsCollection;
class StMuRpsCollection;
class StSpinDbMaker;
class TClonesArray;
class StMuFmsCluster;
class StMuFmsPoint;
class StMuFmsHit;
class StMuRpsTrack;
class StMuRpsTrackPoint;
class TH2F;
//class TMatrixT;

const int MAX_POUT = 1264; // max size for p_out tree arrays
const int n_tracks_max = 1000;
const int n_ang = 3; // should be equal to StMuRpsTrack::mNumberOfAngleTypes; double-checked in Init() since I can't seem to access it here

class StFmsAnalysisMaker : public StMaker{
  public: 
    StFmsAnalysisMaker(StMuDstMaker * maker, 
                       StSpinDbMaker * spindb,
                       const Char_t* name="FmsOffQA");
    ~StFmsAnalysisMaker();
    Int_t Init();
    Int_t Make();
    Int_t Finish();

    void setFileName(char* file){mFilename=file;} 
    void setPrint(int v) {mPrint=v;}
    void ResetBranches(Int_t i_);


    static int NTracksMax() { return n_tracks_max; };

  private:
    StMuDstMaker * mMuDstMaker;
    StSpinDbMaker * mSpinDb;
    StFmsDbMaker* mFmsDbMaker;
    StMuFmsCollection* mMuFmsColl;
    StMuRpsCollection* mMuRpsColl;

    Int_t Rnum_tmp;

    char* mFilename;
    TFile* mFile;
    int mPrint;

    // tree for events
    //TTree * evtr;
    UInt_t bc[2]; // bxing counter
    Int_t evid; // event id
    Int_t nhits,nclusters,npoints; // number of hits, clusters, points, pairs
    Int_t nTriggers;
    Int_t trigids[128];


    // tree for clusters
    TTree * clustr;
    Int_t det_ctr;
    Int_t category_ctr,ntow_ctr,nphot_ctr;
    Float_t en_ctr,pt_ctr,x_ctr,y_ctr;
    Float_t sigmax_ctr,sigmin_ctr;
    Float_t chi2ndf_onephot_ctr;
    Float_t chi2ndf_twophot_ctr;
    Int_t id_ctr;

    // tree for point pairs
    TTree * pairtr;
    Int_t llss; // 0=large, 1=small
    Int_t nclust_in_pair; // number of clusters in pair
    Int_t npoints_in_pair; // number of points in pair
    Float_t en_tr,pt_tr,eta_tr,phi_tr,mass_tr,dgg_tr,zgg_tr,x_tr,y_tr; // pair kinematics
    Float_t cone_rad_tr,cone_en_tr,cone_enfrac_tr; // cone kinematics
    Float_t en_p_tr[2]; // pair point energy
    Float_t x_p_tr[2]; // pair point x
    Float_t y_p_tr[2]; // pair point y
    Float_t dist_from_edge_p_tr[2]; // distance of pair point from edge
    Int_t edge_type_p_tr[2]; // edge type nearby pair point
    // cluster variables (if one cluster, first entry of array == second)
    Int_t det_p_tr[2]; // detector of cluster associated to pair point
    Float_t sigmax_p_tr[2]; // sigma max of cluster(s)
    Float_t sigmin_p_tr[2]; // sigma min of cluster(s)
    Float_t chi2ndf_onephot_p_tr[2]; // chi2/ndf for 1 photon fit
    Float_t chi2ndf_twophot_p_tr[2]; // chi2/ndf for 2 photon fit
    Int_t clust_id_p_tr[2]; // cluster id?
    

    // OFile tree
    TTree * p_out;
    Int_t spin; //-------------- using v124 bits until spinqa complete
    Int_t nphotons; 
    Int_t nwrds; 
    Int_t tpes[MAX_POUT]; 
    Float_t pxyzt[MAX_POUT]; 
    Int_t Rnum; 
    Int_t Bunchid7bit; 
    UChar_t BBcSums[5]; 
    Float_t BBcVertex[7]; 
    Int_t EventN; 
    Int_t ievt; 
    Int_t nSavedHits; 
    Int_t SavedHits[MAX_POUT]; 
    Int_t TrigBits; 
    Int_t nCluster; 
    Int_t nPhotonClu; 
    Int_t SCIndex[MAX_POUT]; 
    Int_t SPCIndex[MAX_POUT]; 
    Float_t SPCEnergy[MAX_POUT]; 
    UInt_t L2sum[2]; 
    UInt_t lastdsm[8]; 
    UInt_t Fpde[8]; 

    // OFile tree RP branches
    Int_t n_tracks,n_trackpoints;
    // track variables -- prefixed with t_
    Int_t t_index[n_tracks_max];
    Int_t t_branch[n_tracks_max];
    Int_t t_type[n_tracks_max];
    UInt_t t_planesUsed[n_tracks_max];
    Double_t t_p[n_tracks_max];
    Double_t t_pt[n_tracks_max];
    Double_t t_eta[n_tracks_max];
    Double_t t_time[n_tracks_max];
    Double_t t_theta[n_tracks_max][n_ang]; // [angle (X,Y,full)]
    Double_t t_thetaRP[n_tracks_max][n_ang]; // [angle (X,Y,full)]
    Double_t t_phi[n_tracks_max];
    Double_t t_phiRP[n_tracks_max];
    Double_t t_t[n_tracks_max];
    Double_t t_xi[n_tracks_max];
    Bool_t t_gold[n_tracks_max];
    Bool_t t_isBad[n_tracks_max];
    Double_t t_qualHash[n_tracks_max];
    // trackpoint variables -- prefixed with p0_ and p1_
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
    Double_t p_time_pmt1[2][n_tracks_max];
    Double_t p_time_pmt2[2][n_tracks_max];
    

    


    
    // other variables
    unsigned long long L2sum_full;
    Int_t spinbit;
    Int_t bNib,yNib;
    Int_t spinbyte,spinbyte_db;
    Int_t spinbyte_curr,spinbyte_db_curr;
    Bool_t spinbyte_same;
    char spc_b,spc_y;
    Int_t sp_b,sp_y;
    Bool_t dump_spinbits;
    Bool_t verbose,verbose_clu,verbose_pho,verbose_rp;
    unsigned short nstb;
    unsigned short chan;
    unsigned short hit_adc;
    float hit_en;
    int hit_r,hit_c,hit_n;

    Int_t trkSortIdx[n_tracks_max];
    Double_t trkSortQual[n_tracks_max];
    Double_t trkMom;
    Int_t trkPlanes;
    Int_t trkGlobal;
    Int_t trkInGeom;
    Double_t trkThetaX;
    Double_t trkThetaY;



    // object arrays
    TClonesArray * mHits;
    TClonesArray * mClustersAll;
    TObjArray * mClusters[4];
    TClonesArray * mPoints;

    // looping objects
    StMuFmsCluster * clu;
    StMuFmsPoint * pho;
    StMuFmsHit * hit;
    StMuRpsTrack * trk;
    Int_t clu_id,pho_id;
    Float_t clu_en,clu_x,clu_y,clu_sigmin,clu_sigmax,clu_csqn1,clu_csqn2;
    Float_t pho_en,pho_x,pho_y;
    Int_t ent[4];
    Float_t en_sum_all;
    Float_t en_sum[4];
    Float_t clu_en_sortval[MAX_POUT];
    Int_t clu_en_sortind[MAX_POUT];
    StLorentzVectorF vecs[MAX_POUT];
    Int_t num_null;


    // from root12fms/FpdRoot::Geom class
    Float_t FpdTowWid[2][4]; // [ew][nstb]
    Float_t ZFPD[2][4]; // [ew][nstb]
    Float_t xOffset[2][4]; // [ew][nstb]
    Float_t yOffset[2][4]; // [ew][nstb]
    Float_t EnergySum,Energy;

    Int_t nRows[4]; // [nstb]
    Int_t nCols[4]; // [nstb]

    TMatrix * eMat[4]; // [nstb]

    unsigned int p1,p2,p3,p4,pa;
    

    virtual const char *GetCVS() const
    {static const char cvs[]="Tag $Name:  $ $Id: StFmsAnalysisMaker.h,v 1.5 2015/12/08 17:00:03 akio Exp $ built " __DATE__ " " __TIME__ ; return cvs;}

    ClassDef(StFmsAnalysisMaker,0);
};

#endif
