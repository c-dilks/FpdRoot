// \class StFmsAnalysisMaker
// \author Akio Ogawa
//
//  $Id: StFmsAnalysisMaker.cxx,v 1.6 2015/12/08 17:00:03 akio Exp $
//

#include "StFmsAnalysisMaker.h"

#include "StMessMgr.h"
#include "Stypes.h"

#include "StFmsDbMaker/StFmsDbMaker.h"
#include "StEnumerations.h"
#include "StEventTypes.h"

#include "StMuDSTMaker/COMMON/StMuTypes.hh"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuFmsCollection.h"
#include "StMuDSTMaker/COMMON/StMuFmsCluster.h"
#include "StMuDSTMaker/COMMON/StMuFmsHit.h"
#include "StMuDSTMaker/COMMON/StMuFmsPoint.h"

#include "StMuDSTMaker/COMMON/StMuRpsCollection.h"
#include "StMuDSTMaker/COMMON/StMuRpsTrackPoint.h"
#include "StMuDSTMaker/COMMON/StMuRpsTrack.h"

#include "StEvent/StFmsCollection.h"
#include "StEvent/StFmsHit.h"
#include "StEvent/StFmsPoint.h"
#include "StEvent/StFmsPointPair.h"
#include "StEvent/StEnumerations.h"

#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h"

#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TList.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TMath.h"
#include "TMatrixF.h"
#include "TVector3.h"
#include "TLorentzVector.h"


#include <bitset>

ClassImp(StFmsAnalysisMaker);



/*
StFmsAnalysisMaker::StFmsAnalysisMaker(const Char_t* name):
  StMaker(name),mFilename((char *)"fmsOffQa.root"),mPrint(1)
{}
*/
StFmsAnalysisMaker::StFmsAnalysisMaker(StMuDstMaker * maker,
                                       StSpinDbMaker * spindb,
                                       const Char_t* name):
  StMaker(name),mFilename((char *)"fmsOffQa.root"),mPrint(1) {
    LOG_INFO << "[++++++] begin StFmsAnalysisMaker instantiation..." << endm;
    mMuDstMaker = maker;
    mSpinDb = spindb;

    // logging switches
    dump_spinbits = true; // if true, dumps spin pattern during first event
    verbose = false; // if true, prints out event-by-event, cluster-by-cluster, hit-by-hit info
    verbose_clu = true; // if true, prints out kinematics for each cluster (if !verbose, then verbose_clu forced to be false)
    verbose_pho = true; // if true, prints out kinematics for each photon along with types[nwrds] and pxyzt[nwrds] (if !verbose, forced to be false)
    verbose_rp = false; // RP diagnostics

    if(!verbose) {
      verbose_clu = false;
      verbose_pho = false;
    };


    // geometry from root12fms/FpdRoot::Geom class; the following code was generated using 
    // heppellab16f:~/15tran/root12fms/PRINT_GEOM_TABLE.C; numbers are from geom.txt
    FpdTowWid[0][0] = 3.822000;
    FpdTowWid[0][1] = 3.822000;
    FpdTowWid[0][2] = 0.000000;
    FpdTowWid[0][3] = 0.000000;
    FpdTowWid[1][0] = 5.812000;
    FpdTowWid[1][1] = 5.812000;
    FpdTowWid[1][2] = 3.822000;
    FpdTowWid[1][3] = 3.822000;

    ZFPD[0][0] = -811.000000;
    ZFPD[0][1] = -811.000000;
    ZFPD[0][2] = 0.000000;
    ZFPD[0][3] = 0.000000;
    ZFPD[1][0] = 734.099976;
    ZFPD[1][1] = 734.099976;
    ZFPD[1][2] = 729.700012;
    ZFPD[1][3] = 729.700012;

    xOffset[0][0] = -20.400000;
    xOffset[0][1] = 20.400000;
    xOffset[0][2] = 0.000000;
    xOffset[0][3] = 0.000000;
    xOffset[1][0] = -0.300000;
    xOffset[1][1] = 0.300000;
    xOffset[1][2] = -0.930000;
    xOffset[1][3] = 0.930000;

    yOffset[0][0] = 13.380000;
    yOffset[0][1] = 13.380000;
    yOffset[0][2] = 0.000000;
    yOffset[0][3] = 0.000000;
    yOffset[1][0] = 98.800003;
    yOffset[1][1] = 98.800003;
    yOffset[1][2] = 46.500000;
    yOffset[1][3] = 46.500000;
    // -------------- end generated code
    
    nRows[0] = 34;
    nRows[1] = 34;
    nRows[2] = 24;
    nRows[3] = 24;

    nCols[0] = 17;
    nCols[1] = 17;
    nCols[2] = 12;
    nCols[3] = 12;

    for(int nn=0; nn<4; nn++) eMat[nn] = new TMatrix(nRows[nn],nCols[nn]);

    if(verbose) {
      LOG_INFO << "[--------]  rpsAngleThetaX = " << 
                 StMuRpsTrack::rpsAngleThetaX << endm;
      LOG_INFO << "[--------]  rpsAngleThetaY = " << 
                 StMuRpsTrack::rpsAngleThetaY << endm;
      LOG_INFO << "[--------]  rpsAngleTheta = " << 
                 StMuRpsTrack::rpsAngleTheta << endm;
    };
    LOG_INFO << "[++++++] ....instantiated." << endm;
};



StFmsAnalysisMaker::~StFmsAnalysisMaker(){}



Int_t StFmsAnalysisMaker::Init(){  
  mFmsDbMaker=static_cast<StFmsDbMaker*>(GetMaker("fmsDb"));  
  if(!mFmsDbMaker){
    LOG_ERROR  << "StFmsAnalysisMaker::InitRun Failed to get StFmsDbMaker" << endm;
    return kStFatal;
  }    
  mFile=new TFile(mFilename,"RECREATE");

  
  // event tree
  /*
  evtr = new TTree("evtr","evtr");
  evtr->Branch("bc",bc,"bc[2]/i"); // bunch crossing counter
  evtr->Branch("evid",&evid,"evid/I"); // event id
  evtr->Branch("nhits",&nhits,"nhits/I"); // number of hits
  evtr->Branch("nclusters",&nclusters,"nclusters/I"); // number of clusters
  evtr->Branch("npoints",&npoints,"npoints/I"); // number of points
  evtr->Branch("nTriggers",&nTriggers,"nTriggers/I"); // number of triggers fired
  evtr->Branch("trigids",trigids,"trigids[nTriggers]/I"); // list of triggers id fired
  */



  // OFile p_out tree
  p_out = new TTree("p_out","p_out");
  p_out->Branch("spin",&spin,"spin/I"); // spin of protons in bXing (0=--, 1=-+, 2=+-, 3=++ (notation (B)(Y)), 40=otherwise)
  p_out->Branch("nphotons",&nphotons,"nphotons/I"); //  number of photons found by yiqun
  p_out->Branch("br_nwrds",&nwrds,"nwrds/I"); // number of "words"
                                                // incremented by 1 for each NSTB which passes "esum-cut"
                                                // incremented also by 1 for each photon found 
                                                // finally, multiplied by 4 after all increments above have occurred for all 4 NSTBs
  p_out->Branch("br_types",tpes,"tpes[nwrds]/I");   // "types array" -- the type of int[4]-length "word" (either new NSTB or new photon... see below)
                                                    // the way the code loops below after the "event-cut" is that it loops through the 4 NSTBs, and then
                                                    // through the photons found after the reconstruction
                                                    // tpes is filled with "5+NSTB" at entry nwrds*4, then nwrds is incremented by 1
                                                    // tpes is then filled with "305+NSTB" at entry nwrds+4 for each new photon, then nwrds is incremented by 1
                                                    // all other entries are 0
                                                    // example -- say 2 photons in nstb0 and 1 in nstb1, zero in the others --> nwrds = [ (4 nstbs) + (3 total photons) ] * 4 = 7*4 = 28
                                                    // index   entry
                                                    // 0       5 + 0 = 5
                                                    // 1       0
                                                    // 2       0
                                                    // 3       0
                                                    // 4       305 + 0 = 305
                                                    // 5       0
                                                    // 6       0
                                                    // 7       0
                                                    // 8       305 + 0 = 305 
                                                    // 9       0
                                                    // 10      0
                                                    // 11      0
                                                    // 12      5 + 1 = 6
                                                    // 13      0
                                                    // 14      0
                                                    // 15      0
                                                    // 16      305 + 1 = 306
                                                    // 17      0
                                                    // 18      0
                                                    // 19      0
                                                    // 20      5 + 2 = 7
                                                    // 21      0
                                                    // 22      0
                                                    // 23      0
                                                    // 24      5 + 3 = 8
                                                    // 25      0
                                                    // 26      0
                                                    // 27      0
  p_out->Branch("br_pxyzt",pxyzt,"pxyzt[nwrds]/F"); // momentum 4-vectors for each word
                                                      // if type<300, then the 4-vector is total 4-momentum of all of the hits in the NSTB  (see AnalTools::FourMom)
                                                      // if type>300, then the 4-vector is the 4-momentum  of the photon
                                                      // the reason the types are filled every 4 entries above is because we store 4-vector info in this array: 
                                                      // example -- say 2 photons in nstb0 and 1 in nstb1, zero in the others --> nwrds = [ (4 nstbs) + (3 total photons) ] * 4 = 7*4 = 28
                                                      // index   entry
                                                      // 0       p_x of NSTB 0
                                                      // 1       p_y of NSTB 0
                                                      // 2       p_z of NSTB 0
                                                      // 3       Esum of NSTB 0
                                                      // 4       p_x of 1st photon in NSTB 0
                                                      // 5       p_y of 1st photon in NSTB 0
                                                      // 6       p_z of 1st photon in NSTB 0
                                                      // 7       E of 1st photon in NSTB 0
                                                      // 8       p_x of 2nd photon in NSTB 0
                                                      // 9       p_y of 2nd photon in NSTB 0
                                                      // 10      p_z of 2nd photon in NSTB 0
                                                      // 11      E of 2nd photon in NSTB 0
                                                      // 12      p_x of NSTB 1
                                                      // 13      p_y of NSTB 1
                                                      // 14      p_z of NSTB 1
                                                      // 15      E of NSTB 1
                                                      // 16      p_x of 1st photon in NSTB 1
                                                      // 17      p_y of 1st photon in NSTB 1
                                                      // 18      p_z of 1st photon in NSTB 1
                                                      // 19      E of 1st photon in NSTB 1
                                                      // 20      p_x of NSTB 2
                                                      // 21      p_y of NSTB 2
                                                      // 22      p_z of NSTB 2
                                                      // 23      E of NSTB 2
                                                      // 24      p_x of NSTB 3
                                                      // 25      p_y of NSTB 3
                                                      // 26      p_z of NSTB 3
                                                      // 27      E of NSTB 3
  p_out->Branch("br_Rnum",&Rnum,"Rnum/I"); // run number 
  p_out->Branch("br_Bunchid7bit",&Bunchid7bit,"Bunchid7bit/I"); // bXing number
  p_out->Branch("br_BBcSums",BBcSums,"BBcSums[5]/b"); // not needed for TwoTr creation; just passed along
  p_out->Branch("br_BBcVertex",BBcVertex,"BBcVertex[7]/F"); // not needed for TwoTr creation; just passed along


  p_out->Branch("br_EventN",&EventN,"EventN/I"); // deprecated ?; comes from dataSet::event, which is from branch "event" from much older version h111 trees
  p_out->Branch("br_ievt",&ievt,"ievt/I"); // event number (assuming each segment has 10,000 events
  p_out->Branch("br_nSavedHits",&nSavedHits,"nSavedHits/I"); // number of saved hits, see AnalTools::storeCluster, which calls AnalTools::SaveClHit, which increments this number
                                                               // for call to AnalTools::storeCluser, search for "STORE-CLUSTER" --> minimum ADC = 1; min cluster energy = 0.5 
  p_out->Branch("br_SavedHits",SavedHits,"SavedHits[nSavedHits]/I"); // encoding of a single hit: encodes east/west, nstb, row, column, and ADC
                                                                       // see AnalTools::storeCluster which calls AnalTools::SaveClHit, which does this encoding for each hit tower 
                                                                       //
                                                                       // the 4-byte encoding is layed out in the following way:
                                                                       //  - E = EW
                                                                       //  - N = NSTB
                                                                       //  - R = Row
                                                                       //  - C = Col
                                                                       //  - A = ADC
                                                                       //
                                                                       //   |  E N N N  |  R R R R  |  R R C C  |  C C C C  |
                                                                       //   |  A A A A  |  A A A A  |  A A A A  |  A A A A  | 
                                                                       //  
                                                                       // decode it by using the following (decimal places put in hex numbers for readability!)
                                                                       // ( based on AnalTools::PrClHit )
                                                                       // --> EW =   ( SavedHits & 0x.8000.0000 ) / 0x.8000.0000
                                                                       // --> NSTB = ( SavedHits & 0x.7000.0000 ) / 0x.1000.0000
                                                                       // --> Row  = ( SavedHits & 0x.0F90.0000 ) / 0x.0040.0000
                                                                       // --> Col =  ( SavedHits & 0x.003F.0000 ) / 0x.0001.0000
                                                                       // --> ADC =    SavedHits & 0x.0000.FFFF 
                                                                       //
                                                                       // ARE THESE HITS ORDERED?? --> orderd by highest-energy cluster
  p_out->Branch("br_TrigBits",&TrigBits,"TrigBits/I");    // FMS layer 2 (DSMs) output bits (i.e., input bits to TCU, known as lastdsm[5] in the h111 tree)
                                                       // in fact this is exactly lastdsm[5] in this tree
                                                       // to filter, use (TrigBits >> (L2_output_bit)) & 0x1
                                                       // 2015 pp output bits:
                                                       // 0 - smbs3
                                                       // 1 - smbs2
                                                       // 2 - smbs1
                                                       // 3 - lgbs3
                                                       // 4 - lgbs2
                                                       // 5 - lgbs1
                                                       // 6 - dibs
                                                       // 7 - jp2
                                                       // 8 - jp1
                                                       // 9 - jp0
                                                       // 10 - dijet
                                                       // 11 - unused
                                                       // 12 - unused
                                                       // 13 - unused
                                                       // 14 - unused
                                                       // 15 - unused

  p_out->Branch("br_nCluster",&nCluster,"nCluster/I"); // number of clusters; incremented every time AnalTools::storeCluster is called
  p_out->Branch("br_nPhotonClu",&nPhotonClu,"nPhotonClu/I"); // total number of photons within clusters (equivalent to StMuFmsPoints?)
                                                               // n.b. nPhotonClu-nCluster >= 0 for all events
                                                               // this number is only incremented as we loop through the clusters and look at
                                                               // their photon contents
                                                               // BUT... in a sample OFile, I found nPhotonClu==nphotons for all events
  p_out->Branch("br_SCIndex",SCIndex,"SCIndex[nCluster]/I"); // number of saved hits (nSavedHits) for current cluster, but... 
                                                                        // nSavedHits is incremented as we loop through each cluster, so basically:
                                                                        // to find out how many hits are in cluster number 3, we do: SCIndex[3]-SCIndex[2]
                                                                        // to find out how many hits are in the first cluster, we do: SCIndex[1]-SCindex[0] (n.b. SCindex[0]==0 for all events)
                                                                        // to find out how many hits are in the last cluster, we do: nSavedHits-SCindex[nCluster-1]
  p_out->Branch("br_SPCIndex",SPCIndex,"SPCIndex[nPhotonClu]/I"); // cluster index (starting at 0) of nth photon (ordered by highest-energy clusters)
  p_out->Branch("br_SPCEnergy",SPCEnergy,"SPCEnergy[nPhotonClu]/F"); // photon energy
  p_out->Branch("br_L2sum",L2sum,"L2sum[2]/i"); // essentially the trigger ID
  p_out->Branch("br_lastdsm",lastdsm,"lastdsm[8]/i"); // full TCU input bitstream
  //p_out->Branch("br_Fpde",Fpde,"Fpde[8]/i"); // unecessary for run 15++ ???


  // Roman Pot Branches
  p_out->Branch("RP_n_tracks",&n_tracks,"RP_n_tracks/I");
  p_out->Branch("RP_n_trackpoints",&n_trackpoints,"RP_n_trackpoints/I");

  // tracks
  p_out->Branch("RP_t_index",t_index,"RP_t_index[RP_n_tracks]/I"); // track index number (my definition, since RP tracks are sorted differently)
  p_out->Branch("RP_t_branch",t_branch,"RP_t_branch[RP_n_tracks]/I"); // RP branch (0=EU 1=ED 2=WU 3=WD)
  p_out->Branch("RP_t_type",t_type,"RP_t_type[RP_n_tracks]/I"); // track type (see next line)
    /* 0=rpsLocal -- 1 track point
     * 1=rpsGlobal -- 2 track points
     * 2=rpsUndefined -- track not defined
     */
  p_out->Branch("RP_t_planesUsed",t_planesUsed,"RP_t_planesUsed[RP_n_tracks]/i"); // number of SSD planes hit by track points in track
  p_out->Branch("RP_t_p",t_p,"RP_t_p[RP_n_tracks]/D"); // momentum
  p_out->Branch("RP_t_pt",t_pt,"RP_t_pt[RP_n_tracks]/D"); // transverse momentum
  p_out->Branch("RP_t_eta",t_eta,"RP_t_eta[RP_n_tracks]/D"); // pseudorapidity
  p_out->Branch("RP_t_time",t_time,"RP_t_time[RP_n_tracks]/D"); // time of track detection

  if(n_ang != StMuRpsTrack::mNumberOfAngleTypes) {
    LOG_ERROR << "StFmsAnalysisMaker::Init -- n_ang not equal to StMuRpsTrack::mNumberOfAngleTypes" << endm;
    return kStErr;
  };
  TString t_theta_str = Form("RP_t_theta[RP_n_tracks][%d]/D",StMuRpsTrack::mNumberOfAngleTypes);
  TString t_thetaRP_str = Form("RP_t_thetaRP[RP_n_tracks][%d]/D",StMuRpsTrack::mNumberOfAngleTypes);
  p_out->Branch("RP_t_theta",t_theta,t_theta_str.Data()); // polar angle at RP according to STAR coord sys
  p_out->Branch("RP_t_thetaRP",t_thetaRP,t_thetaRP_str.Data()); // polar angle at RP according to STAR survey

  p_out->Branch("RP_t_phi",t_phi,"RP_t_phi[RP_n_tracks]/D"); // azimuth at RP according to STAR coord sys
  p_out->Branch("RP_t_phiRP",t_phiRP,"RP_t_phiRP[RP_n_tracks]/D"); // azimuth at RP according to STAR survey
  p_out->Branch("RP_t_t",t_t,"RP_t_t[RP_n_tracks]/D"); // squared 4-momentum transfer
  p_out->Branch("RP_t_xi",t_xi,"RP_t_xi[RP_n_tracks]/D"); // fractional momentum loss (pbeam-p)/pbeam
  p_out->Branch("RP_t_gold",t_gold,"RP_t_gold[RP_n_tracks]/O"); // my track quality variable (2 track points in all 2x4=8 Si planes)
  p_out->Branch("RP_t_isBad",t_isBad,"RP_t_isBad[RP_n_tracks]/O"); // true if any kinematics have nonsense values
  p_out->Branch("RP_t_qualHash",t_qualHash,"RP_t_qualHash[RP_n_tracks]/D"); // hash track info into a "quality" number, for sorting tracks

  // track point 0
  p_out->Branch("RP_p0_tpExists",p_tpExists[0],"RP_p0_tpExists[RP_n_tracks]/O"); // true if track point 0 exists
  p_out->Branch("RP_p0_RPid",p_RPid[0],"RP_p0_RPid[RP_n_tracks]/I"); // RP id (see next line)
    /* 0=E1U  1=E1D  2=E2U  3=E2D
     * 4=W1U  5=W1D  6=W2U  7=W2D
     */
  p_out->Branch("RP_p0_clustid_s1",p_clustid_s1[0],"RP_p0_clustid_s1[RP_n_tracks]/I"); // track point cluster IDs, indexed by SSD plane (1-4)
                                                                                       // =-1 if plane not used in track point
  p_out->Branch("RP_p0_clustid_s2",p_clustid_s2[0],"RP_p0_clustid_s2[RP_n_tracks]/I");
  p_out->Branch("RP_p0_clustid_s3",p_clustid_s3[0],"RP_p0_clustid_s3[RP_n_tracks]/I");
  p_out->Branch("RP_p0_clustid_s4",p_clustid_s4[0],"RP_p0_clustid_s4[RP_n_tracks]/I");
  p_out->Branch("RP_p0_quality",p_quality[0],"RP_p0_quality[RP_n_tracks]/I"); // track point quality (see next line)
    /* 0=rpsNormal -- not golden and not undefined
     * 1=rpsGolden -- single cluster in all 4 SSD planes
     * 2=rpsNotSet -- undefined track point
     */
  p_out->Branch("RP_p0_planesUsed",p_planesUsed[0],"RP_p0_planesUsed[RP_n_tracks]/i"); // number of SSD planes with valid cluster for track point
  p_out->Branch("RP_p0_x",p_x[0],"RP_p0_x[RP_n_tracks]/D"); // STAR survey coords x-position
  p_out->Branch("RP_p0_y",p_y[0],"RP_p0_y[RP_n_tracks]/D"); // STAR survey coords y-position
  p_out->Branch("RP_p0_z",p_z[0],"RP_p0_z[RP_n_tracks]/D"); // STAR survey coords z-position
  p_out->Branch("RP_p0_time_pmt1",p_time_pmt1[0],"RP_p0_time_pmt1[RP_n_tracks]/D"); // TAC (in ns) of PMT 1 hit (2 PMT's per vessel)
  p_out->Branch("RP_p0_time_pmt2",p_time_pmt2[0],"RP_p0_time_pmt2[RP_n_tracks]/D"); // TAC (in ns) of PMT 1 hit (2 PMT's per vessel)
  
  // track point 1
  p_out->Branch("RP_p1_tpExists",p_tpExists[1],"RP_p1_tpExists[RP_n_tracks]/O");
  p_out->Branch("RP_p1_RPid",p_RPid[1],"RP_p1_RPid[RP_n_tracks]/I");
  p_out->Branch("RP_p1_clustid_s1",p_clustid_s1[1],"RP_p1_clustid_s1[RP_n_tracks]/I");
  p_out->Branch("RP_p1_clustid_s2",p_clustid_s2[1],"RP_p1_clustid_s2[RP_n_tracks]/I");
  p_out->Branch("RP_p1_clustid_s3",p_clustid_s3[1],"RP_p1_clustid_s3[RP_n_tracks]/I");
  p_out->Branch("RP_p1_clustid_s4",p_clustid_s4[1],"RP_p1_clustid_s4[RP_n_tracks]/I");
  p_out->Branch("RP_p1_quality",p_quality[1],"RP_p1_quality[RP_n_tracks]/I");
  p_out->Branch("RP_p1_planesUsed",p_planesUsed[1],"RP_p1_planesUsed[RP_n_tracks]/i");
  p_out->Branch("RP_p1_x",p_x[1],"RP_p1_x[RP_n_tracks]/D");
  p_out->Branch("RP_p1_y",p_y[1],"RP_p1_y[RP_n_tracks]/D");
  p_out->Branch("RP_p1_z",p_z[1],"RP_p1_z[RP_n_tracks]/D");
  p_out->Branch("RP_p1_time_pmt1",p_time_pmt1[1],"RP_p1_time_pmt1[RP_n_tracks]/D");
  p_out->Branch("RP_p1_time_pmt2",p_time_pmt2[1],"RP_p1_time_pmt2[RP_n_tracks]/D");


  Rnum_tmp = 0;

  return kStOK;
}

Int_t StFmsAnalysisMaker::Finish(){
  LOG_INFO << Form("Writing and closing %s",mFilename) << endm;
  mFile->Write();
  mFile->Close();
  return kStOK;
}

Int_t StFmsAnalysisMaker::Make(){
  // obtain event instance
  StMuEvent * muEvent = mMuDstMaker->muDst()->event();
  if(!muEvent) {
    LOG_ERROR << "StFmsAnalysisMaker::Make did not find StMuEvent" << endm;
    return kStErr;
  };


  // obtain collection and arrays
  mMuFmsColl = mMuDstMaker->muDst()->muFmsCollection();
  if(!mMuFmsColl) {LOG_ERROR << "StFmsAnalysisMaker::Make did not find MuFmsCollection in StMuEvent"<<endm; return kStErr;}
  mMuRpsColl = mMuDstMaker->muDst()->RpsCollection();
  if(!mMuRpsColl) {LOG_ERROR << "StFmsAnalysisMaker::Make did not find MuRpsCollection in StMuEvent"<<endm; return kStErr;}
  StRunInfo runInfo = muEvent->runInfo();
  
  
  // get FMS object arrays
  mHits = mMuFmsColl->getHitArray();
  mClustersAll = mMuFmsColl->getClusterArray();
  mPoints = mMuFmsColl->getPointArray();

  nhits = mMuFmsColl->numberOfHits();
  nclusters = mMuFmsColl->numberOfClusters();
  npoints = mMuFmsColl->numberOfPoints();

  if(verbose) LOG_INFO << "[+] nhits=" << nhits << "  nclusters=" << nclusters << "  npoints=" << npoints << endm;


  // event and run numbers
  evid = muEvent->eventId();
  Rnum = muEvent->runNumber(); // --> p_out
  ievt = evid; // --> p_out
  EventN = 0; // --> p_out  // unused in OFile->Output chain?
  nTriggers = muEvent->triggerIdCollection().nominal().triggerIds().size();
  if(verbose) LOG_INFO << " >>>>>>>>>>>>>>>> nTriggers=" << nTriggers << endm;

  for(int u=0; u<nTriggers; ++u) {
    trigids[u] = (muEvent->triggerIdCollection().nominal().triggerIds())[u];
    if(verbose) LOG_INFO << "                triggerid " << trigids[u] << " fired" << endm;
  };


  // evtr->Fill();


  // get lastdsm (TCU inputs) and L2sum (trigger IDs satisfied)
  for(int dd=0; dd<8; dd++) 
    lastdsm[dd] = muEvent->triggerData()->lastDSM(dd); // --> p_out
  TrigBits = lastdsm[5]; // --> p_out

  L2sum_full = muEvent->triggerData()->l2sum();
  for(int ld=0; ld<2; ld++)
    L2sum[ld] = (L2sum_full >> ld*32) & 0xFFFFFFFF; // --> p_out

  // bXing number
  Bunchid7bit = muEvent->triggerData()->bunchId7Bit(); // --> p_out
  for(int b=0; b<2; b++) bc[b]=(muEvent->eventInfo()).bunchCrossingNumber(b);




  // obtain spinbit from spindb 
  if(Rnum_tmp!=Rnum) {
    mSpinDb->InitRun(Rnum);
    //assert(mSpinDb->isValid() && "mSpinDb->isValid = false");
    //mSpinDb->print(1);
    Rnum_tmp = Rnum;
  };
  spinbyte_db = (mSpinDb->getSpin8bits())[Bunchid7bit];
  spinbit = mSpinDb->spin4usingBX7(Bunchid7bit);

  // obtain spinbyte from v124bits from spindbmaker
  if(verbose) LOG_INFO << "[+] Bunchid7bit=" << Bunchid7bit << endm;
  yNib = (mSpinDb->getRawV124bits())[ (Bunchid7bit*3 + 240) % 360 ] & 0x0f;
  bNib = (mSpinDb->getRawV124bits())[ Bunchid7bit*3 ] >> 4;
  spinbyte = (bNib<<4) + yNib;


  if( (spinbyte >> 0 ) & 0x1) {
    if( (spinbyte >> 2) & 0x1) sp_y=1;
    else if( (spinbyte >> 1) & 0x1) sp_y=0;
  }
  else sp_y=-1;

  if( (spinbyte >> 4 ) & 0x1) {
    if( (spinbyte >> 6) & 0x1) sp_b=1;
    else if( (spinbyte >> 5) & 0x1) sp_b=0;
  }
  else sp_b=-1;

  spin = (sp_y>=0 && sp_b>=0) ? (sp_b<<1)+sp_y : 40; // --> p_out


  if(verbose) {
    LOG_INFO << "[+] spinbyte=" << std::bitset<8>(spinbyte)
             << "  spinbyte_db=" << std::bitset<8>(spinbyte_db)
             << "  spinbit=" << spinbit
             << "  spin=" << spin <<  endm;
  };


  // dump spinbytes (only for first event)
  if(dump_spinbits) {
    LOG_INFO << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS" << endm;
    LOG_INFO << "--- BEGIN SPINBYTE DUMP ---" << endm;
    for(int bx=0; bx<120; bx++) {
      yNib = (mSpinDb->getRawV124bits())[ (bx*3 + 240) % 360 ] & 0x0f;
      bNib = (mSpinDb->getRawV124bits())[ bx*3 ] >> 4;
      spinbyte_curr = (bNib<<4) + yNib;
      spinbyte_db_curr = (mSpinDb->getSpin8bits())[bx];
      spinbyte_same = spinbyte_curr==spinbyte_db_curr;


      if( (spinbyte_curr >> 0 ) & 0x1) {
        if( (spinbyte_curr >> 2) & 0x1) spc_y='+';
        else if( (spinbyte_curr >> 1) & 0x1) spc_y='-';
      }
      else spc_y='0';

      if( (spinbyte_curr >> 4 ) & 0x1) {
        if( (spinbyte_curr >> 6) & 0x1) spc_b='+';
        else if( (spinbyte_curr >> 5) & 0x1) spc_b='-';
      }
      else spc_b='0';

      LOG_INFO << std::bitset<4>(spinbyte_curr>>4) << 
                  "." << 
                  std::bitset<4>(spinbyte_curr) << 
                  "  " << 
                  std::bitset<4>(spinbyte_db_curr>>4) << 
                  "." << 
                  std::bitset<4>(spinbyte_db_curr) << 
                  "  " << 
                  "same=" << spinbyte_same << " " <<  
                  spc_b << " " << spc_y << " bx=" << bx << endm;
    };
    LOG_INFO << "--- END SPINBYTE DUMP ---" << endm;
    LOG_INFO << "SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS" << endm;
    dump_spinbits = false;
  }; // eo if(dump_spinbits)



  // BBc* variables for p_out
  for(int bb=0; bb<5; bb++) BBcSums[bb]=0; // --> p_out
  for(int bb=0; bb<7; bb++) BBcVertex[bb]=0; // --> p_out


  // FPDE variable
  for(int bb=0; bb<8; bb++) 
    Fpde[bb] = muEvent->triggerData()->fpdLayer1DSMRaw(StBeamDirection::west,bb); // --> p_out
 
  
  // debugging diagnostics for trigger info
  if(verbose) {
    LOG_INFO << Form("TRIGGER -- bx=%d  |  spinbyte=0x%X=",
      Bunchid7bit,spinbyte) << std::bitset<8>(spinbyte) << endm;
    //LOG_INFO << Form("TRIGGER -- Bunchid7bit=%d  spinbit=%d  bc[1]=%d bc[0]=%d  (b7bittest=%d)",
      //Bunchid7bit,spinbit,bc[1],bc[0],((bc[0]-7)%120)==Bunchid7bit) << endm;
    //for(int dd=0; dd<8; dd++) LOG_INFO << Form("           lastdsm[%d]=0x%x",dd,lastdsm[dd]) << endm;
    LOG_INFO << Form("           TrigBits=0x%x",TrigBits) << endm;
    for(int ld=0; ld<2; ld++) LOG_INFO << Form("           l2sum[%d]=0x%x",ld,L2sum[ld]) << endm;
    LOG_INFO << Form("           l2sum_full=0x%llx",L2sum_full) << endm;
  };



  //////////////////////////////////////////////////////////////////////////////////
  //
  //    FMS Sector 
  //
  //////////////////////////////////////////////////////////////////////////////////
  

  
  // reorganize clusters array into arrays for each NSTB and compute cluster sums for each NSTB
  en_sum_all=0;
  for(int nn=0; nn<4; nn++) {
    mClusters[nn] = new TObjArray();
    en_sum[nn]=0;
  };

  for(int ci=0; ci<nclusters; ci++) {
    if(mMuFmsColl->getCluster(ci)) { 
      nstb = mMuFmsColl->getCluster(ci)->detectorId() - 7; // NSTB = detectorId-7 (see StEnumeerations.h)
      if(nstb>=1 && nstb<=4) {
        mClusters[nstb-1]->AddLast(mMuFmsColl->getCluster(ci));
        clu_en = mMuFmsColl->getCluster(ci)->energy();
        en_sum[nstb-1] += clu_en;
        en_sum_all += clu_en;
      }
      else LOG_ERROR << "NSTB of cluster out of range; IGNORING THIS CLUSTER" << endm;
    };
  };

  for(int nn=0; nn<4; nn++) ent[nn] = mClusters[nn]->GetEntries();
  if(verbose) {
    LOG_INFO << "CLUSTER LOOP --- " << "nclusters=" << nclusters
             << " n(1)=" << ent[0] << " n(2)=" << ent[1] << " n(3)=" << ent[2] << " n(4)=" << ent[3] << endm;
  };


  // fill energy matrix
  for(int nn=0; nn<4; nn++) {
    for(int nr=0; nr<nRows[nn]; nr++) {
      for(int nc=0; nc<nCols[nn]; nc++) {
        (*(eMat[nn]))(nr,nc)=0.0;
      };
    };
  };
  for(int hi=0; hi<nhits; hi++) {
    if(mMuFmsColl->getHit(hi)) {
      nstb = mMuFmsColl->getHit(hi)->detectorId() - 7;
      chan = mMuFmsColl->getHit(hi)->channel();
      hit_en = mMuFmsColl->getHit(hi)->energy();

      if(verbose) LOG_INFO << "[++++++++++++] nstb=" << nstb << endm;


      if(nstb>=1 && nstb<=4 && hit_en>0) {

        hit_r = (chan-1)/nCols[nstb-1];
        hit_c = (chan-1)%nCols[nstb-1];

        if(verbose) {
          LOG_INFO << " ENERGY MATRIX LOOP :: hi=" << hi
                   << " nstb=" << nstb
                   << " chan=" << chan
                   << " (r" << hit_r << ",c" << hit_c << ")"
                   << " en=" << hit_en
                   << " adc=" << mMuFmsColl->getHit(hi)->adc()
                   << endm;
        };

        if(hit_r>=0 && hit_r<nRows[nstb-1] && hit_c>=0 && hit_c<nCols[nstb-1]) {
          (*(eMat[nstb-1]))(hit_r,hit_c) = hit_en;
        };
      };
    };
  };
  

  // reset cluster & hit counters & arrays for p_out
  nwrds = 0; 
  nSavedHits = 0;
  nCluster = 0;
  nPhotonClu = 0;
  nphotons = 0;
  num_null = 0;
  for(int ii=0; ii<MAX_POUT; ii++) {
    tpes[ii] = 0;
    pxyzt[ii] = 0.0;
    SavedHits[ii] = 0;
    SCIndex[ii] = 0;
    SPCIndex[ii] = 0;
    SPCEnergy[ii] = 0.0;
  };


  // BEGIN MAIN CLUSTER LOOP
  // -----------------------
  // loop through each nstb, in turn, looping through their clusters
  for(int nn=0; nn<4; nn++) {
    nstb=nn+1;
    if(verbose_clu) LOG_INFO << "NSTB " << nn+1 << "-----------------------------------" << endm;
    for(int ii=0; ii<MAX_POUT; ii++) clu_en_sortval[ii]=-1; // reset sortval array

    // "esum-cut" for each NSTB
    if(en_sum[nn]>1.0 && en_sum[nn]<700.0) {
      
      // compute "average" momentum 4-vector for this nstb
      // (following AnalTools->FourMom->mom(TMatrix,...) )
      TVector3 off(0,0,0);
      TVector3 momsum(0.,0.,0.);
      TVector3 xyz(0.,0.,0.);
      Float_t width = FpdTowWid[1][nn];
      off[2] = ZFPD[1][nn];
      off[1] = yOffset[1][nn];
      off[0] = xOffset[1][nn];
      Int_t signx = (nn==0||nn==2) ? 1:-1;
      EnergySum = 0;
      for(int ir=0; ir<nRows[nn]; ir++) {
        for(int ic=0; ic<nCols[nn]; ic++) {
          xyz[0] = off[0] - signx*ic*width;
          xyz[1] = off[1] - ir*width;
          xyz[2] = off[2];
          xyz = xyz * (1/xyz.Mag());
          Energy = TMath::Max((*eMat[nn])(ir,ic),(Float_t)0.0);
          xyz = xyz * Energy;
          EnergySum += Energy; // 4-momentum energy
          momsum = momsum + xyz; // 4-momentum momenta
        };
      };


      // add NSTB blocks to types and 4-momenta arrays
      tpes[nwrds*4] = nn+5; // --> p_out
      vecs[nwrds] = StLorentzVectorF(momsum[0],momsum[1],momsum[2],EnergySum);
      nwrds++;


      // sort clusters by energy
      if(verbose_clu) LOG_INFO << "     unsorted clusters:" << endm;
      for(int ci=0; ci<ent[nn]; ci++)
      {
        clu_en_sortval[ci] = ((StMuFmsCluster*)(mClusters[nn]->At(ci)))->energy();
        if(verbose_clu)  LOG_INFO << "          en=" << clu_en_sortval[ci] << endm;
      }; 
      TMath::Sort(MAX_POUT,clu_en_sortval,clu_en_sortind);
      

      // sorted cluster loop
      if(verbose_clu) LOG_INFO << "     sorted clusters:" << endm;
      for(int ci=0; ci<ent[nn]; ci++) {
        clu = (StMuFmsCluster*)(mClusters[nn]->At(clu_en_sortind[ci]));
        clu_en = clu->energy();
        clu_x = clu->x();
        clu_y = clu->y();
        clu_sigmin = clu->sigmaMin();
        clu_sigmax = clu->sigmaMax();
        clu_csqn1 = clu->chi2Ndf1Photon();
        clu_csqn2 = clu->chi2Ndf2Photon();
        clu_id = clu->id();
        if(verbose_clu) {
          LOG_INFO << "          id=" << clu_id
                   << " en=" << clu_en
                   << " (x=" << clu_x << ",y=" << clu_y << ")"
                   << " nph=" << clu->photons()->GetEntries()
                   << endm;
        };

        // loop through this cluster's photons
        for(int pi=0; pi<clu->photons()->GetEntries(); pi++) {
          pho = (StMuFmsPoint*)(clu->photons()->At(pi));
          if(pho) {

            vecs[nwrds] = pho->fourMomentum(); // NOTE: default arguments means fms at 735.8 cm from IP; could later put in actual zvertex correction...
            tpes[nwrds*4] = 305 + nn; // --> p_out

            pho_id = pho->id();
            pho_en = pho->energy();
            pho_x = pho->x();
            pho_y = pho->y();

            SPCIndex[nPhotonClu] = nCluster; // --> p_out
            SPCEnergy[nPhotonClu] = pho_en; // --> p_out
            nPhotonClu++;

            if(verbose_pho) {
              LOG_INFO << "            photon " << pi+1
                       << " id=" << pho_id
                       << " en=" << pho_en
                       << " (x=" << pho_x << ",y=" << pho_y << ")"
                       << " 4-mom=(" 
                          << vecs[nwrds].px() << ","
                          << vecs[nwrds].py() << ","
                          << vecs[nwrds].pz() << ","
                          << vecs[nwrds].e() << ")"
                       << "  parent_clu_id=" << pho->cluster()->id()
                       << " xyz=("
                          << (pho->xyz()).x() << ","
                          << (pho->xyz()).y() << ","
                          << (pho->xyz()).z() << ")"
                       << endm;
              LOG_INFO << " ---- SPCIndex[" << nPhotonClu-1
                       << "] = " << SPCIndex[nPhotonClu-1] << endm;
              LOG_INFO << " ---- SPCEnergy[" << nPhotonClu-1
                       << "] = " << SPCEnergy[nPhotonClu-1] << endm;
            };

            nwrds++;
            nphotons++; // --> p_out

          }
          else {
            num_null++;
            if(verbose_pho) LOG_INFO << "            photon " << pi+1 << " -- NULL POINTER... potential issue??" << endm;
          };

        }; // eo cluster's photon loop


        // append to SCIndex
        SCIndex[nCluster] = nSavedHits; // --> p_out
        if(verbose_pho) LOG_INFO << "      SCIndex[" << nCluster << "]=" << SCIndex[nCluster] << endm;
        nCluster++; // --> p_out


        // loop through this cluster's hits
        for(int hi=0; hi<clu->hits()->GetEntries(); hi++) {
          hit = (StMuFmsHit*)(clu->hits()->At(hi));
          if(hit) {
            hit_n = hit->detectorId()-7;
            hit_r = (hit->channel()-1) / nCols[hit_n-1];
            hit_c = (hit->channel()-1) % nCols[hit_n-1];
            hit_adc = hit->adc();

            p1 =                    0x80000000;
            p2 = ((hit_n-1)&7)    * 0x10000000;
            p3 = ((hit_r-1)&0x3F) * 0x00400000;
            p4 = ((hit_c-1)&0x3F) * 0x00010000;

            SavedHits[nSavedHits] = (hit_adc&0xFFFF) + p1 + p2 + p3 + p4; // --> p_out
            nSavedHits++; // --> p_out

            if(verbose_pho) {
              LOG_INFO << " SavedHits loop: (n" << hit_n << ",r" << hit_r << ",c" << hit_c << ")  adc="
                       << hit_adc << " SavedHits[" << nSavedHits-1 << "]="
                       << Form("0x%8x",SavedHits[nSavedHits-1])
                       << endm;
            };
          }
          else {
            if(verbose_pho) LOG_INFO << "            hit " << hi+1 << " -- NULL POINTER... potential issue??" << endm;
          };
        }; // eo cluster's hit loop
        



      }; // eo sorted cluster loop
    }; // eo if ent[nn]>0
  }; // eo NSTB loop
 
  
  // fill 4-momentum array
  int wcnt = 0;
  for(int kk=0; kk<nwrds; kk++) {
    pxyzt[wcnt++] = vecs[kk].px(); // --> p_out
    pxyzt[wcnt++] = vecs[kk].py(); // --> p_out
    pxyzt[wcnt++] = vecs[kk].pz(); // --> p_out
    pxyzt[wcnt++] = vecs[kk].e(); // --> p_out
  };

  // nwrds = (3+1 dimensions) * (number of valid photons + number of NSTBs which pass esum-cut)
  nwrds *= 4; // --> p_out

  


  //////////////////////////////////////////////////////////////////////////////////
  //
  //    RP Sector 
  //
  //////////////////////////////////////////////////////////////////////////////////

  Double_t thetaRpLimits[2][2] = { 
    {-1.5e-3, 5.0e-3},
    { 1.0e-3, 5.5e-3} 
  }; // [rad]; theta limits from rafal

  // get number of tracks/trackpoints
  n_trackpoints = 0;
  n_tracks = 0;
  while(mMuRpsColl->trackPoint(n_trackpoints)) ++n_trackpoints;
  while(mMuRpsColl->track(n_tracks)) ++n_tracks;


  // reset track and track point variables for this event
  for(int i=0; i<n_tracks_max; i++) ResetBranches(i);
  

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  // pre-track loop, for sorting tracks 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(n_tracks>n_tracks_max) {
    LOG_INFO << "WARNING WARNING WARNING -- n_tracks > n_tracks_max; only analysing first " 
            << n_tracks_max << " tracks (also will set n_trackpoints to -1)" << endm;
    n_tracks = n_tracks_max;
    n_trackpoints = -1;
  }
  if(n_trackpoints>2*n_tracks_max) {
    LOG_ERROR << "WARNING WARNING WARNING -- n_trackpoints > 2 * n_tracks_max; setting to -1" << endm;
    n_trackpoints = -1;
  };

  if(n_tracks>0) {
    for(int i = 0; i<n_tracks; i++) {
      trk = mMuRpsColl->track(i);
      trkPlanes = trk->planesUsed();
      trkMom = trk->p();
      trkThetaX = trk->thetaRp(StMuRpsTrack::rpsAngleThetaX);
      trkThetaY = trk->thetaRp(StMuRpsTrack::rpsAngleThetaY);
      trkGlobal = (trk->type()==StMuRpsTrack::rpsGlobal) ? 1:0;
      if(trkThetaX > thetaRpLimits[0][0] && trkThetaX < thetaRpLimits[0][1] &&
         TMath::Abs(trkThetaY) > thetaRpLimits[1][0] && TMath::Abs(trkThetaY) < thetaRpLimits[1][1])
        trkInGeom = 1;
      else trkInGeom=0;

      // bound momentum in order to get reasonable track hash for sorting
      if(trkMom<0) {
        LOG_ERROR << "WARNING WARNING WARNING -- trkMom<0" << endm;
        trkMom=0;
      }

      // suppress tracks with nonsensical momenta (c.f. stronger cut later, where t_isBad is set)
      if(trkMom>300 || trkMom<0) {
        trkMom=0;
        trkPlanes=0;
        trkGlobal=0;
        trkInGeom=0;
      };

      // hash track kinematics etc. into a number for sorting
      // hash = Min(trkMom,300) + 1,000*trkPlanes + 10,000*trkGlobal + 20,000*trkInGeom
      // -- trkGlobal=0 && trkInGeom=0  --> hash in range  2,000 -  8,300
      // -- trkGlobal=0 && trkInGeom=1  --> hash in range 12,000 - 18,300
      // -- trkGlobal=1 && trkInGeom=0  --> hash in range 22,000 - 28,300
      // -- trkGlobal=1 && trkInGeom=1  --> hash in range 32,000 - 38,300
      //  -- trkMom bound between 0 and 300 (to prevent bad hash of nonsensical momenta)
      //  -- nplanes in a single track point is >=2; nplanes in a single track is <=8
      trkSortQual[i] = trkMom + 
                      1000*trkPlanes + 
                      10000*trkInGeom +
                      20000*trkGlobal;
    };
    TMath::Sort(n_tracks_max,trkSortQual,trkSortIdx);

    // dump sort table (DEBUGGING)
    if(false) {
      for(int i = 0; i<n_tracks; i++) {
        LOG_INFO << std::setprecision(20) << "trkSortIdx[" << i << "]=" << trkSortIdx[i] <<
                    "  trkSortQual[" << trkSortIdx[i] << "]=" << trkSortQual[trkSortIdx[i]] <<
                    endm;
      };
    };
  };


  

  // ----------------------------------------------
  // main track loop
  // ----------------------------------------------
  if(n_tracks>0) {
    for(int i = 0; i<n_tracks; i++) {
      if(verbose_rp) LOG_INFO << "for i="<<i<<endm;

      if(trkSortIdx[i]<0) LOG_ERROR << "ERROR: negative trkSortIdx" << endm;
      if(trkSortIdx[i]>n_tracks_max) LOG_ERROR << "ERROR: trkSortIdx > n_tracks_max" << endm;

      trk = mMuRpsColl->track(trkSortIdx[i]);

      t_qualHash[i] = trkSortQual[trkSortIdx[i]];

      // fill track leaves
      t_index[i] = i;
      t_branch[i] = trk->branch();
      t_type[i] = (Int_t)(trk->type());
      t_planesUsed[i] = trk->planesUsed();
      t_p[i] = trk->p();
      t_pt[i] = trk->pt();
      t_eta[i] = trk->eta();
      t_time[i] = trk->time();
      

      for(int a=0; a<n_ang; a++) {
        t_theta[i][a] = trk->theta(a);
        t_thetaRP[i][a] = trk->thetaRp(a);
      };


      t_phi[i] = trk->phi();
      t_phiRP[i] = trk->phiRp();
      t_t[i] = trk->t( runInfo.beamEnergy(trk->branch()<2 ? StBeamDirection::east : StBeamDirection::west) );
      t_xi[i] = trk->xi( runInfo.beamEnergy(trk->branch()<2 ? StBeamDirection::east : StBeamDirection::west) );

      // fill track point branches
      for(int j=0; j<2; j++) {
        if(trk->trackPoint(j)!=NULL) {
          p_tpExists[j][i] = true;
          p_RPid[j][i] = trk->trackPoint(j)->rpId();
          p_clustid_s1[j][i] = trk->trackPoint(j)->clusterId(0);
          p_clustid_s2[j][i] = trk->trackPoint(j)->clusterId(1);
          p_clustid_s3[j][i] = trk->trackPoint(j)->clusterId(2);
          p_clustid_s4[j][i] = trk->trackPoint(j)->clusterId(3);
          p_quality[j][i] = (Int_t)(trk->trackPoint(j)->quality());
          p_planesUsed[j][i] = trk->trackPoint(j)->planesUsed();
          p_x[j][i] = trk->trackPoint(j)->x();
          p_y[j][i] = trk->trackPoint(j)->y();
          p_z[j][i] = trk->trackPoint(j)->z();
          p_time_pmt1[j][i] = trk->trackPoint(j)->time(0);
          p_time_pmt2[j][i] = trk->trackPoint(j)->time(1);
        }
        else {
          p_tpExists[j][i] = false;
          p_RPid[j][i] = -1;
          p_clustid_s1[j][i]=-1;
          p_clustid_s2[j][i]=-1;
          p_clustid_s3[j][i]=-1;
          p_clustid_s4[j][i]=-1;
          p_quality[j][i] = -1;
          p_x[j][i] = -1000;
          p_y[j][i] = -1000;
          p_z[j][i] = -1000;
          p_time_pmt1[j][i]=-1;
          p_time_pmt2[j][i]=-1;
        };
      };


      // deprecated 'golden track' quality cut
      if(t_type[i]==StMuRpsTrack::rpsGlobal && 
         p_quality[0][i]==StMuRpsTrackPoint::rpsGolden &&
         p_quality[1][i]==StMuRpsTrackPoint::rpsGolden) t_gold[i] = true;
      else t_gold[i] = false;


      // suppress possible nonsense values and mark track as t_isBad=true
      t_isBad[i]=false;
      if(t_branch[i]<0 || t_branch[i]>3) t_isBad[i]=true;
      if(t_type[i]<0 || t_type[i]>2) t_isBad[i]=true;
      if(t_planesUsed[i]<0 || t_planesUsed[i]>8) t_isBad[i]=true;
      if(t_p[i]<0 || t_p[i]>300) t_isBad[i]=true;
      if(t_pt[i]<0 || t_pt[i]>300) t_isBad[i]=true;
      if(t_eta[i]<-30 || t_eta[i]>30) t_isBad[i]=true;
      if(t_xi[i]<-3 || t_xi[i]>3) t_isBad[i]=true;
      if(t_t[i]<-1000 || t_t[i]>0) t_isBad[i]=true;

      if(t_isBad[i]) ResetBranches(i);


    }; // eo track loop
  }; // eo if n_tracks>0



  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////



  // fill p_out tree and print arrays
  if(nwrds>0) {
    p_out->Fill();
    if(verbose_pho) {
      LOG_INFO << "nSavedHits=" << nSavedHits << endm;
      LOG_INFO << "nPhotonClu=" << nPhotonClu << endm;
      LOG_INFO << "TYPES AND PXYZT ARRAYS -----------------------" << endm;
      for(int kk=0; kk<nwrds; kk++) {
        LOG_INFO << "     ent = " << kk
                 << "  types = " << tpes[kk]
                 << "  pxyzt = " << pxyzt[kk]
                 << endm;
      };
    };
  };



  // print some statistics
  if(verbose) {
    // note: diff is npoints-nphotons, where npoints is the total number of points
    // and nphotons is the number of photons added to p_out; this number can differ
    // because of the above esum-cut
    LOG_INFO << ">>>>>   "
             << " nhits=" << nhits
             << " nclusters=" << nclusters
             << " npoints=" << npoints
             << " nphotons=" << nphotons
             << " diff=" << npoints-nphotons
             << " num_null=" << num_null
             << endm;
  };


  // print list of photons (for comparing to photons in clusters)
  // -- sometimes you'll see photons in this list which aren't in the list of photons
  //    within each cluster, collected above; these are just the photons in clusters which
  //    did not pass esum-cut
  if(verbose_pho) {
    LOG_INFO << "POINT LOOP" << endm;
    for(int pi=0; pi<npoints; pi++) {
      pho = mMuFmsColl->getPoint(pi);
      if(pho) {
        LOG_INFO << " *** found photon: id=" << pho->id() 
                 << " en=" << pho->energy() 
                 << " nstb=" << pho->detectorId()-7
                 << " clu_id=" << pho->cluster()->id()
                 << endm;
      }
      else LOG_ERROR << "*** WARNING: null pointer to photon in mMuFmsColl, entry " << pi << endm;
    };
  };




  // clean up
  for(int nn=0; nn<4; nn++) if(mClusters[nn]) delete mClusters[nn];

  return kStOK;
};

void StFmsAnalysisMaker::ResetBranches(Int_t i_) {
  t_index[i_]=-1;
  t_branch[i_]=-1;
  t_type[i_]=-1;
  t_planesUsed[i_]=0;
  t_p[i_]=-1;
  t_pt[i_]=-1;
  t_eta[i_]=-100;
  t_time[i_]=-1;

  for(int a=0; a<n_ang; a++) {
    t_theta[i_][a] = -1;
    t_thetaRP[i_][a] = -1;
  };

  t_phi[i_]=-100;
  t_phiRP[i_]=-100;
  t_t[i_]=-1000;
  t_xi[i_]=-1000;
  t_gold[i_]=false;
  //t_isBad[i_]=true; // DO NOT RESET isBad HERE
  t_qualHash[i_]=-1000;

  trkSortQual[i_]=-1000;
  trkSortIdx[i_]=-1;

  for(int j=0; j<2; j++) {
    p_tpExists[j][i_]=false;
    p_RPid[j][i_]=-1;
    p_clustid_s1[j][i_]=-1;
    p_clustid_s2[j][i_]=-1;
    p_clustid_s3[j][i_]=-1;
    p_clustid_s4[j][i_]=-1;
    p_quality[j][i_]=-1;
    p_planesUsed[j][i_]=0; 
    p_x[j][i_]=-1000;
    p_y[j][i_]=-1000;
    p_z[j][i_]=-1000;
    p_time_pmt1[j][i_]=-1;
    p_time_pmt2[j][i_]=-1;
  };
};
