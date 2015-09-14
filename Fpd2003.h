///////////////////////////////////////////////
//
// 2003-10-03
//
// Increase "MAX_NUMB_PHOTONS" from 4 to 7
// Increase "MAX_NUMER_CLUSTERS" from 4 to 6
//
////////////////////////////////////////////////
//
//  2003-07-02
//
//  Define common FPD FY2003 setup
//  Detector channel mapping is in 
//    "FPD2003Map.h"
//
///////////////////////////////////////////////

#ifndef FPD2003_H
#define FPD2003_H

#include "TROOT.h"

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
const Int_t nNSTow = NS_ETA_NUM * NS_PHI_NUM ;
const Int_t nTBTow = TB_ETA_NUM * TB_PHI_NUM ;

#define MAX_NUMB_PHOTONS 7

// maximum number of clusters that will can be handled
//
#define MAX_NUMER_CLUSTERS 6


#endif


