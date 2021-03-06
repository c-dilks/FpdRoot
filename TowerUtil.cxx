#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "TText.h"
#include "TCutG.h"
#include "TCanvas.h"

#include "TowerUtil.h"
using namespace std;


ClassImp(TowerUtil);

  
void  TowerUtil::TowerUtilCommon(WasExternal* Pwe)
{
  pwe=Pwe;
  h2=NULL;
  arrValley=NULL;
  neighbor=NULL;
  pPaveText=NULL;

  nNSTow = NS_ETA_NUM * NS_PHI_NUM ;
  nTBTow = TB_ETA_NUM * TB_PHI_NUM ;
  maxDistanceFromPeak=0.3;
  
  minTowerCatag02=5;
  
  cutEcSigma[0][0]=2.1;
  cutEcSigma[0][1]=7.0;
  cutEcSigma[1][0]=2.1;
  cutEcSigma[1][1]=2.0;
  
  
  minEcSigma2Ph=35.;
  maxEcSigma1Ph=10.;
  minTowerEnergy=0.01;
  minRatioPeakTower=1.6;
  SetMomentEcutoff();
};

TowerUtil::~TowerUtil()
{
  if(h2)delete h2;
  if(arrValley)delete arrValley;
  if(neighbor)delete neighbor;
  if(pPaveText)delete pPaveText;
};

Int_t TowerUtil::FindTowerCluster(TObjArray *arrTow, HitCluster *clust)
{
  
  
  // the neighbor TObjArray
  //
  neighbor = new TObjArray(20);
  neighbor->Clear();
  
  // the "valley" TObjArray
  //
  arrValley = new TObjArray(16);
  arrValley->Clear();
  
  // number of "peaks" that has the same shortest distance to a "valley" tower
  //
  Int_t nPeakSameDist[nNSTow];
  
  // store which "peaks" have the same shortest distance to a "valley" tower
  //
  Int_t peaksToValley[nNSTow][MAX_NUMER_CLUSTERS];
  
  
  Int_t nClusts = 0 ;
  
  // We assume that "TObjArray *arrTow" is already sorted, however, to be safe,
  // sort tower by energy, if not done already
  //
  if( !(arrTow->IsSorted()) )
    arrTow->Sort();
  
  TowerFPD *high;
  
  // "TObjArray::Sort()" sorts the array from small to big ( [0]<=[1]<=...<=[48] )
  // need to take care of that
  //
  // have to go last object first!
  //
  // The algoriths is such, first get the highest tower (which is ALWAYS the last one),
  // and it and all its neighbors to the next cluster. Then repeat the process over the
  // remaining towers.
  //
  while( !arrTow->IsEmpty() ) {
    
    // By design, this tower is the highest tower in "arrTow", but it could be lower
    // than a tower in "neighbor"
    //
    high = (TowerFPD *) arrTow->Last() ;
    
    // when tower energy is less than minTowerEnergy, break out the loop!
    //
    if( high->energy < minTowerEnergy )
      break;

		// 2003-08-15
		// Fix a logical loop hole in deciding if a tower is
		//    a peak; Need to first compare the highest tower
		//    remained in "arrTow" to all towers in "neighbor" which is its neighbor,
		//    and if it is lower than any of those, it is
		//    a neighbor. Move it to "neighbor" and continue to
		//    the next tower in "arrTow".
		//
    TowerFPD *nbTow;
    Bool_t isPeak;
    isPeak = true ;
    for(Int_t jn=0; jn<neighbor->GetEntriesFast(); jn++) 
      {
	nbTow = (TowerFPD *) neighbor->At(jn);
	if( high->IsNeighbor(nbTow) ) 
	  {
				//
				// 2003-09-13
				// Now require "peak" tower energy be at least "minRatioPeakTower" times
				//   that of "neighbor" tower
				//
	    if( high->energy < (minRatioPeakTower * nbTow->energy) ) 
	      {
		isPeak = false ;
		break;
	      }
	  }
      } // loop over towers in "neighbor"


    TowerFPD *resTow;
    Int_t nT;

    // if "high" is not a peak, move it to "neighbor"
		//
    if( !isPeak ) 
      {
	arrTow->Remove(high);
	neighbor->Add(high);
      }
    //
    // else if "high" is a peak
    // remove the high tower from the original TObjArray, and add it to the next cluster
    //
    else 
      {
	high->cluster = nClusts ;
	arrTow->Remove(high);
	(clust[nClusts].tow)->Add(high);
	// 		arrTow->Compress();


	// loop over rest of original TObjArray, and move any towers neighboring "high"
	// to "neighbor"
	//
	nT = arrTow->GetEntriesFast();
	for(Int_t i=nT-1; i>=0; i--) 
	  {
	    resTow = (TowerFPD *) arrTow->At(i);
	    
	    // do NOT add tower with "zero" energy to cluster!
	    // when tower energy is less than minTowerEnergy, break out the loop over "neighbor"!
	    //
	    if( resTow->energy < minTowerEnergy )
	      break;
	    
	    // if "resTow" is the immediate neighbor of "high", add it to the "neighbor"
	    //
	    if( resTow->IsNeighbor(high) ) 
	      {
		arrTow->Remove(resTow);
		neighbor->Add(resTow);
	      }
	    
	  } // loop over the rest of towers in "arrTow"
	
      } // when "high" is a "peak"
    
    
    // 2003-09-27
    // Do the follow, no matter "high" is a "peak" or "neighbor"!
    //
    // To close the logic loop hole that a tower which is seperated (by towers of the same energy) from
    //   the "neighbor" towers becomes a peak, just because it happens to be in front of those towers
    //   of the same energy (since the sorting can not distinguish that).
    //
    // We need to again loop over "arrTow", move any tower (that is neighboring any of the "neighbor"
    //   towers and also has lower (or equal) energy than that "neighbor" tower) to "neighbor"
    //
    // Every time we remove a tower from "arrTow" (except we just simply go over all items in TObjArray
    //   sequentially without worrying the relative order), we need to compress "arrTow"!
    //   Since we assume that it has no emty slot!
    //
    arrTow->Compress();
    nT = arrTow->GetEntriesFast();
    for(Int_t ii=nT-1; ii>=0; ii--) 
      {
	resTow = (TowerFPD *) arrTow->At(ii);
      
	// do NOT add tower with "zero" energy to cluster!
	// when tower energy is less than minTowerEnergy, break out the loop over "arrTow"!
	//
	if( resTow->energy < minTowerEnergy )
	  break;
      
	for(Int_t kn=0; kn<neighbor->GetEntriesFast(); kn++) 
	  {
	    nbTow = (TowerFPD *) neighbor->At(kn);
	    if( resTow->IsNeighbor(nbTow) ) 
	      {
		//
		// 2003-09-27
		// unless "resTow" is a peak ("minRatioPeakTower" times of energy of "nbTow"), it is a neighbor!
		//
		if( resTow->energy < (minRatioPeakTower * nbTow->energy) ) 
		  {
		    arrTow->Remove( resTow);
		    neighbor->Add(resTow);
		    break;
		  }
	      }
	  }
      }
    
    // compress "arrTow" again
    //
    // Every time we remove a tower from "arrTow" (except we just simply go over all items in TObjArray
    //   sequentially without worrying the relative order), we need to compress "arrTow"!
    //   Since we assume that it has no emty slot!
    //
    arrTow->Compress();


    // increment "nClusts" when we find a "peak"
    //
    if( isPeak ) 
      {
	nClusts++ ;
	if(pwe){if(pwe->YPrintLevel>0)printf("adding to cluster nClusts++ = %d \n",nClusts);};
	
	// when reaching maximum number of clusters, break out the loop!
	//
	if( nClusts >= MAX_NUMER_CLUSTERS )
	  break;
      }
    
  } // loop over "arrTow"
  
  
  
  // now that we know all the peaks, let's decide the affiliation of
  // those remote neighbor-towers in "neighbor" TObjArray
  //
  // First, we need to sort the "neighbor" TObjArray, because we want to consider the "neighbor" towers
  //    from higher towers to lower towers
  //
  neighbor->UnSort();
  neighbor->Sort();
  // 	neighbor->Compress();

  // extremely faraway distance (no distance between towers can be this large!)
  //
  const Float_t ExtremelyFaraway = 99999 ;

  // distance to peak of cluster
  //
  Float_t distToClust[MAX_NUMER_CLUSTERS] ;

  TowerFPD *nbT;
  TowerFPD *pkT;

  // All towers in "neighbor" must belong to a cluster.
  // Loop over them, check if it is bordering a cluster.
  //   If yes, move it to the nearest cluster, and move on to the next tower.
  //   If no, move on to the next tower.
  //   When reach the end of the array, start from the beginning.
  //
  Int_t jjn = neighbor->GetEntriesFast() - 1 ;

  while( !neighbor->IsEmpty() ) 
    {
      
      nbT = (TowerFPD *) neighbor->At(jjn);
      
      // towers in "neighbor" should NEVER be lower than "minTowerEnergy"
      //
      if( nbT->energy < minTowerEnergy ) 
	{
	  cout << "Something is wrong! A tower in \"neighbor\" has energy " << nbT->energy;
	  cout << ". Lower than " << minTowerEnergy << ".\n" << endl;
	}
      
      // which cluster should this tower belong?
      //
      Int_t whichCluster;
      //
      // the smallest distance to the peaks
      //
      // 2003-09-30
      // now the smallest distance to clusters
    //
    Float_t minDist;
    minDist = ExtremelyFaraway ;

    for(Int_t ic=0; ic<nClusts; ic++) 
      {

	// first set the distance to the peak of a cluster to be unreasonably high
	//
	distToClust[ic] = ExtremelyFaraway ;
	
	// peak-tower is always the first one in cluster's array
	//
	pkT = (TowerFPD *) (clust[ic].tow)->First();
	
	
	// do not consider if the peak is lower than the "neighbor" tower
	//
	// 2003-08-15
	// Because of digitization, it is possible that the "neighbor" tower
	// has the exact energy as the peak tower!!!
	// Change the condition from "if( pkT->energy <= nbT->energy )"
	// to "if( pkT->energy < nbT->energy )", so that this "neighbor" tower
	// would not be hung dry!
	//
	// 2003-09-27
	// On rare occasions, this requirement that a "neighbor" must not be higher than
	//    the "peak" could leave the "neighbor" hang on dry (infinite loop). It happens
	//    (run 4126039, EN, for example. Pedestals subtraction not properly done?!) when
	//    "neighbor" is higher than the closest "peak", but towers surrounding it all belongs
	//    to that cluster.
	//
	// The requirement is reasonable. So we have to find a way to break out of the infinite
	//    loop and print out an error message. (Counts the number of towers remaining in
	//    "neighbor". If it is the same as the last iteration, we have an infinite loop!)
	//
	if( pkT->energy < nbT->energy )
	  continue;
	
	// loop over all towers in this cluster
	//
	TowerFPD *towInClust;
	
	for(Int_t jt=0; jt<(clust[ic].tow)->GetEntriesFast(); jt++) 
	  {
	    
	    // check if "nbT" is neigboring any tower in this cluster
	    //
	    towInClust = (TowerFPD *) (clust[ic].tow)->At(jt) ;
	    
	    if( nbT->IsNeighbor( towInClust ) ) {
	      
	      // make sure that "nbT" is not more than "minRatioPeakTower" times of "towInClust"
	      //
	      if( nbT->energy > (minRatioPeakTower * towInClust->energy) ) 
		continue;
	      
	      // calculate distance to the "peak" of the cluster
	      //
	      Float_t delc, delr;
	      
	      // 2003-10-03
	      // Revert to using distance to "peak" tower.
	      //
	      // 2003-09-30
	      // Calculate cluster moment (only need center position) while adding
	      //    "neighbor" towers to clusters. Use the center of cluster to
	      //    ecide the distance from a "neighbor" tower to a cluster (previously
	      //    use the "peak" tower position). This way, a tower just in between
	      //    of two "peaks" will probably go to the cluster of lower "peak".
	      //    Thus, the fitting program will less likely to try to find a 2nd
	      //    peak in the direction of another (lower) cluster.
	      //
	      // 					CalClusterMoment(&clust[ic]);
	      //  					delc = clust[ic].x0 - (nbT->col - 0.5) ;
	      // 					delr = clust[ic].y0 - (nbT->row - 0.5) ;
	      
	      delc = pkT->col - nbT->col ;
	      delr = pkT->row - nbT->row ;
	      distToClust[ic] = sqrt( delc * delc + delr * delr ) ;
	      
	      // once if the tower is found to be a neighbor of (and does not have higher enough energy than) one tower
	      //   in the cluster, it is not necessary to contine
	      //
	      break;
	    }
	  } // loop over all towers in a cluster
	
	// check if the distance to the peak of this cluster is smaller
	//
	if( distToClust[ic] < minDist ) 
	  {
	    minDist = distToClust[ic] ;
	    whichCluster = ic ;
	  }
	
      } // loop over all clusters
    
    // move the tower to the appropriate cluster
    //
    if( minDist < ExtremelyFaraway ) 
      {
	
	// loop over all clusters, and count the number of "peaks" that have the same distance to "nbT"
	//
	Int_t numbValleyTower;
	numbValleyTower = arrValley->GetEntriesFast();
	nPeakSameDist[numbValleyTower] = 0 ;
	for(Int_t llc=0; llc<nClusts; llc++) 
	  {
	    if( distToClust[llc] == minDist ) 
	      {
		peaksToValley[numbValleyTower][ nPeakSameDist[numbValleyTower] ] = llc ;
		nPeakSameDist[numbValleyTower] ++ ;
	      }
	  }
	
	if( nPeakSameDist[numbValleyTower] == 1 ) 
	  {
	    //
	    // Only one "peak" is closest to "nbT". "nbT" belongs to this "peak"!
	    //
	    nbT->cluster = whichCluster ;
	    neighbor->Remove(nbT);
	    (clust[whichCluster].tow)->Add(nbT);
	  }
	else if( nPeakSameDist[numbValleyTower] > 1 ) 
	  {
	    neighbor->Remove(nbT);
	    arrValley->Add(nbT);
	  }
      else 
	{
	  cout << "Something wrong in your logic! nPeakSameDist = " << nPeakSameDist << "! Error!" << endl;
	}
      }
    
    // move forward on to the next tower
    //
    jjn-- ;
    if( jjn == -1 ) 
      {
	
	// Counts number of towers in "neighbor". If the next iteration does not move any tower to a cluster
	//    (same number of towers), we have an infinite loop. Break out and print out error message!
	//
	Int_t numbTowBefore;
	numbTowBefore = neighbor->GetEntriesFast() ;
	neighbor->Compress();
	jjn = neighbor->GetEntriesFast() - 1;
	
	if( numbTowBefore > 0 && (jjn + 1) == numbTowBefore ) 
	  {
	    break;
	  }
      }
    
    } // loop over TObjArray "neighbor"
  
  
  // 2003-10-12
  // calculate the moment of clusters, then decide where the "valley" towers should go
  //
  for(Int_t ic=0; ic<nClusts; ic++) 
    {
      CalClusterMoment(&clust[ic]);
    }
  
  Int_t numbVal = arrValley->GetEntriesFast();
  for(Int_t iVal = 0; iVal<numbVal; iVal++) {
    nbT = (TowerFPD *) arrValley->At(iVal);
    
    
    // which cluster should this tower belong?
    //
    Int_t whichCluster = -1;
    
    Float_t minDist;
    minDist = ExtremelyFaraway ;
    
    for(Int_t lc=0; lc<nPeakSameDist[iVal]; lc++) 
      {
	
	// which cluster is this?
	//
	Int_t jkc;
	jkc = peaksToValley[iVal][lc] ;
	
	Float_t delc, delr;
	delc = clust[jkc].x0 - (nbT->col - 0.5) ;
	delr = clust[jkc].y0 - (nbT->row - 0.5) ;
	distToClust[lc] = sqrt( delc * delc + delr * delr ) ;
	
	// check if the distance to the "center" of this cluster is smaller
	//
	if( distToClust[lc] < minDist ) 
	  {
	    minDist = distToClust[lc] ;
	    whichCluster = jkc ;
	  }
      }
    
    // move the tower to the appropriate cluster
    //
    if( minDist < ExtremelyFaraway ) 
      {
	nbT->cluster = whichCluster ;
	neighbor->Remove(nbT);
	(clust[whichCluster].tow)->Add(nbT);
      }	
    else 
      {
	cout << "Something is wrong! The following \"Valley\" tower does not belong to any cluster! Error!" << endl;
	nbT->Print();
	cout << "!!!!!!!!\n" << endl;
      }
  }
  
  // 2003-10-12
  // If there are still towers left in "neighbor", distribute them to clusters
  //
  neighbor->Compress();
  jjn = neighbor->GetEntriesFast() - 1 ;
  
  while( !neighbor->IsEmpty() ) 
    {
      
      nbT = (TowerFPD *) neighbor->At(jjn);
      
      
      // towers in "neighbor" should NEVER be lower than "minTowerEnergy"
      //
      if( nbT->energy < minTowerEnergy ) 
	{
	  cout << "Something is wrong! A tower in \"neighbor\" has energy " << nbT->energy;
	  cout << ". Lower than " << minTowerEnergy << ".\n" << endl;
	}
      
      // which cluster should this tower belong?
      //
      Int_t whichCluster;
      //
      // the smallest distance to the peaks
      //
      // 2003-09-30
      // now the smallest distance to clusters
      //
      Float_t minDist;
      minDist = ExtremelyFaraway ;
      
      for(Int_t ic=0; ic<nClusts; ic++) 
	{
	  
	  // first set the distance to the peak of a cluster to be unreasonably high
	  //
	  distToClust[ic] = ExtremelyFaraway ;
	  
	  // peak-tower is always the first one in cluster's array
	  //
	  pkT = (TowerFPD *) (clust[ic].tow)->First();
	  
	  // do not consider if the peak is lower than the "neighbor" tower
	  //
	  if( pkT->energy < nbT->energy )
	    continue;
	  
	  // loop over all towers in this cluster
	  //
	  TowerFPD *towInClust;
	  
	  for(Int_t jt=0; jt<(clust[ic].tow)->GetEntriesFast(); jt++) 
	    {
	      
	      // check if "nbT" is neigboring any tower in this cluster
	      //
	      towInClust = (TowerFPD *) (clust[ic].tow)->At(jt) ;
	      
	      if( nbT->IsNeighbor( towInClust ) ) 
		{
		  
		  // make sure that "nbT" is not more than "minRatioPeakTower" times of "towInClust"
		  //
		  if( nbT->energy > (minRatioPeakTower * towInClust->energy) ) 
		    continue;
		  
		  // calculate distance to the "peak" of the cluster
		  //
		  Float_t delc, delr;
		  
		  CalClusterMoment(&clust[ic]);
		  delc = clust[ic].x0 - (nbT->col - 0.5) ;
		  delr = clust[ic].y0 - (nbT->row - 0.5) ;
		  distToClust[ic] = sqrt( delc * delc + delr * delr ) ;
		  
		  // once if the tower is found to be a neighbor of (and does not have higher enough energy than) one tower
		  //   in the cluster, it is not necessary to contine
		  //
		  break;
		}
	    } // loop over all towers in a cluster
	  
	  // check if the distance to the peak of this cluster is smaller
	  //
	  if( distToClust[ic] < minDist ) 
	    {
	      minDist = distToClust[ic] ;
	      whichCluster = ic ;
	    }
	  
	} // loop over all clusters
      
      // move the tower to the appropriate cluster
      //
      if( minDist < ExtremelyFaraway ) 
	{
	  nbT->cluster = whichCluster ;
	  neighbor->Remove(nbT);
	  (clust[whichCluster].tow)->Add(nbT);
	}
    
    
      // move forward on to the next tower
      //
      jjn-- ;
      if( jjn == -1 ) 
	{
	  
	  // Counts number of towers in "neighbor". If the next iteration does not move any tower to a cluster
	  //    (same number of towers), we have an infinite loop. Break out and print out error message!
	  //
	  Int_t numbTowBefore;
	  numbTowBefore = neighbor->GetEntriesFast() ;
	  neighbor->Compress();
	  jjn = neighbor->GetEntriesFast() - 1;
	  
	  if( numbTowBefore > 0 && (jjn + 1) == numbTowBefore ) 
	    {
	      cout << "Infinite loop! The following towers are not claimed by any cluster!" << endl;
	      for(Int_t jbt=jjn; jbt>=0; jbt--) 
		{
		  cout << "\t";
		  ( (TowerFPD *) neighbor->At(jbt) )->Print();
		}
	      cout << "\n  Algorithm could not deal with the case when \"neighbor\" is higher than the ";
	      cout << "closest \"peak\", but towers surrounding it all belongs to that cluster!\n" << endl;
	      cout << "Check the event! No (or bad) pedestal subtraction!?? \n" << endl;
	      break;
	    }
	}
      
    } // loop over TObjArray "neighbor"
  
  
  // 2003-08-26
  //
  // Sort towers by energy (descending, higher energy towers first)
  //
  for(Int_t jc=0; jc<nClusts; jc++) 
    {
      (clust[jc].tow)->UnSort();
      (clust[jc].tow)->Sort();
      //
      // TObjArray sort ascending! Not what I want!
      // Do an exchange of towers: 0<-->nTT-1, 1<-->nTT-2, etc.
      //
      TowerFPD *tmp1;
      TowerFPD *tmp2;
      Int_t nTT;
      nTT = (clust[jc].tow)->GetEntriesFast();
      for(Int_t itt=0; itt<nTT/2; itt++) {
	tmp1 = (TowerFPD *) (clust[jc].tow)->RemoveAt(itt) ; 
	tmp2 = (TowerFPD *) (clust[jc].tow)->RemoveAt(nTT-1-itt) ;
	(clust[jc].tow)->AddAt(tmp1, nTT-1-itt);
	(clust[jc].tow)->AddAt(tmp2, itt);
      }
    }
  
  // 2003-08-30
  // put center of towers at 0.5 lgd, because this is the more natural way.
  //
  // calculate various moment of clusters
  //
  for(Int_t ic=0; ic<nClusts; ic++) 
    {
      
      CalClusterMoment(&clust[ic]);
      
      // set initial nPhoton to 0 & catag to be -1
      //
      clust[ic].nPhoton =  0 ;
      clust[ic].catag   = -1 ;
    }
  
  
  // 2003-09-14
  // now add those "zero" towers to the clusters
  // those towers serve the purpose of preventing the creation of bogus peak
  //   (peak where there is no energy deposited at the tower)
  //
  arrTow->Compress();
  for(Int_t jz=0; jz<arrTow->GetEntriesFast(); jz++) 
    {
      nbT = (TowerFPD *) arrTow->At(jz);
      
      // which cluster should this tower belong?
      //
      Int_t whichCluster = 0;
      //
      // the smallest distance to the peaks
      //
      Float_t minDist;
      minDist = ExtremelyFaraway ;
      
      // loop over all clusters
      //
      for(Int_t ijc=0; ijc<nClusts; ijc++) 
	{
	  
	  Float_t dist, delc, delr;
	  //
	  // peak-tower is always the first one in cluster's array
	  //
	  pkT = (TowerFPD *) (clust[ijc].tow)->First();
	  //
	  // distance to this peak
	  //
	  delc = nbT->col - pkT->col;
	  delr = nbT->row - pkT->row;
	  //
	  // 2003-10-11
	  // distance to "center" of cluster is used
	  //
	  // 			delc = clust[ijc].x0 - (nbT->col - 0.5) ;
	  // 			delr = clust[ijc].y0 - (nbT->row - 0.5) ;
	  dist = sqrt( delc*delc + delr*delr ) ;
	  
	  // since the higher-peak cluster is considered first, when "dist" is the same, favor the
	  // higher-peak cluster (naturally)
	  //
	  if( dist < minDist ) 
	    {
	      minDist = dist ;
	      whichCluster = ijc;
	    }
	}
      
      // if the distance is smaller than "maxDistanceFromPeak"
      //    move the tower to the appropriate cluster
      // Do not want to add too many "zero" towers to a cluster!
      //
      if( minDist < maxDistanceFromPeak ) 
	{
	  nbT->cluster = whichCluster ;
	  arrTow->Remove(nbT);
	  (clust[whichCluster].tow)->Add(nbT);
	}
    }
  
  neighbor->Clear();
  delete neighbor;
  neighbor=NULL;
  arrValley->Clear();
  delete arrValley;
  arrValley=NULL;
  return nClusts;
}



// cluster moment calculation is now a seperate function
//
void TowerUtil::CalClusterMoment(HitCluster *clust)
{
  if(clust)
    {
      clust->CalClusterMoment(Ecutoff);
    };
  clust->numbTower = clust->tow->GetEntriesFast() ;
  /*
    Float_t w0, w1, mtmp, mx, my, sigx, sigy, sigXY;
	w0 = w1 = mtmp = mx = my = sigx = sigy = sigXY = 0 ;
	
 	clust->numbTower = clust->tow->GetEntriesFast() ;
	
	
	for(Int_t it=0; it<clust->numbTower; it++) 
	  {
	    
	    TowerFPD * oneTow;
	    oneTow = (TowerFPD *) clust->tow->At(it) ;
	    Float_t xxx, yyy;
	    xxx = oneTow->col - 0.5 ;
	    yyy = oneTow->row - 0.5 ;
	    mtmp = log(oneTow->energy+1.-Ecutoff)>0 ? log(oneTow->energy+1.-Ecutoff) : 0;
	    w1 += mtmp;
	    w0    += oneTow->energy ;
	    mx    += mtmp * xxx ;
	    my    += mtmp * yyy ;
	    sigx  += mtmp * xxx * xxx ;
	    sigy  += mtmp * yyy * yyy ;
	    sigXY += mtmp * xxx * yyy ;
	  }
	
	clust->energy  = w0 ;
	if(w1>0)
	  {
	    clust->x0      = mx / w1 ;
	    clust->y0      = my / w1 ;
	    clust->sigmaX  = sqrt( fabs(sigx / w1 - clust->x0 * clust->x0 ) ) ;
	    clust->sigmaY  = sqrt( fabs(sigy / w1 - clust->y0 * clust->y0 ) ) ;
	    clust->sigmaXY = sigXY / w1 - clust->x0 * clust->y0 ;
	  }
	else
	  {
	    clust->x0      = 0 ;
	    clust->y0      = 0 ;
	    clust->sigmaX  = 0 ;
	    clust->sigmaY  = 0 ;
	    clust->sigmaXY = 0 ;
	  };
  */
	//cout<<"ene="<<clust->energy<<" x0="<<clust->x0<<" y0="<<clust->y0
	//<<" sigmaX="<<clust->sigmaX<<" sigmaY="<<clust->sigmaY
	//<<" sigmaXY="<<clust->sigmaXY<<endl;
}


// 2003-08-28
// a new way of catagorize clusters
//
// 2003-09-13
// new parameters
//
// lines of seperation:
// 1. y=1.825(x-6), below this line, (and x>19) is almost certainly 2-photon (1 exception)
// 2. y=1.825(x-0), above this, (and x<25) is almost certainly 1-photon (some exceptions, but probably very weak 2nd photon)
//
// catagozie cluster (1-photon, 2-photon, or need-to-fit-for-both-and-see)
//
Int_t TowerUtil::CatagBySigmXY(HitCluster *clust)
{
  // only need to be called once "ReadParamters()"
  //
  
  
  // if the number of towers in a cluster is less than "minTowerCatag02"
  //    consider the cluster a catag-1 cluster (with only 1 gamma)
  //
  if( clust->numbTower < minTowerCatag02 ) 
    {
      clust->catag = 1 ;
      return 1;
    }
  
  Float_t sMaxEc = clust->sigmaMax * clust->energy;

  if( clust->energy < cutEcSigma[0][0] * ( sMaxEc - cutEcSigma[0][1] ) ) 
    {
      if( sMaxEc > minEcSigma2Ph ) 
	{
	  clust->catag = 2 ;
	  return 2;
	}
      else 
	{
	  clust->catag = 0 ;
	  return 0;
	}
    }
  else if( clust->energy > cutEcSigma[1][0] * ( sMaxEc - cutEcSigma[1][1] ) ) 
    {
      if( sMaxEc < maxEcSigma1Ph ) 
	{
	  clust->catag = 1 ;
	  return 1;
	}
      else 
	{
	  clust->catag = 0 ;
	  return 0;
	}
    }
  else 
    {
      clust->catag = 0 ;
    }
  
  return clust->catag ;
}



void TowerUtil::PrintTowers(const TowerFPD *tows)
{
  
  for(Int_t j=0; j<nNSTow; j++) 
    {
      if( j%7 == 0 ) printf("\n");
      printf("%3d", tows[j].cluster+1);
    }
  
  for(Int_t j=0; j<nNSTow; j++) 
    {
      if( j%7 == 0 ) printf("\n");
      printf("%10.4f", tows[j].energy);
    }
  
  printf("\n\n");
  
}


void TowerUtil::VisualizeFit(const Int_t iEv, const Int_t nRow, const Float_t widthLG, const TowerFPD *tows, const Int_t nPh, const Float_t *xPh, const Float_t *yPh, const Float_t *ePh, const Float_t *chi2, TF2 *f2, char path[64])
{
  
  TCanvas c1("c1", "test", 700, 500);
  
  // 	cout << "Event #" << iEv << ":" << endl;
  
  // first fill the tower energy hitogram
  //
  Float_t boundary[nRow+1];
  for(Int_t i=0; i<=nRow; i++) 
    {
    boundary[i] = i * widthLG ;
    }
  
  h2 = new TH2F("h2", Form("Event #%d, %d photons", iEv, nPh), nRow, boundary, nRow, boundary);
  mySetHistTitleLabelProperty(h2);
  h2->SetXTitle("x (cm)");
  h2->SetYTitle("y (cm)");
  h2->SetZTitle("E (GeV)");
  
  for(Int_t j=0; j<nRow*nRow; j++) 
    {
      h2->SetBinContent(tows[j].col, tows[j].row, tows[j].energy);
    }
  
  Double_t xyMax = nRow * widthLG ;
  Double_t zMax = h2->GetMaximum() * 1.2;
  
  h2->SetStats(false);
  h2->SetMaximum(zMax);
  h2->Draw("lego2");
  
  pPaveText = new TPaveText(0.72,0.80,0.99,0.99,"NDC");
  pPaveText->SetTextSizePixels(8);
  
  f2->SetMaximum(zMax);
  
  TCutG * cutg[nPh];
  
  // now draw the shape of fitted photons
  //
  for(Int_t kp=0; kp<nPh; kp++) 
    {
      
      // try to avoid color "5" (yellow! hard to see the numbers on canvas)
      //
      Int_t myColor;
      myColor = kp+2;
      if( myColor == 5 ) 
	{
	  myColor = 2 + nPh;
	}
      
      TText *pText;
      pText = pPaveText->AddText(Form("#gamma #%d: (%3.1f,%3.1f), E=%4.1f, #chi^{2}=%4.1f", kp, xPh[kp], yPh[kp], ePh[kp], chi2[kp]));
      pText->SetTextColor(myColor);
      
      Double_t xcutga[4], ycutga[4];
      xcutga[0] = xPh[kp]-4 ;
      if( xcutga[0] < 0 )
	xcutga[0]	= 0 ;
      
      xcutga[1]	= xcutga[0];
      
      xcutga[3] = xPh[kp]+4 ;
      if( xcutga[3] > xyMax )
	xcutga[3] = xyMax ;
      
      xcutga[2] = xcutga[3];
      
      ycutga[0] = yPh[kp]-4 ;
      if( ycutga[0] < 0 )
	ycutga[0]	= 0 ;

      ycutga[3]	= ycutga[0];
      
      ycutga[2] = yPh[kp]+4 ;
      if( ycutga[2] > xyMax )
	ycutga[2] = xyMax ;
      
      ycutga[1] = ycutga[2];
      
      // cuts that will limit the drawing range of "f2"
      //
      
      cutg[kp] = new TCutG(Form("ev%d_cutg%d", iEv, kp), 4, xcutga, ycutga);
      
      f2->SetParameter(7, xPh[kp]);
      f2->SetParameter(8, yPh[kp]);
      f2->SetParameter(9, ePh[kp]);
      f2->SetLineColor(myColor);
      
      TF2 *ff = (TF2 *) f2->Clone( Form("f2%d", kp) );
      ff->Draw( Form("surf same [ev%d_cutg%d]", iEv, kp) );
      
      // 		cutg->Draw("C");
      
      // 		Double_t param[10];
      // 		ff->GetParameters(param);
      // 		for(Int_t iip=0; iip<10; iip++) {
      // 			cout << param[iip] << ",\t";
      // 		}
      // 		cout << "\n" << endl;
    }
  
  pPaveText->AddText(Form("global fit #chi^{2}=%5.2f", chi2[nPh]) );
  pPaveText->SetAllWith(":", "align", 11);
  // 	pPaveText->Draw("SAME");
  
  
  for(Int_t jp=0; jp<nPh; jp++) 
    {
      cutg[jp]->Draw("C");
      // 		cutg->Write();
    }
  
  c1.SaveAs(Form("%s/event_%d.root", path, iEv));
  
  for(Int_t jp=0; jp<nPh; jp++) 
    {
      delete cutg[jp];
    }
  
  // 	f.Close();
  
  delete h2;
  h2=0;
  
};

void TowerUtil::mySetHistTitleLabelProperty(TH1* h)
{
  
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelOffset(0.002);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTickLength(0.03);
  
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelOffset(0.002);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.03);
  
  h->GetZaxis()->SetTitleOffset(1.0);
  h->GetZaxis()->SetTitleSize(0.045);
  h->GetZaxis()->SetLabelOffset(0.002);
  h->GetZaxis()->SetLabelSize(0.04);
  h->GetZaxis()->SetTickLength(0.03);
};

