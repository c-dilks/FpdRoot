#include "Yiqun.h"
using namespace std;


ClassImp(Yiqun)


// the real pointers are declared here
//
//-----------------------------static

  WasExternal FitTower::we;

Float_t Yiqun::FitOnePhoton(HitCluster* p_clust)
{
  if(pwe->YPrintLevel>1)cout<<"Call FitPePhoton\n";


	Double_t chiSq;

	// fit parameters(starting positions), errors, and gradients of function
	//
	Double_t param[4];
	Double_t error[4];
	Double_t gradient[4];

	// starting position, lower and upper limit of parameters
	//
	Double_t start[4], lowLim[4], upLim[4];

	// fit status, and flag
	//
	Int_t status, iflag=1;


	// constant parameter: 1 photon
	//
	start[0] = 1;
	lowLim[0] = 0.5;
	lowLim[0] = 1.5;

	start[1] = widLG[0] * p_clust->x0;
	start[2] = widLG[1] * p_clust->y0;
	start[3] = p_clust->energy;

	if( pwe->useLimitForPara == 1 ) {
		//
		// use limit
		//
		lowLim[1] = start[1] - posDif_1PC * widLG[0] ;
		lowLim[2] = start[2] - posDif_1PC * widLG[1] ;
		lowLim[3] = start[3] * (1 - eneRat_1PC) ;
	
		upLim[1] = start[1] + posDif_1PC * widLG[0] ;
		upLim[2] = start[2] + posDif_1PC * widLG[1] ;
		upLim[3] = start[3] * (1 + eneRat_1PC) ;
	}
	else if( pwe->useLimitForPara == 2 ) {
		//
		// just limit the energy
		//
		lowLim[1] = upLim[1] = lowLim[2] = upLim[2] = 0 ;
		lowLim[3] = 0;
		upLim[3]  = 100;
	}
	else {
		// no limit when "useLimitForPara==0"
		//
		for(Int_t jkp=1; jkp<=3; jkp++) {
			lowLim[jkp] = upLim[jkp] = 0 ;
		}
	}

	// call fitter->Fit(...), return status
	//
	status = fitter->Fit(start, pwe->step, lowLim, upLim);


	// 2003-10-09
	// check return status
	//
	if( status != 0 ) {
		std::cout << "Minuit fit returns " << status << "!" << "\n";
	}


	// get the fit result
	//
	fitter->fMn->GetParameter(0, param[0], error[0]);
	Int_t nPar = 3*((Int_t) param[0])+1;
	for(Int_t ipar=1; ipar<nPar; ipar++) {
		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
	}


	// put the fit result back in "clust"
	//
	p_clust->photon[0].xPos    = param[1];
	p_clust->photon[0].errXPos = error[1];
	p_clust->photon[0].yPos    = param[2];
	p_clust->photon[0].errYPos = error[2];
	p_clust->photon[0].energy  = param[3];
	p_clust->photon[0].errEne  = error[3];


	// evaluate the Chi-square function
	//
	fitter->fMn->Eval(4, gradient, chiSq, param, iflag);

	if(PRINT_FIT_1_RESULT)
	  {
	    printf("chiSq = %f\n", chiSq);
	    printf(" start    step   lowLim   upLim     par    error \n");
	    for(Int_t jpar=0; jpar<nPar; jpar++) {
	      printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", start[jpar],
		     pwe->step[jpar],lowLim[jpar],upLim[jpar],
		     param[jpar],error[jpar]);
	    };
	  };
	return chiSq;
}


Float_t Yiqun::FitTwoPhoton(HitCluster* p_clust)
{
  if(pwe->YPrintLevel>1)cout<<"Call FitTwoPhoton\n";
	Double_t chiSq;

	// fit parameters(starting positions), errors, and gradients of function
	//
	Double_t param[7];
	Double_t error[7];
	Double_t gradient[7];

	// starting position, lower and upper limit of parameters
	//
	Double_t start[7], lowLim[7], upLim[7];

	// fit status, and flag
	//
	Int_t status, iflag=1;

	if(pwe->YPrintLevel>0)
	  { printf("New 2Phot fit on cluster\n");
	    p_clust->Print();
	  }
	// constant parameter: 2 photon
	//
	start[0] = 2;
	lowLim[0] = 1.5;
	lowLim[0] = 2.5;

	start[1]  = widLG[0] * p_clust->photon[0].xPos;
	start[2]  = widLG[1] * p_clust->photon[0].yPos;
	start[3]  = p_clust->photon[0].energy;

	start[4]  = widLG[0] * p_clust->photon[0].xPos;
	start[5]  = widLG[1] * p_clust->photon[0].yPos;
	start[6]  = p_clust->photon[0].energy;

	if( pwe->useLimitForPara == 1 ) {
		lowLim[1] = start[1] - 0.7 * widLG[0] ;
		lowLim[2] = start[2] - 0.7 * widLG[1] ;
		lowLim[3] = start[3] * 0.6 ;
	
		upLim[1]  = start[1] + 0.7 * widLG[0] ;
		upLim[2]  = start[2] + 0.7 * widLG[1] ;
		upLim[3]  = start[3] * 1.6 ;

		lowLim[4] = start[4] - 0.6 * widLG[0];
		lowLim[5] = start[5] - 0.6 * widLG[1];
		lowLim[6] = start[6] * 0.6 ;
	
		upLim[4]  = start[4] + 0.6 * widLG[0];
		upLim[5]  = start[5] + 0.6 * widLG[1];
		upLim[6]  = start[6] * 1.7 ;
	}
	else if( pwe->useLimitForPara == 2 ) {
		//
		// just limit the energy
		//
		lowLim[1] = upLim[1] = lowLim[2] = upLim[2] = lowLim[4] = upLim[4] = lowLim[5] = upLim[5] = 0 ;
		lowLim[3] = lowLim[6] = 0;
		upLim[3]  = upLim[6]   = 100;
	}
	else {
		// no limit when "useLimitForPara==0"
		//
		for(Int_t jkp=1; jkp<=6; jkp++) {
			lowLim[jkp] = upLim[jkp] = 0 ;
		}
	}

	// call fitter->Fit(...), return status
	//
	status = fitter->Fit(start, pwe->step, lowLim, upLim);


	// 2003-10-09
	// check return status
	//
	if( status != 0 ) {
		std::cout << "Minuit fit returns " << status << "!" << "\n";
	}


	// get the fit result
	//
	fitter->fMn->GetParameter(0, param[0], error[0]);
	Int_t nPar = 3*((Int_t) param[0])+1;
	for(Int_t ipar=1; ipar<nPar; ipar++) {
		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
	}

	// put the fit result back in "clust"
	//
	p_clust->photon[0].xPos    = param[1];
	p_clust->photon[0].errXPos = error[1];
	p_clust->photon[0].yPos    = param[2];
	p_clust->photon[0].errYPos = error[2];
	p_clust->photon[0].energy  = param[3];
	p_clust->photon[0].errEne  = error[3];

	p_clust->photon[1].xPos    = param[4];
	p_clust->photon[1].errXPos = error[4];
	p_clust->photon[1].yPos    = param[5];
	p_clust->photon[1].errYPos = error[5];
	p_clust->photon[1].energy  = param[6];
	p_clust->photon[1].errEne  = error[6];

	// evaluate the Chi-square function
	//
	fitter->fMn->Eval(7, gradient, chiSq, param, iflag);

	if(PRINT_FIT_2_RESULT){
	  printf("chiSq = %f\n", chiSq);
	  printf(" start    step   lowLim   upLim     par    error \n");
	  for(Int_t jpar=0; jpar<nPar; jpar++) {
	    printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", start[jpar],
		   pwe->step[jpar],lowLim[jpar],upLim[jpar],param[jpar],
		   error[jpar]);
	  }
	};
	return chiSq;
}




// 2003-10-10
// Global fit, return x (or y) positions (or errors) in "cm"!
//
// fit "nPh" of photons, using "*fitter" on "nCl" of clusters "*clust"
//
// Check "nCl>=1 && nPh>=2", this is a global fit. Restrict the range accordingly.
//

Float_t Yiqun::GlobalFit(const Int_t nPh, const Int_t nCl, HitCluster *p_clust)
{
  if(pwe->YPrintLevel>1)cout<<"Call GlobalFit\n";

	// by design, we can only fit up to "MAX_NUMB_PHOTONS" (currently 4) photons
	//
	if( nPh > MAX_NUMB_PHOTONS || nPh < 2 ) {
	  std::cout << "Global fit! Can not fit " << nPh << " photons! ERROR!" << "\n";
		return -9999;
	}

	// check that there is at least one cluster
	//
	if( nCl < 1 ) {
		std::cout << nCl << " clusters! Global fit will NOT work! ERROR!" << "\n";
		return -9999;
	}


	// fit parameters(starting positions), errors, and gradients of function
	//
	Double_t param[3*MAX_NUMB_PHOTONS+1];
	Double_t error[3*MAX_NUMB_PHOTONS+1];
	Double_t gradient[3*MAX_NUMB_PHOTONS+1];

	// starting position, lower and upper limit of parameters
	//
	Double_t start[3*MAX_NUMB_PHOTONS+1], lowLim[3*MAX_NUMB_PHOTONS+1], upLim[3*MAX_NUMB_PHOTONS+1];


	// at least 2 photons: global fit!
	// The positions (such as "clust[ic].photon[jp].xPos") are already in unit of "cm"!!!
	//
	// here all clusters have all their fields properly filled (for example cluster[].photon[0] should NOT be NULL!)
	//
	Int_t totPh;
	totPh = 0;

	// loop over all clusters
	//
	for(Int_t ic=0; ic< nCl; ic++) {

		// loop over all photons in cluster
		//
		for(Int_t jp=0; jp<p_clust[ic].nPhoton; jp++)
		  {
		    
		    if( totPh > MAX_NUMB_PHOTONS ) {
		      std::cout << "Total # of photons in " << nCl << " clusters is at least " << totPh << "! I can NOT do fit! ERROR!" << "\n";
		      return -9999;
		    }
		    
		    Int_t kpar;
		    kpar = 3 * totPh + 1 ;
		    
		    start[ kpar] = p_clust[ic].photon[jp].xPos;
		    lowLim[kpar] = start[kpar] - posDif_Gl ;        // in "cm", not in unit of lgd!
		    upLim[ kpar] = start[kpar] + posDif_Gl ;        // in "cm", not in unit of lgd!
		    
		    kpar++ ;
		    start[ kpar] = p_clust[ic].photon[jp].yPos;
		    lowLim[kpar] = start[kpar] - posDif_Gl ;        // in "cm", not in unit of lgd!
		    upLim[ kpar] = start[kpar] + posDif_Gl ;        // in "cm", not in unit of lgd!
		    
		    kpar++ ;
		    start[ kpar] = p_clust[ic].photon[jp].energy;
		    lowLim[kpar] = start[kpar] * (1 - eneRat_Gl) ;
		    upLim[ kpar] = start[kpar] * (1 + eneRat_Gl) ;

		    totPh++ ;
		  }			
	}
	
	if( totPh != nPh ) {
		std::cout << "WARNING! Total # of photons in " << nCl << " clusters is at least " << totPh << "! Not the same as the nPh = ";
		std::cout  << nPh << "! I will try " << totPh << " instead!" << "\n";
	}
	start[0] = totPh ;

	lowLim[0] = 0.5;
	lowLim[0] = MAX_NUMB_PHOTONS + 0.5 ;


// 	std::cout << "useLimitForPara = " << pwe->useLimitForPara << "\n";


	// 2003-09-23
	// no limit when "useLimitForPara==0" !!!
	//
	if( pwe->useLimitForPara == 0 ) {
		for(Int_t jkp=1; jkp<=3*nPh; jkp++) {
			lowLim[jkp] = upLim[jkp] = 0 ;
		}
	}
	else if( pwe->useLimitForPara == 2 ) {
		//
		// just limit the energy
		//
		for(Int_t jkp=1; jkp<=3*nPh; jkp++) {
			lowLim[jkp] = upLim[jkp] = 0 ;
			if( jkp%3 == 0 ){
				lowLim[jkp] = 0;
				upLim[jkp]  = 100;
			}
		}
	}


	// fit status, and flag
	//
	Int_t status, iflag=1;

	// call fitter->Fit(...), return status
	//
	status = fitter->Fit(start, pwe->step, lowLim, upLim);


	// 2003-10-09
	// check return status
	//
	if( status != 0 ) {
		std::cout << "Minuit fit returns " << status << "!" << "\n";
	}


	// get the fit result
	//
	fitter->fMn->GetParameter(0, param[0], error[0]);
	Int_t nPar = 3*((Int_t) param[0])+1;
	for(Int_t ipar=1; ipar<nPar; ipar++) {
		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
	}


	// put the fit result back in "clust"
	//
	// loop over all clusters
	//
	Int_t tPh = 0 ;
	for(Int_t ic=0; ic< nCl; ic++) {

		// loop over all photons in cluster
		//
		for(Int_t jp=0; jp<p_clust[ic].nPhoton; jp++) {

			Int_t kpar;
			kpar = 3 * tPh + 1 ;

			p_clust[ic].photon[jp].xPos    = param[kpar] ;
			p_clust[ic].photon[jp].errXPos = error[kpar] ;

			p_clust[ic].photon[jp].yPos    = param[kpar+1] ;
			p_clust[ic].photon[jp].errYPos = error[kpar+1] ;

			p_clust[ic].photon[jp].energy  = param[kpar+2] ;
			p_clust[ic].photon[jp].errEne  = error[kpar+2] ;

			tPh ++;
		}

	}

	// evaluate the Chi-square function
	//
	Double_t chiSq;
	fitter->fMn->Eval(7, gradient, chiSq, param, iflag);

	if(PRINT_FIT_ALL_RESULT){
	  printf("chiSq = %f\n", chiSq);
	  printf(" start    step   lowLim   upLim     par    error \n");
	  for(Int_t jpar=0; jpar<nPar; jpar++) {
	    printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", start[jpar],
		   pwe->step[jpar],lowLim[jpar],upLim[jpar],param[jpar],
		   error[jpar]);
	  }
	};
	return chiSq;
	
}



// special initial step for 2-photon-cluster
//

// 2003-08-30
// special 2-photon cluster fitting
//
Float_t Yiqun::Fit2PhotonClust(HitCluster* p_clust)
{
  const Double_t step2[7] = {0, 0.02, 0.02, 0.01, 0.01, 0.01, 0.1} ;
    if(pwe->YPrintLevel>1)cout<<"Call Fit2PhotonClust\n";
  // random number (for Z_gg)
  //
  if( !(pwe->myRand) ) {
    std::cout << "Initialize random generator (for Z_gg and d_gg)!" << "\n";
    (pwe->myRand) = new TRandom();
  }
  
  // sigma of cluster
  //
  // 	Double_t sigma;
  // 	sigma  = clust->sigmaX * clust->sigmaX ;
  // 	sigma += clust->sigmaY * clust->sigmaY ;
  // 	sigma = sqrt(sigma);
  // ratio of sigmaMin/sigmaMax
  //
  Double_t ratioSigma = p_clust->sigmaMin / p_clust->sigmaMax ;
  
  // possible maximum angle
  //
  // 	Double_t maxTheta = asin( thetaPara * ratioSigma) ;
  //
  // 2003-10-11
  // tigher limits
  //
  Double_t maxTheta = ratioSigma / thetaPara ;
  if( maxTheta > (TMath::Pi()/2.0) ) {
    maxTheta = TMath::Pi() / 2.0 ;
  }
  // use for restrict d_gg
  //
  Double_t EcSigmaMax = p_clust->energy * p_clust->sigmaMax ;
  
  
  Double_t chiSq;
  
  // fit parameters(starting positions), errors, and gradients of function
  //
  Double_t param[7];
  Double_t error[7];
  Double_t gradient[7];
  
  // starting position, lower and upper limit of parameters
  //
  Double_t start[7], lowLim[7], upLim[7];
  
  // fit status, and flag
  //
  Int_t status, iflag=1;
  
  
  // constant parameter: 2 photon
  //
  start[0] = 2;
  lowLim[0] = 1.5;
  lowLim[0] = 2.5;
  
  
  // parameter starting points and limits are determined by looking at cluster information
  //
  //  xPi and yPi:   rarely do they go beyond 0.3 unit of lgd
  //  theta:         have a narrow theta range (for r=sigmaMin/sigmaMax, |theta|<0.5*r/0.65
  //                      when r<0.65, and linear increase from 0.5 to Pi/2 for 0.65<r<1)
  //  E_gg:          given by Ec (+/- 20% or less)
  //  z_gg:          should just let it vary from -1 to 1.
  //  d_gg:          a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
  //                      d_gg > Max( 2.5*(r-0.6), 0.5 )
  //
  start[1]  = widLG[0] * p_clust->x0;
  start[2]  = widLG[1] * p_clust->y0;
  start[6]  = p_clust->energy;
  
  start[4]  = p_clust->thetaAxis;
  
  // 2003-10-10
  // d_gg should be more closely related to sigma_max!
  //
  // 2003-10-09
  // Change start[3] (related to  2 * sigma)
  //
  // 	Double_t dggPi0Peak;
  // 	dggPi0Peak = towerWidthCM * sigma * 2.0;
  //
  start[3] = dggPara[1] * widLG[0] * p_clust->sigmaMax ;
  
  
  // randomize the starting point of Z_gg (from -0.1 to 0.1)
  //
  start[5]  = 0.1 * (2 * (pwe->myRand)->Rndm() - 1) ;
  if(pwe->YPrintLevel>1)
    {
      cout<<"start values:";
      for(int ipa=0;ipa<6;ipa++)
	{
	  cout<<start[ipa]<<" ";
	};
      cout<<"\n";
    }
  
  if( pwe->useLimitForPara == 1 ) {
    lowLim[1] = start[1] - posDif_2PC * widLG[0] ;
    lowLim[2] = start[2] - posDif_2PC * widLG[1] ;
    lowLim[6] = start[6] * (1 - eneRat_2PC) ;
    
    upLim[1]  = start[1] + posDif_2PC * widLG[0] ;
    upLim[2]  = start[2] + posDif_2PC * widLG[1] ;
    upLim[6]  = start[6] * (1 + eneRat_2PC) ;
    
    lowLim[4] = start[4] - maxTheta ;
    lowLim[5] = - 1.0 ;
    
    // 		// old limit
    // 		//
    // 		lowLim[3] = dggPara[1] * (dggPara[0] - EcSigmaMax ) ;
    //
    // 2003-10-12
    // new limit
    //
    lowLim[3] = dggPara[0] / pow(EcSigmaMax, 0.8) ;
    if( lowLim[3] < dggPara[2] ) {
      lowLim[3] = dggPara[2] ;
    }
    lowLim[3] *= widLG[0];
    if( lowLim[3] >= start[3] ) {
      lowLim[3] = start[3] * 0.9 ;
    }
    
    upLim[3] = dggPara[4] * (dggPara[3] - EcSigmaMax ) ;
    if( upLim[3] > dggPara[5] ) {
      upLim[3] = dggPara[5] ;
    }
    upLim[3] *= widLG[0] ;
    if( upLim[3]  <= start[3] ) {
      upLim[3] = start[3] * 1.1 ;
    }
    
    upLim[4]  = start[4] + maxTheta ;
    upLim[5]  =   1.0 ;
  }
  else if( pwe->useLimitForPara == 2 ) {
    //
    // limit the energy (and Z_gg), d_gg (in "cm")
    //
    
    lowLim[1] = upLim[1] = lowLim[2] = upLim[2] =  lowLim[4] = upLim[4] = 0 ;
    lowLim[3] = 0.4 * widLG[0] ;
    upLim[3]  = 7.0 * widLG[0] ;
    lowLim[5] = -1.0;
    upLim[5]  =  1.0;
    lowLim[6] = 0;
    upLim[6]  = 100;
  }
  else {
    // no limit when "useLimitForPara==0"
    //
    for(Int_t jkp=1; jkp<=6; jkp++) {
      lowLim[jkp] = upLim[jkp] = 0 ;
    }
  }
  
  
  // 	// 2003-10-13
  // 	// fix E_total and theta
  // 	//
  // 	fitter->fMn->FixParameter(6);
  // 	fitter->fMn->FixParameter(4);
  
  
  // call special 2-photon-cluster fitter
  //
  status = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim);
  
  // 2003-10-09
  // check return status
  //
  if( status != 0 ) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }
  
  
  // get the fit result
  //
  fitter->fMn->GetParameter(0, param[0], error[0]);
  Int_t nPar = 3*((Int_t) param[0])+1;
  for(Int_t ipar=1; ipar<nPar; ipar++) {
    fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
  }
  //  float mm=param[6]*sqrt(1-param[5]*param[5])*param[3]/730./2.;
  // printf("mm=%f fitx=%f fity=%f dgg=%f zgg=%f egg=%f \n",mm,
  //	 param[1],param[2],param[3],param[5],param[6]);

  // evaluate the Chi-square function
  //
  fitter->fMn->Eval(7, gradient, chiSq, param, iflag);

  if(PRINT_FIT_2_RESULT){
    printf("1st fit:\n");
    printf("chiSq = %f\n", chiSq);
    printf(" start    step   lowLim   upLim     par    error \n");
    for(Int_t jpar=0; jpar<nPar; jpar++) {
      printf("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f\n", start[jpar],step2[jpar],lowLim[jpar],upLim[jpar],param[jpar],error[jpar]);
    }
  };
  
  
  // 	// 2003-10-13
  // 	// 2nd fit: use result of 1st fit as starting point, and free all fixed parameter
  // 	//
  // 	for(Int_t ipar=1; ipar<nPar; ipar++) {
  // 		start[ipar] = param[ipar];
  // 	}
  // 	fitter->fMn->mnfree(0);
  
  // 	// call special 2-photon-cluster fitter
  // 	//
  // 	status = fitter->Fit2Pin1Clust(start, step2, lowLim, upLim);
  
  
  // 2003-10-09
  // check return status
  //
  if( status != 0 ) {
    std::cout << "Minuit fit returns " << status << "!" << "\n";
  }
  
  
  // 	// get the fit result
  // 	//
  // 	fitter->fMn->GetParameter(0, param[0], error[0]);
  // 	for(Int_t ipar=1; ipar<nPar; ipar++) {
  // 		fitter->fMn->GetParameter(ipar, param[ipar], error[ipar]);
  // 	}
  
  // 	// evaluate the Chi-square function
  // 	//
  // 	fitter->fMn->Eval(7, gradient, chiSq, param, iflag);
  
  
  // put the fit result back in "clust"
  //
  p_clust->photon[0].xPos    = param[1] + cos(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errXPos = error[1] + (cos(param[4])*error[3]-error[4]*sin(param[4])*param[3])*(1-param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].yPos    = param[2] + sin(param[4]) * param[3] * (1 - param[5]) / 2.0 ;
  p_clust->photon[0].errYPos = error[2] + (sin(param[4])*error[3]+error[4]*cos(param[4])*param[3])*(1-param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[0].energy  = param[6] * (1 + param[5]) / 2.0 ;
  p_clust->photon[0].errEne  = error[6]*(1+param[5])/2.0 + param[6]*error[5]/2.0;
  
  p_clust->photon[1].xPos    = param[1] - cos(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errXPos = error[1] + (-cos(param[4])*error[3]+error[4]*sin(param[4])*param[3])*(1+param[5])/2 - cos(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].yPos    = param[2] - sin(param[4]) * param[3] * (1 + param[5]) / 2.0 ;
  p_clust->photon[1].errYPos = error[2] + (sin(param[4])*error[3]-error[4]*cos(param[4])*param[3])*(1+param[5])/2 - sin(param[4])*param[3]*error[5]/2.0;
  p_clust->photon[1].energy  = param[6] * (1 - param[5]) / 2.0 ;
  p_clust->photon[1].errEne  = error[6]*(1-param[5])/2.0 - param[6]*error[5]/2.0;
  
  
  
  
  // 2003-10-13
  // do a global fit, using result of 1st fit as starting point
  //
  // need to set "nPhoton" before calling "GlobalFit(..)"
  //
  p_clust->nPhoton = 2;
  if(pwe->DoGlobal) chiSq = GlobalFit(2, 1, p_clust);
  
  
  if(PRINT_FIT_2_RESULT){
    printf("2nd fit:\n");
    printf("chiSq = %f\n", chiSq);
    printf("\n\n");
  };
  return chiSq;
};


Int_t Yiqun::FitEvent(Int_t nTows, Int_t &nClusts, Int_t &nRealClusts, Double_t &chiSqG, Bool_t &junkyEvent)
{  
  if(pwe->YPrintLevel>1)cout<<"Call FitEvent\n";
  if(pwe->YPrintLevel>0)
    {PRINT_FIT_2_RESULT=true;} else {PRINT_FIT_2_RESULT=false;}
	// possible alternative clusters for 1-photon fit: for catagory 0
	//

	HitCluster altClu ;
	

	// number of photons found by fit
	//
	Int_t nPh = nClusts = 0 ;

	TObjArray arrTow(nTows);

	// at the beginning of each event, first need to clear the TObjArray
	//
	
	arrTow.Clear();
	for(Int_t ic=0; ic<MAX_NUMER_CLUSTERS; ic++) {
	  if(clust[ic].tow)clust[ic].tow->Clear();
		clust[ic].Clear();
	}
	// re-add the TowerFPD objects to "arrTow" for each event!
	//
	for(Int_t itow=0; itow<nTows; itow++) {
	  arrTow.AddAt(&towers[itow], itow);
	}

	// find clusters
	//
	// "TObjArray::Sort()" sort array from small to large ( [0]<[1]<...<[48] )
	//
	// Need to first mark "TObjArray *arrTow" as unsorted! Then sort!
	//
	// calling "TObjArray::GetAbsLast()" will restore "fLast" to a value other than -2 (-2 is setho after "Sort()" is called) ?!
	//
	arrTow.UnSort();
	arrTow.Sort();
	if(pwe->YPrintLevel>0)printf("Yiqun:Call FindTowerCluster\n");
	nClusts = pTowerUtil->FindTowerCluster(&arrTow, clust);
	if(pwe->YPrintLevel>0)
	  {
	    printf("Yiqun:Back from FindToweCluster nClust=%d \n",nClusts);
	  };

	if(pwe->SubClu2)
	  {
	    /*
	      try alternate clusters
	    */

	    int nclu2=Clu2->GetEntries();
	    TIter clnext(Clu2);
	    HitCluster* cltmp;
	    int cnt=0;
	    /*	if(nClusts>0)printf("was ecluster0=%f  nClusts = %d ",clust[0].energy,
		nClusts);*/
	    if(pwe->YPrintLevel>0)
	      {
		printf("Yiqun:alternate clusters\n");
	      };
	    for(int j=0;j<nClusts;j++)
	      {
		clust[j].FindClusterAxis(pTowerUtil->GetMomentEcutoff());
	      };
	    while(cltmp=(HitCluster*) clnext())
	      {
		if(cnt<=MAX_NUMER_CLUSTERS)
		  {
		    if(cnt>0)
		      {
			if(sqrt(pow(cltmp->x0-clust[0].x0,2)+
				pow(cltmp->y0-clust[0].y0,2))<30)
			  {
			    clust[0].AddCluster(cltmp);
			    clust[0].CalClusterMoment(pTowerUtil->GetMomentEcutoff());
			  }
			else
			  {
			    clust[cnt]=*cltmp;
			    cnt++;
			  };
		      } 
		    else 
		      {
			clust[cnt]=*cltmp;
			cnt++;
		      };
		  };
		nClusts=cnt;
	      };
	    for(int j=0;j<nClusts;j++)
	      {
		Float_t saveEcutoff=pTowerUtil->GetMomentEcutoff();
		Float_t cutoff[5]={.2,.4,.5,.6,1.0};
		for(int jf=0;jf<5;jf++)
		  {
		    HitCluster* scl =new HitCluster();
		    pTowerUtil->SetMomentEcutoff(cutoff[jf]);  
		    *scl=clust[j];		
		    pTowerUtil->CalClusterMoment(&(clust[j]));
		    scl->FindClusterAxis(pTowerUtil->GetMomentEcutoff());
		    ShadowClusters->Add(scl);
		    clust[j].LinkShadowPointer(scl);
		  };
		
		pTowerUtil->SetMomentEcutoff(saveEcutoff);  
		pTowerUtil->CalClusterMoment(&(clust[j]));
		clust[j].FindClusterAxis(pTowerUtil->GetMomentEcutoff());
	      };
	    
	  };
	
	


	// 2003-09-12
	//
	// Cluster energy should be at least 2 GeV (parameter "minRealClusterEne")
	//
	// 2003-09-16
	// New algorithm to move all clusters above "minRealClusterEne" to front, and
	//   clusters below "minRealClusterEne" to the end. This way, even if a low
	//   energy cluster were before a high energy one, we would not miscount!
	//
	// "ifront" movers from the 1st clusters, "jend" moves from the last cluster
	// 2008-03-19 SH generalizes to "move to end" clusters with too low energy OR to many hits
	Int_t ifront = 0;
	Int_t jend   = nClusts - 1 ;
	while(1) 
	  {
	    
	    // break out the loop when "ifront" meets "jend"
	    //
	    if( ifront > jend ) 
	      {
		break;
	      }
	    
	    if( clust[ifront].energy >= minRealClusterEne 
		&& clust[ifront].tow->GetEntries()<= maxHitsInRealCluster ) 
	      {
		ifront ++;
	      }
	    else 
	      {
		// search from "jend" and forward, until find a clusters that is higher than "minRealClusterEne"
		//
		while(jend>=ifront && ( clust[jend].energy < minRealClusterEne 
							     || clust[jend].tow->GetEntries()> maxHitsInRealCluster )) {
		  jend --;
		}
		//
		// check that "jend" is in front of "iend", we have exhausted the array. Stop!
		//
		if( jend < ifront ) {
		  break;
		}
		//
		// exchange the clusters: "ifront"<-->"jend", then increase "ifront" by 1 and decrease "jend" by 1;
		//
		else {
		  altClu        = clust[ifront] ;
		  clust[ifront] = clust[jend]   ;
		  clust[jend]   = altClu        ;
		  ifront ++;
		  jend   --;		
		}
		
	      }
	  }
	//
	// at the end of the above code, "ifront" counts the number of clusters above "minRealClusterEne"
	//
	nRealClusts = ifront ;

	// do mement analysis before catagozie!
	//
	//cout<<"test nRealClusts="<<nRealClusts<<endl;
	if(pwe->YPrintLevel>0)
	  {
	    printf("Yiqun:Post Increment    ifront=%d jend=%d minclE=%f maxHitsClu=%d \n",ifront,jend,minRealClusterEne,maxHitsInRealCluster );
	    for(int jcc=0;jcc<nClusts;jcc++)printf("clust[%d]=%d \n",jcc,clust[jcc].numbTower);
	  };
	
	for(Int_t iic=0; iic<nRealClusts; iic++) {

	  clust[iic].FindClusterAxis(pTowerUtil->GetMomentEcutoff());
	}


	// loop over clusters, catagorize, guess the photon locations for cat 0 or 2 clusters
	// then fit, compare, and choose the best fit
	//
	Int_t clustCatag;

	// bad events?
	//
	junkyEvent = false;

	if(pwe->YPrintLevel>0)
	  {
	    printf("Yiqun:nRealClusts=%d\n",nRealClusts);
	  };

	for(Int_t icc=0; icc<nRealClusts; icc++) {

		clustCatag = pTowerUtil->CatagBySigmXY( &clust[icc] );

		// point to the real TObjArray that contains the towers to be fitted
		// it is the same tower array for the cluster or all alternative clusters
		//
		pwe->tow2Fit = clust[icc].tow ;

		// store the goodness of fit
		//
		Double_t chiSq1, chiSq2;

		// Number of Degree of Freedom for the fit
		//
		Int_t ndf;
		ndf = clust[icc].tow->GetEntriesFast();

		// for catag-1 cluster, do 1-photon fit
		if(pwe->NoCatag)clustCatag=pwe->ForceCatag;//from SH	
		//
		if( clustCatag == 1 ) {
			clust[icc].nPhoton   = 1 ;
			chiSq1 = FitOnePhoton(&clust[icc]);

			photons[nPh] = clust[icc].photon[0] ;
			photons[nPh].pHitCluster=&(clust[icc]);
// 			photons[nPh].xPos    = clust[icc].photon[0].xPos;
// 			photons[nPh].yPos    = clust[icc].photon[0].yPos;
// 			photons[nPh].energy  = clust[icc].photon[0].energy;

// 			photons[nPh].errXPos = clust[icc].photon[0].errXPos ;
// 			photons[nPh].errYPos = clust[icc].photon[0].errYPos ;
// 			photons[nPh].errEne  = clust[icc].photon[0].errEne  ;

			ndf -= 3 ;
			if( ndf<=0 )
				ndf = 1;

			clust[icc].chiSquare = chiSq1 / ndf ;

			nPh ++ ;

		}
		//
		// for catagory-2 cluster, do 2-photon fit!
		// If the fit is good enough, it is 1-photon. Else also
		// try 2-photon fit, and find the best fit (including 1-photon fit).
		//
		else if( clustCatag == 2 ) {
			clust[icc].nPhoton   = 2 ;
			chiSq2 = Fit2PhotonClust(&clust[icc]);

			ndf -= 6 ;
			if( ndf<=0 ) 
				ndf = 1;

			clust[icc].chiSquare = chiSq2 / ndf ;

			if( clust[icc].chiSquare > MaxChi2Catag2 ) {
				junkyEvent = true;
			}

			// we need to fill in the information anyway
			//
			clust[icc].nPhoton   = 2 ;
			for(Int_t kp=0; kp<2; kp++) 
			  {
			    
			    photons[nPh] = clust[icc].photon[kp] ;
			    photons[nPh].pHitCluster=&(clust[icc]);
// 				photons[nPh].xPos    = clust[icc].photon[kp].xPos   ;
// 				photons[nPh].yPos    = clust[icc].photon[kp].yPos   ;
// 				photons[nPh].energy  = clust[icc].photon[kp].energy ;
// 				photons[nPh].errXPos = clust[icc].photon[kp].errXPos ;
// 				photons[nPh].errYPos = clust[icc].photon[kp].errYPos ;
// 				photons[nPh].errEne  = clust[icc].photon[kp].errEne  ;
				nPh ++ ;
			}
		}
		//
		// for catagory-0 cluster, first try 1-photon fit!
		// If the fit is good enough, it is 1-photon. Else also
		// try 2-photon fit, and find the best fit (including 1-photon fit).
		//
		else if( clustCatag == 0 ) {
			// 2-photon or 1-photon
			//
			Bool_t is2Photon;
			is2Photon = true;

			altClu = clust[icc] ;
			altClu.nPhoton   = 1 ;
			chiSq1 = FitOnePhoton(&altClu);
			

			Int_t NDF1, NDF2;

			NDF1= ndf - 3 ;
			if( NDF1 <= 0 )
			NDF1 = 1;

			NDF2= ndf - 6 ;
			if( NDF2 <= 0 )
				NDF2 = 1;

			// decide if this 1-photon fit is good enough
			//
			if( chiSq1 < maxGood1PhChi2NDF * NDF1 ) {
				is2Photon = false ;
			}
			//
			// else try 2-photon fit
			//
			else {

				// 2003-09-15
				//
       // New logic for deciding if a 2-photon fit is good!
       // If one photon peak lies on top of low (compared to photon energy) or even zero tower,
       //    this photon is definitely bogus. This could happen if a nearby cluster took away
       //    towers that might make the fit having a large Chi-square (deviation).
       //
       // So for 2-photon fit, first check that the fitted photon of lower energy (higher energy
       //    photon should be fine) is over one of the non-zero tower of the cluster. First of all,
       //    this ensures that we don't have an "outside of cluster" bogus photon, i.e. a bogus
       //    photon could be the result of minimizing the chi-square over towers that do not include
       //    the supposed peak tower. 

			  chiSq2 = Fit2PhotonClust(&clust[icc]);
       // check that the fitted photon of lower energy is contained within one of the non-zero towers
       //    of the cluster
       //
			  Bool_t isBogus2ndPhoton;
			  isBogus2ndPhoton = true;
			  Int_t secondPhoton;
			  if( clust[icc].photon[0].energy < clust[icc].photon[1].energy ) 
			    {
			      secondPhoton = 0;
			    }
			  else 
			    {
			      secondPhoton = 1;
			    }
			  
			  // tower where the fitted photon of lower energy should hit
			  //
			  Int_t ix, iy;
			  ix = 1 + (Int_t) (clust[icc].photon[secondPhoton].xPos / widLG[0]) ;
			  iy = 1 + (Int_t) (clust[icc].photon[secondPhoton].yPos / widLG[1]) ;
			  TowerFPD tow2ndPhoton(clust[icc].photon[secondPhoton].energy, ix, iy, icc) ;
			  
       // now check whether this tower is one of the non-zero towers of the cluster
       //
			  for(Int_t itNZ=0; itNZ<clust[icc].numbTower; itNZ++) {
			    TowerFPD * nZTow = (TowerFPD *) clust[icc].tow->At(itNZ) ;
			    if( tow2ndPhoton.IsEqual( nZTow ) ) {
			      isBogus2ndPhoton = false;
			      tow2ndPhoton.energy = nZTow->energy ;
			      
			      break;
			    }
			  }
			  
			  // and the fitted energy is too large compared to the energy of the tower
			  //
			  if( !isBogus2ndPhoton ) {
			    if( tow2ndPhoton.energy < minHTEneOverPhoton * clust[icc].photon[secondPhoton].energy ) {
			      isBogus2ndPhoton = true;
			    }
			  }
			  
	// 2003-09-29
        // Check that the 2nd photon's "High-Tower" enery is too large compared to its fitted energy.
	//    If so, it is probably splitting one photon into two!
	//
			  if( !isBogus2ndPhoton ) 
			    {
			      Double_t eSS;

			      eSS = EnergyInTowerByPhoton(widLG[0], &tow2ndPhoton, 
							  &( clust[icc].photon[secondPhoton] ) );
			      if( tow2ndPhoton.energy > maxHTEneOverPhoton * eSS ) 
				{
				  isBogus2ndPhoton = true;
				}
			    }
			  
			  // 2003-09-29
			  // Check that the 2nd photon is not near the edge of another cluster.
			  // Namely, we check what would be the energy deposited in other clusters by this photon vs
			  //    energy deposited in its own cluster.
	// If the ratio is too high, this fitted photon is probably a bogus one!
	//

				if( !isBogus2ndPhoton ) {
					Double_t eneWInCluster = EnergyInClusterByPhoton(widLG[0], 
					       &clust[icc], &( clust[icc].photon[secondPhoton] ));

	// loop over all clusters, but skip its own cluster
	//
					for(Int_t kkc=0; kkc<nRealClusts; kkc++) 
					  {
					    if( kkc == icc )
					      continue;
					    
					    if( EnergyInClusterByPhoton(widLG[0], &clust[kkc], 
									&(clust[icc].photon[secondPhoton])) > 
						(maxRatioSpill * eneWInCluster) ) 
					      {
						isBogus2ndPhoton = true;
						break ;
					      }
					  }
				}

				
				if( isBogus2ndPhoton ) 
				  {
				    is2Photon = false ;
				  }
				else 
				  {
				    
				    // 2003-09-12
				    // Need to change!
				    //
				    // What is the criteria of choosing 2-photon fit over 1-photon fit?
				    // How about "chiSq2/NDF > chiSq1/NDF"?
				    //
				    // But there are cases when NDF<=0 !!!
				    //
				    // for catagory 0 cluster and when 1-photon fit is better
				    //
				    if( (chiSq1/NDF1) < (chiSq2/NDF2) ) {
				      is2Photon = false ;
				    }
				    else {
				      is2Photon = true ;
				    }
				    
				  } // isBogus2ndPhoton == true
				
			} // try 2-photon fit
			
			// now fill in the fit result
			//
			
			if( !is2Photon ) 
			  {
			    
			    // 1-photon fit
			    //			
			    clust[icc].nPhoton   = altClu.nPhoton   = 1 ;
			    clust[icc].chiSquare = altClu.chiSquare = chiSq1 / NDF1 ;
			    
			    photons[nPh] = clust[icc].photon[0] =  altClu.photon[0] ;
			    photons[nPh].pHitCluster=&(clust[icc]);
			    
//    			photons[nPh].xPos    = clust[icc].photon[0].xPos    = altClu.photon[0].xPos   ;
//   			photons[nPh].yPos    = clust[icc].photon[0].yPos    = altClu.photon[0].yPos   ;
// 			photons[nPh].energy  = clust[icc].photon[0].energy  = altClu.photon[0].energy ;
// 			photons[nPh].errXPos = clust[icc].photon[0].errXPos = altClu.photon[0].errXPos;
//			photons[nPh].errYPos = clust[icc].photon[0].errYPos = altClu.photon[0].errYPos;
// 			photons[nPh].errEne  = clust[icc].photon[0].errEne  = altClu.photon[0].errEne ;
			    nPh ++ ;
			  }
			else 
			  {
			    // 2-photon fit is better
			    //
			    clust[icc].nPhoton   = 2 ;
			    clust[icc].chiSquare = chiSq2 / NDF2 ;
			    
			    // check if this cluster has too much stuff in it!
			    //
			    if( clust[icc].chiSquare > MaxChi2Catag2 ) {
			      junkyEvent = true;
			    }
			    
			    
			    for(Int_t kp=0; kp<2; kp++) {
			      
			      photons[nPh] = clust[icc].photon[kp] ;
			      photons[nPh].pHitCluster=&(clust[icc]);
// 					photons[nPh].xPos    = clust[icc].photon[kp].xPos   ;
// 					photons[nPh].yPos    = clust[icc].photon[kp].yPos   ;		      
// 					photons[nPh].energy  = clust[icc].photon[kp].energy ;
// 					photons[nPh].errXPos = clust[icc].photon[kp].errXPos ;
// 					photons[nPh].errYPos = clust[icc].photon[kp].errYPos ;
// 					photons[nPh].errEne  = clust[icc].photon[kp].errEne  ;
					nPh ++ ;
			    }
			} // else ( 2-photon fit is better)
			
		} // catag-0
		//
		// should not happen!
		//
		else {
		  std::cout << "Your logic of catagory is wrong! Something impossible happens! This a catagory-" << clustCatag;
		  std::cout << " clusters! Don't know how to fit it!\n" << "\n";
		}
		
	} // loop over all real clusters

	// 2003-10-04
	// Do a global fit for all events (even for 1 real cluster). Return
	//    Chi2/NDF instead of Chi2 in "FitEvent(...)".
	//
	// 2003-09-08
	// if there are more than 1 cluster, do a global fit!
	//
	if( nPh > MAX_NUMB_PHOTONS ) 
	  {
	    //
	    // myFitter can only do up to "MAX_NUMB_PHOTONS"-photon fit
	    //
	    std::cout << "Can not fit " << nPh << " (more than " << MAX_NUMB_PHOTONS << " photons!" << "\n";
	    chiSqG = -1;
	    return nPh;
	  }
	
	// 2003-09-08
	// for global fit, add all towers from all clusters
	//
	// first clear the TObjArray!
	//
	TObjArray allTow(nTows);
	allTow.Clear();
	Int_t ndfg = 0 ;
	for(Int_t jjc=0; jjc<nRealClusts; jjc++) 
	  {
	    allTow.AddAll(clust[jjc].tow);
	    ndfg += (clust[jjc].tow->GetEntriesFast() - 3 * clust[jjc].nPhoton);
	  }
	
	if( ndfg <= 0 ) 
	  {
	    ndfg = 1;
	  }
	
	pwe->tow2Fit = &allTow ;

	// 2003-10-13
	// only do global fit for 2 or more clusters (2-photon fit for one cluster already has global fit)
	//
	if(pwe->DoGlobal && nRealClusts > 1 ) 
	  {
	    chiSqG = GlobalFit(nPh, nRealClusts, clust) / ndfg ;
	    Int_t iph=0;
	    for(Int_t icl=0;icl<nRealClusts;icl++)
	      {
		
		for(Int_t iclph=0;iclph<clust[icl].nPhoton;iclph++)
		  {
		    
		    photons[iph] = clust[icl].photon[iclph];  
		    photons[iph].pHitCluster=&(clust[icl]);
		    iph++;
		    if(iph>nPh){printf("ERROR total nPh=%d iph=%d \n",nPh,iph);break;};
		  }
	      };
	    if(iph!=nPh){printf("ERROR total nPh=%d iph=%d \n",nPh,iph);};

	  }
	else if( nRealClusts == 1 ) 
	  {
	    chiSqG = clust[0].chiSquare ;
	  }
	else 
	  {
	    // temp revove comment
	    //std::cout << "Something wrong! Should NOT have " << nRealClusts << " real clusters! Error!" << "\n";
	    chiSqG = -1 ;
	  };
	
	return nPh;
}



// Calculate the energy deposit in a cluster by a photon (from shower-shape function).
// In this case, we only need to consider non-zero towers.
//
Double_t Yiqun::EnergyInClusterByPhoton(Double_t widthLG, HitCluster *p_clust, PhotonHitFPD *p_photon)
{
	Double_t eSS = 0;
	for(Int_t it=0; it<p_clust->numbTower; it++) {
		eSS += EnergyInTowerByPhoton(widthLG, (TowerFPD *) p_clust->tow->At(it), p_photon);
	};

	return eSS ;
};


// Calculate the energy deposit in a tower by a photon (from shower-shape function).
//
Double_t Yiqun::EnergyInTowerByPhoton(Double_t widthLG, TowerFPD *p_tower, PhotonHitFPD* p_photon)
{
	Double_t xx, yy, eSS;
	//widthLG no longer used SH

	xx = ( (Double_t) p_tower->col - 0.5 ) * widLG[0] - p_photon->xPos ;
	yy = ( (Double_t) p_tower->row - 0.5 ) * widLG[1] - p_photon->yPos ;
	eSS = p_photon->energy * (pwe->fcnSS)->Eval( xx, yy, 0 );

	return eSS;
};
Yiqun::Yiqun(TMatrix* pEm,Geom* pgeom,Int_t iew,Int_t nstb)
{
  p_geom=pgeom;
  p_gain=0;
  p_gaincorr=0;
  EW=iew;
  NSTB=nstb;
  Y(pEm);
};
Yiqun::Yiqun(TMatrix* pEm,Geom* pgeom,CalibStr* pgain,CalibStr* pgaincorr,Int_t iew,Int_t nstb)
{  
  p_geom=pgeom;
  p_gain=pgain;
  p_gaincorr=pgaincorr;
  EW=iew;
  NSTB=nstb;
  Y(pEm);
};
void Yiqun::Y(TMatrix* pEm)
{   
  NPh=0;
  NTower=0;
  ChiSqG=10000.;
  NClusts=0;
  NRealClusts=0;
  pwe=&(FitTower::we);
  
  ShadowClusters=new TObjArray(20,0);
  //  pTowerUtil=new TowerUtil();
  pTowerUtil=new TowerUtil(pwe);
  pTowerUtil->SetMomentEcutoff(.5);  
  Clu2=0;
  tow_Arr=0;
  p_LastFit=NULL;
  choiceChi2=2;  
  if(pwe->dev)delete pwe->dev;
  if(pwe->dchi2)delete pwe->dchi2;
  widLG[0]=widLG[1]=(p_geom->FpdTowWid(EW,NSTB))[0];
  if(p_geom->FMSGeom)widLG[1]=(p_geom->FpdTowWid(EW,NSTB))[1];
  //pwe->myRand=NULL;

  /*
  if( !(pwe->myRand) ) {
    std::cout << "Initialize random generator (for Z_gg and d_gg)!" << "\n";
    (pwe->myRand) = new TRandom();
  }
  */

  pwe->choiceChi2= choiceChi2;
  NTower=pEm->GetNoElements();

  if(NTower>578)
    {
      printf("Too many towers for Fit\n");
      exit(-1);
    };

  nrows=pEm->GetNrows();
  ncols=pEm->GetNcols();
  pwe->dev=new TMatrix(*pEm);
  pwe->dchi2=new TMatrix(*pEm);
  *(pwe->dev)=0.;
  *(pwe->dchi2)=0.;

  Int_t cnt=-1;
  tow_Arr=new TObjArray(NTower);

  for(Int_t ir=0;ir<nrows;ir++)
    {
      for(Int_t ic=0;ic<ncols;ic++)
	{
	  cnt++;  
	  towers[cnt].energy  = (*pEm)(ir,ic);
	  towers[cnt].col     = ic+1 ;
	  towers[cnt].row     = ir+1 ;
	  towers[cnt].cluster = -1 ; 
	  
	  if((p_gain!=0) && (p_gaincorr!=0))
	    {
	      towers[cnt].adc_over_ped=0;
	      Float_t fullgain=p_gain->GetValue(EW,NSTB,ir,ic)*(p_gaincorr->GetValue(EW,NSTB,ir,ic));
	      if(fullgain>0)towers[cnt].adc_over_ped=(Int_t) (towers[cnt].energy/fullgain+.5);
	      if(towers[cnt].adc_over_ped<1.)towers[cnt].energy=0.;
	    };
	  if(towers[cnt].energy>.001)tow_Arr->Add(&towers[cnt]);
	};
    }

  MakeClu2();  // get second cluster code
 
  /* Seems not used (SH)
  for(Int_t ic=0; ic<MAX_NUMER_CLUSTERS; ic++) {
    clust[ic].tow = &tclu[ic] ;
  }
  */
  PRINT_FIT_1_RESULT=false;
  PRINT_FIT_2_RESULT=false;
  PRINT_FIT_ALL_RESULT=false;
  PRINT_FIT_PARA=false;

  // choice of Chi-square function
  // 1:  Steve's
  // 2:  Larisa's
  //
  
  // clusters with just 1 photon
  //

  posDif_1PC=0.5;
  
  eneRat_1PC=0.15;
  posDif_2PC=0.2;
  eneRat_2PC=0.05;
  dggPara[0]=18.0;
  dggPara[1]=2.2;
  dggPara[2]=0.5;
  dggPara[3]=60.0;
  dggPara[4]=0.085;
  dggPara[5]=3.5;
  
  thetaPara=2.8;
  
  // global fit (more than 1 cluster)
  //

  posDif_Gl=1.25;
  eneRat_Gl=0.3;
 if(!pwe->UseThis_Err)
    {
      pwe->errFactor=0.03;
      pwe->errQ=0.01;
    };
  maxGood1PhChi2NDF=5;
  minHTEneOverPhoton=0.25;
  maxHTEneOverPhoton=1.5;
  maxRatioSpill=0.2;
  MaxChi2Catag2=10.;
  minRealClusterEne=2.0;
  maxHitsInRealCluster=49;
  if(EW==2 && (NSTB==1 || NSTB==2))
    {
      maxHitsInRealCluster=10;
      minRealClusterEne=.75;
      //      maxHitsInRealCluster=4;
      maxHitsInRealCluster=25;
    };
  pwe->useLimitForPara=1;
  pwe->showerShapeFunc=0; 
  NClusts=0;
  NRealClusts=0;
  ChiSqG=0.;
  JunkyEvent=false;
  fitter=0;
  if(pwe->BlockFit)return;

  fitter = new FitTower(pEm,p_geom,EW,NSTB);
  pwe->showerWidthX=1.0;
  pwe->showerWidthY=1.0;

  fitter->SetTWidthCM(*(p_geom->FpdTowWid(EW,NSTB)));
  if(p_geom->FMSGeom)fitter->SetXYTWidthCM(p_geom->FpdTowWid(EW,NSTB));
  
 NPh = FitEvent(NTower, NClusts, NRealClusts, ChiSqG, JunkyEvent);
};
Int_t Yiqun::TryGlobal()
{
  HitCluster* twoc;
  // loop through cluster list and find 2 photon clusters
  for(int icl=0;icl<NRealClusts;icl++)
    {
      twoc=&clust[icl];
      if(twoc->nPhoton==2)
	{
	  printf("Before Global FiQt for cluster %d.0: E=%f x=%f y=%f \n",icl,
		 twoc->photon[0].energy,twoc->photon[0].xPos,twoc->photon[0].yPos);
	  printf("                                %d.1: E=%f x=%f y=%f \n",icl,
		 twoc->photon[1].energy,twoc->photon[1].xPos,twoc->photon[1].yPos);
	  //	  GlobalFit(2,1,twoc);
	  printf("After Global Fit for cluster %d.1: E=%f x=%f y=%f \n",icl,
		 twoc->photon[0].energy,twoc->photon[0].xPos,twoc->photon[0].yPos);
	  printf("                                %d.1: E=%f x=%f y=%f \n",icl,
		 twoc->photon[1].energy,twoc->photon[1].xPos,twoc->photon[1].yPos);
	};
    };
  return 0;
};
void Yiqun::MakeClu2()
{
  if(!tow_Arr)return;  
  ((TowerFPD*) tow_Arr->First())->ClearClusterInfo(tow_Arr);
  TIter next(tow_Arr);
  while(TowerFPD* tow=(TowerFPD*) next())
    {
      tow->SetContext(tow_Arr,p_gain,p_gaincorr,EW,NSTB);
    };
  
  Int_t NC=0;
  Int_t NT=0;
  if((NT=tow_Arr->GetEntries())>0)
    {
      NC=((TowerFPD*) tow_Arr->First())->CreateContiguous();
    };
  
  Clu2=new TObjArray(NC,0);
  
  if(NT>0)
    {
      TowerFPD* FirstTow=(TowerFPD*) tow_Arr->First();
      Int_t clustercnt=0;
      for(int j=0;j<NC;j++)
	{	  
	  
	  TObjArray* tar=FirstTow->CreateClusterList(j);
	  if(tar==0)continue;
	  HitCluster* hcl=new HitCluster();
	  
	  hcl->tow=tar;
	  hcl->energy=0;
	  Int_t clusterlength=hcl->tow->GetEntries();
	  Float_t Energysum=0;
	  if(clusterlength>0)
	    {
	      hcl->numbTower=clusterlength;
	      hcl->catag=((TowerFPD*) tow_Arr->First())->ClusterCategory(hcl->tow);
	      Clu2->AddAt(hcl,clustercnt);
	      hcl->cluster2=clustercnt;
	      hcl->UpdateEnergy();
	      hcl->Setcluster2Number(clustercnt);
	      clustercnt++;
	    };
	};
    };
  
  next.Reset();
  
  /*  ------------------------
      while(TowerFPD* tow=(TowerFPD*) next())
      {
      if(tow->LocalHighTower())
      {
      printf("High Tower NSTB=%d,row=%d col=%d ,E=%f, ClaimedID=%d\n",
      tow->NSTB,tow->row,tow->col,tow->energy,tow->ClaimedID);
      }
      else
      {
      printf("Other Tower NSTB=%d,row=%d col=%d ,E=%f, ClaimedID=%d\n",
      tow->NSTB,tow->row,tow->col,tow->energy,tow->ClaimedID);
      };
      };
      
      TIter nextC(Clu2);
      Int_t knt=0;
      while( HitCluster* cl=(HitCluster*) nextC())
      {
      printf("print cluster %d \n",knt);
      cl->Print();
      knt++;
      }
            -------------- */
  
  
};

void Yiqun::PrintClu()
{
  printf("PrintClu()----------- %d photons   %d clusters\n ",NPh,NClusts); 
  for(Int_t i=0;i<NPh;i++)
    {
      printf("photons[%d]:\n",i);
      photons[i].Print();
    };
  for(Int_t i=0;i<NClusts;i++)
    {
      printf("clust[%d]:\n",i);
      clust[i].Print();
      if(clust[i].nPhoton>1)
	{
	  Float_t m12=(mom(&clust[i].photon[0])+mom(&clust[i].photon[1])).Mag();
	  printf("Mass ph0-ph1=%f\n",m12);
	};
    };
};

TMatrix Yiqun::FittedMat()
{
  TMatrix mat(nrows,ncols);
  mat=0.;
  TF2* SSFunct=fitter->GetFunctShowShape();
  Double_t lgwidth;
  fitter->GetTWidthCM(lgwidth);
  for(Int_t ir=0;ir<nrows;ir++)
    {
      for(Int_t ic=0;ic<ncols;ic++)
	{
	  Float_t energySS=0;
	  for(Int_t inp=0;inp<NPh;inp++)
	   {
	     //must change + back to - ....
	     Double_t xx=(((Double_t) ic)+.5)*widLG[0] - photons[inp].xPos;
	     Double_t yy=(((Double_t) ir)+.5)*widLG[1] - photons[inp].yPos;
	     
	     energySS+= photons[inp].energy * (*fitter->GetFunctShowShape()).Eval( xx, yy,0 );
	   };
	  mat(ir,ic)=0;
	  if((p_gain!=0) && (p_gaincorr!=0))
	    {
	      Float_t g1=p_gain->GetValue(EW,NSTB,ir,ic);
	      Float_t g2=p_gaincorr->GetValue(EW,NSTB,ir,ic);
	      if(g1*g2<=0)
		{
		  mat(ir,ic)=.001;
		}
	      else
		{
		  mat(ir,ic)=energySS-fmod(energySS,g1*g2);
		}	      
	    }
	  else
	    {
	      if(energySS>.01)mat(ir,ic)=energySS;
	    };
	};
    };
  return mat;
};

void Yiqun::FittedMat(TMatrix* mat)
{

  TF2* SSFunct=fitter->GetFunctShowShape();
  Double_t lgwidth;
  fitter->GetTWidthCM(lgwidth);
  for(Int_t ir=0;ir<nrows;ir++)
    {
      for(Int_t ic=0;ic<ncols;ic++)
	{
	  Float_t energySS=0;
	  for(Int_t inp=0;inp<NPh;inp++)
	   {
	     //must change + back to - ....
	     Double_t xx=(((Double_t) ic)+.5)*widLG[0] - photons[inp].xPos;
	     Double_t yy=(((Double_t) ir)+.5)*widLG[1] - photons[inp].yPos;
	     
	     energySS+= photons[inp].energy * (*fitter->GetFunctShowShape()).Eval( xx, yy,0 );
	   };
	  (*mat)(ir,ic)=0;
	  if((p_gain!=0) && (p_gaincorr!=0))
	    {
	      Float_t g1=p_gain->GetValue(EW,NSTB,ir,ic);
	      Float_t g2=p_gaincorr->GetValue(EW,NSTB,ir,ic);
	      if(g1*g2<=0)
		{
		  (*mat)(ir,ic)=.001;
		}
	      else
		{
		  (*mat)(ir,ic)=energySS-fmod(energySS,g1*g2);
		}	      
	    }
	  else
	    {
	      if(energySS>.01) (*mat)(ir,ic)=energySS;
	    };
	};
    };

};
TVector3 Yiqun::ph_coord_lab(PhotonHitFPD* phot)
{
  TVector3 xyz;
  xyz[2]=*(p_geom->ZFPD(EW,NSTB));
  xyz[0]=phot->xPos;
  xyz[1]=phot->yPos;
  if(NSTB==1 || (EW==2 && NSTB==3))
    {
      xyz[0]=(*(p_geom->xOffset(EW,NSTB)))-xyz[0];
    }
  else
    {
      xyz[0]=(*(p_geom->xOffset(EW,NSTB)))+xyz[0];
    };
  xyz[1]=*(p_geom->yOffset(EW,NSTB))-xyz[1];
  return xyz;
};
TVector3 Yiqun::ph_coord_lab(Int_t ph_num)
{
  if(ph_num>= NPh)
    {
      printf("request more photons than are defined \n");
      TVector3 tmp(0,0,0);
      return tmp;
    };
  TVector3 xyz;
  xyz[2]=*(p_geom->ZFPD(EW,NSTB));
  xyz[0]=photons[ph_num].xPos;
  xyz[1]=photons[ph_num].yPos;
  if(NSTB==1 || (EW==2 && NSTB==3))
    {
      xyz[0]=(*(p_geom->xOffset(EW,NSTB)))-xyz[0];
    }
  else
    {
      xyz[0]=(*(p_geom->xOffset(EW,NSTB)))+xyz[0];
    };
  xyz[1]=*(p_geom->yOffset(EW,NSTB))-xyz[1];
  return xyz;
};

TLorentzVector  Yiqun::mom(Int_t ph_num)
{
  if(ph_num>= NPh)
    {
      printf("request more photons than are defined \n");
      TLorentzVector tmp(0.,0.,0.,0.);
      return tmp;
    };  
  TVector3 vvec=ph_coord_lab(ph_num);
  Double_t dist=vvec.Mag();
  TVector3 uvec(0.,0.,0.);
  if(dist!=0)uvec=(1./dist)*vvec;
  Float_t phEnergy=photons[ph_num].energy;
  if(pwe->UseEDepCorrection)
    {
      if(NSTB<3)
	{
	  phEnergy=phEnergy*pwe->EDepCorrection->Eval(20)/pwe->EDepCorrection->Eval(phEnergy);
	}
      else
	{
	  phEnergy=phEnergy*pwe->EDepCorrection->Eval(30.)/pwe->EDepCorrection->Eval(phEnergy);
	};
    };
  TVector3  mom3;
  mom3=phEnergy*uvec;
  TLorentzVector mom4(mom3,phEnergy);
  return mom4;
};
TLorentzVector  Yiqun::mom(PhotonHitFPD* phot)
{
  TVector3 vvec=ph_coord_lab(phot);
  Double_t dist=vvec.Mag();
  TVector3 uvec(0.,0.,0.);
  if(dist!=0)uvec=(1./dist)*vvec;
  Float_t phEnergy=phot->energy;
  TVector3  mom3;
  mom3=phEnergy*uvec;
  TLorentzVector mom4(mom3,phEnergy);
  return mom4;
};
void Yiqun::Print()
{
  printf("ShowerShape=%d\n",pwe->showerShapeFunc);
  printf("ShowerWidthX=%f\n",pwe->showerWidthX);
  printf("ShowerWidthY=%f\n",pwe->showerWidthY);
  
  
  
  std::cout << "choiceChi2 = " << choiceChi2 << "\n";
  std::cout << "posDif_1PC = " << posDif_1PC << "\n";
  std::cout << "eneRat_1PC = " << eneRat_1PC << "\n";
  
  std::cout << "posDif_2PC = " << posDif_2PC << "\n";
  std::cout << "eneRat_2PC = " << eneRat_2PC << "\n";
  for(Int_t ii=0; ii<6; ii++) 
    {
      std::cout << "dggPara[" << ii << "] = " << dggPara[ii] << "\n";
    }
  std::cout << "thetaPara = " << thetaPara << "\n";
  
  std::cout << "posDif_Gl = " << posDif_Gl << "\n";
  std::cout << "eneRat_Gl = " << eneRat_Gl << "\n";
  
  std::cout << "errFactor = " << pwe->errFactor << "\n";
  std::cout << "errQ = " << pwe->errQ << "\n";
  std::cout << "maxGood1PhChi2NDF = " << maxGood1PhChi2NDF << "\n";
  std::cout << "minHTEneOverPhoton = " << minHTEneOverPhoton << "\n";
  std::cout << "maxHTEneOverPhoton = " << maxHTEneOverPhoton << "\n";
  std::cout << "maxRatioSpill = " << maxRatioSpill << "\n";
  std::cout << "MaxChi2Catag2 = " << MaxChi2Catag2 << "\n";
  std::cout << "minRealClusterEne = " << minRealClusterEne << "\n";
  std::cout << "useLimitForPara = " << pwe->useLimitForPara << "\n";
};

Yiqun::~Yiqun()
{
  if(fitter!=0)delete fitter;

  if(p_LastFit!=0)delete p_LastFit;

  
  if(pTowerUtil)delete pTowerUtil;
 
  TIter nxt(Clu2);
  while(HitCluster* ccl=((HitCluster*) nxt()))
    {
      delete (ccl->tow);
      (ccl->tow)=0;
    }
  if(Clu2)delete Clu2;
  if(tow_Arr)delete tow_Arr;
  if(ShadowClusters)
    {
      TIter nextsh(ShadowClusters);
      HitCluster* hc;
      while(hc=(HitCluster*) nextsh())delete hc;
      delete ShadowClusters;
    };
};
