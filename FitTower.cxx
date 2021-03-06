#include "FitTower.h"

ClassImp(FitTower)

FitTower::FitTower(TMatrix* pEm,Geom* pgeom,Int_t iew,Int_t nstb)
{
  SetStep();
  pTowerUtil=new TowerUtil();
  fCol = pEm->GetNcols() ;
  fRow = pEm->GetNrows() ;
  fTWidthCM=*(pgeom->FpdTowWid(iew,nstb));
  if(!we.UseThis_ab)
    {
      
      Setwe_ab(.8,.3,-.1,.8,.2,7.6);
    };

  Int_t numbPara = 10;
  Double_t para[numbPara];
  
  para[0] = fTWidthCM ;
  para[1] =  we.a1 ;
  para[2] =  we.a2 ;
  para[3] = we.a3 ;
  //SH preferred
  //  para[2] =  0.36 ;
  //  para[3] = -0.16 ;
  para[4] =  we.b1 ;
  para[5] =  we.b2 ;
  para[6] =  we.b3 ;
  para[7] =  0.0 ;
  para[8] =  0.0 ;
  para[9] =  1.0 ;
  we.showerWidthX=1.0;
  we.showerWidthY=1.0;
  if(we.fcnSS==0)
    {
      we.fcnSS = new TF2("fFunctSS", &GGams, -25.0, 25.0, -25.0, 25.0, numbPara);

      we.fcnSS->SetParameters(para); 
      fFunctSS=new TF2(*(we.fcnSS));
    }
  else
    {
      we.fcnSS->SetParameters(para); 
      fFunctSS=new TF2(*(we.fcnSS));
    };
  
  // copy the global Shower-shape function pointer to the local function
  //
  // create a Minuit instance
  //

  fMn = new TMinuit(3*MAX_NUMB_PHOTONS+1);

};

void FitTower::SetSubclu2(Bool_t useclu2)
{
  we.SubClu2=useclu2;
};
void FitTower::SetFreeGlobals(Double_t* fglobals,int dim)
{
  if(dim>0 && dim<=10)
    {
      for(int j=0;j<dim;j++)we.FreeGlobals[j]=fglobals[j];
    };
};
Bool_t FitTower::Setwe_EDepCor(Bool_t useEdepCor)
{
  we.UseEDepCorrection=useEdepCor;
  return true;
};
Bool_t  FitTower::SetForceMass(Float_t fmass)
{
  we.Force2Mass=fmass;
  return true;
}  
Bool_t FitTower::SetDoGlobal(Bool_t sdg)
{
  //bool to true "do" or false "dont do" global fit
  we.DoGlobal=sdg;
  return we.DoGlobal;
};
void FitTower::SetYPrintLevel(int PLevel)
{
  we.YPrintLevel=PLevel;
}
Bool_t FitTower::SetNoCatag(Bool_t value,Int_t cat)
{
  //NoCatag=false implies no effect 
  //NoCat==true implies catag fixed to cat
  we.NoCatag=value;
  we.ForceCatag=cat;
  return value;
}
Bool_t FitTower::Setwe_ErrFactors(float errQ,float errFactor,float p1,float p2,float energy_study)
{
  we.errQ=errQ;
  we.errFactor=errFactor;
  we.UseThis_Err=true;
  we.Power1=p1;
  we.Power2=p2;
  we.Energy_study=energy_study;
};

void FitTower::SetBlockFit(Bool_t bl)
{
  we.BlockFit=bl;
}

Bool_t FitTower::Setwe_ab(Float_t a1,Float_t a2,Float_t a3,Float_t b1,Float_t b2,Float_t b3)
{
  we.a1=a1;
  we.a2=a2;
  we.a3=a3;
  we.b1=b1;
  we.b2=b2;
  we.b3=b3;
  we.UseThis_ab=true;
  return true;
};

/*
FitTower::FitTower(const Int_t dim[2], const Double_t wd, TF2 *ssFunct)
{
  
  SetStep();
  pTowerUtil=new TowerUtil();
  fCol = dim[0] ;
  fRow = dim[1] ;
  fTWidthCM = wd ;
  
  
  fFunctSS = ssFunct ;

  we.fcnSS = fFunctSS ;
  
  // create a Minuit instance
  //
  fMn = new TMinuit(3*MAX_NUMB_PHOTONS+1);  
}
*/
FitTower::~FitTower()
{
  
  if( fMn ) 
    {
      delete fMn;
    };
  if(fFunctSS)delete fFunctSS;
  
  delete pTowerUtil;
};

void FitTower::SetStep()
{
  fMn=0;

  for(int j=0;j<3*MAX_NUMB_PHOTONS+1;j++)step[j]= we.step[j];
};

Double_t FitTower::FGams(Double_t *x, Double_t *para)
{
  
  Double_t f=0;
  Double_t xx=x[0]/we.showerWidthX;
  Double_t yy=x[1]/we.showerWidthY;
  for(Int_t i=1;i<=3;i++) {
    Int_t j;
    j = i + 3 ;
    f += para[i]*( atan( xx*yy/ (para[j]* sqrt(para[j]*para[j]+xx*xx+yy*yy) ) ) );
  };
  return f/(2 * TMath::Pi() ) ;
}

Double_t FitTower::GGams(Double_t *x, Double_t *para)
{
  Double_t gg, s[2];
  gg = 0 ;  
  for(Int_t ix=0; ix<2; ix++) {
    for(Int_t iy=0; iy<2; iy++) {
      Double_t ax, ay;
      ax = pow(-1.0, ix);
      ay = pow(-1.0, iy);
      s[0] = x[0] - para[7] + ax * para[0] / 2.0 ;
      s[1] = x[1] - para[8] + ay * para[0] / 2.0 ;
      gg += ax * ay * FGams(s, para) ;
    }
  }
  return gg*para[9];
}

void FitTower::Fcn1(Int_t& npara, Double_t* grad,  Double_t& fval, Double_t* para, Int_t iflag)
{
  *(we.dev)=0;
  *(we.dchi2)=0;
  // number of expected photons
  // should ALWAYS be the first parameter "para[0]"
  //
  Int_t numbPh;
  numbPh = (Int_t) para[0] ;
  
  fval = 0 ;
  Double_t err ;
  Double_t xx, yy;
  Double_t eSS, eMeas;
  Double_t dev;
  
  // we are in GeV, not ADC count
  //
  //   Double_t we.errFactor = 0.15 ;
  
  TowerFPD * oneTow;
  
  
  // first get cluster sum
  //
  Double_t sumCl = 0;
  
  if( we.choiceChi2 == 2 ) 
    {
      TIter next(we.tow2Fit);
      while(oneTow=(TowerFPD*) next())sumCl+=oneTow->energy;
    }
  
  // loop over all towers that are involved in the fit
  //
  TIter next(we.tow2Fit);
  while(oneTow=(TowerFPD*) next())
    {
      
    // center of tower in unit of "cm"
    //
    // my towers are center at 0.5 to 6.5, as Steve Heppelmann
    //
    //Note from SFH need more clever position for FMS

    xx = ( (Double_t) oneTow->col - 0.5 ) * we.widLG[0] ;
    yy = ( (Double_t) oneTow->row - 0.5 ) * we.widLG[1] ;
    
    // measured energy
    //
    eMeas = oneTow->energy;
    
    // expected energy from Shower-Shape
    //
    eSS = 0 ;

    for(Int_t iph=0; iph<numbPh; iph++) {
      Int_t j;
      j = 3 * iph ;
      
      //2004 feb 25, added by akio
      //If showerShapeFunc>0, it means we choose energy dependent
      //shower shape. "measured" photon energy is in para[j+3].
      //all we need to do is set the parameters in fncSS:
      if(we.showerShapeFunc==1) 
	{
	  double e=para[j+3];			  
	  //guess Egamma from Emeasured
	  //  double e_measured=e;
	  //  double e_gamma=0.6022+1.192*e_measured;
	  //  double e=egamma;
	  //get shower shape parameters as function of Egamma			  
	  double loge=log(e);
	  double a0 = 0.5728 + 0.07085*loge;
	  double a1 = 0.01169 + 0.09739*loge;
	  double a2 = 0.6311 - 0.02977*loge;
	  double a3 = a0-a1-a2;
	  double b1 = 2.512 - 0.4304*loge;
	  double b2 = 0.8271 - 0.1500*loge;
	  double b3 = 19.95 + 2.282*loge;
	  //re-normarize to a1+a2+a3=1
	  a1=a1/a0;
	  a2=a2/a0;
	  a3=a3/a0;
	  //set parameters
	  we.fcnSS->SetParameter(1,a1);
	  we.fcnSS->SetParameter(2,a2);
	  we.fcnSS->SetParameter(3,a3);
	  we.fcnSS->SetParameter(4,b1);
	  we.fcnSS->SetParameter(5,b2);
	  we.fcnSS->SetParameter(6,b3);
	  //printf("SS: E=%f a=%f %f %f %f b=%f %f %f\n",e,a0,a1,a2,a3,b1,b2,b3);
	}
      
      //
      // shower-shape function calculate the fraction of energy
      // in coords of center of tower relative to photon
      //
      
      Double_t Eshape = para[j+3] * we.fcnSS->Eval(xx-para[j+1], yy-para[j+2], 0);
      //	    LastShape.Fill(oneTow->col - 0.5,oneTow->row - 0.5,Eshape);
      eSS+=Eshape;
    }
    
    
    dev = eMeas - eSS ;
    //		printf("inFcn  col=%d row=%d eMeas=%f: ",oneTow->col,oneTow->row,eMeas);
    //		printf("fitted eSS= %f \n",eSS);
    Double_t dchi2;
    (*(we.dev))[oneTow->row-1][oneTow->col-1]=dev;

    if( we.choiceChi2 == 2 ) {
      //
      // Larisa'e Chi2 function
      //
      
      if(!we.UseThis_Err)
	{
	  err = we.errFactor * eMeas * (1 - eMeas/sumCl) + we.errQ ;
	}
      else
	{
	  Float_t power1=we.Power1;
	  Float_t power2=we.Power2;
	  err = (we.errFactor * pow(eMeas/sumCl,power1-.001*sumCl) * 
		 pow(1 - eMeas/sumCl,power2-.007*sumCl))*sumCl
	    +we.errQ;
	};
      dchi2 = dev * dev / err ;
      float dsign=1.;
      if(dev<0)dsign=-1.;
      (*(we.dchi2))[oneTow->row-1][oneTow->col-1]=dsign*dchi2;
    }
    else if( we.choiceChi2 == 1 ) {
      //
      // Steve's Chi2 function
      //
      // 2003-09-11
      //
      // translate into unit of GeV
      //
      err = we.errFactor * sqrt( eMeas ) ;
      if( err < 0.25 )
	err = sqrt( err * err + 0.25 ) ;
      
      // 2003-09-11
      // do not understand yet!
      //
      if( eSS < 1.25 )
	err = sqrt( eSS + eMeas + 3 / 4 );
      
      dchi2 = dev * dev / err / err ;
    }
    else {
      std::cout << "Your \"choiceChi2\" = " << we.choiceChi2 << " is invalid! Quit!!!" << "\n";
      exit(-1);
    }
    
    fval += dchi2 ;
    
  };
  
  
  // 2003-09-14
  // require that the fraction be positive!
  //
  if( fval < 0 )
    fval = 0;
  
  
  //std::cout<<f<<"\n";
  
};


// 2003-08-29
// a different set of parameters for 2-photon clusters only:
//
//    param[0]:  still a constant parameter, should be set to 2 for 2-photon fitting
//    param[1]:  xPi      (x-position of pi^0)
//    param[2]:  yPi      (y-position of pi^0)
//    param[3]:  d_gg     (distance between 2 photons)
//    param[4]:  theta    (theta angle of displacement vector from photon 2 to photon 1)
//    param[5]:  z_gg     (this z_gg can go from -1 to +1, so we do not set E1>E2)
//    param[6]:  E_gg     (total energy of two photons)
//
// Thus, in more conventional parameterization: x1, y1, E1, x2, y2, E2:
//
//     E1 = E_gg * (1 + z_gg)/2
//     E2 = E_gg * (1 - z_gg)/2
//     x1 = xPi + cos(theta) * d_gg * (1 - z_gg)/2
//     y1 = yPi + sin(theta) * d_gg * (1 - z_gg)/2
//     x2 = xPi - cos(theta) * d_gg * (1 + z_gg)/2
//     y2 = yPi - sin(theta) * d_gg * (1 + z_gg)/2
//
// The advantage of the new parameterization is that for 2-photon cluster fitting, we can
//    ensure that the two photons never get to close. The old parameterization certainly
//    suffers from this shortcoming if we let the parameters vary freely.
//
// What we already know about the limits of the new parameters:
//
//    xPi and yPi:   rarely do they go beyond 0.3 unit of lgd
//    theta:         have a narrow theta range (for r=sigmaMin/sigmaMax, |theta|<0.5*r/0.65
//                      when r<0.65, and linear increase from 0.5 to Pi/2 for 0.65<r<1)
//    E_gg:          given by Ec (+/- 20% or less)
//    z_gg:          should just let it vary from -1 to 1.
//    d_gg:          a lower bound is given by r=sqrt(sigmaX^2+sigmaY^2). 
//                      d_gg > Max( 2.5*(r-0.6), 0.5 )
//
//

void FitTower::SetFCN(void (*fcn)(Int_t &, Double_t *, Double_t &, Double_t *, Int_t))
{
  fMn->SetFCN(fcn);
}

void FitTower::SetNumberPhoton(const Int_t nP)
{
  if( nP < 1 || nP > MAX_NUMB_PHOTONS ) {
    fNumbPhotons = 1;
    std::cerr << "nP = " << nP << "! Number of photons must be between 1 and " << MAX_NUMB_PHOTONS << "! Set it to be 1!" << "\n";
  }
  else {
    fNumbPhotons = nP ;
  }
  
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter!
  //
  Int_t ierflg = 0;
  fMn->mnparm(0, "nph", fNumbPhotons, 1.0, 1.0, 4.0, ierflg);
  fMn->FixParameter(0);
}


Int_t FitTower::Fit(const Double_t *para, const Double_t *step, const Double_t *low, const Double_t *up)
{
  
  // check that there is a pointer to TObjArray of towers
  //
  if( !(we.tow2Fit) ) {
    std::cerr << "no tower data available! return -1!" << "\n";
    return -1;
  }
  
  
  // 	std::cout << "do fit!" << "\n";
  if(we.YPrintLevel>1)
    {
      fMn->SetPrintLevel(we.YPrintLevel-1);}
  else  fMn->SetPrintLevel(-1);
  fMn->fLwarn = false ;
  // 	fMn->fLwarn = true ;
  
  // must set the function to "Fcn1"!
  //
  fMn->SetFCN(Fcn1);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  
  fMn->mnexcm("SET ERR", arglist, 1, ierflg);
  
  Int_t nPh = (Int_t) para[0];
  
  if( nPh < 1 || nPh > MAX_NUMB_PHOTONS ) {
    fNumbPhotons = 1;
    std::cerr << "nPh = " << nPh << "! Number of photons must be between 1 and " << MAX_NUMB_PHOTONS << "! Set it to be 1!" << "\n";
  }
  else {
    fNumbPhotons = nPh ;
  }
  
  // clear old parameters, so we can define the new parameters
  //
  fMn->mncler();
  
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter!
  //
  ierflg = 0;
  fMn->mnparm(0, "nph", fNumbPhotons, 0, 0.5, 4.5, ierflg);
  // 	fMn->FixParameter(0);
  
  // set the rest of parameters, 3 parameters for 1 photon
  //
  for(Int_t i=0; i<fNumbPhotons; i++) {
    Int_t j ;
    j = 3*i+1 ;
    fMn->mnparm(j, Form("x%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++ ;
    fMn->mnparm(j, Form("y%d", i+1), para[j], step[j], low[j], up[j], ierflg);
    j++ ;
    fMn->mnparm(j, Form("E%d", i+1), para[j], step[j], low[j], up[j], ierflg);
  }
  
  arglist[0] = 1000;
  arglist[1] = 1.;
  ierflg = 0;
  fMn->mnexcm("MIGRAD", arglist ,2,ierflg);
  return fMn->GetStatus() ;
}

Int_t FitTower::Fit2Pin1Clust(const Double_t *para, const Double_t *step, const Double_t *low, const Double_t *up)
{
  fMn->SetPrintLevel(-1);

if(we.YPrintLevel>1)
    {
      printf("Fit2Pin1Clust called\n");
      fMn->SetPrintLevel(+1);
    }
  // check that there is a pointer to TObjArray of towers
  //
  if( !(we.tow2Fit) ) {
    std::cerr << "no tower data available! return -1!" << "\n";
    return -1;
  }
  
    fMn->fLwarn = false ;
  
  // must set the function to "Fcn2"!
  //
  fMn->SetFCN(Fcn2);
  
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  fMn->mnexcm("SET ERR", arglist, 1, ierflg);
  
  Int_t nPh = (Int_t) para[0];
  
  fNumbPhotons = 2;
  
  if( nPh != 2 ) {
    std::cerr << "number of photons must be 2 for special 2-photon cluster fitter \"Int_t FitTower::Fit2Pin1Clust(...)\"!";
    std::cerr << " Set it to be 2!" << "\n";
    fNumbPhotons = 2;
  }
  
  // clear old parameters, so we can define the new parameters
  //
  fMn->mncler();
  
  // The first parameter tells Minuit how many photons to fit!
  // It should be a fixed parameter!
  //
  ierflg = 0;
  fMn->mnparm(0, "nph", fNumbPhotons, 0, 1.5, 2.5, ierflg);
  
  fMn->mnparm(1, "xPi"  , para[1], step[1], low[1], up[1], ierflg);
  fMn->mnparm(2, "yPi"  , para[2], step[2], low[2], up[2], ierflg);
  fMn->mnparm(3, "d_gg" , para[3], step[3], low[3], up[3], ierflg);
  fMn->mnparm(4, "theta", para[4], step[4], low[4], up[4], ierflg);
  fMn->mnparm(5, "z_gg" , para[5], step[5], low[5], up[5], ierflg);
  fMn->mnparm(6, "E_gg" , para[6], step[6], low[6], up[6], ierflg);
  
  arglist[0] = 1000;
  arglist[1] = 1.;
  arglist[2] = 5;
  ierflg = 0;
  
  
  // 2003-10-13
  // fix E_total and theta
  //
  fMn->FixParameter(6);
  fMn->FixParameter(4);
  fMn->mnexcm("CALl", &arglist[2] ,1,ierflg);
  fMn->mnexcm("MIGRAD", arglist ,2,ierflg);
  fMn->mnexcm("IMP",arglist,2,ierflg);
  
  fMn->mnfree(0);
  return fMn->GetStatus() ;
}

void FitTower::Fcn2(Int_t & nparam, Double_t *grad, Double_t &fval, Double_t *param, Int_t iflag)
{

  // only need to translate into the old parameterization
  //
  Double_t oldParam[7];
  //// playing
  
  float dd0=we.Force2Mass*2./param[6]/(sqrt(1-param[5]*param[5]))*730.;
  float dd=param[3];
  ///
  oldParam[0] = param[0] ;
  oldParam[1] = param[1] + cos(param[4]) * dd * (1 - param[5]) / 2.0 ;
  oldParam[2] = param[2] + sin(param[4]) * dd * (1 - param[5]) / 2.0 ;
  oldParam[3] = param[6] * (1 + param[5]) / 2.0 ;
  oldParam[4] = param[1] - cos(param[4]) * dd * (1 + param[5]) / 2.0 ;
  oldParam[5] = param[2] - sin(param[4]) * dd * (1 + param[5]) / 2.0 ;
  oldParam[6] = param[6] * (1 - param[5]) / 2.0 ;
  
  // then just call "Fcn1(...)"
  //
  Fcn1(nparam, grad, fval, oldParam, iflag);
  Double_t dfval=pow(fabs(dd-dd0)/.002,2)-1;
  if(we.Force2Mass>0&& (dfval>0))
    {
      fval=fval+dfval;
      float e0=param[6];
      float zz=param[5];
      float mm=e0*sqrt(1-zz*zz)*dd/730./2.;
      if(we.YPrintLevel>1)
	{      printf("fcn with E0=%f zgg=%f m=%f  dd=%f dd0=%f val=%f \n",e0,zz,mm,dd,dd0,fval);
	  printf("x0=%f y0=%f x1=%f y1=%f \n",oldParam[1],oldParam[2],oldParam[4],oldParam[5]);
	};
    };
};


