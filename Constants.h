#ifndef Constants_
#define COnstants_
TArrayI* Vorder(TVector&);
//Float_t GevPerADC=.25;
Float_t GevPerADC=1.;
Float_t TowerWidthCM=3.81;
Float_t DistanceToInteraction=750.;
Double_t a[3]={.8,.3,-.1};
//Double_t a[3]={.8,.36,-.16};
Double_t b[3]={.8,.2,7.6};
Double_t Pi=3.1415926;
TH2F* FitHist;
TH2F* chi2Hist;
#endif
