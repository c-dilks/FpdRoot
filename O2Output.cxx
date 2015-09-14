#include "O2Output.h"
ClassImp(O2Output)

O2Output::O2Output() {
  OutFileName = "Output.root"; //the default name of the output file.
  NumberEventsToProcess = 0; //0 means process all events.
  keyIndex = kSingleKey; //for OutputMycells to output to a single key: TObject::kSingleKey = 1;
  maxsep = 0.07;
  minsep = 0.0;
  ECutMin = 60.0;
  ECutMax = 90.0;
  ZMax = 0.75;
};

O2Output::~O2Output() {

};
  
Int_t O2Output::OutputTwoTr(char* filelist, FilesSet* p_files) {
  // creates a small Output.root file (minimum) from a OFile -- similar to readq1 but with much fewer branches in the output file - by Saroj Adhikari - for learning purposes. Only one tree TwoTr (similar to the TwoTr of readq1) is created with only some properties thought of to be useful in preliminary spin analysis//
  TFile* out = new TFile(OutFileName,"recreate");  
  TTree TwoTr("TwoTr", "consists of pair events and their various properties");
  
  typedef struct{
    Int_t Bunchid7bit;
    Int_t spin;
    Int_t N12;
    Float_t M12;
    Float_t E12;
    Float_t Pt;
    Float_t Eta;
    Float_t Ntracks;
    Float_t Phiaway;
    Float_t Ptaway;
    Float_t Maway;
    Float_t Yaway;
    Float_t Y1;
    Float_t Phi;   
    Float_t Z;
  } Properties_t;
  Properties_t p2;

  TwoTr.Branch("spin", &p2.spin, "spin/I");
  TwoTr.Branch("Bunchid7bit", &p2.Bunchid7bit, "Bunchid7bit/I");
  TwoTr.Branch("N12", &p2.N12, "N12/I");
  TwoTr.Branch("M12", &p2.M12, "M12/F");
  TwoTr.Branch("E12", &p2.E12, "E12/F");
  TwoTr.Branch("Pt", &p2.Pt, "Pt/F");
//TwoTr.Branch("Pttest", &p2.Pttest, "Pt/F");
  TwoTr.Branch("Eta", &p2.Eta, "Eta/F");
  TwoTr.Branch("Ntracks", &p2.Ntracks, "Ntracks/F");
  TwoTr.Branch("Phiaway", &p2.Phiaway, "Phiaway/F");
  TwoTr.Branch("Ptaway", &p2.Ptaway, "Ptaway/F");
  TwoTr.Branch("Maway", &p2.Maway, "Maway/F");
  TwoTr.Branch("Yaway", &p2.Yaway, "Yaway/F");
  TwoTr.Branch("Y1", &p2.Y1, "Y1/F");
  TwoTr.Branch("Phi", &p2.Phi, "Phi/F");
  TwoTr.Branch("Z", &p2.Z, "Z/F"); //Z is a measure of the energy differences between two photons in case of N12==2; if one takes most of the energy then, Z is large and it is more likely that the event is not a real one as it is very easy to get a low energy photon in the background. a value of 1 is assigned for N12!=2.
  
  pout = new poutTree(filelist);
  pout->EnableEdepCorr=0; // do not enable Edep corrections

  //set energy threshold to include photons
  for(Int_t kk=1;kk<9;kk++) {
    pout->MinEnergy[kk]=6.;
    //Large Cell threshold
    if(kk==5 || kk==6)pout->MinEnergy[kk]=4.;
    };

  Int_t nentries = pout->nentries;
  if (NumberEventsToProcess>0) nentries = NumberEventsToProcess;  

  for (Int_t ievent=0; ievent<nentries; ievent++) {
      pout->GetEntry(ievent);
      p2.Bunchid7bit = pout->Bunchid7bit;
      
      pout->ClearScratch();
      pout->AllToScratch(false); // make a scratch list of all hard photons -- if the argument were true then the scratch list would also include soft photons.
      pout->ClusterwithYPhi(false);// cluster for angular size

      TObjArray* Clust=pout->ClusterScratch(maxsep);
      TLorentzVector Vscr=pout->SumScratch();
      TIter nxt(Clust);
      TObjArray* ctmp;

      if((ievent%100000)==0) printf("%d Events Read \n", ievent);
      while(ctmp=(TObjArray*) nxt())
      	{
      	  Int_t ClusterIndexOf=Clust->IndexOf(ctmp);
      	        
	        if (ctmp->GetEntries()>0)
	          {
	            p2.Z = 1.0;
	            p2.N12=ctmp->GetEntries();
	            p2.spin=pout->spin;
	            
	            TLorentzVector TotalVector = pout->SumList(ctmp);
	            if (TotalVector.Pt()==0) continue;
	            
	            p2.E12 = TotalVector.E();
	            p2.M12 = TotalVector.M();
	            //p2.Pttest = TotalVector.Pt();
	            p2.Eta = TotalVector.PseudoRapidity();
	            p2.Ntracks=pout->scratchlist->GetEntries();
	            
	            LVec* tv1=(LVec*) ctmp->First();
	            
	            //p2.Phi1=tv1->Phi();
	            //p2.Esoft=pout->ClusterSoftE(ClusterIndexOf,.9).E();
	            
	            TVector3 uv1=tv1->Vect();
	            uv1.SetMag(1.);
	            //p2.Y1=0;
	            
	            if(tv1->Pt()==0)continue;
	            
	            p2.Y1=tv1->PseudoRapidity();
	            p2.Phi=TotalVector.Phi();
	            
	            //if (ctmp->GetEntries()==2) p2.Pt = TotalVector.Pt();
	            p2.Pt = TotalVector.Pt();
	            
	            if(p2.Ntracks>p2.N12) // otherwise (i.e. if Ntracks=N12) these should be zero anyway. 
	            /* added later: but it turns out that what gets written is not zero but the last value set to the variable; so, its necessary to set the values to 0 explicitly. */
		            {
		              p2.Ptaway=(Vscr-TotalVector).Pt();
		              p2.Phiaway=(Vscr-TotalVector).Phi();
		              if(p2.Ptaway==0)continue;
		              p2.Yaway=(Vscr-TotalVector).PseudoRapidity();
		              p2.Maway=(Vscr-TotalVector).Mag();
		            }
		          else
		            {
		              // if these are not set to zero explicitly here then they get the values of what was last assigned to them.
		              p2.Ptaway = 0;
		              p2.Phiaway = 0;
		              p2.Yaway=0;
		              p2.Maway=0;
		            }
		          if (p2.N12==2) {
		            LVec* tv2 = (LVec*) ctmp->Last();
		            p2.Z = fabs(tv1->E() - tv2->E())/p2.E12;
		            }
		            
		          TwoTr.Fill(); //Fill the tree
		        }
		    }   
    };
  TwoTr.Write();	/* write to the file OutFileName that is  open */
  return 0;
};

Int_t O2Output::OutputMycells(char* filelist, FilesSet* p_files, Bool_t doSingle) {
  // creates a small output file (minimum) from a OFile -- similar to readqM but with much fewer branches in the output file - by Saroj Adhikari - for learning purposes; the output file contains Cell instances for each of the cells which can then be used to create _Cellr%_c%_%.root files needed for calibration analysis; this code will only create an output file with minimum information necessary for calibration and will be expanded later as needed //
  TFile* out = new TFile(OutFileName, "recreate");
  Qt qt(p_files);
  Geom* p_geom=new Geom(p_files);
  pout = new poutTree(filelist);
  pout->GetEntry(0);
  Int_t RunNum = pout->Rnum;
  p_files->Print(); // necessary to initialize ->path below or else use Path()
  CalibStr gcor(RunNum, (const char*) p_files->p_fpdgaincorr()->path);
  CalibStr gain(RunNum, (const char*) p_files->p_fpdgain()->path);
  
  Cell* cell[34][17][4];
  TObjArray mcells(2400,0);  
  Bool_t ProcessAll = false;
  
  for (Int_t i=0; i<4; i++) {
    Int_t up=34;
    if (i>1) up=24;
    for (Int_t j=0; j<34; j++) {
      for (Int_t k=0; k<17; k++) {
        cell[j][k][i] = 0;
        if (j>=up || k>=up/2) continue;
        if ((i<2) && (k<8) && (j>8) && (j<25)) continue;
        if ((i>1) && (k<5) && (j>6) && (j<17)) continue;
        
        char nam[40];
	      sprintf(nam,"Cellr%d_c%d_%d",j,k,i);
	      std::cout<<"about to make: "<<nam<<"\n";        
        out->cd();
        Cell* pcell = new Cell(2, i+1, j+1, k+1, RunNum,p_geom, &gain, &gcor);
        pcell->SetName(nam);
        cell[j][k][i] = pcell;
        mcells.Add(pcell);
      
        sprintf(nam,"CellD%d_c%d_%d",j,k,i);
	      pcell->InitTree(nam);
      };
    };
  };
  
  TH2F* qtHist[4];
  TH2F* qtHistLed[4];
  TFile adcroot("adc.root");
  if(adcroot.IsOpen()) {
    qtHist[0]=(TH2F*) adcroot.GetKey("qtHistLN")->ReadObj();
    qtHist[1]=(TH2F*) adcroot.GetKey("qtHistLS")->ReadObj();
    qtHist[2]=(TH2F*) adcroot.GetKey("qtHistSN")->ReadObj();
    qtHist[3]=(TH2F*) adcroot.GetKey("qtHistSS")->ReadObj();

    qtHistLed[0]=(TH2F*) adcroot.GetKey("qtHistLedLN")->ReadObj();
    qtHistLed[1]=(TH2F*) adcroot.GetKey("qtHistLedLS")->ReadObj();
    qtHistLed[2]=(TH2F*) adcroot.GetKey("qtHistLedSN")->ReadObj();
    qtHistLed[3]=(TH2F*) adcroot.GetKey("qtHistLedSS")->ReadObj();
  }

  remove ("data.txt");
  FILE* fp=fopen("data.txt","a");  
  ProcessAll=true;
  out->cd();

  for(Int_t kk=1;kk<9;kk++)pout->MinEnergy[kk]=4;
  pout->MinEnergy[5]=2.;
  pout->MinEnergy[6]=2.;
  
  Int_t nentries = pout->nentries;
  if(NumberEventsToProcess>0)nentries=NumberEventsToProcess;  
  if(NumberEventsToProcess>nentries)nentries=pout->nentries;
  p_geom->FMSGeom=true;
  Float_t EnergyNSTB[4];
  printf("Number of Events to Process: %d \n", nentries);
  for (Int_t ievent=0; ievent<nentries; ievent++) {
    if((ievent%100000)==0) printf("%d Events Read \n", ievent);    
    pout->GetEntry(ievent);
    RunNum = pout->Rnum;
    Float_t mass;
    Float_t energy;
    Float_t Z;
    //printf("here0 \n");
    pout->ClearScratch();
    pout->AllToScratch(false);
    pout->ClusterwithYPhi(false);
    TObjArray* Clust = pout->ClusterScratch(maxsep); //maxseparation = .045 for pi0 reconstruction.
    //TLorentzVector Vscr = pout->SumScratch(); //not used
    TIter nxt(Clust);
    TObjArray* ctmp;
    //printf("here-1 \n");
    while(ctmp=(TObjArray*) nxt()) {
  	  Int_t ClusterIndexOf = Clust->IndexOf(ctmp);
  	  /* the claibration works by looking at pi0 masses, which is detected as two photons close together on the detector (since pi0 decays into two photons separated by a small angle from its path of motion); so, we just look for two photons. The small angular separation is characterized by the .045 anglular separation above in pout->ClusterScratch */
      //printf("here \n");
  	  if (!doSingle && ctmp->GetEntries()==2) {
	      LVec* photon[2];
	      photon[0] = (LVec*) ctmp->First(); //first photon
	      photon[1] = (LVec*) ctmp->After(photon[0]); //other photon
	      TLorentzVector twophotons = *(photon[0]) + *(photon[1]);
	      mass = twophotons.Mag();
	      energy = twophotons.E();
        Float_t y2 = (*photon[0] + *photon[1]).PseudoRapidity();

// eta selection based on energy and pseudorapidity
        if (energy<ECutMin || energy>ECutMax) continue;
        if (fabs(y2-3.7)>0.2) continue;
        /*
        if (fabs(energy-85.0)<5.0) {
          if (fabs(y2-3.75)>0.15) continue;
        }
        if (fabs(energy-95.0)<5.0) {
          if (fabs(y2-3.75)>0.25) continue;
        }
        if (fabs(energy-110.0)<10.0) {
          if (fabs(y2-3.85)>0.35) continue;
        }
        if (fabs(energy-135.0)<15.0) {
          if (fabs(y2-3.9)>0.3) continue;
        }
        */
                
	      if (photon[0]->Nstb==0) continue;
	      if (photon[1]->Nstb==0) continue;
	      if (photon[0]->Nstb != photon[1]->Nstb) continue; // remove events crossing boundary
	      TVector3 phpos0 = pout->PosInDet(photon[0], p_geom, 0.);
	      TVector3 phpos1 = pout->PosInDet(photon[1], p_geom, 0.);
	      TVector3 xyab[2];
	      xyab[0] = p_geom->LocalXYZ(2, photon[0]->Nstb, phpos0, true);
	      xyab[1] = p_geom->LocalXYZ(2, photon[1]->Nstb, phpos1, true);
	      Int_t phRow[2], phCol[2];
	      phRow[0] = (Int_t) xyab[0].Y();
	      phRow[1] = (Int_t) xyab[1].Y();
	      phCol[0] = (Int_t) xyab[0].X();
	      phCol[1] = (Int_t) xyab[1].X();
	      
	      for (Int_t iz=0; iz<2; iz++) {
	        /* check if the Row and Column numbers are within the detector range */
	        Bool_t good = true;

	        Int_t iz2 = 1 - iz;
	        if (phRow[iz]<0 || phCol[iz]<0 || phRow[iz2]<0 || phCol[iz2]<0) good = false;
	        if (phRow[iz]>=qt.ROW_NUM[photon[iz]->Nstb-1] || phCol[iz]>=qt.COL_NUM[photon[iz]->Nstb-1]) good = false;
	        if (phRow[iz2]>=qt.ROW_NUM[photon[iz2]->Nstb-1] || phCol[iz2]>=qt.COL_NUM[photon[iz2]->Nstb-1]) good = false;
	        
	        if (photon[iz]->Iew==2) {
	          if (photon[iz]->Nstb<3) {
	            if (phCol[iz]<8 && phRow[iz]>8 && phRow[iz]<25) good = false;
	          }
	          else {
	            if (phCol[iz]<5 && phRow[iz]>6 && phRow[iz]<17) good = false;
	          }
	        }   

	        if (good) {

	          Bool_t FillIt = false; 
	          Cell* ocl = cell[phRow[iz2]][phCol[iz2]][photon[iz2]->Nstb-1];
	          Cell* tcl = cell[phRow[iz]][phCol[iz]][photon[iz]->Nstb-1];

	          if (ocl!=0 && tcl!=0) {

	            if ((ocl->Masspeakfraction>.15 && ocl->Massmaxcontents>10) || ProcessAll) {
	              FillIt = true;
	            }
	          }

	          if (tcl!=0 && ProcessAll) FillIt = true;
            
            Z = (fabs(photon[0]->E()-photon[1]->E()))/(photon[0]->E()+photon[1]->E());
		        if (Z>ZMax) FillIt = false;

	          if (FillIt) {  
	            //printf("fillit is true \n");          
	            
	            tpCellDat cd;
	            cd.Mcell = mass;
	            //printf("mass: %f \n", mass);
	            cd.Epair = photon[iz]->E() + photon[iz2]->E();
	            cd.Ypair = y2;
	            cd.E1 = photon[iz]->E();
	            cd.E2 = photon[iz2]->E();
	            cd.X1 = xyab[iz].X();
	            cd.X2 = xyab[iz2].X();
	            cd.Y1 = xyab[iz].Y();
	            cd.Y2 = xyab[iz2].Y();
	            cd.NSTB = photon[iz]->Nstb; // detector number for the first photon
	            cd.ievt = pout->ievt;
	            cd.Rnum = pout->Rnum;
	            cd.nSavedHits = pout->nSavedHits;
	            
	            if (cd.nSavedHits > 100) cd.nSavedHits = 100;
	            for (int in = 0; in < cd.nSavedHits; in++) {
	              cd.SavedHits[in] = pout->SavedHits[in];
	            };

	            tcl->FillTree(cd);
              char setvalchar[100];
	            sprintf(setvalchar,"%d,%d,%d\n",phRow[iz],phCol[iz],cd.NSTB);
              fputs(setvalchar,fp);
              
	            //printf("here filling tree %d %d \n", phRow[iz], phCol[iz]);
	          };
	        };
	      };      
  	  }
  	  else {
  	    if (doSingle && ctmp->GetEntries()==1) {
  	    // create Output.root for single photon events Ntracks==1 (and then obviously N12==1;
  	      LVec* photon; 
	        photon = (LVec*) ctmp->First(); //first and only photon 	    
  	      mass = photon->Mag();
  	      energy = photon->E();
  	      if (photon->Nstb==0) continue;
  	      TVector3 phpos = pout->PosInDet(photon, p_geom, 0.);
  	      TVector3 xyab;
  	      xyab = p_geom->LocalXYZ(2, photon->Nstb, phpos, true);
  	      Int_t phRow, phCol;
  	      phRow = (Int_t) xyab.Y();
  	      phCol = (Int_t) xyab.X();
  	      
          TLorentzVector TotalVector = pout->SumList(ctmp);
	        if (TotalVector.Pt()==0) continue;  	      
  	
  	      Bool_t good = true;
  	      //check if the row and cols and within the detecor range
  	      if (phRow<0 || phCol<0 || phRow>=qt.ROW_NUM[photon->Nstb-1] || phCol >= qt.COL_NUM[photon->Nstb-1] ) good = false;
         
          if (photon->Iew==2) {
	          if (photon->Nstb<3) {
	            if (phCol<8 && phRow>8 && phRow<25) good = false;
	          }
	          else {
	            if (phCol<5 && phRow>6 && phRow<17) good = false;
	          }
	        }
	        
          if (good) {
	          Bool_t FillIt = false; 
	          Cell* ocl = cell[phRow][phCol][photon->Nstb-1];
	          Cell* tcl = cell[phRow][phCol][photon->Nstb-1];

	          if (ocl!=0 && tcl!=0) {

	            if ((ocl->Masspeakfraction>.15 && ocl->Massmaxcontents>10) || ProcessAll) {
	              FillIt = true;
	            }
	          }

	          if (tcl!=0 && ProcessAll) FillIt = true;

	          if (FillIt) {  
	            //printf("fillit is true \n");          
	            
	            tpCellDat cd;
	            cd.Mcell = mass;
	            //printf("mass: %f \n", mass);
	            cd.Epair = energy;
	            cd.Ypair = photon->PseudoRapidity();
	            cd.E1 = photon->E();
	            cd.E2 = pout->scratchlist->GetEntries(); //the Ntracks variable set to E2 which is of no use here.
	            cd.X1 = xyab.X();
	            cd.X2 = TotalVector.Pt();
	            cd.Y1 = xyab.Y();
	            cd.Y2 = pout->spin;
	            cd.NSTB = photon->Nstb; // detector number for the first photon
	            cd.ievt = pout->ievt;
	            cd.Rnum = pout->Rnum;
	            cd.nSavedHits = pout->nSavedHits;
	            
	            if (cd.nSavedHits > 100) cd.nSavedHits = 100;
	            for (int in = 0; in < cd.nSavedHits; in++) {
	              cd.SavedHits[in] = pout->SavedHits[in];
	            };

	            tcl->FillTree(cd);
              char setvalchar[100];
	            sprintf(setvalchar,"%d,%d,%d\n",phRow,phCol,cd.NSTB);
              fputs(setvalchar,fp);
              
	            //printf("here filling tree %d %d \n", phRow[iz], phCol[iz]);
	          };
	        };	          	      
  	    }
  	  }
  	};    
  };
  fclose(fp);
  // out of nentries loop
  printf("out of nentries loop");
  
  TIter next2(&mcells);
  Cell* pcl;
  while(pcl=(Cell*) next2()) {
    Int_t idet=1;
    if(pcl->Iew==2&& ((idet=pcl->Instb)<5)) {
      Int_t chan=(pcl->Row1-1)*qt.COL_NUM[idet-1] + pcl->Col1;
  	  if(pcl->p_adc)delete pcl->p_adc;
  	  if(pcl->p_adcLed)delete pcl->p_adcLed;
  	  char nam[40];
  	  sprintf(nam,"adc_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
  	  printf("set adc %s \n",nam);
  	  pcl->p_adc=qtHist[idet-1]->ProjectionY(nam,chan,chan);
  	  sprintf(nam,"adcLed_r%d_c%d_%d",pcl->Row1-1,pcl->Col1-1,pcl->Instb-1);
  	  printf("set adcLed %s \n",nam);
  	  pcl->p_adcLed=qtHistLed[idet-1]->ProjectionY(nam,chan,chan);
    };
  };
  std::cout<<"Finish adc setup\n";
  out->cd();
  for (Int_t j=0; j<4; j++) {
    qtHist[j]->Write();
    qtHist[j]->GetNbinsX(), qtHist[j]->GetNbinsY();
  }

  mcells.Print();
  printf("%d", keyIndex);
  mcells.Write("Mycells", keyIndex);
};


