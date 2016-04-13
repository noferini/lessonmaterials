TF1 *fresp;
TF1 *fseparation;

Float_t getCoef(Int_t iref,Int_t icont,Float_t pr[3],Float_t pt);

void priorsCorr(Int_t step){
  TMatrix  A_C(3,3);
  TMatrix B(3,1);
  TMatrix Delta(3,1);

  fresp = new TF1("fresp","gaus(0)*gaus(3)*3.98942280401432703e-01/(gaus(6)+gaus(9)+gaus(12))",-20,20);

  TFile *f = new TFile(Form("step%i.root",step));
  TFile *f2 = new TFile(Form("step%i.root",step-1));
  
  TH1D *hpr[3];
  TH1D *hmeas[3];
  TH1D *htr[3];
  TH1D *hdelta[3];

  hpr[0] = (TH1D *) f2->Get("priorsPtPiP");
  hpr[1] = (TH1D *) f2->Get("priorsPtKaP");
  hpr[2] = (TH1D *) f2->Get("priorsPtPrP");

  hmeas[0] = (TH1D *) f->Get("priorsPtPiP");
  hmeas[1] = (TH1D *) f->Get("priorsPtKaP");
  hmeas[2] = (TH1D *) f->Get("priorsPtPrP");
  
  htr[0] = (TH1D *) f->Get("truePtPiP");
  htr[1] = (TH1D *) f->Get("truePtKaP");
  htr[2] = (TH1D *) f->Get("truePtPrP");

  hdelta[0] = new TH1D(*htr[0]);
  hdelta[0]->SetName("deltaPtPiP");
  hdelta[1] = new TH1D(*htr[2]);
  hdelta[1]->SetName("deltaPtKaP");
  hdelta[2] = new TH1D(*htr[1]);
  hdelta[2]->SetName("deltaPtPrP");

  fseparation = new TF1("fsep","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.);
  fseparation->SetParameter(1,7.);

  Int_t nbin = hpr[0]->GetNbinsX();

  Float_t pt;
  Float_t alfa[3][3],pr[3],meas[3];


  for(Int_t i=1;i<=nbin;i++){
    pt = hpr[0]->GetBinCenter(i);

    for(Int_t ii=0;ii<3;ii++){
      pr[ii] = hpr[ii]->GetBinContent(i);
      meas[ii] = hmeas[ii]->GetBinContent(i);
    }

    // prepare matrix
    for(Int_t ii=0;ii<3;ii++)
      for(Int_t kk=0;kk<3;kk++){
	alfa[ii][kk] = getCoef(ii,kk,pr,pt);
	A_C(ii,kk) = -alfa[ii][kk];
	if(ii==kk){
	  for(Int_t jj=0;jj<3;jj++)
	    A_C(ii,kk) = A_C(ii,kk) + pr[jj]/pr[ii]*alfa[ii][jj];
	}
      }
    // end matrix preparation
    if(A_C.Determinant() != 0) A_C.Invert();
    else{
      printf("Determinant equal to zero\n");
      for(Int_t ii=0;ii<3;ii++)
	for(Int_t jj=0;jj<3;jj++)
	  A_C(ii,jj)=0;
    }

 
    // prepare Vector
    for(Int_t ii=0;ii<3;ii++){
      B(ii,0) = 0;
      for(Int_t jj=0;jj<3;jj++)
	B(ii,0) = B(ii,0)  + alfa[ii][jj]*(meas[jj]/meas[ii] - pr[jj]/pr[ii]);
      B(ii,0) = B(ii,0) * meas[ii];
    }
    // end vector preparation


    // compute errors
    A_C.Print();
    B.Print();
    Delta = A_C * B;

    Delta.Print();

    for(Int_t ii=0;ii<3;ii++)
      hdelta[ii]->SetBinContent(i,Delta(ii,0));
  }
  hmeas[0]->Add(htr[0],-1);
  hmeas[1]->Add(htr[1],-1);
  hmeas[2]->Add(htr[2],-1);
  
  hmeas[0]->SetLineColor(2);
  hmeas[1]->SetLineColor(2);
  hmeas[2]->SetLineColor(2);
  
  hdelta[0]->Divide(htr[0]);
  hmeas[0]->Divide(htr[0]);
  hdelta[1]->Divide(htr[0]);
  hmeas[1]->Divide(htr[0]);
  hdelta[2]->Divide(htr[0]);
  hmeas[2]->Divide(htr[0]);

  
  for(Int_t i=1;i<=nbin;i++){
    hmeas[0]->SetBinError(i,0);
    hmeas[1]->SetBinError(i,0);
    hmeas[2]->SetBinError(i,0);
  }

  hmeas[1]->Draw();
  hdelta[1]->Draw("SAME");
}



Float_t getCoef(Int_t iref,Int_t icont,Float_t pr[3],Float_t pt){
  Float_t nsigma[3];
  nsigma[0] = -fseparation->Eval(pt);
  nsigma[1] = 0;
  nsigma[2] = fseparation->Eval(pt);

  // set coef
  fresp->SetParameter(0,1);
  fresp->SetParameter(3,1);
  fresp->SetParameter(6,pr[0]/pr[iref]);
  fresp->SetParameter(9,pr[1]/pr[iref]);
  fresp->SetParameter(12,pr[2]/pr[iref]);

  // set mean
  fresp->SetParameter(1,nsigma[icont]-nsigma[icont]);
  fresp->SetParameter(4,nsigma[iref]-nsigma[icont]);
  fresp->SetParameter(7,nsigma[0]-nsigma[icont]);
  fresp->SetParameter(10,nsigma[1]-nsigma[icont]);
  fresp->SetParameter(13,nsigma[2]-nsigma[icont]);

  // set sigma to 1
  fresp->SetParameter(2,1);
  fresp->SetParameter(5,1);
  fresp->SetParameter(8,1);
  fresp->SetParameter(11,1);
  fresp->SetParameter(14,1);

  return fresp->Integral(-20,20);
}

