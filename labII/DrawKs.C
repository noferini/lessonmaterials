DrawKs(Int_t step=1,Float_t ptmin=0,Float_t ptmax=10){
  TFile *fcur = TFile::Open(Form("step%i.root",step));
  TFile *fpre = TFile::Open(Form("step%i.root",step-1));

  const char *comb = "PiKa";

  TH2D *hcurT = (TH2D *) fcur->Get(Form("truePidKs%s",comb));
  TH2D *hcurP = (TH2D *) fcur->Get(Form("priorsKs%s",comb));

  TH2D *hpreT = (TH2D *) fpre->Get(Form("truePidKs%s",comb));
  TH2D *hpreP = (TH2D *) fpre->Get(Form("priorsKs%s",comb));

  TH2D *hPion = (TH2D *) fcur->Get("priorsKsPiPi");
//   hPion->Add((TH2D *) fcur->Get("priorsKsPiKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPiPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaPi"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrPi"));


  TH1D *hprior = hcurP->ProjectionX("meas",hcurP->GetYaxis()->FindBin(ptmin+0.001),hcurP->GetYaxis()->FindBin(ptmax-0.001));
  TH1D *htrue = hcurT->ProjectionX("true",hcurT->GetYaxis()->FindBin(ptmin+0.001),hcurT->GetYaxis()->FindBin(ptmax-0.001));

  TH1D *hpriorOld = hpreP->ProjectionX("measOld",hpreP->GetYaxis()->FindBin(ptmin+0.001),hpreP->GetYaxis()->FindBin(ptmax-0.001));

  TH1D *hnorm = hPion->ProjectionX("norm",hPion->GetYaxis()->FindBin(ptmin+0.001),hPion->GetYaxis()->FindBin(ptmax-0.001));


  for(Int_t i=1;i <= hprior->GetNbinsX();i++){
    Float_t err = TMath::Abs(hprior->GetBinContent(i)-hpriorOld->GetBinContent(i));

    if(hprior->GetBinContent(i)){
      err *= hnorm->GetBinContent(i);
      err /= hprior->GetBinContent(i);
    }
    hprior->SetBinError(i,err);
    


  }

  TF1 *bw = new TF1("bw","[0]/((x-[1])*(x-[1]) +([2]*[2]*0.25)) + pol1(3)",0.8,1);
  TF1 *bwsig = new TF1("bw","[0]/((x-[1])*(x-[1]) +([2]*[2]*0.25))",0.8,1);

  bw->SetLineColor(2);
  bw->SetParameter(0,100);
  bw->FixParameter(1,0.896);
  bw->FixParameter(2,0.0505);
  bw->SetParameter(3,0);
  bw->SetParameter(4,0);
  bw->SetParameter(5,0);

  TCanvas *c = new TCanvas();
  c->Divide(1,2);
  c->cd(1);
  hprior->Draw("ERR3");
  hprior->SetMinimum(0);
  htrue->SetLineColor(2);
  htrue->Draw("SAME");
  htrue->Fit(bw);

  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));
  
  Float_t truesig = bwsig->Integral(0.8,1) / htrue->GetBinWidth(1);
					       
						
//   hpriorOld->SetLineColor(4);
//   hpriorOld->Draw("SAME");

//   for(Int_t i=1;i<= htrue->GetNbinsX();i++)
//     htrue->SetBinError(i,0);

  c->cd(2);
  bw->SetParameter(0,100);
  bw->SetParameter(3,0);
  bw->SetParameter(4,0);
  bw->SetParameter(5,0);

  TH1D *hr=new TH1D(*hprior);
  hr->Add(htrue,-1);
  hr->Draw();
  hr->GetYaxis()->UnZoom();
  hr->Fit(bw);

  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));

  Float_t deltasig = bwsig->Integral(0.8,1) / hr->GetBinWidth(1);

  printf("signal = %f -- delta = %f --> Precision = %f%c\n",truesig,deltasig,deltasig/truesig*100,'%');

}
