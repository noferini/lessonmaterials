void DrawLc3D(Int_t step=1,Float_t ptmin=0,Float_t ptmax=10){
  TFile *fcur = TFile::Open(Form("step%i.root",step));
  TFile *fpre = TFile::Open(Form("step%i.root",step-1));

  if(! fpre) fpre = fcur;

  const char *comb = "PiKaPr";

  Int_t rebin=1;

  TH3D *hcurT = (TH3D *) fcur->Get(Form("truePidLc%s",comb));
  TH3D *hcurP = (TH3D *) fcur->Get(Form("priorsLc%s",comb));
  TH3D *htrueLc = (TH3D *) fcur->Get(Form("trueLc"));
  TH3D *hmypidLc = (TH3D *) fcur->Get(Form("mypidLc"));

  TH3D *hpreT = (TH3D *) fpre->Get(Form("truePidLc%s",comb));
  TH3D *hpreP = (TH3D *) fpre->Get(Form("priorsLc%s",comb));

  TH3D *hPion = (TH3D *) fcur->Get("priorsLcPiPiPi");
//   hPion->Add((TH2D *) fcur->Get("priorsKsPiKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPiPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsKaPi"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrKa"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrPr"));
//   hPion->Add((TH2D *) fcur->Get("priorsKsPrPi"));

  hcurT->RebinX(rebin);
  hcurP->RebinX(rebin);
  hpreT->RebinX(rebin);
  hpreP->RebinX(rebin);
  hPion->RebinX(rebin);
  htrueLc->RebinX(rebin);
  if(hmypidLc) hmypidLc->RebinX(rebin);

  TH1D *hprior = hcurP->ProjectionX("meas",hcurP->GetYaxis()->FindBin(ptmin+0.001),hcurP->GetYaxis()->FindBin(ptmax-0.001),0,-1);

  TH1D *htrue = hcurT->ProjectionX("true",hcurT->GetYaxis()->FindBin(ptmin+0.001),hcurT->GetYaxis()->FindBin(ptmax-0.001),0,-1);

  TH1D *hreal = htrueLc->ProjectionX("real",hcurT->GetYaxis()->FindBin(ptmin+0.001),hcurT->GetYaxis()->FindBin(ptmax-0.001),0,-1);

  TH1D *hmypid = NULL;
  if(hmypidLc) hmypid = hmypidLc->ProjectionX("mypid",hcurT->GetYaxis()->FindBin(ptmin+0.001),hcurT->GetYaxis()->FindBin(ptmax-0.001),0,-1);

  TH1D *hpriorOld = hpreP->ProjectionX("measOld",hpreP->GetYaxis()->FindBin(ptmin+0.001),hpreP->GetYaxis()->FindBin(ptmax-0.001),0,-1);

  TH1D *hnorm = hPion->ProjectionX("norm",hPion->GetYaxis()->FindBin(ptmin+0.001),hPion->GetYaxis()->FindBin(ptmax-0.001),0,-1);



  for(Int_t i=1;i <= hprior->GetNbinsX();i++){
    Float_t err = TMath::Abs(hprior->GetBinContent(i)-hpriorOld->GetBinContent(i));

    if(hnorm->GetBinContent(i) > hprior->GetBinContent(i)){
      err *= hnorm->GetBinContent(i);
      err /= hprior->GetBinContent(i);
    }

    err = sqrt(err*err + hprior->GetBinContent(i));

    if(err < 1) err = 1;

    if(htrue->GetBinError(i)<1)htrue->SetBinError(i,1);

    hprior->SetBinError(i,err);
    


  }

  TF1 *bw = new TF1("bw","gaus(0) + pol1(3)",2.23,2.35);
  TF1 *bwsig = new TF1("bw","gaus(0)",0.8,1);

  bw->SetLineColor(2);
  bw->SetParameter(0,100);
  bw->FixParameter(1,2.28646);
  bw->FixParameter(2,0.008);
  bw->SetParameter(3,0);
  bw->SetParameter(4,0);
  bw->SetParameter(5,0);

  TCanvas *c = new TCanvas();
  c->Divide(1,2);
  c->cd(1);
  hprior->Draw("ERR3");
  hprior->SetMinimum(0);
  htrue->Draw("SAME");
  htrue->SetLineColor(2);
  htrue->Fit(bw,"EI","",2.25,2.32);

  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));
  
  Float_t truesig = bwsig->Integral(2.28646-3*0.008,2.28646+3*0.008) / htrue->GetBinWidth(1);
  Float_t backgrd = bw->Integral(2.28646-3*0.008,2.28646+3*0.008) / htrue->GetBinWidth(1);

  hprior->Draw("SAME");
  hprior->Fit(bw,"EI","",2.25,2.32);

  Float_t fiterror = 0;
  if(bw->GetParameter(0)) fiterror = TMath::Abs(bw->GetParError(0)/bw->GetParameter(0));

  hreal->Draw("SAME");

  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));
  
  Float_t meassig = bwsig->Integral(2.28646-3*0.008,2.28646+3*0.008) / htrue->GetBinWidth(1);
					       
						
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
  hr->GetYaxis()->SetTitle("(meas - true)/true");
  hr->Add(htrue,-1);
  hr->Divide(htrue);
  hr->Draw();
  hr->GetYaxis()->UnZoom();
  //  hr->Fit(bw,"","",0.8,1.);

  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));

  Float_t deltasig = bwsig->Integral(2.28646-3*0.008,2.28646+3*0.008) / hr->GetBinWidth(1);

  printf("signal = %f(true) %f(meas) -- delta  = %f (%f) --> Relative Error = %f(stat=%f, fit=%f)%c\n",truesig,meassig,meassig-truesig,deltasig,(meassig-truesig)/truesig*100,100./sqrt(TMath::Abs(truesig)),fiterror*100,'%');

  printf("significance = %f\n",hreal->Integral(1,hreal->GetNbinsX())/sqrt(backgrd));
  printf("significance measured = %f\n",meassig/sqrt(backgrd));

  printf("reco Lambdac = %i\n",hreal->Integral(1,hreal->GetNbinsX()));

  if(!hmypidLc) return;

  new TCanvas();
  hmypid->Draw();
  hmypid->Fit(bw,"EI","",2.25,2.32);
  for(Int_t i=0;i < 3;i++) bwsig->SetParameter(i,bw->GetParameter(i));

  Float_t mypidsig = bwsig->Integral(2.28646-3*0.008,2.28646+3*0.008) / hmypid->GetBinWidth(1);
  printf("my pid signal = %f\n",mypidsig);
  printf("my pid significance = %f\n",mypidsig/sqrt(backgrd));
}
