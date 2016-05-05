#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

#include "particle.h"

TF1 *fseparation;
TF1 *fseparationPiKa;
TF1 *fseparationKaPr;

TF2 *fTOFpi;
TF2 *fTOFka;
TF2 *fTOFpr;
TF1 *fTPCpi;
TF1 *fTPCka;
TF1 *fTPCpr;

// bethe block parameter
Float_t fKp1=0.0283086;
Float_t fKp2=2.63394e+01;
Float_t fKp3=5.04114e-11;
Float_t fKp4=2.12543;
Float_t fKp5=4.88663;
Double_t BetheBlochAleph(Double_t *x,Double_t *par);

TTree *t;
Float_t pt,pz,phi,eta;
Int_t id,reco;
Float_t sig;
Float_t sigM;
Float_t sigTOF;
Float_t sigTPC;
Int_t mother;
Float_t vxy;
Int_t tof;
Int_t ev=0;

void FillTree(particle part);
void FillKine(particle &part,TH1D *h);

const char *filein="4050";

Bool_t kALICEseparation=kTRUE;


TF1 *fEfficiencyPi;
TF1 *fEfficiencyKa;
TF1 *fEfficiencyPr;
TF1 *fEfficiencyPiTOF;
TF1 *fEfficiencyKaTOF;
TF1 *fEfficiencyPrTOF;

int main(){

  fEfficiencyPi = new TF1("fEfficiencyPi","(x > 0.2)*(x-0.2)*(x<0.5)*3 + (x>0.5)*0.9",0,10);
  fEfficiencyKa = new TF1("fEfficiencyPi","(x > 0.2)*(x-0.2)*(x<0.5)*3 + (x>0.5)*0.9",0,10);
  fEfficiencyPr = new TF1("fEfficiencyPi","(x > 0.2)*(x-0.2)*(x<0.5)*3 + (x>0.5)*0.9",0,10);
  fEfficiencyPiTOF = new TF1("fEfficiencyPiTOF","(x > 0.3)*0.7",0,10);
  fEfficiencyKaTOF = new TF1("fEfficiencyKaTOF","(x > 0.3)*(x-0.3)*(x<1) + (x>1)*0.7",0,10);
  fEfficiencyPrTOF = new TF1("fEfficiencyPrTOF","(x > 0.3)*(x-0.3)*(x<1) + (x>1)*0.7",0,10);


  // x=p, y=pt/p (normalized at the number of sigma assuming 80 ps resolution)
  fTOFpi = new TF2("fTOFpi","3.7/y*(sqrt(x*x+0.0193210)/x-1)*37.47405725",0.3,10,0.5,1);
  fTOFka = new TF2("fTOFka","3.7/y*(sqrt(x*x+0.243049)/x-1)*37.47405725",0.3,10,0.5,1);
  fTOFpr = new TF2("fTOFpr","3.7/y*(sqrt(x*x+0.879844)/x-1)*37.47405725",0.3,10,0.5,1);


  // x=p, already normalized in number of sigma (sigma assumed 3.5=7% of the MIP)
  fTPCpi = new TF1("fTPCpi",BetheBlochAleph,0,10,6);
  fTPCpi->SetParameter(0,fKp1);
  fTPCpi->SetParameter(1,fKp2);
  fTPCpi->SetParameter(2,fKp3);
  fTPCpi->SetParameter(3,fKp4);
  fTPCpi->SetParameter(4,fKp5);
  fTPCpi->SetParameter(5,0.139);
  fTPCka = new TF1("fTPCka",BetheBlochAleph,0,10,6);
  fTPCka->SetParameter(0,fKp1);
  fTPCka->SetParameter(1,fKp2);
  fTPCka->SetParameter(2,fKp3);
  fTPCka->SetParameter(3,fKp4);
  fTPCka->SetParameter(4,fKp5);
  fTPCka->SetParameter(5,0.493);
  fTPCpr = new TF1("fTPCpr",BetheBlochAleph,0,10,6);
  fTPCpr->SetParameter(0,fKp1);
  fTPCpr->SetParameter(1,fKp2);
  fTPCpr->SetParameter(2,fKp3);
  fTPCpr->SetParameter(3,fKp4);
  fTPCpr->SetParameter(4,fKp5);
  fTPCpr->SetParameter(5,0.938);

  // simulation parameters (tuned on PbPb 10-20% @ 2.76 ATeV)
  Int_t nevents = 1000;

  TH1D *hspectra[7];
  Float_t npartPerY[7] = {0,0,0,0,0,0,0};

  TFile *f = new TFile(Form("spectra%s.root",filein));
  hspectra[0] = (TH1D *) f->Get(Form("hpi%s",filein));
  hspectra[1] = (TH1D *) f->Get(Form("hka%s",filein));
  hspectra[2] = (TH1D *) f->Get(Form("hpr%s",filein));
  hspectra[3] = (TH1D *) f->Get(Form("hks%s",filein));
  hspectra[4] = (TH1D *) f->Get(Form("hph%s",filein));
  hspectra[5] = (TH1D *) f->Get(Form("hde%s",filein));
  hspectra[6] = (TH1D *) f->Get(Form("hlc%s",filein));

  // integrated yields per unit of rapidity
  for(Int_t isp=0;isp<7;isp++){
    for(Int_t ib=1;ib <= hspectra[isp]->GetNbinsX();ib++){
      npartPerY[isp] += hspectra[isp]->GetBinContent(ib)* hspectra[isp]->GetBinWidth(ib);
    }
  }


  Int_t npion = npartPerY[0]*2*2;
  Int_t nkaon = npartPerY[1]*2*2;
  Int_t nproton = npartPerY[2]*2*2;
  Int_t nk0star = npartPerY[3]*2*2;
  Int_t nphi = npartPerY[4]*2;
  Int_t ndelta = npartPerY[5]*2*2; // not measured
  Float_t nlambdac = npartPerY[6]*2*2; // not measured
  
  // branching ratio for lambda_c
  Float_t br1 = 0.028; // Lambda_c -> pi K p no resonant
  Float_t br2 = 0.016; // Lambda_c -> K0* p -> pi K p
  Float_t br3 = 0.0086; // Lambda_c -> K Delta -> pi K p 

  gRandom->SetSeed(0);

  fseparation = new TF1("f","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.);
  fseparation->SetParameter(1,7.);

  fseparationPiKa = new TF1("fPiKa","[0]+[1]/TMath::Power(x,2.5)",0,100);
  fseparationPiKa->SetParameter(0,2.34);
  fseparationPiKa->SetParameter(1,1);

  fseparationKaPr = new TF1("fKaPr","[0]+[1]/TMath::Power(x,2.5)",0,100);
  fseparationKaPr->SetParameter(0,1);
  fseparationKaPr->SetParameter(1,56);

  TFile *fout = new TFile("out.root","RECREATE");
  t = new TTree("tree","tree");
  t->Branch("ev",&ev,"ev/I"); // number of event
  t->Branch("id",&id,"id/I"); // id particle = integer corresponding to the position in the particle type array
  t->Branch("pt",&pt,"pt/F"); // transverse momentum
  t->Branch("pz",&pz,"pz/F"); // longitudinal momentum
  t->Branch("eta",&eta,"eta/F"); // pseudorapidity
  t->Branch("phi",&phi,"phi/F"); // azhimuthal angle
  t->Branch("signal",&sigM,"signal/F"); // PID signal defined using fseparation
  t->Branch("signalTPC",&sigTPC,"signalTPC/F"); // PID signal defined using fseparation
  t->Branch("signalTOF",&sigTOF,"signalTOF/F"); // PID signal defined using fseparation
  t->Branch("mother",&mother,"mother/I"); // -1=primary, otherwise id particle of the mother
  t->Branch("reco",&reco,"reco/I"); // -1=primary, otherwise id particle of the mother
  // t->Branch("vxy",&vxy,"vxy/F");
  // t->Branch("tof",&tof,"tof/I");


  // define particle types (particle type array)
  particle::AddParticleType("pi+",0.139,1); // 0
  particle::AddParticleType("pi-",0.139,-1); // 1
  particle::AddParticleType("K+",0.493,1); // 2 
  particle::AddParticleType("K-",0.493,-1); // 3
  particle::AddParticleType("p+",0.938,1); // 4 
  particle::AddParticleType("p-",0.938,-1); // 5
  particle::AddParticleType("K0*",0.896,0,5.05e-02); // 6
  particle::AddParticleType("K0bar*",0.896,0,5.05e-02); // 7
  particle::AddParticleType("Phi",1.02,0,0.00426); // 8
  particle::AddParticleType("Delta++",1.232,2,0.118); // 9 
  particle::AddParticleType("Delta--",1.232,-2,0.118);  // 10
  particle::AddParticleType("Lambdac+",2.28646,1,0.08); // 11
  particle::AddParticleType("Lambdacbar-",2.28646,-1,0.08); // 12


  // define same dummy particles useful to perform decays
  particle res1("K0*");
  particle res2("K0bar*");
  particle res3("Delta++");
  particle res4("Delta--");
  particle d1("pi+");
  particle d2("K+");
  particle d3("p+");
  particle d4("pi-");
  particle d5("K-");
  particle d6("p-");

  particle part;

  for(ev=0;ev < nevents;ev++){ // event loop
    //    t->Reset();

    mother = -1;
    vxy = 0;
    tof = 1;

    // pions
    for(Int_t j=0;j < npion;j++){
      part.ChangeParticleType(gRandom->Rndm() > 0.5);
      FillKine(part,hspectra[0]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
      else sig = -fseparation->Eval(pt);
      FillTree(part);
    }

    // kaons
    for(Int_t j=0;j < nkaon;j++){
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+2);
      FillKine(part,hspectra[1]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = 0;
      FillTree(part);
    }

    // protons
    for(Int_t j=0;j < nproton;j++){
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+4);
      FillKine(part,hspectra[2]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
      else  sig = fseparation->Eval(pt);
      FillTree(part);
    }

    // K0*
    for(Int_t j=0;j < nk0star;j++){
      mother = -1;

      part.ChangeParticleType((gRandom->Rndm() > 0.5)+6);
      FillKine(part,hspectra[3]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);
      
      if(part.GetParticleType() == 7){
	if(gRandom->Rndm() < 0.5){
	  mother = 7;
	  part.Decay2body(d1,d5);
	  pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	  sig = 0;
	  FillTree(d5);
	}
      }
      else{
	if(gRandom->Rndm() < 0.5){
	  mother = 6;
	  part.Decay2body(d4,d2);
	  pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	  sig = 0;
	  FillTree(d2);
	}
      }

    }

    // Phi
    for(Int_t j=0;j < nphi;j++){
      mother = -1;
      part.ChangeParticleType(8);
      FillKine(part,hspectra[4]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);
      if(gRandom->Rndm() < 0.492){
	mother = 8;
	part.Decay2body(d2,d5);
	pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	sig = 0;
	FillTree(d2);
	pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	sig = 0;
	FillTree(d5);
      }
    }

    // Delta
    for(Int_t j=0;j < ndelta;j++){
      mother = -1;
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+9);
      FillKine(part,hspectra[5]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);

      if(part.GetParticleType() == 9){
	mother = 9;
	part.Decay2body(d1,d3);
	pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	else 	sig = -fseparation->Eval(pt);
	FillTree(d1);
	pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	else sig = fseparation->Eval(pt);
	FillTree(d3);
      }
      else{
	mother = 10;
	part.Decay2body(d4,d6);
	pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	else sig = -fseparation->Eval(pt);
	FillTree(d4);
	pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	else sig = fseparation->Eval(pt);
	FillTree(d6);
      }
    }

    // Lambdac
    // Int_t n = 0;
    // if(gRandom->Rndm() < nlambdac) n = 1;
    for(Int_t j=0;j < nlambdac;j++){
      mother = -1;

      part.ChangeParticleType((gRandom->Rndm() > 0.5)+11);
      FillKine(part,hspectra[6]);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);

      Float_t xvar = gRandom->Rndm();

      mother = part.GetParticleType();

      if(part.GetParticleType() == 11){
	if(xvar < br1){ // no resonant
	  part.Decay3body(d1,d5,d3);
	  pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	  sig = 0;
	  FillTree(d5);
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d3);
	}
	else if(xvar < br1+br2){ // k0*
	  part.Decay2body(d3,res2);
	  
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d3);
	  pt = TMath::Sqrt(res2.GetPx()*res2.GetPx() + res2.GetPy()*res2.GetPy());
	  sig = -999;
	  FillTree(res2);
	  if(gRandom->Rndm() < 0.5){
	    res2.Decay2body(d1,d5);
	    mother = res2.GetParticleType();
	    pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	    if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	    else sig = -fseparation->Eval(pt);
	    FillTree(d1);
	    pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	    sig = 0;
	    FillTree(d5);
	  }
	}
	else if(xvar < br1+br2+br3){ // Delta++
	  part.Decay2body(d5,res3);
	  res3.Decay2body(d1,d3);
	  
	  pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	  sig = 0;
	  FillTree(d5);
	  pt = TMath::Sqrt(res3.GetPx()*res3.GetPx() + res3.GetPy()*res3.GetPy());
	  sig = -999;
	  FillTree(res3);
	  mother = res3.GetParticleType();
	  pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d3);
	}
      }
      else{
	if(xvar < br1){ // no resonant
	  part.Decay3body(d4,d2,d6);
	  pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	  sig = 0;
	  FillTree(d2);
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d6);
	}
	else if(xvar < br1+br2){ // k0s
	  part.Decay2body(d6,res1);
	  
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d6);
	  pt = TMath::Sqrt(res1.GetPx()*res1.GetPx() + res1.GetPy()*res1.GetPy());
	  sig = -999;
	  FillTree(res1);
	  if(gRandom->Rndm() < 0.5){
	    res1.Decay2body(d4,d2);
	    mother = res1.GetParticleType();
	    pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	    if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	    else sig = -fseparation->Eval(pt);
	    FillTree(d4);
	    pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	    sig = 0;
	    FillTree(d2);
	  }
	}
	else if(xvar < br1+br2+br3){ // Delta--
	  part.Decay2body(d2,res4);
	  res4.Decay2body(d4,d6);
	  
	  pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	  sig = 0;
	  FillTree(d2);
	  pt = TMath::Sqrt(res4.GetPx()*res4.GetPx() + res4.GetPy()*res4.GetPy());
	  sig = -999;
	  FillTree(res4);
	  mother = res4.GetParticleType();
	  pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	  if(kALICEseparation) sig = -fseparationPiKa->Eval(pt);
	  else sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  if(kALICEseparation) sig = fseparationKaPr->Eval(pt);
	  else sig = fseparation->Eval(pt);
	  FillTree(d6);
	}
      }
    }


  }

  t->Write();
  fout->Close();

  return 0;
}

void FillTree(particle part){
  sigTPC=-999;
  sigTOF=-999;
  reco = 0;
  id = part.GetParticleType();

  if(id < 2){
    reco = gRandom->Rndm() < fEfficiencyPi->Eval(pt);
    if(reco) reco += gRandom->Rndm() < fEfficiencyPiTOF->Eval(pt);
  }
  else if(id < 4){
    reco = gRandom->Rndm() < fEfficiencyKa->Eval(pt);
    if(reco) reco += gRandom->Rndm() < fEfficiencyKaTOF->Eval(pt);
  }
  else if(id < 6){
    reco = gRandom->Rndm() < fEfficiencyPr->Eval(pt);
    if(reco) reco += gRandom->Rndm() < fEfficiencyPrTOF->Eval(pt);
  }

  sigM = gRandom->Gaus(sig,1);
  pz = part.GetPz();
  Float_t p = TMath::Sqrt(pz*pz + pt*pt);
  eta = 0.5*TMath::Log((p+pz)/(p-pz));


  if(TMath::Abs(eta) > 0.8) reco = 0;

  if(id < 2 && reco){
    sigTPC = fTPCpi->Eval(p);
    if(reco==2)
      sigTOF = fTOFpi->Eval(p,pt/p);
  }
  else if(id < 4 && reco){
    sigTPC = fTPCka->Eval(p);
    if(reco==2)
      sigTOF = fTOFka->Eval(p,pt/p);
  }
  else if(id < 6 && reco){
    sigTPC = fTPCpr->Eval(p);
    if(reco==2)
      sigTOF = fTOFpr->Eval(p,pt/p);
  }
  
  if(reco){
    sigTPC += gRandom->Gaus(0,1);
    if(reco==2)
      sigTOF += gRandom->Gaus(0,1);
  }

  phi = TMath::ATan2(part.GetPy(),part.GetPx());
  if(pt > 0.) t->Fill();
}

void FillKine(particle &part,TH1D *h){
  Float_t phit = gRandom->Rndm()*2*TMath::Pi();
  Float_t pt=h->GetRandom();//-TMath::Log(gRandom->Rndm()) * ptav;

  Float_t y = gRandom->Rndm()*2-1;
  Float_t var = TMath::Exp(2*y);
  var = (var + 1)/(var -1);
  var = 1./(var*var - 1);
  Float_t m = part.GetMass();

  if(y > 0) part.SetP(pt*cos(phit),pt*sin(phit),TMath::Sqrt(var*(m*m + pt*pt)));
  else part.SetP(pt*cos(phit),pt*sin(phit),-TMath::Sqrt(var*(m*m + pt*pt)));

}

Double_t BetheBlochAleph(Double_t *x,Double_t *par) {
  //
  // This is the empirical ALEPH parameterization of the Bethe-Bloch formula.
  // It is normalized to 1 at the minimum.
  //
  // bg - beta*gamma
  // 
  // The default values for the kp* parameters are for ALICE TPC.
  // The returned value is in MIP units
  //

  Double_t bg = x[0]/par[5];
  Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);

  Double_t aa = TMath::Power(beta,par[3]);
  Double_t bb = TMath::Power(1./bg,par[4]);

  bb=TMath::Log(par[2]+bb);
  
  return 14.29*(par[1]-aa-bb)*par[0]/aa;
}
