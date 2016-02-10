#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TTree.h"

#include "particle.h"

TF1 *fseparation;

TTree *t;
Float_t pt,pz,phi;
Int_t id;
Float_t sig;
Float_t sigM;
Int_t mother;
Float_t vxy;
Int_t tof;
Int_t ev=0;

void FillTree(particle part);
void FillKine(particle &part);

int main(){
  /*
  gROOT->LoadMacro("particleType.cxx++");
  gROOT->LoadMacro("resonanceType.cxx++");
  gROOT->LoadMacro("particle.cxx++");
  */

  gRandom->SetSeed(0);

  fseparation = new TF1("f","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.5);
  fseparation->SetParameter(1,10);

  t = new TTree("tree","tree");
  t->Branch("ev",&ev,"ev/I");
  t->Branch("id",&id,"id/I");
  t->Branch("pt",&pt,"pt/F");
  t->Branch("pz",&pz,"pz/F");
  t->Branch("phi",&phi,"phi/F");
  t->Branch("signal",&sigM,"signal/F");
  t->Branch("mother",&mother,"mother/I");
  t->Branch("vxy",&vxy,"vxy/F");
  t->Branch("tof",&tof,"tof/I");

  particle::AddParticleType("pi+",0.139,1);
  particle::AddParticleType("pi-",0.139,-1);
  particle::AddParticleType("K+",0.493,1);  
  particle::AddParticleType("K-",0.493,-1);
  particle::AddParticleType("p+",0.938,1);  
  particle::AddParticleType("p-",0.938,-1);
  particle::AddParticleType("K0*",0.896,0,5.05e-02);
  particle::AddParticleType("K0bar*",0.896,0,5.05e-02);
  particle::AddParticleType("Phi",1.02,0,0.00426);
  particle::AddParticleType("Delta++",1.232,2,0.118);
  particle::AddParticleType("Delta--",1.232,-2,0.118);
  particle::AddParticleType("Lambdac+",2.28646,1,0.08);
  particle::AddParticleType("Lambdacbar+",2.28646,1,0.08);


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

  Float_t br1 = 0.028; // Lambda_c -> pi K p no resonant
  Float_t br2 = 0.016; // Lambda_c -> K0* p -> pi K p
  Float_t br3 = 0.0086; // Lambda_c -> K Delta -> pi K p 

  particle part;

  for(ev=0;ev < 10000;ev++){ // event loop
    //    t->Reset();

    mother = -1;
    vxy = 0;
    tof = 1;

    // pions
    for(Int_t j=0;j < 1000;j++){
      part.ChangeParticleType(gRandom->Rndm() > 0.5);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -fseparation->Eval(pt);
      FillTree(part);
    }

    // kaons
    for(Int_t j=0;j < 100;j++){
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+2);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = 0;
      FillTree(part);
    }

    // protons
    for(Int_t j=0;j < 75;j++){
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+4);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = fseparation->Eval(pt);
      FillTree(part);
    }

    // K0s
    for(Int_t j=0;j < 7;j++){
      mother = -1;

      part.ChangeParticleType((gRandom->Rndm() > 0.5)+6);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);
      
      if(part.GetParticleType() == 7){
	mother = 7;
	part.Decay2body(d1,d5);
	pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	sig = -fseparation->Eval(pt);
	FillTree(d1);
	pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	sig = 0;
	FillTree(d5);
      }
      else{
	mother = 6;
	part.Decay2body(d4,d2);
	pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	sig = -fseparation->Eval(pt);
	FillTree(d4);
	pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	sig = 0;
	FillTree(d2);
      }

    }

    // Phi
    for(Int_t j=0;j < 2;j++){
      mother = -1;
      part.ChangeParticleType(8);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = -999;
      FillTree(part);
      Float_t xran = gRandom->Rndm();
      if(xran < 0.492){
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
    for(Int_t j=0;j < 1;j++){
      mother = -1;
      part.ChangeParticleType((gRandom->Rndm() > 0.5)+9);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = fseparation->Eval(pt);
      FillTree(part);

      if(part.GetParticleType() == 9){
	mother = 9;
	part.Decay2body(d1,d3);
	pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	sig = -fseparation->Eval(pt);
	FillTree(d1);
	pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	sig = fseparation->Eval(pt);
	FillTree(d3);
      }
      else{
	mother = 10;
	part.Decay2body(d4,d6);
	pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	sig = -fseparation->Eval(pt);
	FillTree(d4);
	pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	sig = fseparation->Eval(pt);
	FillTree(d6);
      }
    }

    // Lambdac
    Int_t n = 0;
    if(gRandom->Rndm() < 0.25) n = 1;
    for(Int_t j=0;j < n;j++){
      mother = -1;

      part.ChangeParticleType((gRandom->Rndm() > 0.5)+11);
      FillKine(part);
      pt = TMath::Sqrt(part.GetPx()*part.GetPx() + part.GetPy()*part.GetPy());
      sig = fseparation->Eval(pt);
      FillTree(part);

      Float_t xvar = gRandom->Rndm();

      mother = part.GetParticleType();

      if(part.GetParticleType() == 11){
	if(xvar < br1){ // no resonant
	  part.Decay3body(d1,d5,d3);
	  pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	  sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	  sig = 0;
	  FillTree(d5);
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d3);
	}
	else if(xvar < br1+br2){ // k0s
	  part.Decay2body(d3,res2);
	  res2.Decay2body(d1,d5);
	  
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d3);
	  pt = TMath::Sqrt(res2.GetPx()*res2.GetPx() + res2.GetPy()*res2.GetPy());
	  sig = -999;
	  FillTree(res2);
	  mother = res2.GetParticleType();
	  pt = TMath::Sqrt(d1.GetPx()*d1.GetPx() + d1.GetPy()*d1.GetPy());
	  sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d5.GetPx()*d5.GetPx() + d5.GetPy()*d5.GetPy());
	  sig = 0;
	  FillTree(d5);
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
	  sig = -fseparation->Eval(pt);
	  FillTree(d1);
	  pt = TMath::Sqrt(d3.GetPx()*d3.GetPx() + d3.GetPy()*d3.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d3);
	}
      }
      else{
	if(xvar < br1){ // no resonant
	  part.Decay3body(d4,d2,d6);
	  pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	  sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	  sig = 0;
	  FillTree(d2);
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d6);
	}
	else if(xvar < br1+br2){ // k0s
	  part.Decay2body(d6,res1);
	  res1.Decay2body(d4,d2);
	  
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d6);
	  pt = TMath::Sqrt(res1.GetPx()*res1.GetPx() + res1.GetPy()*res1.GetPy());
	  sig = -999;
	  FillTree(res1);
	  mother = res1.GetParticleType();
	  pt = TMath::Sqrt(d4.GetPx()*d4.GetPx() + d4.GetPy()*d4.GetPy());
	  sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d2.GetPx()*d2.GetPx() + d2.GetPy()*d2.GetPy());
	  sig = 0;
	  FillTree(d2);
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
	  sig = -fseparation->Eval(pt);
	  FillTree(d4);
	  pt = TMath::Sqrt(d6.GetPx()*d6.GetPx() + d6.GetPy()*d6.GetPy());
	  sig = fseparation->Eval(pt);
	  FillTree(d6);
	}
      }
    }


  }

  TFile *fout = new TFile("out.root","RECREATE");
  t->Write();
  fout->Close();

  return 0;
}

void FillTree(particle part){
  id = part.GetParticleType();
  sigM = gRandom->Gaus(sig,1);
  pz = part.GetPz();
  phi = TMath::ATan2(part.GetPy(),part.GetPx());
  if(pt > 0.) t->Fill();
}

void FillKine(particle &part){
  Float_t phit = gRandom->Rndm()*2*TMath::Pi();
  Float_t pt=-5*TMath::Log(gRandom->Rndm());

  Float_t y = gRandom->Rndm()*2-1;
  Float_t var = TMath::Exp(2*y);
  var = (var + 1)/(var -1);
  var = 1./(var*var - 1);
  Float_t m = part.GetMass();

  part.SetP(pt*cos(phit),pt*sin(phit),var * TMath::Sqrt(m*m + pt*pt));
}
