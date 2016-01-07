{ 

  gROOT->LoadMacro("particleType.cxx++");
  gROOT->LoadMacro("resonanceType.cxx++");
  gROOT->LoadMacro("particle.cxx++");
  gROOT->LoadMacro("event.cxx++");


  // massa in GeV/c^2
  // Carica in unita'di cariche elementari (1 => 1.6*10E-19 C) 
  particle::AddParticleType("pi+",0.139,1);
  particle::AddParticleType("pi-",0.139,-1);
  particle::AddParticleType("K+",0.493,1);  
  particle::AddParticleType("K-",0.493,-1);
  particle::AddParticleType("p+",0.938,1);  
  particle::AddParticleType("p-",0.938,-1);
  particle::AddParticleType("K0*",0.896,0,5.05e-02);
  
  // stampa le informazioni dei tipi di particella
  particle::PrintParticleType();


  Int_t isp;
  Float_t p[3];

  event *eve = new event();
  
  TFile *f = new TFile("treecomplex.root");
  TTree *tree = (TTree *) f->Get("part");
  tree->SetBranchAddress("event",&eve); // indirizzo di puntatore

  Int_t nev = tree->GetEntries();

  TH1D *h = new TH1D("hptot","p_{tot};p (GeV/c);entries",100,0,4);

  Float_t ptot,phi;
  for(Int_t i=0;i < nev;i++){
    tree->GetEvent(i);
    for(Int_t j=0;j < eve->GetEntries();j++){
      const particle *part = eve->GetParticle(j);
      ptot = TMath::Sqrt(part->GetPx()*part->GetPx() + part->GetPy()*part->GetPy() + part->GetPz()*part->GetPz());
      phi = TMath::ATan2(part->GetPy(),part->GetPx());
      
      h->Fill(ptot);
    }
  }

  h->Draw();

}


