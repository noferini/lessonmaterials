{ 

  gROOT->LoadMacro("particleType.cxx++");
  gROOT->LoadMacro("resonanceType.cxx++");
  gROOT->LoadMacro("particle.cxx++");


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

  particle *part = new particle();
  
  TTree *tree = new TTree("part","part");
  // tree->Branch("isp",&isp,"isp/I");
  // tree->Branch("p",p,"p[3]/F");
  tree->Branch("part","particle",part);

  TF1 *fexp = new TF1("fexp","expo",0,4);
  fexp->SetParameter(0,1);
  fexp->SetParameter(1,-1);
  fexp->SetNpx(1000);

  Float_t rr;
  Float_t ptot,phi,theta;
  for(Int_t i=0;i < 10000;i++){
    rr = gRandom->Rndm();
    if(rr < 0.5)isp=0;
    else isp=1;

    ptot = fexp->GetRandom(0,4);
    phi = gRandom->Rndm()*TMath::Pi()*2;
    theta = (gRandom->Rndm()-0.5)*TMath::Pi();
    
    p[0] = ptot * TMath::Cos(phi) * TMath::Cos(theta);
    p[1] = ptot * TMath::Sin(phi) * TMath::Cos(theta);
    p[2] = ptot * TMath::Sin(theta);

    part->ChangeParticleType(isp);
    part->SetP(p[0],p[1],p[2]);
    
    tree->Fill();
  }


  TFile *fout = new TFile("treeour.root","RECREATE");
  tree->Write();
  fout->Close();

}

