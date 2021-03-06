{
  gROOT->LoadMacro("particleType.cxx++");
  gROOT->LoadMacro("resonanceType.cxx++");
  gROOT->LoadMacro("particle.cxx++");

  particle::AddParticleType("pi+",0.139,1);
  particle::AddParticleType("pi-",0.139,-1);
  particle::AddParticleType("K+",0.493,1);  
  particle::AddParticleType("K-",0.493,-1);
  particle::AddParticleType("p+",0.938,1);  
  particle::AddParticleType("p-",0.938,-1);
  particle::AddParticleType("K0*",0.896,0,5.05e-02/2.355);
  particle::AddParticleType("Delta++",1.232,2,0.118/2.355);
  particle::AddParticleType("Lambdac+",2.28646,1,0.08/2.355);
 
  particle moth("Lambdac+");
  particle d1("pi+");
  particle d2("K-");
  particle d3("p+");
  particle res1("K0*");
  particle res2("Delta++");

  TH2D *hdalitz = new TH2D("hdalitz",";m_{#piK}^2;m_{#pip}^2",200,0,4,200,1,9);
  TH2D *hdalitz2 = new TH2D("hdalitz2",";m_{K#pi}^2;m_{Kp}^2",200,0,4,200,1,9);
  TH2D *hdalitz3 = new TH2D("hdalitz3",";m_{p#pi}^2;m_{pK}^2",200,0,4,200,1,9);
  TH1D *hmass = new TH1D("hmass","m_{\piKp}",100,1,3);

  Float_t br1 = 2.8/5.26;
  Float_t br2 = 1.6/5.26;
  Float_t br3 = 0.86/5.26;

  for(Int_t i=0;i < 10000*br1;i++){

    moth->Decay3body(d1,d2,d3);

    hdalitz->Fill(d1.InvMass(d2)**2,d1.InvMass(d3)**2);
    hdalitz2->Fill(d2.InvMass(d1)**2,d2.InvMass(d3)**2);
    hdalitz3->Fill(d3.InvMass(d1)**2,d3.InvMass(d2)**2);
    hmass->Fill(d1.InvMass(d2,d3));
  }


  for(Int_t i=0;i < 10000*br2;i++){

    moth->Decay2body(d3,res1);
    res1->Decay2body(d1,d2);

    hdalitz->Fill(d1.InvMass(d2)**2,d1.InvMass(d3)**2);
    hdalitz2->Fill(d2.InvMass(d1)**2,d2.InvMass(d3)**2);
    hdalitz3->Fill(d3.InvMass(d1)**2,d3.InvMass(d2)**2);
    hmass->Fill(d1.InvMass(d2,d3));
  }

  for(Int_t i=0;i < 10000*br3;i++){

    moth->Decay2body(d2,res2);
    res2->Decay2body(d1,d3);

    hdalitz->Fill(d1.InvMass(d2)**2,d1.InvMass(d3)**2);
    hdalitz2->Fill(d2.InvMass(d1)**2,d2.InvMass(d3)**2);
    hdalitz3->Fill(d3.InvMass(d1)**2,d3.InvMass(d2)**2);
    hmass->Fill(d1.InvMass(d2,d3));
  }

  TCanvas *c = new TCanvas();
  c->Divide(3,1);
  c->cd(1);
  hdalitz->Draw("colz");
  c->cd(2);
  hdalitz2->Draw("colz");
  c->cd(3);
  hdalitz3->Draw("colz");
  new TCanvas();
  hmass->Draw();
}
