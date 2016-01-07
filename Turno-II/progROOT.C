{ 

  gROOT->LoadMacro("particleType.cxx++g");
  gROOT->LoadMacro("resonanceType.cxx++g");
  gROOT->LoadMacro("particle.cxx++g");


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

  // esempio di particelle
  particle p1("K0*",1.0,0.0,0.0); 
  particle p2("pi-");
  particle p3("K+");

  p1.Decay2body(p2,p3);

  p1.Print(); 
  p2.Print(); 
  p3.Print(); 

  printf("Inv mass pi-,K+ = %f\n",p2.InvMass(p3));

  p1.Decay2body(p2,p3);

  p1.Print(); 
  p2.Print(); 
  p3.Print(); 

  printf("Inv mass pi-,K+ = %f\n",p2.InvMass(p3));

}

