TF1 *fseparation;

analyze(Int_t step=1){
  gROOT->LoadMacro("../lessonmaterials/labII/particleType.cxx++");
  gROOT->LoadMacro("../lessonmaterials/labII/resonanceType.cxx++");
  gROOT->LoadMacro("../lessonmaterials/labII/particle.cxx++");

  fseparation = new TF1("f","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.5);
  fseparation->SetParameter(1,10);

  TH1D *priorsPt[6];
  TH1D *newpriorsPt[6];
  TH1D *truePt[6];

  if(step==0){
    priorsPt[0] = new TH1D("oldpriorsPtPiP","Pion (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
    for(Int_t i=1;i<=100;i++)
      priorsPt[0]->SetBinContent(i,1);
    priorsPt[1] = new TH1D("oldpriorsPtKaP","Kaon (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
    priorsPt[2] = new TH1D("oldpriorsPtPrP","Proton (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
    priorsPt[3] = new TH1D("oldpriorsPtPiM","Pion (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
    priorsPt[4] = new TH1D("oldpriorsPtKaM","Kaon (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
    priorsPt[5] = new TH1D("oldpriorsPtPrM","Proton (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
    priorsPt[1]->Add(priorsPt[0]);
    priorsPt[2]->Add(priorsPt[0]);
    priorsPt[3]->Add(priorsPt[0]);
    priorsPt[4]->Add(priorsPt[0]);
    priorsPt[5]->Add(priorsPt[0]);
  }
  else{
    TFile *fin = new TFile(Form("step%i.root",step));
    priorsPt[0] = (TH1D *) fin->Get("priorsPtPiP"); 
    priorsPt[0]->SetName("oldpriorsPtPiP");
    priorsPt[1] = (TH1D *) fin->Get("priorsPtKaP"); 
    priorsPt[1]->SetName("oldpriorsPtPiP");
    priorsPt[2] = (TH1D *) fin->Get("priorsPtPrP"); 
    priorsPt[2]->SetName("oldpriorsPtPiP");
    priorsPt[3] = (TH1D *) fin->Get("priorsPtPiM"); 
    priorsPt[3]->SetName("oldpriorsPtPiM");
    priorsPt[4] = (TH1D *) fin->Get("priorsPtKaM"); 
    priorsPt[4]->SetName("oldpriorsPtKaM");
    priorsPt[5] = (TH1D *) fin->Get("priorsPtPrM"); 
    priorsPt[5]->SetName("oldpriorsPtPrM");
  }
  
  newpriorsPt[0] = new TH1D("priorsPtPiP","Pion (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[1] = new TH1D("priorsPtKaP","Kaon (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[2] = new TH1D("priorsPtPrP","Proton (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[3] = new TH1D("priorsPtPiM","Pion (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[4] = new TH1D("priorsPtKaM","Kaon (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[5] = new TH1D("priorsPtPrM","Proton (-) priors;p_{T} (GeV/#it{c});N",100,0,10);

  truePt[0]  = new TH1D("truePtPiP","Pion (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[1]  = new TH1D("truePtKaP","Kaon (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[2]  = new TH1D("truePtPrP","Proton (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[3]  = new TH1D("truePtPiM","Pion (-) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[4]  = new TH1D("truePtKaM","Kaon (-) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[5]  = new TH1D("truePtPrM","Proton (-) truth;p_{T} (GeV/#it{c});N",100,0,10);

  particle::AddParticleType("pi+",0.139,1);
  particle::AddParticleType("pi-",0.139,-1);
  particle::AddParticleType("K+",0.493,1);  
  particle::AddParticleType("K-",0.493,-1);
  particle::AddParticleType("p+",0.938,1);  
  particle::AddParticleType("p-",0.938,-1);
  particle::AddParticleType("K0*",0.896,0,5.05e-02/2.355);
  particle::AddParticleType("K0bar*",0.896,0,5.05e-02/2.355);
  particle::AddParticleType("Phi",1.02,0,0.00426/2.355);
  particle::AddParticleType("Delta++",1.232,2,0.118/2.355);
  particle::AddParticleType("Delta--",1.232,-2,0.118/2.355);
  particle::AddParticleType("Lambdac+",2.28646,1,0.08/2.355);
  particle::AddParticleType("Lambdacbar+",2.28646,-1,0.08/2.355);

  Int_t charge[] = {1,-1,1,-1,1,-1,0,0,0,2,-2,1,-1};

  particle ppos[20000];
  particle pneg[20000];

  Float_t weightsPos[20000][3];
  Float_t weightsNeg[20000][3];

  Int_t npos=0;
  Int_t nneg=0;

  TFile *f = new TFile("out.root");
  TTree *t = f->Get("tree");
  Int_t n = t->GetEntries()/10;

  Float_t signal,pt,pz,phi;

  Int_t iev=-1,id;
  Int_t cev;
  for(Int_t i=0;i < n;i++){
    t->GetEvent(i);
    cev = t->GetLeaf("ev")->GetValue();
    id = t->GetLeaf("id")->GetValue();

    if(cev != iev){
      //printf("%i\n",cev);
      iev = cev;

      // start the analysis on the combinatorial


      // reset variables
      npos=0;
      nneg=0;
      i--;
    }
    else if(id < 6){
      signal = t->GetLeaf("signal")->GetValue();
      pt = t->GetLeaf("pt")->GetValue();
      pz = t->GetLeaf("pz")->GetValue();
      phi = t->GetLeaf("phi")->GetValue();

      Float_t priors[3],prob[3];

      if(charge[id] > 0){
	ComputeWeights(weightsPos[npos],signal,pt);
	priors[0] = priorsPt[0]->Interpolate(pt);
	priors[1] = priorsPt[1]->Interpolate(pt);
	priors[2] = priorsPt[2]->Interpolate(pt);
	GetProb1(weightsPos[npos],priors,prob);
	newpriorsPt[0]->Fill(pt,prob[0]);
	newpriorsPt[1]->Fill(pt,prob[1]);
	newpriorsPt[2]->Fill(pt,prob[2]);
	if(id==0) truePt[0]->Fill(pt);
	else if(id==2) truePt[1]->Fill(pt);
	else if(id==4) truePt[2]->Fill(pt);
	ppos[npos].ChangeParticleType(id);
	ppos[npos].SetP(pt*cos(phi),pt*sin(phi),pz);
	npos++;
      }
      else{
	ComputeWeights(weightsNeg[nneg],signal,pt);
	priors[0] = priorsPt[3]->Interpolate(pt);
	priors[1] = priorsPt[4]->Interpolate(pt);
	priors[2] = priorsPt[5]->Interpolate(pt);
	GetProb1(weightsNeg[nneg],priors,prob);
	newpriorsPt[3]->Fill(pt,prob[0]);
	newpriorsPt[4]->Fill(pt,prob[1]);
	newpriorsPt[5]->Fill(pt,prob[2]);
	if(id==1) truePt[0]->Fill(pt);
	else if(id==3) truePt[1]->Fill(pt);
	else if(id==5) truePt[2]->Fill(pt);
	pneg[npos].ChangeParticleType(id);
	pneg[npos].SetP(pt*cos(phi),pt*sin(phi),pz);
	nneg++;
      }
      // perform the analysis on the single species
    }
  }

  TFile *fout = new TFile(Form("step%i.root",step+1),"RECREATE");
  for(Int_t i=0;i<6;i++){
    newpriorsPt[i]->Write();
    truePt[i]->Write();
  }

  fout->Close();
}

void ComputeWeights(Float_t weights[3],Float_t signal,Float_t pt){
  Float_t expect = fseparation->Eval(pt);

  weights[0] = TMath::Exp(-(signal+expect)*(signal+expect)*0.5);
  weights[1] = TMath::Exp(-signal*signal*0.5);
  weights[2] = TMath::Exp(-(signal-expect)*(signal-expect)*0.5);
}

void GetProb1(Float_t weights[3],Float_t priors[3],Float_t prob[3]){
    prob[0] = weights[0]*priors[0];
    prob[1] = weights[1]*priors[1];
    prob[2] = weights[2]*priors[2];
    Float_t R= prob[0]+prob[1]+prob[2];
    R = 1./R;
    prob[0] *= R;
    prob[1] *= R;
    prob[2] *= R;
}
