#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TMath.h"
#include "TLeaf.h"
#include "../particle.h"

TF1 *fseparation;
void ComputeWeights(Float_t weights[3],Float_t signal,Float_t pt);
void GetProb1(Float_t weights[3],Float_t priors[3],Float_t prob[3]);
void GetProb2(Float_t weights1[3],Float_t weights2[3],Float_t priors[3][3],Float_t prob[3][3]);
void analyze(Int_t step=0);

int main(int argc, char* argv[]){
  Int_t nstep=0;
  if(argc > 1){
    sscanf(argv[1],"%d",&nstep);
    printf("Run for step %i\n",nstep);
  }
  analyze(nstep);
  return 0;
}

void analyze(Int_t step){
  /*
  gROOT->LoadMacro("../lessonmaterials/labII/particleType.cxx++");
  gROOT->LoadMacro("../lessonmaterials/labII/resonanceType.cxx++");
  gROOT->LoadMacro("../lessonmaterials/labII/particle.cxx++");
  */
  fseparation = new TF1("f","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.5);
  fseparation->SetParameter(1,10);

  TH1D *priorsPt[6];
  TH1D *newpriorsPt[6];
  TH1D *truePt[6];
  TH1D *allPtPos = new TH1D("allPtP","All positive;p_{T} (GeV/#it{c});N",100,0,10);
  TH1D *allPtNeg = new TH1D("allPtN","All negative;p_{T} (GeV/#it{c});N",100,0,10);

  TH2D *priorsKs[3][3];
  TH2D *priorsPhi[3][3];
  TH2D *newpriorsKs[3][3];
  TH2D *newpriorsPhi[3][3];
  TH2D *truePidKs[3][3];
  TH2D *truePidPhi[3][3];
  TH2D *trueKs;
  TH2D *truePhi;

  const char *spec[3] = {"Pi","Ka","Pr"};

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

    // Ks and phi priors
    for(Int_t i=0; i< 3;i++){
      for(Int_t j=0; j< 3;j++){
	priorsKs[i][j] =  new TH2D(Form("oldpriorsKs%s%s",spec[i],spec[j]),Form("K^{0}_{s} priors for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.8,1,40,0,10);
	priorsPhi[i][j] =  new TH2D(Form("oldpriorsPhi%s%s",spec[i],spec[j]),Form("#phi priors for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);


	if(i==0 && j==0){
	  for(Int_t ibx=1;ibx<=100;ibx++)
	    for(Int_t iby=1;iby<=40;iby++){
	      priorsKs[i][j]->SetBinContent(ibx,iby,1);
	      priorsPhi[i][j]->SetBinContent(ibx,iby,1);
	    }
	}
	else{
	  priorsKs[i][j]->Add(priorsKs[0][0]);
	  priorsPhi[i][j]->Add(priorsPhi[0][0]);
	}
      }
    }
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

    // Ks and phi priors
    for(Int_t i=0; i< 3;i++){
      for(Int_t j=0; j< 3;j++){
	priorsKs[i][j] =  (TH2D *) fin->Get(Form("priorsKs%s%s",spec[i],spec[j]));
	priorsKs[i][j]->SetName(Form("oldpriorsKs%s%s",spec[i],spec[j]));
	priorsPhi[i][j] =  (TH2D *) fin->Get(Form("priorsPhi%s%s",spec[i],spec[j]));
	priorsPhi[i][j]->SetName(Form("oldpriorsPhi%s%s",spec[i],spec[j]));
      }
    }
  }
  
  newpriorsPt[0] = new TH1D("priorsPtPiP","Pion (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[1] = new TH1D("priorsPtKaP","Kaon (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[2] = new TH1D("priorsPtPrP","Proton (+) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[3] = new TH1D("priorsPtPiM","Pion (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[4] = new TH1D("priorsPtKaM","Kaon (-) priors;p_{T} (GeV/#it{c});N",100,0,10);
  newpriorsPt[5] = new TH1D("priorsPtPrM","Proton (-) priors;p_{T} (GeV/#it{c});N",100,0,10);

  // Ks and phi priors distributions
  for(Int_t i=0; i< 3;i++){
    for(Int_t j=0; j< 3;j++){
      newpriorsKs[i][j] =  new TH2D(Form("priorsKs%s%s",spec[i],spec[j]),Form("K^{0}_{s} priors for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.8,1,40,0,10);
      newpriorsPhi[i][j] =  new TH2D(Form("priorsPhi%s%s",spec[i],spec[j]),Form("#phi priors for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);
    }
  }
  
  truePt[0]  = new TH1D("truePtPiP","Pion (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[1]  = new TH1D("truePtKaP","Kaon (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[2]  = new TH1D("truePtPrP","Proton (+) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[3]  = new TH1D("truePtPiM","Pion (-) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[4]  = new TH1D("truePtKaM","Kaon (-) truth;p_{T} (GeV/#it{c});N",100,0,10);
  truePt[5]  = new TH1D("truePtPrM","Proton (-) truth;p_{T} (GeV/#it{c});N",100,0,10);

  // Ks and phi truePid distributions
  for(Int_t i=0; i< 3;i++){
    for(Int_t j=0; j< 3;j++){
      truePidKs[i][j] =  new TH2D(Form("truePidKs%s%s",spec[i],spec[j]),Form("K^{0}_{s} truePid for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.8,1,40,0,10);
      truePidPhi[i][j] =  new TH2D(Form("truePidPhi%s%s",spec[i],spec[j]),Form("#phi truePid for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);
    }
  }

  trueKs =  new TH2D(Form("trueKs"),Form("K^{0}_{s} true;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),100,0.8,1,40,0,10);
  truePhi =  new TH2D(Form("truePhi"),Form("#phi true;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),100,0.98,1.05,40,0,10);
  
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

  particle d1("pi+");
  particle d2("K+");
  particle d3("p+");
  particle d4("pi-");
  particle d5("K-");
  particle d6("p-");

  Int_t charge[] = {1,-1,1,-1,1,-1,0,0,0,2,-2,1,-1};

  particle ppos[20000];
  particle pneg[20000];

  Float_t weightsPos[20000][3];
  Float_t weightsNeg[20000][3];

  Int_t npos=0;
  Int_t nneg=0;

  TFile *f = new TFile("out.root");
  TTree *t = (TTree *) f->Get("tree");
  Int_t n = t->GetEntries()/1000;

  Float_t signal,pt,pz,phi,ptComb,invmass;
  Float_t ptd,pzd,phid;

  Float_t priors[3],prob[3];
  Float_t priors2[3][3],prob2[3][3];
      
  Int_t iev=-1,id;
  Int_t cev;
  for(Int_t i=0;i < n;i++){
    t->GetEvent(i);
    cev = t->GetLeaf("ev")->GetValue();
    id = t->GetLeaf("id")->GetValue();
    pt = t->GetLeaf("pt")->GetValue();

    if(cev != iev){
      //printf("%i\n",cev);
      iev = cev;
      // start the analysis on the combinatorial
      for(Int_t ip=0;ip < npos;ip++){
	for(Int_t jn=0;jn < nneg;jn++){

	  ptComb = (ppos[ip].GetPx() + pneg[jn].GetPx())*(ppos[ip].GetPx() + pneg[jn].GetPx());
	  ptComb += (ppos[ip].GetPy() + pneg[jn].GetPy())*(ppos[ip].GetPy() + pneg[jn].GetPy());
	  ptComb = sqrt(ptComb);

	  Int_t isp1=ppos[ip].GetParticleType(),isp2=pneg[jn].GetParticleType();
	  if(isp1 == 0 || isp1 == 1) isp1 = 0;
	  else if(isp1 == 2 || isp1 == 3) isp1 = 1;
	  else isp1 = 2;
	  if(isp2 == 0 || isp2 == 1) isp2 = 0;
	  else if(isp2 == 2 || isp2 == 3) isp2 = 1;
	  else isp2 = 2;

	  // Ks case
	  d1.SetP(ppos[ip].GetPx(),ppos[ip].GetPy(),ppos[ip].GetPz());
	  d2.SetP(pneg[jn].GetPx(),pneg[jn].GetPy(),pneg[jn].GetPz());
	  invmass = d1.InvMass(d2);

	  if(invmass > 0.8 && invmass < 1 && ptComb < 10){ 
	    for(Int_t ipr=0;ipr<3;ipr++)
	      for(Int_t jpr=0;jpr<3;jpr++)
		priors2[ipr][jpr] = priorsKs[ipr][jpr]->Interpolate(invmass,ptComb);

	    truePidKs[isp1][isp2]->Fill(invmass,ptComb);

	    GetProb2(weightsPos[ip],weightsNeg[jn],priors2,prob2);
	    for(Int_t ipr=0;ipr<3;ipr++)
	      for(Int_t jpr=0;jpr<3;jpr++)
		newpriorsKs[ipr][jpr]->Fill(invmass,ptComb,prob2[ipr][jpr]);
	  }

	  // Phi case
	  d2.SetP(ppos[ip].GetPx(),ppos[ip].GetPy(),ppos[ip].GetPz());
	  d5.SetP(pneg[jn].GetPx(),pneg[jn].GetPy(),pneg[jn].GetPz());
	  invmass = d2.InvMass(d5);

	  if(invmass > 0.98 && invmass < 1.05 && ptComb < 10){ 

// 	    if(ptComb < 0.5){
// 	      printf("%i %i -> %f\n",ip,jn,invmass);
// 	      ppos[ip].Print();
// 	      pneg[jn].Print();
// 	      getchar();
// 	    }

	    for(Int_t ipr=0;ipr<3;ipr++)
	      for(Int_t jpr=0;jpr<3;jpr++){
		priors2[ipr][jpr] = priorsPhi[ipr][jpr]->Interpolate(invmass,ptComb);
		//	printf("%f\n",priors2[ipr][jpr] );
	      }
	    truePidPhi[isp1][isp2]->Fill(invmass,ptComb);

	    GetProb2(weightsPos[ip],weightsNeg[jn],priors2,prob2);
	    
	    for(Int_t ipr=0;ipr<3;ipr++)
	      for(Int_t jpr=0;jpr<3;jpr++)
		newpriorsPhi[ipr][jpr]->Fill(invmass,ptComb,prob2[ipr][jpr]);
	  }
	}
      }

      // reset variables
      npos=0;
      nneg=0;
      i--;
    }
    else if(id < 6){
      signal = t->GetLeaf("signal")->GetValue();
      pz = t->GetLeaf("pz")->GetValue();
      phi = t->GetLeaf("phi")->GetValue();

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
	allPtPos->Fill(pt);
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
	if(id==1) truePt[3]->Fill(pt);
	else if(id==3) truePt[4]->Fill(pt);
	else if(id==5) truePt[5]->Fill(pt);
	allPtNeg->Fill(pt);
	pneg[nneg].ChangeParticleType(id);
	pneg[nneg].SetP(pt*cos(phi),pt*sin(phi),pz);
	nneg++;
      }
      // perform the analysis on the single species
    }
    else if(id == 6 || id == 7){ // is K0s or a K0sb
      t->GetEvent(i+1);
      ptd = t->GetLeaf("pt")->GetValue();
      pzd = t->GetLeaf("pz")->GetValue();
      phid = t->GetLeaf("phi")->GetValue();
      d1.SetP(ptd*cos(phid),ptd*sin(phid),pzd);
      t->GetEvent(i+2);
      ptd = t->GetLeaf("pt")->GetValue();
      pzd = t->GetLeaf("pz")->GetValue();
      phid = t->GetLeaf("phi")->GetValue();
      d2.SetP(ptd*cos(phid),ptd*sin(phid),pzd);
      trueKs->Fill(d1.InvMass(d2),pt);
    }
    else if(id == 8){ // is a Phi
      t->GetEvent(i+1);
      ptd = t->GetLeaf("pt")->GetValue();
      pzd = t->GetLeaf("pz")->GetValue();
      phid = t->GetLeaf("phi")->GetValue();
      d2.SetP(ptd*cos(phid),ptd*sin(phid),pzd);
      t->GetEvent(i+2);
      ptd = t->GetLeaf("pt")->GetValue();
      pzd = t->GetLeaf("pz")->GetValue();
      phid = t->GetLeaf("phi")->GetValue();
      d5.SetP(ptd*cos(phid),ptd*sin(phid),pzd);
      truePhi->Fill(d2.InvMass(d5),pt);
    }
  }

  for(Int_t i=0;i<6;i++){
    if(step == 0){
      priorsPt[i]->Reset();
      if(i < 3) priorsPt[i]->Add(allPtPos);
      else priorsPt[i]->Add(allPtNeg);
      priorsPt[i]->Scale(1./3);      
    }
    
    priorsPt[i]->Add(newpriorsPt[i],-1);

    for(Int_t j=1;j<=100;j++){
      if(newpriorsPt[i]->GetBinContent(j)-priorsPt[i]->GetBinContent(j)*0.5 > 0)
	newpriorsPt[i]->SetBinContent(j,newpriorsPt[i]->GetBinContent(j)-priorsPt[i]->GetBinContent(j)*0.5);
      if(i < 3) newpriorsPt[i]->SetBinError(j,priorsPt[i]->GetBinContent(j)/newpriorsPt[i]->GetBinContent(j)*allPtPos->GetBinContent(j));
      else       newpriorsPt[i]->SetBinError(j,priorsPt[i]->GetBinContent(j)/newpriorsPt[i]->GetBinContent(j)*allPtNeg->GetBinContent(j));

    }
  }

  TFile *fout = new TFile(Form("step%i.root",step+1),"RECREATE");
  for(Int_t i=0;i<6;i++){
    newpriorsPt[i]->Write();
    truePt[i]->Write();
  }
  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<3;j++){
      priorsPhi[i][j]->Write();
      newpriorsKs[i][j]->Write();
      newpriorsPhi[i][j]->Write();
      truePidKs[i][j]->Write();
      truePidPhi[i][j]->Write();
    }
  }
  trueKs->Write();
  truePhi->Write();
      
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

void GetProb2(Float_t weights1[3],Float_t weights2[3],Float_t priors[3][3],Float_t prob[3][3]){
  Float_t R= 0;
  for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
    prob[i][j] = weights1[i]*weights2[j]*priors[i][j];
    R += prob[i][j];
    }
  }
  R = 1./R;
  Float_t R2=0;
   for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
      prob[i][j] *= R;
      R2 += prob[i][j];
      //if(prob[i][j] > 1) printf("why? %i-%i = %f\n",i,j,prob[i][j]);
      //if(prob[i][j] < 0.001) prob[i][j] = 0;//printf("why so low? %i-%i = %f\n",i,j,prob[i][j]);
    }
   }
}
