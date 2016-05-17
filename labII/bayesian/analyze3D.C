#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TRandom.h"
#include "TList.h"

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

Float_t ptminPi = 0.3;
Float_t ptminKa = 0.4;
Float_t ptminPr = 0.6;

// bethe block parameter
Float_t fKp1=0.0283086;
Float_t fKp2=2.63394e+01;
Float_t fKp3=5.04114e-11;
Float_t fKp4=2.12543;
Float_t fKp5=4.88663;
Double_t BetheBlochAleph(Double_t *x,Double_t *par);

Int_t lambdacGood(Float_t ptLc,Float_t ptPi,Float_t ptKa,Float_t ptPr);

Bool_t kALICEseparation=kTRUE;

// return the weight array for pion, kaon and proton hypoteses
// conditioned probability: a given particles releases the measured signal
void ComputeWeights(Float_t weights[3],Float_t signal,Float_t pt);
void ComputeWeightsALICE(Float_t weights[3],Float_t signalTPC,Float_t signalTOF,Float_t pt,Float_t p);

// return the probability array for pion, kaon and proton hypoteses
void GetProb1(Float_t weights[3],Float_t priors[3],Float_t prob[3]);

// return the probability array for all the pairs 3x3: (pi,pi), (pi,ka), ... , (pr,ka), (pr,pr)
void GetProb2(Float_t weights1[3],Float_t weights2[3],Float_t priors[3][3],Float_t prob[3][3]);

void GetProb3(Float_t weights1[3],Float_t weights2[3],Float_t weights3[3],Float_t priors[3][3][3],Float_t prob[3][3][3]);

void analyze(Int_t step=0);

Float_t GetPolariz(particle moth,particle dau);

// method you can use to manage your own PID
Int_t passMyPID(Int_t isp,Float_t p,Float_t pt,Float_t TPCsign,Float_t TOFsign);
Int_t 




















passMyPID2(Int_t isp,Float_t p,Float_t pt,Float_t TPCsign,Float_t TOFsign);

int main(int argc, char* argv[]){
  Int_t nstep=0;
  if(argc > 1){
    sscanf(argv[1],"%d",&nstep);
    printf("Run for step %i\n",nstep);
  }
  analyze(nstep);
  return 0;
}

Float_t invwidth;
Float_t addshift;
Float_t invwidthTOF;
Float_t addshiftTOF;

TF1 *fEfficiencyPiTOF;
TF1 *fEfficiencyKaTOF;
TF1 *fEfficiencyPrTOF;

void analyze(Int_t step){

  // TOF propagation factors (TOF efficiencies)
  fEfficiencyPiTOF = new TF1("fEfficiencyPiTOF","(x > 0.3)*0.7",0,10);
  fEfficiencyKaTOF = new TF1("fEfficiencyKaTOF","(x > 0.3)*(x-0.3)*(x<1) + (x>1)*0.7",0,10);
  fEfficiencyPrTOF = new TF1("fEfficiencyPrTOF","(x > 0.3)*(x-0.3)*(x<1) + (x>1)*0.7",0,10);

  // teoretical separation (perfect if equal to the one simualted in sim.C)
  fseparation = new TF1("f","[0]+[1]/x",0,100);
  fseparation->SetParameter(0,0.);
  fseparation->SetParameter(1,7.);
 
  fseparationPiKa = new TF1("fPiKa","[0]+[1]/TMath::Power(x,2.5)",0,100);
  fseparationPiKa->SetParameter(0,2.34);
  fseparationPiKa->SetParameter(1,10);

  fseparationKaPr = new TF1("fKaPr","[0]+[1]/TMath::Power(x,2.5)",0,100);
  fseparationKaPr->SetParameter(0,1);
  fseparationKaPr->SetParameter(1,56);

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

  Float_t width = 1.0;
  addshift =0;

  invwidth = 1./width;

  Float_t widthTOF = 1.0;
  addshiftTOF =0;

  invwidthTOF = 1./widthTOF;

  TH1D *priorsPt[6];
  TH1D *newpriorsPt[6];
  TH1D *truePt[6];
  TH1D *allPtPos = new TH1D("allPtP","All positive;p_{T} (GeV/#it{c});N",100,0,10);
  TH1D *allPtNeg = new TH1D("allPtN","All negative;p_{T} (GeV/#it{c});N",100,0,10);

  TH3D *priorsKs[3][3];
  TH3D *newpriorsKs[3][3];
  TH3D *truePidKs[3][3];
  TH3D *trueKs;

  TH2D *priorsPhi[3][3];
  TH2D *newpriorsPhi[3][3];
  TH2D *truePidPhi[3][3];
  TH2D *truePhi;

  TH3D *priorsLc[3][3][3];
  TH3D *newpriorsLc[3][3][3];
  TH3D *truePidLc[3][3][3];
  TH3D *trueLc,*mypidLc;

  TH3D *priorsLcbar[3][3][3];
  TH3D *newpriorsLcbar[3][3][3];
  TH3D *truePidLcbar[3][3][3];
  TH3D *trueLcbar,*mypidLcbar;

  Int_t nbinPtFrKa = 8;
  Int_t nbinPtFrPi = 8;
  Int_t nbinY = 1;
  Int_t nbinpol=nbinPtFrKa*nbinPtFrPi*nbinY;
  Double_t normbin = 1./nbinpol;
  Int_t nbinmlc = 100;
  Int_t nbinptlc = 10;

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
	priorsKs[i][j] =  new TH3D(Form("oldpriorsKs%s%s",spec[i],spec[j]),Form("K^{0*} priors for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),200,0.4,1.4,40,0,10,nbinpol,-1.001,1.001);
	priorsPhi[i][j] =  new TH2D(Form("oldpriorsPhi%s%s",spec[i],spec[j]),Form("#phi priors for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);


	if(i==0 && j==0){
	  for(Int_t ibx=1;ibx<=200;ibx++)
	    for(Int_t iby=1;iby<=40;iby++){
	      for(Int_t ibz=1;ibz <= nbinpol;ibz++)
		priorsKs[i][j]->SetBinContent(ibx,iby,ibz,1);
	      priorsPhi[i][j]->SetBinContent(ibx,iby,1);
	    }
	}
	else{
	  priorsKs[i][j]->Add(priorsKs[0][0]);
	  priorsPhi[i][j]->Add(priorsPhi[0][0]);
	}

	for(Int_t k=0; k< 3;k++){
	  priorsLc[i][j][k] =  new TH3D(Form("oldpriorsLc%s%s%s",spec[i],spec[j],spec[k]),Form("#Lambda_{c}^{+} priors for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
	  priorsLcbar[i][j][k] =  new TH3D(Form("oldpriorsLcbar%s%s%s",spec[i],spec[j],spec[k]),Form("#overline{#Lambda}_{c}^{-} priors for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
	  if(i==0 && j==0 && k==0){
	    for(Int_t ibx=1;ibx<=nbinmlc;ibx++)
	      for(Int_t iby=1;iby<=nbinptlc;iby++){
		for(Int_t ibz=1;ibz <= nbinpol;ibz++){
		  priorsLc[i][j][k]->SetBinContent(ibx,iby,ibz,1);
		  priorsLcbar[i][j][k]->SetBinContent(ibx,iby,ibz,1);
		}
	      }
	  }
	  else{
	    priorsLc[i][j][k]->Add(priorsLc[0][0][0]);
	    priorsLcbar[i][j][k]->Add(priorsLc[0][0][0]);
	  }
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
	priorsKs[i][j] =  (TH3D *) fin->Get(Form("priorsKs%s%s",spec[i],spec[j]));
	priorsKs[i][j]->SetName(Form("oldpriorsKs%s%s",spec[i],spec[j]));
	priorsPhi[i][j] =  (TH2D *) fin->Get(Form("priorsPhi%s%s",spec[i],spec[j]));
	priorsPhi[i][j]->SetName(Form("oldpriorsPhi%s%s",spec[i],spec[j]));
	for(Int_t k=0; k< 3;k++){
	  priorsLc[i][j][k] =  (TH3D *) fin->Get(Form("priorsLc%s%s%s",spec[i],spec[j],spec[k]));
	  priorsLc[i][j][k]->SetName(Form("oldpriorsLc%s%s%s",spec[i],spec[j],spec[k]));
	  priorsLcbar[i][j][k] =  (TH3D *) fin->Get(Form("priorsLcbar%s%s%s",spec[i],spec[j],spec[k]));
	  priorsLcbar[i][j][k]->SetName(Form("oldpriorsLcbar%s%s%s",spec[i],spec[j],spec[k]));
	}
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
      newpriorsKs[i][j] =  new TH3D(Form("priorsKs%s%s",spec[i],spec[j]),Form("K^{0*} priors for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),200,0.4,1.4,40,0,10,nbinpol,-1.001,1.001);
      newpriorsPhi[i][j] =  new TH2D(Form("priorsPhi%s%s",spec[i],spec[j]),Form("#phi priors for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);
      for(Int_t k=0; k< 3;k++){
	newpriorsLc[i][j][k] =  new TH3D(Form("priorsLc%s%s%s",spec[i],spec[j],spec[k]),Form("#Lambda_{c}^{+} priors for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
	newpriorsLcbar[i][j][k] =  new TH3D(Form("priorsLcbar%s%s%s",spec[i],spec[j],spec[k]),Form("#overline{#Lambda}_{c}^{-} priors for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
      }
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
      truePidKs[i][j] =  new TH3D(Form("truePidKs%s%s",spec[i],spec[j]),Form("K^{0*} truePid for %s-%s;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),200,0.4,1.4,40,0,10,nbinpol,-1.001,1.001);
      truePidPhi[i][j] =  new TH2D(Form("truePidPhi%s%s",spec[i],spec[j]),Form("#phi truePid for %s-%s;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j]),100,0.98,1.05,40,0,10);

      for(Int_t k=0; k< 3;k++){
	truePidLc[i][j][k] =  new TH3D(Form("truePidLc%s%s%s",spec[i],spec[j],spec[k]),Form("#Lambda_{c}^{+} truePid for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
	truePidLcbar[i][j][k] =  new TH3D(Form("truePidLcbar%s%s%s",spec[i],spec[j],spec[k]),Form("#overline{#Lambda}_{c}^{-} truePid for %s-%s-%s;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N",spec[i],spec[j],spec[k]),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
      }

    }
  }

  trueKs =  new TH3D(Form("trueKs"),Form("K^{0*} true;m_{#piK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),200,0.4,1.4,40,0,10,nbinpol,-1.001,1.001);
  truePhi =  new TH2D(Form("truePhi"),Form("#phi true;m_{KK} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),100,0.98,1.05,40,0,10);

  trueLc =  new TH3D(Form("trueLc"),Form("#Lambda_{c}^{+} true;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
  trueLcbar =  new TH3D(Form("trueLcbar"),Form("#overline{#Lambda}_{c}^{-} true;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1);

  mypidLc =  new TH3D(Form("mypidLc"),Form("#Lambda_{c}^{+} true;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);
  mypidLcbar =  new TH3D(Form("mypidLcbar"),Form("#Lambda_{c}^{+} true;m_{#piKp} (GeV/#it{c}^2);p_{T} (GeV/#it{c});N"),nbinmlc,2.1,2.5,nbinptlc,7,27,nbinpol,0,1.);

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
  particle::AddParticleType("Lambdac+",2.28646,1,0.008); // 11
  particle::AddParticleType("Lambdacbar-",2.28646,-1,0.008); // 12
  particle::AddParticleType("Lambda1520",1.5195,-1,0.0000156); // 12

  particle d1("pi+");
  particle d2("K+");
  particle d3("p+");
  particle d4("pi-");
  particle d5("K-");
  particle d6("p-");

  particle prong1;
  particle prong2;
  particle prong3;

  particle polarKs("K0*");
  particle polarLc("Lambdac+");

  Int_t charge[] = {1,-1,1,-1,1,-1,0,0,0,2,-2,1,-1};

  particle ppos[20000];
  particle pneg[20000];

  Float_t weightsPos[20000][3];
  Float_t weightsNeg[20000][3];
  Int_t passMyPIDPos[20000][3];
  Int_t passMyPIDNeg[20000][3];

  Int_t npos=0;
  Int_t nneg=0;

  Float_t signal,signalTOF,signalTPC,pt,pz,phi,ptComb,ptComb3prong,invmass;
  Float_t ptd,pzd,phid;

  Float_t priors[3],prob[3];
  Float_t priors2[3][3],prob2[3][3];
  Float_t priors3[3][3][3],prob3[3][3][3];
      
  Int_t iev=-1,id,mother;
  Int_t cev;
  
  TFile *fout = new TFile(Form("step%i.root",step+1),"RECREATE");
  // TTree *treeKs = new TTree("treeKs","treeKs");
  // Float_t ptPair,massPair,ptD1,ptD2,weightD1[3],weightD2[3],weightFill;
  // Int_t isTrue,isTruePid;
  // treeKs->Branch("ptPair",&ptPair,"ptPair/F");
  // treeKs->Branch("massPair",&massPair,"massPair/F");
  // treeKs->Branch("ptPi",&ptD1,"ptPi/F");
  // treeKs->Branch("ptKa",&ptD2,"ptKa/F");
  // treeKs->Branch("weightPi",weightD1,"wightPi[3]/F");
  // treeKs->Branch("weightKa",weightD2,"wightKa[3]/F");
  // treeKs->Branch("weightFill",&weightFill,"wightFill/F");
  // treeKs->Branch("isTruePid",&isTruePid,"isTruePid/I");
  // treeKs->Branch("isTrue",&isTrue,"isTrue/I");
  TH1F *hcentr = new TH1F("hcentr","",100,0,100);

  FILE *flist = fopen("lista","r");
  char namefile[100];
  Float_t weight1[3],weight2[3],weight3[3];
  Float_t ptPi,ptKa,ptPr;

  TH1F *htemp;

  while(fscanf(flist,"%s",namefile)==1){
    TFile *fin = new TFile(namefile);
    printf("file = %s\n",namefile);
    TList *l = (TList *) fin->Get("TOFpid");
    htemp = (TH1F *) l->At(0);
    if(!hcentr) hcentr = new TH1F(*htemp);
    else hcentr->Add(htemp);
    TTree *t = (TTree *) l->At(1);
    Int_t n = t->GetEntries();
    
    for(Int_t i=0;i < n;i++){
      t->GetEvent(i);

      ptPi = t->GetLeaf("ptPi")->GetValue();
      ptKa = t->GetLeaf("ptPi")->GetValue();
      ptPr = t->GetLeaf("ptPi")->GetValue();
      pt = t->GetLeaf("pt")->GetValue();
      invmass = t->GetLeaf("mass")->GetValue();
      weight1[0] = t->GetLeaf("weightPi")->GetValue(0);
      weight1[1] = t->GetLeaf("weightPi")->GetValue(1);
      weight1[2] = t->GetLeaf("weightPi")->GetValue(2);
      weight2[0] = t->GetLeaf("weightKa")->GetValue(0);
      weight2[1] = t->GetLeaf("weightKa")->GetValue(1);
      weight2[2] = t->GetLeaf("weightKa")->GetValue(2);
      weight3[0] = t->GetLeaf("weightPr")->GetValue(0);
      weight3[1] = t->GetLeaf("weightPr")->GetValue(1);
      weight3[2] = t->GetLeaf("weightPr")->GetValue(2);

      Float_t pt1 = Int_t(ptPi/(ptPi+ptKa+ptPr)*nbinPtFrPi);
      Float_t pt2 = Int_t(ptKa/(ptPi+ptKa+ptPr)*nbinPtFrKa);
      Float_t polar = 0;//TMath::Abs(polarLc.GetY());//ptComb3prong/ptot;//(pt2*nbinpol + pt1)*invpollc;
      polar = ((pt1*nbinPtFrKa + pt2 + polar)*nbinY)*normbin;

      Int_t ibinx = priorsLc[0][0][0]->GetXaxis()->FindBin(invmass);
      Int_t ibiny = priorsLc[0][0][0]->GetYaxis()->FindBin(pt);
      Int_t ibinz = priorsLc[0][0][0]->GetZaxis()->FindBin(polar);
      
      for(Int_t ipr=0;ipr<3;ipr++)
	for(Int_t jpr=0;jpr<3;jpr++)
	  for(Int_t kpr=0;kpr<3;kpr++)
	    priors3[ipr][jpr][kpr] = priorsLc[ipr][jpr][kpr]->GetBinContent(ibinx,ibiny,ibinz);
      
      GetProb3(weight1,weight2,weight3,priors3,prob3);

      for(Int_t ipr=0;ipr<3;ipr++)
	for(Int_t jpr=0;jpr<3;jpr++)
	  for(Int_t kpr=0;kpr<3;kpr++){
	    newpriorsLc[ipr][jpr][kpr]->Fill(invmass,pt,polar,prob3[ipr][jpr][kpr]);
	  }
    }
    t->Delete();
    fin->Close();

  }



  printf("Write output\n");
  fout->cd();
  hcentr->Write();
  //if(step==0) treeKs->Write();
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
      for(Int_t k=0;k<3;k++){
       	newpriorsLc[i][j][k]->Write();
	truePidLc[i][j][k]->Write();
	newpriorsLcbar[i][j][k]->Write();
	truePidLcbar[i][j][k]->Write();
      }
    }
  }
  trueKs->Write();
  truePhi->Write();
  trueLc->Write();
  trueLcbar->Write();
  mypidLc->Write();
  mypidLcbar->Write();
  fout->Close();
}

// return the weight array for pion, kaon and proton hypoteses
// conditioned probability: a given particles releases the measured signal
void ComputeWeights(Float_t weights[3],Float_t signal,Float_t pt){
  Float_t expectPi = -fseparation->Eval(pt);
  Float_t expectPr = fseparation->Eval(pt);
  
  if(kALICEseparation){
    expectPi = -fseparationPiKa->Eval(pt);
    expectPr = fseparationKaPr->Eval(pt);
  }
  
  signal -= addshift;

  weights[0] = TMath::Exp(-(signal-expectPi)*(signal-expectPi)*0.5*invwidth*invwidth);
  weights[1] = TMath::Exp(-signal*signal*0.5*invwidth*invwidth);
  weights[2] = TMath::Exp(-(signal-expectPr)*(signal-expectPr)*0.5*invwidth*invwidth);
}

void ComputeWeightsALICE(Float_t weights[3],Float_t signalTPC,Float_t signalTOF,Float_t pt,Float_t p){
  Float_t expectPiTPC = fTPCpi->Eval(p);
  Float_t expectKaTPC = fTPCka->Eval(p);
  Float_t expectPrTPC = fTPCpr->Eval(p);
  
  signalTPC -= addshift;

  weights[0] = TMath::Exp(-(signalTPC-expectPiTPC)*(signalTPC-expectPiTPC)*0.5*invwidth*invwidth);
  weights[1] = TMath::Exp(-(signalTPC-expectKaTPC)*(signalTPC-expectKaTPC)*0.5*invwidth*invwidth);
  weights[2] = TMath::Exp(-(signalTPC-expectPrTPC)*(signalTPC-expectPrTPC)*0.5*invwidth*invwidth);


  Float_t propFactorPi = 1 - fEfficiencyPiTOF->Eval(pt);
  Float_t propFactorKa = 1 - fEfficiencyKaTOF->Eval(pt);
  Float_t propFactorPr = 1 - fEfficiencyPrTOF->Eval(pt);

  if(signalTOF > -999){
    Float_t expectPiTOF = fTOFpi->Eval(p,pt/p);
    Float_t expectKaTOF = fTOFka->Eval(p,pt/p);
    Float_t expectPrTOF = fTOFpr->Eval(p,pt/p);
    
    signalTOF -= addshiftTOF;

    // apply the TOF propagation factor here
    weights[0] *= TMath::Exp(-(signalTOF-expectPiTOF)*(signalTOF-expectPiTOF)*0.5*invwidthTOF*invwidthTOF);
    weights[1] *= TMath::Exp(-(signalTOF-expectKaTOF)*(signalTOF-expectKaTOF)*0.5*invwidthTOF*invwidthTOF);
    weights[2] *= TMath::Exp(-(signalTOF-expectPrTOF)*(signalTOF-expectPrTOF)*0.5*invwidthTOF*invwidthTOF);

    propFactorPi = 1 - propFactorPi;
    propFactorKa = 1 - propFactorKa;
    propFactorPr = 1 - propFactorPr;
  }

  // correct for propagation factor
  weights[0] *= propFactorPi;
  weights[1] *= propFactorKa;
  weights[2] *= propFactorPr;

  // 5 sigma truncation
  weights[0] += 1E-6;
  weights[1] += 1E-6;
  weights[2] += 1E-6;

}


// return the probability array for pion, kaon and proton hypoteses
void GetProb1(Float_t weights[3],Float_t priors[3],Float_t prob[3]){
  prob[0] = weights[0]*(priors[0]+0.000001);
  prob[1] = weights[1]*(priors[1]+0.000001);
  prob[2] = weights[2]*(priors[2]+0.000001);
    Float_t R= prob[0]+prob[1]+prob[2];
    R = 1./R;
    prob[0] *= R;
    prob[1] *= R;
    prob[2] *= R;
}

// return the probability array for all the pairs 3x3: (pi,pi), (pi,ka), ... , (pr,ka), (pr,pr)
void GetProb2(Float_t weights1[3],Float_t weights2[3],Float_t priors[3][3],Float_t prob[3][3]){
  Float_t R= 0;
  for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
      prob[i][j] = weights1[i]*weights2[j]*(priors[i][j]+0.000001);
    R += prob[i][j];
    }
  }
  R = 1./R;
   for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
      prob[i][j] *= R;
    }
   }
}

// return the probability array for all the pairs 3x3: (pi,pi), (pi,ka), ... , (pr,ka), (pr,pr)
void GetProb3(Float_t weights1[3],Float_t weights2[3],Float_t weights3[3],Float_t priors[3][3][3],Float_t prob[3][3][3]){
  Float_t R= 0;
  for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
      for(Int_t k=0;k < 3;k++){
	prob[i][j][k] = weights1[i]*weights2[j]*weights3[k]*(priors[i][j][k]+0.000001);
	R += prob[i][j][k];
      }
    }
  }

  R = 1./R;
   for(Int_t i=0;i < 3;i++){
    for(Int_t j=0;j < 3;j++){
      for(Int_t k=0;k < 3;k++){
	prob[i][j][k] *= R;
      }
    }
   }
}

Float_t GetPolariz(particle moth,particle dau){
  double bx = -moth.GetPx()/moth.GetEnergy();
  double by = -moth.GetPy()/moth.GetEnergy();
  double bz = -moth.GetPz()/moth.GetEnergy();

  Float_t momMoth = TMath::Sqrt(moth.GetPx()*moth.GetPx() + moth.GetPy()*moth.GetPy() + moth.GetPz()*moth.GetPz());

  if(momMoth == 0) return 0;

  dau.Boost(bx,by,bz);

  Float_t momDau = TMath::Sqrt(dau.GetPx()*dau.GetPx() + dau.GetPy()*dau.GetPy() + dau.GetPz()*dau.GetPz());

  if(momDau == 0) return 0;

  Float_t polar = dau.GetPx()*moth.GetPx() + dau.GetPy()*moth.GetPy() + dau.GetPz()*moth.GetPz();
  polar /= momMoth*momDau;

  return polar;
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

Int_t lambdacGood(Float_t ptLc,Float_t ptPi,Float_t ptKa,Float_t ptPr){

  //  printf("%f %f %f\n",ptPi,ptKa,ptPr);
//   if(ptPi < ptminPi) return 0;
//   if(ptKa < ptminKa) return 0;
//   if(ptPr < ptminPr) return 0;

  if(ptPi < ptLc*0.06 + 0.15) return 0;
  if(ptKa < ptLc*0.12 + 0.18) return 0;
  if(ptPr < ptLc*0.3) return 0;

  return 1;
}

Int_t passMyPID(Int_t isp,Float_t p,Float_t pt,Float_t TPCsign,Float_t TOFsign){
  // isp=0(pi), 1(ka), 2(pr)
  Float_t expectTPC[3];
  expectTPC[0]= fTPCpi->Eval(p);
  expectTPC[1]= fTPCka->Eval(p);
  expectTPC[2]= fTPCpr->Eval(p);
  
  if(TOFsign > -999){
    Float_t expectTOF[3];
    expectTOF[0] = fTOFpi->Eval(p,pt/p);
    expectTOF[1] = fTOFka->Eval(p,pt/p);
    expectTOF[2] = fTOFpr->Eval(p,pt/p);
  }

  // for example TPC Nsigma variable for pions is: TPCsign - expectTPC[0]
  // for example TOF Nsigma variable for pions is: TOFsign - expectTOF[0] (please check the TOF is available: TOFsign > -999)


  return 1;
}

Int_t passMyPID2(Int_t isp,Float_t p,Float_t pt,Float_t TPCsign,Float_t TOFsign){
  Float_t weights[3],prob[3];
  ComputeWeightsALICE(weights,TPCsign,TOFsign,pt,p);

  // define your priors
  Float_t priors[3] = {5,3,2};

  GetProb1(weights,priors,prob);

  // for example to require a probability > 50%
  if(prob[isp] < 0.3) return 0;

  return 1;
}
