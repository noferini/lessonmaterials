#include "AliAnalysisTaskPid.h"

// ROOT includes
#include <TMath.h>

// AliRoot includes
#include "AliInputEventHandler.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliCentrality.h"
#include "AliVHeader.h"
#include "AliMCEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TChain.h"

#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTOFHeader.h"
#include "TRandom.h"
#include "AliFlowBayesianPID.h"

ClassImp(AliAnalysisTaskPid)

//_____________________________________________________________________________
AliAnalysisTaskPid::AliAnalysisTaskPid():
  AliAnalysisTaskSE(),
  fIsMC(kFALSE),
  fList(new TList()),
  fPIDResponse(NULL),
  fPIDCombined(NULL),
  fIsFirstEvent(kFALSE),
  fCentrality(0.),
  fVertexZ(0.),
  fCentrTOF(NULL),
  fTreePr(NULL),
  fPt(0),
  fMismWeight(0.0),
  fTOF(0),
  fTPC(0),
  fNtofCl(0),
  fPiTOFsigma(100),
  fBayesPID(NULL)
{
  // Default constructor (should not be used)
  fList->SetName("TOFpid");

  fList->SetOwner(kTRUE); 
}

//______________________________________________________________________________
AliAnalysisTaskPid::AliAnalysisTaskPid(const char *name):
  AliAnalysisTaskSE(name),
  fIsMC(kFALSE),
  fList(new TList()),
  fPIDResponse(NULL),
  fPIDCombined(NULL),
  fIsFirstEvent(kFALSE),
  fCentrality(0.),
  fVertexZ(0.),
  fCentrTOF(NULL),
  fTreePr(NULL),
  fPt(0),
  fMismWeight(0.0),
  fTOF(0),
  fTPC(0),
  fNtofCl(0),
  fPiTOFsigma(100),
  fBayesPID(NULL)
{

  DefineOutput(1, TList::Class());

  // Output slot #1 writes into a TTree
  fList->SetName("TOFpid");

  fList->SetOwner(kTRUE); 
}
//_____________________________________________________________________________
AliAnalysisTaskPid::~AliAnalysisTaskPid()
{

}

//______________________________________________________________________________
void AliAnalysisTaskPid::UserCreateOutputObjects()
{
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC+AliPIDResponse::kDetTOF);

  fCentrTOF = new TH1F("hCentrTOF","p_{T} disribution for TOF tracks;centrality;N",100,0,100);
  fList->Add(fCentrTOF);
  fTreePr = new TTree("treePr","treePr");
  fTreePr->Branch("centr",&fCentrality,"centr/F");
  fTreePr->Branch("pt",&fPt,"pt/F");
  fTreePr->Branch("pz",&fPz,"pz/F"); // longitudinal momentum
  fTreePr->Branch("phi",&fPhi,"phi/F"); // azhimuthal angle
  fTreePr->Branch("weightPi",fWeight1,"weightPi[3]/F");
  fTreePr->Branch("weightKa",fWeight2,"weightKa[3]/F");
  fTreePr->Branch("weightPr",fWeight3,"weightPr[3]/F");
  fTreePr->Branch("mass",&fMass,"mass/F");
  fTreePr->Branch("ptPi",&fPtPi,"ptPi/F");
  fTreePr->Branch("ptKa",&fPtKa,"ptKa/F");
  fTreePr->Branch("ptPr",&fPtPr,"ptPr/F");
//   fTreePr->Branch("nTOFcl",&fNtofCl,"nTOFcl/I");
//   fTreePr->Branch("nTOFpiSigma",&fPiTOFsigma,"nTOFpiSigma/F");

  fList->Add(fTreePr);

  fBayesPID = new AliFlowBayesianPID();
  fBayesPID->SetNewTrackParam();
  fBayesPID->ForceOldDedx();

  // Post output data.
  PostData(1, fList);

}

//______________________________________________________________________________
void AliAnalysisTaskPid::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
  
  // reset variables;
  fVertexZ = 999.;

  Float_t weightsPos[20000][3];
  Float_t weightsNeg[20000][3];
  Float_t pxPos[20000],pyPos[20000],pzPos[20000];
  Float_t pxNeg[20000],pyNeg[20000],pzNeg[20000];

  Int_t npos=0;
  Int_t nneg=0;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");


  fOutputAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if((!fOutputAOD)){
    Printf("%s:%d AODEvent/ESDEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }
 
  const AliTOFHeader *tofHeader = fOutputAOD->GetTOFHeader();
  
//   Double_t probTPCTOF[AliPID::kSPECIESC];
  
//   Int_t run = fOutputAOD->GetRunNumber();
    
  const AliAODVertex* vtxAOD = fOutputAOD->GetPrimaryVertex();
  if(vtxAOD->GetNContributors()>0)
    fVertexZ = vtxAOD->GetZ();

  AliAODMCHeader *mcHeader = NULL;
  TClonesArray *mcArray = NULL;
  
  //Get the MC object
  if(fIsMC){
    mcHeader = dynamic_cast<AliAODMCHeader*>(fOutputAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      AliError("Could not find MC Header in AOD");
      return;
    }
    if(mcHeader)
      mcArray = (TClonesArray*)fOutputAOD->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  }
  
  AliCentrality *centrality = fOutputAOD->GetCentrality();
  Float_t v0Centr  = 100.;
  Float_t trkCentr  = 100.;
  if(centrality){
    trkCentr  = centrality->GetCentralityPercentile("TRK");
    v0Centr = centrality->GetCentralityPercentile("V0M"); 
  }
  else printf("no centrality available\n\n");
  
  fCentrality =v0Centr;
  if(centrality){
    if(fCentrality < 0 || fCentrality > 99.999) fCentrality = 99.999;
  }

  if(!fIsFirstEvent){    
    fIsFirstEvent = kTRUE;
  }

  fBayesPID->SetDetResponse(fOutputAOD, fCentrality,AliESDpid::kTOF_T0);

  Int_t nTPConlytrack=0;
   
  if(tofHeader) fNtofCl=tofHeader->GetNumberOfTOFclusters(); 

  fCentrTOF->Fill(fCentrality);

  Int_t ntracks=fOutputAOD->GetNumberOfTracks();
  for (Int_t itrack=0; itrack<ntracks; ++itrack){
    AliAODTrack *track = (AliAODTrack *) fOutputAOD->GetTrack(itrack);

    Bool_t isTPConly = kFALSE;
    if(track->TestFilterBit(1)) isTPConly = kTRUE;

    Bool_t isGlobal = kFALSE;
    if(track->TestFilterBit(32)) isGlobal = kTRUE;

    if(isTPConly && TMath::Abs(track->Eta()) < 0.8 && track->Pt() > 0.2)
      nTPConlytrack++;

    if(!isTPConly || TMath::Abs(track->Eta()) > 0.8) continue;
    

    fTPC = track->GetTPCsignal() > 10;

    fTOF = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME) && (track->GetStatus() & AliVTrack::kITSrefit);
    
    if(!(fTOF || fTPC)) continue;

    fPt = track->Pt();
    fPz = track->Pz();
    fPhi = track->Phi();

    Double_t tpcW[9] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
    Double_t tofW[9] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};

    fPiTOFsigma = 100;

    fMismWeight = 0;

    if(fTPC){
      // official framework
      fPIDResponse->ComputeTPCProbability(track,9,tpcW);

      // flow framework to use my own TPC spline
      //fBayesPID->ComputeProb(track);
      //const Float_t *tpcWeight = fBayesPID->GetWeights(0);
      //for(Int_t ii=0;ii<9;ii++) tpcW[ii] = tpcWeight[ii];
    }
    if(fTOF){
      Float_t oldMism;
      fPIDResponse->ComputeTOFProbability(track,9,tofW);
      oldMism = fPIDResponse->GetTOFMismatchProbability();

      //      printf("mism = %f -- pi = %f -- centr = %f\n",fMismWeight,tofW[2],fCentrality);

      fMismWeight = fPIDResponse->GetTOFResponse().GetMismatchProbability(track->GetTOFsignal(),track->Eta());

      fPiTOFsigma = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kPion);

      for(Int_t is=0;is < 9;is++) tofW[is] -= oldMism;
    }


    for(Int_t is=0;is < 9;is++) fWeightTPC[is] = tpcW[is],fWeightTOF[is] = tofW[is];

    //    UInt_t detUsed = fPIDCombined->ComputeProbabilities(track, fPIDResponse, probTPCTOF);

    Float_t ranv = 0;
    if(fPt < 2 && fCentrality < 99 && (tofW[5]+tofW[6]+tofW[7]+tofW[8] < 1E-8)) ranv = gRandom->Rndm();

    if(fPt > 0.6){
      if(track->Charge()>0){
	weightsPos[npos][0] = tpcW[2]*(tofW[2]+fMismWeight*0.1);
	weightsPos[npos][1] = tpcW[3]*(tofW[3]+fMismWeight*0.1);
	weightsPos[npos][2] = tpcW[4]*(tofW[4]+fMismWeight*0.1);

	pxPos[npos] = fPt*TMath::Cos(fPhi);
	pyPos[npos] = fPt*TMath::Sin(fPhi);
	pzPos[npos] = fPz;
	npos++;
      }
      else{
	weightsNeg[nneg][0] = tpcW[2]*(tofW[2]+fMismWeight*0.1);
	weightsNeg[nneg][1] = tpcW[3]*(tofW[3]+fMismWeight*0.1);
	weightsNeg[nneg][2] = tpcW[4]*(tofW[4]+fMismWeight*0.1);

	pxNeg[nneg] = fPt*TMath::Cos(fPhi);
	pyNeg[nneg] = fPt*TMath::Sin(fPhi);
	pzNeg[nneg] = fPz;
	nneg++;
      }


    }
    

  }
  

  Float_t ptPi,ptKa,ptPr,pLambdac[3],ptLambdac,momLambdac,eLambdac,mLambdac;


  for(Int_t i=0;i<npos;i++){
    fPtPr = TMath::Sqrt(pxPos[i]*pxPos[i] + pyPos[i]*pyPos[i]);
    if(fPtPr < 2.1) continue;
    for(Int_t j=0;j<nneg;j++){
      fPtKa = TMath::Sqrt(pxNeg[j]*pxNeg[j] + pyNeg[j]*pyNeg[j]);
      if(fPtKa < 1) continue;
      for(Int_t k=0;k<npos;k++){
	if(k==i) continue;
	fPtPi = TMath::Sqrt(pxPos[k]*pxPos[k] + pyPos[k]*pyPos[k]);

	pLambdac[0] = pxPos[i] + pxNeg[j] + pxPos[k];
	pLambdac[1] = pyPos[i] + pyNeg[j] + pyPos[k];
	pLambdac[2] = pzPos[i] + pzNeg[j] + pzPos[k];
	ptLambdac = TMath::Sqrt(pLambdac[0]*pLambdac[0] + pLambdac[1]*pLambdac[1]); 
	if(ptLambdac < 7 || ptLambdac > 20 || 0.06*ptLambdac+0.15>fPtPi || 0.12*ptLambdac+0.18>fPtKa || 0.3*ptLambdac>fPtPr) continue;
	momLambdac = TMath::Sqrt(ptLambdac*ptLambdac + pLambdac[2]*pLambdac[2]); 
	eLambdac = TMath::Sqrt(pxPos[i]*pxPos[i]+pyPos[i]*pyPos[i]+pzPos[i]*pzPos[i]+0.938*0.938);
	eLambdac += TMath::Sqrt(pxNeg[j]*pxNeg[j]+pyNeg[j]*pyNeg[j]+pzNeg[j]*pzNeg[j]+0.493*0.493);
	eLambdac += TMath::Sqrt(pxPos[k]*pxPos[k]+pyPos[k]*pyPos[k]+pzPos[k]*pzPos[k]+0.139*0.139);

	mLambdac = TMath::Sqrt(eLambdac*eLambdac - momLambdac*momLambdac);

	if(mLambdac > 2.5|| mLambdac < 2.1) continue;

	fPt = ptLambdac;
	fPhi = TMath::ATan2(pLambdac[1],pLambdac[0]);
	fPz = pLambdac[2];
	fWeight1[0] = weightsPos[k][0];
	fWeight1[1] = weightsPos[k][1];
	fWeight1[2] = weightsPos[k][2];
	fWeight2[0] = weightsNeg[j][0];
	fWeight2[1] = weightsNeg[j][1];
	fWeight2[2] = weightsNeg[j][2];
	fWeight3[0] = weightsPos[i][0];
	fWeight3[1] = weightsPos[i][1];
	fWeight3[2] = weightsPos[i][2];	
	fMass = mLambdac;
	fTreePr->Fill();


      }
    }
  }

  for(Int_t i=0;i<nneg;i++){
    fPtPr = TMath::Sqrt(pxNeg[i]*pxNeg[i] + pyNeg[i]*pyNeg[i]);
    if(fPtPr < 2.1) continue;
    for(Int_t j=0;j<npos;j++){
      fPtKa = TMath::Sqrt(pxPos[j]*pxPos[j] + pyPos[j]*pyPos[j]);
      if(fPtKa < 1) continue;
      for(Int_t k=0;k<nneg;k++){
	if(k==i) continue;
	fPtPi = TMath::Sqrt(pxNeg[k]*pxNeg[k] + pyNeg[k]*pyNeg[k]);

	pLambdac[0] = pxNeg[i] + pxPos[j] + pxNeg[k];
	pLambdac[1] = pyNeg[i] + pyPos[j] + pyNeg[k];
	pLambdac[2] = pzNeg[i] + pzPos[j] + pzNeg[k];
	ptLambdac = TMath::Sqrt(pLambdac[0]*pLambdac[0] + pLambdac[1]*pLambdac[1]); 
	if(ptLambdac < 7 || ptLambdac > 20 || 0.06*ptLambdac+0.15>fPtPi || 0.12*ptLambdac+0.18>fPtKa || 0.3*ptLambdac>fPtPr) continue;
	momLambdac = TMath::Sqrt(ptLambdac*ptLambdac + pLambdac[2]*pLambdac[2]); 
	eLambdac = TMath::Sqrt(pxPos[j]*pxPos[j]+pyPos[j]*pyPos[j]+pzPos[j]*pzPos[j]+0.493*0.493);
	eLambdac += TMath::Sqrt(pxNeg[i]*pxNeg[i]+pyNeg[i]*pyNeg[i]+pzNeg[i]*pzNeg[i]+0.938*0.933);
	eLambdac += TMath::Sqrt(pxNeg[k]*pxNeg[k]+pyNeg[k]*pyNeg[k]+pzNeg[k]*pzNeg[k]+0.193*0.193);

	mLambdac = TMath::Sqrt(eLambdac*eLambdac - momLambdac*momLambdac);

	if(mLambdac > 2.5|| mLambdac < 2.1) continue;

	fPt = ptLambdac;
	fPhi = TMath::ATan2(pLambdac[1],pLambdac[0]);
	fPz = pLambdac[2];
	fWeight1[0] = weightsNeg[k][0];
	fWeight1[1] = weightsNeg[k][1];
	fWeight1[2] = weightsNeg[k][2];
	fWeight2[0] = weightsPos[j][0];
	fWeight2[1] = weightsPos[j][1];
	fWeight2[2] = weightsPos[j][2];
	fWeight3[0] = weightsNeg[i][0];
	fWeight3[1] = weightsNeg[i][1];
	fWeight3[2] = weightsNeg[i][2];	
	fMass = -mLambdac;

	fTreePr->Fill();

      }
    }
  }





  // Fill once per event
}

//_____________________________________________________________________________
void AliAnalysisTaskPid::Terminate(Option_t *)
{ 
  // Terminate loop
  Printf("Terminate()");
}
//_____________________________________________________________________________
