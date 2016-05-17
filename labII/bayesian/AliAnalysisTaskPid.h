#ifndef ALIANALYSISTASKPID_H
#define ALIANALYSISTASKPID_H

// ROOT includes
#include <TObject.h>
#include <TClonesArray.h>
#include <TList.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>

// output object
#include <TTree.h>
#include "AliVEvent.h"
#include "TH1F.h"

class AliPIDResponse;
class AliPIDCombined;
class AliFlowBayesianPID;

class AliAnalysisTaskPid : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPid();
  AliAnalysisTaskPid(const char *name);

  virtual ~AliAnalysisTaskPid();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 
  
  virtual void SetMC(Bool_t flag = kTRUE){fIsMC = flag;};

 private:
  AliAnalysisTaskPid(const AliAnalysisTaskPid &old); 
  AliAnalysisTaskPid& operator=(const AliAnalysisTaskPid &source); 

  Bool_t fIsMC;                 // if MC
  TList *fList;                 //! List for output objects
  AliPIDResponse *fPIDResponse; //! PID response object
  AliPIDCombined *fPIDCombined; //! PID combined object

  Bool_t fIsFirstEvent;         //! to execute some routines only for first event
  Float_t fCentrality;          //! centrality if defined
  Float_t fVertexZ;             //! primary vertex (z) in cm

  TH1F *fCentrTOF;              //! centrality distribtion of TOF tracks
  TTree *fTreePr;               //! tree for priors
  Float_t fPt;                  //! pt track
  Float_t fPz;                  //! pz track
  Float_t fPhi;                 //! phi track
  Float_t fWeightTOF[9];        //! array weights
  Float_t fWeightTPC[9];        //! array weights
  Float_t fWeight1[3];          //! array weights pi
  Float_t fWeight2[3];          //! array weights ka
  Float_t fWeight3[3];          //! array weights pr
  Float_t fMass;                //! mass lambdac
  Float_t fPtPi;                //! pt pion
  Float_t fPtKa;                //! pt kaon
  Float_t fPtPr;                //! pt proton
  Float_t fMismWeight;          //! mismatch weight
  Int_t fTOF,fTPC;              //! pid flags
  Int_t fNtofCl;                //! number of TOF clusters
  Float_t fPiTOFsigma;          //! tof n sigma
  AliFlowBayesianPID *fBayesPID;//! Bayesian pid object flow package

  ClassDef(AliAnalysisTaskPid, 1);    //Analysis task v2 and v3 analysis on AOD
};

#endif
