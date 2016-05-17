void makeAna
(
 const char *kCollectionFile="wn.xml",            // XML file containing tags
 Long64_t    nentries=TChain::kBigNumber
 )
{
  // Connecting to the PROOF cluster
  Bool_t grid = 1; // run on grid?
  Bool_t connect = kFALSE; // need the grid?
  if(grid) connect = kTRUE;

  Bool_t isMCdata=kFALSE; // if the data are MC
  Bool_t isMC=kFALSE; // if you have the kinematics information
  if(isMC) isMCdata=kTRUE;

  Bool_t tuneOnData=kTRUE;
  if(!isMC) tuneOnData=kFALSE;

  // include the path you need to compile and the library you need
  gSystem->SetIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");   

  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libPWGPPpid.so");
  gSystem->Load("libPWGflowBase.so");
  gSystem->Load("libPWGflowTasks.so");

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");

  gROOT->LoadMacro("AliAnalysisTaskPid.cxx++");
  gROOT->LoadMacro("AddTaskPid.C");

  // Connect to the GRID
  if(connect){
    TGrid::Connect("alien://");
  }

  // Create the analysis manager
  mgr = new AliAnalysisManager("testAnalysis");

  // Add AOD handler
  AliAODInputHandler* aodH = new AliAODInputHandler;
  mgr->SetInputEventHandler(aodH);
  // for MC
  if(isMC){
 //    AliMCEventHandler* mcH = new AliMCEventHandler();
//     mgr->SetMCtruthEventHandler(mcH);
  }

  // create chain of files to read
  TChain* analysisChain = new TChain("aodTree");

  TAlienCollection *myCollection = TAlienCollection::Open(kCollectionFile);
  if (!myCollection) {
    Error("AliRsnReadTaskAlien", Form("Cannot create an AliEn collection from %s", kCollectionFile));
    return;
  }
  myCollection->Reset();    
  
  // loop on the entries of the XML input file
  Int_t ifile = 0;
  Int_t nmaxfile =30; // to limit the number of ESD files
    
  // read the list from a xml collection "wn.xml" (grid case)
  while (grid && myCollection->Next() && ifile < nmaxfile) {
    ifile++;
    char esdFile[255];
    sprintf(esdFile, "%s", myCollection->GetTURL(""));
    Info("AliTaskAlienForTOF", Form("Adding %s", esdFile));
    analysisChain->Add(esdFile);
  }
  // read the list from a local text file "lista" (local case)
  if(! grid){
    char nomefile[100];
    FILE *fl = fopen("lista","r");
    while(fscanf(fl,"%s",nomefile) == 1 && ifile < nmaxfile){
      analysisChain->Add(nomefile);
      ifile++;
    }
    fclose(fl);
  }
  Info("Task", Form("CHAIN HAS %d ENTRIES", (Int_t)analysisChain->GetEntries()));

  // Quality tasks
  // compile my own code


//  AddTaskCentrality(1,1);
  AddTaskPIDResponse(isMC,kTRUE,tuneOnData,2,0);
  AddTaskPid(isMC);

  printf("Tasks added\n");

  // Enable debug printouts
  // mgr->SetDebugLevel(2);

  // Run analysis
  //	AliLog::SetGlobalLogLevel(AliLog::kError);
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("testAnalysis", analysisChain);
  
}

