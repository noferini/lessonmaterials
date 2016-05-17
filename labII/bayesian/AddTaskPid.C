AliAnalysisTask *AddTaskPid(Bool_t ismc=kFALSE){

  //get the current analysis manager
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("No manager found in AddTaskPid. Why?");
    return 0;
  }

  //========= Add tender to the ANALYSIS manager and set default storage =====
  char mytaskName[100];
  snprintf(mytaskName,100,"AliAnalysisTaskPid.cxx"); 

  AliAnalysisTaskPid *task = new AliAnalysisTaskPid(mytaskName);
  if(ismc) task->SetMC(ismc);

  mgr->AddTask(task);

  //Attach input to my tasks
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(),AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());

  // Attach output to my tasks
  AliAnalysisDataContainer *cOutputL= mgr->CreateContainer("TOFpid",TList::Class(), AliAnalysisManager::kOutputContainer, AliAnalysisManager::GetCommonFileName());
  mgr->ConnectOutput(task, 1, cOutputL);

  printf("task really added\n");

  return task;
}

