#include "runAction.hh"

RunAction::RunAction(){
	// Create Ntuples when program is activated, not for every run
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->CreateNtuple("Hits", "Hits");
	man->CreateNtupleIColumn("fEventID"); // 0
	man->CreateNtupleIColumn("fParticleTrackID"); 
	man->CreateNtupleIColumn("fParentTrackID"); // 2
	man->CreateNtupleIColumn("fPDGCode");
	man->CreateNtupleIColumn("fDetectorCopyNumber"); // 4
	man->CreateNtupleDColumn("fXDetector"); 
	man->CreateNtupleDColumn("fYDetector"); // 6
	man->CreateNtupleDColumn("fZDetector");
	man->CreateNtupleDColumn("fXParticle"); // 8
	man->CreateNtupleDColumn("fYParticle"); 
	man->CreateNtupleDColumn("fZParticle"); // 10
	man->CreateNtupleDColumn("fPx");
	man->CreateNtupleDColumn("fPy"); // 12
	man->CreateNtupleDColumn("fPz");
	man->CreateNtupleDColumn("fEnergy"); // 14
	man->FinishNtuple(0); // NTuple number 0. If you want additional NTuples, you have to increment the number

}
RunAction::~RunAction(){}

void RunAction::BeginOfRunAction(const G4Run* run){
	// Create a file every run
	G4AnalysisManager* man = G4AnalysisManager::Instance();

	G4int runID = run->GetRunID();

	std::stringstream strRunID;
	strRunID << runID;
	man->OpenFile("/home/cas/phd/FCT_toymodel/GEANT4FCT_v2/build/output/O2event_"+strRunID.str()+".root");

}

void RunAction::EndOfRunAction(const G4Run*){
	G4AnalysisManager* man = G4AnalysisManager::Instance();
	man->Write();
	man->CloseFile();
}