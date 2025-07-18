#include "sensitiveDetector.hh"

// In FCT known as TrackerSD

SensitiveDetector::SensitiveDetector(G4String name): G4VSensitiveDetector(name){}
SensitiveDetector::~SensitiveDetector(){}

G4bool SensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist){
	G4Track* track = aStep->GetTrack();
	G4int trackID = track->GetTrackID();

	G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	G4StepPoint* preStepPoint = aStep->GetPreStepPoint(); // When the particle enters the detector
	G4StepPoint* postStepPoint = aStep->GetPostStepPoint(); // When the particle leaves the detector
	// For charged particles this is more difficult 
	// becasue every step of the trajectory has a pre and post step point, not only where it entered and exited the volume

	G4ThreeVector posParticle = preStepPoint->GetPosition();


	// Get the copy number of the volume the photon entered
	const G4VTouchable *touchable = aStep->GetPreStepPoint()->GetTouchable();
	G4int copyNo = touchable->GetCopyNumber();

	// Get position of the physical volume the photon entered
	G4VPhysicalVolume *physVol = touchable->GetVolume();
	G4ThreeVector posDetector = physVol->GetTranslation();

	G4ThreeVector momParticle = track->GetMomentum();
	// if(preStepPoint->GetTotalEnergy() > 0.1){
	// 	G4cout << "------- Particle detected. Event: " << eventID << " -------" << G4endl;
	// 	G4cout << "TrackiD: " << trackID << G4endl;
	// 	G4cout << "ParentID: " << track->GetParentID() << G4endl;
	// 	G4cout << "PDGCode: " << track->GetParticleDefinition()->GetPDGEncoding() << G4endl;
	// 	G4cout << "Particle position: " << posParticle << G4endl;
	// 	G4cout << "Copy number: " << copyNo << G4endl;
	// 	G4cout << "Particle momentum: " << momParticle << G4endl;
	// 	G4cout << "Detector position: " << posDetector << G4endl;
	// 	G4cout << "E: " << preStepPoint->GetTotalEnergy() << G4endl;
	// 	G4cout << "---------------------------------" << G4endl;
	// }

	// Now put the data into the Ntuples
	G4AnalysisManager* man = G4AnalysisManager::Instance();

	man->FillNtupleIColumn(0, eventID);
	man->FillNtupleIColumn(1, trackID);
	man->FillNtupleIColumn(2, track->GetParentID());
	man->FillNtupleIColumn(3, track->GetParticleDefinition()->GetPDGEncoding());
	man->FillNtupleIColumn(4, copyNo);
	man->FillNtupleDColumn(5, posDetector[0]);
	man->FillNtupleDColumn(6, posDetector[1]);
	man->FillNtupleDColumn(7, posDetector[2]);
	man->FillNtupleDColumn(8, posParticle[0]);
	man->FillNtupleDColumn(9, posParticle[1]);
	man->FillNtupleDColumn(10, posParticle[2]);
	man->FillNtupleDColumn(11, momParticle[0]);
	man->FillNtupleDColumn(12, momParticle[1]);
	man->FillNtupleDColumn(13, momParticle[2]);
	man->FillNtupleDColumn(14, preStepPoint->GetTotalEnergy());
	man->AddNtupleRow(0); // Select which Ntuple. 0'th in this case		
	
}