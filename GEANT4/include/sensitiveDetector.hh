#ifndef SENSITIVEDETECTOR_HH
#define SENSITIVEDETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"

class SensitiveDetector : public G4VSensitiveDetector{
	public:
		SensitiveDetector(G4String);
		~SensitiveDetector();
	private:
		G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
};



#endif