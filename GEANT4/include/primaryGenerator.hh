#ifndef PRIMARYGENERATOR_HH
#define PRIMARYGENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4GenericMessenger.hh"
#include "G4String.hh"
#include "Randomize.hh"
#include <math.h>

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
	public:
		PrimaryGenerator(std::vector<G4double> firstLayer, G4String particleFile);
		~PrimaryGenerator();

		void GeneratePrimaries(G4Event* anEvent) override;
	
	private:
		G4GenericMessenger *fMessenger;
		G4int mEventType; // 0: single particle, 1: particles from file
		G4double mEnergy; // energy is ignored for event type 0
		G4int mParticlePDG; // particle pdg code
		G4String mParticleFile; // file with particles if provided
		G4int mEventNumber; // event number
		G4int mNumberOfEvents; // if event number > mNumberOfEvents, return an empty event

		G4ParticleGun* fParticleGun;
		G4int layerType;
		G4double zPos;
		G4double xPosHole;
		G4double yPosHole;
		G4double holeRadius;
		G4double holeR2;

		// If first layer is type 0
		G4double xCenter;
		G4double yCenter;
		G4double outerRadius;

		// If first layer is type 1
		G4double xLeftBoundary;
		G4double xRightBoundary;
		G4double yDownBoundary;
		G4double yUpBoundary;

		//
		void loadFirstLayerParams(std::vector<G4double> firstLayer);
		void loadParticlesFromFile();
		void generateSingleParticle(G4Event *anEvent);
		void generateParticlesFromFile(G4Event *anEvent);

		struct particle{
			G4int pdgCode;
			G4double px, py, pz;
			G4double x, y, z;
			particle(G4int pdg, G4double momx, G4double momy, G4double momz, G4double posx, G4double posy, G4double posz):
				pdgCode(pdg),
				px(momx),
				py(momy),
				pz(momz),
				x(posx),
				y(posy),
				z(posz)
			{}
		};
		std::map<G4int, std::vector<particle>> mEventMap;

};

#endif