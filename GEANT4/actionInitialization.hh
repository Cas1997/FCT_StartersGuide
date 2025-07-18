#ifndef ACTIONINITIALIZATION_HH
#define ACTIONINITIALIZATION_HH

#include "G4VUserActionInitialization.hh"

#include "primaryGenerator.hh"
#include "runAction.hh"
#include "globals.hh"
#include "G4String.hh"

class ActionInitialization : public G4VUserActionInitialization {
	public:
		ActionInitialization(std::vector<G4double> firstLayer, G4String particleFile);
		~ActionInitialization();

		void Build() const override;
		void BuildForMaster() const override;
	private:
		std::vector<G4double> mFirstLayerBoundaries;
		G4String mParticleFile;

};


#endif