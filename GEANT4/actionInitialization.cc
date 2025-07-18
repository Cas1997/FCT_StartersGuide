#include "actionInitialization.hh"

ActionInitialization::ActionInitialization(std::vector<G4double> firstLayer, G4String particleFile) : 
	mFirstLayerBoundaries(firstLayer),
	mParticleFile(particleFile)
	{}

ActionInitialization::~ActionInitialization(){}

void ActionInitialization::BuildForMaster() const{

}

void ActionInitialization::Build() const{
	std::cout << "Attempting to construct generator" << std::endl;
	PrimaryGenerator* generator = new PrimaryGenerator(mFirstLayerBoundaries, mParticleFile);
	std::cout << "Setting the generator" << std::endl;
	SetUserAction(generator);
	std::cout << "Generator okay!" << std::endl;

	RunAction* runAction = new RunAction();
	SetUserAction(runAction);
	std::cout << "RunAction okay!" << std::endl;
}