#include <iostream>

#include "G4String.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "actionInitialization.hh"
#include "constructDetector.hh"

int main(int argc, char** argv){
	// Fetch arguments
	G4String layout("FCT_circular.cfg");
	G4String executableFile("");
	G4String visualizationFile("vis.mac");
	G4String particleFile("");
	G4bool visualizationToggle(true);
	for(G4int i = 1; i<argc; i=i+2){
		if(G4String(argv[i]) == "-lf") layout = argv[i+1];
		if(G4String(argv[i]) == "-ef") executableFile = argv[i+1];
		if(G4String(argv[i]) == "-pf") particleFile = argv[i+1];
		if(G4String(argv[i]) == "-v"){
			if(G4String(argv[i+1]) == "n" || G4String(argv[i+1]) == "N" || G4String(argv[i+1]) == "no" || G4String(argv[i+1]) == "false"){
				visualizationToggle = false;
			}
		}
		if(G4String(argv[i]) == "-vf") visualizationFile = argv[i+1];
	}

	// Run
	G4RunManager * runManager = new G4RunManager();

	std::cout << layout << std::endl;
	DetectorConstruction* detector = new DetectorConstruction(layout, true);
	runManager->SetUserInitialization(detector);

	G4PhysListFactory factory;
	auto* physicsList = factory.GetReferencePhysList("FTFP_BERT_EMV");
	// Biasing
	G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
	biasingPhysics->Bias("gamma");
	physicsList->RegisterPhysics(biasingPhysics);
	physicsList->RegisterPhysics(new G4OpticalPhysics());
    // physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

	runManager->SetUserInitialization(new ActionInitialization(detector->getFirstLayerBoundaries(), particleFile));
	runManager->Initialize();
	
	G4UIExecutive* ui = 0;
	if(visualizationToggle){ // If no commands are given, run in interactive mode
		ui = new G4UIExecutive(argc, argv);
	}
	
	G4VisManager *visManager = new G4VisExecutive();
	visManager->Initialize();

	G4UImanager* UIManager = G4UImanager::GetUIpointer();

	G4String command = "/control/execute ";
	if(ui){
		UIManager->ApplyCommand(command+visualizationFile);
		ui->SessionStart();
	} else if(executableFile != ""){ // If a command is given, run that macro file, e.g. run.mac
		UIManager->ApplyCommand(command+executableFile);
	} else {
		G4cout << "Either toggle the visualization or provide an executable file." << G4endl;
	}



	// These commands will be moved to a seperate file, but are left here for instructional purposes
	// UIManager->ApplyCommand("/vis/open OGL"); // Open GL
	// UIManager->ApplyCommand("/vis/view/set/viewpointVector 1 1 1"); // Set initial position
	// UIManager->ApplyCommand("/vis/drawVolume"); // Draw all volumes
	// UIManager->ApplyCommand("/vis/viewer/set/autoRefresh"); // Refresh the scene for every run
	// UIManager->ApplyCommand("/vis/scene/add/trajectories smooth"); // Draw the trajectories
	// UIManager->ApplyCommand("/vis/scene/endOfEventAction accumulate"); // Accumulate all events in one run
	// Use command /random/setSeed X Y to set the seed. Very nice

	// Execute the file - CMakeLists gets modified to include all macro files (copies macro files to build folder)

	return 0;
}
