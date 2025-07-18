#include "primaryGenerator.hh"

PrimaryGenerator::PrimaryGenerator(std::vector<G4double> firstLayer, G4String particleFile) :
	mParticleFile(particleFile)
{
	// Set boundary conditions for first layer
	std::cout << "Constructor here!" << std::endl;
	fMessenger = new G4GenericMessenger(this, "/particleGun/", "Particle Gun Settings");
	fMessenger->DeclareProperty("particlePDG", mParticlePDG, "Select your particle pdg code");
	fMessenger->DeclareProperty("energy", mEnergy, "Particle energy for event type 1-2 (MeV)");

	PrimaryGenerator::loadFirstLayerParams(firstLayer);
	if(mParticleFile != ""){
		PrimaryGenerator::loadParticlesFromFile();
		mEventType = 1; // 0: single particle, 1: particles from file
		mEventNumber = 0;
	} else {
		// For single particles
		mParticlePDG = 22; // default is photon
		mEnergy = 100;	// default is 100 MeV
		mEventType = 0; // 0: single particle, 1: particles from file
	}
	
	// // Initiate particle gun
	fParticleGun = new G4ParticleGun();

	std::cout << "Boundaries set. Particle gun made" << std::endl;
}

void PrimaryGenerator::loadFirstLayerParams(std::vector<G4double> firstLayer)
{
	layerType = firstLayer[0];
	zPos = firstLayer[3];
	holeRadius = firstLayer[4];
	holeR2 = holeRadius * holeRadius;
	xPosHole = firstLayer[5];
	yPosHole = firstLayer[6];

	if (layerType == 0)
	{
		xCenter = firstLayer[1];
		yCenter = firstLayer[2];
		outerRadius = firstLayer[7];
	}
	else if (layerType == 1)
	{
		xLeftBoundary = firstLayer[1] - firstLayer[7];
		xRightBoundary = firstLayer[1] + firstLayer[7];
		yDownBoundary = firstLayer[2] - firstLayer[8];
		yUpBoundary = firstLayer[2] + firstLayer[8];
	}
}

void PrimaryGenerator::loadParticlesFromFile(){
	// Read in file. Layers must start after FT3 (so after abs(z) = 361 cm)
	std::ifstream ifs(mParticleFile.c_str());
	if (!ifs.good()) {
    	std::cout << " Invalid Particle File!" << std::endl;
	}
	std::string tempstr;
	G4int eventNumber, nHits, pdgCode;
	G4double px, py, pz, x, y, z;
	G4int nEvents = 0;

	while (std::getline(ifs, tempstr)) {
		std::istringstream iss(tempstr);
		iss >> eventNumber;
		iss >> nHits;
		iss >> pdgCode;
		iss >> px;
		iss >> py;
		iss >> pz;
		iss >> x;
		iss >> y;
		iss >> z;
		// G4cout << px << " " << py << " " << pz << " " << x << " " << y << " " << z << G4endl;
		mEventMap[eventNumber].push_back(particle(pdgCode, px, py, pz, x, y, z));
		if(eventNumber + 1 > nEvents){
			nEvents = eventNumber + 1;
		}
	}
	mNumberOfEvents = nEvents;
	ifs.close();
}

PrimaryGenerator::~PrimaryGenerator()
{
	delete fParticleGun;
}

void PrimaryGenerator::generateSingleParticle(G4Event *anEvent)
{
	G4double xDir;
	G4double yDir;
	G4double zDir;

	fParticleGun->SetNumberOfParticles(1);

	if (layerType == 0)
	{
		G4double outerRadiusR2 = outerRadius * outerRadius;
		G4bool unchosen = true;
		while (unchosen)
		{
			G4double randR = outerRadius * G4UniformRand();
			G4double randTheta = 2. * M_PI * G4UniformRand();
			xDir = xCenter + randR * std::cos(randTheta);
			yDir = yCenter + randR * std::sin(randTheta);
			if ((xDir - xPosHole) * (xDir - xPosHole) + (yDir - yPosHole) * (yDir - yPosHole) > holeR2)
			{
				unchosen = false;
			}
		}
	}
	else if (layerType == 1)
	{
		G4bool unchosen = true;
		while (unchosen)
		{
			xDir = xLeftBoundary + (xRightBoundary - xLeftBoundary) * G4UniformRand();
			yDir = yDownBoundary + (yUpBoundary - yDownBoundary) * G4UniformRand();
			if ((xDir - xPosHole) * (xDir - xPosHole) + (yDir - yPosHole) * (yDir - yPosHole) > holeR2)
			{
				unchosen = false;
			}
		}
	}
	zDir = zPos;

	// default particle kinematic
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	if (mParticlePDG == 22)
	{
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "gamma");
		fParticleGun->SetParticleDefinition(particle);

		mEnergy = 100 + (400 - 100) * G4UniformRand();
	}
	else if (mParticlePDG == 11)
	{
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "e-");
		fParticleGun->SetParticleDefinition(particle);
	}
	else if (mParticlePDG == -11)
	{
		G4ParticleDefinition *particle = particleTable->FindParticle(particleName = "e+");
		fParticleGun->SetParticleDefinition(particle);
	}

	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(xDir * cm, yDir * cm, zDir * cm));
	fParticleGun->SetParticleEnergy(mEnergy * MeV);
	fParticleGun->SetParticlePosition(G4ThreeVector(0. * cm, 0. * cm, 0. * cm));

	// Fire particle
	fParticleGun->GeneratePrimaryVertex(anEvent);
}

void PrimaryGenerator::generateParticlesFromFile(G4Event *anEvent){
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	fParticleGun->SetNumberOfParticles(1);
	G4cout << "-------------------------------------------------------" << G4endl;
	G4cout << "Event number: " << mEventNumber << G4endl;
	G4cout << "-------------------------------------------------------" << G4endl;

	if(mEventNumber < mNumberOfEvents){
		// Loop over eventMap
		std::vector<particle>::iterator particle_iter = mEventMap[mEventNumber].begin();
		while(particle_iter != mEventMap[mEventNumber].end()){
			G4ParticleDefinition *particle = particleTable->FindParticle((*particle_iter).pdgCode);
			if(!particle){
				particle_iter++;
				continue;
			}
			fParticleGun->SetParticleDefinition(particle);
			fParticleGun->SetParticleMomentum(G4ThreeVector((*particle_iter).px * MeV, (*particle_iter).py * MeV, (*particle_iter).pz * MeV));
			fParticleGun->SetParticlePosition(G4ThreeVector((*particle_iter).x * cm, (*particle_iter).y * cm, (*particle_iter).z * cm));
			fParticleGun->GeneratePrimaryVertex(anEvent);
			
			// G4cout << "Shooting particle " << (*particle_iter).pdgCode << G4endl;
			G4cout << "Px: " << (*particle_iter).px << " Py: " << (*particle_iter).py << " Pz: " << (*particle_iter).pz << G4endl;
			G4cout << "X: " << (*particle_iter).x << " Y: " << (*particle_iter).y << " Z: " << (*particle_iter).z << G4endl;
			G4cout << "-------------------------------------------------------" << G4endl;
			particle_iter++;
		}
	} else { // load default event. Very soft photon. Should only happen when the user input of number of event mismatches the number of events in the file
		G4ParticleDefinition *particle = particleTable->FindParticle(22);
		fParticleGun->SetParticleDefinition(particle);
		fParticleGun->SetParticleMomentum(G4ThreeVector(0 * MeV, 0 * MeV, 10 * MeV));
		fParticleGun->SetParticlePosition(G4ThreeVector(0 * cm, 0 * cm, 0 * cm));
		fParticleGun->GeneratePrimaryVertex(anEvent);
	}
	mEventNumber++;
}

void PrimaryGenerator::GeneratePrimaries(G4Event *anEvent)
{
	if(mEventType == 0){
		PrimaryGenerator::generateSingleParticle(anEvent);
	} else if(mEventType == 1){
		PrimaryGenerator::generateParticlesFromFile(anEvent);
	}
}