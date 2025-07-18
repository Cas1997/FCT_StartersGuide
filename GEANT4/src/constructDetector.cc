#include "constructDetector.hh"

DetectorConstruction::DetectorConstruction(const std::string& FCTCfg, G4bool addFT3) : mIncludeFT3(addFT3){
	createMaterials();
	loadFCTGeometry(FCTCfg);
	calculateWorldAndBoxPlacement();
	mCheckOverlaps = true;
}

DetectorConstruction::~DetectorConstruction(){}

std::vector<G4double> DetectorConstruction::getFirstLayerBoundaries(){
	return mBoundaryFirstLayer;
}

void DetectorConstruction::createMaterials(){
	// https://www.fe.infn.it/u/paterno/Geant4_tutorial/slides_further/Geometry/G4_Nist_Materials.pdf
	G4NistManager *nist = G4NistManager::Instance();
	mWorldMat = nist->FindOrBuildMaterial("G4_Galactic");
	mFCTVolMat = nist->FindOrBuildMaterial("G4_AIR");
	mCherVolMat = nist->FindOrBuildMaterial("G4_AIR");
	mLayerMat = nist->FindOrBuildMaterial("G4_Si");
	mSiPMMat = nist->FindOrBuildMaterial("G4_Si");
	mPipeMat = nist->FindOrBuildMaterial("G4_Be");
	mCherHe = nist->FindOrBuildMaterial("G4_He");
	mCherNe = nist->FindOrBuildMaterial("G4_Ne");
	mAluminium = nist->FindOrBuildMaterial("G4_Al");

	// https://www.youtube.com/watch?v=-H-QEIphAsA&list=PLLybgCU6QCGWgzNYOV0SKen9vqg4KXeVL&index=6 13:25
	// Aerogel
	// From https://pdg.lbl.gov/2023/AtomicNuclearProperties/HTML/silica_aerogel.html
	G4Material *SiO2 = new G4Material("SiO2", 2.201*g/cm3, 2);
	SiO2->AddElement(nist->FindOrBuildElement("Si"), 1);
	SiO2->AddElement(nist->FindOrBuildElement("O"), 2);

	G4Material *H2O = new G4Material("H2O", 1.000*g/cm3, 2);
	H2O->AddElement(nist->FindOrBuildElement("H"), 2);
	H2O->AddElement(nist->FindOrBuildElement("O"), 1);

	// Check O2 for percentages
	mCherAerogel = new G4Material("Aerogel", 0.200*g/cm3, 2);
	mCherAerogel->AddMaterial(SiO2, 97.0*perCent);
	mCherAerogel->AddMaterial(H2O, 3.0*perCent);

	// https://pdg.lbl.gov/2023/reviews/rpp2022-rev-atomic-nuclear-prop.pdf
	mCherCO2 = new G4Material("CO2", 0.001842*g/cm3, 2);
	mCherCO2->AddElement(nist->FindOrBuildElement("C"), 1);
	mCherCO2->AddElement(nist->FindOrBuildElement("O"), 2);

	// Set refractive index
	G4double c = 299792458; // speed of light in m/s
	G4double h = 6.62607015e-34; // planck constant in joules
	G4double e = 1.602176634e-19; // elementary charge in coulomb
	G4double photWavelengthToEnergy = h * c / e * 1e6; // Photon wavelength in micrometer to energy conversion factor ~ 1.239841939

	G4double energy[2] = {photWavelengthToEnergy * eV / 0.75, photWavelengthToEnergy * eV / 0.3};
	// Refractive indices
	// CO2: 1.00045
	// He: 1.000036
	G4double rindexHe[2] = {1.000036, 1.000036};
	G4double rindexNe[2] = {1.000063, 1.000063};
	G4double rindexCO2[2] = {1.00045, 1.00045};
	G4double rindexAerogel[2] = {1.03, 1.03};

	G4MaterialPropertiesTable *mptAerogel = new G4MaterialPropertiesTable();
	mptAerogel->AddProperty("RINDEX", energy, rindexAerogel, 2);
	mCherAerogel->SetMaterialPropertiesTable(mptAerogel);

	G4MaterialPropertiesTable *mptCO2 = new G4MaterialPropertiesTable();
	mptCO2->AddProperty("RINDEX", energy, rindexCO2, 2);
	mCherCO2->SetMaterialPropertiesTable(mptCO2);

	G4MaterialPropertiesTable *mptHe = new G4MaterialPropertiesTable();
	mptHe->AddProperty("RINDEX", energy, rindexHe, 2);
	mCherHe->SetMaterialPropertiesTable(mptHe);

	G4MaterialPropertiesTable *mptNe = new G4MaterialPropertiesTable();
	mptNe->AddProperty("RINDEX", energy, rindexNe, 2);
	mCherNe->SetMaterialPropertiesTable(mptNe);

	// Now define refractive index for world material nad SiPMs
	G4double rindexCherVolume[2] = {1.0, 1.0};
	G4MaterialPropertiesTable *mptCherVolume = new G4MaterialPropertiesTable();
	mptCherVolume->AddProperty("RINDEX", energy, rindexCherVolume, 2);
	mCherVolMat->SetMaterialPropertiesTable(mptCherVolume);

	G4double rindexSiPM[2] = {1.0, 1.0};
	G4MaterialPropertiesTable *mptSiPM = new G4MaterialPropertiesTable();
	mptSiPM->AddProperty("RINDEX", energy, rindexSiPM, 2);
	// mptSiPM->AddProperty("EFFICIENCY", energy, rindexSiPM, 2);
	mSiPMMat->SetMaterialPropertiesTable(mptSiPM);

	// Define absorption for aluminium, otherwise some photons will be eternally reflected
	G4MaterialPropertiesTable *mptAl = new G4MaterialPropertiesTable();
	G4double efficiencyAl[2] = {1., 1.};
	mptAl->AddProperty("EFFICIENCY", energy, efficiencyAl, 2);
	mAluminium->SetMaterialPropertiesTable(mptAl);

}

void DetectorConstruction::loadFCTGeometry(const std::string& FCTCfg){
	// Read in file. Layers must start after FT3 (so after abs(z) = 361 cm)
	std::ifstream ifs(FCTCfg.c_str());
	if (!ifs.good()) {
    	std::cout << " Invalid FCT configFile!" << std::endl;
	}
	std::string tempstr;
	int layerNumber = 0;
	int layerType;
	double xPos, yPos, zPos, holeRadius, xPosHole, yPosHole, thickness, outerRadius, xLen, yLen;

	double fctMinZ;
	double fctMaxZ;
	double fctMaxXY;
	double fctMinThickness;

	G4double radLenSi = 9.37;
	int firstLayerNumber = 0;
	while (std::getline(ifs, tempstr)) {
		if(tempstr[0] == '#'){
			int loc_FT3_toggle = tempstr.find("FT3");
      		if (loc_FT3_toggle != -1){mAddFT3 = true;}
			continue;
		}
		std::istringstream iss(tempstr);
		iss >> layerType;
		iss >> xPos;
		iss >> yPos;
		iss >> zPos;
		iss >> holeRadius;
		iss >> xPosHole;
		iss >> yPosHole;
		iss >> thickness;
		
		if(layerType == 0){
			iss >> outerRadius;
			mFCTLayerParams.emplace_back(layerType, layerNumber, xPos, yPos, zPos, holeRadius, xPosHole, yPosHole, thickness * radLenSi, outerRadius);
			std::cout << "Thickness: " << thickness * radLenSi << " cm" << std::endl;
		} else if(layerType == 1) {
			iss >> xLen;
			iss >> yLen;
			mFCTLayerParams.emplace_back(layerType, layerNumber, xPos, yPos, zPos, holeRadius, xPosHole, yPosHole, thickness * radLenSi, xLen, yLen);
		}
		layerNumber++;
	}
}

void DetectorConstruction::calculateWorldAndBoxPlacement(){
	std::vector<structLayer>::iterator iter = mFCTLayerParams.begin();

	G4double xyLenMax = 0.;
	G4double zPosMin = 0.; // Closest to z = 0
	G4double zPosMax = 0.; // Furthest away from z = 0
	G4int direction; // 1: forward, -1 backward (neg z)

	bool firstLayer = true;
	G4int firstLayerNumber;
	for(iter; iter < mFCTLayerParams.end(); iter++){
		if(firstLayer){
			zPosMin = (*iter).zPos;
			direction = (G4int)((*iter).zPos / std::abs((*iter).zPos)); // Assuming all layers are in the same direction
			firstLayer = false;
			firstLayerNumber = (*iter).layerNumber;
		}
		// Check for xyLen
		if((*iter).type == 0){
			if(xyLenMax < std::abs((*iter).xPos) + (*iter).outerRadius){
				xyLenMax = std::abs((*iter).xPos) + (*iter).outerRadius;
			}
			if(xyLenMax < std::abs((*iter).yPos) + (*iter).outerRadius){
				xyLenMax = std::abs((*iter).yPos) + (*iter).outerRadius;
			}
		} else if((*iter).type == 1){
			if(xyLenMax < std::abs((*iter).xPos) + (*iter).xLen){
				xyLenMax = std::abs((*iter).xPos) + (*iter).xLen;
			}
			if(xyLenMax < std::abs((*iter).yPos) + (*iter).yLen){
				xyLenMax = std::abs((*iter).yPos) + (*iter).yLen;
			}
		}
		// Check for z len
		if(zPosMax < std::abs((*iter).zPos)){
			zPosMax = std::abs((*iter).zPos);
		}
		if(zPosMin > std::abs((*iter).zPos)){
			zPosMin = std::abs((*iter).zPos);
			firstLayerNumber = (*iter).layerNumber;
		}

	}

	// Make them slightly bigger so everything is nicely inside their own box
	zPosMax += 10.;
	zPosMin -= 10.;
	xyLenMax += 10.;
	mFCTBoxXYLen = 2 * xyLenMax;
	mFCTBoxZLen = zPosMax - zPosMin;
	mFCTBoxZPos = 0.5 * (zPosMax + zPosMin) * direction;

	if(mFCTBoxXYLen > 160.){
		mWorldXY = mFCTBoxXYLen;
		mFT3BoxXYLen = mFCTBoxXYLen;
	} else {
		mWorldXY = 160.;
		mFT3BoxXYLen = 160.;
	}
	
	mWorldZ = 2. * zPosMax;
	mFT3BoxZLen = mWorldZ - 2 * mFCTBoxZLen;

	// Now add the Cher, which is about 100 cm (currently this includes the box that conatins the gas as well, 
	// so effective photon production length is a little less, but in the order of mm)
	mCherBoxZLen = 100;
	mCherBoxXYLen = mWorldXY;
	mWorldZ += 2 * mCherBoxZLen;
	mCherBoxZPos = mWorldZ / 2. - mCherBoxZLen / 2.;

	std::cout << "###################################" << std::endl;
	std::cout << "###################################" << std::endl;
	std::cout << "FCT Z Len:  " << mFCTBoxZLen << std::endl;
	std::cout << "FCT Z Pos:  " << mFCTBoxZPos << std::endl;
	std::cout << "FT3 Z Len:  " << mFT3BoxZLen << std::endl;
	std::cout << "Cher Z Len: " << mCherBoxZLen << std::endl;
	std::cout << "Cher Z Pos: " << mCherBoxZPos << std::endl;
	std::cout << "World Z Len: " << mWorldZ << std::endl;
	std::cout << "World XYLen: " << mWorldXY << std::endl;
	std::cout << "###################################" << std::endl;
	std::cout << "###################################" << std::endl;

	// Save first layer params - Used to determine where to fire the photons in PrimaryGeneratorAction
	// The "first layer" is the one closest to z = 0
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].type);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].xPos);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].yPos);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].zPos);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].holeRadius);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].xPosHole);
	mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].yPosHole);
	if(mFCTLayerParams[firstLayerNumber].type == 0){
		mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].outerRadius);
	} else {
		mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].xLen);
		mBoundaryFirstLayer.push_back(mFCTLayerParams[firstLayerNumber].yLen);
	}

	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(mWorldZ);

}

void DetectorConstruction::constructWorld(){
	mSolidWorld = new G4Box("World", 0.5 * mWorldXY * cm, 0.5 * mWorldXY * cm, 0.5 * mWorldZ * cm);
	mLogicWorld = new G4LogicalVolume(mSolidWorld, mWorldMat, "World");
	mPhysWorld = new G4PVPlacement(nullptr, G4ThreeVector(), mLogicWorld, "World", nullptr, false, 0, mCheckOverlaps);
}

void DetectorConstruction::constructFT3(){
	mSolidFT3Box = new G4Box("FT3SolidBox", 0.5 * mFT3BoxXYLen*cm, 0.5 * mFT3BoxXYLen*cm, 0.5 * mFT3BoxZLen*cm);
	mLogicFT3Box = new G4LogicalVolume(mSolidFT3Box, mWorldMat, "FT3LogicVol");
	mPhysFT3 = new G4PVPlacement(nullptr,
							     G4ThreeVector(),
								 mLogicFT3Box,
								 "FT3PhysVol",
								 mLogicWorld,
								 false,
								 0,
								 true);
	
	int nLayers = 12;

	// FT3 Layers are the same in the positive as negative direction
	std::vector<std::array<G4double, 4>> FT3LayersConfig{
    	{26., .5, 2.5, 0.1f * 0.01f * mRadLenSi}, // {z_layer, r_in, r_out, Layerx2X0}
    	{30., .5, 2.5, 0.1f * 0.01f * mRadLenSi},
    	{34., .5, 2.5, 0.1f * 0.01f * mRadLenSi},
    	{77., 5.0, 35., 0.01f * mRadLenSi},
    	{100., 5.0, 35., 0.01f * mRadLenSi},
    	{122., 5.0, 35., 0.01f * mRadLenSi},
    	{150., 5.0, 68.f, 0.01f * mRadLenSi},
    	{180., 5.0, 68.f, 0.01f * mRadLenSi},
    	{220., 5.0, 68.f, 0.01f * mRadLenSi},
    	{260., 5.0, 68.f, 0.01f * mRadLenSi},
    	{300., 5.0, 68.f, 0.01f * mRadLenSi},
    	{350., 5.0, 68.f, 0.01f * mRadLenSi}
	};

	for(int layerNumber = 0; layerNumber < nLayers; layerNumber++){
		std::ostringstream strs_solidLayer;
		strs_solidLayer << "FT3SolidLayer_" << layerNumber;
		G4String solidLayerName = strs_solidLayer.str();
		G4Tubs* solidLayer = new G4Tubs(solidLayerName,
								  FT3LayersConfig[layerNumber][1]*cm,
								  FT3LayersConfig[layerNumber][2]*cm,
								  0.5 * FT3LayersConfig[layerNumber][3]*cm,
								  0.*deg,
								  360.*deg);
		
		std::ostringstream strs_logicLayer;
		strs_logicLayer << "FT3LogicLayer_" << layerNumber;
		G4String logicLayerName = strs_logicLayer.str();		
		G4LogicalVolume* logicLayer = new G4LogicalVolume(solidLayer, mLayerMat, logicLayerName);
		mFT3LogicLayers.emplace_back(logicLayer);

		std::ostringstream strs_physLayerPos;
		strs_physLayerPos << "FT3PhysLayer_posZ_" << layerNumber;
		G4String physLayerPosName = strs_physLayerPos.str();				
		G4VPhysicalVolume* physLayerPos = new G4PVPlacement(nullptr, 
																G4ThreeVector(0.*cm, 0.*cm, FT3LayersConfig[layerNumber][0]*cm),
																logicLayer,
																physLayerPosName,
																mLogicFT3Box,
																false,
																layerNumber * 2,
																mCheckOverlaps);
		mFT3PhysLayers.emplace_back(physLayerPos);
		
		std::ostringstream strs_physLayerNeg;
		strs_physLayerNeg << "FT3PhysLayer_negZ_" << layerNumber;
		G4String physLayerNegName = strs_physLayerNeg.str();			
		G4VPhysicalVolume* physLayerNeg = new G4PVPlacement(nullptr, 
																G4ThreeVector(0.*cm, 0.*cm, -FT3LayersConfig[layerNumber][0]*cm),
																logicLayer,
																physLayerNegName,
																mLogicFT3Box,
																false,
																layerNumber * 2 + 1,
																mCheckOverlaps);
		mFT3PhysLayers.emplace_back(physLayerNeg);

	}

	mSolidFT3Pipe = new G4Tubs("SolidPipeFT3",
							   (3.7)*cm,
							   (3.7 + 0.08)*cm,
							   0.5*mFT3BoxZLen*cm,
							   0.*deg,
  							   360.*deg);
	
	mLogicFT3Pipe = new G4LogicalVolume(mSolidFT3Pipe, mPipeMat, "LogicPipeFT3");

	mPhysFT3Pipe = new G4PVPlacement(nullptr,
									 G4ThreeVector(),
									 mLogicFT3Pipe,
									 "PhysPipeFT3",
									 mLogicFT3Box,
									 false,
									 0,
									 mCheckOverlaps);


}

void DetectorConstruction::constructFCT(){
	mSolidFCTBox = new G4Box("SolidFCTBox", 0.5 * mFCTBoxXYLen * cm, 0.5 * mFCTBoxXYLen * cm, 0.5 * mFCTBoxZLen *cm);
	mLogicFCTBox = new G4LogicalVolume(mSolidFCTBox, mFCTVolMat, "LogicFCTBox");
	mPhysFCT = new G4PVPlacement(nullptr,
								 G4ThreeVector(0. * cm, 0. * cm, mFCTBoxZPos * cm),
								 mLogicFCTBox,
								 "FCT",
								 mLogicWorld,
								 false,
								 0,
								 mCheckOverlaps);


    // mLogicWorld->SetVisAttributes(boxVisAtt);
	// mFCTLogicBox->SetVisAttributes(boxVisAtt);

	std::vector<structLayer>::iterator iter = mFCTLayerParams.begin();
	for(iter; iter < mFCTLayerParams.end(); iter++){
		if((*iter).type == 0){
			std::ostringstream strs_solidLayerDisc;
			strs_solidLayerDisc << "FCTSolidLayerDisc_" << (*iter).layerNumber;
			G4String solidLayerDiscName = strs_solidLayerDisc.str();
			G4Tubs* solidLayerDisc = new G4Tubs(solidLayerDiscName,
									  0.*cm,
									  (*iter).outerRadius*cm,
									  0.5*(*iter).thickness*cm,
									  0.*deg,
									  360.*deg);
			
			std::ostringstream strs_solidLayerHole;
			strs_solidLayerHole << "FCTSolidLayerHole_" << (*iter).layerNumber;
			G4String solidLayerHoleName = strs_solidLayerHole.str();
			G4Tubs* solidLayerHole = new G4Tubs(solidLayerHoleName,
									  0.*cm,
									  (*iter).holeRadius*cm,
									  0.5*(*iter).thickness*cm,
									  0.*deg,
									  360.*deg);
			
			std::ostringstream strs_solidLayer;
			strs_solidLayer << "FCTSolidLayer_" << (*iter).layerNumber;
			G4String solidLayerName = strs_solidLayer.str();			
			G4SubtractionSolid* solidLayer = new G4SubtractionSolid(solidLayerName, solidLayerDisc, solidLayerHole, nullptr, G4ThreeVector((*iter).xPosHole*cm, (*iter).yPosHole*cm, 0.*cm));
			
			std::ostringstream strs_logicLayer;
			strs_logicLayer << "FCTLogicLayer_" << (*iter).layerNumber << "_Log";
			G4String logicLayerName = strs_logicLayer.str();
			G4LogicalVolume* logicLayer = new G4LogicalVolume(solidLayer, mLayerMat, logicLayerName);
			mFCTLogicLayers.emplace_back(logicLayer);

			std::ostringstream strs_physLayer;
			strs_physLayer << "FCTPhysLayer_" << (*iter).layerNumber;
			G4String physLayerName = strs_physLayer.str();

			G4VPhysicalVolume* physLayer = new G4PVPlacement(nullptr, 
															 G4ThreeVector((*iter).xPos*cm, (*iter).yPos*cm, ((*iter).zPos - mFCTBoxZPos)*cm),
															 logicLayer,
															 physLayerName,
															 mLogicFCTBox,
															 false,
															 (*iter).layerNumber,
															 mCheckOverlaps);
			mFCTPhysLayers.emplace_back(physLayer);
		} else if((*iter).type == 1){
			std::ostringstream strs_solidLayerSquare;
			strs_solidLayerSquare << "FCTSolidLayerSquare_" << (*iter).layerNumber;
			G4String solidLayerSquareName = strs_solidLayerSquare.str();
			G4Box* solidLayerSquare = new G4Box(solidLayerSquareName,
									   (*iter).xLen*cm,
									   (*iter).yLen*cm,
									   (*iter).thickness*cm);
			
			std::ostringstream strs_solidLayerHole;
			strs_solidLayerHole << "FCTSolidLayerHole_" << (*iter).layerNumber;
			G4String solidLayerHoleName = strs_solidLayerHole.str();
			G4Tubs* solidLayerHole = new G4Tubs(solidLayerHoleName,
									  0.*cm,
									  (*iter).holeRadius*cm,
									  0.5*(*iter).thickness*cm,
									  0.*deg,
									  360.*deg);
			
			std::ostringstream strs_solidLayer;
			strs_solidLayer << "FCTSolidLayer_" << (*iter).layerNumber;
			G4String solidLayerName = strs_solidLayer.str();			
			G4SubtractionSolid* solidLayer = new G4SubtractionSolid(solidLayerName, solidLayerSquare, solidLayerHole, nullptr, G4ThreeVector((*iter).xPosHole*cm, (*iter).yPosHole*cm, 0.*cm));
			
			std::ostringstream strs_logicLayer;
			strs_logicLayer << "FCTLogicLayer_" << (*iter).layerNumber;
			G4String logicLayerName = strs_logicLayer.str();
			G4LogicalVolume* logicLayer = new G4LogicalVolume(solidLayer, mLayerMat, logicLayerName);
			mFCTLogicLayers.emplace_back(logicLayer);

			std::ostringstream strs_physLayer;
			strs_physLayer << "FCTPhysLayer_" << (*iter).layerNumber;
			G4String physLayerName = strs_physLayer.str();
			G4VPhysicalVolume* physLayer = new G4PVPlacement(nullptr, 
															 G4ThreeVector((*iter).xPos*cm, (*iter).yPos*cm, ((*iter).zPos - mFCTBoxZPos)*cm),
															 logicLayer,
															 physLayerName,
															 mLogicFCTBox,
															 false,
															 (*iter).layerNumber,
															 mCheckOverlaps);
			mFCTPhysLayers.emplace_back(physLayer);
		}
	}

	// Make the FCT part of the pipe
	mSolidFCTPipe = new G4Tubs("SolidPipeFCT",
							   (3.7)*cm,
							   (3.7 + 0.08)*cm,
							   0.5*mFCTBoxZLen*cm,
							   0.*deg,
  							   360.*deg);
	
	mLogicFCTPipe = new G4LogicalVolume(mSolidFCTPipe, mPipeMat, "LogicPipeFCT");

	mPhysFCTPipe = new G4PVPlacement(nullptr,
									 G4ThreeVector(),
									 mLogicFCTPipe,
									 "PhysPipeFCT",
									 mLogicFCTBox,
									 false,
									 0,
									 mCheckOverlaps);


}

void DetectorConstruction::constructCher(){
	mSolidCherBox = new G4Box("SolidCherBox", 0.5 * mCherBoxXYLen * cm, 0.5 * mCherBoxXYLen * cm, 0.5 * mCherBoxZLen *cm);
	mLogicCherBox = new G4LogicalVolume(mSolidCherBox, mCherVolMat, "LogicCherBox");
	mPhysCher = new G4PVPlacement(nullptr,
								 G4ThreeVector(0. * cm, 0. * cm, mCherBoxZPos * cm),
								 mLogicCherBox,
								 "Cher",
								 mLogicWorld,
								 false,
								 0,
								 mCheckOverlaps);
	
	// Make hollow cube
	G4double wallThickness = 0.01; // cm
	G4Box *cherGasCube = new G4Box("CherGasCube", 0.5 * (mCherBoxXYLen - 2. * wallThickness) * cm, 0.5 * (mCherBoxXYLen - 2. * wallThickness) * cm, 0.5 * (mCherBoxZLen - 2. * wallThickness) *cm);
	G4SubtractionSolid *cherCubeWalls = new G4SubtractionSolid("HollowBox", mSolidCherBox, cherGasCube, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	G4Tubs *beamPipeHole = new G4Tubs("beamPipeHole", 0.*cm, (5. + wallThickness)*cm, 0.5*mCherBoxZLen*cm, 0.*deg, 360.*deg);
	
	mSolidCherBoxWalls = new G4SubtractionSolid("SolidCherDetBoxWalls", cherCubeWalls, beamPipeHole, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	mLogicCherBoxWalls = new G4LogicalVolume(mSolidCherBoxWalls, mAluminium, "LogicCherDetBoxWalls");
	mPhysCherBoxWalls = new G4PVPlacement(nullptr,
								      	  G4ThreeVector(0. * cm, 0. * cm, 0. * cm),
								          mLogicCherBoxWalls,
								          "PhysCherDetBoxWalls",
								          mLogicCherBox,
								          false,
								          0,
								          mCheckOverlaps);

	mSolidCherGasVol = new G4SubtractionSolid("SolidCherGasVol", cherGasCube, beamPipeHole, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	mLogicCherGasVol = new G4LogicalVolume(mSolidCherGasVol, mCherNe, "LogicCherGasVol");
	mPhysCherGasVol = new G4PVPlacement(nullptr,
								      	G4ThreeVector(0. * cm, 0. * cm, 0. * cm),
								        mLogicCherGasVol,
								        "PhysCherGasVol",
								        mLogicCherBox,
								        false,
								        0,
								        mCheckOverlaps);

	mSolidCherCylinderWall = new G4Tubs("SolidCherDetCylinderWalls", 5.*cm, (5. + wallThickness)*cm, 0.5*mCherBoxZLen*cm, 0.*deg, 360.*deg);
	mLogicCherCylinderWall = new G4LogicalVolume(mSolidCherCylinderWall, mAluminium, "LogicCherDetCylinderWalls");
	mPhysCherCylinderWall = new G4PVPlacement(nullptr,
								      		  G4ThreeVector(0. * cm, 0. * cm, 0. * cm),
								      		  mLogicCherCylinderWall,
								      		  "PhysCherCylinderWall",
								      		  mLogicCherBox,
								      		  false,
								      		  0,
								      		  mCheckOverlaps);

	G4double SiPMThickness = 0.2; // cm
	G4Box *SiPMSquare = new G4Box("SiPMSquare", 0.5 * (mCherBoxXYLen - 2. * wallThickness) * cm, 0.5 * (mCherBoxXYLen - 2. * wallThickness) * cm, 0.5 * SiPMThickness * cm);
	G4Tubs* SiPMHole = new G4Tubs("SiPMHole",
								  0.*cm,
								  (5 + wallThickness) * cm,
								  0.5 * SiPMThickness * cm,
								  0. * deg,
								  360. * deg);	
	mCherSolidSiPM = new G4SubtractionSolid("SolidCherSiPM", SiPMSquare, SiPMHole, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	mCherLogicSiPM = new G4LogicalVolume(mCherSolidSiPM, mSiPMMat, "LogicCherSiPM");
	mCherPhysSiPM = new G4PVPlacement(nullptr,
								      G4ThreeVector(0. * cm, 0. * cm, 0.5 * (mCherBoxZLen - 2. * wallThickness - 1. * SiPMThickness) * cm),
								      mCherLogicSiPM,
								      "Cher_SiPM",
								      mLogicCherGasVol,
								      false,
								      100,
								      mCheckOverlaps);

	// G4double aerogelThickness = 2;
	// G4double radiationDistance = 20;
	// G4double relAerogelZPos = -mCherBoxZLen / 2. + aerogelThickness / 2. + 1.;
	// G4Box *aerogelSquare = new G4Box("AerogelSquare", 0.5 * mCherBoxXYLen * cm, 0.5 * mCherBoxXYLen * cm, 0.5 * aerogelThickness * cm);
	// G4Tubs* aerogelHole = new G4Tubs("AerogelHole",
	// 								 0.*cm,
	// 								 5 * cm,
	// 								 0.5 * aerogelThickness * cm,
	// 								 0. * deg,
	// 								 360. * deg);
	// mSolidCherAerogel = new G4SubtractionSolid("SolidCherAerogel", aerogelSquare, aerogelHole, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	// mLogicCherAerogel = new G4LogicalVolume(mSolidCherAerogel, mCherAerogel, "LogicRICGAerogel");
	// mPhysCherAerogel = new G4PVPlacement(nullptr,
	// 									 G4ThreeVector(0. * cm, 0. * cm, relAerogelZPos * cm),
	// 									 mLogicCherAerogel,
	// 									 "Cher_Aerogel",
	// 									 mLogicCherBox,
	// 									 false,
	// 									 99,
	// 									 mCheckOverlaps);
	
	// // Place sensitive volume
	// G4double radiationDistance = 20;
	// G4double SiPMThickness = aerogelThickness;
	// G4Box *SiPMSquare = new G4Box("SiPMSquare", 0.5 * mCherBoxXYLen * cm, 0.5 * mCherBoxXYLen * cm, 0.5 * SiPMThickness * cm);
	// G4Tubs* SiPMHole = new G4Tubs("SiPMHole",
	// 							  0.*cm,
	// 							  5 * cm,
	// 							  0.5 * SiPMThickness * cm,
	// 							  0. * deg,
	// 							  360. * deg);	
	// mCherSolidSiPM = new G4SubtractionSolid("SolidCherSiPM", SiPMSquare, SiPMHole, nullptr, G4ThreeVector(0.*cm, 0.*cm, 0.*cm));
	// mCherLogicSiPM = new G4LogicalVolume(mCherSolidSiPM, mCherVolMat, "LogicCherSiPM");
	// mCherPhysSiPM = new G4PVPlacement(nullptr,
	// 							      G4ThreeVector(0. * cm, 0. * cm, (relAerogelZPos + radiationDistance + 0.5 * SiPMThickness) * cm),
	// 							      mCherLogicSiPM,
	// 							      "Cher_SiPM",
	// 							      mLogicCherBox,
	// 							      false,
	// 							      100,
	// 							      mCheckOverlaps);
	
	// Make the Cher part of the pipe
	mSolidCherPipe = new G4Tubs("SolidPipeCher",
							   (3.7)*cm,
							   (3.7 + 0.08)*cm,
							   0.5*mCherBoxZLen*cm,
							   0.*deg,
  							   360.*deg);
	
	mLogicCherPipe = new G4LogicalVolume(mSolidCherPipe, mPipeMat, "LogicPipeCher");

	mPhysCherPipe = new G4PVPlacement(nullptr,
									 G4ThreeVector(),
									 mLogicCherPipe,
									 "PhysPipeCher",
									 mLogicCherBox,
									 false,
									 0,
									 mCheckOverlaps);
}

void DetectorConstruction::constructPassive(){

}

G4VPhysicalVolume* DetectorConstruction::Construct(){
	constructWorld(); // done
	constructFT3();  // done
	constructFCT(); // done
	constructCher();
	constructPassive();
	return mPhysWorld;
}

void DetectorConstruction::ConstructSDandField(){
	SensitiveDetector *sensDet = new SensitiveDetector("SensitiveDetector");
	std::vector<G4LogicalVolume*>::iterator iterFCT = mFCTLogicLayers.begin();
	for(iterFCT; iterFCT < mFCTLogicLayers.end(); iterFCT++){
		(*iterFCT)->SetSensitiveDetector(sensDet);
	}

	OptrMultiParticleChangeCrossSection *bias = 
		new OptrMultiParticleChangeCrossSection();

	bias->AddParticle("gamma");
	bias->AttachTo(mFCTLogicLayers[0]);
	std::cout << "Attaching biasing operator " << bias->GetName()
			  << "to logical volume " << mFCTLogicLayers[0]->GetName() << std::endl;

	// std::vector<G4LogicalVolume*>::iterator iterFT3 = mFT3LogicLayers.begin();
	// for(iterFT3; iterFT3 < mFT3LogicLayers.end(); iterFT3++){
	// 	(*iterFT3)->SetSensitiveDetector(sensDet);
	// }

	mCherLogicSiPM->SetSensitiveDetector(sensDet);
	// mLogicCherAerogel->SetSensitiveDetector(sensDet);

	// Set the magnetic fields
	G4MagneticField *magFieldFCT;
	magFieldFCT = new G4UniformMagField(G4ThreeVector(0., 0.25 * tesla, 0.));
	G4FieldManager *magFieldFCTManager = new G4FieldManager();
	magFieldFCTManager->SetDetectorField(magFieldFCT);
	magFieldFCTManager->CreateChordFinder(magFieldFCT);
	mLogicFCTBox->SetFieldManager(magFieldFCTManager, true);
	G4AutoDelete::Register(magFieldFCT);
	G4AutoDelete::Register(magFieldFCTManager);

	G4MagneticField *magFieldFT3;
	magFieldFT3 = new G4UniformMagField(G4ThreeVector(0., 0., 2.0 * tesla));
	G4FieldManager *magFieldFT3Manager = new G4FieldManager();
	magFieldFT3Manager->SetDetectorField(magFieldFT3);
	magFieldFT3Manager->CreateChordFinder(magFieldFT3);
	mLogicFT3Box->SetFieldManager(magFieldFT3Manager, true);
	G4AutoDelete::Register(magFieldFT3);
	G4AutoDelete::Register(magFieldFT3Manager);
}