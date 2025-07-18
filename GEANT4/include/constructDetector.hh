#ifndef CONSTRUCTDETECTOR_HH
#define CONSTRUCTDETECTOR_HH

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4AutoDelete.hh"

#include "sensitiveDetector.hh"
#include "optrMultiParticleChangeCrossSection.hh"

#include <vector>

class DetectorConstruction : public G4VUserDetectorConstruction {
	public:
		DetectorConstruction(const std::string& FCTCfg, G4bool inclFT3);
		~DetectorConstruction();
		
		G4VPhysicalVolume* Construct() override;
		std::vector<G4double> getFirstLayerBoundaries();

	private:
		void ConstructSDandField() override;
		void createMaterials();
		void loadFCTGeometry(const std::string& FCTCfg);
		void calculateWorldAndBoxPlacement();
		void constructWorld();
		void constructFT3();
		void constructFCT();
		void constructCher();
		void constructPassive();

		// world and mother volumes
		G4Box *mSolidWorld, *mSolidFCTBox, *mSolidFT3Box, *mSolidCherBox;
		G4LogicalVolume *mLogicWorld, *mLogicFCTBox, *mLogicFT3Box, *mLogicCherBox;
		G4VPhysicalVolume *mPhysWorld, *mPhysFCT, *mPhysFT3, *mPhysCher;
		// passive volumes
		G4SubtractionSolid *mSolidCherBoxWalls, *mSolidCherGasVol;
		G4Tubs *mSolidFT3Pipe, *mSolidFCTPipe, *mSolidCherPipe, *mSolidCherCylinderWall;
		G4LogicalVolume *mLogicFT3Pipe, *mLogicFCTPipe, *mLogicCherPipe;
		G4LogicalVolume *mLogicCherBoxWalls, *mLogicCherGasVol, *mLogicCherCylinderWall;
		G4VPhysicalVolume *mPhysFT3Pipe, *mPhysFCTPipe, *mPhysCherPipe;
		G4VPhysicalVolume *mPhysCherBoxWalls, *mPhysCherGasVol, *mPhysCherCylinderWall;
		// list of (sensitive) FCT and FT3 logical and physical layers
		std::vector<G4LogicalVolume*> mFCTLogicLayers;
		std::vector<G4VPhysicalVolume*> mFCTPhysLayers;
		std::vector<G4LogicalVolume*> mFT3LogicLayers;
		std::vector<G4VPhysicalVolume*> mFT3PhysLayers;
		// Cherenkov detector volumes

		
		G4SubtractionSolid *mCherSolidSiPM;
		G4LogicalVolume *mCherLogicSiPM;
		G4VPhysicalVolume *mCherPhysSiPM;

		// world and mother sizes
		G4double mWorldXY;
		G4double mWorldZ;

		G4double mFCTBoxXYLen;
		G4double mFCTBoxZLen;
		G4double mFCTBoxZPos;

		G4bool mIncludeFT3;
		G4double mFT3BoxXYLen; // FT3 is placed at 0., 0., 0.
		G4double mFT3BoxZLen;

		G4double mCherBoxXYLen;
		G4double mCherBoxZLen;
		G4double mCherBoxZPos;
		// Materials
		G4Material *mWorldMat;
		G4Material *mFCTVolMat;
		G4Material *mLayerMat;
		G4Material *mPipeMat;
		G4Material *mCherVolMat;
		G4Material *mSiPMMat;
		G4Material *mCherAerogel, *mCherCO2, *mCherHe, *mCherNe, *mAluminium;
		const G4double mRadLenSi = 9.7; // radiation length silicon in cm
		// Settings
		G4bool mCheckOverlaps;
		G4bool mAddFT3;
		// First layer boundaries
		std::vector<G4double> mBoundaryFirstLayer;


		struct structLayer{
			structLayer(G4int t, G4int layerN, G4double x, G4double y, G4double z, G4double holeR, G4double xHole, G4double yHole, G4double thick, G4double radius):
				type(t),
				layerNumber(layerN),
				xPos(x),
				yPos(y),
				zPos(z),
				holeRadius(holeR),
				xPosHole(xHole),
				yPosHole(yHole),
				thickness(thick),
				outerRadius(radius)
			{}
			structLayer(G4int t, G4int layerN, G4double x, G4double y, G4double z, G4double holeR, G4double xHole, G4double yHole, G4double thick, G4double xL, G4double yL):
				type(t),
				layerNumber(layerN),
				xPos(x),
				yPos(y),
				zPos(z),
				holeRadius(holeR),
				xPosHole(xHole),
				yPosHole(yHole),
				thickness(thick),
				xLen(xL),
				yLen(yL)
			{}
			G4int type; // 0 circular, 1 rectangular
			G4int layerNumber;
			G4double xPos; // cm
			G4double yPos; // cm
			G4double zPos; // cm
			G4double holeRadius; // cm
			G4double xPosHole; // with respect to layer center (cm)
			G4double yPosHole; // with respect to layer center (cm)
			G4double thickness; // cm

			// circular layer attributes
			G4double outerRadius; // cm

			// rectangular layer attributes
			G4double xLen; // x (cm)
			G4double yLen; // y (cm) 
		};

		// list of FCT layers parameters. FT3 layers doesn't get this list, since they don't change
		std::vector<structLayer> mFCTLayerParams;

};



#endif