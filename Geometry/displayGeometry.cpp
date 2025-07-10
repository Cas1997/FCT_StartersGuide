#include <string>
#include <unordered_set>

#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "FCTSimulation/Detector.h"
#include "FCTBase/FCTBaseParam.h"
#include "ITSMFTSimulation/Hit.h"
#include "SimulationDataFormat/BaseHits.h"
#include "Steer/MCKinematicsReader.h"
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCEventHeader.h"
#include "SimulationDataFormat/MCUtils.h"
#endif

#include "TCanvas.h"
#include "TEveGeoNode.h"
#include "TEveGeoShape.h"
#include "TEveGeoShapeExtract.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGLOrthoCamera.h"
#include "TGLViewer.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"

#include <fmt/format.h>


TEveTrackPropagator* g_prop = 0;

void displayTracks(){
	// Define target directory - should contain 
	// - o2sim_geometry.root
	// - o2sim_HitsFCT.root
	// - o2sim_Kine.root
	TString directory = "/home/cas/phd/FCTO2_Code/Nazar_track_vis/input/test_folder";

	//__________________________________________
	// event manager
	fmt::print("Creating TEveManager\n");
	TEveManager::Create();

	//__________________________________________
	// geometry
	fmt::print("Importing geometry\n");
	gStyle->SetCanvasPreferGL(true);
	gGeoManager->Import(Form("%s/o2sim_geometry.root", directory.Data()));

	auto* topv = gGeoManager->GetTopVolume();
	fmt::print("Top volume:{}\n", topv->GetName());

	gGeoManager->SetVisLevel(4);
	auto layerColor = kBlack;

	// 	if(strcmp(thisNode->GetName(), "TRKV_2") == 0){
	// 		TEveGeoTopNode* thisEveNode = new TEveGeoTopNode(gGeoManager, thisNode);
	// 		gEve->AddGlobalElement(thisEveNode);
	// 		std::cout << "Added " << thisNode->GetName() << " to viewer" << std::endl;
	// 	}

	// 	// if(strcmp(thisNode->GetName(), "TRKV_2") != 0){continue;}
	// 	// if(strcmp(thisNode->GetName(), "VACUUM_BASE_1") == 0){continue;}
	// 	// thisEveNode->SetVisLevel(1);
	// }

	int nNodes = gGeoManager->GetVolume("barrel")->GetNdaughters();
	for(int i = 0; i < nNodes; i++){
		TGeoNode* thisNode = gGeoManager->GetVolume("barrel")->GetNode(i);
		TEveGeoTopNode* thisEveNode = new TEveGeoTopNode(gGeoManager, thisNode);
		gEve->AddGlobalElement(thisEveNode);
		std::cout << "Added " << thisNode->GetName() << " to viewer" << std::endl;
	}

	// gGeoManager->CheckOverlaps(0.01);
	// gGeoManager->GetListOfOverlaps()->Print();

	gEve->GetDefaultGLViewer()->UseLightColorSet();
	// gEve->GetDefaultGLViewer()->UpdateScene();

	double c[3] = {0, 0, -495};
	gEve->GetDefaultGLViewer()->CurrentCamera().Configure(9., 800., c, -20. / 180. * 3.14, 150. / 180. * 3.14);

	// gEve->GetDefaultGLViewer()->SetGuideState(TGLUtil::kAxesEdge, true, false, nullptr);
	gEve->GetDefaultGLViewer()->DrawGuides();

	gEve->Redraw3D(false);

}
