// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Original authors: M. Sitta, F. Grosa
// Author: M. Concas

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TText.h>
#include <THStack.h>

// #include "FairGeoParSet.h"     // for FairGeoParSet
#include <fairlogger/Logger.h> // for LOG, LOG_IF
#include <stdio.h>

#include <vector>
#include <string>
#include <algorithm>
#endif

constexpr int nBinsPhiScan = 90;
constexpr int nBinsEtaScan = 200;
constexpr int nBinsZvtxScan = 300;
constexpr float maxEtaScan = 5.;
constexpr int n = 1e6;      // testtracks
constexpr float len = 1000; // cm

void vacuumFormMaterial(TGeoMaterial* mat)
{
	constexpr double kAirA = 14.00674;
	constexpr double kAirZ = 7.;
	constexpr double kAirDensity = 0.0; // set for vacuum
	constexpr double kAirRadLen = std::numeric_limits<double>::max();
	// std::cout << "Vacuum forming {} ...", mat->GetName());
	// std::cout << "\t A = " << mat->GetA() << " -> " << kAirA << std::endl;
	// std::cout << "\t Z = " << mat->GetZ() << " -> " << kAirZ << std::endl;
	// std::cout << "\t Density = " << mat->GetDensity() << " -> " << kAirDensity << std::endl;
	// std::cout << "\t RadLen = " << mat->GetRadLen() << " -> " << kAirRadLen << std::endl;

	// Make this material air-like
	mat->SetA(kAirA);
	mat->SetZ(kAirZ);
	mat->SetDensity(kAirDensity);
	mat->SetRadLen(kAirRadLen);
}



std::vector<std::string> printMaterialDefinitions(TGeoManager* gman)
{
	std::vector<std::string> materialNames;
	TGeoMedium* med;
	TGeoMaterial* mat;
	char mediaName[50], matName[50], shortName[50];

	int nMedia = gman->GetListOfMedia()->GetEntries();

	std::cout << " =================== ALICE 3 Material Properties ================= " << std::endl;
	std::cout << "    A      Z   d (g/cm3)  RadLen (cm)  IntLen (cm)   Name  " << std::endl;

	for (int i = 0; i < nMedia; i++) {
		med = (TGeoMedium*)(gman->GetListOfMedia()->At(i));
		mat = med->GetMaterial();
		printf("%5.1f %6.1f %8.3f %13.1f %11.1f\t\t%s\n", mat->GetA(), mat->GetZ(), mat->GetDensity(), mat->GetRadLen(), mat->GetIntLen(), mat->GetName());
		std::vector<std::string> tokens;
		std::string matNameStr(mat->GetName());
		if (matNameStr.back() == '$') {
			matNameStr.pop_back();
		}
		size_t pos = 0;
		while ((pos = matNameStr.find("_")) != std::string::npos) {
			std::string part = matNameStr.substr(0, pos);
			tokens.push_back(part);
			matNameStr.erase(0, pos + 1);
		}
		std::transform(matNameStr.begin(), matNameStr.end(), matNameStr.begin(), ::toupper);
		tokens.push_back(matNameStr);
		if (tokens.back() == "NF") { // Manually manage air_NF
			continue;
		}
		if (std::find(materialNames.begin(), materialNames.end(), tokens.back()) == materialNames.end()) {
			materialNames.push_back(tokens.back());
		}
	}

	// print material names for debug
	for (auto& name : materialNames) {
		std::cout << "Unique material name: " << name << std::endl;
	}

	return materialNames;
}

void listMaterials(const float rmax = 350, const float rmin = 0.1, const std::string OnlyMat = "all", const std::string fileName = "o2sim_geometry.root", const string path = "./")
{
	double etaPos = 1.;

	TGeoManager::Import((path + fileName).c_str());
	auto materials = printMaterialDefinitions(gGeoManager);

	const double phiMin = 0;
	const double phiMax = 2 * TMath::Pi();
	const double len = 1200.;

	int count = 2;
	float maxPhiHist = 0.f, maxEtaHist = 0.f, maxZHist = 0.f;
	for (size_t iMaterial{0}; iMaterial < materials.size(); ++iMaterial) {
		if (materials[iMaterial] == "VACUUM") {
			continue;
		}
		if (OnlyMat != "all" && materials[iMaterial] != OnlyMat) {
			continue;
		}
		TGeoManager::Import((path + fileName).c_str());
		std::cout << " ********* Processing material: {} ********* " << materials[iMaterial] << std::endl;
		auto nMedia = gGeoManager->GetListOfMedia()->GetEntries();
		for (int i = 0; i < nMedia; i++) {
			auto* med = (TGeoMedium*)(gGeoManager->GetListOfMedia()->At(i));
			auto* mat = med->GetMaterial();
			std::string matname{mat->GetName()};
			std::transform(matname.begin(), matname.end(), matname.begin(), ::toupper);
			if (matname.find(materials[iMaterial]) == std::string::npos) {
				vacuumFormMaterial(mat);
			} else {
				std::cout << "Original name: " << matname << " Name now: " << materials[iMaterial] << std::endl;
				std::cout << materials[iMaterial] << " found as element " << iMaterial << std::endl;
			}
		}
		delete gGeoManager;
	}
}