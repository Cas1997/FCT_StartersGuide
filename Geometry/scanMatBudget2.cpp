#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "DetectorsBase/GeometryManager.h"
#include "ITSBase/GeometryTGeo.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TGeoManager.h>
#include <TH2F.h>
#include <TH1.h>
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

// #include "FairGeoParSet.h"   // for FairGeoParSet
#include <fairlogger/Logger.h> // for LOG, LOG_IF

#include <vector>
#include <string>
#include <algorithm>
#endif

constexpr int nBinsEtaScan = 200;
constexpr float minEtaScan = 3;
constexpr float maxEtaScan = 6;

void vacuumFormMaterial(TGeoMaterial* mat)
{
	constexpr double kAirA = 14.00674;
	constexpr double kAirZ = 7.;
	constexpr double kAirDensity = 0.0; // set for vacuum
	constexpr double kAirRadLen = std::numeric_limits<double>::max();
	// LOGP(info, "Vacuum forming {} ...", mat->GetName());
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

void ComputeMaterialBudget(double rmin, double zEnd, double phiMin, double phiMax, TH1F* xOverX0VsEta){
	// Ancillary function to compute material budget between rmin and rmax

	double x1, y1, z1, x2, y2, z2;
	x1 = 0.;
	y1 = 0.;
	z1 = 0.;
	double rmax;
	
	int n = nBinsEtaScan;
	double binWidthEta = (maxEtaScan - minEtaScan)/((double)n);
	int nPhiBins = 500;
	double binWidthPhi = (phiMax - phiMin)/((double)nPhiBins);
	for (int it = 0; it < n; it++) { // we simulate flat in phi and eta, from Zvtx=0
		double eta = minEtaScan + binWidthEta * (it + 0.5);
		double theta = TMath::ATan(TMath::Exp(-eta)) * 2.;

		if(eta == 0){
			rmax = 50;
		} else {
			rmax = (eta/abs(eta)) * zEnd * TMath::Tan(theta);
		}
		// PHI VS ETA
		double meanPhiContribution = 0.;
		for(int phiBin = 0; phiBin < nPhiBins; ++phiBin){
			double phi = phiMin + binWidthPhi * (phiBin + 0.5);

			// x1 = rmin * TMath::Cos(phi);
			// y1 = rmin * TMath::Sin(phi);
			// z1 = rmin / TMath::Tan(theta);
			x2 = rmax * TMath::Cos(phi);
			y2 = rmax * TMath::Sin(phi);
			z2 = rmax / TMath::Tan(theta);

			auto mparam = o2::base::GeometryManager::meanMaterialBudget(x1, y1, z1, x2, y2, z2);
			meanPhiContribution += mparam.meanX2X0;
			if(eta > -4.5 && eta < -4.4){
				std::cout << "Phi: " << phi << " Eta: " << eta << " Theta: " << theta << std::endl;
				std::cout << "x2: " << x2 << " y2: " << y2 << " z2: " << z2 << std::endl;
				std::cout << "meanX2X0: " << mparam.meanX2X0 << std::endl;
				std::cout << std::endl;
			}
		}

		xOverX0VsEta->Fill(eta, meanPhiContribution/nPhiBins);

	}

}

std::vector<std::string> printMaterialDefinitions(TGeoManager* gman)
{
	std::vector<std::string> materialNames;
	TGeoMedium* med;
	TGeoMaterial* mat;
	char mediaName[50], matName[50], shortName[50];

	int nMedia = gman->GetListOfMedia()->GetEntries();

	LOGP(info, " =================== ALICE 3 Material Properties ================= ");
	LOGP(info, "    A      Z   d (g/cm3)  RadLen (cm)  IntLen (cm)\t Name\n");

	for (int i = 0; i < nMedia; i++) {
		med = (TGeoMedium*)(gman->GetListOfMedia()->At(i));
		mat = med->GetMaterial();
		LOGP(info, "{:5.1f} {:6.1f} {:8.3f} {:13.2f} {:11.1f}\t {}", mat->GetA(), mat->GetZ(), mat->GetDensity(), mat->GetRadLen(), mat->GetIntLen(), mat->GetName());

		std::vector<std::string> tokens;
		std::string matNameStr(mat->GetName());
		if (matNameStr.back() == '$') {
			matNameStr.pop_back();
		}
		std::transform(matNameStr.begin(), matNameStr.end(), matNameStr.begin(), ::toupper);
		// size_t pos = 0;
		// while ((pos = matNameStr.find("_")) != std::string::npos) {
		// 	std::string part = matNameStr.substr(0, pos);
		// 	tokens.push_back(part);
		// 	matNameStr.erase(0, pos + 1);
		// }
		if (matNameStr == "CAVE_AIR_NF") { // Manually manage air_NF
			continue;
		}
		if(matNameStr.find("AIR") != std::string::npos){matNameStr = "AIR";}
		tokens.push_back(matNameStr);

		if (std::find(materialNames.begin(), materialNames.end(), tokens.back()) == materialNames.end()) {
			materialNames.push_back(tokens.back());
		}
		std::cout << "New name: " << matNameStr << std::endl;
	}

	// print material names for debug
	for (auto& name : materialNames) {
		LOGP(info, "Unique material name: {}", name);
	}

	return materialNames;
}

void scanMatBudget2(){
	// Input file
	const string path = "o2sim_geometry.root";
	// Coverage
	const float rmax = 50.;
	const float rmin = 0.1;
	const double firstLayerFCTZ = 440;
	const double phiMin = 0;
	const double phiMax = 2 * TMath::Pi();
	// Visual settings
	gStyle->SetPadTopMargin(0.035);
	gStyle->SetPadRightMargin(0.035);
	gStyle->SetPadBottomMargin(0.14);
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetTitleOffset(1.4, "y");
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetPalette(kSolar);

	std::map<const string, const string> matContributions = {{"AIR", "Air"}, {"ALICE3_PIPE_BERILLIUM", "Beam pipe"}, {"TRK_SILICON", "IRIS tracker"}, {"ALICE3_TRKSERVICES_POLYETHYLENE", "Services polyethylene"}, {"ALICE3_TRKSERVICES_BERYLLIUM", "IRIS vacuum vessel"}, {"RCH_AEROGEL", "RICH aerogel"}, {"RCH_ARGON", "RICH argon"}, {"RCH_SILICON", "RICH silicon"}, {"ALICE3_TRKSERVICES_CERAMIC", "Cold plate"}, {"FT3_SILICON", "Tracker disks"}, {"TF3_SILICON", "TOF"}};
	std::map<const string, int> matColors = {{"Air", 1}, {"Beam pipe", 2}, {"IRIS tracker", 3}, {"IRIS vacuum vessel", 4}, {"RICH aerogel", 5}, {"RICH argon", 6}, {"Cold plate", 7}, {"Services polyethylene", 8}, {"Tracker disks", 95}, {"RICH aerogel", 10}, {"RICH silicon", 11}, {"TOF", 12}};

	TCanvas* canv = new TCanvas("canv", "canv", 2400, 1400);
	
	TLegend* legVsEta = new TLegend(0.25, 0.55, 0.85, 0.95);
	legVsEta->SetFillColor(kWhite);
	legVsEta->SetTextSize(0.025);
	legVsEta->SetHeader("ALICE 3, 0 < #it{#varphi} < #pi, #it{Z}_{vtx} = 0");

	TLegend *leg = new TLegend(0.5, 0.55, 0.75, 0.95);
    leg->SetBorderSize(0);
	leg->SetTextSize(0.04);

	THStack *xOverX0VsEtaStack = new THStack("xOverX0VsEta", ";#eta;#it{X/X0}");

	std::vector<TH1F*> xOverX0VsEta;

	TGeoManager::Import(path.c_str());
	auto materials = printMaterialDefinitions(gGeoManager);

	std::vector<int> colors = {kAzure + 4, kRed + 1};

	// delete gGeoManager; // We re-import the geometry at each iteration 
	// since we set all other materials except the one under consideration to be vacuum
	auto cols = TColor::GetPalette();
	int stackCounter = 0;
	for (size_t iMaterial{0}; iMaterial < materials.size(); ++iMaterial) {
		if (materials[iMaterial] == "ALICE3_PIPE_VACUUM") {
			continue;
		}
		// if (materials[iMaterial] != "BERILLIUM"){// && materials[iMaterial] != "BERYLLIUM") {
		// 	continue;
		// }
		TGeoManager::Import(path.c_str());
		LOGP(info, " ********* Processing material: {} ********* ", materials[iMaterial]);
		auto nMedia = gGeoManager->GetListOfMedia()->GetEntries();
		for (int i = 0; i < nMedia; i++) {
			auto* med = (TGeoMedium*)(gGeoManager->GetListOfMedia()->At(i));
			auto* mat = med->GetMaterial();
			std::string matname{mat->GetName()};
			std::transform(matname.begin(), matname.end(), matname.begin(), ::toupper);
			if (matname.find(materials[iMaterial]) == std::string::npos) {
				// std::cout << "Transforming " << matname << " to vacuum" << std::endl;
				vacuumFormMaterial(mat);
			} else {
				// LOGP(info, "\t {} found as {} element.", materials[iMaterial], iMaterial);
				// std::cout << "Material under consideration: " << matname << std::endl;
			}
		}
		xOverX0VsEta.emplace_back(new TH1F(Form("xOverX0VsEta_step%zu", iMaterial), "", nBinsEtaScan, minEtaScan, maxEtaScan));
		ComputeMaterialBudget(rmin, firstLayerFCTZ, phiMin, phiMax, xOverX0VsEta.back());

		double meanX0vsEta = 0;
		for (int ix = 1; ix <= xOverX0VsEta.back()->GetNbinsX(); ix++) {
    		meanX0vsEta += xOverX0VsEta.back()->GetBinContent(ix);
    	}
    	meanX0vsEta /= xOverX0VsEta.back()->GetNbinsX();
    	
		if (!meanX0vsEta) {
    		LOGP(info, "Material {} not present in this range. Skipping it.", materials[iMaterial]);
			std::cout << std::endl;
    		continue;
    	}
		std::cout << std::endl;
		stackCounter++;
		// xOverX0VsEta.back()->SetFillColor(gROOT->GetColor(stackCounter) - 9);
		xOverX0VsEta.back()->SetFillColor(matColors[matContributions[materials[iMaterial].c_str()]]);

		xOverX0VsEtaStack->Add(xOverX0VsEta.back());

		TString componentName;
		if (matContributions.find(materials[iMaterial].c_str()) == matContributions.end()) {
			// not found
			componentName = materials[iMaterial].c_str();
		} else {
  			// found
			componentName = matContributions[materials[iMaterial].c_str()];
		}
		leg->AddEntry(xOverX0VsEta.back(), componentName, "f");
	
		delete gGeoManager;
	}

	xOverX0VsEtaStack->Draw("hist");
	leg->Draw();
	canv->Print("matBudget_2.pdf");
	canv->Close();
	delete canv;

	// Put combined THStack in one TH1D
	TH1F *histToSave = new TH1F("matBudgetHist", "", nBinsEtaScan, minEtaScan, maxEtaScan);

	std::vector<TH1F*>::iterator histIter = xOverX0VsEta.begin();
	while(histIter != xOverX0VsEta.end()){
		histToSave->Add((*histIter));
		histIter++;
	}
	
	for(int bin = 1; bin <= histToSave->GetNbinsX(); ++bin){
		std::cout << "Bin: " << bin << " content: " << histToSave->GetBinContent(bin) << std::endl;
	}
	// Save in root file
	TFile *file = new TFile("matBudget_2.root", "RECREATE");
	histToSave->Write();
	file->Close();
	delete file;
}