#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCEventHeader.h"
#endif

#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TVector3.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TRandom3.h"

using o2::MCTrackT;
using o2::itsmft::Hit;

TRandom* myRand = new TRandom3();

struct productionProcesses{
	TString TMCProcessName[50] = {
		"Primary particle emission", // 0
		"Multiple scattering",
		"Energy loss",
		"Bending in magnetic field",
		"Secondary decay", // 4. Called "Decay" in TMCProcess, but that name is confusing
		"Lepton pair production",
		"Compton scattering",
		"Photoelectric effect",
		"Bremstrahlung", // 8
		"Delta ray",
		"Positron annihilation",
		"Positron annihilation at rest",
		"Positron annihilation in flight",
		"Hadronic interaction",
		"Nuclear evaporation",
		"Nuclear fission",
		"Nuclear absorbtion",
		"Antiproton annihilation",
		"Antineutron annihilation",
		"Neutron capture",
		"Hadronic elastic",
		"Hadronic incoherent elastic",
		"Hadronic coherent elastic",
		"Hadronic inelastic",
		"Photon inelastic",
		"Muon nuclear interaction",
		"Electron nuclear interaction",
		"Positron nuclear interaction",
		"Time of flight limit",
		"Nuclear photofission",
		"Rayleigh effect",
		"No active process",
		"Energy threshold",
		"Light absorption",
		"Light detection",
		"Light scattering",
		"Maximum allowed step",
		"Cerenkov production",
		"Cerenkov feed back photon",
		"Cerenkov photon reflection",
		"Cerenkov photon refraction",
		"Synchrotron radiation",
		"Scintillation",
		"Transportation",
		"Unknown process",
		"Coulomb scattering",
		"Photo nuclear interaction",
		"User defined process",
		"Optical photon wavelength shifting",
		"Transition radiation"
	};

	std::map<string, int> colorMap;
	productionProcesses(){
		colorMap.insert(pair<string, int>("Primary particle emission", 1));
		colorMap.insert(pair<string, int>("Secondary decay", 2));
		colorMap.insert(pair<string, int>("Bremstrahlung", 3));
		colorMap.insert(pair<string, int>("Other sources", 4));
		colorMap.insert(pair<string, int>("Positron annihilation", 4));
		colorMap.insert(pair<string, int>("Low photon", 6));
		colorMap.insert(pair<string, int>("Neutron capture", 7));
		colorMap.insert(pair<string, int>("Hadronic inelastic", 8));
		colorMap.insert(pair<string, int>("Electron nuclear interaction", 9));
		colorMap.insert(pair<string, int>("Positron nuclear interaction", 29));
		colorMap.insert(pair<string, int>("Hadronic interaction", 41));
	}
};

struct chargedParticle{
	chargedParticle(Double_t x, Double_t y, Double_t z):
		mX(x), mY(y), mZ(z){};
	Double_t mX, mY, mZ;
};

double calcEta(double x, double y, double z){
	double r = TMath::Sqrt(x*x + y*y);
	double theta;
	if(z > 0.){
		theta = TMath::ATan(r/z);
	} else {
		theta = TMath::Pi() - TMath::ATan(r/abs(z));
	}
	double eta = -TMath::Log(TMath::Tan(theta/2.));
	return eta;
}

double smearPT(double pT, double smearFactor){
	double smearedPT = myRand->Gaus(pT, smearFactor * pT);
	return smearedPT;
}

double calcRecPT(double eta, double P){
	double recPT = P / TMath::CosH(eta);
	return recPT;
}

double calcDistance(double photX, double photY, double chargedX, double chargedY){
	double distance = TMath::Sqrt((photX - chargedX)*(photX - chargedX) + (photY - chargedY)*(photY - chargedY));
	return distance;
}

void photSpectrumClosestDistanceCut(TString jobDirectory, TString saveDirectory, Int_t nJobs = 100){

	Double_t paCutVal = 0.025;
	Double_t smearFactor = 0.15;
	Float_t firstLayerZ = 442.; // cm
	TString outputFileName = "closestDistanceCut_zoomOut";

	// For the histograms
	Int_t nBins = 100;
	Double_t histBegin = 0.;
	Double_t histEnd = 100.;
	Double_t yLowRange = 1e-5;
	Double_t yHighRange = 1e2;

	TChain *hitsChain = new TChain("o2sim");
	for(int iFile = 0; iFile < nJobs; iFile++){
		TString input_path_hits = Form("%s/job%i/o2sim_HitsFCT.root", jobDirectory.Data(), iFile);
		hitsChain->Add(input_path_hits);
	}

	TChain *kineChain = new TChain("o2sim");
	for(int iFile = 0; iFile < nJobs; iFile++){
		TString input_path_kine = Form("%s/job%i/o2sim_Kine.root", jobDirectory.Data(), iFile);
		kineChain->Add(input_path_kine);
	}

	// Hits from FCT
	std::vector<Hit>* hitFCT = nullptr;
	hitsChain->SetBranchAddress("FCTHit", &hitFCT);
	Int_t nFCTEvents = hitsChain->GetEntries();
	std::cout << "checking " << nFCTEvents << " events" << std::endl;

	// Kine file
	std::vector<o2::MCTrack>* mcTr = nullptr;
	kineChain->SetBranchAddress("MCTrack", &mcTr);
	o2::dataformats::MCEventHeader* eventHeader = nullptr;
	kineChain->SetBranchAddress("MCEventHeader.", &eventHeader);
	Int_t nFCTEventsKine = kineChain->GetEntries();

	if(nFCTEvents != nFCTEventsKine){
		std::cout << "Number of events in the kine and hits files do not match!" << std::endl;
		std::cout << "Hits: " << nFCTEvents << std::endl;
		std::cout << "Kine: " << nFCTEventsKine << std::endl;
		std::cout << "Closing event analysis" << std::endl;
		return;
	}

	Int_t nCuts = 7;
	Double_t cdCuts[7] = {1., 1.5, 2., 2.5, 3., 4., 5.};

	std::vector<std::vector<TH1D*>> photonSpectras;
	std::vector<std::vector<Int_t>> allPPPs; // All photon production processes. A histogram for every process will be made
	std::vector<std::map<int, int>> photonSpectraIndices;
	std::vector<TLegend*> legends;
	for(int i = 0; i < nCuts; ++i){
		photonSpectras.push_back({});
		allPPPs.push_back({});
		photonSpectraIndices.push_back({});
		legends.push_back(new TLegend(0.65, 0.65, 0.85, 0.85));
	}

	productionProcesses prodNames;
	TDatabasePDG* pdg = o2::O2DatabasePDG::Instance();

	for(Int_t event = 0; event < nFCTEvents; event++){
		hitsChain->GetEvent(event);
		kineChain->GetEvent(event);
		Int_t nFCTHits = hitFCT->size();

		std::vector<chargedParticle> chargedParticles;
		// Collect charged particles
		for(Int_t iHit = 0; iHit < nFCTHits; ++iHit){
			// Hit
			Hit* thisHit = &(*hitFCT)[iHit];
			Float_t thisX = thisHit->GetStartX();
			Float_t thisY = thisHit->GetStartY();
			Float_t thisZ = thisHit->GetStartZ();
			if(!(thisZ > firstLayerZ - 0.5 && thisZ < firstLayerZ + 0.5)){continue;}
			// Track
			Int_t thisTrackID = thisHit->GetTrackID();
			MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrackID);
			Int_t thisPID = thisTrack->GetPdgCode();
			TParticlePDG* thisPart = pdg->GetParticle(thisPID);
			if(!thisPart){continue;} // Happens for rare particles in excited modes that people do not bother with logging in the pdg database
			
			Double_t charge = 0.;
			charge = thisPart->Charge();
			if(charge != 0.){
				chargedParticles.emplace_back(thisX, thisY, thisZ);
			}
		}

		// Apply pointing angle cut and save pT (smeared) of photons
		for(Int_t iHit = 0; iHit < nFCTHits; ++iHit){
			Hit* thisHit = &(*hitFCT)[iHit];
			Float_t thisX = thisHit->GetStartX();
			Float_t thisY = thisHit->GetStartY();
			Float_t thisZ = thisHit->GetStartZ();
			Float_t thisE = thisHit->GetE();

			Int_t thisTrackID = thisHit->GetTrackID();
			MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrackID);
			Int_t thisPID = thisTrack->GetPdgCode();

			// Rejections
			if(!(thisPID == 22)){continue;} // Check if particle is a photon
			if(!(thisZ > firstLayerZ - 0.5 && thisZ < firstLayerZ + 0.5)){continue;} // Check if inside first layer of FCT
			if(!(thisE > 0.05)){continue;} // Check if energy is above 50 MeV

			Double_t thisRecEta = calcEta(thisX, thisY, thisZ);
			Double_t thisRecPT = 1000. * calcRecPT(thisRecEta, thisE); // in MeV/c
			Double_t smearedRecpT = smearPT(thisRecPT, smearFactor); // in MeV/c
			// if(smearedRecpT > 10){continue;} // To keep it in the range of interest

			// Pointing angle cut
			TVector3 momVecMC(thisHit->GetPx(), thisHit->GetPy(), thisHit->GetPz());
			TVector3 posVecMC(thisHit->GetStartX(), thisHit->GetStartY(), thisHit->GetStartZ());
			Double_t pointingAngleMC = momVecMC.Angle(posVecMC);
			if(pointingAngleMC > paCutVal){continue;}

			// Calculate the distance to the closest charged particle
			Double_t closestDistance = 2.*17.;
			for(int chargedIndex = 0; chargedIndex < chargedParticles.size(); ++chargedIndex){
				Double_t distance = calcDistance(thisX, thisY, chargedParticles[chargedIndex].mX, chargedParticles[chargedIndex].mY);
				if(closestDistance > distance){
					closestDistance = distance;
				}
			}

			// Production process. 201 is from Low's theorem
			Int_t thisProdProcess = thisTrack->getProcess();
			o2::mcgenstatus::MCGenStatusEncoding statusCodes = thisTrack->getStatusCode(); // If a photon has flag 201, this is a low photon
			Int_t genStatusCode = statusCodes.gen;
			if(genStatusCode == 201){
				thisProdProcess = 201;
			}

			for(int i = 0; i < nCuts; ++i){
				if(closestDistance < cdCuts[i]){continue;}
				if(!(std::find(allPPPs[i].begin(), allPPPs[i].end(), thisProdProcess) != allPPPs[i].end())){
					allPPPs[i].push_back(thisProdProcess);
					if(thisProdProcess == 201){
						photonSpectraIndices[i][thisProdProcess] = photonSpectras[i].size();
						photonSpectras[i].push_back(new TH1D(Form("Low photon_%d", i), "", nBins, histBegin, histEnd));
						photonSpectras[i].back()->SetLineColor(prodNames.colorMap["Low photon"]);
						legends[i]->AddEntry(photonSpectras[i].back(), "Low photon", "l");
					} else {
						photonSpectraIndices[i][thisProdProcess] = photonSpectras[i].size();
						photonSpectras[i].push_back(new TH1D(Form("%s_%d", prodNames.TMCProcessName[thisProdProcess].Data(), i), "", nBins, histBegin, histEnd));
						photonSpectras[i].back()->SetLineColor(prodNames.colorMap[prodNames.TMCProcessName[thisProdProcess].Data()]);
						legends[i]->AddEntry(photonSpectras[i].back(), prodNames.TMCProcessName[thisProdProcess], "l");
					}
				}
				photonSpectras[i][photonSpectraIndices[i][thisProdProcess]]->Fill(smearedRecpT);
			}
		}
	}

	for(int i = 0; i < nCuts; ++i){
		// Draw
		TCanvas *c = new TCanvas(Form("c%d", i), "My Root Plots", 2000, 1000);
		c->DrawFrame(0, 0, 1, 1);
		c->SetGrid();
		c->SetLogy();
		c->cd();

		photonSpectras[i][0]->SetStats(0);
		photonSpectras[i][0]->GetXaxis()->SetTitle("k_{T}");
		photonSpectras[i][0]->GetYaxis()->SetTitle("1/N_{evt} dN/dk_{T} (MeV/c)^{-1}");
		legends[i]->SetBorderSize(0);

		double scaleFactor = 1. / (Double_t) nFCTEvents * nBins / (histEnd - histBegin);

		photonSpectras[i][0]->Scale(scaleFactor);
		photonSpectras[i][0]->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
		photonSpectras[i][0]->Draw("hist");
		for(int j = 1; j < photonSpectras[i].size(); ++j){
			photonSpectras[i][j]->Scale(scaleFactor);
			photonSpectras[i][j]->Draw("hist same");
		}
		legends[i]->Draw();

		TString pdfFileName = Form("%s/%s_%1.1f.pdf", saveDirectory.Data(), outputFileName.Data(), cdCuts[i]);
		c->Print(pdfFileName);
		c->Close();
		delete c;
	
		// Save
		TString rootFileName = Form("%s/%s_%1.1f.root", saveDirectory.Data(), outputFileName.Data(), cdCuts[i]);
		TFile *saveFile = new TFile(rootFileName, "RECREATE");
		for(int j = 0; j < photonSpectras[i].size(); ++j){
			photonSpectras[i][j]->Write();
		}
		saveFile->Close();
		delete saveFile;
	}
}