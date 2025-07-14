#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCEventHeader.h"
#endif

#include <iostream>
#include <map>

#include "TFile.h"
#include "TTree.h"
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
        colorMap.insert(pair<string, int>("Positron annihilation", 4));
        colorMap.insert(pair<string, int>("Low photon", 6));
        colorMap.insert(pair<string, int>("Neutron capture", 7));
        colorMap.insert(pair<string, int>("Hadronic inelastic", 8));
        colorMap.insert(pair<string, int>("Electron nuclear interaction", 9));
        colorMap.insert(pair<string, int>("Positron nuclear interaction", 29));
        colorMap.insert(pair<string, int>("Hadronic interaction", 41));
    }
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

struct pi0{
	public:
		pi0(Int_t motherID, Double_t motherE, Double_t motherP, Double_t motherPT, Double_t motherRap, Int_t phot1ID, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t E):
			mMotherPartID(motherID), mMotherPartE(motherE), mMotherPartP(motherP), mMotherPartPT(motherPT), mMotherPartRap(motherRap),
			mPhoton1ID(phot1ID), mPhoton1x(x), mPhoton1y(y), mPhoton1z(z), mPhoton1px(px), mPhoton1py(py), mPhoton1pz(pz), mPhoton1E(E),
			mHasBothPhotons(false){
		}

		Double_t reconstruct_pi0_mass(){
			ROOT::Math::PxPyPzEVector phot1(mPhoton1px, mPhoton1py, mPhoton1pz, mPhoton1E);
			ROOT::Math::PxPyPzEVector phot2(mPhoton2px, mPhoton2py, mPhoton2pz, mPhoton2E);
			
			Double_t pi0Mass = TMath::Sqrt((phot1 + phot2).Dot(phot1 + phot2));

			return pi0Mass;
		}

		Double_t reconstruct_pi0_mass_smeared(){
			ROOT::Math::PxPyPzEVector phot1(mPhoton1px, mPhoton1py, mPhoton1pz, mPhoton1E);
			Double_t smearedE1 = myRand->Gaus(mPhoton1E, mSmearingFactor * mPhoton1E);
			phot1 *= smearedE1/mPhoton1E;
			
			ROOT::Math::PxPyPzEVector phot2(mPhoton2px, mPhoton2py, mPhoton2pz, mPhoton2E);
			Double_t smearedE2 = myRand->Gaus(mPhoton2E, mSmearingFactor * mPhoton2E);
			phot2 *= smearedE2/mPhoton2E;

			Double_t pi0Mass = TMath::Sqrt((phot1 + phot2).Dot(phot1 + phot2));
			return pi0Mass;
		}

		void addPhoton(Int_t phot2ID, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t E){
			if(phot2ID == mPhoton1ID){
				std::cout << "Error. Same photon. Not adding it" << std::endl;
				return;
			}
			mPhoton2ID = phot2ID;
			mPhoton2x = x;
			mPhoton2y = y;
			mPhoton2z = z;
			mPhoton2px = px;
			mPhoton2py = py;
			mPhoton2pz = pz;
			mPhoton2E = E;
			mHasBothPhotons = true;
		}

		Int_t getMotherID(){return mMotherPartID;}
		Double_t getMotherPT(){return mMotherPartPT;}
		Double_t getMotherRap(){return mMotherPartRap;}
		Bool_t hasBothPhotons(){return mHasBothPhotons;}

	private:
		Int_t mMotherPartID;
		Double_t mMotherPartE;
		Double_t mMotherPartP;
		Double_t mMotherPartPT;
		Double_t mMotherPartRap;

		Int_t mPhoton1ID;
		Double_t mPhoton1E;
		Double_t mPhoton1x, mPhoton1y, mPhoton1z;
		Double_t mPhoton1px, mPhoton1py, mPhoton1pz;
		
		Int_t mPhoton2ID;
		Double_t mPhoton2E;
		Double_t mPhoton2x, mPhoton2y, mPhoton2z;
		Double_t mPhoton2px, mPhoton2py, mPhoton2pz;

		Bool_t mHasBothPhotons;
		Double_t mSmearingFactor = 0.15;

};

Double_t reconstruct_pi0_mass(ROOT::Math::PxPyPzEVector phot1, ROOT::Math::PxPyPzEVector phot2, Double_t smearingFactor = 0.15){
	Double_t smearedE1 = myRand->Gaus(phot1.E(), smearingFactor * phot1.E());
	phot1 *= smearedE1/phot1.E();
	
	Double_t smearedE2 = myRand->Gaus(phot2.E(), smearingFactor * phot2.E());
	phot2 *= smearedE2/phot2.E();

	Double_t pi0Mass = TMath::Sqrt((phot1 + phot2).Dot(phot1 + phot2));
	return pi0Mass;
}

void pi0_pTRapidityDist(TString jobDirectory, TString saveDirectory, Int_t nJobs = 100){

	double smearFactor = 0.15; // pT resolution simulation
	Float_t firstLayerZ = 442; // cm

	// For the histograms
	Int_t nBins = 100;
	Double_t histBegin = 0; // GeV/c
	Double_t histEnd = 5; // GeV/c
	Double_t yLowRange = 1e-5;
	Double_t yHighRange = 1e2;

	Double_t rapBegin = 4.;
	Double_t rapEnd = 5.;

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

	productionProcesses prodNames;
    TDatabasePDG* pdg = o2::O2DatabasePDG::Instance();

	TH2D *pi0pTRapidityDist = new TH2D("pi0pTy", "", nBins, 0., 6., nBins, 0., 5.);
	TH1D *pi0pTDist = new TH1D("pi0pT", "", nBins, histBegin, histEnd);

	for(int event = 0; event < nFCTEvents; event++){
		hitsChain->GetEvent(event);
        kineChain->GetEvent(event);
		int nFCTHits = hitFCT->size();

		std::vector<ROOT::Math::PxPyPzEVector> photons;
		std::vector<pi0> pi0s;

		for(Int_t iHit = 0; iHit < nFCTHits; iHit++){ // Apply other cuts
			// Hit
            Hit* thisHit = &(*hitFCT)[iHit];
			Float_t thisX = thisHit->GetStartX();
            Float_t thisY = thisHit->GetStartY();
			Float_t thisZ = thisHit->GetStartZ();
			Float_t thisR = TMath::Sqrt(thisX * thisX + thisY * thisY);
			Float_t thisE = thisHit->GetE();
			// Track
            Int_t thisTrackID = thisHit->GetTrackID();
            MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrackID);
            int thisPID = thisTrack->GetPdgCode();

			// Rejections
			if(!(thisPID == 22)){continue;} // Check if particle is a photon
			if(!(thisZ > firstLayerZ - 0.5 && thisZ < firstLayerZ + 0.5)){continue;} // Check if inside first layer of FCT
            if(!(thisE > 0.05)){continue;} // Check if energy is above 50 MeV

			// Production process. 201 is from Low's theorem
            Int_t thisProdProcess = thisTrack->getProcess();
            // o2::mcgenstatus::MCGenStatusEncoding statusCodes = thisTrack->getStatusCode(); // If a photon has flag 201, this is a low photon
            // Int_t genStatusCode = statusCodes.gen;
            // if(genStatusCode == 201){
            //     thisProdProcess = 201;
            // }

			// True pi0s
			if(thisProdProcess != 0){continue;}
            Int_t motherTrackID = thisTrack->getMotherTrackId();
            if(motherTrackID == -1){ // Check if mother particle exists
                continue;
            }
            MCTrackT<float>* thisMother = &(*mcTr).at(motherTrackID);
            if(abs(thisMother->GetPdgCode()) != 111){continue;}

			std::vector<pi0>::iterator pi0_iter = pi0s.begin();
			Bool_t pi0_found = false;

			while(pi0_iter != pi0s.end()){
				if((*pi0_iter).getMotherID() == motherTrackID){
					(*pi0_iter).addPhoton(thisTrackID, thisX, thisY, thisZ, thisHit->GetPx(), thisHit->GetPy(), thisHit->GetPz(), thisHit->GetE());
					pi0_found = true;
				}
				++pi0_iter;
			}
			if(!pi0_found){
				pi0s.emplace_back(motherTrackID, thisMother->GetP(), thisMother->GetEnergy(), thisMother->GetPt(), thisMother->GetRapidity(), thisTrackID, thisX, thisY, thisZ, thisHit->GetPx(), thisHit->GetPy(), thisHit->GetPz(), thisHit->GetE());
			}
		}

		// True pi0s
		std::vector<pi0>::iterator pi0_iter = pi0s.begin();
		while(pi0_iter != pi0s.end()){
			if((*pi0_iter).hasBothPhotons()){
				pi0pTRapidityDist->Fill((*pi0_iter).getMotherRap(), (*pi0_iter).getMotherPT());
				if((*pi0_iter).getMotherRap() > rapBegin && (*pi0_iter).getMotherRap() < rapEnd){
					pi0pTDist->Fill((*pi0_iter).getMotherPT());
				}
			}
			++pi0_iter;
		}
	}

	// Draw
	TCanvas *c = new TCanvas("c", "My Root Plots", 2000, 1000);
	c->DrawFrame(0, 0, 1, 1);
	c->SetGrid();
	c->cd();

	pi0pTRapidityDist->SetStats(0);
	pi0pTRapidityDist->GetXaxis()->SetTitle("y");
	pi0pTRapidityDist->GetYaxis()->SetTitle("p_{T} (GeV/c)");
	pi0pTRapidityDist->Draw("colz");

	TString pdfFileName = saveDirectory + "/pi0_pT_Rapidity_dist.pdf";
	c->Print(pdfFileName);
	c->Close();
	delete c;

	// Save
	TString rootFileName = saveDirectory + "/pi0_pT_Rapidity_dist.root";
	TFile *saveFile = new TFile(rootFileName, "RECREATE");
	pi0pTRapidityDist->Write();
	saveFile->Close();
	delete saveFile;



	// Draw
	TCanvas *c2 = new TCanvas("c2", "My Root Plots", 2000, 1000);
	c2->DrawFrame(0, 0, 1, 1);
	c2->SetGrid();
	c2->cd();

	Double_t scaleFactor = 1. / (Double_t) nFCTEvents * nBins / (histEnd - histBegin);
	pi0pTDist->SetStats(0);
	pi0pTDist->Scale(scaleFactor);
	pi0pTDist->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	pi0pTDist->GetYaxis()->SetTitle("1/N_{evt} dN/dp_{T} (GeV/c)^{-1}");
	pi0pTDist->Draw("hist");

	TString pdfFileName2 = saveDirectory + "/pi0_pT_dist.pdf";
	c2->Print(pdfFileName2);
	c2->Close();
	delete c2;

	// Save
	TString rootFileName2 = saveDirectory + "/pi0_pT_dist.root";
	TFile *saveFile2 = new TFile(rootFileName2, "RECREATE");
	pi0pTDist->Write();
	saveFile2->Close();
	delete saveFile2;
}