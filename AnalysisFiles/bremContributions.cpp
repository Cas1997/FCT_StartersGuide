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
        "Positron annihilation in flight", // 12
        "Hadronic interaction",
        "Nuclear evaporation",
        "Nuclear fission",
        "Nuclear absorbtion", // 16
        "Antiproton annihilation",
        "Antineutron annihilation",
        "Neutron capture",
        "Hadronic elastic", // 20
        "Hadronic incoherent elastic",
        "Hadronic coherent elastic",
        "Hadronic inelastic",
        "Photon inelastic", // 24
        "Muon nuclear interaction",
        "Electron nuclear interaction",
        "Positron nuclear interaction",
        "Time of flight limit", // 28
        "Nuclear photofission",
        "Rayleigh effect",
        "No active process",
        "Energy threshold", // 32
        "Light absorption",
        "Light detection",
        "Light scattering",
        "Maximum allowed step", // 36
        "Cerenkov production",
        "Cerenkov feed back photon",
        "Cerenkov photon reflection",
        "Cerenkov photon refraction", // 40
        "Synchrotron radiation",
        "Scintillation",
        "Transportation",
        "Unknown process", // 44
        "Coulomb scattering",
        "Photo nuclear interaction",
        "User defined process",
        "Optical photon wavelength shifting", // 48
        "Transition radiation" // 49
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

void bremContributions(TString jobDirectory, TString saveDirectory, Int_t nJobs = 100){

	Float_t firstLayerZ = 442; // cm

	// For the histograms
	Int_t nBins = 100;
	Double_t histBegin = 0.;
	Double_t histEnd = 10.;
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

	// std::vector<TH1D*> bremOrigin;
    TH1D *bremOrigin = new TH1D("BremContributHist", "", nBins, histBegin, histEnd);

    TDatabasePDG* pdg = o2::O2DatabasePDG::Instance();

	TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);

    Int_t electronBremPhotons = 0;
    Int_t chargedParticleBremPhotons = 0;
    Int_t motherDoesNotExist = 0;
    Int_t grandmotherDoesNotExist = 0;
    Int_t grandmother2DoesNotExist = 0;

    // Electrons that produced bremstrahlung
    Int_t pairProductionElectrons = 0;
    Int_t otherProductionElectrons = 0;

    // Electrons that were produced via pair production
    Int_t photonMotherParticleElectrons = 0;
    Int_t otherMotherParticleElectrons = 0;

    // Photons that pair produced electrons
    Int_t piEtaOriginPhotons = 0;
    Int_t otherOriginPhotons = 0;

    // Photons from pi or eta origin that are within the eta acceptance of the FCT
    Int_t withinEtaAcceptance = 0;
    Int_t outsideEtaAcceptance = 0;

    // Pions and Etas coming from primary collision
    Int_t piEtaPrimaryVertex = 0;
    Int_t piEtaOtherOrigin = 0;

	for(int event = 0; event < nFCTEvents; event++){
		hitsChain->GetEvent(event);
        kineChain->GetEvent(event);
		int nFCTHits = hitFCT->size();

		for(Int_t iHit = 0; iHit < nFCTHits; iHit++){ // Apply other cuts
			// Hit
            Hit* thisHit = &(*hitFCT)[iHit];
			Float_t thisX = thisHit->GetStartX();
            Float_t thisY = thisHit->GetStartY();
			Float_t thisZ = thisHit->GetStartZ();
			Float_t thisE = thisHit->GetE();
            Float_t thisPT = TMath::Sqrt(thisHit->GetPx()*thisHit->GetPx() + thisHit->GetPy()*thisHit->GetPy());
			// Track
            Int_t thisTrackID = thisHit->GetTrackID();
            MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrackID);
            int thisPID = thisTrack->GetPdgCode();

			// Rejections
			if(!(thisPID == 22)){continue;} // Check if particle is a photon
			if(!(thisZ > firstLayerZ - 0.5 && thisZ < firstLayerZ + 0.5)){continue;} // Check if inside first layer of FCT
            if(!(thisE > 0.05)){continue;} // Check if energy is above 50 MeV
            if(!(thisPT < 0.01)){continue;} // Check if transverse momentum is below 10 MeV/c
            Int_t thisProdProcess = thisTrack->getProcess();
            if(thisProdProcess != 8){continue;} // Only allow photons which came from bremstrahlung

            
            // Check mother particle
            Int_t motherTrackID = thisTrack->getMotherTrackId();
            if(motherTrackID == -1){ // Check if mother particle exists
                motherDoesNotExist++;
                continue;
            }
            MCTrackT<float>* thisMother = &(*mcTr).at(motherTrackID);
            // if(abs(thisMother->GetPdgCode()) == 11){
            //     electronBremPhotons++;
            //     if(thisMother->getProcess() == 5){
            //         pairProductionElectrons++;
                    
            //     } else {
            //         otherProductionElectrons++;
            //     }
            // } else {
            //     chargedParticleBremPhotons++;
            // }

            if(abs(thisMother->GetPdgCode()) != 11){
                chargedParticleBremPhotons++;
                continue;
            }
            electronBremPhotons++;

            if(thisMother->getProcess() != 5){
                otherProductionElectrons++;
                continue;
            }
            pairProductionElectrons++;

            Int_t grandmotherTrackID = thisMother->getMotherTrackId();
            if(grandmotherTrackID == -1){ // Check if grandmother particle exists
                grandmotherDoesNotExist++;
                continue;
            }
            MCTrackT<float>* thisGrandmother = &(*mcTr).at(grandmotherTrackID);
            if(abs(thisGrandmother->GetPdgCode()) != 22){
                otherMotherParticleElectrons++;
                continue;
            }
            photonMotherParticleElectrons++;

            Int_t grandmother2TrackID = thisGrandmother->getMotherTrackId();
            if(grandmother2TrackID == -1){ // Check if grandmother second order particle exists
                grandmother2DoesNotExist++;
                continue;
            }
            MCTrackT<float>* thisGrandmother2 = &(*mcTr).at(grandmother2TrackID);
            if(abs(thisGrandmother2->GetPdgCode()) == 111 || abs(thisGrandmother2->GetPdgCode()) == 221){
                piEtaOriginPhotons++;
            } else {
                otherOriginPhotons++;
                continue;
            }

            Double_t eta = thisGrandmother->GetEta();
            if(!(eta > 4 && eta < 5.)){
                outsideEtaAcceptance++;
                
                continue;
            }
            withinEtaAcceptance++;

            if(thisGrandmother2->getProcess() != 0){
                piEtaOtherOrigin++;
                continue;
            }
            piEtaPrimaryVertex++;

            bremOrigin->Fill(1000. * thisPT);
		}
	}

    Int_t totalBremPhotons = motherDoesNotExist + electronBremPhotons + chargedParticleBremPhotons;
    std::cout << "Of the " << totalBremPhotons << " bremstrahulng photons" << std::endl;
    std::cout << electronBremPhotons << " came from electrons. " << 100. * ((Double_t) electronBremPhotons)/((Double_t) totalBremPhotons) << " percent" << std::endl;
    std::cout << chargedParticleBremPhotons << " came from charged particles. " << 100. * ((Double_t) chargedParticleBremPhotons)/((Double_t) totalBremPhotons) << " percent" << std::endl;
    std::cout << motherDoesNotExist << " had no mother particle. " << 100. * ((Double_t) motherDoesNotExist)/((Double_t) totalBremPhotons) << " percent" << std::endl;
	std::cout << "############################################" << std::endl;
    std::cout << "Of the electrons that produced bremstrahlung photons" << std::endl;
    std::cout << pairProductionElectrons << " came from pair production. " << 100. * ((Double_t) pairProductionElectrons)/((Double_t) electronBremPhotons) << " percent" << std::endl;
    std::cout << otherProductionElectrons << " came from other production processes. " << 100. * ((Double_t) otherProductionElectrons)/((Double_t) electronBremPhotons) << " percent" << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << "Of the electrons that were produced via pair production" << std::endl;
    std::cout << photonMotherParticleElectrons << " had photons as mother particles. " << 100. * ((Double_t) photonMotherParticleElectrons)/((Double_t) pairProductionElectrons) << " percent" << std::endl;
    std::cout << otherMotherParticleElectrons << " had other particles as mother particles. " << 100. * ((Double_t) otherMotherParticleElectrons)/((Double_t) pairProductionElectrons) << " percent" << std::endl;
    std::cout << grandmotherDoesNotExist << " had no mother particle. " << 100. * ((Double_t) grandmotherDoesNotExist)/((Double_t) pairProductionElectrons) << " percent" << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << "Of the photons that pair produced electrons" << std::endl;
    std::cout << piEtaOriginPhotons << " had pions or etas as mother particles. " << 100. * ((Double_t) piEtaOriginPhotons)/((Double_t) photonMotherParticleElectrons) << " percent" << std::endl;
    std::cout << otherOriginPhotons << " had other particles as mother particles. " << 100. * ((Double_t) otherOriginPhotons)/((Double_t) photonMotherParticleElectrons) << " percent" << std::endl;
    std::cout << grandmother2DoesNotExist << " had no mother particle. " << 100. * ((Double_t) grandmother2DoesNotExist)/((Double_t) photonMotherParticleElectrons) << " percent" << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << "Of the photons with a pion or eta as mother particle" << std::endl;
    std::cout << withinEtaAcceptance << " came from within the FCT acceptance. " << 100. * ((Double_t) withinEtaAcceptance)/((Double_t) piEtaOriginPhotons) << " percent" << std::endl;
    std::cout << outsideEtaAcceptance << " came from outside the FCT acceptance. " << 100. * ((Double_t) outsideEtaAcceptance)/((Double_t) piEtaOriginPhotons) << " percent" << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << "Of the pions and eta that produced a photon in the FCT acceptance" << std::endl;
    std::cout << piEtaPrimaryVertex << " came from the primary vertex. " << 100. * ((Double_t) piEtaPrimaryVertex)/((Double_t) withinEtaAcceptance) << " percent" << std::endl;
    std::cout << piEtaOtherOrigin << " came from other production processes. " << 100. * ((Double_t) piEtaOtherOrigin)/((Double_t) withinEtaAcceptance) << " percent" << std::endl;
    std::cout << "############################################" << std::endl;
    std::cout << 100. * ((Double_t) piEtaPrimaryVertex)/((Double_t) totalBremPhotons) << " percent of the photons find their origin from these pions and etas" << std::endl;
    
    // Draw
	TCanvas *c = new TCanvas("c", "My Root Plots", 2000, 1000);
	c->DrawFrame(0, 0, 1, 1);
	c->SetGrid();
	c->SetLogy();
	c->cd();

	bremOrigin->SetStats(0);
	bremOrigin->GetXaxis()->SetTitle("k_{T} (MeV/c)");
	bremOrigin->GetYaxis()->SetTitle("1/N_{evt} dN/dk_{T} (MeV/c)^{-1}");

	double scaleFactor = 1. / (Double_t) nFCTEvents * nBins / (histEnd - histBegin);

	bremOrigin->Scale(scaleFactor);
	bremOrigin->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	bremOrigin->Draw("hist");

	TString pdfFileName = saveDirectory + "/bremContributions_inclEta.pdf";
	c->Print(pdfFileName);
	c->Close();
	delete c;

	// Save
	TString rootFileName = saveDirectory + "/bremContributions_inclEta.root";
	TFile *saveFile = new TFile(rootFileName, "RECREATE");
	bremOrigin->Write();
	saveFile->Close();
	delete saveFile;
}