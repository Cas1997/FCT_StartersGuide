#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "SimulationDataFormat/MCTrack.h"
#include "SimulationDataFormat/MCEventHeader.h"
#endif

#include <iostream>
#include <map>

#include "TFile.h"
#include "TLine.h"
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
        "Primary particle emission",
        "Multiple scattering",
        "Energy loss",
        "Bending in magnetic field",
        "Secondary decay", // Called "Decay" in TMCProcess, but that name is confusing
        "Lepton pair production",
        "Compton scattering",
        "Photoelectric effect",
        "Bremstrahlung",
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

std::map<int, string> particles = {{11, "e^{#pm}"}, {13, "#mu^{#pm}"}, {211, "#pi^{#pm}"}, {321, "K^{#pm}"}, {0, "other"}, {2212, "p^{#pm}"}};

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

void chargedParticleSpectrum(TString jobDirectory, TString saveDirectory, Int_t nJobs = 100){

	TString outputFileName = "chargedParticleSpectrum";

	// For the histograms
	Int_t nBins = 200;
	Double_t histBegin = 0.;
	Double_t histEnd = 50.;
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

    Float_t firstLayerZ = 442; // cm
    Float_t lastLayerZ = 500.; // cm
	std::vector<TH1D*> chargedParticleSpectra;
	std::map<Int_t, Bool_t> abundantParticles = {{11, true}, {13, true}, {211, true}, {321, true}, {2212, true}};
	std::map<int, int> chargedParticleSpectraIndex = {{11, 1}, {13, 2}, {211, 3}, {321, 4}, {2212, 5}};
	
	TLegend* legend = new TLegend(0.65, 0.65, 0.85, 0.85);
	
	chargedParticleSpectra.emplace_back(new TH1D("Other", "", nBins, histBegin , histEnd));
	chargedParticleSpectra.emplace_back(new TH1D("Electron", "", nBins, histBegin , histEnd));
	chargedParticleSpectra.emplace_back(new TH1D("Muon", "", nBins, histBegin , histEnd));
	chargedParticleSpectra.emplace_back(new TH1D("Pion", "", nBins, histBegin , histEnd));
	chargedParticleSpectra.emplace_back(new TH1D("Kaon", "", nBins, histBegin , histEnd));
	chargedParticleSpectra.emplace_back(new TH1D("Proton", "", nBins, histBegin , histEnd));

	for(int i = 0; i < 6; ++i){
		chargedParticleSpectra[i]->SetStats(0);
		chargedParticleSpectra[i]->SetLineColor(i + 1);
		legend->AddEntry(chargedParticleSpectra[i], chargedParticleSpectra[i]->GetName(), "l");
	}

	TH1D *electronSpectrumDetail = new TH1D("Electron detail", "", 500, 0, 25);
	Int_t nElectrons = 0;

    TDatabasePDG* pdg = o2::O2DatabasePDG::Instance();

	for(int event = 0; event < nFCTEvents; event++){
		hitsChain->GetEvent(event);
        kineChain->GetEvent(event);
		int nFCTHits = hitFCT->size();

		for(Int_t iHit = 0; iHit < nFCTHits; iHit++){
			// Hit
            Hit* thisHit = &(*hitFCT)[iHit];
			Float_t thisX = thisHit->GetStartX();
            Float_t thisY = thisHit->GetStartY();
            Float_t thisZ = thisHit->GetStartZ();
            Float_t mcP = TMath::Sqrt(thisHit->GetPx() * thisHit->GetPx() + thisHit->GetPy() * thisHit->GetPy() + thisHit->GetPz() * thisHit->GetPz());
			// Track
            Int_t thisTrackID = thisHit->GetTrackID();
            MCTrackT<float>* thisTrack = &(*mcTr).at(thisTrackID);
            int thisPID = abs(thisTrack->GetPdgCode());

			// Rejections
			if(!(thisZ > firstLayerZ - 0.5 && thisZ < firstLayerZ + 0.5)){continue;} // Check if inside first layer of FCT

            // Check if charged particle
            TParticlePDG* thisPartPDG = pdg->GetParticle(thisPID);
            if(!thisPartPDG){continue;} // Happens for rare particles in excited modes that people do not bother with logging in the pdg database
            
            Double_t charge = 0.;
            charge = thisPartPDG->Charge();
            if(charge == 0){
				continue;
            }

			if(abundantParticles[thisPID]){
				chargedParticleSpectra[chargedParticleSpectraIndex[thisPID]]->Fill(mcP);
			} else {
				chargedParticleSpectra[0]->Fill(mcP);
			}

			if(thisPID == 11){
				electronSpectrumDetail->Fill(mcP);
				nElectrons++;
			}
		}
	}

	// Draw
	TCanvas *c = new TCanvas("c", "My Root Plots", 2000, 1000);
	c->DrawFrame(0, 0, 1, 1);
	c->SetGrid();
	c->SetLogy();
	c->cd();

	chargedParticleSpectra[0]->SetStats(0);
	chargedParticleSpectra[0]->GetXaxis()->SetTitle("k_{T} (MeV/c)");
	chargedParticleSpectra[0]->GetYaxis()->SetTitle("1/N_{evt} dN/dk_{T} (MeV/c)^{-1}");
	legend->SetBorderSize(0);

	double scaleFactor = 1. / (Double_t) nFCTEvents * nBins / (histEnd - histBegin);

	chargedParticleSpectra[0]->Scale(scaleFactor);
	chargedParticleSpectra[0]->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	chargedParticleSpectra[0]->Draw("hist");
	for(int i = 1; i < chargedParticleSpectra.size(); ++i){
		chargedParticleSpectra[i]->Scale(scaleFactor);
		chargedParticleSpectra[i]->Draw("hist same");
	}
	legend->Draw();

	TString pdfFileName = saveDirectory + "/" + outputFileName + ".pdf";
	c->Print(pdfFileName);
	c->Close();
	delete c;

	// Save
	TString rootFileName = saveDirectory + "/" + outputFileName + ".root";
	TFile *saveFile = new TFile(rootFileName, "RECREATE");
	for(int i = 0; i < chargedParticleSpectra.size(); ++i){
		chargedParticleSpectra[i]->Write();
	}
	saveFile->Close();
	delete saveFile;


	// Draw & save cumulative electron cut
	TH1D *cumElectronCut = new TH1D("Cumulative electron cut", "", 500, 0, 25);
	cumElectronCut->GetXaxis()->SetTitle("p (GeV/c)");
	cumElectronCut->GetYaxis()->SetTitle("% e^{#pm} cut");
	cumElectronCut->SetStats(0);
	Int_t cumElectron = 0;
	for(int bin = 1; bin <= 500; ++bin){
		cumElectron += electronSpectrumDetail->GetBinContent(bin);
		cumElectronCut->SetBinContent(bin, ((Double_t)cumElectron)/((Double_t)nElectrons));
	}
	// Momentum thresholds
	Double_t n[7] = {1.003, 1.0014, 1.0008, 1.0005, 1.00045, 1.000063, 1.000036};
	Double_t pionMass = 0.13957039; // GeV
	Double_t pth[7] = {0, 0, 0, 0, 0, 0, 0};
	for(int i = 0; i < 7; ++i){
		pth[i] = pionMass / TMath::Sqrt(n[i]*n[i] - 1.);
		std::cout << pth[i] << std::endl;
	}

	TCanvas *c2 = new TCanvas("c2", "My Root Plots", 1600, 1000);
	c2->DrawFrame(0, 0, 1, 1);
	c2->SetGrid();
	c2->cd();
	
	cumElectronCut->GetYaxis()->SetRangeUser(0.2, 1.05);
	cumElectronCut->Draw("hist");
	for(int i = 0; i < 7; ++i){
		TLine *line = new TLine(pth[i], 0.2, pth[i], 1.05);
		line->SetLineColor(kRed);
		line->Draw("same");
	}

	TString pdfFileName2 = saveDirectory + "/cumulativeElectronCut.pdf";
	c2->Print(pdfFileName2);
	c2->Close();
	delete c2;

	// Save
	TString rootFileName2 = saveDirectory + "/cumulativeElectronCut.root";
	TFile *saveFile2 = new TFile(rootFileName2, "RECREATE");
	cumElectronCut->Write();
	saveFile2->Close();
	delete saveFile2;	


}

