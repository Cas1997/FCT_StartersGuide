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
#include "TRandom3.h"

using o2::MCTrackT;
using o2::itsmft::Hit;

TRandom* myRand = new TRandom3();

struct pos{
	pos(Double_t radius, Double_t zPos):
		r(radius),
		z(zPos){}
	Double_t r;
	Double_t z; 
};

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
		// Universal production processes
        colorMap.insert(pair<string, int>("Primary particle emission", 1));
		colorMap.insert(pair<string, int>("Secondary decay", 2));
		colorMap.insert(pair<string, int>("Hadronic interaction", 41));
		colorMap.insert(pair<string, int>("Hadronic inelastic", 8));

		// Photon production processes
        colorMap.insert(pair<string, int>("Bremstrahlung", 3));
        colorMap.insert(pair<string, int>("Positron annihilation", 4));
        colorMap.insert(pair<string, int>("Low photon", 6));
        colorMap.insert(pair<string, int>("Neutron capture", 7));
        colorMap.insert(pair<string, int>("Electron nuclear interaction", 9));
        colorMap.insert(pair<string, int>("Positron nuclear interaction", 29));

		// Electron production processes
        colorMap.insert(pair<string, int>("Lepton pair production", 3));
        colorMap.insert(pair<string, int>("Compton scattering", 4));
        colorMap.insert(pair<string, int>("Delta ray", 6));
        colorMap.insert(pair<string, int>("Photoelectric effect", 7));
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

double calcRecPT(double eta, double P){
    double recPT = P / TMath::CosH(eta);
    return recPT;
}

double smearPT(double pT, double smearFactor){
	double smearedPT = myRand->Gaus(pT, smearFactor * pT);
	return smearedPT;
}

double calcDistance(double photX, double photY, double chargedX, double chargedY){
	double distance = TMath::Sqrt((photX - chargedX)*(photX - chargedX) + (photY - chargedY)*(photY - chargedY));
	return distance;
}

struct chargedParticle{
	chargedParticle(Double_t x, Double_t y, Double_t z):
		mX(x), mY(y), mZ(z){};
	Double_t mX, mY, mZ;
};

struct regionRZ{
	public:
		regionRZ(double rMin, double rMax, double zMin, double zMax, TString name):
			mRMin(rMin), mRMax(rMax), mZMin(zMin), mZMax(zMax), mName(name){
				mNSecondaryInteractionPhotons = 0;
				mNExternalBremsstrahlungPhotons = 0;
				mNPositronAnnihilationPhotons = 0;
				mNOtherPhotons = 0;
			}

		bool withinRegion(double r, double z){
			if(r < mRMin || r > mRMax){return false;}
			if(z < mZMin || z > mZMax){return false;}
			return true;
		}

		bool withinRegion(double x, double y, double z){
			double r = TMath::Sqrt(x*x + y*y);
			if(r < mRMin || r > mRMax){return false;}
			if(z < mZMin || z > mZMax){return false;}
			return true;
		}

		void addPhoton(int prodProcess){
			if(prodProcess == 4){
				mNSecondaryInteractionPhotons++;
			} else if(prodProcess == 8){
				mNExternalBremsstrahlungPhotons++;
			} else if(prodProcess == 10){
				mNPositronAnnihilationPhotons++;
			} else {
				mNOtherPhotons++;
			}
			return;
		}
		int getNPhotons(int prodProcess){
			if(prodProcess == 4){
				return mNSecondaryInteractionPhotons;
			} else if(prodProcess == 8){
				return mNExternalBremsstrahlungPhotons;
			} else if(prodProcess == 10){
				return mNPositronAnnihilationPhotons;
			} else {
				return mNOtherPhotons;
			}
		}

		double getRMin(){return mRMin;}
		double getRMax(){return mRMax;}
		double getZMin(){return mZMin;}
		double getZMax(){return mZMax;}
		TString mName;

	private:
		double mRMin, mRMax, mZMin, mZMax;
		int mNSecondaryInteractionPhotons;
		int mNExternalBremsstrahlungPhotons;
		int mNPositronAnnihilationPhotons;
		int mNOtherPhotons;
};

void photonOrigin_improved(TString jobDirectory, TString saveDirectory, Int_t nJobs){
	// Set the cut values
	Double_t paCutVal = 0.025;
	Double_t smearFactor = 0.15;
	Float_t firstLayerZ = 442.; // cm
	Double_t cdCut = 1.; // cm
	
	// Construct the regions
	std::vector<regionRZ> regionsRZ;

	// FT3
	int nFT3Disks = 9;
	double zFT3[9] = {77, 100, 122, 150, 180, 220, 260, 300, 350};
	double rMinFT3[9] = {5., 5., 5., 5., 5., 5., 5., 5., 5.};
	double rMaxFT3[9] = {35., 35., 35., 68., 68., 68., 68., 68., 68.};
	double thicknessFT3 = 0.1; // cm
	for(int i = 0; i < nFT3Disks; ++i){
		regionsRZ.emplace_back(rMinFT3[i], rMaxFT3[i], zFT3[i] - thicknessFT3/2., zFT3[i] + thicknessFT3/2., Form("FT3 disk %d", i));
	}

	// TRK
	int nTRKLayers = 8;
	double zMinTRK[8] = {-62., -62., -62., -62., -62., -129., -129., -129.};
	double zMaxTRK[8] = {62., 62., 62., 62., 62., 129., 129., 129.};
	double rTRK[8] = {7.05, 9.05, 12.05, 20.05, 30.05, 45.05, 60.05, 80.05};
	double thicknessTRK = 0.1;
	for(int i = 0; i < nTRKLayers; ++i){
		regionsRZ.emplace_back(rTRK[i] - thicknessTRK/2., rTRK[i] + thicknessTRK/2., zMinTRK[i], zMaxTRK[i], Form("TRK layer %d", i));
	}

	// Beam pipe (BP) sections
	int nBPSections = 9;
	double zMinBP[9] = {600., 500., 400., 350., 300., 250., 200., 150., 38.5};
	double zMaxBP[9] = {700., 600., 500., 400., 350., 300., 250., 200., 150.};
	double rBP = 1.85; // 
	double thicknessBP = 0.1; // Also for vacuum vessel. Actual thickness is 0.08, but something is going on with the step size.
	for(int i = 0; i < nBPSections; ++i){
		regionsRZ.emplace_back(rBP - thicknessBP/2., rBP + thicknessBP/2., zMinBP[i], zMaxBP[i], Form("Beam pipe section %d", i));
	}

	// Vacuum vessel sections
	double forwardWallRMin = 1.8;
	double forwardWallRMax = 5.6;
	double forwardWallZ = 38.;
	regionsRZ.emplace_back(forwardWallRMin, forwardWallRMax, forwardWallZ - .5, forwardWallZ + .5, "Forward vacuum vessel wall");

	// IRIS tracker
	double zMinIRIS = -35.5;
	double zMaxIRIS = 35.5;
	double rMinIRIS = 0.48;
	double rMaxIRIS = 3.2;
	regionsRZ.emplace_back(rMinIRIS, rMaxIRIS, zMinIRIS, zMaxIRIS, "IRIS");

	// Interaction point
	double zMinIP = -0.5;
	double zMaxIP = 0.5;
	double rMinIP = 0.;
	double rMaxIP = 0.47;
	regionsRZ.emplace_back(rMinIP, rMaxIP, zMinIP, zMaxIP, "Interaction point");

	// Inside the beam pipe vacuum (BPV) (how & why I don't know)
	double zMinBPV[2] = {0., 35.5};
	double zMaxBPV[2] = {35.5, 700.};
	double rMinBPV[2] = {0., 0.};
	double rMaxBPV[2] = {0.48, 1.8};
	for(int i = 0; i < 2; ++i){
		regionsRZ.emplace_back(rMinBPV[i], rMaxBPV[i], zMinBPV[i], zMaxBPV[i], Form("Inside beam pipe vacuum %d", i));
	}

	// Big sections (BS) to narrow down where photons come from if not from any of the above specific sections
	int nBigSections = 10;
	double zMinBS[10] = {500., 300., 100., 100., 100., 100., 100., 100., 0., -500.};
	double zMaxBS[10] = {700., 500., 300., 300., 300., 300., 300., 300., 100., 0.};
	double rMinBS[10] = {1.9, 1.9, 1.9, 5., 10., 15., 20., 30., 1.9, 1.9};
	double rMaxBS[10] = {50., 50., 5., 10., 15., 20., 30., 50., 50., 50.};

	for(int i = 0; i < nBigSections; ++i){
		regionsRZ.emplace_back(rMinBS[i], rMaxBS[i], zMinBS[i], zMaxBS[i], Form("Big section %d", i));
	}

	// Any photons that originate form outside of any of the sections
	int nPhotonsOutside = 0; 
	int totalNPhotons = 0;

	// Load in trees
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
			if(smearedRecpT > 10){continue;} // To keep it in the range of interest
			
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
			if(closestDistance < cdCut){continue;}

			// Production process. 201 is from Low's theorem
			Int_t thisProdProcess = thisTrack->getProcess();
			o2::mcgenstatus::MCGenStatusEncoding statusCodes = thisTrack->getStatusCode(); // If a photon has flag 201, this is a low photon
			Int_t genStatusCode = statusCodes.gen;
			if(genStatusCode == 201){
				thisProdProcess = 201;
			}

			Double_t thisXOrigin = thisTrack->Vx();
			Double_t thisYOrigin = thisTrack->Vy();
			Double_t thisZOrigin = thisTrack->Vz();

			bool inside = false;
			for(int i = 0; i < regionsRZ.size(); ++i){
				inside = regionsRZ[i].withinRegion(thisXOrigin, thisYOrigin, thisZOrigin);
				if(inside){
					regionsRZ[i].addPhoton(thisProdProcess);
					break;
				}
			}
			if(!inside){
				nPhotonsOutside++;
				double r = TMath::Sqrt(thisXOrigin*thisXOrigin + thisYOrigin*thisYOrigin);
				std::cout << "R: " << r << " Z: " << thisZOrigin << " prodProcess: " << thisProdProcess << std::endl;
			}

			totalNPhotons++;
		}
	}

	int nPhotonsInsideRegions = 0;
	int nBremSecIntPosAnihiPhotons = 0;
	for(int i = 0; i < regionsRZ.size(); ++i){
		nBremSecIntPosAnihiPhotons += regionsRZ[i].getNPhotons(4) + regionsRZ[i].getNPhotons(8) + regionsRZ[i].getNPhotons(10);
		nPhotonsInsideRegions += regionsRZ[i].getNPhotons(0) + regionsRZ[i].getNPhotons(4) + regionsRZ[i].getNPhotons(8) + regionsRZ[i].getNPhotons(10);
	}

	for(int i = 0; i < regionsRZ.size(); ++i){
		std::cout << "#################################################" << std::endl;
		std::cout << "Region: " << regionsRZ[i].mName << std::endl;
		std::cout << "rMin: " << regionsRZ[i].getRMin() << " rMax: " << regionsRZ[i].getRMax() << std::endl;
		std::cout << "zMin: " << regionsRZ[i].getZMin() << " zMax: " << regionsRZ[i].getZMax() << std::endl;
		std::cout << "Number of Ext. Brems. photons: " << regionsRZ[i].getNPhotons(8) << std::endl;
		std::cout << "Number of Sec. Inter. photons: " << regionsRZ[i].getNPhotons(4) << std::endl;
		std::cout << "Number of Pos. Anihi. photons: " << regionsRZ[i].getNPhotons(10) << std::endl;
		std::cout << "Number of Other proc. photons: " << regionsRZ[i].getNPhotons(0) << std::endl;
		std::cout << "% of phot from Ext. Brem, Sec. Int., Pos. Anihi.: " << 100. * (regionsRZ[i].getNPhotons(4) + regionsRZ[i].getNPhotons(8) + regionsRZ[i].getNPhotons(10))/(nBremSecIntPosAnihiPhotons + nPhotonsOutside) << std::endl;
		std::cout << "% of phot from all sources   : " << 100. * (regionsRZ[i].getNPhotons(0) + regionsRZ[i].getNPhotons(4) + regionsRZ[i].getNPhotons(8) + regionsRZ[i].getNPhotons(10))/(nPhotonsOutside + nPhotonsInsideRegions) << std::endl;
	}
	std::cout << "#################################################" << std::endl;
	std::cout << "N. Photons originating outside the regions: " << nPhotonsOutside << std::endl;
	std::cout << "N. Photons originating inside the regions : " << nPhotonsInsideRegions << std::endl;
	std::cout << "Sum                                       : " << nPhotonsOutside + nPhotonsInsideRegions << std::endl;
	std::cout << "Total number of photons recorded          : " << totalNPhotons << std::endl;
	std::cout << "% of photons outside of regions			: " << 100. * nPhotonsOutside/totalNPhotons << std::endl;
	std::cout << "#################################################" << std::endl;
}