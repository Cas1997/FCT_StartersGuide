#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>

void remakePlot(){

	//TString fileLocation = "<Insert your path here>"; // Include .root
	// New text to be put in the figure
	std::vector<TString> texts = {"ALICE 3 study",
								  "pp,#sqrt{s} = 14 TeV", 
                                  "FCT: 4 < #eta < 5",
                                  "E_{#gamma} > 50 MeV",
								  "10M Events"};
	// std::vector<TString> texts = {};
	Double_t textSize = 0.04;

	// New y-axis range
	double yLowRange = 1e-5;
	double yHighRange = 1e1;
	Bool_t yAxisLogScale = true;

	// Hists to load. Use the names you gave the histograms
	std::vector<TString> histNames = {"Primary particle emission",
									  "Secondary decay",
									  "Bremstrahlung",
									  "Low photon",
									  "Positron annihilation"};
	// std::vector<TString> histNames = {"pi0MassCBG",
	// 								  "pi0MassTrue"};
	// std::vector<TString> histNames = {"significance"};
	// Legend names - same order and length as histNames. Leave spot empty if same as histNames
	std::vector<TString> legendNames = {"Hadronic decays and direct #gamma's",
									    "Secondary interactions/decays",
									    "External bremsstrahlung",
									    "Internal bremsstrahlung",
										"Other sources"};
	// std::vector<TString> legendNames = {"m_{#pi^{0}} + CBG",
	// 									"m_{#pi^{0}}"};
	Bool_t displayLegend = true;

	std::vector<Int_t> colors = {}; // Leave empty if unchanged, otherwise same number of entries as histNames
	Double_t legendTextSize = 0.04;

	// Rename and resize x and y-axis (leave empty if unchanged)
	TString xAxisTitle = "k_{T} (MeV/c)";
	// TString xAxisTitle = "Closest distance (cm)";
	// TString xAxisTitle = "z (cm)";
	// TString xAxisTitle = "m_{#pi^{0}} (MeV)";
	// TString xAxisTitle = "p_{T} (GeV/c)";
	Double_t xAxisTitleSize = 0.06;
	Double_t xAxisLabelSize = 0.05;
	TString yAxisTitle = "1/N_{evt} dN/dk_{T} (MeV/c)^{-1}";
	// TString yAxisTitle = "1/N_{evt} dN/d(c.d.) (cm)^{-1}";
	// TString yAxisTitle = "1/N_{evt} dN/dz (cm)^{-1}";
	// TString yAxisTitle = "1/N_{evt} dN/dm (MeV)^{-1}";
	// TString yAxisTitle = "#frac{S}{#sqrt{S+B}}";
	Double_t yAxisTitleSize = 0.06;
	Double_t yAxisLabelSize = 0.05;

	// Set title - Empty if no title
	TString title = "";

	// Aspect ratio of canvas
	Int_t xAspect = 1000;
	Int_t yAspect = 800;
	// Left and bottom margin canvas
	// Double_t leftMargin = 0.16;
	Double_t leftMargin = 0.18;
	Double_t rightMargin = 0.05;
	Double_t topMargin = 0.05;
	Double_t bottomMargin = 0.14;
	// Left and bottom margins of the legend box
	Double_t legendLeftMargin = 0.40;
	Double_t legendRightMargin = 0.94;
	Double_t legendBottomMargin = 0.73;
	Double_t legendTopMargin = 0.94;

	// Double_t legendLeftMargin = 0.65;
	// Double_t legendRightMargin = 0.94;
	// Double_t legendBottomMargin = 0.73;
	// Double_t legendTopMargin = 0.94;

	// ###########################################################################
	TFile* file = new TFile(fileLocation.Data(), "READ");
	TH1D** histArray = new TH1D*[histNames.size()];
	TLegend* legend = new TLegend(legendLeftMargin, legendBottomMargin, legendRightMargin, legendTopMargin);
	legend->SetBorderSize(0);
	legend->SetTextSize(legendTextSize);
	for(int i = 0; i < histNames.size(); i++){
		histArray[i] = (TH1D*)file->Get(histNames[i].Data());
		if(!colors.empty()){	
			histArray[i]->SetLineColor(colors[i]);
		}
		if(legendNames[i].IsNull()){
			legend->AddEntry(histArray[i], histNames[i].Data(), "l");
		} else {
			legend->AddEntry(histArray[i], legendNames[i].Data(), "l");
		}
	}

	if(!xAxisTitle.IsNull()){
		histArray[0]->GetXaxis()->SetTitle(xAxisTitle);
	}
	histArray[0]->GetXaxis()->SetTitleSize(xAxisTitleSize);
	histArray[0]->GetXaxis()->SetLabelSize(xAxisLabelSize);
	if(!yAxisTitle.IsNull()){
		histArray[0]->GetYaxis()->SetTitle(yAxisTitle);
	}
	histArray[0]->GetYaxis()->SetTitleSize(yAxisTitleSize);
	histArray[0]->GetYaxis()->SetLabelSize(yAxisLabelSize);

	histArray[0]->SetTitle(title.Data());

	TCanvas *c = new TCanvas("c", "My Root Plots", xAspect, yAspect);
    c->DrawFrame(0, 0, 1, 1);
    if(yAxisLogScale){c->SetLogy();}
    c->SetGrid();
	c->SetLeftMargin(leftMargin);
	c->SetRightMargin(rightMargin);
	c->SetBottomMargin(bottomMargin);
	c->SetTopMargin(topMargin);

	histArray[0]->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	histArray[0]->SetStats(0);
	histArray[0]->Draw("hist");
	for(int i = 1; i < histNames.size(); i++){
		histArray[i]->Draw("same hist");
	}
	if(displayLegend){legend->Draw();}

	TLatex t;
    t.SetTextSize(textSize);
	double starting_value = .91;
	double step_size = textSize;
	auto tLatexText_iter = texts.begin();
	int textCounter = 0;
	while(tLatexText_iter != texts.end()){
		t.DrawLatexNDC(0.20, starting_value - textCounter * step_size, (*tLatexText_iter));
		textCounter++;
		tLatexText_iter++;
	}

	fileLocation.ReplaceAll(".root", "_Remade.pdf");
	c->Print(fileLocation);
	file->Close();
}
