#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>

void significance(TString fileLocation = ""){

	fileLocation = "<Insert your path here>";

	// New y-axis range
	double yLowRange = 0;
	double yHighRange = 20;

	double backgroundUncertainty = 0.05;
	double enhancement = 4.;

	// Hists to load
	std::vector<TString> histNames = {"Primary particle emission",
									  "Secondary decay",
									  "Bremstrahlung",
									  "Low photon",
									  "Positron annihilation"};

	// std::vector<TString> histNames = {"Hadronic decays",
	// 								  "Secondary interactions",
	// 								  "External bremsstrahlung",
	// 								  "Internal bremsstrahlung",
	// 								  "Other sources"};

	// Rename and resize x and y-axis (leave empty if unchanged)
	TString xAxisTitle = "k_{T} (MeV/c)";
	Double_t xAxisTitleSize = 0.06;
	Double_t xAxisLabelSize = 0.05;
	TString yAxisTitle = "s/#sigma_{s}";
	Double_t yAxisTitleSize = 0.06;
	Double_t yAxisLabelSize = 0.05;

	// Set title - Empty if no title
	TString title = "";

	// Aspect ratio of canvas
	Int_t xAspect = 1000;
	Int_t yAspect = 800;
	// Left and bottom margin canvas
	Double_t leftMargin = 0.16;
	Double_t rightMargin = 0.05;
	Double_t topMargin = 0.05;
	Double_t bottomMargin = 0.14;

	// ###########################################################################3
	TFile* file = new TFile(fileLocation.Data(), "READ");
	TH1D* signalHist;
	TH1D* backgroundHist = new TH1D("Background", "", 100, 0, 10);
	TH1D* fullSpectrumHist = new TH1D("FullSpectrum", "", 100, 0, 10);
	for(int i = 0; i < histNames.size(); i++){
		if(histNames[i] == histNames[3]){
			signalHist = (TH1D*)file->Get(histNames[i].Data());
			std::cout << "Done Low" << std::endl;
		} else {
			backgroundHist->Add((TH1D*)file->Get(histNames[i].Data()));
			std::cout << "Done" << std::endl;
		}
		fullSpectrumHist->Add((TH1D*)file->Get(histNames[i].Data()));
	}
	for(int i = 1; i <= fullSpectrumHist->GetNbinsX(); ++i){
		fullSpectrumHist->SetBinError(i, backgroundUncertainty);
	}


	// TH1D* backgroundSubtracted = new TH1D("backgroundSubtracted", "", 100, 0, 10);
	// backgroundSubtracted->Add(fullSpectrumHist);
	// backgroundSubtracted->Add(backgroundHist, -1.);
	// for(int i = 1; i <= backgroundSubtracted->GetNbinsX(); ++i){
	// 	double errorThisbin = backgroundUncertainty * fullSpectrumHist->GetBinContent(i);
	// 	backgroundSubtracted->SetBinError(i, errorThisbin);
	// }
	// backgroundSubtracted->SetStats(0);

	// if(!xAxisTitle.IsNull()){
	// 	backgroundSubtracted->GetXaxis()->SetTitle(xAxisTitle);
	// }
	// backgroundSubtracted->GetXaxis()->SetTitleSize(xAxisTitleSize);
	// backgroundSubtracted->GetXaxis()->SetLabelSize(xAxisLabelSize);
	// if(!yAxisTitle.IsNull()){
	// 	backgroundSubtracted->GetYaxis()->SetTitle(yAxisTitle);
	// }
	// backgroundSubtracted->GetYaxis()->SetTitleSize(yAxisTitleSize);
	// backgroundSubtracted->GetYaxis()->SetLabelSize(yAxisLabelSize);

	TH1D* significance = new TH1D("significance", "", 100, 0, 10);
	significance->Add(fullSpectrumHist);
	significance->Add(backgroundHist, -1.);
	TH1D *enhancedSignificance = new TH1D("enhanced", "", 100, 0, 10);
	enhancedSignificance->Add(significance, 4);
	for(int i = 1; i <= significance->GetNbinsX(); ++i){
		// double errorThisbin = backgroundUncertainty * fullSpectrumHist->GetBinContent(i);
		double errorThisbin = backgroundUncertainty * backgroundHist->GetBinContent(i);
		if(errorThisbin == 0){
			significance->SetBinContent(i, 0);
			enhancedSignificance->SetBinContent(i, 0);
		} else {
			significance->SetBinContent(i, significance->GetBinContent(i)/errorThisbin);
			enhancedSignificance->SetBinContent(i, enhancedSignificance->GetBinContent(i)/errorThisbin);
		}
	}
	significance->SetStats(0);
	significance->SetLineColor(kBlue);
	enhancedSignificance->SetStats(0);

	if(!xAxisTitle.IsNull()){
		enhancedSignificance->GetXaxis()->SetTitle(xAxisTitle);
	}
	enhancedSignificance->GetXaxis()->SetTitleSize(xAxisTitleSize);
	enhancedSignificance->GetXaxis()->SetLabelSize(xAxisLabelSize);
	if(!yAxisTitle.IsNull()){
		enhancedSignificance->GetYaxis()->SetTitle(yAxisTitle);
	}
	enhancedSignificance->GetYaxis()->SetTitleSize(yAxisTitleSize);
	enhancedSignificance->GetYaxis()->SetLabelSize(yAxisLabelSize);

	enhancedSignificance->SetLineColor(kRed);

	// Define TLegends
	TLegend* legend = new TLegend(0.55, 0.75, 0.94, 0.94);
	legend->SetBorderSize(0);
	legend->SetTextSize(0.04);
	legend->AddEntry(enhancedSignificance, "Enhancement", "l");
	legend->AddEntry(significance, "No enhancement", "l");

	TCanvas *c = new TCanvas("c", "My Root Plots", xAspect, yAspect);
    c->DrawFrame(0, 0, 1, 1);
    // c->SetLogy();
    c->SetGrid();
	c->SetLeftMargin(leftMargin);
	c->SetRightMargin(rightMargin);
	c->SetBottomMargin(bottomMargin);
	c->SetTopMargin(topMargin);

	significance->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	enhancedSignificance->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	enhancedSignificance->Draw("hist");
	significance->Draw("hist same");
	legend->Draw();

	// TLatex t;
    // t.SetTextSize(textSize);
	// double starting_value = .90;
	// double step_size = textSize;
	// auto tLatexText_iter = texts.begin();
	// int textCounter = 0;
	// while(tLatexText_iter != texts.end()){
	// 	t.DrawLatexNDC(.20, starting_value - textCounter * step_size, (*tLatexText_iter));
	// 	textCounter++;
	// 	tLatexText_iter++;
	// }

	fileLocation.ReplaceAll(".root", "_significance.pdf");
	c->Print(fileLocation);
	file->Close();
}
