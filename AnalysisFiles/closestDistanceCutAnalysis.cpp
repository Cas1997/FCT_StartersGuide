#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>

void closestDistanceCutAnalysis(TString fileLocation = ""){
	fileLocation = "<Insert your path here>/";
	TString inputFileName = "closestDistance.root";
	
	TString outputFileName = "closestDistanceCutAnalysis.pdf";
	TString outputFileNameSignificance = "closestDistanceSignificance.pdf";

	// New y-axis range
	double yLowRange = 1e-2;
	double yHighRange = 1e1;

	double errorBackground = 0.05;

	// Hists to load
	std::vector<TString> histNames = {"Primary particle emission",
									  "Secondary decay",
									  "Bremstrahlung",
									  "Low photon",
									  "Positron annihilation"};

	// Rename and resize x and y-axis (leave empty if unchanged)
	TString xAxisTitle = "p_{T} (MeV/c)";
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
	TFile* file = new TFile((fileLocation + inputFileName).Data(), "READ");
	TH1D** histArray = new TH1D*[histNames.size()];
	for(int i = 0; i < histNames.size(); i++){
		histArray[i] = (TH1D*)file->Get(histNames[i].Data());
	}

	TH1D *percentageExtBremCut = new TH1D("percentageExtBremCut", "", 100, 0, 20.);
	percentageExtBremCut->GetXaxis()->SetTitle("Closest distance cut (cm)");
	percentageExtBremCut->GetYaxis()->SetTitle("Percentage cut");
	percentageExtBremCut->GetXaxis()->SetTitleSize(xAxisTitleSize);
	percentageExtBremCut->GetXaxis()->SetLabelSize(xAxisLabelSize);
	percentageExtBremCut->GetYaxis()->SetTitleSize(yAxisTitleSize);
	percentageExtBremCut->GetYaxis()->SetLabelSize(yAxisLabelSize);
	percentageExtBremCut->SetStats(0);
	percentageExtBremCut->SetLineColor(kGreen);

	TH1D *percentageIntBremCut = new TH1D("percentageIntBremCut", "", 100, 0, 20.);
	percentageIntBremCut->GetXaxis()->SetTitle("Closest distance cut (cm)");
	percentageIntBremCut->GetYaxis()->SetTitle("Percentage cut");
	percentageIntBremCut->SetStats(0);
	percentageIntBremCut->SetLineColor(kViolet);

	TH1D *significance = new TH1D("significance", "", 100, 0, 20.);
	significance->GetXaxis()->SetTitle("Closest distance cut (cm)");
	significance->GetYaxis()->SetTitle("S/#sigma_{S}");
	significance->GetXaxis()->SetTitleSize(xAxisTitleSize);
	significance->GetXaxis()->SetLabelSize(xAxisLabelSize);
	significance->GetYaxis()->SetTitleSize(yAxisTitleSize);
	significance->GetYaxis()->SetLabelSize(yAxisLabelSize);

	Double_t totalPrimPart = 0.;
	Double_t totalSecDec = 0.;
	Double_t totalExtBrem = 0;
	Double_t totalIntBrem = 0;
	Double_t totalPosAn = 0.;
	for(int bin = 1; bin <= 200; ++bin){
		totalPrimPart += histArray[0]->GetBinContent(bin);
		totalSecDec += histArray[1]->GetBinContent(bin);
		totalExtBrem += histArray[2]->GetBinContent(bin);
		totalIntBrem += histArray[3]->GetBinContent(bin);
		totalPosAn += histArray[4]->GetBinContent(bin);
	}
	Double_t totalBackground = totalPrimPart + totalSecDec + totalExtBrem + totalPosAn;

	Double_t cumPrimPart = 0.;
	Double_t cumSecDec = 0.;
	Double_t cumExtBrem = 0;
	Double_t cumIntBrem = 0;
	Double_t cumPosAn = 0.;
	Double_t cumBackground = 0.;
	for(int bin = 1; bin <= 100; ++bin){
		cumPrimPart += histArray[0]->GetBinContent(bin);
		cumSecDec += histArray[1]->GetBinContent(bin);
		cumExtBrem += histArray[2]->GetBinContent(bin);
		cumIntBrem += histArray[3]->GetBinContent(bin);
		cumPosAn += histArray[4]->GetBinContent(bin);
		cumBackground = cumPrimPart + cumSecDec + cumExtBrem + cumPosAn;

		percentageExtBremCut->SetBinContent(bin, cumExtBrem/totalExtBrem);
		percentageIntBremCut->SetBinContent(bin, cumIntBrem/totalIntBrem);

		if(totalBackground - cumBackground == 0.){
			significance->SetBinContent(bin, 0.);
		} else {
			significance->SetBinContent(bin, (totalIntBrem - cumIntBrem)/(errorBackground * (totalBackground - cumBackground)));
		}
	}

	TCanvas *c = new TCanvas("c", "My Root Plots", xAspect, yAspect);
    c->DrawFrame(0, 0, 1, 1);
    c->SetGrid();
	c->SetLeftMargin(leftMargin);
	c->SetRightMargin(rightMargin);
	c->SetBottomMargin(bottomMargin);
	c->SetTopMargin(topMargin);

	TLegend* legend = new TLegend(0.50, 0.23, 0.94, 0.43);
	legend->SetBorderSize(0);
	legend->AddEntry(percentageExtBremCut, "External bremsstrahlung", "l");
	legend->AddEntry(percentageIntBremCut, "Internal bremsstrahlung", "l");
	legend->SetTextSize(0.04);

	percentageExtBremCut->GetYaxis()->SetRangeUser(0, 1.05);

	percentageExtBremCut->Draw("hist");
	percentageIntBremCut->Draw("same hist");
	legend->Draw();

	// // TLatex t;
    // // t.SetTextSize(textSize);
	// // double starting_value = .90;
	// // double step_size = textSize;
	// // auto tLatexText_iter = texts.begin();
	// // int textCounter = 0;
	// // while(tLatexText_iter != texts.end()){
	// // 	t.DrawLatexNDC(.20, starting_value - textCounter * step_size, (*tLatexText_iter));
	// // 	textCounter++;
	// // 	tLatexText_iter++;
	// // }

	c->Print(fileLocation + outputFileName);
	c->Close();
	delete c;

	TCanvas *c2 = new TCanvas("c2", "My Root Plots", xAspect, yAspect);
    c2->DrawFrame(0, 0, 1, 1);
    c2->SetGrid();
	c2->SetLeftMargin(leftMargin);
	c2->SetRightMargin(rightMargin);
	c2->SetBottomMargin(bottomMargin);
	c2->SetTopMargin(topMargin);

	significance->SetStats(0);
	significance->Draw("hist");

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

	c2->Print(fileLocation + outputFileNameSignificance);	

	file->Close();
}
