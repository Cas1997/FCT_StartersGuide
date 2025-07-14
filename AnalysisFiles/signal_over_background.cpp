#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"

#include <iostream>

void signal_over_background(){

	TString fileLocation = "<Insert your path here>";

	// New y-axis range
	double yLowRange = 1e-2;
	double yHighRange = 1e1;

	// Hists to load
	std::vector<TString> histNames = {"Hadronic decays",
									  "Secondary interactions",
									  "External bremsstrahlung",
									  "Internal bremsstrahlung",
									  "Other sources"};

	// Rename and resize x and y-axis (leave empty if unchanged)
	TString xAxisTitle = "k_{T} (MeV/c)";
	Double_t xAxisTitleSize = 0.06;
	Double_t xAxisLabelSize = 0.05;
	TString yAxisTitle = "S/B";
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
	for(int i = 0; i < histNames.size(); i++){
		if(histNames[i] == histNames[3]){
			signalHist = (TH1D*)file->Get(histNames[i].Data());
		} else {
			backgroundHist->Add((TH1D*)file->Get(histNames[i].Data()));
		}
	}

	TH1D* s_over_b_hist = new TH1D("s_over_bg_Hist", "", 100, 0, 10);
	s_over_b_hist->Divide(signalHist, backgroundHist);
	s_over_b_hist->SetStats(0);

	if(!xAxisTitle.IsNull()){
		s_over_b_hist->GetXaxis()->SetTitle(xAxisTitle);
	}
	s_over_b_hist->GetXaxis()->SetTitleSize(xAxisTitleSize);
	s_over_b_hist->GetXaxis()->SetLabelSize(xAxisLabelSize);
	if(!yAxisTitle.IsNull()){
		s_over_b_hist->GetYaxis()->SetTitle(yAxisTitle);
	}
	s_over_b_hist->GetYaxis()->SetTitleSize(yAxisTitleSize);
	s_over_b_hist->GetYaxis()->SetLabelSize(yAxisLabelSize);



	TCanvas *c = new TCanvas("c", "My Root Plots", xAspect, yAspect);
    c->DrawFrame(0, 0, 1, 1);
    c->SetLogy();
    c->SetGrid();
	c->SetLeftMargin(leftMargin);
	c->SetRightMargin(rightMargin);
	c->SetBottomMargin(bottomMargin);
	c->SetTopMargin(topMargin);

	s_over_b_hist->GetYaxis()->SetRangeUser(yLowRange, yHighRange);
	s_over_b_hist->Draw("hist");

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

	fileLocation.ReplaceAll(".root", "_s_over_bg.pdf");
	c->Print(fileLocation);
	file->Close();
}