/*
Author: Daniel Kikoła, 2024 - 2025.
This macro uses CERN ROOT and RooUnfold frameworks. 
*/

#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TRandom.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

TStyle* setWindowDressing();

void unfoldHadrSmearing(TString outFileName = "Unfolding-results.root", bool useExternalPseudoData = false){

    int N_DDbar_pairs = 1000;

    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);


    setWindowDressing();
    gROOT->ForceStyle();

    TString inFileName = "Response-Matrix-Smearing-2fm-T150MeV.root"; 
    TString inFileNameWithResponseMatrix = "Response-Matrix-Smearing-3fm-T150MeV.root";  
    TString inFileNamePseudoData =   "Unfolding-results.root";


    TFile* fIn = new TFile(inFileName);    
    TFile* fInResponseMatrix = new TFile(inFileNameWithResponseMatrix);    
    TFile* fInPseudoData = new TFile(inFileNamePseudoData);    



    TH1D* hTrue = (TH1D*)fInResponseMatrix->Get("truth");
    TH1D* hMeasured = (TH1D*)fInResponseMatrix->Get("measured");
    TH2D* hResponse = (TH2D*)fInResponseMatrix->Get("response");

    TH1D* hDeltaPhiGaussianSmearing_orig = (TH1D*)fIn->Get("truth");
    TH1D* hMeasured_true = (TH1D*)fIn->Get("measured");
    hDeltaPhiGaussianSmearing_orig->SetName("hDeltaPhiGaussianSmearing_orig");
    TH1D* hDeltaPhiGaussianSmearing_unfolded = (TH1D*)hDeltaPhiGaussianSmearing_orig->Clone("hDeltaPhiGaussianSmearing_unfolded");
    
    TH1D* hDeltaPhiGaussianSmearing_exp = 0;
    TH1D* hDeltaPhiGaussianSmearing_exp_original_binning = 0;
    
    if (useExternalPseudoData)
    {
        hDeltaPhiGaussianSmearing_exp = (TH1D*)fInPseudoData->Get("hDeltaPhiGaussianSmearing_exp_original_binning");
        hDeltaPhiGaussianSmearing_exp->SetName("hDeltaPhiGaussianSmearing_exp");
    } else {
        // create new set of pseudo-data
        hDeltaPhiGaussianSmearing_exp = (TH1D*)hDeltaPhiGaussianSmearing_orig->Clone("hDeltaPhiGaussianSmearing_exp");
        hDeltaPhiGaussianSmearing_exp->Reset();
        hDeltaPhiGaussianSmearing_exp->Sumw2();
        hDeltaPhiGaussianSmearing_exp->FillRandom(hMeasured_true, N_DDbar_pairs);
    }

    hDeltaPhiGaussianSmearing_exp_original_binning = (TH1D*)hDeltaPhiGaussianSmearing_exp->Clone("hDeltaPhiGaussianSmearing_exp_original_binning");


    hDeltaPhiGaussianSmearing_unfolded->Reset();
    hDeltaPhiGaussianSmearing_orig->Sumw2();
    hDeltaPhiGaussianSmearing_unfolded->Sumw2();


    // create response matrix
   //Parameters: const TH1* measured, const TH1* truth, const TH2* response,
    RooUnfoldResponse response (hMeasured,hTrue,hResponse);

    TCanvas* c1 = new TCanvas("c","c",1200,600);
    c1->Divide(2,1);

    c1->cd(1);


    TH1D* hTrueCrossCheck = (TH1D*)response.Htruth();
    hTrueCrossCheck->Draw("hist");
    
    TH1D* hMeasuredCrosscheck = (TH1D*)response.Hmeasured();
    hMeasuredCrosscheck->SetLineColor(kRed);
    hMeasuredCrosscheck->SetMarkerColor(kRed);
    hMeasuredCrosscheck->Draw("hist same");

    hTrueCrossCheck->SetXTitle("#Phi");
    hTrueCrossCheck->SetYTitle("dN/d#Phi");


    TLegend* legendResponse= new TLegend(0.35,0.75,0.9,0.95);
    legendResponse->SetTextSize(0.04);
    legendResponse->AddEntry(hTrueCrossCheck,"Response: True dist.","l");
    legendResponse->AddEntry(hMeasuredCrosscheck,"Response: Measured dist.","l");
    legendResponse->Draw("same");


   
    cout << "==================================== UNFOLD ===================================" << endl;
    RooUnfoldBayes  unfold (&response, hDeltaPhiGaussianSmearing_exp,10);

    
    hDeltaPhiGaussianSmearing_unfolded = (TH1D*) unfold.Hunfold();
    hDeltaPhiGaussianSmearing_unfolded->SetLineColor(kRed);

    hDeltaPhiGaussianSmearing_unfolded->SetLineColor(kBlack);
    hDeltaPhiGaussianSmearing_exp->SetLineColor(kBlue);

    hDeltaPhiGaussianSmearing_exp->Rebin(10);
    hDeltaPhiGaussianSmearing_unfolded->Rebin(10);
    hDeltaPhiGaussianSmearing_orig->Rebin(10);

    // correct the stat. uncertties to include the stat. ones and those from unfolding
    for(int i = 1; i<= hDeltaPhiGaussianSmearing_unfolded->GetNbinsX();i++){
        double err2 = TMath::Power(hDeltaPhiGaussianSmearing_unfolded->GetBinError(i),2) + hDeltaPhiGaussianSmearing_unfolded->GetBinContent(i);
        hDeltaPhiGaussianSmearing_unfolded->SetBinError(i,TMath::Sqrt(err2));
    } 

    c1->cd(2);

    hDeltaPhiGaussianSmearing_exp->SetLineColor(kRed);
    hDeltaPhiGaussianSmearing_unfolded->SetLineColor(kBlue);
    hDeltaPhiGaussianSmearing_exp->SetMarkerColor(kRed);
    hDeltaPhiGaussianSmearing_unfolded->SetMarkerColor(kBlue);
    hDeltaPhiGaussianSmearing_exp->Scale(1./hDeltaPhiGaussianSmearing_exp->Integral("width"));
    hDeltaPhiGaussianSmearing_orig->Scale(1./hDeltaPhiGaussianSmearing_orig->Integral("width"));
    hDeltaPhiGaussianSmearing_unfolded->Scale(1./hDeltaPhiGaussianSmearing_unfolded->Integral("width"));

    hDeltaPhiGaussianSmearing_orig->SetMaximum(1.1);
    hDeltaPhiGaussianSmearing_orig->SetXTitle("#Phi");
    hDeltaPhiGaussianSmearing_orig->SetYTitle("dN/d#Phi");    

    hDeltaPhiGaussianSmearing_orig->Draw("hist e");
    hDeltaPhiGaussianSmearing_exp->Draw("hist same e");
    hDeltaPhiGaussianSmearing_unfolded->Draw("hist same e");

    TLegend* legend= new TLegend(0.35,0.75,0.9,0.95);
    legend->SetTextSize(0.04);
    legend->AddEntry(hDeltaPhiGaussianSmearing_orig,"Gauss","l");
    legend->AddEntry(hDeltaPhiGaussianSmearing_exp,"Experiment (Gauss + Hadr.)","l");
    legend->AddEntry(hDeltaPhiGaussianSmearing_unfolded,"Experiment unfolded","l");
    legend->Draw("same");


    hDeltaPhiGaussianSmearing_unfolded->SetName("hDeltaPhiGaussianSmearing_unfolded");


   TFile* fOut = new TFile(outFileName,"recreate");  
   hDeltaPhiGaussianSmearing_orig->Write();
   hDeltaPhiGaussianSmearing_exp->Write();
   hDeltaPhiGaussianSmearing_unfolded->Write();
   hDeltaPhiGaussianSmearing_exp_original_binning->Write();

   fOut->Close();


}

TStyle* setWindowDressing(){

    float LegFontSize = 0.09;
    float LegMargin  = 0.25;
    float MarkerSize  = 2;
    float LineWidth  = 1.0;

    TStyle *plotsStyle= new TStyle("plotsStyle","plain plots style");
    //
    // use plain black on white colors

    plotsStyle->SetFrameBorderMode(0);
    plotsStyle->SetCanvasBorderMode(0);
    plotsStyle->SetPadBorderMode(0);
    plotsStyle->SetPadColor(0);
    plotsStyle->SetCanvasColor(0);
    plotsStyle->SetStatColor(0);
    plotsStyle->SetFrameFillColor(0);
    plotsStyle->SetEndErrorSize(1.);
    plotsStyle->SetTitleFillColor(0);

    // set the paper & margin sizes
    plotsStyle->SetPaperSize(20,26);
    plotsStyle->SetPadTopMargin(0.02);
    plotsStyle->SetPadRightMargin(0.025);
    plotsStyle->SetPadBottomMargin(0.15);
    plotsStyle->SetPadLeftMargin(0.16);

    // use large Helvetica fonts
    plotsStyle->SetTextFont(42);
    plotsStyle->SetTextSize(0.08);
    plotsStyle->SetLabelFont(42,"x");
    plotsStyle->SetTitleFont(42,"x");
    plotsStyle->SetTitleFont(42,"y");
    plotsStyle->SetLabelFont(42,"y");
    plotsStyle->SetLabelFont(42,"z");
    plotsStyle->SetLabelSize(0.05,"x");
    plotsStyle->SetTitleSize(0.06,"x");
    plotsStyle->SetLabelSize(0.05,"y");
    plotsStyle->SetTitleSize(0.06,"y");
    plotsStyle->SetLabelSize(0.04,"z");
    plotsStyle->SetTitleSize(0.07,"z");
    plotsStyle->SetTitleOffset(1.3,"y");
    plotsStyle->SetTitleOffset(1.1,"x");

    plotsStyle->SetLegendBorderSize(0.0);

    // use bold lines and markers
    plotsStyle->SetMarkerStyle(20);
    plotsStyle->SetMarkerSize(MarkerSize);
    plotsStyle->SetLineWidth(LineWidth);
    plotsStyle->SetHistLineWidth(3);
    plotsStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
    //
    // get rid of X error bars and y error bar caps
    //plotsStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
    plotsStyle->SetOptTitle(0);
    plotsStyle->SetOptStat(0);
    plotsStyle->SetOptFit(0);

    // put tick marks on top and RHS of plots
    plotsStyle->SetPadTickX(1);
    plotsStyle->SetPadTickY(1);

    // set large legend fonts

    plotsStyle->cd();
    return plotsStyle;
}