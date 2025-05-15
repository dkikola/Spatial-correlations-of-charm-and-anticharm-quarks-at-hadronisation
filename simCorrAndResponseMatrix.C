/*
Author: Daniel Kikoła, 2024 - 2025.
This macro uses CERN ROOT and RooUnfold frameworks. 
*/


#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <TH1D.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TLatex.h>
#include "RooUnfoldResponse.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif

auto D0mass = 1.869; // D0 mass
TF1 fBoltzman("fBoltzman", "x*x*exp(-TMath::Sqrt(1.869*1.869+x*x)/[0])", 0.0, 10);
TF1 fGammaWeight("fGammaWeight", "1./TMath::Sqrt(1+x/([0]))", 0, 10); // 1/gamma factor vs x^2, par0 is t_fo^2

TRandom3 *randGen;

auto tauFreezeOut = 9.0;  // freeze-out time
auto Rmax = 6.0;          // max distance from the c.m. of the expanding system; assumed 6 fm for Pb+Pb at SPS energies
auto Rmax2 = Rmax * Rmax; // max distance from the c.m. of the expanding system

auto c = 1.0;
auto sigmaBrownianMotion = 0.0;

void getThermalProduction(TLorentzVector &DMesonMom, TVector3 &fluidVelocity)
{

    double phiC = randGen->Uniform(0, TMath::TwoPi());
    double cosTheta = randGen->Uniform(-1., 1.);
    double thetaC = TMath::ACos(cosTheta);
    double p = fBoltzman.GetRandom();

    double px = p * TMath::Sin(thetaC) * TMath::Cos(phiC);
    double py = p * TMath::Sin(thetaC) * TMath::Sin(phiC);
    double pz = p * TMath::Cos(thetaC);

    DMesonMom.SetXYZM(px, py, pz, D0mass);
    DMesonMom.Boost(fluidVelocity);
}

double gamma(double v)
{
    return 1. / TMath::Sqrt(1 - v * v);
}

double uToV(double u)
{
    return u / TMath::Sqrt(1 + u * u);
}

double getDeltaPhi(double phiC, double phiCbar)
{
    double dphi = phiC - phiCbar;
    if (dphi < 0)
        dphi = dphi + TMath::TwoPi();
    if (dphi > TMath::Pi())
        dphi = TMath::TwoPi() - dphi;
    return dphi;
}

void addBrownianMotion3D(double &x, double &y, double &z)
{

    bool accepted = false;

    while (accepted == false)
    {
        double dx = randGen->Gaus(0, sigmaBrownianMotion);
        double dy = randGen->Gaus(0, sigmaBrownianMotion);
        double dz = randGen->Gaus(0, sigmaBrownianMotion);

        if (((x + dx) * (x + dx) + (y + dy) * (y + dy) + (z + dz) * (z + dz)) < Rmax2)
        {
            x = x + dx;
            y = y + dy;
            z = z + dz;
            return;
        }
    }
}

bool acceptPointWithLorentzContaction(double rSquare)
{

    float test = randGen->Uniform(0, 1.);
    if (test < fGammaWeight.Eval(rSquare))
        return true;
    else
        return false;
}

TStyle *setWindowDressing();

void simCorrAndResponseMatrix(TString outFileName = "Test-Results-dEta-dPhi-corr-T-150MeV.root", double sigmaBrownian = 2.0, int nMax = 10000000, double Temp = 150 /*MeV*/)
{

    setWindowDressing();

    fBoltzman.SetParameter(0, Temp / 1000.);                   // freeze-out temperature in GeV
    fGammaWeight.SetParameter(0, tauFreezeOut * tauFreezeOut); // weight factor to account for the Lorentz contraction

    sigmaBrownianMotion = sigmaBrownian;

    TH1D *hDeltaPhiNonLocalWithHadr = new TH1D("hDeltaPhiNonLocalWithHadr", "dNdPhi;#Delta#phi;dN/d#Delta#phi", 100, 0, TMath::Pi());
    TH1D *hDeltaPhiLocalWithHadr = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiLocalWithHadr");

    TH1D *hDeltaPhiNonLocalFlowOnly = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiNonLocalFlowOnly");
    TH1D *hDeltaPhiLocalFlowOnly = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiLocalFlowOnly");
    TH1D *hDeltaPhiLocalFlowAndSmearing = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiLocalFlowAndSmearing");

    TH1D *hDeltaPhiLocalWithGaussSmearingAndHadr = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiLocalWithGaussSmearingAndHadr");
    TH1D *hDeltaPhiNonLocalWithGaussSmearingAndHadr = (TH1D *)hDeltaPhiNonLocalWithHadr->Clone("hDeltaPhiNonLocalWithGaussSmearingAndHadr");

    RooUnfoldResponse response(hDeltaPhiLocalWithHadr->GetNbinsX(), 0, TMath::Pi());

    randGen = new TRandom3(0);

    TVector3 fluidVelocity, fluidVelocityC, fluidVelocityCbar;

    TLorentzVector DMesonMom, DbarMesonMom;
    DMesonMom.SetPtEtaPhiM(1, 0, 0, D0mass);
    DbarMesonMom.SetPtEtaPhiM(1, 0, 0, D0mass);

    TLatex label;
    label.SetTextSize(0.055);

    TString txt;
    txt.Form("#sigma_{Gaussian smearing} = %1.2lf fm", sigmaBrownian);

    // local production (same production point of c and c-bar quarks)
    auto nPts = 0;
    while (nPts < nMax)
    {

        double x = randGen->Uniform(-Rmax, Rmax);
        double y = randGen->Uniform(-Rmax, Rmax);
        double z = randGen->Uniform(-Rmax, Rmax);
        double r2 = (x * x + y * y + z * z);

        // model the density change due to Lorentz contraction
        if (acceptPointWithLorentzContaction(r2) == false)
            continue;

        if (r2 < Rmax2)
        {

            auto r = TMath::Sqrt(r2);
            auto u = r / tauFreezeOut;

            double uToVConvFact = 1. / TMath::Sqrt(1 + u * u);

            double vx = x / tauFreezeOut * uToVConvFact;
            double vy = y / tauFreezeOut * uToVConvFact;
            double vz = z / tauFreezeOut * uToVConvFact;

            double v = u * uToVConvFact;

            nPts++;

            fluidVelocity.SetXYZ(vx, vy, vz);

            getThermalProduction(DMesonMom, fluidVelocity);
            getThermalProduction(DbarMesonMom, fluidVelocity);

            double dphi = getDeltaPhi(DMesonMom.Phi(), DbarMesonMom.Phi());

            hDeltaPhiLocalWithHadr->Fill(dphi);
            hDeltaPhiLocalFlowOnly->Fill(0);

            double xCBrownian, yCBrownian, zCBrownian;
            double xCbarBrownian, yCbarBrownian, zCbarBrownian;

            xCBrownian = x;
            yCBrownian = y;
            zCBrownian = z;

            xCbarBrownian = x;
            yCbarBrownian = y;
            zCbarBrownian = z;

            addBrownianMotion3D(xCBrownian, yCBrownian, zCBrownian);
            addBrownianMotion3D(xCbarBrownian, yCbarBrownian, zCbarBrownian);

            double rCBrownian = TMath::Sqrt(xCBrownian * xCBrownian + yCBrownian * yCBrownian + zCBrownian * zCBrownian);
            double rCbarBrownian = TMath::Sqrt(xCbarBrownian * xCbarBrownian + yCbarBrownian * yCbarBrownian + zCbarBrownian * zCbarBrownian);

            double uC = rCBrownian / tauFreezeOut;
            double uCbar = rCbarBrownian / tauFreezeOut;

            double uToVConvFactC = 1. / TMath::Sqrt(1 + uC * uC);
            double uToVConvFactCbar = 1. / TMath::Sqrt(1 + uCbar * uCbar);

            fluidVelocityC.SetXYZ(xCBrownian / tauFreezeOut * uToVConvFactC,
                                  yCBrownian / tauFreezeOut * uToVConvFactC,
                                  zCBrownian / tauFreezeOut * uToVConvFactC);

            fluidVelocityCbar.SetXYZ(xCbarBrownian / tauFreezeOut * uToVConvFactCbar,
                                     yCbarBrownian / tauFreezeOut * uToVConvFactCbar,
                                     zCbarBrownian / tauFreezeOut * uToVConvFactCbar);

            dphi = getDeltaPhi(fluidVelocityC.Phi(), fluidVelocityCbar.Phi());

            hDeltaPhiLocalFlowAndSmearing->Fill(dphi);

            float dphiBeforeHadronisation = dphi;

            getThermalProduction(DMesonMom, fluidVelocityC);
            getThermalProduction(DbarMesonMom, fluidVelocityCbar);

            dphi = getDeltaPhi(DMesonMom.Phi(), DbarMesonMom.Phi());
            response.Fill(dphi, dphiBeforeHadronisation);

            hDeltaPhiLocalWithGaussSmearingAndHadr->Fill(dphi);
        }
    }

    nPts = 0;

    // non-local production
    while (nPts < nMax)
    {
        double xC = randGen->Uniform(-Rmax, Rmax);
        double yC = randGen->Uniform(-Rmax, Rmax);
        double zC = randGen->Uniform(-Rmax, Rmax);

        double rC2 = (xC * xC + yC * yC + zC * zC);

        if (acceptPointWithLorentzContaction(rC2) == false)
            continue;

        double xCbar = randGen->Uniform(-Rmax, Rmax);
        double yCbar = randGen->Uniform(-Rmax, Rmax);
        double zCbar = randGen->Uniform(-Rmax, Rmax);

        double rCbar2 = (xCbar * xCbar + yCbar * yCbar + zCbar * zCbar);

        if (acceptPointWithLorentzContaction(rCbar2) == false)
            continue;

        if ((rC2 < Rmax2) && (rCbar2 < Rmax2))
        {
            nPts++;

            auto rC = TMath::Sqrt(rC2);
            auto rCbar = TMath::Sqrt(rCbar2);

            auto uC = rC / tauFreezeOut;
            auto uCbar = rCbar / tauFreezeOut;

            double uToVConvFactC = 1. / TMath::Sqrt(1 + uC * uC);
            double uToVConvFactCbar = 1. / TMath::Sqrt(1 + uCbar * uCbar);

            fluidVelocityC.SetXYZ(xC / tauFreezeOut * uToVConvFactC, yC / tauFreezeOut * uToVConvFactC, zC / tauFreezeOut * uToVConvFactC);
            fluidVelocityCbar.SetXYZ(xCbar / tauFreezeOut * uToVConvFactCbar, yCbar / tauFreezeOut * uToVConvFactCbar, zCbar / tauFreezeOut * uToVConvFactCbar);

            double vC = fluidVelocityC.Mag();
            double vCbar = fluidVelocityCbar.Mag();

            double dphi = getDeltaPhi(fluidVelocityC.Phi(), fluidVelocityCbar.Phi());

            hDeltaPhiNonLocalFlowOnly->Fill(dphi);

            getThermalProduction(DMesonMom, fluidVelocityC);
            getThermalProduction(DbarMesonMom, fluidVelocityCbar);

            dphi = getDeltaPhi(DMesonMom.Phi(), DbarMesonMom.Phi());

            hDeltaPhiNonLocalWithHadr->Fill(dphi);

            double xCBrownian, yCBrownian, zCBrownian;
            double xCbarBrownian, yCbarBrownian, zCbarBrownian;

            xCBrownian = xC;
            yCBrownian = yC;
            zCBrownian = zC;

            xCbarBrownian = xCbar;
            yCbarBrownian = yCbar;
            zCbarBrownian = zCbar;

            addBrownianMotion3D(xCBrownian, yCBrownian, zCBrownian);
            addBrownianMotion3D(xCbarBrownian, yCbarBrownian, zCbarBrownian);

            // recalculate u and v for c and c-bar

            rC = TMath::Sqrt(xCBrownian * xCBrownian + yCBrownian * yCBrownian + zCBrownian * zCBrownian);
            rCbar = TMath::Sqrt(xCbarBrownian * xCbarBrownian + yCbarBrownian * yCbarBrownian + zCbarBrownian * zCbarBrownian);

            uC = rC / tauFreezeOut;
            uCbar = rCbar / tauFreezeOut;

            uToVConvFactC = 1. / TMath::Sqrt(1 + uC * uC);
            uToVConvFactCbar = 1. / TMath::Sqrt(1 + uCbar * uCbar);

            fluidVelocityC.SetXYZ(xCBrownian / tauFreezeOut * uToVConvFactC, yCBrownian / tauFreezeOut * uToVConvFactC, zCBrownian / tauFreezeOut * uToVConvFactC);
            fluidVelocityCbar.SetXYZ(xCbarBrownian / tauFreezeOut * uToVConvFactCbar, yCbarBrownian / tauFreezeOut * uToVConvFactCbar, zCbarBrownian / tauFreezeOut * uToVConvFactCbar);

            getThermalProduction(DMesonMom, fluidVelocityC);
            getThermalProduction(DbarMesonMom, fluidVelocityCbar);

            dphi = getDeltaPhi(DMesonMom.Phi(), DbarMesonMom.Phi());
            hDeltaPhiNonLocalWithGaussSmearingAndHadr->Fill(dphi);
        }
    }

    TCanvas *c1 = new TCanvas("Corr", "Corr", 1000, 600);
    hDeltaPhiLocalWithHadr->SetLineColor(kRed);
    hDeltaPhiLocalFlowOnly->SetLineColor(kGreen + 2);
    hDeltaPhiNonLocalFlowOnly->SetLineColor(kBlue);
    hDeltaPhiLocalWithGaussSmearingAndHadr->SetLineColor(kOrange);
    hDeltaPhiNonLocalWithGaussSmearingAndHadr->SetLineColor(kCyan - 6);
    hDeltaPhiLocalWithHadr->Draw();
    hDeltaPhiNonLocalWithHadr->Draw("same");
    hDeltaPhiLocalFlowOnly->Draw("same");
    hDeltaPhiNonLocalFlowOnly->Draw("same");
    hDeltaPhiLocalWithGaussSmearingAndHadr->Draw("same");
    hDeltaPhiNonLocalWithGaussSmearingAndHadr->Draw("same");

    auto legend = new TLegend(0.4, 0.6, 0.9, 0.9);
    legend->AddEntry(hDeltaPhiLocalFlowOnly, "Local prod.", "l");
    legend->AddEntry(hDeltaPhiLocalWithHadr, "Local + Hadr.", "l");
    legend->AddEntry(hDeltaPhiLocalWithGaussSmearingAndHadr, "Local + Smearing + Hadr.", "l");
    legend->AddEntry(hDeltaPhiNonLocalFlowOnly, "Non-local prod.", "l");
    legend->AddEntry(hDeltaPhiNonLocalWithHadr, "Non-local + Hadr.", "l");
    legend->AddEntry(hDeltaPhiNonLocalWithGaussSmearingAndHadr, "Non-local + Smearing + Hadr.", "l");
    legend->Draw();
    label.DrawLatexNDC(0.5, 0.52, txt);

    TFile *fOut = new TFile(outFileName, "recreate");
    c1->Write();

    hDeltaPhiLocalWithHadr->Write();
    hDeltaPhiNonLocalWithHadr->Write();
    hDeltaPhiLocalFlowOnly->Write();
    hDeltaPhiNonLocalFlowOnly->Write();
    hDeltaPhiLocalWithGaussSmearingAndHadr->Write();
    hDeltaPhiNonLocalWithGaussSmearingAndHadr->Write();
    hDeltaPhiLocalFlowAndSmearing->Write();

    TH1D *hTrue = (TH1D *)response.Htruth();
    hTrue->Write();

    TH1D *hMeasured = (TH1D *)response.Hmeasured();
    hMeasured->Write();

    TH2D *hResponse = (TH2D *)response.Hresponse();
    hResponse->Write();

    fOut->Close();
}

float LegFontSize = 0.09;
float LegMargin = 0.25;
float MarkerSize = 2;
float LineWidth = 2.0;

TStyle *setWindowDressing()
{

    TStyle *plotsStyle = new TStyle("plotsStyle", "plain plots style");
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

    // set the paper & margin sizes
    plotsStyle->SetPaperSize(20, 26);
    plotsStyle->SetPadTopMargin(0.05);
    plotsStyle->SetPadRightMargin(0.03);
    plotsStyle->SetPadBottomMargin(0.16);
    plotsStyle->SetPadLeftMargin(0.12);
    //
    // use large Times-Roman fonts
    plotsStyle->SetTextFont(132);
    plotsStyle->SetTextSize(0.08);
    plotsStyle->SetLabelFont(132, "x");
    plotsStyle->SetTitleFont(132, "x");
    plotsStyle->SetTitleFont(132, "y");
    plotsStyle->SetLabelFont(132, "y");
    plotsStyle->SetLabelFont(132, "z");
    plotsStyle->SetLabelSize(0.06, "x");
    plotsStyle->SetTitleSize(0.07, "x");
    plotsStyle->SetLabelSize(0.06, "y");
    plotsStyle->SetTitleSize(0.07, "y");
    plotsStyle->SetLabelSize(0.06, "z");
    plotsStyle->SetTitleSize(0.07, "z");
    plotsStyle->SetTitleOffset(0.85, "y");
    plotsStyle->SetTitleOffset(0.9, "x");

    plotsStyle->SetLegendBorderSize(0.0);

    // use bold lines and markers
    plotsStyle->SetMarkerStyle(20);
    plotsStyle->SetMarkerSize(MarkerSize);
    plotsStyle->SetLineWidth(LineWidth);
    plotsStyle->SetHistLineWidth(LineWidth);
    plotsStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

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
