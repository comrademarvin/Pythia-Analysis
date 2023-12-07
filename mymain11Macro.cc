#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain11Macro() {
    TFile *infile_W_plus = TFile::Open("/mnt/d/Pythia_Results/mymain11_W+_100k.root", "READ");
    TFile *infile_W_minus = TFile::Open("/mnt/d/Pythia_Results/mymain11_W-_100k.root", "READ");

    TFile* outFile = new TFile("mymain11Macro.root", "RECREATE");

    // Ratio plot of pT cross-section
    TH1F* W_plus_pt = (TH1F*)infile_W_plus->Get("W_pt_sigma");
    TH1F* W_minus_pt = (TH1F*)infile_W_minus->Get("W_pt_sigma");

    TCanvas *canvasWPt = new TCanvas("W_pt_sigma","W_pt_sigma");

    gPad->SetLogy();

    W_plus_pt->SetLineColor(1);
    W_minus_pt->SetLineColor(2);

    auto rp1 = new TRatioPlot(W_plus_pt, W_minus_pt);
    rp1->SetH1DrawOpt("E");
    rp1->Draw();

    auto legendW = new TLegend();
    legendW->AddEntry(W_plus_pt,"W+ -> #mu","l");
    legendW->AddEntry(W_minus_pt,"W- -> #mu","l");
    //legendW->AddEntry(muonFONLL,"FONLL","l");
    legendW->Draw("SAME");

    canvasWPt->Write();

    // eta distributions
    TH1F* W_plus_eta = (TH1F*)infile_W_plus->Get("W_eta_sigma");
    TH1F* W_minus_eta = (TH1F*)infile_W_minus->Get("W_eta_sigma");

    TCanvas *canvasWEta = new TCanvas("W_eta_sigma","W_eta_sigma");

    gPad->SetLogy();

    W_plus_eta->SetLineColor(1);
    W_minus_eta->SetLineColor(2);

    auto rp2 = new TRatioPlot(W_plus_eta, W_minus_eta);
    rp2->SetH1DrawOpt("E");
    rp2->Draw();

    auto legendWEta = new TLegend();
    legendWEta->AddEntry(W_plus_eta,"W+ -> #mu","l");
    legendWEta->AddEntry(W_minus_eta,"W- -> #mu","l");
    //legendWEta->AddEntry(muonFONLL,"FONLL","l");
    legendWEta->Draw("SAME");

    canvasWEta->Write();

    delete outFile;
}