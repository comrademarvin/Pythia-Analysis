#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>


void mymain04Macro() {
    // Read in NTuple and info from generation output
    TFile *infile = TFile::Open("mymain04_136_100k.root", "READ");

    TNtuple *hardTuple = (TNtuple*)infile->Get("hardQCD");
    TNtuple *softTuple = (TNtuple*)infile->Get("softQCD");

    std::vector<double> *genInfoHard;
    infile->GetObject("genInfoHard", genInfoHard);

    std::vector<double> *genInfoSoft;
    infile->GetObject("genInfoSoft", genInfoSoft);

    // histograms
    const unsigned int nBins = 20;
    TH1F *hardQCDpTHat = new TH1F("hard_QCD_pTHat","Contribution to Hardest Process;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", nBins, 0, 20);
    hardTuple->Draw("pTHat>>hard_QCD_pTHat", "pTHat<20 && pTHat>2");

    TH1F *softQCDpTHat = new TH1F("soft_QCD_pTHat", "Contribution to hardest process;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", nBins, 0, 20);
    softTuple->Draw("pTHat>>soft_QCD_pTHat", "pTHat<20", "SAME");

    // normalising to cross-section
    hardQCDpTHat->Scale(1/((*genInfoHard)[0]), "width"); // normalise by pythia weightSum and bin width
    TH1F *hardQCDsigmaGen = new TH1F("hard_QCD_sigma_gen","", nBins, 0, 20);
    for (int iBin = 1; iBin <= nBins; iBin++) {
        hardQCDsigmaGen->SetBinContent(iBin, (*genInfoHard)[1]);
        hardQCDsigmaGen->SetBinError(iBin, (*genInfoHard)[2]);
    }
    hardQCDpTHat->Multiply(hardQCDsigmaGen);

    softQCDpTHat->Scale(1/((*genInfoSoft)[0]), "width"); // normalise by pythia weightSum and bin width
    TH1F *softQCDsigmaGen = new TH1F("soft_QCD_sigma_gen","", nBins, 0, 20);
    for (int iBin = 1; iBin <= nBins; iBin++) {
        softQCDsigmaGen->SetBinContent(iBin, (*genInfoSoft)[1]);
        softQCDsigmaGen->SetBinError(iBin, (*genInfoSoft)[2]);
    }
    softQCDpTHat->Multiply(softQCDsigmaGen);

    //Plotting
    TFile* outFile = new TFile("mymain04Macro.root", "RECREATE");

    hardQCDsigmaGen->Write();

    TCanvas *canvasPtHat = new TCanvas("pT_hat_contributions","pT_hat_contributions"); // contributions from soft/hard QCD to pT-hat
    gPad->SetLogy();

    hardQCDpTHat->SetMinimum(0.01);
    hardQCDpTHat->SetStats(0);
    hardQCDpTHat->SetLineColor(1);
    hardQCDpTHat->SetLineWidth(3);
    hardQCDpTHat->SetMarkerStyle(20);
    hardQCDpTHat->SetMarkerColor(1);
    hardQCDpTHat->SetMarkerSize(1.5);
    hardQCDpTHat->Draw("SAME");

    softQCDpTHat->SetStats(0);
    softQCDpTHat->SetLineColor(2);
    softQCDpTHat->SetLineWidth(3);
    softQCDpTHat->SetMarkerStyle(21);
    softQCDpTHat->SetMarkerColor(2);
    softQCDpTHat->SetMarkerSize(1.5);
    softQCDpTHat->Draw("SAME");

    auto legend = new TLegend();
    legend->AddEntry(hardQCDpTHat,"HardQCD:all","p");
    legend->AddEntry(softQCDpTHat,"SoftQCD:nonDiffractive","p");
    legend->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "This Work");
    labelCuts->DrawLatex(0.0, 0.0, "Pythia8 pp @ #sqrt{s} = 13.6 TeV");
    labelCuts->DrawLatex(0.0, 0.0, "Full Monash Tune");
    labelCuts->Draw("SAME");

    canvasPtHat->Write();

    delete outFile;
}