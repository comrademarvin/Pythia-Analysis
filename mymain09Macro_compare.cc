#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_compare() {
    TFile *infile = TFile::Open("mymain09.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    TH1F *muonPtTotal = new TH1F("muon_full","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", 9, 2.0, 20.0);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", 9, 2.0, 20.0);

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "(y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);

        iBin++;
    }

    TFile* outFile = new TFile("mymain09Hist_compare.root", "RECREATE");

    // Decay Status Contributions
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    muonPtTotal->SetLineColor(1);
    muonPtTotal->Draw();

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.02 TeV, 2.5 < y < 4");
    labelCuts->DrawLatex(0.0, 0.0, "2 < p_{T} < 20");
    labelCuts->Draw("SAME");

    canvasMuon->Write();

    delete outFile;
}