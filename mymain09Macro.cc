#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TNtuple.h>

void mymain09Macro() {
    TFile *infile = TFile::Open("mymain09.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    TH1F *muonPtTotal = new TH1F("muon_full","Produced Muon Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", 35, 0.0, 70.0);

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);

        iBin++;
    }

    TFile* outFile = new TFile("mymain09Hist.root", "RECREATE");

    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    //muonPtTotal->SetLineColor(1);
    //muonPtTotal->SetStats(0);
    muonPtTotal->Draw();

    canvasMuon->Write();

    delete outFile;
}