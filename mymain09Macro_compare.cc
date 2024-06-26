#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_compare() {
    TFile *infile = TFile::Open("mymain09.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    TH1F *muonPtTotal = new TH1F("muon_full","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", NBINS, edges);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", NBINS, edges);

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonPtPart->Reset();
        if (iBin == 0) {
            muonTuple->Draw("pt>>muon_pt_part", "(pAbs > 4.0) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        } else {
            muonTuple->Draw("pt>>muon_pt_part", "(pAbs > 4.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        }
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);

        iBin++;
    }

    // Data to compare with
    TFile *datafile = TFile::Open("Joyful Data/inclusive_muoncross_data.root", "READ");
    TH1F *muonData = new TH1F("muon_data","", NBINS, edges);
    muonData = (TH1F*)datafile->Get("crosssection_graph");

    TFile *FONLLfile = TFile::Open("Joyful Data/total_cross_section_FONNL.root", "READ");
    TH1F *muonFONLL = new TH1F("muon_fonll","", NBINS, edges);
    muonFONLL = (TH1F*)FONLLfile->Get("FONLL total");

    // Output file
    TFile* outFile = new TFile("mymain09Hist_compare.root", "RECREATE");

    // Just results
    muonPtTotal->Write("Muon_sigma_results");

    // Decay Status Contributions
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    gPad->SetLogy();

    muonPtTotal->SetLineColor(1);
    // muonPtTotal->Draw();

    muonData->SetLineColor(2);
    // muonData->Draw("SAME");

    auto rp = new TRatioPlot(muonPtTotal, muonData);
    rp->SetH1DrawOpt("E");
    rp->Draw();

    rp->GetUpperPad()->cd();
    muonFONLL->Draw("SAME");

    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtTotal,"Pythia","l");
    legendMuon->AddEntry(muonData,"Joyful Data","l");
    //legendMuon->AddEntry(muonFONLL,"FONLL","l");
    legendMuon->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.02 TeV, 2.5 < y < 4");
    labelCuts->DrawLatex(0.0, 0.0, "2 < p_{T} < 20");
    labelCuts->DrawLatex(0.0, 0.0, "|p| > 4");
    labelCuts->Draw("SAME");

    canvasMuon->Write();

    delete outFile;
}