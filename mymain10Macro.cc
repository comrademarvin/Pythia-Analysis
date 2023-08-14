#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TNtuple.h>

void mymain10Macro() {
    TFile *infile = TFile::Open("results/mymain10_100k.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    TH1F *HFEventsTotal = new TH1F("HF_events_full","HF Production Contributions;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 50, 0.0, 100.0);
    TH1F *HFEventsPart = new TH1F("HF_events_pt_hat_part","", 50, 0.0, 100.0);

    TH1F *HFPair = new TH1F("pair_creation","", 50, 0.0, 100.0);
    TH1F *HFExcite = new TH1F("flavour_excitation","", 50, 0.0, 100.0);
    TH1F *HFSplit = new TH1F("gluon_splitting","", 50, 0.0, 100.0);

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *HFEventsTuple = (TNtuple*)infile->Get(Form("HF_events%d", iBin));

        HFEventsPart->Reset();
        HFEventsTuple->Draw("ptHat>>HF_events_pt_hat_part");
        HFEventsPart->Scale(1/(*it), "width");
        HFEventsTotal->Add(HFEventsPart);

        HFEventsPart->Reset();
        HFEventsTuple->Draw("ptHat>>HF_events_pt_hat_part", "prodCount == 2");
        HFEventsPart->Scale(1/(*it), "width");
        HFPair->Add(HFEventsPart);

        HFEventsPart->Reset();
        HFEventsTuple->Draw("ptHat>>HF_events_pt_hat_part", "prodCount == 1");
        HFEventsPart->Scale(1/(*it), "width");
        HFExcite->Add(HFEventsPart);

        HFEventsPart->Reset();
        HFEventsTuple->Draw("ptHat>>HF_events_pt_hat_part", "prodCount == 0");
        HFEventsPart->Scale(1/(*it), "width");
        HFSplit->Add(HFEventsPart);

        iBin++;
    }

    TFile* outFile = new TFile("mymain10_production_contributions.root", "RECREATE");

    // Decay Status Contributions
    TCanvas *canvasHFEvents = new TCanvas("HF_events_sigma","HF_events_sigma");

    gPad->SetLogy();

    HFEventsTotal->SetLineColor(1);
    HFEventsTotal->Draw();

    HFPair->SetLineColor(2);
    HFPair->SetStats(0);
    HFPair->Draw("SAME");

    HFExcite->SetLineColor(3);
    HFExcite->SetStats(0);
    HFExcite->Draw("SAME");

    HFSplit->SetLineColor(4);
    HFSplit->SetStats(0);
    HFSplit->Draw("SAME");

    auto legendHFProductions = new TLegend();
    legendHFProductions->AddEntry(HFEventsTotal,"HF Events","l");
    legendHFProductions->AddEntry(HFPair,"Pair Creation Events","l");
    legendHFProductions->AddEntry(HFExcite,"Flavour Excitation","l");
    legendHFProductions->AddEntry(HFSplit,"Gluon Splitting","l");
    legendHFProductions->Draw("SAME");

    canvasHFEvents->Write();

    delete outFile;
}