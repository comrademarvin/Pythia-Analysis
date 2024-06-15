#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_HPC_compare_forced() {
    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    TH1F *muonPt = new TH1F("muon_pt","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", 40, 0, 80);
    TH1F *muonPtCompare = new TH1F("muon_pt_compare","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", NBINS, edges);

    for(int iBin = 0; iBin < 7; iBin++) {
        TFile *infile = TFile::Open(Form("mymain09_HPC_root/mymain09_%d.root", iBin), "READ");

        std::vector<double> *binLuminocity;
        infile->GetObject("luminocity", binLuminocity);

        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        TH1F *muonPtD = new TH1F("muon_pt_D","", 40, 0, 80);
        TH1F *muonPtD0 = new TH1F("muon_pt_D0","", 40, 0, 80);
        TH1F *muonPtRest = new TH1F("muon_pt_rest","", 40, 0, 80);

        TH1F *muonPtDCompare = new TH1F("muon_pt_D_compare","", NBINS, edges);
        TH1F *muonPtD0Compare = new TH1F("muon_pt_D0_compare","", NBINS, edges);
        TH1F *muonPtRestCompare = new TH1F("muon_pt_rest_compare","", NBINS, edges);

        if (iBin == 0) {
            // pt cross-section
            muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            // compare to data
            muonTuple->Draw("pt>>muon_pt_D_compare", "(firstMother == 411) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0_compare", "(firstMother == 421) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest_compare", "(firstMother != 411) && (firstMother != 421) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        } else {
            // pt cross-section
            muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            // compare to data
            muonTuple->Draw("pt>>muon_pt_D_compare", "(firstMother == 411) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0_compare", "(firstMother == 421) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest_compare", "(firstMother != 411) && (firstMother != 421) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        }

        // Scale by Branching Ratio
        double D_BR = 0.176;
        double D0_BR = 0.067;

        muonPtD->Scale(D_BR);
        muonPtD0->Scale(D0_BR);

        muonPtDCompare->Scale(D_BR);
        muonPtD0Compare->Scale(D0_BR);

        muonPtRest->Add(muonPtD);
        muonPtRest->Add(muonPtD0);
        muonPtRest->Scale(1/(*binLuminocity)[iBin], "width");
        muonPt->Add(muonPtRest);

        muonPtRestCompare->Add(muonPtDCompare);
        muonPtRestCompare->Add(muonPtD0Compare);
        muonPtRestCompare->Scale(1/(*binLuminocity)[iBin], "width");
        muonPtCompare->Add(muonPtRestCompare);
    }

    // Data to compare with
    // TFile *datafile = TFile::Open("Joyful Data/muon_pt_differential_published.root", "READ");
    // TH1F *muonData = new TH1F("muon_data","", NBINS, edges);
    // muonData = (TH1F*)datafile->Get("crosssection_graph");

    // TFile *FONLLfile = TFile::Open("Joyful Data/total_cross_section_FONNL.root", "READ");
    // TH1F *muonFONLL = new TH1F("muon_fonll","", NBINS, edges);
    // muonFONLL = (TH1F*)FONLLfile->Get("FONLL total");

    // Output file
    TFile* outFile = new TFile("mymain09Hist_HPC_compare_forced.root", "RECREATE");

    // Pt yield results
    TCanvas *canvasMuonPt = new TCanvas("Muon_sigma_pt","Muon_sigma_pt");
    gPad->SetLogy();
    muonPt->SetMinimum(1);
    muonPt->Draw();
    canvasMuonPt->Write();

    // // Decay Status Contributions
    // TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    // gPad->SetLogy();

    // muonPtCompare->SetLineColor(1);
    // muonPtCompare->Draw();

    // muonData->SetLineColor(2);
    // // muonData->Draw("SAME");

    // auto rp = new TRatioPlot(muonPtCompare, muonData);
    // rp->SetH1DrawOpt("E");
    // rp->Draw();

    // rp->GetUpperPad()->cd();
    // //muonFONLL->Draw("SAME");

    // auto legendMuon = new TLegend();
    // legendMuon->AddEntry(muonPtCompare,"Pythia (CTEQ66.00, NLO)","l");
    // legendMuon->AddEntry(muonData,"Published","l");
    // //legendMuon->AddEntry(muonFONLL,"FONLL","l");
    // legendMuon->Draw("SAME");

    // auto labelCuts = new TLatex();
    // labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.36 TeV, 2.5 < y < 4");
    // labelCuts->DrawLatex(0.0, 0.0, "2 < p_{T} < 20");
    // labelCuts->DrawLatex(0.0, 0.0, "|p| > 4");
    // labelCuts->Draw("SAME");

    // canvasMuon->Write();

    delete outFile;
}