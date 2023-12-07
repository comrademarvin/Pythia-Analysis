#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_HPC_compare_PDF() {
    int pdf_count = 4;
    int pdf_runs[] = {8, 9, 13, 19};
    const char* pdf_labels[] = {"CTEQ6L1, LO", "CTEQ66.00, NLO", "NNPDF2.3 QCD+QED LO (Monash)", "NNPDF3.1 QCD+LUXQED NLO"};

    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    vector<TH1F*> muonPtTotal(pdf_count);

    for (int iPDF = 0; iPDF < pdf_count; iPDF++) {
        muonPtTotal[iPDF] = new TH1F(Form("muon_full_pdf_%d", pdf_runs[iPDF]),"HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", NBINS, edges);

        for(int iBin = 0; iBin < 7; iBin++) {
            TFile *infile = TFile::Open(Form("HPC_4M_PDF%d/mymain09_%d.root", pdf_runs[iPDF], iBin), "READ");

            std::vector<double> *binLuminocity;
            infile->GetObject("luminocity", binLuminocity);

            TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

            TH1F *muonPtD = new TH1F("muon_pt_D","", NBINS, edges);
            TH1F *muonPtD0 = new TH1F("muon_pt_D0","", NBINS, edges);
            TH1F *muonPtRest = new TH1F("muon_pt_rest","", NBINS, edges);

            if (iBin == 0) {
                muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
                muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
                muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (pt < 10.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            } else {
                muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
                muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
                muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            }

            // Scale by Branching Ratio
            double D_BR = 0.176;
            double D0_BR = 0.067;
            muonPtD->Scale(D_BR);
            muonPtD0->Scale(D0_BR);

            muonPtRest->Add(muonPtD);
            muonPtRest->Add(muonPtD0);
            muonPtRest->Scale(1/(*binLuminocity)[iBin], "width");
            muonPtTotal[iPDF]->Add(muonPtRest);
        }
    }

    // Output file
    TFile* outFile = new TFile("mymain09Hist_HPC_compare_PDF.root", "RECREATE");

    // Decay Status Contributions
    TCanvas *canvasMuon = new TCanvas("Muon_sigma_PDF","Muon_sigma_PDF");

    gPad->SetLogy();

    auto legendMuon = new TLegend();

    for (int iPDF = 0; iPDF < pdf_count; iPDF++) {
        muonPtTotal[iPDF]->SetLineColor(iPDF+1);
        muonPtTotal[iPDF]->SetStats(0);
        muonPtTotal[iPDF]->Draw("SAME");

        legendMuon->AddEntry(muonPtTotal[iPDF],pdf_labels[iPDF],"l");
    }

    legendMuon->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.02 TeV, 2.5 < y < 4");
    labelCuts->DrawLatex(0.0, 0.0, "2 < p_{T} < 20");
    labelCuts->DrawLatex(0.0, 0.0, "|p| > 4");
    labelCuts->Draw("SAME");

    canvasMuon->Write();

    delete outFile;
}