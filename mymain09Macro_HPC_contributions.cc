#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_HPC_contributions() {
    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    vector<TH1F*> binHist(7);

    for(int iBin = 0; iBin < 7; iBin++) {
        TFile *infile = TFile::Open(Form("mymain09_HPC_root/mymain09_%d.root", iBin), "READ");

        std::vector<double> *binLuminocity;
        infile->GetObject("luminocity", binLuminocity);

        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        TH1F *muonPtPart = new TH1F(Form("muon_pt_part_%d", iBin),"", 40, 0, 80);

        if (iBin == 0) {
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iBin), "(pt < 8.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        } else {
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iBin), "(pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        }

        muonPtPart->Scale(1/(*binLuminocity)[iBin], "width");

        binHist[iBin] = muonPtPart;
    }

    TFile* outFile = new TFile("mymain09Hist_HPC_contribution.root", "RECREATE");

    for(int iBin = 0; iBin < 7; iBin++) {
        TCanvas *canvasMuon = new TCanvas(Form("muon_sigma_%d", iBin), Form("muon_sigma_%d", iBin));

        gPad->SetLogy();
        binHist[iBin]->Draw();

        canvasMuon->Write();
    }

    delete outFile;
}