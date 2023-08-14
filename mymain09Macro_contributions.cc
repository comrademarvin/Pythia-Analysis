#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_contributions() {
    TFile *infile = TFile::Open("mymain09.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    TFile* outFile = new TFile("mymain09Hist_contribution.root", "RECREATE");

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        TH1F *muonPtPart = new TH1F(Form("muon_pt_part_%d", iBin),"", NBINS, edges);

        if (iBin == 0) {
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iBin), "(pt < 8.0) && (pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        } else {
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iBin), "(pt > 2.0) && (y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        }

        muonPtPart->Scale(1/(*it), "width");
        
        TCanvas *canvasMuon = new TCanvas(Form("muon_sigma_%d", iBin), Form("muon_sigma_%d", iBin));

        gPad->SetLogy();
        muonPtPart->Draw();

        canvasMuon->Write();

        iBin++;
    }

    delete outFile;
}