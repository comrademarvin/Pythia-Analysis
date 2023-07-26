#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_contributions() {
    TFile *infile = TFile::Open("results/mymain09_test_500k.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    TFile* outFile = new TFile("mymain09Hist_contribution.root", "RECREATE");

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        TH1F *muonPtPart = new TH1F(Form("muon_pt_part_%d", iBin),"", 9, 2.0, 20.0);
        
        muonTuple->Draw(Form("pt>>muon_pt_part_%d", iBin), "(y > 2.5) && (y < 4) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonPtPart->Scale(1/(*it), "width");
        
        TCanvas *canvasMuon = new TCanvas(Form("muon_sigma_%d", iBin), Form("muon_sigma_%d", iBin));

        gPad->SetLogy();
        muonPtPart->Draw();

        canvasMuon->Write();

        iBin++;
    }

    delete outFile;
}