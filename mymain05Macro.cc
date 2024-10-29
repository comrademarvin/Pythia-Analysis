#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain05Macro() {
    TFile *infile = TFile::Open("mymain05.root", "READ");

    TNtuple *genInfo = (TNtuple*)infile->Get("genInfo");

    TFile* outFile = new TFile("mymain05Macro.root", "RECREATE");

    // Output histograms
    int N_bins_muon_pt = 20;
    float lowerPt = 0.0;
    float upperPt = 20.0;
    TH1F *muonPt = new TH1F("muon_pt","", N_bins_muon_pt, lowerPt, upperPt);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", N_bins_muon_pt, lowerPt, upperPt);

    // iterate over bins and access generated cross-section info
    float weightSum, sigmaGen, sigmaErr;
    genInfo->SetBranchAddress("weightSum",&weightSum);
    genInfo->SetBranchAddress("sigmaGen",&sigmaGen);
    genInfo->SetBranchAddress("sigmaErr",&sigmaErr);

    for(int iBin = 0; iBin < genInfo->GetEntries(); iBin++){
        genInfo->GetEntry(iBin);
        muonPtPart->Reset();
        
        // access the associated bin forward muons
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonTuple->Draw("pt>>muon_pt_part"); // pT of forward muons

        // normalise to cross-section
        muonPtPart->Scale(1/weightSum, "width");
        TH1F *sigmaGenHistPt = new TH1F("sigma_gen_pt","", N_bins_muon_pt, lowerPt, upperPt);
        for (int iBin = 0; iBin < N_bins_muon_pt; iBin++) {
            sigmaGenHistPt->SetBinContent(iBin, sigmaGen);
            sigmaGenHistPt->SetBinError(iBin, sigmaErr);
        }
        muonPtPart->Multiply(sigmaGenHistPt);
        muonPtPart->Write(Form("muon_pt_part_%d", iBin));

        muonPt->Add(muonPtPart);
    }

    TCanvas *canvasPt = new TCanvas("muon_pt","muon_pt");
    gPad->SetLogy();
    muonPt->SetMinimum(0.000000001);
    muonPt->Draw();
    canvasPt->Write();

    delete outFile;
}