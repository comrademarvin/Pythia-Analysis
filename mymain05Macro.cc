#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain05Macro() {
    TFile *infile = TFile::Open("mymain05_1M.root", "READ");

    TNtuple *genInfo = (TNtuple*)infile->Get("genInfo");

    TFile* outFile = new TFile("mymain05Macro.root", "RECREATE");

    // Output histograms
    int N_bins_muon_pt = 20;
    float lowerPt = 0.0;
    float upperPt = 20.0;
    const unsigned int N_cont = 8;
    const string contrib[N_cont] = {"total", "lepton", "boson", "light meson", "charm meson", "bottom meson", "charmonium", "baryon"};
    const unsigned int contrib_lower_pdg[N_cont] = {0, 10, 20, 110, 410, 510, 440, 1110};
    const unsigned int contrib_upper_pdg[N_cont] = {10000000, 19, 38, 340, 440, 550, 450, 5555};

    vector<TH1F*> muonPt(N_cont);
    vector<TH1F*> muonPtPart(N_cont);
    for (int iCont = 0; iCont < N_cont; iCont++) {
        muonPt[iCont] = new TH1F(Form("muon_pt_%d", iCont),"Contributions to Forward Muon Production;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_muon_pt, lowerPt, upperPt);
        muonPtPart[iCont] = new TH1F(Form("muon_pt_part_%d", iCont),"", N_bins_muon_pt, lowerPt, upperPt);
    }

    // iterate over bins and access generated cross-section info
    float weightSum, sigmaGen, sigmaErr;
    genInfo->SetBranchAddress("weightSum",&weightSum);
    genInfo->SetBranchAddress("sigmaGen",&sigmaGen);
    genInfo->SetBranchAddress("sigmaErr",&sigmaErr);

    for(int iBin = 0; iBin < genInfo->GetEntries(); iBin++){
        genInfo->GetEntry(iBin);
        for (int iCont = 0; iCont < N_cont; iCont++) muonPtPart[iCont]->Reset();
        
        // access the associated bin forward muons
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        for (int iCont = 0; iCont < N_cont; iCont++) {
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iCont), 
                            Form("firstMother > %d && firstMother < %d", contrib_lower_pdg[iCont], contrib_upper_pdg[iCont]));
        }

        // normalise to cross-section
        TH1F *sigmaGenHistPt = new TH1F("sigma_gen_pt","", N_bins_muon_pt, lowerPt, upperPt);
        for (int iBin = 1; iBin <= N_bins_muon_pt; iBin++) {
            sigmaGenHistPt->SetBinContent(iBin, sigmaGen);
            sigmaGenHistPt->SetBinError(iBin, sigmaErr);
        }
        for (int iCont = 0; iCont < N_cont; iCont++) {
            muonPtPart[iCont]->Scale(1/weightSum, "width");
            muonPtPart[iCont]->Multiply(sigmaGenHistPt);
            if (iCont == 0) muonPtPart[iCont]->Write(Form("muon_pt_part_%d", iBin));
            muonPt[iCont]->Add(muonPtPart[iCont]);
        }
    }

    int lineColours[N_cont] = {1,2,3,4,6,7,8,9};
    TCanvas *canvasPt = new TCanvas("muon_pt","muon_pt");
    gPad->SetLogy();
    auto legend = new TLegend();
    for (int iCont = 0; iCont < N_cont; iCont++) {
        muonPt[iCont]->SetMinimum(0.000000001);
        muonPt[iCont]->SetStats(0);
        muonPt[iCont]->SetLineColor(lineColours[iCont]);
        muonPt[iCont]->SetMarkerStyle(20+iCont);
        muonPt[iCont]->SetMarkerColor(lineColours[iCont]);
        muonPt[iCont]->Draw("SAME");
        legend->AddEntry(muonPt[iCont],(contrib[iCont]).c_str(),"p");
    }
    legend->Draw("SAME");
    canvasPt->Write();

    // also plot total cross-section
    TH1F* total_cs_pt_hat = (TH1F*)infile->Get("pT_hat");
    TCanvas *canvasPtHat = new TCanvas("total_cs_pt_hat","total_cs_pt_hat");
    gPad->SetLogy();

    total_cs_pt_hat->SetMaximum(100);
    total_cs_pt_hat->SetMinimum(0.0000001);
    total_cs_pt_hat->SetStats(0);
    total_cs_pt_hat->SetLineColor(3);
    total_cs_pt_hat->SetLineWidth(3);
    total_cs_pt_hat->SetMarkerStyle(20);
    total_cs_pt_hat->SetMarkerColor(3);
    total_cs_pt_hat->SetMarkerSize(1.3);
    total_cs_pt_hat->Draw("SAME");

    auto labelPtHat = new TLatex();
    labelPtHat->DrawLatex(0.0, 0.0, "Pythia8 pp @ #sqrt{s} = 5.36 TeV, Monash Tune");
    labelPtHat->DrawLatex(0.0, 0.0, "Minimum Bias (QCD)");
    labelPtHat->Draw("SAME");

    canvasPtHat->Write();

    delete outFile;
}