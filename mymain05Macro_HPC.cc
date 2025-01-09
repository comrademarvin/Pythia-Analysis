#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain05Macro_HPC() {
    // pT-hat bin used to generate events
    const int nBins = 7;
    static const double tempArray[8] = {0.0, 10.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0};

    // Output histograms
    int N_bins_muon_pt = 20;
    float lowerPt = 0.0;
    float upperPt = 20.0;
    const unsigned int N_cont = 8;
    const string contrib[N_cont] = {"total", "lepton", "boson", "light meson", "charm meson", "bottom meson", "charmonium", "baryon"};
    const unsigned int contrib_lower_pdg[N_cont] = {0, 10, 20, 110, 410, 510, 440, 1110};
    const unsigned int contrib_upper_pdg[N_cont] = {10000000, 19, 38, 340, 440, 550, 450, 5555};

    vector<TH1F*> muonPt(N_cont);
    for (int iCont = 0; iCont < N_cont; iCont++) {
        muonPt[iCont] = new TH1F(Form("muon_pt_%d", iCont),"Contributions to Forward Muon Production;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_muon_pt, lowerPt, upperPt);
    }

    // read in root output per pT-hat bin
    for(int iBin = 0; iBin < nBins; iBin++) {
        TFile *infile = TFile::Open(Form("mymain05_HPC_root_1M_536/mymain05_%d.root", iBin), "READ");

        std::vector<double> *genInfoNorm;
        infile->GetObject("genInfoNorm", genInfoNorm);

        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        vector<TH1F*> muonPtPart(N_cont);

        for (int iCont = 0; iCont < N_cont; iCont++) {
            muonPtPart[iCont] = new TH1F(Form("muon_pt_part_%d", iCont),"", N_bins_muon_pt, lowerPt, upperPt);
            muonTuple->Draw(Form("pt>>muon_pt_part_%d", iCont), 
                            Form("firstMother > %d && firstMother < %d", contrib_lower_pdg[iCont], contrib_upper_pdg[iCont]));
        }

        // normalise
        for (int iCont = 0; iCont < N_cont; iCont++) {
            muonPtPart[iCont]->Scale((*genInfoNorm)[iBin], "width");
            muonPt[iCont]->Add(muonPtPart[iCont]);
        }
    };

    // output root file
    TFile* outFile = new TFile("mymain05Macro_HPC.root", "RECREATE");

    // plotting
    int lineColours[N_cont] = {1,2,3,4,6,28,8,9};
    int markerStyles[N_cont] = {20,21,22,23,29,33,34,47};
    TCanvas *canvasPt = new TCanvas("muon_pt","muon_pt");
    gPad->SetLogy();
    auto legend = new TLegend();
    for (int iCont = 0; iCont < N_cont; iCont++) {
        muonPt[iCont]->GetYaxis()->SetTitle("#frac{d#sigma_{X#rightarrow#mu}}{dp_{T}} (mb/GeV/c)");
        muonPt[iCont]->SetMaximum(0);
        muonPt[iCont]->SetMinimum(0.0000000001);
        muonPt[iCont]->SetStats(0);
        muonPt[iCont]->SetLineColor(lineColours[iCont]);
        muonPt[iCont]->SetLineWidth(3);
        muonPt[iCont]->SetMarkerStyle(markerStyles[iCont]);
        muonPt[iCont]->SetMarkerColor(lineColours[iCont]);
        muonPt[iCont]->SetMarkerSize(1.7);
        muonPt[iCont]->Draw("SAME");
        legend->AddEntry(muonPt[iCont],(contrib[iCont]).c_str(),"p");
    }
    legend->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "Pythia8 pp @ #sqrt{s} = 5.36 TeV");
    labelCuts->DrawLatex(0.0, 0.0, "Monash Tune, Minimum Bias (QCD)");
    labelCuts->DrawLatex(0.0, 0.0, "2.5 < #eta_{muon} < 4");
    labelCuts->Draw("SAME");

    canvasPt->Write();

    delete outFile;
}