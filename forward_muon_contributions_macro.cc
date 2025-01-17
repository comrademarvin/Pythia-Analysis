#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"


void forward_muon_contributions_macro() {
    // background from mb
    TFile* infile_mb = TFile::Open("mymain09Hist_HPC_compare_forced_20M_536.root", "READ");
    TH1F* HF_muon_pt = (TH1F*) infile_mb->Get("HF_muon_pt_contribution");

    // W
    TFile* infile_W = TFile::Open("mymain11Hist_joined_central.root", "READ");
    TH1F* W_muon_pt = (TH1F*) infile_W->Get("W_muon_pt_forward");

    // plot contributions together
    TFile* outFile = new TFile("forward_muon_contributions_hist.root", "RECREATE");

    TCanvas *canvasMuonPt = new TCanvas("W_muon_pt","W_muon_pt");
    gPad->SetLogy();

    HF_muon_pt->SetTitle("Contributions to Forward Region Semileptonic Muon Yield");
    HF_muon_pt->GetYaxis()->SetTitle("#frac{d#sigma_{#mu}}{dp_{T}} (pb/GeV/c)");
    HF_muon_pt->SetMaximum(100000000);
    HF_muon_pt->SetMinimum(0.0001);
    HF_muon_pt->SetStats(0);
    HF_muon_pt->SetLineColor(4);
    HF_muon_pt->SetLineWidth(3);
    HF_muon_pt->SetMarkerStyle(20);
    HF_muon_pt->SetMarkerColor(4);
    HF_muon_pt->SetMarkerSize(1.5);
    HF_muon_pt->Draw("SAME");

    W_muon_pt->SetLineColor(3);
    W_muon_pt->SetLineWidth(3);
    W_muon_pt->SetMarkerStyle(21);
    W_muon_pt->SetMarkerColor(3);
    W_muon_pt->SetMarkerSize(1.5);
    W_muon_pt->SetStats(0);
    W_muon_pt->Draw("SAME");

    auto legendMuon = new TLegend();
    legendMuon->AddEntry(HF_muon_pt,"c,b #rightarrow #mu (Pythia)","p");
    legendMuon->AddEntry(W_muon_pt,"W #rightarrow #mu (Powheg+Pythia)","p");
    legendMuon->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.36 TeV, 2.5 < #eta < 4");
    labelCuts->DrawLatex(0.0, 0.0, "Monash Tune (NNPDF2.3 QCD+QED LO)");
    labelCuts->Draw("SAME");

    canvasMuonPt->Write();

    delete outFile;
};