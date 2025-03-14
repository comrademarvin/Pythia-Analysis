#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_HPC_compare_forced() {
    // generated pT-hat bins info
    const int nBins = 7;
    static const float pT_hat_binEdges[nBins+1] = {0.0, 10.0, 40.0, 70.0, 100.0, 150.0, 200.0, 300.0};

    // published data bins
    const Int_t NBINS_PUB = 26;
    Double_t edges_pub[NBINS_PUB + 1] = {2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16.5, 18, 20};

    // different binning for comparison with data and W->mu cross-section
    TH1F *muonPt = new TH1F("muon_pt","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b#rightarrow#mu}}{dp_{T}} (pb/GeV/c)", 40, 0, 80);
    TH1F *muonPtCompare = new TH1F("muon_pt_compare","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b#rightarrow#mu}}{dp_{T}} (pb/GeV/c)", NBINS_PUB, edges_pub);

    vector<TH1F*> muonPtContribs(nBins);

    for(int iBin = 0; iBin < nBins; iBin++) {
        TFile *infile = TFile::Open(Form("mymain09_HPC_root_20M_502_tune2/mymain09_%d.root", iBin), "READ");

        std::vector<double> *genInfoNorm;
        infile->GetObject("genInfoNorm", genInfoNorm);

        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        TH1F *muonPtD = new TH1F("muon_pt_D","", 40, 0, 80);
        TH1F *muonPtD0 = new TH1F("muon_pt_D0","", 40, 0, 80);
        TH1F *muonPtRest = new TH1F("muon_pt_rest","", 40, 0, 80);

        TH1F *muonPtDCompare = new TH1F("muon_pt_D_compare","", NBINS_PUB, edges_pub);
        TH1F *muonPtD0Compare = new TH1F("muon_pt_D0_compare","", NBINS_PUB, edges_pub);
        TH1F *muonPtRestCompare = new TH1F("muon_pt_rest_compare","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", NBINS_PUB, edges_pub);

        // pt cross-section
        muonTuple->Draw("pt>>muon_pt_D", "(pAbs > 4) && (firstMother == 411) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonTuple->Draw("pt>>muon_pt_D0", "(pAbs > 4) && (firstMother == 421) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonTuple->Draw("pt>>muon_pt_rest", "(pAbs > 4) && (firstMother != 411) && (firstMother != 421) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");

        muonTuple->Draw("pt>>muon_pt_D_compare", "(pAbs > 4) && (firstMother == 411) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonTuple->Draw("pt>>muon_pt_D0_compare", "(pAbs > 4) && (firstMother == 421) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        muonTuple->Draw("pt>>muon_pt_rest_compare", "(pAbs > 4) && (firstMother != 411) && (firstMother != 421) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");

        // Scale by Branching Ratio
        double D_BR = 0.176;
        double D0_BR = 0.067;

        muonPtD->Scale(D_BR);
        muonPtD0->Scale(D0_BR);

        muonPtDCompare->Scale(D_BR);
        muonPtD0Compare->Scale(D0_BR);

        muonPtRest->Add(muonPtD);
        muonPtRest->Add(muonPtD0);
        muonPtRest->Scale((*genInfoNorm)[iBin]*pow(10,9), "width");  // with mb to pb conversion
        muonPt->Add(muonPtRest);

        muonPtRestCompare->Add(muonPtDCompare);
        muonPtRestCompare->Add(muonPtD0Compare);
        muonPtRestCompare->Scale((*genInfoNorm)[iBin]*pow(10,9), "width");
        muonPtCompare->Add(muonPtRestCompare);

        muonPtContribs[iBin] = new TH1F(*muonPtRestCompare);
    }

    // Data to compare with
    TH1F* pub_data_hist = new TH1F("pub_data_hist", "", NBINS_PUB, edges_pub);
    float pub_data[NBINS_PUB] = {5010950, 2422650, 1280790, 719302, 415048, 251516, 164166, 105080, 67632.7, 50979, 35201.8, 23867.8, 19428.1, 12946.5, 9452.1, 7511.8, 5730.3, 3931.54, 3296.76, 2725.1, 1819.58, 1194.78, 693.276, 557.351, 333.684, 163.149};
    float pub_data_stat_err[NBINS_PUB] = {10816.73689, 7987.452823, 6046.238161, 4660.875555, 3583.329359, 2824.952257, 2327.709714, 1858.728596, 1501.46623, 1319.351814, 876.764192, 720.296789, 649.048136, 522.053371, 442.150334, 394.811194, 345.55027, 283.734524, 259.820952, 236.296146, 136.383344, 109.068001, 80.05258, 59.754158, 44.993951, 26.821532};
    float pub_data_sys_err[NBINS_PUB] = {387765.3504, 138291.1609, 54371.84092, 23641.58656, 12017.09227, 6622.542038, 4047.200815, 2431.172912, 1473.202524, 1040.002187, 1190.584719, 813.722519, 653.168836, 444.356246, 330.711965, 266.299319, 208.36918, 151.560867, 133.733729, 117.958406, 88.561688, 69.314803, 48.367995, 48.794464, 37.9362, 24.562735};

    for (int iBin = 0; iBin < NBINS_PUB; iBin++) {
        float pub_data_err_norm = sqrt(pow(pub_data_stat_err[iBin], 2) + pow(pub_data_sys_err[iBin], 2));
        pub_data_hist->SetBinContent(iBin+1, pub_data[iBin]);
        pub_data_hist->SetBinError(iBin+1, pub_data_err_norm);

    };

    // Output file
    TFile* outFile = new TFile("mymain09Hist_HPC_compare_forced_20M_502_tune2.root", "RECREATE");

    // Contribution to HF->mu from each pT-hat bin
    int lineColours[nBins] = {1,2,3,4,6,8,9};
    int markerStyles[nBins] = {20,21,22,23,29,33,34};
    TCanvas *canvasMuonPtContrib = new TCanvas("Muon_sigma_pt_contrib","Muon_sigma_pt_contrib");
    gPad->SetLogy();
    auto legendContrib = new TLegend();
    for(int iBin = 0; iBin < nBins; iBin++) {
        muonPtContribs[iBin]->SetMaximum(10000000);
        muonPtContribs[iBin]->SetMinimum(0.001);
        muonPtContribs[iBin]->SetStats(0);
        muonPtContribs[iBin]->SetLineColor(lineColours[iBin]);
        muonPtContribs[iBin]->SetLineWidth(3);
        muonPtContribs[iBin]->SetMarkerStyle(markerStyles[iBin]);
        muonPtContribs[iBin]->SetMarkerColor(lineColours[iBin]);
        muonPtContribs[iBin]->SetMarkerSize(1.7);
        muonPtContribs[iBin]->Draw("SAME");
        legendContrib->AddEntry(muonPtContribs[iBin],Form("%.1f <= #hat{p}_{T} < %.1f", pT_hat_binEdges[iBin], pT_hat_binEdges[iBin+1]),"p");
    }
    legendContrib->Draw("SAME");

    auto labelContrib = new TLatex();
    labelContrib->DrawLatex(0.0, 0.0, "This Work");
    labelContrib->DrawLatex(0.0, 0.0, "Pythia8 pp #sqrt{s} = 5.02 TeV, Full Monash Tune");
    labelContrib->DrawLatex(0.0, 0.0, "Minimum Bias (QCD), 2.5 < #eta_{#mu} < 4");
    labelContrib->Draw("SAME");

    canvasMuonPtContrib->Write();

    // Pt yield results
    muonPt->SetMinimum(0.01);
    muonPt->Write("HF_muon_pt_contribution");

    // Decay Status Contributions
    TCanvas *canvasMuonPtCompare = new TCanvas("Muon_sigma_pt_compare","Muon_sigma_pt_compare");
    gPad->SetLogy();

    muonPtCompare->SetStats(0);
    muonPtCompare->SetLineColor(1);
    muonPtCompare->SetLineWidth(3);
    muonPtCompare->SetMarkerStyle(20);
    muonPtCompare->SetMarkerColor(1);
    muonPtCompare->SetMarkerSize(1.7);
    muonPtCompare->SetMinimum(100);
    //muonPtCompare->Draw("SAME");
    pub_data_hist->SetStats(0);
    pub_data_hist->SetLineColor(2);
    pub_data_hist->SetLineWidth(3);
    pub_data_hist->SetMarkerStyle(21);
    pub_data_hist->SetMarkerColor(2);
    pub_data_hist->SetMarkerSize(1.7);
    pub_data_hist->SetMinimum(100);
    //pub_data_hist->Draw("SAME");

    auto rp = new TRatioPlot(muonPtCompare, pub_data_hist);
    rp->SetH1DrawOpt("E");
    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("ratio");

    rp->GetUpperPad()->cd();
    // //muonFONLL->Draw("SAME");

    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtCompare,"Pythia8 (Tune 1)","p");
    legendMuon->AddEntry(pub_data_hist,"ALICE Data","p");
    legendMuon->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "This Work");
    labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.02 TeV");
    labelCuts->DrawLatex(0.0, 0.0, "2.5 < #eta < 4, 2 < p_{T} < 20");
    labelCuts->Draw("SAME");

    canvasMuonPtCompare->Write();

    delete outFile;
}