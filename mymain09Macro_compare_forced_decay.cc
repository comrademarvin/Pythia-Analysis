#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_compare_forced_decay() {
    TFile *infile = TFile::Open("mymain09_10k.root", "READ");

    TNtuple *genInfo = (TNtuple*)infile->Get("genInfo");

    // branching ratios for forced decays
    const double D_BR = 0.176;
    const double D0_BR = 0.067;

    TFile* outFile = new TFile("mymain09Hist_compare_forced.root", "RECREATE");

    // published data bins
    const Int_t NBINS_PUB = 26;
    Double_t edges_pub[NBINS_PUB + 1] = {2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13, 14, 15, 16.5, 18, 20};

    // HF contribution histograms
    TH1F *muonPtTotal = new TH1F("muon_full","HF Muon Decay Cross-Section;p_{T} (GeV/c);#frac{d#sigma_{c,b->#mu}}{dp_{T}} (pb/GeV/c)", NBINS_PUB, edges_pub);
    TH1F *muonPtD = new TH1F("muon_pt_D","", NBINS_PUB, edges_pub);
    TH1F *muonPtD0 = new TH1F("muon_pt_D0","", NBINS_PUB, edges_pub);
    TH1F *muonPtRest = new TH1F("muon_pt_rest","", NBINS_PUB, edges_pub);

    // iterate over bins and access generated cross-section info
    float weightSum, sigmaGen, sigmaErr;
    genInfo->SetBranchAddress("weightSum",&weightSum);
    genInfo->SetBranchAddress("sigmaGen",&sigmaGen);
    genInfo->SetBranchAddress("sigmaErr",&sigmaErr);

    for(int iBin = 0; iBin < genInfo->GetEntries(); iBin++) {
        genInfo->GetEntry(iBin);

        // access the associated bin forward muons
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonPtD->Reset();
        muonPtD0->Reset();
        muonPtRest->Reset();

        if (iBin == 0) {
            muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (pt < 10.0) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (pt < 10.0) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (pt < 10.0) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        } else {
            muonTuple->Draw("pt>>muon_pt_D", "(firstMother == 411) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_D0", "(firstMother == 421) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
            muonTuple->Draw("pt>>muon_pt_rest", "(firstMother != 411) && (firstMother != 421) && (pt > 2.0) && (decayStatus == 0 || decayStatus == 1 || decayStatus == 2)");
        }

        // Scale by Branching Ratio
        muonPtD->Scale(D_BR);
        muonPtD0->Scale(D0_BR);

        muonPtRest->Add(muonPtD);
        muonPtRest->Add(muonPtD0);

        // normalise to cross-section
        TH1F *sigmaGenHistPt = new TH1F("sigma_gen_pt","", NBINS_PUB, edges_pub);
        for (int iBin = 0; iBin < NBINS_PUB; iBin++) {
            sigmaGenHistPt->SetBinContent(iBin, sigmaGen);
            sigmaGenHistPt->SetBinError(iBin, sigmaErr);
        }
        sigmaGenHistPt->Scale(pow(10,9)); // mb to pb conversion

        muonPtRest->Scale(1/weightSum, "width");
        muonPtRest->Multiply(sigmaGenHistPt);
        muonPtRest->Write(Form("hf_muon_pt_bin_%d", iBin));
        muonPtTotal->Add(muonPtRest);
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

    // compare simulated HF to experimental data
    TCanvas *canvasMuonPtCompare = new TCanvas("Muon_sigma_pt_compare","Muon_sigma_pt_compare");
    gPad->SetLogy();

    muonPtTotal->SetStats(0);
    muonPtTotal->SetLineColor(1);
    muonPtTotal->SetMarkerStyle(33);
    muonPtTotal->SetMarkerColor(1);
    muonPtTotal->Draw("SAME");
    pub_data_hist->SetStats(0);
    pub_data_hist->SetLineColor(2);
    pub_data_hist->SetMarkerStyle(21);
    pub_data_hist->SetMarkerColor(2);
    pub_data_hist->Draw("SAME");

    // auto rp = new TRatioPlot(muonPtTotal, pub_data_hist);
    // rp->SetH1DrawOpt("E");
    // rp->Draw();

    // rp->GetUpperPad()->cd();
    // // //muonFONLL->Draw("SAME");

    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtTotal,"Pythia (Monash Tune)","p");
    legendMuon->AddEntry(pub_data_hist,"Published","p");
    //legendMuon->AddEntry(muonFONLL,"FONLL","l");
    legendMuon->Draw("SAME");

    // auto labelCuts = new TLatex();
    // labelCuts->DrawLatex(0.0, 0.0, "pp #sqrt{s} = 5.02 TeV, 2.5 < y < 4");
    // labelCuts->DrawLatex(0.0, 0.0, "2 < p_{T} < 20");
    // labelCuts->DrawLatex(0.0, 0.0, "|p| > 4");
    // labelCuts->Draw("SAME");

    canvasMuonPtCompare->Write();

    delete outFile;
}