#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"

void mymain11Macro_join() {
    // bin constants
    // const Int_t multBinCount = 5;
    // Double_t multBins[multBinCount+1] = {1, 10, 20, 30, 40, 50};

    const Int_t multBinCount = 7;
    Double_t multBins[multBinCount+1] = {1, 5, 10, 15, 20, 25, 30, 35};

    const Int_t pTBinCount = 50;
    const Int_t pTBinMin = 10;
    const Int_t pTBinMax = 80;

    Double_t multBinsAverage[multBinCount+1];

    // Access Minimum Bias data
    TFile* infile_mb = TFile::Open("mymain01Macro.root", "READ");
    TH1D* mb_mult_central = (TH1D*) infile_mb->Get("multiplicity_events_central");
    TH1D* mb_mult_central_raw = (TH1D*) infile_mb->Get("multiplicity_events_central_raw");

    std::cout << "Mean from mult binned: " << mb_mult_central->GetMean() << std::endl;
    std::cout << "Mean from raw binned: " << mb_mult_central_raw->GetMean()/2 << std::endl;

    // normalise charged particle multiplicity by the average
    float mb_mult_average = mb_mult_central_raw->GetMean()/2;
    for (int iMultBin = 0; iMultBin < multBinCount+1; iMultBin++) {
        multBinsAverage[iMultBin] = multBins[iMultBin]/mb_mult_average;
    }

    // W+/- output files from mymain11Macro_multiplicity
    TFile *infile_W_plus = TFile::Open("mymain11Hist_mult_join_plus.root", "READ");
    TFile *infile_W_minus = TFile::Open("mymain11Hist_mult_join_minus.root", "READ");

    // kinematics
    TH1F* W_muon_pt_plus = (TH1F*)infile_W_plus->Get("W_muon_pt");
    W_muon_pt_plus->Scale(pow(10,9));
    TH1F* W_muon_pt_minus = (TH1F*)infile_W_minus->Get("W_muon_pt");
    W_muon_pt_minus->Scale(pow(10,9));
    TH1F* W_muon_pt = new TH1F(*W_muon_pt_plus);
    W_muon_pt->Add(W_muon_pt_minus);

    // multiplicity
    // full cross-section not normalised multiplicity, display bins
    TH1D* W_muon_mult_raw_plus = (TH1D*)infile_W_plus->Get("W_muon_mult_raw");
    TH1D* W_muon_mult_raw_minus = (TH1D*)infile_W_minus->Get("W_muon_mult_raw");
    TH1D* W_muon_mult_raw = new TH1D(*W_muon_mult_raw_plus);
    W_muon_mult_raw->Add(W_muon_mult_raw_minus);

    // full cross-section
    TH1D* W_muon_mult_plus = (TH1D*)infile_W_plus->Get("W_muon_mult");
    TH1D* W_muon_mult_minus = (TH1D*)infile_W_minus->Get("W_muon_mult");
    TH1D* W_muon_mult = new TH1D(*W_muon_mult_plus);
    W_muon_mult->Add(W_muon_mult_minus);

    // cross-section ratio
    TH1D* cs_ratio_plus = (TH1D*)infile_W_plus->Get("cs_ratio");
    TH1D* cs_ratio_minus = (TH1D*)infile_W_minus->Get("cs_ratio");
    TH1D* cs_ratio = new TH1D(*cs_ratio_plus);
    cs_ratio->Add(cs_ratio_minus);

    // forward region yield cross-section
    TH1D* W_muon_yield_mult_plus = (TH1D*)infile_W_plus->Get("W_muon_yield_mult");
    TH1D* W_muon_yield_mult_minus = (TH1D*)infile_W_minus->Get("W_muon_yield_mult");
    TH1D* W_muon_yield_mult = new TH1D(*W_muon_yield_mult_plus);
    W_muon_yield_mult->Add(W_muon_yield_mult_minus);

    // minimum bias yield normalised
    TH1D* W_muon_yield_mult_mb = new TH1D();
    *W_muon_yield_mult_mb = (*W_muon_yield_mult)*(*cs_ratio); // normalise yield by mb ratio
    W_muon_yield_mult_mb->Scale(1/W_muon_yield_mult_mb->Integral());

    // pt yield mult binned normalised
    vector<TH1D*> W_muon_pt_mult_binned(multBinCount);
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        TH1D* W_muon_pt_mult_binned_plus = (TH1D*)infile_W_plus->Get(Form("W_muon_pt_%d", iBin));
        TH1D* W_muon_pt_mult_binned_minus = (TH1D*)infile_W_minus->Get(Form("W_muon_pt_%d", iBin));
        W_muon_pt_mult_binned[iBin] = new TH1D(*W_muon_pt_mult_binned_plus);
        W_muon_pt_mult_binned[iBin]->Add(W_muon_pt_mult_binned_minus);
        W_muon_pt_mult_binned[iBin]->Scale(1/W_muon_yield_mult->Integral());
    }

    // profile for average per mult bin
    TProfile* multProfile = new TProfile("mult_prof","Average of W->#mu Yield per Multiplicity bin;dN_{ch}/d#eta_{|#eta|<1}/<dN_{ch}/d#eta>;<#frac{d^{2}N_{W}}{dp_{T}dy}>", multBinCount, multBinsAverage);

    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        double multBinCenter = multBinsAverage[iMultBin] + ((multBinsAverage[iMultBin+1] - multBinsAverage[iMultBin])/2);
        for (int iPtBin = 1; iPtBin <= pTBinCount; iPtBin++) {
            double yield_pt_bin = (W_muon_pt_mult_binned[iMultBin])->GetBinContent(iPtBin);
            multProfile->Fill(multBinCenter, yield_pt_bin/((4-2.5)*((pTBinMax-pTBinMin)/pTBinCount)), 1);
        }
    }

    TH1D* W_muon_yield_average = multProfile->ProjectionX("yield_average");

    TH1D* W_muon_norm_yield_mult = new TH1D();
    *W_muon_norm_yield_mult = (*W_muon_yield_mult_mb)/(*W_muon_yield_average); // normalise yield by pt average

    // plotting
    TFile* outFile = new TFile("mymain11Hist_joined.root", "RECREATE");

    // kinematics
    TCanvas *canvasMuonPt = new TCanvas("W_muon_pt","W_muon_pt");
    gPad->SetLogy();
    W_muon_pt->SetTitle("Vector Boson Muon Decay p_{T}-differential Cross-Section;p_{T} GeV/c;#frac{d#sigma}{dp_{T}} (pb/GeV/c)");
    W_muon_pt->SetLineColor(1);
    W_muon_pt->SetMarkerStyle(4);
    W_muon_pt->SetMarkerColor(1);
    W_muon_pt->SetStats(0);
    W_muon_pt->Draw("SAME");
    W_muon_pt_plus->SetLineColor(3);
    W_muon_pt_plus->SetMarkerStyle(26);
    W_muon_pt_plus->SetMarkerColor(3);
    W_muon_pt_plus->SetStats(0);
    W_muon_pt_plus->Draw("SAME");
    W_muon_pt_minus->SetLineColor(4);
    W_muon_pt_minus->SetMarkerStyle(27);
    W_muon_pt_minus->SetMarkerColor(4);
    W_muon_pt_minus->SetStats(0);
    W_muon_pt_minus->Draw("SAME");
    auto legendPt = new TLegend();
    legendPt->AddEntry(W_muon_pt,"W #rightarrow #mu","p");
    legendPt->AddEntry(W_muon_pt_plus,"W+ #rightarrow #mu^{+}","p");
    legendPt->AddEntry(W_muon_pt_minus,"W- #rightarrow #mu^{-}","p");
    legendPt->Draw("SAME");
    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "POWHEG pp @ #sqrt{s} = 5.36 TeV");
    labelCuts->DrawLatex(0.0, 0.0, "2.5 < y < 4");
    labelCuts->DrawLatex(0.0, 0.0, "|p| > 4");
    labelCuts->Draw("SAME");
    canvasMuonPt->Write();

    // multiplicity
    // full cross-section not normalised multiplicity, display bins
    TCanvas *canvasMuonMultRaw = new TCanvas("W_muon_mult_raw","W_muon_mult_raw");
    gPad->SetLogy();
    mb_mult_central_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    mb_mult_central_raw->SetMinimum(0.000000000001);
    mb_mult_central_raw->SetStats(0);
    mb_mult_central_raw->SetLineColor(1);
    mb_mult_central_raw->SetMarkerStyle(25);
    mb_mult_central_raw->SetMarkerColor(1);
    mb_mult_central_raw->Draw("SAME");
    W_muon_mult_raw->SetLineColor(3);
    W_muon_mult_raw->SetMarkerStyle(4);
    W_muon_mult_raw->SetMarkerColor(3);
    W_muon_mult_raw->SetStats(0);
    W_muon_mult_raw->Draw("SAME");
    auto legendMultRaw = new TLegend();
    legendMultRaw->AddEntry(mb_mult_central_raw,"Minimum Bias (|#eta| < 1.0)","p");
    legendMultRaw->AddEntry(W_muon_mult_raw,"W #rightarrow #mu","p");
    legendMultRaw->Draw("SAME");
    canvasMuonMultRaw->Write();

    // full cross-section
    TCanvas *canvasMuonMult = new TCanvas("W_muon_mult","W_muon_mult");
    W_muon_mult->Draw("SAME");
    // W_muon_mult_plus->Draw("SAME");
    // W_muon_mult_minus->Draw("SAME");
    canvasMuonMult->Write();

    // cross-section ratio
    TCanvas *canvasMuonRatio = new TCanvas("cs_ratio","cs_ratio");
    cs_ratio->SetLineColor(4);
    cs_ratio->SetMarkerStyle(8);
    cs_ratio->SetMarkerColor(4);
    cs_ratio->Draw("SAME");
    // cs_ratio_plus->Draw("SAME");
    // cs_ratio_minus->Draw("SAME");
    canvasMuonRatio->Write();

    // minimum bias yield normalised
    TCanvas *canvasMuonYield = new TCanvas("W_muon_yield_mult_mb","W_muon_yield_mult_mb");
    W_muon_yield_mult_mb->Draw("SAME");
    canvasMuonYield->Write();

    // pt yield mult binned normalised
    TCanvas *canvasMuonPtBins = new TCanvas("W_muon_pt_mult_binned","W_muon_pt_mult_binned");
    gPad->SetLogy();
    auto legendPtBins = new TLegend();
    int lineColors[multBinCount] = {1,2,3,4,6,28,7};
    int lineMarkers[multBinCount] = {24,25,26,27,28,30,32};
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        W_muon_pt_mult_binned[iBin]->SetStats(0);
        W_muon_pt_mult_binned[iBin]->SetMinimum(0.00001);
        W_muon_pt_mult_binned[iBin]->SetLineColor(lineColors[iBin]);
        W_muon_pt_mult_binned[iBin]->SetMarkerStyle(lineMarkers[iBin]);
        W_muon_pt_mult_binned[iBin]->SetMarkerColor(lineColors[iBin]);
        W_muon_pt_mult_binned[iBin]->Draw("SAME");
        legendPtBins->AddEntry(W_muon_pt_mult_binned[iBin], Form("%d <= dN_{ch}/d#eta < %d", (int)multBins[iBin], (int)multBins[iBin+1]),"p");
    }
    legendPtBins->Draw("SAME");
    canvasMuonPtBins->Write();

    // profile average per mult bin
    TCanvas *canvasMuonAverage = new TCanvas("W_muon_average_mult","W_muon_average_mult");
    W_muon_yield_average->Draw();
    canvasMuonAverage->Write();

    // self-normalised mult dependence
    TCanvas *canvasYieldMultDep = new TCanvas("yield_mult_dep","yield_mult_dep");
    W_muon_norm_yield_mult->SetTitle("Self Normalised Charged Particle Dependence of W_{#mu}");
    W_muon_norm_yield_mult->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}dy}/<#frac{d^{2}#sigma}{dp_{T}dy}>");
    W_muon_norm_yield_mult->SetStats(0);
    W_muon_norm_yield_mult->SetLineColor(4);
    W_muon_norm_yield_mult->SetMarkerStyle(8);
    W_muon_norm_yield_mult->SetMarkerColor(4);
    // fit straight line
    TF1* fitLinear = new TF1("fitLinear", "pol1");
    W_muon_norm_yield_mult->Fit("fitLinear");
    W_muon_norm_yield_mult->Draw();
    // add plot of y=x
    double y[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    double x[8] = {0, 1, 2, 3, 4, 5, 6, 7};
    auto gLine = new TGraph(8, x, y);
    gLine->Draw("SAME");
    auto legendLinear = new TLegend();
    legendLinear->AddEntry(fitLinear, "Linear Best Fit");
    legendLinear->AddEntry(gLine, "y = x", "l");
    legendLinear->Draw("SAME");
    canvasYieldMultDep->Write();

    // write muon pt for semileptonic muon cross-section
    //W_muon_pt->Scale(pow(10,9));
    W_muon_pt->Write("W_muon_pt_forward");

    delete outFile;
}