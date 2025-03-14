#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"
#include <TNtuple.h>

void mymain11Macro_multiplicity() {
    // estimate multiplicity in different regions
    const int nRegions = 7;
    const int selectedRegion = 0; // select the desired multiplicity estimation region here
    const string region_label[nRegions] = {"central", "forward", "V0C", "central_CR_off", "central_MPI_off", "central_CR_mode_2", "central_MPI_level_0"};
    const float region_eta_min[nRegions] = {-1.0, 2.5, 1.7, -1.0, -1.0, -1.0, -1.0};
    const float region_eta_max[nRegions] = {1.0, 4.0, 3.7, 1.0, 1.0, 1.0, 1.0};
    const float region_plot_max[nRegions] = {100.0, 70.0, 90.00, 120.0, 60.0, 100.0, 100.0};
    const int region_plot_bins[nRegions] = {50, 35, 45, 60, 30, 50, 50};
    const float region_eta_width = region_eta_max[selectedRegion] - region_eta_min[selectedRegion];

    // multiplicity analysis bins
    const Int_t multBinCount = 7;
    Double_t multBins[nRegions][multBinCount+1] = {
        {0, 4, 8, 12, 16, 20, 24, 28},
        {0, 4, 8, 12, 16, 20, 24, 28},
        {0, 4, 8, 12, 16, 20, 24, 28},
        {0, 5, 10, 15, 20, 25, 30, 35},
        {0, 2, 4, 6, 8, 10, 12, 14},
        {0, 4, 8, 12, 16, 20, 24, 28},
        {0, 4, 8, 12, 16, 20, 24, 28}
    };

    // pT-bins for mult dependence
    const Int_t pTBinCount = 60;
    const float pTBinMin = 30.0;
    const float pTBinMax = 60.0;

    // Access Minimum Bias data
    TFile* infile_mb = TFile::Open("mymain01Macro_central.root", "READ");
    TH1D* mb_mult = (TH1D*) infile_mb->Get("multiplicity_events");
    TH1D* mb_mult_raw = (TH1D*) infile_mb->Get("multiplicity_events_raw");
    //TH1D* mb_mult_raw_norm = (TH1D*) infile_mb->Get("mult_raw_norm");

    // normalise charged particle multiplicity by the average
    float mb_mult_average = mb_mult_raw->GetMean()/region_eta_width; // use mult average before re-binning
    Double_t multBinsAverage[multBinCount+1];
    for (int iMultBin = 0; iMultBin < multBinCount+1; iMultBin++) {
        multBinsAverage[iMultBin] = multBins[selectedRegion][iMultBin]/mb_mult_average;
    } 

    // Access W+/- showered data
    TFile* infile = TFile::Open("mymain11_W+_10M_central.root", "READ");

    // W->muon distributions (pt+y, event multiplicities, multiplicity bins, yield)
    TH1F* W_muon_pt = new TH1F("W_muon_pt", "", 40, 0, 80);
    TH1F* W_muon_y = new TH1F("W_muon_y", "W->Muon parameter space for rapidity;y;#sigma (mb)", 120, -6, 6);
    TH1F* W_muon_eta = new TH1F("W_muon_eta", "W->Muon parameter space for pseudorapidity;#eta;#sigma (mb)", 120, -6, 6);
    TH1D* W_muon_mult_raw = new TH1D("W_muon_mult_raw", "Primary charged particle dependence;N_{ch};#sigma_{W->#mu} (mb)", 
                                        region_plot_bins[selectedRegion], 0, region_plot_max[selectedRegion]);
    TH1D* W_muon_multiplicity = new TH1D("W_muon_multiplicity", "Primary charged particle dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma_{W->#mu} (mb)", 
                                        multBinCount, multBins[selectedRegion]);
    TH1D* W_muon_yield_mult = new TH1D("W_muon_yield_mult", "Primary charged particle dependence;dN_{ch}/d#eta_{|#eta|<1}/<dN_{ch}/d#eta>;#frac{d^{2}N}{dp_{T}dy}", 
                                        multBinCount, multBinsAverage);

    vector<TH1D*> W_muon_pt_mult_binned(multBinCount);
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        W_muon_pt_mult_binned[iBin] = new TH1D(Form("W_muon_pt_%d", iBin), "Muon Yield per Multiplicity bin; p_{T}; #frac{1}{#sigma_{W}}#frac{d^{2}#sigma}{dp_{T}dy}", 
                                                pTBinCount, pTBinMin, pTBinMax);
    }

    // Read in generated cross-section as integrated luminocity
    std::vector<double> *genInfo;
    infile->GetObject("genInfo", genInfo);

    // Read in multiplicities of events
    std::vector<int> *event_mult;
    infile->GetObject("event_multiplicity", event_mult);
    int event_count = event_mult->size();

    // iterate through muons and connect to event multiplicities
    TNtuple *muonTuple = (TNtuple*)infile->Get("W_muon");
    int muons_count = muonTuple->GetEntries();
    Float_t eventIndex, eta, pt, rapidity, pAbs;
    muonTuple->SetBranchAddress("eventTag", &eventIndex);
    muonTuple->SetBranchAddress("y", &rapidity);
    muonTuple->SetBranchAddress("eta", &eta);
    muonTuple->SetBranchAddress("pt", &pt);
    muonTuple->SetBranchAddress("pAbs", &pAbs);

    for (int iMuon = 0; iMuon < muons_count; iMuon++) {
        muonTuple->GetEntry(iMuon);
        W_muon_y->Fill(rapidity);
        W_muon_eta->Fill(eta);
        int muon_event_mult = (*event_mult)[static_cast<int>(eventIndex)];
        double muon_event_mult_norm = muon_event_mult/region_eta_width;
        W_muon_mult_raw->Fill(muon_event_mult);
        W_muon_multiplicity->Fill(muon_event_mult_norm);

        // kinematics
        if ((rapidity >= 2.5) && (rapidity <= 4) && (abs(pAbs) > 4)) {
            W_muon_pt->Fill(pt);
        }

        // multiplicity
        if ((rapidity >= 2.5) && (rapidity <= 4) && (pt >= pTBinMin) && (pt <= pTBinMax) && (abs(pAbs) > 4)) {
            double event_mult_average_norm = muon_event_mult_norm/mb_mult_average;
            W_muon_yield_mult->Fill(event_mult_average_norm);
            for (int iBin = 0; iBin < multBinCount; iBin++) {
                if ((event_mult_average_norm >= multBinsAverage[iBin]) && (event_mult_average_norm < multBinsAverage[iBin+1])) {
                    W_muon_pt_mult_binned[iBin]->Fill(pt);
                    break;
                }
            }
        }
    };

    // Cross Section Scaling + cuts normalisation
    W_muon_mult_raw->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_multiplicity->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_pt->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_y->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_eta->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_yield_mult->Scale((*genInfo)[1]/(*genInfo)[0], "width");
    W_muon_yield_mult->Scale(1/(pTBinMax-pTBinMin)); // normalise by integration width

    vector<TH1D*> W_muon_pt_mult_binned_cs(multBinCount);
    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        W_muon_pt_mult_binned[iMultBin]->Scale((*genInfo)[1]/(*genInfo)[0], "width");
        W_muon_pt_mult_binned[iMultBin]->Scale(1/(4-2.5)); // normalise by integration width
        W_muon_pt_mult_binned_cs[iMultBin] = new TH1D(*(W_muon_pt_mult_binned[iMultBin]));
    }

    // Compute the multiplicity dependant yield

    // Cross section ratio multiplicity binned
    TH1D* cs_ratio = new TH1D();

    *cs_ratio = (*W_muon_multiplicity)/(*mb_mult);

    TH1D* W_muon_yield_mult_mb = new TH1D();
    *W_muon_yield_mult_mb = (*W_muon_yield_mult)*(*cs_ratio); // normalise yield by mb ratio
    W_muon_yield_mult_mb->Scale(1/W_muon_yield_mult_mb->Integral());

    // Determine yield of W->muon as function of multiplicity
    TProfile* multProfile = new TProfile("mult_prof","Average of W->#mu Yield per Multiplicity bin;dN_{ch}/d#eta_{|#eta|<1}/<dN_{ch}/d#eta>;<#frac{d^{2}N_{W}}{dp_{T}dy}>", multBinCount, multBinsAverage);

    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        double multBinCenter = multBinsAverage[iMultBin] + ((multBinsAverage[iMultBin+1] - multBinsAverage[iMultBin])/2);
        W_muon_pt_mult_binned[iMultBin]->Scale(1/(W_muon_yield_mult->Integral()));

        for (int iPtBin = 1; iPtBin <= pTBinCount; iPtBin++) {
            double yield_pt_bin = (W_muon_pt_mult_binned[iMultBin])->GetBinContent(iPtBin);
            multProfile->Fill(multBinCenter, yield_pt_bin, 1);
        }
    }

    TH1D* W_muon_yield_average = multProfile->ProjectionX("yield_average");
    W_muon_yield_average->Scale(1, "width"); // normalise for bin-wi

    TH1D* W_muon_norm_yield_mult = new TH1D();
    *W_muon_norm_yield_mult = (*W_muon_yield_mult_mb)/(*W_muon_yield_average); // normalise yield by pt average

    // Output file
    TFile* outFile = new TFile("mymain11Hist_multiplicity.root", "RECREATE");

    TCanvas *canvasMuonPtY = new TCanvas("Muon_y","Muon_y");
    W_muon_y->Draw();
    canvasMuonPtY->Write();

    TCanvas *canvasMuonMult = new TCanvas("Muon_mult","Muon_mult");
    gPad->SetLogy();
    mb_mult->SetMinimum(0.00000000000001);
    mb_mult->SetLineColor(1);
    mb_mult->Draw("SAME");
    W_muon_multiplicity->SetStats(0);
    W_muon_multiplicity->SetMinimum(0.00000000000001);
    W_muon_multiplicity->Draw("SAME");
    auto legendMult = new TLegend();
    legendMult->AddEntry(mb_mult,"Minimum Bias","l");
    legendMult->AddEntry(W_muon_multiplicity,"W->#mu","l");
    legendMult->Draw("SAME");
    canvasMuonMult->Write();

    TCanvas *canvasMultBinned = new TCanvas("mult_binned_ratio","mult_binned_ratio");
    cs_ratio->SetStats(0);
    //cs_ratio->SetMinimum(0.00000001);
    cs_ratio->SetTitle("Ratio of W muon Cross Section wrt Minimum Bias");
    cs_ratio->GetYaxis()->SetTitle("#sigma_{W#rightarrow#mu}/#sigma_{MB}");
    cs_ratio->Draw();
    canvasMultBinned->Write();

    TCanvas *canvasYield = new TCanvas("yield_mult_binned","yield_mult_binned");
    gPad->SetLogy();
    W_muon_yield_mult_mb->SetTitle("Minimum Bias Scaled W#rightarrow#mu Yield");
    W_muon_yield_mult_mb->GetYaxis()->SetTitle("#frac{1}{sigma_{W_{#mu},MB}}#frac{d^{2}#sigma_{W_{#mu},MB}}{dp_{T}d#eta}");
    //W_muon_yield_mult_mb->SetMinimum(0.00001);
    W_muon_yield_mult_mb->Draw();
    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "2.5 < #eta_{#mu} < 4, 30 < p_{T} < 60");
    labelCuts->Draw("SAME");
    canvasYield->Write();

    TCanvas *canvasYieldBinned = new TCanvas("pt_yield_mult_binned","pt_yield_mult_binned");
    gPad->SetLogy();
    auto legendPtBins = new TLegend();
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        W_muon_pt_mult_binned[iBin]->SetStats(0);
        W_muon_pt_mult_binned[iBin]->SetMinimum(0.00001);
        W_muon_pt_mult_binned[iBin]->SetLineColor(iBin+1);
        W_muon_pt_mult_binned[iBin]->Draw("SAME");
        legendPtBins->AddEntry(W_muon_pt_mult_binned[iBin], Form("%.1f <= N_{ch} < %.1f", multBins[selectedRegion][iBin], multBins[selectedRegion][iBin+1]),"l");
    }
    // auto labelCuts = new TLatex();
    // labelCuts->DrawLatex(0.0, 0.0, "2.5 < y < 4");
    // labelCuts->Draw("SAME");
    legendPtBins->Draw("SAME");
    canvasYieldBinned->Write();

    TCanvas *canvasYieldAverage = new TCanvas("yield_mult_average","yield_mult_average");
    gPad->SetLogy();
    multProfile->SetStats(0);
    multProfile->Draw();
    canvasYieldAverage->Write();

    TCanvas *canvasYieldMultDep = new TCanvas("yield_mult_dep","yield_mult_dep");
    W_muon_norm_yield_mult->SetTitle("Self Normalised Charged Particle Dependence of W->muon");
    W_muon_norm_yield_mult->GetYaxis()->SetTitle("#frac{d^{2}N}{dp_{T}dy}/<#frac{d^{2}N}{dp_{T}dy}>");
    W_muon_norm_yield_mult->SetStats(0);
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
    legendLinear->AddEntry(W_muon_norm_yield_mult, "W #rightarrow #mu", "l");
    legendLinear->AddEntry(fitLinear, "Linear Best Fit");
    legendLinear->AddEntry(gLine, "y=x", "l");
    legendLinear->Draw("SAME");
    canvasYieldMultDep->Write();

    delete outFile;

    // save histograms to join W+/- analysis
    TFile* outFileJoin = new TFile("mymain11Hist_mult_join_plus_central.root", "RECREATE");

    // kinematics
    W_muon_pt->Write();
    W_muon_y->Write();
    W_muon_eta->Write();

    // multiplicity dependence
    W_muon_mult_raw->Write("W_muon_mult_raw");
    W_muon_multiplicity->Write("W_muon_mult");
    cs_ratio->Write("cs_ratio");
    W_muon_yield_mult->Write("W_muon_yield_mult");

    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        W_muon_pt_mult_binned_cs[iMultBin]->Write();
    }

    delete outFileJoin;
};