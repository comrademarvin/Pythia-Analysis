#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain01Macro() {
    // pT-hat bins
    static const int nBins = 6;
    static const double binEdges[nBins+1] = {0.0, 15.0, 30.0, 50.0, 70.0, 100.0, 150.0};

    // estimate multiplicity in different regions
    const int nRegions = 5;
    const int selectedRegion = 4; // select the desired multiplicity estimation region here
    const string region_label[nRegions] = {"central", "forward", "V0C", "central_CR_off", "central_MPI_off"};
    const float region_eta_min[nRegions] = {-1.0, 2.5, 1.7, -1.0, -1.0};
    const float region_eta_max[nRegions] = {1.0, 4.0, 3.7, 1.0, 1.0};
    const float region_plot_max[nRegions] = {100.0, 60.0, 100.00, 120.0, 60.0};
    const int region_plot_bins[nRegions] = {50, 30, 50, 60, 30};
    const float region_eta_width = region_eta_max[selectedRegion] - region_eta_min[selectedRegion];

    // multiplicity analysis bins
    const Int_t multBinCount = 7;
    Double_t multBins[nRegions][multBinCount+1] = {
        {0, 4, 8, 12, 16, 20, 24, 28},
        {0, 4, 8, 12, 16, 20, 24, 32},
        {0, 4, 8, 12, 16, 20, 24, 32},
        {0, 5, 10, 15, 20, 25, 30, 35},
        {0, 4, 8, 12, 16, 20, 24, 28}
    };

    // pThat multiplicity
    // histogram for input for multiplicity analysis
    TH1D* multiplicity_events = new TH1D("multiplicity_events", "Integrated Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", 
                                            multBinCount, multBins[selectedRegion]);
    TH1D* multiplicity_events_part = new TH1D("multiplicity_events_part", "", multBinCount, multBins[selectedRegion]);

    // histograms for output/comparison
    TH1D* multiplicity_events_raw = new TH1D("multiplicity_events_raw", "Primary Charged Particle Multiplicity Dependence;N_{ch};#sigma (mb)", 
                                                region_plot_bins[selectedRegion], 0, region_plot_max[selectedRegion]);
    TH1D* multiplicity_events_raw_part = new TH1D("multiplicity_events_raw_part", "", region_plot_bins[selectedRegion], 0, region_plot_max[selectedRegion]);

    TH1D* multiplicity_events_raw_norm = new TH1D("multiplicity_events_raw_norm", "Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", 
                                                    region_plot_bins[selectedRegion]/2, 0, region_plot_max[selectedRegion]/region_eta_width);
    TH1D* multiplicity_events_raw_norm_part = new TH1D("multiplicity_events_raw_norm_part", "", 
                                                        region_plot_bins[selectedRegion]/2, 0, region_plot_max[selectedRegion]/region_eta_width);

    TFile *infile = TFile::Open("mymain01_10k_536_4pi_MPI_off.root", "READ");

    TNtuple *genInfo = (TNtuple*)infile->Get("genInfo");

    // iterate over bins and access generated cross-section info
    float weightSum, sigmaGen, sigmaErr;
    genInfo->SetBranchAddress("weightSum",&weightSum);
    genInfo->SetBranchAddress("sigmaGen",&sigmaGen);
    genInfo->SetBranchAddress("sigmaErr",&sigmaErr);

    for(int iBin = 0; iBin < nBins; iBin++) {
        genInfo->GetEntry(iBin);

        TNtuple *chargedTuple = (TNtuple*)infile->Get(Form("charged%d", iBin));

        multiplicity_events_part->Reset();
        multiplicity_events_raw_part->Reset();
        multiplicity_events_raw_norm_part->Reset();

        int particle_count = chargedTuple->GetEntries();
        Float_t eventIndex, eta, pAbs, id, pt;
        chargedTuple->SetBranchAddress("eventTag", &eventIndex);
        chargedTuple->SetBranchAddress("eta", &eta);
        chargedTuple->SetBranchAddress("pAbs", &pAbs);
        chargedTuple->SetBranchAddress("id", &id);
        chargedTuple->SetBranchAddress("pt", &pt);

        // iterate through the events
        int event_counter = 0;
        int multiplicity_count = 0;
        for (int iParticle = 0; iParticle < particle_count; iParticle++) {
            chargedTuple->GetEntry(iParticle);
            if (eventIndex != event_counter) {
                if (multiplicity_count > 0) { // exclude events where N_ch = 0
                    multiplicity_events_part->Fill(multiplicity_count/region_eta_width);
                    multiplicity_events_raw_part->Fill(multiplicity_count);
                    multiplicity_events_raw_norm_part->Fill(multiplicity_count/region_eta_width);
                }
                event_counter++;
                multiplicity_count = 0;
            }
            // in region of interest
            if ((region_eta_min[selectedRegion] < eta) && (eta < region_eta_max[selectedRegion])) {
                multiplicity_count++;
            }
        }
        if (multiplicity_count > 0) {
            multiplicity_events_part->Fill(multiplicity_count/region_eta_width);
            multiplicity_events_raw_part->Fill(multiplicity_count);
            multiplicity_events_raw_norm_part->Fill(multiplicity_count/region_eta_width);
        }

        Double_t scaleBin = sigmaGen/weightSum;

        multiplicity_events_part->Scale(scaleBin, "width");
        multiplicity_events->Add(multiplicity_events_part);

        multiplicity_events_raw_part->Scale(scaleBin, "width");
        multiplicity_events_raw->Add(multiplicity_events_raw_part);

        multiplicity_events_raw_norm_part->Scale(scaleBin, "width");
        multiplicity_events_raw_norm->Add(multiplicity_events_raw_norm_part);
    }

    std::cout << "<dN_ch/d_eta> for selected region: " << multiplicity_events_raw_norm->GetMean() << std::endl;
    std::cout << "<dN_ch/d_eta> rebinned for input: " << multiplicity_events->GetMean() << std::endl;

    // Output file
    TFile* outFile = new TFile("mymain01Macro_central_MPI_off.root", "RECREATE");

    multiplicity_events->Draw();
    multiplicity_events->Write();

    multiplicity_events_raw->Write();

    TCanvas *canvasMultRaw = new TCanvas("mult_raw","mult_raw");
    gPad->SetLogy();
    //mb_mult_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_raw->SetMinimum(0.000001);
    //mb_mult_raw->SetStats(0);
    //multiplicity_events_raw->SetLineColor(2);
    // multiplicity_events_raw->SetMarkerStyle(25);
    // multiplicity_events_raw->SetMarkerColor(1);
    multiplicity_events_raw->Draw("SAME");
    canvasMultRaw->Write();

    TCanvas *canvasMultRawNorm = new TCanvas("mult_raw_norm","mult_raw_norm");
    gPad->SetLogy();
    //mb_mult_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_raw_norm->SetMinimum(0.000001);
    //mb_mult_raw->SetStats(0);
    //multiplicity_events_raw_norm->SetLineColor(0);
    // multiplicity_events_raw_norm->SetMarkerStyle(25);
    // multiplicity_events_raw_norm->SetMarkerColor(1);
    multiplicity_events_raw_norm->Draw("SAME");
    canvasMultRawNorm->Write();

    delete outFile;
}