#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain01Macro() {
    // pT-hat bins
    static const int nBins = 6;
    static const double binEdges[nBins+1] = {0.0, 15.0, 30.0, 50.0, 70.0, 100.0, 150.0};

    // W mult analysis bins
    // central mult bins
    const Int_t multBinCount = 7;
    Double_t multBins[multBinCount+1] = {0, 4, 8, 12, 16, 20, 24, 28};
    const float multRawMax = 120.00;

    // // forward mult bins
    // const Int_t multBinCount = 7;
    // Double_t multBins[multBinCount+1] = {0, 3, 6, 9, 12, 15, 18, 21};
    // const float multRawMax = 60.00;

    // either estimate multiplicity in the central or forward region
    bool multCentral = true;
    float multEtaMin, multEtaMax;
    if (multCentral) {
        multEtaMin = -1.0; // range of SPD in ITS
        multEtaMax = 1.0;
    } else {
        multEtaMin = 2.5; // range of V0 (2.8 < eta < 5.1)
        multEtaMax = 4.0;
    }

    // pThat multiplicity
    // histogram for input for multiplicity analysis
    TH1D* multiplicity_events = new TH1D("multiplicity_events", "Integrated Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", multBinCount, multBins);
    TH1D* multiplicity_events_part = new TH1D("multiplicity_events_part", "", multBinCount, multBins);

    // histograms for output/comparison
    TH1D* multiplicity_events_raw = new TH1D("multiplicity_events_raw", "Primary Charged Particle Multiplicity Dependence;N_{ch};#sigma (mb)", 60, 0, multRawMax);
    TH1D* multiplicity_events_raw_part = new TH1D("multiplicity_events_raw_part", "", 60, 0, multRawMax);

    TH1D* multiplicity_events_raw_norm = new TH1D("multiplicity_events_raw_norm", "Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", 30, 0, multRawMax/2);
    TH1D* multiplicity_events_raw_norm_part = new TH1D("multiplicity_events_raw_norm_part", "", 30, 0, multRawMax/2);

    TFile *infile = TFile::Open("mymain01_1M_536_central.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    for(int iBin = 0; iBin < nBins; iBin++) {
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
                if (multiplicity_count > 0) {
                    multiplicity_events_part->Fill(multiplicity_count/2);
                    multiplicity_events_raw_part->Fill(multiplicity_count);
                    multiplicity_events_raw_norm_part->Fill(multiplicity_count/2);
                }
                event_counter++;
                multiplicity_count = 0;
            }
            // in region of interest
            if ((multEtaMin < eta) && (eta < multEtaMax)) {
                multiplicity_count++;
            }
        }
        if (multiplicity_count > 0) {
            multiplicity_events_part->Fill(multiplicity_count/2);
            multiplicity_events_raw_part->Fill(multiplicity_count);
            multiplicity_events_raw_norm_part->Fill(multiplicity_count/2);
        }

        Double_t scaleBin = 1/(*binLuminocity)[iBin];
        //pThat_bins_mult[iBin]->Scale(scaleBin, "width");
        //multiplicity_events->Add(pThat_bins_mult[iBin]);
        multiplicity_events_part->Scale(scaleBin, "width");
        multiplicity_events->Add(multiplicity_events_part);

        multiplicity_events_raw_part->Scale(scaleBin, "width");
        multiplicity_events_raw->Add(multiplicity_events_raw_part);

        multiplicity_events_raw_norm_part->Scale(scaleBin, "width");
        multiplicity_events_raw_norm->Add(multiplicity_events_raw_norm_part);
    }

    std::cout << "<dN_ch/d_eta> for |eta| < 1: " << multiplicity_events_raw_norm->GetMean() << std::endl;
    std::cout << "<dN_ch/d_eta> for |eta| < 1 rebinned for input: " << multiplicity_events->GetMean() << std::endl;

    // Output file
    TFile* outFile = new TFile("mymain01Macro_central.root", "RECREATE");

    multiplicity_events->Draw();
    multiplicity_events->Write();

    multiplicity_events_raw->Write();

    TCanvas *canvasMultRaw = new TCanvas("mult_raw","mult_raw");
    gPad->SetLogy();
    //mb_mult_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_raw->SetMinimum(0.0000001);
    //mb_mult_raw->SetStats(0);
    //multiplicity_events_raw->SetLineColor(2);
    // multiplicity_events_raw->SetMarkerStyle(25);
    // multiplicity_events_raw->SetMarkerColor(1);
    multiplicity_events_raw->Draw("SAME");
    canvasMultRaw->Write();

    TCanvas *canvasMultRawNorm = new TCanvas("mult_raw_norm","mult_raw_norm");
    gPad->SetLogy();
    //mb_mult_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_raw_norm->SetMinimum(0.0000001);
    //mb_mult_raw->SetStats(0);
    //multiplicity_events_raw_norm->SetLineColor(0);
    // multiplicity_events_raw_norm->SetMarkerStyle(25);
    // multiplicity_events_raw_norm->SetMarkerColor(1);
    multiplicity_events_raw_norm->Draw("SAME");
    canvasMultRawNorm->Write();

    delete outFile;
}