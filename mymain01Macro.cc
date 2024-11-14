#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain01Macro() {
    // pT-hat bins
    static const int nBins = 6;
    static const double binEdges[nBins+1] = {0.0, 15.0, 30.0, 50.0, 70.0, 100.0, 150.0};
    // const Int_t multBinCount = 5;
    // Double_t multBins[multBinCount+1] = {1, 10, 20, 30, 40, 50};
    // const Int_t multBinCount = 7;
    // Double_t multBins[multBinCount+1] = {1, 5, 10, 15, 20, 25, 30, 35};
    const Int_t multBinCount = 7;
    Double_t multBins[multBinCount+1] = {0, 4, 8, 12, 16, 20, 24, 28};

    // pThat multiplicity
    // histogram for input for multiplicity analysis
    TH1D* multiplicity_events_central = new TH1D("multiplicity_events_central", "Integrated Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", multBinCount, multBins);
    TH1D* multiplicity_events_central_part = new TH1D("multiplicity_events_central_part", "", multBinCount, multBins);

    // histograms for output/comparison
    TH1D* multiplicity_events_central_raw = new TH1D("multiplicity_events_central_raw", "Primary Charged Particle Multiplicity Dependence;N_{ch};#sigma (mb)", 60, 0, 120);
    TH1D* multiplicity_events_central_raw_part = new TH1D("multiplicity_events_central_raw_part", "", 60, 0, 120);

    TH1D* multiplicity_events_central_raw_norm = new TH1D("multiplicity_events_central_raw_norm", "Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", 30, 0, 60);
    TH1D* multiplicity_events_central_raw_norm_part = new TH1D("multiplicity_events_central_raw_norm_part", "", 30, 0, 60);

    TFile *infile = TFile::Open("mymain01_1M_536.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    for(int iBin = 0; iBin < nBins; iBin++) {
        TNtuple *chargedTuple = (TNtuple*)infile->Get(Form("charged%d", iBin));

        multiplicity_events_central_part->Reset();
        multiplicity_events_central_raw_part->Reset();
        multiplicity_events_central_raw_norm_part->Reset();

        int particle_count = chargedTuple->GetEntries();
        Float_t eventIndex, eta, pAbs, id, pt;
        chargedTuple->SetBranchAddress("eventTag", &eventIndex);
        chargedTuple->SetBranchAddress("eta", &eta);
        chargedTuple->SetBranchAddress("pAbs", &pAbs);
        chargedTuple->SetBranchAddress("id", &id);
        chargedTuple->SetBranchAddress("pt", &pt);

        // iterate through the events
        int event_counter = 0;
        int multiplicity_count_central = 0;
        for (int iParticle = 0; iParticle < particle_count; iParticle++) {
            chargedTuple->GetEntry(iParticle);
            if (eventIndex != event_counter) {
                if (multiplicity_count_central > 0) {
                    multiplicity_events_central_part->Fill(multiplicity_count_central/2);
                    multiplicity_events_central_raw_part->Fill(multiplicity_count_central);
                    multiplicity_events_central_raw_norm_part->Fill(multiplicity_count_central/2);
                }
                event_counter++;
                multiplicity_count_central = 0;
            }
            // Central region
            if (abs(eta) < 1) {
                multiplicity_count_central++;
            }
            // // MFT region
            // if (eta > -3.6 && eta < -2.5) {
            //     multiplicity_count_MFT++;
            //     pt_test->Fill(pt);
            // }
            // // Forward region (muon)
            // if (abs(id) == 13 && eta > -4.0 && eta < -2.5) {
            //     multiplicity_count_Muon++;
            // }
        }
        if (multiplicity_count_central > 0) {
            multiplicity_events_central_part->Fill(multiplicity_count_central/2);
            multiplicity_events_central_raw_part->Fill(multiplicity_count_central);
            multiplicity_events_central_raw_norm_part->Fill(multiplicity_count_central/2);
        }

        Double_t scaleBin = 1/(*binLuminocity)[iBin];
        //pThat_bins_mult[iBin]->Scale(scaleBin, "width");
        //multiplicity_events->Add(pThat_bins_mult[iBin]);
        multiplicity_events_central_part->Scale(scaleBin, "width");
        multiplicity_events_central->Add(multiplicity_events_central_part);

        multiplicity_events_central_raw_part->Scale(scaleBin, "width");
        multiplicity_events_central_raw->Add(multiplicity_events_central_raw_part);

        multiplicity_events_central_raw_norm_part->Scale(scaleBin, "width");
        multiplicity_events_central_raw_norm->Add(multiplicity_events_central_raw_norm_part);
    }

    std::cout << "<dN_ch/d_eta> for |eta| < 1: " << multiplicity_events_central_raw_norm->GetMean() << std::endl;
    std::cout << "<dN_ch/d_eta> for |eta| < 1 rebinned for input: " << multiplicity_events_central->GetMean() << std::endl;

    // Output file
    TFile* outFile = new TFile("mymain01Macro.root", "RECREATE");

    multiplicity_events_central->Draw();
    multiplicity_events_central->Write();

    multiplicity_events_central_raw->Write();

    TCanvas *canvasMultRaw = new TCanvas("mult_raw","mult_raw");
    gPad->SetLogy();
    //mb_mult_central_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_central_raw->SetMinimum(0.0000001);
    //mb_mult_central_raw->SetStats(0);
    //multiplicity_events_central_raw->SetLineColor(2);
    // multiplicity_events_central_raw->SetMarkerStyle(25);
    // multiplicity_events_central_raw->SetMarkerColor(1);
    multiplicity_events_central_raw->Draw("SAME");
    canvasMultRaw->Write();

    TCanvas *canvasMultRawNorm = new TCanvas("mult_raw_norm","mult_raw_norm");
    gPad->SetLogy();
    //mb_mult_central_raw->GetYaxis()->SetTitle("#frac{d#sigma}{dN_{ch}} (mb)");
    multiplicity_events_central_raw_norm->SetMinimum(0.0000001);
    //mb_mult_central_raw->SetStats(0);
    //multiplicity_events_central_raw_norm->SetLineColor(0);
    // multiplicity_events_central_raw_norm->SetMarkerStyle(25);
    // multiplicity_events_central_raw_norm->SetMarkerColor(1);
    multiplicity_events_central_raw_norm->Draw("SAME");
    canvasMultRawNorm->Write();

    delete outFile;
}