#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain01Macro() {
    // pT-hat bins
    const int pThat_bins_count = 1;
    double pThat_bins[pThat_bins_count + 1] = {0.0, 100.0};

    // pThat multiplicity
    TH1F* multiplicity_events_central = new TH1F("multiplicity_events_central", "Multiplicity of all events;N_{ch};#sigma (mb)", 100, 0, 200);
    TH1F* multiplicity_events_MFT = new TH1F("multiplicity_events_MFT", "", 100, 0, 200);

    // vector<TH1F*> pThat_bins_mult(pThat_bins_count);

    // for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
    //     pThat_bins_mult[iBin] = new TH1F(Form("pThat_bin_%d", iBin), "", 150, 0, 300);
    // }

    TFile *infile = TFile::Open("mymain01.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    for(int iBin = 0; iBin < pThat_bins_count; iBin++) {
        TNtuple *chargedTuple = (TNtuple*)infile->Get(Form("charged%d", iBin));

        int particle_count = chargedTuple->GetEntries();
        Float_t eventIndex, eta, pAbs, id;
        chargedTuple->SetBranchAddress("eventTag", &eventIndex);
        chargedTuple->SetBranchAddress("eta", &eta);
        chargedTuple->SetBranchAddress("pAbs", &pAbs);
        chargedTuple->SetBranchAddress("id", &id);

        // iterate through the events
        int event_counter = 0;
        int multiplicity_count_central = 0;
        int multiplicity_count_MFT = 0;
        for (int iParticle = 0; iParticle < particle_count; iParticle++) {
            chargedTuple->GetEntry(iParticle);
            if (eventIndex != event_counter) {
                if (multiplicity_count_central > 0) {
                    multiplicity_events_central->Fill(multiplicity_count_central); 
                }
                if (multiplicity_count_MFT > 0) {
                    multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
                }
                event_counter++;
                multiplicity_count_central = 0;
                multiplicity_count_MFT = 0;
            }
            // Central region
            if (eta > -0.9 && eta < 0.9) {
                multiplicity_count_central++;
            }
            // MFT region
            if (abs(id) == 13) {
                multiplicity_count_MFT++;
            }
        }
        if (multiplicity_count_central > 0) {
            multiplicity_events_central->Fill(multiplicity_count_central); 
        }
        if (multiplicity_count_MFT > 0) {
            multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
        }

        //Double_t scaleBin = 1/(*binLuminocity)[iBin];
        //pThat_bins_mult[iBin]->Scale(scaleBin, "width");
        //multiplicity_events->Add(pThat_bins_mult[iBin]);
    }

    // Output file
    TFile* outFile = new TFile("mymain01Macro.root", "RECREATE");

    TCanvas *canvasEvents = new TCanvas("events_mult","events_mult");

    multiplicity_events_central->Draw();

    multiplicity_events_MFT->SetLineColor(1);
    multiplicity_events_MFT->Draw("SAME");

    auto legendEvents = new TLegend();
    legendEvents->AddEntry(multiplicity_events_central, "Central Barrel Tracks","l");
    legendEvents->AddEntry(multiplicity_events_MFT, "MFT Tracks","l");
    legendEvents->Draw("SAME");

    canvasEvents->Write();

    delete outFile;
}