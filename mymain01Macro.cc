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
    TH1F* multiplicity_events_central = new TH1F("multiplicity_events_central", "Charged Particle Multiplicity for Detector Regions;multiplicity (N_{ch});N_{events}", 140, 0, 140);
    TH1F* multiplicity_events_MFT = new TH1F("multiplicity_events_MFT", "", 140, 0, 140);
    TH1F* multiplicity_events_Muon = new TH1F("multiplicity_events_Muon", "", 140, 0, 140);

    TH1F* pt_test = new TH1F("pt_central", "", 100, 0, 50);

    // vector<TH1F*> pThat_bins_mult(pThat_bins_count);

    // for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
    //     pThat_bins_mult[iBin] = new TH1F(Form("pThat_bin_%d", iBin), "", 150, 0, 300);
    // }

    TFile *infile = TFile::Open("/mnt/d/Pythia_Results/mymain01_13TeV_3_5M.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    for(int iBin = 0; iBin < pThat_bins_count; iBin++) {
        TNtuple *chargedTuple = (TNtuple*)infile->Get(Form("charged%d", iBin));

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
        int multiplicity_count_MFT = 0;
        int multiplicity_count_Muon = 0;
        for (int iParticle = 0; iParticle < particle_count; iParticle++) {
            chargedTuple->GetEntry(iParticle);
            if (eventIndex != event_counter) {
                if (multiplicity_count_central > 0) {
                    multiplicity_events_central->Fill(multiplicity_count_central); 
                }
                if (multiplicity_count_MFT > 0) {
                    multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
                }
                if (multiplicity_count_Muon > 0) {
                    multiplicity_events_Muon->Fill(multiplicity_count_Muon); 
                }
                event_counter++;
                multiplicity_count_central = 0;
                multiplicity_count_MFT = 0;
                multiplicity_count_Muon = 0;
            }
            // Central region
            if (eta > -0.9 && eta < 0.9) {
                multiplicity_count_central++;
            }
            // MFT region
            if (eta > -3.6 && eta < -2.5) {
                multiplicity_count_MFT++;
                pt_test->Fill(pt);
            }
            // Forward region (muon)
            if (abs(id) == 13 && eta > -4.0 && eta < -2.5) {
                multiplicity_count_Muon++;
            }
        }
        if (multiplicity_count_central > 0) {
            multiplicity_events_central->Fill(multiplicity_count_central); 
        }
        if (multiplicity_count_MFT > 0) {
            multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
        }
        if (multiplicity_count_Muon > 0) {
            multiplicity_events_Muon->Fill(multiplicity_count_Muon); 
        }

        //Double_t scaleBin = 1/(*binLuminocity)[iBin];
        //pThat_bins_mult[iBin]->Scale(scaleBin, "width");
        //multiplicity_events->Add(pThat_bins_mult[iBin]);
        // multiplicity_events_central->Scale(scaleBin, "width");
        // multiplicity_events_MFT->Scale(scaleBin, "width");
        // multiplicity_events_Muon->Scale(scaleBin, "width");
    }

    // Output file
    TFile* outFile = new TFile("mymain01Macro.root", "RECREATE");

    TCanvas *canvasEvents = new TCanvas("events_mult","events_mult");

    multiplicity_events_central->Draw();

    multiplicity_events_MFT->SetLineColor(1);
    multiplicity_events_MFT->Draw("SAME");

    multiplicity_events_Muon->SetLineColor(2);
    multiplicity_events_Muon->Draw("SAME");

    auto legendEvents = new TLegend();
    legendEvents->AddEntry(multiplicity_events_central, "Central Barrel (#eta < |0.9|)","l");
    legendEvents->AddEntry(multiplicity_events_MFT, "MFT (-3.6 < #eta < -2.5)","l");
    legendEvents->AddEntry(multiplicity_events_Muon, "Muon Spectrometer (#mu's only & -4 < #eta < -2.5)","l");
    legendEvents->Draw("SAME");

    auto labelCuts = new TLatex();
    labelCuts->DrawLatex(0.0, 0.0, "pp @ #sqrt{s} = 13.6 TeV");
    labelCuts->DrawLatex(0.0, 0.0, "Pythia8 Monash Tune (Tune:pp = 14)");
    labelCuts->DrawLatex(0.0, 0.0, "N_{ch} > 0");
    labelCuts->Draw("SAME");

    // pt
    // pt_test->Draw();
    // pt_test->Write();

    canvasEvents->Write();

    multiplicity_events_central->Draw();
    multiplicity_events_central->Write();

    multiplicity_events_MFT->Draw();
    multiplicity_events_MFT->Write();

    multiplicity_events_Muon->Draw();
    multiplicity_events_Muon->Write();

    delete outFile;
}