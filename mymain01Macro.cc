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
    const Int_t multBinCount = 7;
    Double_t multBins[multBinCount+1] = {1, 5, 10, 15, 20, 25, 30, 35};

    // pThat multiplicity
    TH1D* multiplicity_events_central = new TH1D("multiplicity_events_central", "Integrated Primary Charged Particle Multiplicity Dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma (mb)", multBinCount, multBins);
    TH1D* multiplicity_events_central_part = new TH1D("multiplicity_events_central_part", "", multBinCount, multBins);

    TH1D* multiplicity_events_central_raw = new TH1D("multiplicity_events_central_raw", "Primary Charged Particle Multiplicity Dependence;N_{ch};#sigma (mb)", 50, 0, 100);
    TH1D* multiplicity_events_central_raw_part = new TH1D("multiplicity_events_central_raw_part", "", 50, 0, 100);

    //TH1F* multiplicity_events_MFT = new TH1F("multiplicity_events_MFT", "", 140, 0, 140);
    //TH1F* multiplicity_events_Muon = new TH1F("multiplicity_events_Muon", "", 140, 0, 140);

    //TH1F* pt_test = new TH1F("pt_central", "", 100, 0, 50);

    // vector<TH1F*> pThat_bins_mult(nBins);

    // for (int iBin = 0; iBin < nBins; iBin++) {
    //     pThat_bins_mult[iBin] = new TH1F(Form("pThat_bin_%d", iBin), "", 150, 0, 300);
    // }

    TFile *infile = TFile::Open("mymain01_1M_536.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    for(int iBin = 0; iBin < nBins; iBin++) {
        TNtuple *chargedTuple = (TNtuple*)infile->Get(Form("charged%d", iBin));

        multiplicity_events_central_part->Reset();
        multiplicity_events_central_raw_part->Reset();

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
        // int multiplicity_count_MFT = 0;
        // int multiplicity_count_Muon = 0;
        for (int iParticle = 0; iParticle < particle_count; iParticle++) {
            chargedTuple->GetEntry(iParticle);
            if (eventIndex != event_counter) {
                if (multiplicity_count_central > 0) {
                    multiplicity_events_central_part->Fill(multiplicity_count_central/2);
                    multiplicity_events_central_raw_part->Fill(multiplicity_count_central);
                }
                // if (multiplicity_count_MFT > 0) {
                //     multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
                // }
                // if (multiplicity_count_Muon > 0) {
                //     multiplicity_events_Muon->Fill(multiplicity_count_Muon); 
                // }
                event_counter++;
                multiplicity_count_central = 0;
                // multiplicity_count_MFT = 0;
                // multiplicity_count_Muon = 0;
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
        }
        // if (multiplicity_count_MFT > 0) {
        //     multiplicity_events_MFT->Fill(multiplicity_count_MFT); 
        // }
        // if (multiplicity_count_Muon > 0) {
        //     multiplicity_events_Muon->Fill(multiplicity_count_Muon); 
        // }

        Double_t scaleBin = 1/(*binLuminocity)[iBin];
        //pThat_bins_mult[iBin]->Scale(scaleBin, "width");
        //multiplicity_events->Add(pThat_bins_mult[iBin]);
        multiplicity_events_central_part->Scale(scaleBin, "width");
        multiplicity_events_central->Add(multiplicity_events_central_part);

        multiplicity_events_central_raw_part->Scale(scaleBin, "width");
        multiplicity_events_central_raw->Add(multiplicity_events_central_raw_part);
        // multiplicity_events_MFT->Scale(scaleBin, "width");
        // multiplicity_events_Muon->Scale(scaleBin, "width");
    }

    std::cout << "<dN_ch/d_eta> for |eta| < 1: " << multiplicity_events_central->GetMean() << std::endl;
    std::cout << "RMS of dN_ch/d_eta for |eta| < 1: " << multiplicity_events_central->GetRMS() << std::endl;

    // Output file
    TFile* outFile = new TFile("mymain01Macro.root", "RECREATE");

    multiplicity_events_central->Draw();
    multiplicity_events_central->Write();

    multiplicity_events_central_raw->Draw();
    multiplicity_events_central_raw->Write();

    delete outFile;

    //TCanvas *canvasEvents = new TCanvas("events_mult","events_mult");

    //multiplicity_events_central->Draw();

    // multiplicity_events_MFT->SetLineColor(1);
    // multiplicity_events_MFT->Draw("SAME");

    // multiplicity_events_Muon->SetLineColor(2);
    // multiplicity_events_Muon->Draw("SAME");

    // auto legendEvents = new TLegend();
    // legendEvents->AddEntry(multiplicity_events_central, "Central Barrel (#eta < |0.9|)","l");
    // legendEvents->AddEntry(multiplicity_events_MFT, "MFT (-3.6 < #eta < -2.5)","l");
    // legendEvents->AddEntry(multiplicity_events_Muon, "Muon Spectrometer (#mu's only & -4 < #eta < -2.5)","l");
    // legendEvents->Draw("SAME");

    // auto labelCuts = new TLatex();
    // labelCuts->DrawLatex(0.0, 0.0, "pp @ #sqrt{s} = 13.6 TeV");
    // labelCuts->DrawLatex(0.0, 0.0, "Pythia8 Monash Tune (Tune:pp = 14)");
    // labelCuts->DrawLatex(0.0, 0.0, "N_{ch} > 0");
    // labelCuts->Draw("SAME");

    // pt
    // pt_test->Draw();
    // pt_test->Write();

    //canvasEvents->Write();

    // multiplicity_events_MFT->Draw();
    // multiplicity_events_MFT->Write();

    // multiplicity_events_Muon->Draw();
    // multiplicity_events_Muon->Write();
}