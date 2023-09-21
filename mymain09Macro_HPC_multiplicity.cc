#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain09Macro_HPC_multiplicity() {
    // pT-hat bins
    const int pThat_bins_count = 7;
    double pThat_bins[pThat_bins_count + 1] = {0.0, 10.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0};

    // pT bins
    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};

    // multiplicity bins
    // const int mult_bins_count = 4;
    // int mult_bins[mult_bins_count+1] = {0, 25, 50, 75, 100};

    // vector<TH1F*> muonPt_bins(mult_bins_count);

    // for (int iBin = 0; iBin < mult_bins_count; iBin++) {
    //     muonPt_bins[iBin] = new TH1F(Form("muon_pt_mult_bin_%d", iBin), "", NBINS, edges);
    // }

    // pThat multiplicity
    TH1F* multiplicity_events = new TH1F("", "Multiplicity of all events;N_{ch};N/N_{events}", 300, 0, 300);

    vector<TH1F*> pThat_bins_mult(pThat_bins_count);

    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult[iBin] = new TH1F(Form("pThat_bin_%d", iBin), "", 300, 0, 300);
    }

    // muon multiplicity
    TH1F* multiplicity_total = new TH1F("", "Multiplicity of HF Muon events;N_{ch};N", 150, 0, 300);
    //TH1F* multiplicity_part = new TH1F("", "", 150, 0, 300);

    for(int iBin = 0; iBin < pThat_bins_count; iBin++) {
        //multiplicity_part->Reset();

        TFile *infile = TFile::Open(Form("HPC_4M_PDF9/mymain09_%d.root", iBin), "READ");

        std::vector<double> *binLuminocity;
        infile->GetObject("luminocity", binLuminocity);

        std::vector<int> *mult_events;
        infile->GetObject("multiplicity", mult_events);

        // event multiplicities
        for(std::vector<int>::iterator mult = mult_events->begin(); mult != mult_events->end(); ++mult) {
            if (*mult != -1) {
                pThat_bins_mult[iBin]->Fill(*mult);
            }
        }

        Double_t scale = pThat_bins_mult[iBin]->GetXaxis()->GetBinWidth(1)/(pThat_bins_mult[iBin]->Integral());
        pThat_bins_mult[iBin]->Scale(scale);

        multiplicity_events->Add(pThat_bins_mult[iBin]);

        // muon multiplicities
        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        // iterate through muons and connect to event multiplicities
        int events_count = muonTuple->GetEntries();
        Float_t eventIndex, decayStatus;
        muonTuple->SetBranchAddress("eventTag", &eventIndex);
        muonTuple->SetBranchAddress("decayStatus", &decayStatus);

        for (int iMuon = 0; iMuon < events_count; iMuon++) {
            muonTuple->GetEntry(iMuon);
            int event_mult = (*mult_events)[static_cast<int>(eventIndex)];
            int muon_decay_type = static_cast<int>(decayStatus);

            if (muon_decay_type == 0 || muon_decay_type == 1 || muon_decay_type == 2) {
                multiplicity_total->Fill(event_mult);
            }
        }

        // multiplicity_part->Scale(1/(*binLuminocity)[iBin], "width");
        // multiplicity_total->Add(multiplicity_part);
    }

    // Output file
    TFile* outFile = new TFile("mymain09Hist_HPC_multiplicity.root", "RECREATE");

    // Events mult
    TCanvas *canvasEvents = new TCanvas("Events_mult","Events_mult");

    multiplicity_events->Draw();

    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult[iBin]->Draw("SAME");
    }

    canvasEvents->Write();

    // Muon mult
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    multiplicity_total->Draw();

    canvasMuon->Write();

    delete outFile;
}