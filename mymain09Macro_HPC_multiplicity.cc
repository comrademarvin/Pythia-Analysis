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

    // HF muon pt multiplicity bins
    const int mult_bins_count = 5;
    int mult_bins[mult_bins_count+1] = {0, 50, 100, 150, 200, 300};

    vector<TH1F*> muonPt_bins(mult_bins_count);
    vector<TH1F*> muonPt_bins_temp(mult_bins_count);

    for (int iMultBin = 0; iMultBin < mult_bins_count; iMultBin++) {
        muonPt_bins[iMultBin] = new TH1F(Form("muon_pt_mult_bin_%d", iMultBin), "", NBINS, edges);
        muonPt_bins_temp[iMultBin] = new TH1F(Form("muon_pt_mult_bin_temp_%d", iMultBin), "", NBINS, edges);
    }

    // pThat multiplicity
    TH1F* multiplicity_events = new TH1F("", "Multiplicity of all events;N_{ch};#sigma (mb)", 150, 0, 300);

    vector<TH1F*> pThat_bins_mult(pThat_bins_count);

    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult[iBin] = new TH1F(Form("pThat_bin_%d", iBin), "", 150, 0, 300);
    }

    // muon multiplicity
    TH1F* multiplicity_total = new TH1F("", "Multiplicity of HF Muon events;N_{ch};#frac{1}{N_{events}}N_{HF->#mu}", 150, 0, 300);
    
    vector<TH1F*> pThat_bins_mult_muon(pThat_bins_count);

    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult_muon[iBin] = new TH1F(Form("pThat_muon_bin_%d", iBin), "", 150, 0, 300);
    }

    for(int iBin = 0; iBin < pThat_bins_count; iBin++) {
        //multiplicity_part->Reset();

        TFile *infile = TFile::Open(Form("HPC_4M_PDF9/mymain09_%d.root", iBin), "READ");

        std::vector<double> *binLuminocity;
        infile->GetObject("luminocity", binLuminocity);

        std::vector<int> *mult_events;
        infile->GetObject("multiplicity", mult_events);

        // event multiplicities
        for(std::vector<int>::iterator mult = mult_events->begin(); mult != mult_events->end(); ++mult) {
            if (*mult > 0) {
                pThat_bins_mult[iBin]->Fill(*mult);
            }
        }

        //Double_t scale = pThat_bins_mult[iBin]->GetXaxis()->GetBinWidth(1)/(pThat_bins_mult[iBin]->Integral());
        Double_t scale = pow(10,-9)/(*binLuminocity)[iBin];
        pThat_bins_mult[iBin]->Scale(scale);

        multiplicity_events->Add(pThat_bins_mult[iBin]);

        // muon multiplicities
        TNtuple *muonTuple = (TNtuple*)infile->Get("muon");

        // iterate through muons and connect to event multiplicities
        int events_count = muonTuple->GetEntries();
        Float_t eventIndex, decayStatus, pt;
        muonTuple->SetBranchAddress("eventTag", &eventIndex);
        muonTuple->SetBranchAddress("decayStatus", &decayStatus);
        muonTuple->SetBranchAddress("pt", &pt);

        for (int iMuon = 0; iMuon < events_count; iMuon++) {
            muonTuple->GetEntry(iMuon);
            int event_mult = (*mult_events)[static_cast<int>(eventIndex)];
            int muon_decay_type = static_cast<int>(decayStatus);

            if (muon_decay_type == 0 || muon_decay_type == 1 || muon_decay_type == 2) {
                if (event_mult > 1) {
                    pThat_bins_mult_muon[iBin]->Fill(event_mult);

                    for (int iMultBin = 0; iMultBin < mult_bins_count; iMultBin++) {
                        if ((event_mult>=mult_bins[iMultBin]) && (event_mult<mult_bins[iMultBin+1])) {
                            muonPt_bins_temp[iMultBin]->Fill(pt);
                        }
                    }
                }
            }
        }

        //Double_t scaleMuon = pThat_bins_mult_muon[iBin]->GetXaxis()->GetBinWidth(1)/(pThat_bins_mult_muon[iBin]->Integral());
        Double_t scaleMuon = pow(10,-9)/(*binLuminocity)[iBin];
        pThat_bins_mult_muon[iBin]->Scale(scaleMuon);
        multiplicity_total->Add(pThat_bins_mult_muon[iBin]);

        Double_t scalePt = 1/(*binLuminocity)[iBin];
        for (int iMultBin = 0; iMultBin < mult_bins_count; iMultBin++) {
            muonPt_bins_temp[iMultBin]->Scale(scalePt);
            muonPt_bins[iMultBin]->Add(muonPt_bins_temp[iMultBin]);
        }
    }

    // Output file
    TFile* outFile = new TFile("mymain09Hist_HPC_multiplicity.root", "RECREATE");

    // Events mult
    TCanvas *canvasEvents = new TCanvas("Events_mult","Events_mult");

    multiplicity_events->Draw();

    auto legendEvents = new TLegend();
    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult[iBin]->SetLineColor(iBin+1);
        pThat_bins_mult[iBin]->Draw("SAME");
        legendEvents->AddEntry(pThat_bins_mult[iBin], Form("%d < #hat{p}_{T} < %d", static_cast<int>(pThat_bins[iBin]), static_cast<int>(pThat_bins[iBin+1])),"l");
    }
    legendEvents->Draw("SAME");

    canvasEvents->Write();

    // Muon mult
    TCanvas *canvasMuon = new TCanvas("Muon_mult","Muon_mult");

    multiplicity_total->Draw();

    auto legendMuon = new TLegend();
    for (int iBin = 0; iBin < pThat_bins_count; iBin++) {
        pThat_bins_mult_muon[iBin]->SetLineColor(iBin+1);
        pThat_bins_mult_muon[iBin]->Draw("SAME");
        legendMuon->AddEntry(pThat_bins_mult_muon[iBin], Form("%d < #hat{p}_{T} < %d", static_cast<int>(pThat_bins[iBin]), static_cast<int>(pThat_bins[iBin+1])),"l");
    }
    legendMuon->Draw("SAME");

    canvasMuon->Write();

    // HF muon pt multiplicity bins
    TCanvas *canvasMultPt = new TCanvas("HF_Muon_mult_bins","HF_Muon_mult_bins");

    auto legendMuonPt = new TLegend();
    for (int iMultBin = 0; iMultBin < mult_bins_count; iMultBin++) {
        muonPt_bins[iMultBin]->SetLineColor(iMultBin+1);
        muonPt_bins[iMultBin]->Draw("SAME");
        legendMuonPt->AddEntry(muonPt_bins[iMultBin], Form("%d < N_{ch} < %d", static_cast<int>(mult_bins[iMultBin]), static_cast<int>(mult_bins[iMultBin+1])),"l");
    }
    legendMuonPt->Draw("SAME");

    canvasMultPt->Write();

    delete outFile;
}