#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain11Macro_multiplicity() {
    // Access Minimum Bias data
    

    // Access W+/- showered data
    TFile *infile_plus = TFile::Open("mymain11_W+_500k.root", "READ");
    TFile *infile_minus = TFile::Open("mymain11_W-_500k.root", "READ");

    // // Define multiplicity bins
    // const int mult_bin_count = 4;
    // const int multiplicity_bin_bounds[mult_bin_count+1] = {1, 20, 40, 60, 90};

    // W->muon distributions (pt+y, event multiplicities, multiplicity bins)
    TH2F* W_muon_pt_y = new TH2F("W_muon_pt_y", "W->Muon parameter space for p_{T} and rapidity;p_{T} (GeV/c);y;#sigma_{W->#mu} (mb)", 100, 0, 200, 100, -10, 10);
    TH2F* W_muon_pt_y_part = new TH2F("W_muon_pt_y_part", "", 100, 0, 200, 100, -10, 10);

    TH1F* W_muon_multiplicity = new TH1F("W_muon_multiplicity", "Multiplicity of W->Muon events;N_{ch};#sigma_{W->#mu} (mb)", 30, 0, 60);
    TH1F* W_muon_multiplicity_part = new TH1F("W_muon_multiplicity_part", "", 30, 0, 60);
    // vector<TH1F*> W_muon_multiplicity_bins(mult_bin_count);

    // for (int iBin = 0; iBin < mult_bin_count; iBin++) {
    //     W_muon_multiplicity_bins[iBin] = new TH1F(Form("W_muon_mult_bin_%d", iBin), "", 50, 0, 100);
    // }

    // Read in generated cross-section as integrated luminocity
    std::vector<double> *luminocity_plus;
    infile_plus->GetObject("luminocity", luminocity_plus);

    std::vector<double> *luminocity_minus;
    infile_minus->GetObject("luminocity", luminocity_minus);

    // Read in multiplicities of events
    std::vector<int> *event_mult_plus;
    infile_plus->GetObject("event_multiplicity", event_mult_plus);
    int event_count_plus = event_mult_plus->size();

    std::vector<int> *event_mult_minus;
    infile_minus->GetObject("event_multiplicity", event_mult_minus);
    int event_count_minus = event_mult_minus->size();

    // iterate through muons and connect to event multiplicities
    // For W+ -> muon
    TNtuple *muonTuple_plus = (TNtuple*)infile_plus->Get("W_muon");
    int muons_count = muonTuple_plus->GetEntries();
    Float_t eventIndex, eta, pt, rapidity;
    muonTuple_plus->SetBranchAddress("eventTag", &eventIndex);
    muonTuple_plus->SetBranchAddress("y", &rapidity);
    muonTuple_plus->SetBranchAddress("eta", &eta);
    muonTuple_plus->SetBranchAddress("pt", &pt);

    for (int iMuon = 0; iMuon < muons_count; iMuon++) {
        muonTuple_plus->GetEntry(iMuon);
        W_muon_pt_y_part->Fill(pt, rapidity);
        int muon_event_mult = (*event_mult_plus)[static_cast<int>(eventIndex)];
        W_muon_multiplicity_part->Fill(muon_event_mult/2);
        // if ((eta >= 2.5) && (eta <= 4)) {
        //     int muon_event_mult = (*event_mult)[static_cast<int>(eventIndex)];
        //     W_muon_multiplicity->Fill(muon_event_mult);

        //     // bin in multiplicity classes
        //     for (int iBin = 0; iBin < mult_bin_count; iBin++) {
        //         if ((muon_event_mult > multiplicity_bin_bounds[iBin]) && (muon_event_mult <= multiplicity_bin_bounds[iBin+1])) {
        //             W_muon_multiplicity_bins[iBin]->Fill(pt);
        //             break;
        //         }
        //     }
        // }
    };

    W_muon_multiplicity_part->Scale(1/((*luminocity_plus)[0]), "width");
    W_muon_multiplicity->Add(W_muon_multiplicity_part);

    W_muon_pt_y_part->Scale(1/((*luminocity_plus)[0]), "width");
    W_muon_pt_y->Add(W_muon_pt_y_part);

    // For W- -> muon
    TNtuple *muonTuple_minus = (TNtuple*)infile_minus->Get("W_muon");
    muons_count = muonTuple_minus->GetEntries();
    muonTuple_minus->SetBranchAddress("eventTag", &eventIndex);
    muonTuple_minus->SetBranchAddress("y", &rapidity);
    muonTuple_minus->SetBranchAddress("eta", &eta);
    muonTuple_minus->SetBranchAddress("pt", &pt);
    W_muon_multiplicity_part->Reset();
    W_muon_pt_y_part->Reset();

    for (int iMuon = 0; iMuon < muons_count; iMuon++) {
        muonTuple_minus->GetEntry(iMuon);
        W_muon_pt_y_part->Fill(pt, rapidity);
        int muon_event_mult = (*event_mult_minus)[static_cast<int>(eventIndex)];
        W_muon_multiplicity_part->Fill(muon_event_mult/2);
    };

    W_muon_multiplicity_part->Scale(1/((*luminocity_minus)[0]), "width");
    W_muon_multiplicity->Add(W_muon_multiplicity_part);

    W_muon_pt_y_part->Scale(1/((*luminocity_minus)[0]), "width");
    W_muon_pt_y->Add(W_muon_pt_y_part);

    // Output file
    TFile* outFile = new TFile("mymain11Hist_multiplicity.root", "RECREATE");

    TCanvas *canvasMuonPtY = new TCanvas("Muon_pt_y","Muon_pt_y");
    W_muon_pt_y->Draw("colz");
    canvasMuonPtY->Write();

    TCanvas *canvasMuonMult = new TCanvas("Muon_mult","Muon_mult");
    W_muon_multiplicity->Draw();
    canvasMuonMult->Write();

    // TCanvas *canvasMuonPt = new TCanvas("Muon_pt","Muon_pt");

    // auto legendMuonPt = new TLegend();
    // for (int iBin = 0; iBin < mult_bin_count; iBin++) {
    //     W_muon_multiplicity_bins[iBin]->Scale(1/((*luminocity)[0]), "width");
    //     W_muon_multiplicity_bins[iBin]->SetLineColor(iBin+1);
    //     W_muon_multiplicity_bins[iBin]->Draw("SAME");
    //     legendMuonPt->AddEntry(W_muon_multiplicity_bins[iBin], Form("%d < N_{ch} <= %d", static_cast<int>(multiplicity_bin_bounds[iBin]), static_cast<int>(multiplicity_bin_bounds[iBin+1])),"l");
    // }
    // legendMuonPt->Draw("SAME");

    // canvasMuonPt->Write();

    delete outFile;
};