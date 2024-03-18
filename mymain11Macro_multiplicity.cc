#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain11Macro_multiplicity() {
    TFile *infile = TFile::Open("mymain11_W+_100k.root", "READ");

    // Define multiplicity bins
    const int mult_bin_count = 4;
    const int multiplicity_bin_bounds[mult_bin_count+1] = {1, 20, 40, 60, 90};

    // W->muon event multiplicity
    TH1F* W_muon_multiplicity = new TH1F("", "Multiplicity of W->Muon events;N_{ch};#sigma_{W->#mu} (mb)", 50, 0, 100);

    // Read in generated "luminocity"
    std::vector<double> *luminocity;
    infile->GetObject("luminocity", luminocity);

    // Read in multiplicities of events
    std::vector<int> *event_mult;
    infile->GetObject("event_multiplicity", event_mult);

    // Event multiplicity distributions
    // vector<TH1F*> mult_bin_distributions(mult_bin_count);

    // for (int* iMult = event_mult->begin(); iMult != event_mult->end(); ++iMult) {
        
    // }

    // Read in tuple of W->muon decays
    TNtuple *muonTuple = (TNtuple*)infile->Get("W_muon");

    // iterate through muons and connect to event multiplicities
    int muons_count = muonTuple->GetEntries();
    Float_t eventIndex, eta;
    muonTuple->SetBranchAddress("eventTag", &eventIndex);
    muonTuple->SetBranchAddress("eta", &eta);

    for (int iMuon = 0; iMuon < muons_count; iMuon++) {
        muonTuple->GetEntry(iMuon);
        if ((eta >= 2.5) && (eta <= 4)) {
            int muon_event_mult = (*event_mult)[static_cast<int>(eventIndex)];
            W_muon_multiplicity->Fill(muon_event_mult);
        }
    };

    W_muon_multiplicity->Scale(1/((*luminocity)[0]), "width");

    // Output file
    TFile* outFile = new TFile("mymain11Hist_multiplicity.root", "RECREATE");

    TCanvas *canvasMuonMult = new TCanvas("Muon_mult","Muon_mult");
    W_muon_multiplicity->Draw();
    canvasMuonMult->Write();

    delete outFile;
};