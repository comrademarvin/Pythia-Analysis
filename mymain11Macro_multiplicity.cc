#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

// void processWFile(string filename) {

// }

void mymain11Macro_multiplicity() {
    const Int_t multBinCount = 6;
    Double_t multBins[multBinCount+1] = {0, 10, 20, 30, 40, 50, 60};

    const Int_t pTBinCount = 20;

    // Access Minimum Bias data
    TFile* infile_mb = TFile::Open("mymain01Macro.root", "READ");
    TH1F* mb_mult_central = (TH1F*) infile_mb->Get("multiplicity_events_central");

    // Access W+/- showered data
    TFile* infile = TFile::Open("mymain11_W+_500k.root", "READ");

    // W->muon distributions (pt+y, event multiplicities, multiplicity bins, yield)
    TH2F* W_muon_pt_y = new TH2F("W_muon_pt_y", "W->Muon parameter space for p_{T} and rapidity;p_{T} (GeV/c);y;#sigma_{W->#mu} (mb)", 100, 0, 200, 100, -10, 10);
    TH1F* W_muon_multiplicity = new TH1F("W_muon_multiplicity", "Primary charged particle dependence;dN_{ch}/d#eta_{|#eta|<1};#sigma_{W->#mu} (mb)", multBinCount, multBins);

    vector<TH1F*> W_muon_pt_mult_binned(multBinCount);
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        W_muon_pt_mult_binned[iBin] = new TH1F(Form("W_muon_pt_%d", iBin), "Muon Yield; p_{T}; N", pTBinCount, 10, 60);
    }

    // Read in generated cross-section as integrated luminocity
    std::vector<double> *luminocity;
    infile->GetObject("luminocity", luminocity);

    // Read in multiplicities of events
    std::vector<int> *event_mult;
    infile->GetObject("event_multiplicity", event_mult);
    int event_count = event_mult->size();

    // iterate through muons and connect to event multiplicities
    // For W+ -> muon
    TNtuple *muonTuple = (TNtuple*)infile->Get("W_muon");
    int muons_count = muonTuple->GetEntries();
    Float_t eventIndex, eta, pt, rapidity;
    muonTuple->SetBranchAddress("eventTag", &eventIndex);
    muonTuple->SetBranchAddress("y", &rapidity);
    muonTuple->SetBranchAddress("eta", &eta);
    muonTuple->SetBranchAddress("pt", &pt);

    for (int iMuon = 0; iMuon < muons_count; iMuon++) {
        muonTuple->GetEntry(iMuon);
        W_muon_pt_y->Fill(pt, rapidity);
        int muon_event_mult = (*event_mult)[static_cast<int>(eventIndex)];
        W_muon_multiplicity->Fill(muon_event_mult/2);

        if ((rapidity >= 2.5) && (rapidity <= 4)) {
            for (int iBin = 0; iBin < multBinCount; iBin++) {
                if ((muon_event_mult > multBins[iBin]) && (muon_event_mult < multBins[iBin+1])) {
                    W_muon_pt_mult_binned[iBin]->Fill(pt);
                }
            }
        }
    };

    // Cross Section Scaling
    W_muon_multiplicity->Scale(1/((*luminocity)[0]), "width");
    W_muon_pt_y->Scale(1/((*luminocity)[0]), "width");

    // Compute the multiplicity dependant yield

    // Cross section ratio multiplicity binned
    TH1F* cs_ratio = new TH1F();

    *cs_ratio = (*mb_mult_central)/(*W_muon_multiplicity);

    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        W_muon_pt_mult_binned[iMultBin]->Scale(1/(cs_ratio->GetBinContent(iMultBin+1)));
    }

    // Determine yield of W->muon as function of multiplicity
    // TODO: Scale event count by minimum bias
    TProfile* multProfile = new TProfile("mult_prof","Profile of yield versus multiplicity", multBinCount, multBins);

    for (int iMultBin = 0; iMultBin < multBinCount; iMultBin++) {
        double multBinAverage = multBins[iMultBin] + ((multBins[iMultBin+1] - multBins[iMultBin])/2);
        std::cout << multBinAverage << std::endl;
        for (int iPtBin = 1; iPtBin <= pTBinCount; iPtBin++) {
            double yield_pt_bin = (W_muon_pt_mult_binned[iMultBin])->GetBinContent(iPtBin);
            multProfile->Fill(multBinAverage, yield_pt_bin, 1);
        }
    }

    // Output file
    TFile* outFile = new TFile("mymain11Hist_multiplicity.root", "RECREATE");

    TCanvas *canvasMuonPtY = new TCanvas("Muon_pt_y","Muon_pt_y");
    W_muon_pt_y->Draw("colz");
    canvasMuonPtY->Write();

    TCanvas *canvasMuonMult = new TCanvas("Muon_mult","Muon_mult");
    gPad->SetLogy();
    mb_mult_central->SetMinimum(0.00000000001);
    mb_mult_central->SetLineColor(1);
    mb_mult_central->Draw("SAME");
    W_muon_multiplicity->SetMinimum(0.00000000001);
    W_muon_multiplicity->Draw("SAME");
    auto legendMult = new TLegend();
    legendMult->AddEntry(mb_mult_central,"Minimum Bias","l");
    legendMult->AddEntry(W_muon_multiplicity,"W->#mu","l");
    legendMult->Draw("SAME");
    canvasMuonMult->Write();

    TCanvas *canvasMultBinned = new TCanvas("mult_binned_ratio","mult_binned_ratio");
    gPad->SetLogy();
    cs_ratio->SetStats(0);
    //cs_ratio->SetMinimum(0.00000001);
    cs_ratio->SetTitle("Ratio of W muon Cross Section wrt Minimum Bias");
    cs_ratio->GetYaxis()->SetTitle("#sigma_{MB}/#sigma_{W->#mu}");
    cs_ratio->Draw();
    // mb_mult_central_binned->SetMinimum(0.000000000001);
    // mb_mult_central_binned->Draw("SAME");
    // W_muon_mult_binned->SetMinimum(0.000000000001);
    // W_muon_mult_binned->Draw("SAME");
    canvasMultBinned->Write();

    TCanvas *canvasYieldBinned = new TCanvas("pt_yield_mult_binned","pt_yield_mult_binned");
    gPad->SetLogy();
    auto legendPtBins = new TLegend();
    for (int iBin = 0; iBin < multBinCount; iBin++) {
        W_muon_pt_mult_binned[iBin]->SetLineColor(iBin);
        W_muon_pt_mult_binned[iBin]->Draw("SAME");
        legendPtBins->AddEntry(W_muon_pt_mult_binned[iBin], Form("%f < N_{ch} < %f", multBins[iBin], multBins[iBin+1]),"l");
    }
    // auto labelCuts = new TLatex();
    // labelCuts->DrawLatex(0.0, 0.0, "2.5 < y < 4");
    // labelCuts->Draw("SAME");
    canvasYieldBinned->Write();

    TCanvas *canvasYieldMultDep = new TCanvas("yield_mult_dep","yield_mult_dep");
    multProfile->Draw();
    canvasYieldMultDep->Write();

    delete outFile;
};