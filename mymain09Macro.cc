#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include <TNtuple.h>

void mymain09Macro() {
    TFile *infile = TFile::Open("mymain09.root", "READ");

    std::vector<double> *binLuminocity;
    infile->GetObject("luminocity", binLuminocity);

    TH1F *muonPtTotal = new TH1F("muon_full","Produced Muon Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 40, 0.0, 80.0);
    TH1F *muonPtPart = new TH1F("muon_pt_part","", 40, 0.0, 80.0);

    TH1F *muonPtD = new TH1F("muon_pt_D","", 40, 0.0, 80.0);
    TH1F *muonPtB = new TH1F("muon_pt_B","", 40, 0.0, 80.0);
    TH1F *muonPtBD = new TH1F("muon_pt_BD","", 40, 0.0, 80.0);
    TH1F *muonPtBaryon = new TH1F("muon_pt_baryon","", 40, 0.0, 80.0);
    TH1F *muonPtOnium = new TH1F("muon_pt_onium","", 40, 0.0, 80.0);
    TH1F *muonPtTau = new TH1F("muon_pt_tau","", 40, 0.0, 80.0);

    int iBin = 0;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it) {
        TNtuple *muonTuple = (TNtuple*)infile->Get(Form("muon%d", iBin));

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTotal->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 0");
        muonPtPart->Scale(1/(*it), "width");
        muonPtD->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 1");
        muonPtPart->Scale(1/(*it), "width");
        muonPtB->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 2");
        muonPtPart->Scale(1/(*it), "width");
        muonPtBD->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 3 || decayStatus == 4");
        muonPtPart->Scale(1/(*it), "width");
        muonPtBaryon->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 5 || decayStatus == 6");
        muonPtPart->Scale(1/(*it), "width");
        muonPtOnium->Add(muonPtPart);

        muonPtPart->Reset();
        muonTuple->Draw("pt>>muon_pt_part", "decayStatus == 7");
        muonPtPart->Scale(1/(*it), "width");
        muonPtTau->Add(muonPtPart);

        iBin++;
    }

    TFile* outFile = new TFile("mymain09Hist.root", "RECREATE");

    // Decay Status Contributions
    TCanvas *canvasMuon = new TCanvas("Muon_sigma","Muon_sigma");

    muonPtTotal->SetLineColor(1);
    muonPtTotal->Draw();

    muonPtD->SetLineColor(2);
    muonPtD->SetStats(0);
    muonPtD->Draw("SAME");

    muonPtB->SetLineColor(3);
    muonPtB->SetStats(0);
    muonPtB->Draw("SAME");

    muonPtBD->SetLineColor(4);
    muonPtBD->SetStats(0);
    muonPtBD->Draw("SAME");

    muonPtBaryon->SetLineColor(5);
    muonPtBaryon->SetStats(0);
    muonPtBaryon->Draw("SAME");

    muonPtOnium->SetLineColor(6);
    muonPtOnium->SetStats(0);
    muonPtOnium->Draw("SAME");

    muonPtTau->SetLineColor(7);
    muonPtTau->SetStats(0);
    muonPtTau->Draw("SAME");

    auto legendMuon = new TLegend();
    legendMuon->AddEntry(muonPtTotal,"Inclusive #mu","l");
    legendMuon->AddEntry(muonPtD,"c->D->#mu","l");
    legendMuon->AddEntry(muonPtB,"b->B->#mu","l");
    legendMuon->AddEntry(muonPtBD,"b->B->D->#mu","l");
    legendMuon->AddEntry(muonPtBaryon,"?->Baryon->#mu","l");
    legendMuon->AddEntry(muonPtOnium,"?->Quarkonium->#mu","l");
    legendMuon->AddEntry(muonPtTau,"?->#tau->#mu","l");
    legendMuon->Draw("SAME");

    canvasMuon->Write();

    delete outFile;
}