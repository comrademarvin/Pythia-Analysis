#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

int main() {
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain06.root", "RECREATE");

    // pTHat bins
    int nBins = 4;
    double binEdges[nBins+1] = {16.0, 26.0, 36.0, 50.0, 70.0};
    double binLuminocity[nBins];

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 70, 0.0, 70.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 70, 0.0, 70.0);

    // HF Cross Sections
    vector<TNtuple*> charmTuples(nBins);
    vector<TNtuple*> bottomTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        charmTuples[i] = new TNtuple("charm", "charm", "pt:status");
        bottomTuples[i] = new TNtuple("bottom", "bottom", "pt:status");
    }

    // Number of events to generate per bin.
    int N_events = 1000;

    for (int iBin = 0; iBin < nBins; ++iBin) {
        // set pythia initialization variables
        pythia.readString("HardQCD:all = on");
        // pythia.readString("HardQCD:hardccbar = on");
        // pythia.readString("HardQCD:hardbbbar = on");
        pythia.readString("SoftQCD:nonDiffractive = off");

        pythia.readString("Beams:eCM = 5020.");
        pythia.readString("Tune:pp = 14");
        pythia.readString("411:onMode=off");
        pythia.readString("411:onIfAny=13");
        pythia.readString("411:onMode=off");
        pythia.readString("411:onIfAny=13");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat  = pythia.info.pTHat();
            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);

            for (int i = 0; i < pythia.event.size(); ++i) {
                int particleID = abs(pythia.event[i].id());
                int particleStatus = abs(pythia.event[i].status());
                int particlePt = pythia.event[i].pT();

                if (particleID == 4) { // charm
                    charmTuples[iBin]->Fill(particlePt, particleStatus);
                }

                if (particleID == 5) { // bottom
                    bottomTuples[iBin]->Fill(particlePt, particleStatus);
                }
            }
        }

        // cross-section for the bin
        double luminocity_hard = N_events/(pythia.info.sigmaGen()*pow(10,9));
        binLuminocity[iBin] = luminocity_hard;

        hardPtPart->Scale(1/luminocity_hard, "width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    // Generate histograms from bin tuples
    TH1F *charmPtTotal = new TH1F("charm_full","Charm Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtPart = new TH1F("charm_pt_part","", 35, 0.0, 70.0);
    TH1F *charmPtHardest = new TH1F("charm_pt_hardest","", 35, 0.0, 70.0);
    TH1F *charmPtMPI = new TH1F("charm_pt_mpi","", 35, 0.0, 70.0);
    TH1F *charmPtISS = new TH1F("charm_pt_iss","", 35, 0.0, 70.0);
    TH1F *charmPtFSS = new TH1F("charm_pt_fss","", 35, 0.0, 70.0);
    TH1F *charmPtRemnant = new TH1F("charm_pt_remnant","", 35, 0.0, 70.0);
    TH1F *charmPtHadron = new TH1F("charm_pt_hadron","", 35, 0.0, 70.0);

    TH1F *bottomPtTotal = new TH1F("bottom_full","Bottom Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtPart = new TH1F("bottom_pt_part","", 35, 0.0, 70.0);
    TH1F *bottomPtHardest = new TH1F("bottom_pt_hardest","", 35, 0.0, 70.0);
    TH1F *bottomPtMPI = new TH1F("bottom_pt_mpi","", 35, 0.0, 70.0);
    TH1F *bottomPtISS = new TH1F("bottom_pt_iss","", 35, 0.0, 70.0);
    TH1F *bottomPtFSS = new TH1F("bottom_pt_fss","", 35, 0.0, 70.0);
    TH1F *bottomPtRemnant = new TH1F("bottom_pt_remnant","", 35, 0.0, 70.0);
    TH1F *bottomPtHadron = new TH1F("bottom_pt_hadron","", 35, 0.0, 70.0);

    for (int i = 0; i < nBins; ++i) {
        //charm
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotal->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>20 && status<30");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtHardest->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>30 && status<40");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtMPI->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>40 && status<50");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>50 && status<60");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtFSS->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>60 && status<70");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtRemnant->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>70 && status<80");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtHadron->Add(charmPtPart);

        // bottom
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotal->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>20 && status<30");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtHardest->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>30 && status<40");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtMPI->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>40 && status<50");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>50 && status<60");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtFSS->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>60 && status<70");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtRemnant->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>70 && status<80");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtHadron->Add(bottomPtPart);
    }

    //Plotting
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();

    TCanvas *canvasCharm = new TCanvas("charm_sigma","HF_sigma");

    charmPtTotal->SetLineColor(1);
    charmPtTotal->SetStats(0);
    charmPtTotal->Draw();

    charmPtHardest->SetLineColor(2);
    charmPtHardest->SetStats(0);
    charmPtHardest->Draw("SAME");

    charmPtMPI->SetLineColor(3);
    charmPtMPI->SetStats(0);
    charmPtMPI->Draw("SAME");

    charmPtISS->SetLineColor(4);
    charmPtISS->SetStats(0);
    charmPtISS->Draw("SAME");

    charmPtFSS->SetLineColor(5);
    charmPtFSS->SetStats(0);
    charmPtFSS->Draw("SAME");

    charmPtRemnant->SetLineColor(6);
    charmPtRemnant->SetStats(0);
    charmPtRemnant->Draw("SAME");

    charmPtHadron->SetLineColor(7);
    charmPtHadron->SetStats(0);
    charmPtHadron->Draw("SAME");

    auto legend2 = new TLegend();
    legend2->AddEntry(charmPtTotal,"Charm Total","l");
    legend2->AddEntry(charmPtHardest,"Hardest Subprocess","l");
    legend2->AddEntry(charmPtMPI,"Subsequent Subprocesses","l");
    legend2->AddEntry(charmPtISS,"Initial State Showers","l");
    legend2->AddEntry(charmPtFSS,"Final State Showers","l");
    legend2->AddEntry(charmPtRemnant,"Beam Remnants","l");
    legend2->AddEntry(charmPtHadron,"Hadronization Partons","l");
    legend2->Draw("SAME");

    canvasCharm->Write();

    TCanvas *canvasBottom = new TCanvas("bottom_sigma","HF_sigma");

    bottomPtTotal->SetLineColor(1);
    bottomPtTotal->SetStats(0);
    bottomPtTotal->Draw();

    bottomPtHardest->SetLineColor(2);
    bottomPtHardest->SetStats(0);
    bottomPtHardest->Draw("SAME");

    bottomPtMPI->SetLineColor(3);
    bottomPtMPI->SetStats(0);
    bottomPtMPI->Draw("SAME");

    bottomPtISS->SetLineColor(4);
    bottomPtISS->SetStats(0);
    bottomPtISS->Draw("SAME");

    bottomPtFSS->SetLineColor(5);
    bottomPtFSS->SetStats(0);
    bottomPtFSS->Draw("SAME");

    bottomPtRemnant->SetLineColor(6);
    bottomPtRemnant->SetStats(0);
    bottomPtRemnant->Draw("SAME");

    bottomPtHadron->SetLineColor(7);
    bottomPtHadron->SetStats(0);
    bottomPtHadron->Draw("SAME");

    auto legend3 = new TLegend();
    legend3->AddEntry(bottomPtTotal,"Bottom Total","l");
    legend3->AddEntry(bottomPtHardest,"Hardest Subprocess","l");
    legend3->AddEntry(bottomPtMPI,"Subsequent Subprocesses","l");
    legend3->AddEntry(bottomPtISS,"Initial State Showers","l");
    legend3->AddEntry(bottomPtFSS,"Final State Showers","l");
    legend3->AddEntry(bottomPtRemnant,"Beam Remnants","l");
    legend3->AddEntry(bottomPtHadron,"Hadronization Partons","l");
    legend3->Draw("SAME");

    canvasBottom->Write();
    

    delete outFile;

    return 0;
}