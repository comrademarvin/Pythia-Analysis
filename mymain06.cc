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
    int N_events = 100000;

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
    // Charm
    TH1F *charmPtTotal = new TH1F("charm_full","Charm Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);

    TH1F *charmPtPart = new TH1F("charm_pt_part","", 35, 0.0, 70.0);

    TH1F *charmPtHardest = new TH1F("charm_pt_hardest","Charm Hardest Subprocess Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtHardestInc = new TH1F("charm_pt_hardest_inc","", 35, 0.0, 70.0);
    TH1F *charmPtHardestOut = new TH1F("charm_pt_hardest_out","", 35, 0.0, 70.0);

    TH1F *charmPtMPI = new TH1F("charm_pt_mpi","Charm Subsequent Subprocesses (MPI) Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtMPIInc = new TH1F("charm_pt_mpi_inc","", 35, 0.0, 70.0);
    TH1F *charmPtMPIOut = new TH1F("charm_pt_mpi_out","", 35, 0.0, 70.0);

    TH1F *charmPtISS = new TH1F("charm_pt_iss","Charm Initial-State-Showers Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtISS41 = new TH1F("charm_pt_iss_41","", 35, 0.0, 70.0);
    TH1F *charmPtISS42 = new TH1F("charm_pt_iss_42","", 35, 0.0, 70.0);
    TH1F *charmPtISS43 = new TH1F("charm_pt_iss_43","", 35, 0.0, 70.0);
    TH1F *charmPtISS44 = new TH1F("charm_pt_iss_44","", 35, 0.0, 70.0);

    TH1F *charmPtFSS = new TH1F("charm_pt_fss","Charm Final-State-Showers Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtFSS51 = new TH1F("charm_pt_fss_41","", 35, 0.0, 70.0);
    TH1F *charmPtFSS52 = new TH1F("charm_pt_fss_42","", 35, 0.0, 70.0);
    TH1F *charmPtFSS53 = new TH1F("charm_pt_fss_43","", 35, 0.0, 70.0);

    TH1F *charmPtRemnant = new TH1F("charm_pt_remnant","", 35, 0.0, 70.0);

    // TH1F *charmPtHadron = new TH1F("charm_pt_hadron","", 35, 0.0, 70.0);

    // bottom 
    TH1F *bottomPtTotal = new TH1F("bottom_full","Bottom Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);

    TH1F *bottomPtPart = new TH1F("bottom_pt_part","", 35, 0.0, 70.0);

    TH1F *bottomPtHardest = new TH1F("bottom_pt_hardest","Bottom Hardest Subprocess Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtHardestInc = new TH1F("bottom_pt_hardest_inc","", 35, 0.0, 70.0);
    TH1F *bottomPtHardestOut = new TH1F("bottom_pt_hardest_out","", 35, 0.0, 70.0);

    TH1F *bottomPtMPI = new TH1F("bottom_pt_mpi","Bottom Subsequent Subprocesses (MPI) Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtMPIInc = new TH1F("bottom_pt_mpi_inc","", 35, 0.0, 70.0);
    TH1F *bottomPtMPIOut = new TH1F("bottom_pt_mpi_out","", 35, 0.0, 70.0);

    TH1F *bottomPtISS = new TH1F("bottom_pt_iss","Bottom Initial-State-Showers Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtISS41 = new TH1F("bottom_pt_iss_41","", 35, 0.0, 70.0);
    TH1F *bottomPtISS42 = new TH1F("bottom_pt_iss_42","", 35, 0.0, 70.0);
    TH1F *bottomPtISS43 = new TH1F("bottom_pt_iss_43","", 35, 0.0, 70.0);
    TH1F *bottomPtISS44 = new TH1F("bottom_pt_iss_44","", 35, 0.0, 70.0);

    TH1F *bottomPtFSS = new TH1F("bottom_pt_fss","Bottom Final-State-Showers Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtFSS51 = new TH1F("bottom_pt_fss_41","", 35, 0.0, 70.0);
    TH1F *bottomPtFSS52 = new TH1F("bottom_pt_fss_42","", 35, 0.0, 70.0);
    TH1F *bottomPtFSS53 = new TH1F("bottom_pt_fss_43","", 35, 0.0, 70.0);

    TH1F *bottomPtRemnant = new TH1F("bottom_pt_remnant","", 35, 0.0, 70.0);

    //TH1F *bottomPtHadron = new TH1F("bottom_pt_hadron","", 35, 0.0, 70.0);

    for (int i = 0; i < nBins; ++i) {
        //charm
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotal->Add(charmPtPart);

        // hardest
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>20 && status<25");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtHardest->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==21");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtHardestInc->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==23");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtHardestOut->Add(charmPtPart);

        // MPI
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>30 && status<35");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtMPI->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==31");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtMPIInc->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==33");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtMPIOut->Add(charmPtPart);

        // ISS
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>40 && status<50");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==41");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS41->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==42");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS42->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==43");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS43->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==44");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtISS44->Add(charmPtPart);

        // FSS
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>50 && status<60");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtFSS->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==51");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtFSS51->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==52");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtFSS52->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status==53 || status==54 || status==55");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtFSS53->Add(charmPtPart);

        // Remnant
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status>60 && status<70");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtRemnant->Add(charmPtPart);

        // charmPtPart->Reset();
        // charmTuples[i]->Draw("pt>>charm_pt_part", "status>70 && status<80");
        // charmPtPart->Scale(1/binLuminocity[i], "width");
        // charmPtHadron->Add(charmPtPart);

        // bottom
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotal->Add(bottomPtPart);

        // Hardest
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>20 && status<30");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtHardest->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==21");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtHardestInc->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==23");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtHardestOut->Add(bottomPtPart);

        // MPI
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>30 && status<35");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtMPI->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==31");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtMPIInc->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==33");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtMPIOut->Add(bottomPtPart);

        // ISS
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>40 && status<50");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==41");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS41->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==42");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS42->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==43");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS43->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==44");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtISS44->Add(bottomPtPart);

        // FSS
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>50 && status<60");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtFSS->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==51");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtFSS51->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==52");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtFSS52->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status==53 || status==54 || status==55");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtFSS53->Add(bottomPtPart);

        // Remnant
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>60 && status<70");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtRemnant->Add(bottomPtPart);

        // bottomPtPart->Reset();
        // bottomTuples[i]->Draw("pt>>bottom_pt_part", "status>70 && status<80");
        // bottomPtPart->Scale(1/binLuminocity[i], "width");
        // bottomPtHadron->Add(bottomPtPart);
    }

    //Plotting
    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();

    // Charm Full
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

    // charmPtHadron->SetLineColor(7);
    // charmPtHadron->SetStats(0);
    // charmPtHadron->Draw("SAME");

    auto legend2 = new TLegend();
    legend2->AddEntry(charmPtTotal,"Charm Total","l");
    legend2->AddEntry(charmPtHardest,"Hardest Subprocess","l");
    legend2->AddEntry(charmPtMPI,"Subsequent Subprocesses","l");
    legend2->AddEntry(charmPtISS,"Initial State Showers","l");
    legend2->AddEntry(charmPtFSS,"Final State Showers","l");
    legend2->AddEntry(charmPtRemnant,"Beam Remnants","l");
    // legend2->AddEntry(charmPtHadron,"Hadronization Partons","l");
    legend2->Draw("SAME");

    canvasCharm->Write();

    // Charm Hardest
    TCanvas *canvasCharmHardest = new TCanvas("charm_hardest_sigma","HF_sigma");

    charmPtHardest->SetLineColor(1);
    charmPtHardest->SetStats(0);
    charmPtHardest->Draw();

    charmPtHardestInc->SetLineColor(2);
    //charmPtHardestInc->SetStats(0);
    charmPtHardestInc->Draw("SAME");

    charmPtHardestOut->SetLineColor(4);
    //charmPtHardestOut->SetStats(0);
    charmPtHardestOut->Draw("SAME");

    auto legendCharmHardest = new TLegend();
    legendCharmHardest->AddEntry(charmPtHardest,"Hardest Subprocess","l");
    legendCharmHardest->AddEntry(charmPtHardestInc,"Incoming","l");
    legendCharmHardest->AddEntry(charmPtHardestOut,"Outgoing","l");
    legendCharmHardest->Draw("SAME");

    canvasCharmHardest->Write();

    // Charm MPI
    TCanvas *canvasCharmMPI = new TCanvas("charm_MPI_sigma","HF_sigma");

    charmPtMPI->SetLineColor(1);
    charmPtMPI->SetStats(0);
    charmPtMPI->Draw();

    charmPtMPIInc->SetLineColor(2);
    charmPtMPIInc->SetStats(0);
    charmPtMPIInc->Draw("SAME");

    charmPtMPIOut->SetLineColor(4);
    charmPtMPIOut->SetStats(0);
    charmPtMPIOut->Draw("SAME");

    auto legendCharmMPI = new TLegend();
    legendCharmMPI->AddEntry(charmPtMPI,"Subsequent Subprocesses","l");
    legendCharmMPI->AddEntry(charmPtMPIInc,"Incoming","l");
    legendCharmMPI->AddEntry(charmPtMPIOut,"Outgoing","l");
    legendCharmMPI->Draw("SAME");

    canvasCharmMPI->Write();

    // ISS
    TCanvas *canvasCharmISS = new TCanvas("charm_ISS_sigma","HF_sigma");

    charmPtISS->SetLineColor(1);
    charmPtISS->SetStats(0);
    charmPtISS->Draw();

    charmPtISS41->SetLineColor(2);
    charmPtISS41->SetStats(0);
    charmPtISS41->Draw("SAME");

    charmPtISS42->SetLineColor(3);
    charmPtISS42->SetStats(0);
    charmPtISS42->Draw("SAME");

    charmPtISS43->SetLineColor(4);
    charmPtISS43->SetStats(0);
    charmPtISS43->Draw("SAME");

    charmPtISS44->SetLineColor(5);
    charmPtISS44->SetStats(0);
    charmPtISS44->Draw("SAME");

    auto legendCharmISS = new TLegend();
    legendCharmISS->AddEntry(charmPtISS,"Initial State Showers","l");
    legendCharmISS->AddEntry(charmPtISS41,"Incoming Spacelike Main Branch","l");
    legendCharmISS->AddEntry(charmPtISS42,"Incoming Copy of Recoiler","l");
    legendCharmISS->AddEntry(charmPtISS43,"Outgoing Produced by Branching","l");
    legendCharmISS->AddEntry(charmPtISS44,"Outgoing Shifted by Branching","l");
    legendCharmISS->Draw("SAME");

    canvasCharmISS->Write();

    // FSS
    TCanvas *canvasCharmFSS = new TCanvas("charm_FSS_sigma","HF_sigma");

    charmPtFSS->SetLineColor(1);
    charmPtFSS->SetStats(0);
    charmPtFSS->Draw();

    charmPtFSS51->SetLineColor(2);
    charmPtFSS51->SetStats(0);
    charmPtFSS51->Draw("SAME");

    charmPtFSS52->SetLineColor(3);
    charmPtFSS52->SetStats(0);
    charmPtFSS52->Draw("SAME");

    charmPtFSS53->SetLineColor(4);
    charmPtFSS53->SetStats(0);
    charmPtFSS53->Draw("SAME");

    auto legendCharmFSS = new TLegend();
    legendCharmFSS->AddEntry(charmPtFSS,"Final State Showers","l");
    legendCharmFSS->AddEntry(charmPtFSS51,"Outgoing Produced by Branching","l");
    legendCharmFSS->AddEntry(charmPtFSS52,"Outgoing Copy of Recoiler","l");
    legendCharmFSS->AddEntry(charmPtFSS53,"Copy of Recoiler","l");
    legendCharmFSS->Draw("SAME");

    canvasCharmFSS->Write();

    // Bottom Full
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

    // bottomPtHadron->SetLineColor(7);
    // bottomPtHadron->SetStats(0);
    // bottomPtHadron->Draw("SAME");

    auto legend3 = new TLegend();
    legend3->AddEntry(bottomPtTotal,"Bottom Total","l");
    legend3->AddEntry(bottomPtHardest,"Hardest Subprocess","l");
    legend3->AddEntry(bottomPtMPI,"Subsequent Subprocesses","l");
    legend3->AddEntry(bottomPtISS,"Initial State Showers","l");
    legend3->AddEntry(bottomPtFSS,"Final State Showers","l");
    legend3->AddEntry(bottomPtRemnant,"Beam Remnants","l");
    // legend3->AddEntry(bottomPtHadron,"Hadronization Partons","l");
    legend3->Draw("SAME");

    canvasBottom->Write();

    // Bottom Hardest
    TCanvas *canvasbottomHardest = new TCanvas("bottom_hardest_sigma","HF_sigma");

    bottomPtHardest->SetLineColor(1);
    bottomPtHardest->SetStats(0);
    bottomPtHardest->Draw();

    bottomPtHardestInc->SetLineColor(2);
    //bottomPtHardestInc->SetStats(0);
    bottomPtHardestInc->Draw("SAME");

    bottomPtHardestOut->SetLineColor(4);
    //bottomPtHardestOut->SetStats(0);
    bottomPtHardestOut->Draw("SAME");

    auto legendbottomHardest = new TLegend();
    legendbottomHardest->AddEntry(bottomPtHardest,"Hardest Subprocess","l");
    legendbottomHardest->AddEntry(bottomPtHardestInc,"Incoming","l");
    legendbottomHardest->AddEntry(bottomPtHardestOut,"Outgoing","l");
    legendbottomHardest->Draw("SAME");

    canvasbottomHardest->Write();
    
    // Bottom MPI
    TCanvas *canvasBottomMPI = new TCanvas("bottom_MPI_sigma","HF_sigma");

    bottomPtMPI->SetLineColor(1);
    bottomPtMPI->SetStats(0);
    bottomPtMPI->Draw();

    bottomPtMPIInc->SetLineColor(2);
    bottomPtMPIInc->SetStats(0);
    bottomPtMPIInc->Draw("SAME");

    bottomPtMPIOut->SetLineColor(4);
    bottomPtMPIOut->SetStats(0);
    bottomPtMPIOut->Draw("SAME");

    auto legendbottomMPI = new TLegend();
    legendbottomMPI->AddEntry(bottomPtMPI,"Subsequent Subprocesses","l");
    legendbottomMPI->AddEntry(bottomPtMPIInc,"Incoming","l");
    legendbottomMPI->AddEntry(bottomPtMPIOut,"Outgoing","l");
    legendbottomMPI->Draw("SAME");

    canvasBottomMPI->Write();

    // ISS
    TCanvas *canvasbottomISS = new TCanvas("bottom_ISS_sigma","HF_sigma");

    bottomPtISS->SetLineColor(1);
    bottomPtISS->SetStats(0);
    bottomPtISS->Draw();

    bottomPtISS41->SetLineColor(2);
    bottomPtISS41->SetStats(0);
    bottomPtISS41->Draw("SAME");

    bottomPtISS42->SetLineColor(3);
    bottomPtISS42->SetStats(0);
    bottomPtISS42->Draw("SAME");

    bottomPtISS43->SetLineColor(4);
    bottomPtISS43->SetStats(0);
    bottomPtISS43->Draw("SAME");

    bottomPtISS44->SetLineColor(5);
    bottomPtISS44->SetStats(0);
    bottomPtISS44->Draw("SAME");

    auto legendbottomISS = new TLegend();
    legendbottomISS->AddEntry(bottomPtISS,"Initial State Showers","l");
    legendbottomISS->AddEntry(bottomPtISS41,"Incoming Spacelike Main Branch","l");
    legendbottomISS->AddEntry(bottomPtISS42,"Incoming Copy of Recoiler","l");
    legendbottomISS->AddEntry(bottomPtISS43,"Outgoing Produced by Branching","l");
    legendbottomISS->AddEntry(bottomPtISS44,"Outgoing Shifted by Branching","l");
    legendbottomISS->Draw("SAME");

    canvasbottomISS->Write();

    // FSS
    TCanvas *canvasbottomFSS = new TCanvas("bottom_FSS_sigma","HF_sigma");

    bottomPtFSS->SetLineColor(1);
    bottomPtFSS->SetStats(0);
    bottomPtFSS->Draw();

    bottomPtFSS51->SetLineColor(2);
    bottomPtFSS51->SetStats(0);
    bottomPtFSS51->Draw("SAME");

    bottomPtFSS52->SetLineColor(3);
    bottomPtFSS52->SetStats(0);
    bottomPtFSS52->Draw("SAME");

    bottomPtFSS53->SetLineColor(4);
    bottomPtFSS53->SetStats(0);
    bottomPtFSS53->Draw("SAME");

    auto legendbottomFSS = new TLegend();
    legendbottomFSS->AddEntry(bottomPtFSS,"Final State Showers","l");
    legendbottomFSS->AddEntry(bottomPtFSS51,"Outgoing Produced by Branching","l");
    legendbottomFSS->AddEntry(bottomPtFSS52,"Outgoing Copy of Recoiler","l");
    legendbottomFSS->AddEntry(bottomPtFSS53,"Copy of Recoiler","l");
    legendbottomFSS->Draw("SAME");

    canvasbottomFSS->Write();

    delete outFile;

    return 0;
}