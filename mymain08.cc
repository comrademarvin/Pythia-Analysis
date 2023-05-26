#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

// Prints condensed version of particle mother-tree
void decayStatus(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int momIndex: eventObj[particleIndex].motherList()) {
            cout << " " << eventObj[momIndex].name() << " (" << momIndex << ", " << eventObj[momIndex].id() << ", " << eventObj[momIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int momIndex: eventObj[particleIndex].motherList()) {
            if ((!eventObj[momIndex].isGluon()) && (abs(eventObj[momIndex].id())!=2212)) {
                decayStatus(momIndex, eventObj, decayLevel);
            }
        }
    }
}

// Prints condensed version of particle daughter-tree
void decayStatusDaughters(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            cout << " " << eventObj[daughterIndex].name() << " (" << daughterIndex << ", " << eventObj[daughterIndex].id() << ", " << eventObj[daughterIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            if (!eventObj[daughterIndex].isFinal() && !eventObj[daughterIndex].isGluon()) {
                decayStatusDaughters(daughterIndex, eventObj, decayLevel);
            }
        }
    }
}

// returns index of decay muon in event, returns -1 if no decay muon is found
int muonDecay(int index, Pythia8::Event eventObj) {
    if (abs(eventObj[index].id()) == 13) return index;

    for (int daughterIndex: eventObj[index].daughterList()) {
        int returnIndex = muonDecay(daughterIndex, eventObj);
        if (returnIndex != -1) return returnIndex;
    }

    return -1;
}

int main() {
    // Turn SoftQCD on/off
    bool softQCD = true;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain08.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    if (softQCD) {
        nBins = 6;
        static const double tempArray[7] = {0.0, 8.0, 16.0, 26.0, 36.0, 50.0, 70.0};
        binEdges = &tempArray[0];
    } else {
        nBins = 4;
        static const double tempArray[5] = {16.0, 26.0, 36.0, 50.0, 70.0};
        binEdges = &tempArray[0];
    }

    double binLuminocity[nBins]; // luminocity from generated process sigma to calculate cross-sections

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 70, 0.0, 70.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 70, 0.0, 70.0);

    // HF Cross Sections
    vector<TNtuple*> charmTuples(nBins);
    vector<TNtuple*> bottomTuples(nBins);
    vector<TNtuple*> muonTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        charmTuples[i] = new TNtuple("charm", "charm", "pt:status");
        bottomTuples[i] = new TNtuple("bottom", "bottom", "pt:status");
        muonTuples[i] = new TNtuple("muon", "muon", "pt:decay:hfParent");
    }

    // Number of events to generate per bin.
    int N_events = 1000;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 2) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            // set pythia initialization variables
            pythia.readString("HardQCD:all = on");
            // pythia.readString("HardQCD:hardccbar = on");
            // pythia.readString("HardQCD:hardbbbar = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 5020.");
        pythia.readString("Tune:pp = 14");
        pythia.readString("411:onMode=off");
        pythia.readString("411:onIfAny=13");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat  = pythia.info.pTHat();

            if (softQCD && iBin < 2 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);

            for (int i = 0; i < pythia.event.size(); ++i) {
                int particleID = abs(pythia.event[i].id());
                int particleStatus = abs(pythia.event[i].status());
                int particlePt = pythia.event[i].pT();

                if (particleID == 4) { // charm
                    if (particleStatus == 23 || particleStatus == 33 || particleStatus == 43) { // || particleStatus == 51) {
                        charmTuples[iBin]->Fill(particlePt, particleStatus);
                        decayStatusDaughters(i, pythia.event);
                        cout << endl;
                    }
                }

                if (particleID == 5) { // bottom
                    if (particleStatus == 23 || particleStatus == 33 || particleStatus == 43) { // || particleStatus == 51) {
                        bottomTuples[iBin]->Fill(particlePt, particleStatus);
                        decayStatusDaughters(i, pythia.event);
                        cout << endl;
                    }
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
    TH1F *charmPtTotal = new TH1F("charm_full","Produced Charm Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *charmPtPart = new TH1F("charm_pt_part","", 35, 0.0, 70.0);
    TH1F *charmPtTotalMPI = new TH1F("charm_full_MPI","", 35, 0.0, 70.0);
    TH1F *charmPtTotalMPIShowers = new TH1F("charm_full_MPI_showers","", 35, 0.0, 70.0);

    // bottom 
    TH1F *bottomPtTotal = new TH1F("bottom_full","Produced Bottom Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtPart = new TH1F("bottom_pt_part","", 35, 0.0, 70.0);
    TH1F *bottomPtTotalMPI = new TH1F("bottom_full_MPI","", 35, 0.0, 70.0);
    TH1F *bottomPtTotalMPIShowers = new TH1F("bottom_full_MPI_showers","", 35, 0.0, 70.0);

    for (int i = 0; i < nBins; ++i) {
        //charm
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status<30");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotal->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part", "status<40");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotalMPI->Add(charmPtPart);

        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotalMPIShowers->Add(charmPtPart);

        // bottom
        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status<30");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotal->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part", "status<40");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotalMPI->Add(bottomPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotalMPIShowers->Add(bottomPtPart);
    }

    //Plotting
    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();

    // Produced HF Cross Sections
    TCanvas *canvasHF = new TCanvas("HF_sigma","HF_sigma");

    charmPtTotal->SetLineColor(1);
    charmPtTotal->SetStats(0);
    charmPtTotal->Draw();

    bottomPtTotal->SetLineColor(2);
    bottomPtTotal->SetStats(0);
    bottomPtTotal->Draw("SAME");

    auto legendHF = new TLegend();
    legendHF->AddEntry(charmPtTotal,"Charm","l");
    legendHF->AddEntry(bottomPtTotal,"Bottom","l");
    legendHF->Draw("SAME");

    canvasHF->Write();

    // Charm
    TCanvas *canvasCharm = new TCanvas("Charm_sigma","Charm_sigma");

    charmPtTotal->SetLineColor(1);
    charmPtTotal->SetStats(0);
    charmPtTotal->Draw();

    charmPtTotalMPI->SetLineColor(2);
    charmPtTotalMPI->SetStats(0);
    charmPtTotalMPI->Draw("SAME");

    charmPtTotalMPIShowers->SetLineColor(3);
    charmPtTotalMPIShowers->SetStats(0);
    charmPtTotalMPIShowers->Draw("SAME");

    auto legendCharm = new TLegend();
    legendCharm->AddEntry(charmPtTotal,"Hardest","l");
    legendCharm->AddEntry(charmPtTotalMPI,"Hardest+MPI","l");
    legendCharm->AddEntry(charmPtTotalMPIShowers,"Hardest+MPI+Showers","l");
    legendCharm->Draw("SAME");

    canvasCharm->Write();

    // bottom
    TCanvas *canvasBottom = new TCanvas("Bottom_sigma","Bottom_sigma");

    bottomPtTotal->SetLineColor(1);
    bottomPtTotal->SetStats(0);
    bottomPtTotal->Draw();

    bottomPtTotalMPI->SetLineColor(2);
    bottomPtTotalMPI->SetStats(0);
    bottomPtTotalMPI->Draw("SAME");

    bottomPtTotalMPIShowers->SetLineColor(3);
    bottomPtTotalMPIShowers->SetStats(0);
    bottomPtTotalMPIShowers->Draw("SAME");

    auto legendBottom = new TLegend();
    legendBottom->AddEntry(bottomPtTotal,"Hardest","l");
    legendBottom->AddEntry(bottomPtTotalMPI,"Hardest+MPI","l");
    legendBottom->AddEntry(bottomPtTotalMPIShowers,"Hardest+MPI+Showers","l");
    legendBottom->Draw("SAME");

    canvasBottom->Write();

    delete outFile;

    return 0;
}