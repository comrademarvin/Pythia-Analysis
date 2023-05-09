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
            if (!eventObj[daughterIndex].isFinal()) {
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
    int N_events = 10000;

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

    TH1F *bottomPtTotal = new TH1F("bottom_full","Bottom Cross-Sections;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 35, 0.0, 70.0);
    TH1F *bottomPtPart = new TH1F("bottom_pt_part","", 35, 0.0, 70.0);

    for (int i = 0; i < nBins; ++i) {
        charmPtPart->Reset();
        charmTuples[i]->Draw("pt>>charm_pt_part");
        charmPtPart->Scale(1/binLuminocity[i], "width");
        charmPtTotal->Add(charmPtPart);

        bottomPtPart->Reset();
        bottomTuples[i]->Draw("pt>>bottom_pt_part");
        bottomPtPart->Scale(1/binLuminocity[i], "width");
        bottomPtTotal->Add(bottomPtPart);
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

    auto legend2 = new TLegend();
    legend2->AddEntry(charmPtTotal,"Charm Total","l");
    legend2->Draw("SAME");

    canvasCharm->Write();

    TCanvas *canvasBottom = new TCanvas("bottom_sigma","HF_sigma");

    bottomPtTotal->SetLineColor(1);
    bottomPtTotal->SetStats(0);
    bottomPtTotal->Draw();

    auto legend3 = new TLegend();
    legend3->AddEntry(bottomPtTotal,"Bottom Total","l");
    legend3->Draw("SAME");

    canvasBottom->Write();
    

    delete outFile;

    return 0;
}