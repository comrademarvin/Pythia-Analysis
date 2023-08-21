#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>
#include <string>

using namespace Pythia8;

int main() {
    // Turn SoftQCD on/off
    bool softQCD = true;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain09.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    if (softQCD) {
        nBins = 7;
        static const double tempArray[8] = {0.0, 10.0, 30.0, 50.0, 75.0, 100.0, 150.0, 200.0};
        binEdges = &tempArray[0];
    } else {
        nBins = 4;
        static const double tempArray[5] = {10.0, 30.0, 50.0, 75.0, 100.0};
        binEdges = &tempArray[0];
    }

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 200, 0.0, 200.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 200, 0.0, 200.0);

    // HF Cross Sections
    vector<TNtuple*> muonTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        muonTuples[i] = new TNtuple(Form("muon%d", i), "muon", "binTag:eventTag:pAbs:pt:y:eta:decayStatus:firstMother:lastMother"); // see definition of decayStaus below
    }

    // Number of events to generate per bin.
    int N_events = 10000;

    // Store multiplicity of each event
    // int multBins;
    // int scaleSoftBin = 3;
    // if (softQCD) { 
    //     multBins = nBins+(scaleSoftBin-1);
    // } else {
    //     multBins = nBins;
    // }
    // vector<int> multiplicities(multBins*N_events, -1);

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin == 0) {
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

        // // D+ forced muon decay
        // pythia.readString("411:onMode=off");
        // pythia.readString("411:onIfAny=13");
        // // D0 forced muon decay
        // pythia.readString("421:onMode=off");
        // pythia.readString("421:onIfAny=13");

        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        int eventCount = 0;

        int events_run;
        if (softQCD && iBin == 0) {
            events_run = 10*N_events;
        } else if (iBin == 1) {
            events_run = 5*N_events;
        } else {
            events_run = N_events;
        }

        for (int iEvent = 0; iEvent < events_run; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat  = pythia.info.pTHat();

            if (softQCD && iBin == 0 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            //if (pTHat < binEdges[iBin]) continue;

            eventCount++;

            hardPtPart->Fill(pTHat);

            //cout << "====START OF NEW EVENT====" << endl;

            int multCount = 0;

            for (int i = 0; i < pythia.event.size(); ++i) {
                if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) multCount++;

                int particleID = abs(pythia.event[i].id());
                int particleStatus = abs(pythia.event[i].status());
                double particlePAbs = pythia.event[i].pAbs();
                double particlePt = pythia.event[i].pT();
                double particleRapidity = pythia.event[i].y();
                double particlePseudorapidity = pythia.event[i].eta();

                if (particleID == 13) { // muon
                    int motherIndex = pythia.event[i].mother1();
                    int firstMotherID = abs(pythia.event[motherIndex].id());

                    string decayOutput = "mu";
                    bool motherHadronStatus = pythia.event[motherIndex].isHadron();
                    int prevIndex = -1;
                    while (motherHadronStatus) {
                        string hadronOutput = " -> " + pythia.event[motherIndex].name();
                        decayOutput = decayOutput + hadronOutput;
                        prevIndex = motherIndex;
                        motherIndex = pythia.event[motherIndex].mother1();
                        motherHadronStatus = pythia.event[motherIndex].isHadron();
                    }

                    int lastMotherID = firstMotherID;
                    if (prevIndex != -1) {
                        lastMotherID = abs(pythia.event[prevIndex].id());
                        decayOutput = decayOutput + " -> [";
                        for (int index: pythia.event[prevIndex].motherList()) {
                            if (!pythia.event[index].isGluon()) {
                                decayOutput = decayOutput + pythia.event[index].name();
                                decayOutput = decayOutput + " ";
                            }
                        }
                        decayOutput = decayOutput + "]";
                    } else {
                        string otherOutput = " -> " + pythia.event[motherIndex].name();
                        decayOutput = decayOutput + otherOutput;
                    }

                    // decayStatus definition:
                    // 0: c->D->mu (charm meson)
                    // 1: b->B->mu (bottom meson)
                    // 2: b->B->D->mu (bottom to charm meson)
                    // 3: ?->Baryon_c->mu
                    // 4: ?->Baryon_b->mu
                    // 5: ?->Quarkonium_c->mu
                    // 6: ?->Quarkonium_b->mu
                    // 7: tau->mu

                    int decayStatusTemp = -1;
                    if (firstMotherID >= 411 && firstMotherID <= 435) { // D Meson decay
                        if (lastMotherID >= 411 && lastMotherID <= 435) { // case 0
                            decayStatusTemp = 0;
                        } else if (lastMotherID >= 511 && lastMotherID <= 545) { // case 2
                            decayStatusTemp = 2;
                        }
                    } else if (firstMotherID >= 511 && firstMotherID <= 545) { // B Meson decay
                        if (lastMotherID >= 511 && lastMotherID <= 545) { // case 1
                            decayStatusTemp = 1;
                        }
                    } else if (firstMotherID >= 4122 && firstMotherID <= 4444) { // case 3
                        decayStatusTemp = 3;
                    } else if (firstMotherID >= 5122 && firstMotherID <= 5554) { // case 4
                        decayStatusTemp = 4;
                    } else if (firstMotherID >= 411 && firstMotherID <= 445) { // case 5
                        decayStatusTemp = 5;
                    } else if (firstMotherID >= 551 && firstMotherID <= 557) { // case 5
                        decayStatusTemp = 6;
                    } else if (firstMotherID == 15) {
                        decayStatusTemp = 7;
                    } else {
                        cout << decayOutput << endl;
                    }

                    muonTuples[iBin]->Fill(iBin, iEvent, particlePAbs, particlePt, particleRapidity, particlePseudorapidity, decayStatusTemp, firstMotherID, lastMotherID);
                }
            }

            // if (softQCD) {
            //     if (iBin == 0) {
            //         multiplicities[iEvent] = multCount;
            //     } else {
            //         multiplicities[((iBin+(scaleSoftBin-1))*N_events) + iEvent] = multCount;
            //     }
            // } else {
            //     multiplicities[(iBin*N_events) + iEvent] = multCount;
            // }
        }

        // cross-section for the bin
        double luminocity_hard = events_run/(pythia.info.sigmaGen()*pow(10,9));
        binLuminocity[iBin] = luminocity_hard;

        hardPtPart->Scale(1/luminocity_hard, "width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    hardPt->Write();

    // Read Tuples to output root file
    for (int i = 0; i < nBins; ++i) {
        muonTuples[i]->Write(Form("muon%d", i));
    }

    outFile->WriteObject(&binLuminocity, "luminocity");
    //outFile->WriteObject(&multiplicities, "multiplicity");

    delete outFile;

    return 0;
}