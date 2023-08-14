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
    TFile* outFile = new TFile("mymain10.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    if (softQCD) {
        nBins = 10;
        static const double tempArray[11] = {0.0, 10.0, 15.0, 20.0, 26.0, 32.0, 40.0, 50.0, 62.0, 76.0, 100.0};
        binEdges = &tempArray[0];
    } else {
        nBins = 5;
        static const double tempArray[6] = {16.0, 26.0, 36.0, 50.0, 70.0, 100.0};
        binEdges = &tempArray[0];
    }

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 100, 0.0, 100.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 100, 0.0, 100.0);

    // HF Cross Sections
    vector<TNtuple*> HFEventTuples(nBins);
    vector<TNtuple*> HFProducedTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        HFEventTuples[i] = new TNtuple(Form("HF_events%d", i), "HF_events", "ptHat:prodCount"); // see definition of prodMech below
        HFProducedTuples[i] = new TNtuple(Form("HF_produced%d", i), "HF_produced", "ID:pt:y:prodCount:tau"); // see definition of prodMech below
    }

    // Number of events to generate per bin.
    int N_events = 100000;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin == 0) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            pythia.readString("HardQCD:all = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 5020.");
        pythia.readString("Tune:pp = 14");
        // pythia.readString("PartonLevel:all = on");
        // pythia.readString("HadronLevel:all = off");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        int eventCount = 0;

        int events_run;
        if (softQCD && iBin == 0) {
            events_run = 5*N_events;
        } else {
            events_run = N_events;
        }

        for (int iEvent = 0; iEvent < events_run; ++iEvent) {
            if (!pythia.next()) continue;

            //int subprocessCode = pythia.info.codeSub();

            double pTHat  = pythia.info.pTHat();

            if (softQCD && iBin == 0 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            //if (pTHat < binEdges[iBin]) continue;

            eventCount++;

            hardPtPart->Fill(pTHat);

            //cout << "====START OF NEW EVENT====" << endl;
            int HFCount = 0;
            int tempHFArray[20];

            int hardestHFCount = 0;
            bool HFEvent = false;

            for (int i = 0; i < pythia.event.size(); ++i) {
                int particleID = abs(pythia.event[i].id());
                int particleStatus = abs(pythia.event[i].status());

                if (particleID == 4 || particleID == 5) { // is HF
                    HFEvent = true;
                    if (particleStatus == 23) hardestHFCount++;
                    
                    if (particleStatus > 70 && particleStatus < 80) {
                        tempHFArray[HFCount] = i;
                        HFCount++;
                    } 
                } 
            }

            if (HFEvent) {
                HFEventTuples[iBin]->Fill(pTHat, hardestHFCount);
            }
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
        HFEventTuples[i]->Write(Form("HF_events%d", i));
    }

    outFile->WriteObject(&binLuminocity, "luminocity");

    delete outFile;

    return 0;
}