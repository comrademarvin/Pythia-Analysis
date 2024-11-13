#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>
#include <string>

using namespace Pythia8;

int main(int argc, char *argv[]) {
    // Bin index from command line argument
    int iBin = atoi(argv[1]);

    // Turn SoftQCD on/off
    bool softQCD = true;

    // Pythia object
    PythiaParallel pythia; // generate events in parallel

    // ROOT file for histograms
    TFile* outFile = new TFile(Form("mymain09_HPC_root_20M_502_tune2/mymain09_%d.root", iBin), "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    if (softQCD) {
        nBins = 7;
        static const double tempArray[8] = {0.0, 10.0, 40.0, 70.0, 100.0, 150.0, 200.0, 300.0};
        binEdges = &tempArray[0];
    } else {
        nBins = 4;
        static const double tempArray[5] = {10.0, 30.0, 50.0, 75.0, 100.0};
        binEdges = &tempArray[0];
    }

    // generation info for normalisation
    vector<double> genInfo(nBins, 0.0);

    // Histograms
    // Total Cross Section
    int N_bins_hard_pt = 200;
    TH1F *hardPt = new TH1F("pT_hat","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_hard_pt, 0.0, 200.0);

    // HF Cross Sections
    TNtuple* muonTuple = new TNtuple("muon", "muon", "binTag:eventTag:pAbs:pt:y:eta:firstMother:lastMother:decayStatus");

    // Number of events to generate per bin.
    int N_events = 20000000;

    int events_run;
    if (softQCD && iBin == 0) {
        events_run = 3*N_events;
    } else if (iBin == 1) {
        events_run = 2*N_events;
    } else {
        events_run = N_events;
    }

    // run events for each ptHat bin 
    if (softQCD && iBin == 0) {
        pythia.readString("HardQCD:all = off");
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        // set pythia initialization variables
        pythia.readString("HardQCD:all = on");
        pythia.readString("SoftQCD:nonDiffractive = off");
    }

    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("Tune:pp = 2"); // Tune 1
    pythia.readString("Parallelism:numThreads = 15"); // number of threads assigned per bin
    
    // D+ forced muon decay
    pythia.readString("411:onMode=off");
    pythia.readString("411:onIfAny=13");
    // D0 forced muon decay
    pythia.readString("421:onMode=off");
    pythia.readString("421:onIfAny=13");

    pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
    pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);

    if (!pythia.init()) return 1;	//initiate pythia and output an error if it doesn't initiate.

    int eventCount = 0;

    int iEvent = 0;
    pythia.run(events_run, [&](Pythia* pythiaPtr)
	{
		// giving reference to the instance that generated the event.
		Event& event = pythiaPtr->event; 
		const Info& info = pythiaPtr->info;

        double pTHat  = info.pTHat();

        if (softQCD && iBin == 0 && info.isNonDiffractive()
        && pTHat > binEdges[iBin+1]) return;

        eventCount++;

        hardPt->Fill(pTHat);

        //cout << "====START OF NEW EVENT====" << endl;

        for (int i = 0; i < event.size(); ++i) {
            auto* particle = &event[i];
            int particleID = abs(particle->id());

            if (particleID == 13 && particle->isFinal() && (particle->eta()>2.5) && (particle->eta()<4.0)) { // final state muons in forward rapidity region
                // find first non-copy mother
                auto* muonFinal = particle;
                string decayOutput = "==== Here forward muon decay: mu";
                while ((particle->mother2() != 0) && (particle->mother1() == particle->mother2())) {
                    particle = &event[particle->mother1()];
                    string partOutput = " <- " + particle->name();
                    decayOutput = decayOutput + partOutput;
                }
                auto* firstMother = &event[particle->mother1()];
                decayOutput = decayOutput + " <- ";
                decayOutput = decayOutput + firstMother->name();
                decayOutput = decayOutput + " (firstMother)";

                int firstMotherID = abs(firstMother->id());
                if ((firstMotherID >= 411 && firstMotherID <= 435) || (firstMotherID >= 511 && firstMotherID <= 545)) { // only consider HF->mu decays
                    // then find last mother formed by hadronisation
                    auto* lastMother = firstMother;
                    auto* iterator = &event[firstMother->mother1()];
                    while (iterator->isHadron()) {
                        decayOutput = decayOutput + " <- ";
                        decayOutput = decayOutput + iterator->name();
                        lastMother = iterator;
                        iterator = &event[iterator->mother1()];
                    }
                    decayOutput = decayOutput + " (lastMother)";
                    decayOutput = decayOutput + " <- ";
                    decayOutput = decayOutput + iterator->name();

                    // determine the HF decayStates, defined as:
                    // -1: ? -> D/B -> mu
                    // 0: c->D->mu (charm meson)
                    // 1: b->B->mu (bottom meson)
                    // 2: b->B->D->mu (bottom to charm meson)

                    int decayStatus = -1;
                    int lastMotherID = abs(lastMother->id());
                    if (firstMotherID >= 411 && firstMotherID <= 435) { // muon from D-meson decay
                        if (lastMotherID >= 411 && lastMotherID <= 435) decayStatus = 0;
                        else if (lastMotherID >= 511 && lastMotherID <= 545) decayStatus = 2;
                    } else if (lastMotherID >= 511 && lastMotherID <= 545) decayStatus = 1; // muon from D-meson decay

                    // std::cout << decayOutput << " = decayStatus: " << decayStatus << std::endl;

                    muonTuple->Fill(iBin, iEvent, muonFinal->pAbs(), muonFinal->pT(), muonFinal->y(), muonFinal->eta() , firstMotherID, lastMotherID, decayStatus);
                }
            }
        }

        iEvent += 1;
    });

    // cross-section for the bin
    double genInfoNorm = pythia.sigmaGen()/pythia.weightSum(); // in mb
    genInfo[iBin] = genInfoNorm;
    outFile->WriteObject(&genInfo, "genInfoNorm");

    hardPt->Scale(genInfoNorm, "width");
    hardPt->Write();
    
    muonTuple->Write("muon");

    delete outFile;

    return 0;
}