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

    // Pythia object
    PythiaParallel pythia; // generate events in parallel

    // ROOT file for histograms
    TFile* outFile = new TFile(Form("mymain05_HPC_root_50M_536/mymain05_%d.root", iBin), "RECREATE");

    // Turn SoftQCD on/off
    bool softQCD = true;

    // scale number of events per pT-hat bin
    int N_events = 50000000;
    if (softQCD && iBin == 0) {
        N_events = 3*N_events;
    } else if (iBin == 1) {
        N_events = 2*N_events;
    }

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

    // generation info for normalisation
    vector<double> genInfo(nBins, 0.0);

    // Histograms
    // Total Cross Section
    int N_bins_hard_pt = 200;
    TH1F *hardPt = new TH1F("pT_hat","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_hard_pt, 0.0, 200.0);

    // Generated Muons
    TNtuple* muonTuple = new TNtuple("muon", "muon", "binTag:eventTag:pAbs:pt:y:eta:firstMother");

    // configure and initialize pythia object
    pythia.readString("Beams:eCM = 5360.");
    pythia.readString("Tune:pp = 14"); // Monash Tune
    pythia.readString("Parallelism:numThreads = 15"); // number of threads assigned per bin

    // processes for each ptHat bin 
    if (softQCD && iBin == 0) {
        pythia.readString("HardQCD:all = off");
        pythia.readString("SoftQCD:nonDiffractive = on");
    } else {
        // set pythia initialization variables
        pythia.readString("HardQCD:all = on");
        pythia.readString("SoftQCD:nonDiffractive = off");
    }

    pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
    pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);

    if (!pythia.init()) return 1;	//initiate pythia and output an error if it doesn't initiate.

    // generate events in parallel
    int iEvent = 0;
    pythia.run(N_events, [&](Pythia* pythiaPtr) {
        // giving reference to the instance that generated the event.
		Event& event = pythiaPtr->event; 
		const Info& info = pythiaPtr->info;

        double pTHat  = info.pTHat();

        if (softQCD && iBin == 0 && info.isNonDiffractive()
        && pTHat > binEdges[iBin+1]) return;

        hardPt->Fill(pTHat);

        for (int i = 0; i < event.size(); ++i) {
            // look for forward muons
            auto* particle = &event[i];
            int particleID = abs(particle->id());
            if (particleID == 13 && particle->isFinal() && (particle->eta()>2.5) && (particle->eta()<4.0)) { // final state muons in forward rapidity region
                // find first non-copy mother
                auto* muonFinal = particle;
                //std::cout << "==== Here forward muon decay: mu";
                while ((particle->mother2() != 0) && (particle->mother1() == particle->mother2())) {
                    particle = &event[particle->mother1()];
                    //std::cout << " <- " << abs(particle->id());
                }
                auto* firstMother = &event[particle->mother1()];
                //std::cout << " <- " << abs(firstMother->id()) << " (firstMother)" << std::endl;

                muonTuple->Fill(iBin, iEvent, muonFinal->pAbs(), muonFinal->pT(), muonFinal->y(), muonFinal->eta(), abs(firstMother->id()));
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