#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
//#include "Pythia8Plugins/HepMC2.h"
#include "TFile.h"
#include "TH1.h"
#include <TNtuple.h>

using namespace Pythia8;

int main() {
    // Generator
    Pythia pythia;

    int generatedEvents = 10000;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain11_W+_10k.root", "RECREATE");

    // Total Cross Section
    TH1F *hardPt = new TH1F("SigmaGen","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", 100, 0.0, 100.0);
    vector<double> generatedLuminocity(1);

    // Store muons from W for analysis
    TNtuple* muonTuple = new TNtuple("muon", "muon", "eventTag:pAbs:pt:y:eta");

    // Store charged particles and event charged particle count for multiplicity analysis
    // TNtuple* chargedTuple = new TNtuple("charged", "charged", "eventTag:pAbs:pt:y:eta");
    vector<int> multiplicities(generatedEvents, -1);

    // Store W mother for physics kinematics
    //TNtuple* WMotherTuple = new TNtuple("W_mother", "W_mother", "eventTag:pAbs:pt:y:eta");

    // read events from POWHEG lhe output
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = pwgevents_Josh.lhe");
    pythia.readString("Tune:pp = 14"); // Monash tune
    //pythia.readString("Parallelism:numThreads = 5");

    // Event settings
    pythia.readString("Main:numberOfEvents = 0");
    int nError = 10;
    pythia.readString("Main:timesAllowErrors = 10");

    // POWHEG settings
    pythia.readString("POWHEG:veto = 1");
    pythia.readString("POWHEG:pTdef = 1");
    pythia.readString("POWHEG:emitted = 0");
    pythia.readString("POWHEG:pTemt = 0");
    pythia.readString("POWHEG:pThard = 0");
    pythia.readString("POWHEG:vetoCount = 3");
    pythia.readString("POWHEG:MPIveto = 0");

    // Add in user hooks for shower vetoing.
    shared_ptr<PowhegHooks> powhegHooks;
    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);

    // Initialise and list settings
    pythia.init();

    // Begin event loop; generate until end of LHEF file
    int iEvent = 0, iError = 0, muonFound = 0;
    while (true) {
        // Generate the next event
        if (!pythia.next()) {
            // If failure because reached end of file then exit event loop
            if (pythia.info.atEndOfFile()) break;

            // Otherwise count event failure and continue/exit as necessary
            cout << "Warning: event " << iEvent << " failed" << endl;
            if (++iError == nError) {
                cout << "Error: too many event failures... exiting" << endl;
                break;
            }

            continue;
        }

        // process the event

        double pTHat  = pythia.info.pTHat();
        
        int nCharged = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            int motherIndex;
            if (pythia.event[i].isCharged()) {
                // only consider primary charged particles (According to ALICE's definition of primary)
                double particleLifetime = (pythia.event[i].tau())/10; // Convert mm/c to cm/c
                bool isPrimary = true;
                if (particleLifetime < 1) {
                    isPrimary = false;
                } else {
                    motherIndex = pythia.event[i].mother1();
                    double motherLifetime;
                
                    while (!pythia.event[motherIndex].isQuark() && !pythia.event[motherIndex].isGluon())  {
                        motherLifetime = (pythia.event[motherIndex].tau())/10;
                        if (motherLifetime > 1) {
                            isPrimary = false;
                            break;
                        }
                        motherIndex = pythia.event[motherIndex].mother1();
                    }
                }

                // only consider primary charged particles in the central barrel region
                if (isPrimary && (abs(pythia.event[i].eta()) < 1)) {
                    nCharged++;
                }
            }

            int particleID = abs(pythia.event[i].id());
            if (particleID == 13) { // is muon
                motherIndex = pythia.event[i].mother1();
                int firstMotherID = abs(pythia.event[motherIndex].id());
                if (firstMotherID == 24) { // from W decay
                    // For muon
                    double particlePAbs = pythia.event[i].pAbs();
                    double particlePt = pythia.event[i].pT();
                    double particleRapidity = pythia.event[i].y();
                    double particlePseudorapidity = pythia.event[i].eta();

                    muonTuple->Fill(iEvent, particlePAbs, particlePt, particleRapidity, particlePseudorapidity);
                    muonFound++;

                    // For W mother
                    // double W_mother_rapidity = pythia.event[motherIndex].y();
                    // double W_mother_eta = pythia.event[motherIndex].eta();
                }
            }
        }

        hardPt->Fill(pTHat);
        multiplicities[iEvent] = nCharged;

        ++iEvent;
    }

    //pythia.stat();

    std::cout << "Count of Events: " << iEvent << std::endl;
    std::cout << "Count of W->mu found: " << muonFound << std::endl;

    double luminocity = iEvent/(pythia.info.sigmaGen());
    generatedLuminocity[0] = luminocity;

    hardPt->Scale(1/luminocity, "width");
    hardPt->Write();

    muonTuple->Write("W_muon");
    outFile->WriteObject(&generatedLuminocity, "luminocity");
    outFile->WriteObject(&multiplicities, "event_multiplicity");


    delete outFile;

    return 0;
}