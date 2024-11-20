#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
#include "TFile.h"
#include "TH1.h"
#include <TNtuple.h>

using namespace Pythia8;

int main() {
    // Generator
    Pythia pythia;

    int generatedEvents = 10000000;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain11_W+_10M_forward.root", "RECREATE");

    // Total Cross Section
    TH1F *hardPt = new TH1F("SigmaGen","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", 100, 0.0, 100.0);
    vector<double> genInfo(3);

    // Store muons from W for analysis
    TNtuple* muonTuple = new TNtuple("muon", "muon", "eventTag:pAbs:pt:y:eta");

    // Store charged particles and event charged particle count for multiplicity analysis
    vector<int> multiplicities(generatedEvents, -1);
    bool multCentral = false; // either estimate multiplicity in the central or forward region
    float multEtaMin, multEtaMax;
    if (multCentral) {
        multEtaMin = -1.0; // range of SPD in ITS
        multEtaMax = 1.0;
    } else {
        multEtaMin = 2.5; // range of V0 (2.8 < eta < 5.1)
        multEtaMax = 4.0;
    }

    // Store W mother for physics kinematics
    //TNtuple* WMotherTuple = new TNtuple("W_mother", "W_mother", "eventTag:pAbs:pt:y:eta");

    // read events from POWHEG lhe output
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = pwgevents_W+_10M_536.lhe");
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
            auto* particle = &pythia.event[i];
            int particleID = abs(particle->id());

            int motherIndex;
            if (particle->isCharged() && (particle->eta() > multEtaMin) && (particle->eta() < multEtaMax)) {
                // only consider primary charged particles (According to ALICE's definition of primary)
                double particleLifetime = (particle->tau())/10; // Convert mm/c to cm/c
                bool isPrimary = true;
                if (particleLifetime < 1) {
                    isPrimary = false;
                } else {
                    motherIndex = particle->mother1();
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
                if (isPrimary) nCharged++;
            }

            if (particleID == 13 && particle->isFinal()) { // final state muons in forward rapidity region
                // find first non-copy or muon mother
                auto* muonFinal = particle;
                auto* firstMother = particle;
                string decayOutput = "==== Here forward muon decay: mu";
                while (((firstMother->mother2() != 0) && (firstMother->mother1() == firstMother->mother2())) || (abs(firstMother->id()) == 13)) {
                    firstMother = &pythia.event[firstMother->mother1()];
                    string partOutput = " <- " + firstMother->name();
                    decayOutput = decayOutput + partOutput;
                }
                decayOutput = decayOutput + " (firstMother)";
                decayOutput = decayOutput + " <- ";
                decayOutput = decayOutput + pythia.event[firstMother->mother1()].name();
                //std::cout << decayOutput << std::endl;

                int firstMotherID = abs(firstMother->id());
                if (firstMotherID == 24) { // from W decay
                    muonTuple->Fill(iEvent, muonFinal->pAbs(), muonFinal->pT(), muonFinal->y(), muonFinal->eta());
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

    std::cout << "Count of Events: " << iEvent << std::endl;
    std::cout << "Count of W->mu found: " << muonFound << std::endl;

    // Luminocity for normalization
    genInfo[0] = pythia.info.weightSum();
    genInfo[1] = pythia.info.sigmaGen();
    genInfo[2] = pythia.info.sigmaErr();

    hardPt->Scale(genInfo[1]/genInfo[0], "width");
    hardPt->Write();

    muonTuple->Write("W_muon");
    outFile->WriteObject(&genInfo, "genInfo");
    outFile->WriteObject(&multiplicities, "event_multiplicity");

    delete outFile;

    return 0;
}