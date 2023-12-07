#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
#include "TFile.h"
#include "TH1.h"

using namespace Pythia8;

int main() {
    // Generator
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain11_W+_100k.root", "RECREATE");

    // Total Cross Section
    TH1F *hardPt = new TH1F("SigmaGen","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", 100, 0.0, 100.0);

    TH1F *W_rapidity = new TH1F("W_mother_rapidity_sigma","W Mother of Muons Rapidity;y;#frac{d#sigma}{dy} (mb)", 100, -10.0, 10.0);
    TH1F *W_eta = new TH1F("W_mother_eta_sigma","W Mother of Muons Eta;y;#frac{d#sigma}{d#eta} (mb)", 100, -10.0, 10.0);

    TH1F *muonWPt = new TH1F("W_pt_sigma","W+ -> mu;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", 100, 0.0, 100.0);
    TH1F *muonWEta = new TH1F("W_eta_sigma","W+ -> mu;#eta;#frac{d#sigma}{d#eta} (mb)", 100, -10.0, 10.0);

    // read events from POWHEG lhe output
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = pwgevents_W_plus_100k.lhe");

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
    int iEvent = 0, iError = 0;
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

        for (int i = 0; i < pythia.event.size(); ++i) {
            int particleID = abs(pythia.event[i].id());
            if (particleID == 13) { // is muon
                int motherIndex = pythia.event[i].mother1();
                int firstMotherID = abs(pythia.event[motherIndex].id());
                if (firstMotherID == 24) { // from W+ decay
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    muonWPt->Fill(particlePt);
                    muonWEta->Fill(particleEta);
                    // For W mother
                    double W_mother_rapidity = pythia.event[motherIndex].y();
                    double W_mother_eta = pythia.event[motherIndex].eta();
                    W_rapidity->Fill(W_mother_rapidity);
                    W_eta->Fill(W_mother_eta);
                }
            }
        }

        hardPt->Fill(pTHat);

        ++iEvent;
    }

    pythia.stat();

    double luminocity_hard = iEvent/(pythia.info.sigmaGen());
    hardPt->Scale(1/luminocity_hard, "width");
    muonWPt->Scale(1/luminocity_hard, "width");
    muonWEta->Scale(1/luminocity_hard, "width");
    W_rapidity->Scale(1/luminocity_hard, "width");
    W_eta->Scale(1/luminocity_hard, "width");

    hardPt->Write();
    muonWPt->Write();
    muonWEta->Write();
    W_rapidity->Write();
    W_eta->Write();

    delete outFile;

    return 0;
}