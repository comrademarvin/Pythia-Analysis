#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

bool motherHF(int motherID, Pythia8::Event eventObj) {
    if (motherID == 0) return false;

    if (eventObj[motherID].mother1() == 0) return false;

    if ((eventObj[motherID].id() == 4) || (eventObj[motherID].id() == -4)) return true; // is charm quark

    if ((eventObj[motherID].id() == 5) || (eventObj[motherID].id() == -5)) return true; // is bottom quark

    if (eventObj[motherID].mother1() != eventObj[motherID].mother2()) {
        return (motherHF(eventObj[motherID].mother1(), eventObj) || motherHF(eventObj[motherID].mother2(), eventObj));
    } else {
        return motherHF(eventObj[motherID].mother1(), eventObj);
    }
}

int main() {
    Pythia pythia;

    // Collision settings/parameters
    pythia.readString("HardQCD:all = on");
    // pythia.readString("HardQCD:hardccbar = on");
    // pythia.readString("HardQCD:hardbbbar = on");

    // Collision energy
    pythia.readString("Beams:eCM = 5020.");

    // // PDF
    // pythia.readString("PDF:pSet = 13"); // 9 - CTEQ66.00 (NLO, central member); 13 - NNPDF2.3 QCD+QED LO

    // // MPI
    // pythia.readString("MultipartonInteractions:processLevel = 3"); // also charmonium and bottomonium production, via colour singlet and colour octet channels.
    // //pythia.readString("MultipartonInteractions:pT0Ref = 2.30"); // used by most examples
    // //pythia.readString("MultipartonInteractions:pTmaxMatch = 2"); // always allow multiparton interactions up to the kinematical limit.

    // // Colour reconnection
    // pythia.readString("ColourReconnection:mode = 0"); // default
    // //pythia.readString("ColourReconnection:range = 1.8"); // The higher this number is the more reconnections can occur.

    // Monash 13
    pythia.readString("Tune:pp = 14");

    // initialize
    pythia.init();

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain02.root", "RECREATE");

    // Event Multiplicity counts
    TH1F *mult = new TH1F("mult","charged multiplicity", 200, -1, 399);

    // TNtuple appraoch for cuts and histograms
    TNtuple* muTuple = new TNtuple("hf_muons", "mu", "pt:y:eta");

    // Generate events
    int N_events = 10000;
    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythia.next()) continue;

        int multCount = 0;
        int muCount = 0;
        int muHFCount = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) multCount++;

            if ((pythia.event[i].id() == 13) || (pythia.event[i].id() == -13) && (pythia.event[i].isFinal())) { // is muon
                muCount++;
                cout << "Muon Parents: " << pythia.event[pythia.event[i].mother1()].name() << ", " << pythia.event[pythia.event[i].mother2()].name() << endl;
                // check if it comes from a HF decay
                bool HF_parent = motherHF(i, pythia.event);
                if (HF_parent) {
                    muTuple->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    //cout << "Found HF decay Muon!" << endl;
                    muHFCount++;
                } 
            }
        }

        if (muCount > 0) {
           cout << "Number of final state muons in event: " << muCount << endl;
           pythia.event.list();
           cout << endl;
        }
        // if (muCount != muHFCount) {
        //     cout << "FOUND non-HF MUON!" << endl;
        //     pythia.event.list();
        // }

        mult->Fill( multCount );
    }

    pythia.stat();

    const Int_t NBINS = 17;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40};
    TH1F *muPt = new TH1F("mu_pt","HF muon decay pt;Momentum (GeV/c);d#sigma/dpT (pb/GeV/c)", NBINS, edges);
    //muTuple->Draw("pt>>mu_pt", "pt<20 && pt>2 && y>2.5 && y<4");
    muTuple->Draw("pt>>mu_pt", "pt<40");

    TH1F *muEta = new TH1F("mu_eta","HF muon decay eta;#eta;d#sigma/d#eta", 40, -10, 10);   
    //muTuple->Draw("eta>>mu_eta", "pt<20 && pt>2 && y>2.5 && y<4");
    muTuple->Draw("eta>>mu_eta", "pt<40");

    TH1F *muY = new TH1F("mu_y","HF muon decay rapidity;y;d#sigma/dy (pb)", 24, -6, 6);   
    muTuple->Draw("y>>mu_y", "pt<40");

    // get process cross-section from counts
    muPt->Scale((pythia.info.sigmaGen()*pow(10,9))/N_events);
    muEta->Scale((pythia.info.sigmaGen()*pow(10,9))/N_events);
    muY->Scale((pythia.info.sigmaGen()*pow(10,9))/N_events);

    mult->Write();
    muPt->Write();
    muEta->Write();
    muY->Write();
    //muTuple->Write();

    delete outFile;

    return 0;
}