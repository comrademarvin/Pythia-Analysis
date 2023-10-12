#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>

using namespace Pythia8;

int main() {
    static const int nBins = 5;
    static const double binEdges[nBins+1] = {0.0, 10.0, 30.0, 50.0, 75.0, 100.0};

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 100, 0.0, 100.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 100, 0.0, 100.0);

    // HF Cross Sections
    vector<TNtuple*> chargedTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        chargedTuples[i] = new TNtuple(Form("charged%d", i), "charged", "binTag:eventTag:pAbs:pt:y:eta:id:charge"); // see definition of decayStaus below
    }

    Pythia pythia;

    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("Tune:pp = 14");
    pythia.readString("PDF:pSet = 9");

    int N_events = 50000;

    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (iBin == 0) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            pythia.readString("HardQCD:all = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        // event loop
        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat  = pythia.info.pTHat();

            if (iBin == 0 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            hardPtPart->Fill(pTHat);

            // iterate through event record
            for (int i = 0; i < pythia.event.size(); ++i) {
                if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) {
                    int particleID = pythia.event[i].id();
                    double particlePAbs = pythia.event[i].pAbs();
                    double particlePt = pythia.event[i].pT();
                    double particleRapidity = pythia.event[i].y();
                    double particlePseudorapidity = pythia.event[i].eta();
                    double particleCharge = pythia.event[i].charge();

                    chargedTuples[iBin]->Fill(iBin, iEvent, particlePAbs, particlePt, particleRapidity, particlePseudorapidity, particleID, particleCharge);
                }
            }

        }

        double luminocity = N_events/(pythia.info.sigmaGen());
        binLuminocity[iBin] = luminocity;

        hardPtPart->Scale(1/luminocity, "width");
        hardPt->Add(hardPtPart);
    }

    TFile* outFile = new TFile("mymain01.root", "RECREATE");

    hardPt->Write();

    // Read Tuples to output root file
    for (int iBin = 0; iBin < nBins; ++iBin) {
        chargedTuples[iBin]->Write(Form("charged%d", iBin));
    }

    outFile->WriteObject(&binLuminocity, "luminocity");

    delete outFile;

    return 0;
}