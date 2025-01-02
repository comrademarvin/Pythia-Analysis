#include "Pythia8/Pythia.h"
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include <TNtuple.h>

using namespace Pythia8;

int main() {
    static const int nBins = 6;
    static const double binEdges[nBins+1] = {0.0, 15.0, 30.0, 50.0, 70.0, 100.0, 150.0};

    // generated process info for cross-section normalisation
    TNtuple* genInfo = new TNtuple("genInfo", "genInfo", "weightSum:sigmaGen:sigmaErr");

    // Histograms
    // Total Cross Section
    int N_bins_hard_pt = 150;
    double hard_pt_max = 150.0;
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_hard_pt, 0.0, hard_pt_max);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", N_bins_hard_pt, 0.0, hard_pt_max);

    // HF Cross Sections
    vector<TNtuple*> chargedTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        chargedTuples[i] = new TNtuple(Form("charged%d", i), "charged", "binTag:eventTag:pAbs:pt:y:eta:id:charge"); // see definition of decayStaus below
    }

    Pythia pythia;

    pythia.readString("Beams:eCM = 5360.");
    pythia.readString("Tune:pp = 14"); // Monash Tune

    // MPI
    pythia.readString("PartonLevel:MPI = off");

    // colour reconnection
    //pythia.readString("ColourReconnection:reconnect = off");

    int N_events = 1000000;

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

        int Event_count = 0;

        // event loop
        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat  = pythia.info.pTHat();

            if (iBin == 0 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            hardPtPart->Fill(pTHat);

            Event_count++;

            // iterate through event record
            for (int i = 0; i < pythia.event.size(); ++i) {
                auto* particle = &pythia.event[i];
                if (particle->isCharged()) { // in 4pi region
                    // only consider primary charged particles (According to ALICE's definition of primary)
                    double particleLifetime = (particle->tau())/10; // Convert mm/c to cm/c
                    bool isPrimary = true;
                    if (particleLifetime < 1) {
                        isPrimary = false;
                    } else {
                        int motherIndex = particle->mother1();
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

                    // only consider primary charged particles in the region of interest
                    if (isPrimary) {
                        int particleID = particle->id();
                        double particlePseudorapidity = particle->eta();
                        double particlePAbs = particle->pAbs();
                        double particlePt = particle->pT();
                        double particleRapidity = particle->y();
                        double particleCharge = particle->charge();

                        chargedTuples[iBin]->Fill(iBin, iEvent, particlePAbs, particlePt, particleRapidity, particlePseudorapidity, particleID, particleCharge);
                    }
                }
            }

        }

        // generated cross-section info for the bin
        genInfo->Fill(pythia.info.weightSum(), pythia.info.sigmaGen(), pythia.info.sigmaErr());

        // normalise pT-hat for the bin
        hardPtPart->Scale(1/(pythia.info.weightSum()), "width");
        TH1F *sigmaGenHist = new TH1F("hard_QCD_sigma_gen","", N_bins_hard_pt, 0.0, hard_pt_max);
        for (int iBin = 0; iBin < N_bins_hard_pt; iBin++) {
            sigmaGenHist->SetBinContent(iBin, pythia.info.sigmaGen());
            sigmaGenHist->SetBinError(iBin, pythia.info.sigmaErr());
        }
        hardPtPart->Multiply(sigmaGenHist);
        hardPt->Add(hardPtPart);
    }

    TFile* outFile = new TFile("mymain01_1M_536_4pi_MPI_off.root", "RECREATE");

    hardPt->Write();

    // Read Tuples to output root file
    for (int iBin = 0; iBin < nBins; ++iBin) {
        chargedTuples[iBin]->Write(Form("charged%d", iBin));
    }

    genInfo->Write();

    delete outFile;

    return 0;
}