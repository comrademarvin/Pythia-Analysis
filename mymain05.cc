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
    // ROOT output file
    TFile* outFile = new TFile("mymain05_200k.root", "RECREATE");

    // Number of events to generate per bin.
    int N_events = 200000;

    // Turn SoftQCD on/off
    bool softQCD = true;
    int softEventsScale = 1; // factor to generate more/less softQCD events

    // Pythia object
    Pythia pythia;

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

    TNtuple* genInfo = new TNtuple("genInfo", "genInfo", "weightSum:sigmaGen:sigmaErr"); // generated process info for cross-section normalisation

    // Histograms
    // Total Cross Section
    int N_bins_hard_pt = 50;
    TH1F *hardPt = new TH1F("pT_hat","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (mb/GeV/c)", N_bins_hard_pt, 0.0, 200.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", N_bins_hard_pt, 0.0, 200.0);

    // store forward muons for analysis
    vector<TNtuple*> muonTuples(nBins);

    for (int i = 0; i < nBins; ++i) {
        muonTuples[i] = new TNtuple(Form("muon%d", i), "muon", "binTag:eventTag:pAbs:pt:y:eta:firstMother");
    }

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin == 0) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
        } else {
            pythia.readString("HardQCD:all = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 5360.");
        pythia.readString("Tune:pp = 14");

        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        hardPtPart->Reset();

        int eventCount = 0;

        int events_run = N_events; // run more events for SoftQCD and first bin for improved statistics
        if (softQCD && iBin == 0) {
            events_run = softEventsScale*N_events;
        }

        for (int iEvent = 0; iEvent < events_run; ++iEvent) {
            if (!pythia.next()) continue;

            double pTHat = pythia.info.pTHat();

            if (softQCD && iBin == 0 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            eventCount++;

            hardPtPart->Fill(pTHat);

            for (int i = 0; i < pythia.event.size(); ++i) {
                // look for forward muons
                auto* particle = &pythia.event[i];
                int particleID = abs(particle->id());
                if (particleID == 13 && particle->isFinal() && (particle->eta()>2.5) && (particle->eta()<4.0)) { // final state muons in forward rapidity region
                    // find first non-copy mother
                    auto* muonFinal = particle;
                    std::cout << "==== Here forward muon decay: mu";
                    while ((particle->mother2() != 0) && (particle->mother1() == particle->mother2())) {
                        particle = &pythia.event[particle->mother1()];
                        std::cout << " <- " << abs(particle->id());
                    }
                    auto* firstMother = &pythia.event[particle->mother1()];
                    std::cout << " <- " << abs(firstMother->id()) << " (firstMother)" << std::endl;

                    muonTuples[iBin]->Fill(iBin, iEvent, muonFinal->pAbs(), muonFinal->pT(), muonFinal->y(), muonFinal->eta(), abs(firstMother->id()));
                }
            }
        }
        
        // generated cross-section info for the bin
        genInfo->Fill(pythia.info.weightSum(), pythia.info.sigmaGen(), pythia.info.sigmaErr());

        // normalise pT-hat for the bin
        hardPtPart->Scale(1/(pythia.info.weightSum()), "width");
        TH1F *sigmaGenHist = new TH1F("hard_QCD_sigma_gen","", N_bins_hard_pt, 0.0, 200.0);
        for (int iBin = 0; iBin < N_bins_hard_pt; iBin++) {
            sigmaGenHist->SetBinContent(iBin, pythia.info.sigmaGen());
            sigmaGenHist->SetBinError(iBin, pythia.info.sigmaErr());
        }
        hardPtPart->Multiply(sigmaGenHist);

        hardPtPart->Write(Form("hard_pt_part_%d", iBin));

        hardPt->Add(hardPtPart);
    }

    genInfo->Write();

    // TCanvas *canvasPtHat = new TCanvas("pT_hat_full","pT_hat_full"); // contributions from soft/hard QCD to pT-hat
    // gPad->SetLogy();
    // hardPt->SetMinimum(0.000000001);
    // hardPt->Draw();
    // canvasPtHat->Write();

    hardPt->Write();

    // Read Tuples to output root file
    for (int i = 0; i < nBins; ++i) {
        muonTuples[i]->Write(Form("muon%d", i));
    }

    delete outFile;

    return 0;
}