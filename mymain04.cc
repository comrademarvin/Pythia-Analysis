#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

int main() {
    // Hard QCD
    Pythia pythiaHard;
    pythiaHard.readString("Beams:eCM = 13600.");
    pythiaHard.readString("Tune:pp = 14");
    pythiaHard.readString("HardQCD:all = on");
    pythiaHard.readString("SoftQCD:nonDiffractive = off");
    // pythiaHard.readString("PartonLevel:all = off");
    // pythiaHard.readString("PhaseSpace:pTHatMinDiverge = 8.0");
    pythiaHard.readString("PhaseSpace:pTHatMin = 2.0");
    pythiaHard.init();

    // Soft QCD
    Pythia pythiaSoft;
    pythiaSoft.readString("Beams:eCM = 13600.");
    pythiaSoft.readString("Tune:pp = 14");
    pythiaSoft.readString("HardQCD:all = off");
    pythiaSoft.readString("SoftQCD:nonDiffractive = on");
    // pythiaSoft.readString("PartonLevel:all = on");
    // pythiaSoft.readString("PartonLevel:ISR = off");
    // pythiaSoft.readString("PartonLevel:FSR = off");
    // pythiaSoft.readString("HadronLevel:all = off");
    pythiaSoft.init();

    // TNtuple appraoch for cuts and histograms
    TNtuple* hardestHardQCD = new TNtuple("hardQCD", "hardQCD", "pTHat:thetaHat:phiHat");
    TNtuple* hardestSoftQCD = new TNtuple("softQCD", "softQCD", "pTHat:thetaHat:phiHat");

    // generation info for normalisation
    vector<double> genInfoHard(3);
    vector<double> genInfoSoft(3);

    int N_events = 1000000;

    int eventCounter = 0;

    // generate events for Hard QCD
    eventCounter = 0;
    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythiaHard.next()) continue;
        
        hardestHardQCD->Fill(pythiaHard.info.pTHat(), pythiaHard.info.thetaHat(), pythiaHard.info.phiHat());
        eventCounter++;
    }

    // generate events for Soft QCD
    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythiaSoft.next()) continue;

        hardestSoftQCD->Fill(pythiaSoft.info.pTHat(), pythiaSoft.info.thetaHat(), pythiaSoft.info.phiHat());
        eventCounter++;
    }

    // pythiaHard.stat();
    // pythiaSoft.stat();

    // Luminocity for normalization
    genInfoHard[0] = pythiaHard.info.weightSum();
    genInfoHard[1] = pythiaHard.info.sigmaGen();
    genInfoHard[2] = pythiaHard.info.sigmaErr();

    genInfoSoft[0] = pythiaSoft.info.weightSum();
    genInfoSoft[1] = pythiaSoft.info.sigmaGen();
    genInfoSoft[2] = pythiaSoft.info.sigmaErr();

    // float luminocity_hard = pythiaHard.info.weightSum()/(pythiaHard.info.sigmaGen());
    // float luminocity_soft = pythiaSoft.info.weightSum()/(pythiaSoft.info.sigmaGen());

    // check the error from the cross-section
    std::cout << "Generated cross-section for Soft QCD: " << pythiaSoft.info.sigmaGen() << " +/- " << pythiaSoft.info.sigmaErr() << std::endl;
    std::cout << "Generated cross-section for Hard QCD: " << pythiaHard.info.sigmaGen() << " +/- " << pythiaHard.info.sigmaErr() << std::endl;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain04.root", "RECREATE");

    hardestHardQCD->Write();
    hardestSoftQCD->Write();

    outFile->WriteObject(&genInfoHard, "genInfoHard");
    outFile->WriteObject(&genInfoSoft, "genInfoSoft");

    delete outFile;

    return 0;
}