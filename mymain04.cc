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
    pythiaHard.readString("Beams:eCM = 13700.");
    pythiaHard.readString("Tune:pp = 14");
    pythiaHard.readString("HardQCD:all = on");
    pythiaHard.readString("SoftQCD:nonDiffractive = off");
    pythiaHard.init();

    // Soft QCD
    Pythia pythiaSoft;
    pythiaSoft.readString("Beams:eCM = 13700.");
    pythiaSoft.readString("Tune:pp = 14");
    pythiaSoft.readString("HardQCD:all = off");
    pythiaSoft.readString("SoftQCD:nonDiffractive = on");
    pythiaSoft.init();

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain04.root", "RECREATE");

    // TNtuple appraoch for cuts and histograms
    TNtuple* hardestHardQCD = new TNtuple("hardQCD", "hardQCD", "pTHat:thetaHat:phiHat");
    TNtuple* hardestSoftQCD = new TNtuple("softQCD", "softQCD", "pTHat:thetaHat:phiHat");

    int N_events = 200000;

    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythiaHard.next()) continue;
        if (!pythiaSoft.next()) continue;

        hardestHardQCD->Fill(pythiaHard.info.pTHat(), pythiaHard.info.thetaHat(), pythiaHard.info.phiHat());
        hardestSoftQCD->Fill(pythiaSoft.info.pTHat(), pythiaSoft.info.thetaHat(), pythiaSoft.info.phiHat());
    }

    pythiaHard.stat();
    pythiaSoft.stat();

    // Luminocity for normalization
    float luminocity_hard = N_events/(pythiaHard.info.sigmaGen()*pow(10,9));
    float luminocity_soft = N_events/(pythiaSoft.info.sigmaGen()*pow(10,9));

    //Plotting
    TCanvas *c1 = new TCanvas("c1","c1");

    // Hard QCD hist
    TH1F *hardQCDpTHat = new TH1F("hard_QCD_pTHat","Contribution to Hardest Process;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 40, 0, 40);
    hardQCDpTHat->SetLineColor(1);
    hardQCDpTHat->SetStats(0);
    hardestHardQCD->Draw("pTHat>>hard_QCD_pTHat", "pTHat<40");
    hardQCDpTHat->Scale(1/luminocity_hard, "width");

    // Soft QCD hist
    TH1F *softQCDpTHat = new TH1F("soft_QCD_pTHat", "Contribution to hardest process;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 40, 0, 40);
    softQCDpTHat->SetLineColor(2);
    softQCDpTHat->SetStats(0);
    hardestSoftQCD->Draw("pTHat>>soft_QCD_pTHat", "pTHat<40", "SAME");
    softQCDpTHat->Scale(1/luminocity_soft, "width");

    // Legend
    auto legend = new TLegend();
    legend->AddEntry(hardQCDpTHat,"HardQCD:all","l");
    legend->AddEntry(softQCDpTHat,"SoftQCD:nonDiffractive","l");
    legend->Draw();

    c1->Write();

    delete outFile;

    return 0;
}