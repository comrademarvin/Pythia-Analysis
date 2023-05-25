#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TThread.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

// Thread Handle function, from docs definition
void* handle(void* ptr) {
    int ith = (long)ptr;

    TFile* outFile = new TFile(Form("mymain07_%d.root", ith), "RECREATE");

    TTree* tree = new TTree(Form("tree_%d", ith), "tree");

    double pTHat;

    tree->Branch("pTHat", &pTHat, "pTHat/D");

    int N_events = 100;

    Pythia pythia;
    pythia.readString("HardQCD:all = on");
    pythia.readString("Beams:eCM = 5020.");
    pythia.readString("Tune:pp = 14");
    pythia.settings.parm("PhaseSpace:pTHatMin", (ith+1)*16.0);
    pythia.settings.parm("PhaseSpace:pTHatMax", (ith+2)*16.0);
    pythia.init();

    //TH1F *pTHard = new TH1F(Form("pt_hard_%d", ith),"", 10, 0.0, 20.0);

    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythia.next()) continue;

        pTHat = pythia.info.pTHat();
        
        tree->Fill();
    }

    // TCanvas *canvasTotal = new TCanvas(Form("total_sigma_%d", ith),"total_sigma");

    // pTHard->SetLineColor(1);
    // pTHard->Draw();

    // canvasTotal->Write();
    // canvasTotal->Close();

    outFile->Write();
    outFile->Close();

    delete tree;
    delete outFile;
}

int main() {
    // pT hard bins, spawn one thread for each bin    
    int nBins = 4;
    double binEdges[nBins+1] = {16.0, 26.0, 36.0, 50.0, 70.0};

    // Create and Spawn Threads
    TThread* th[nBins];

    for (int i = 0; i < nBins; i++) {
        th[i] = new TThread(Form("th%d", i), handle, (void*)i);
        th[i]->Run();
    }

    // while (true) {
    //     bool running = false;
    //     for (int i = 0; i < nBins; i++) {
    //         if (th[i]->GetState() == 2) {
    //             cout << "RUNNING\n";
    //             running = true;
    //         }
    //     }

    //     if (!running) break;
    //     cout << "REPEATING\n";
    // }

    for (int i = 0; i < nBins; i++) {
        th[i]->Join();
    }

    // for (int i = 0; i < nBins; i++) {
    //     delete th[i];
    // }

    // delete[] th;

    return 0;
}