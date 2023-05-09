#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

// Prints condensed version of particle mother-tree
void decayStatus(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int momIndex: eventObj[particleIndex].motherList()) {
            cout << " " << eventObj[momIndex].name() << " (" << momIndex << ", " << eventObj[momIndex].id() << ", " << eventObj[momIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int momIndex: eventObj[particleIndex].motherList()) {
            if ((!eventObj[momIndex].isGluon()) && (abs(eventObj[momIndex].id())!=2212)) {
                decayStatus(momIndex, eventObj, decayLevel);
            }
        }
    }
}

// Prints condensed version of particle daughter-tree
void decayStatusDaughters(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            cout << " " << eventObj[daughterIndex].name() << " (" << daughterIndex << ", " << eventObj[daughterIndex].id() << ", " << eventObj[daughterIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            if (!eventObj[daughterIndex].isFinal()) {
                decayStatusDaughters(daughterIndex, eventObj, decayLevel);
            }
        }
    }
}

// returns the status of the earliest HF mother of hadron
int hadronHFStatus(int hadronIndex, Pythia8::Event eventObj, int hfID = 4) {
    int hfIndex = -1;  
    for (int momIndex: eventObj[hadronIndex].motherList()) {
        if (abs(eventObj[momIndex].id()) == hfID) {
            hfIndex = momIndex;
            continue;
        }
    }

    if (hfIndex == -1) return 0; // hadron does not contain requested hf

    while (true) {
        bool charmFound = false;
        for (int momIndex: eventObj[hfIndex].motherList()) {
            if (abs(eventObj[momIndex].id()) == hfID) {
                hfIndex = momIndex;
                charmFound = true;
                continue;
            }
        }

        if (!charmFound) {
            return eventObj[hfIndex].status();
        }
    }
}

// returns index of decay muon in event, returns -1 if no decay muon is found
int muonDecay(int index, Pythia8::Event eventObj) {
    if (abs(eventObj[index].id()) == 13) return index;

    for (int daughterIndex: eventObj[index].daughterList()) {
        int returnIndex = muonDecay(daughterIndex, eventObj);
        if (returnIndex != -1) return returnIndex;
    }

    return -1;
}

int main() {
    Pythia pythiaHard;

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain05.root", "RECREATE");

    // pTHat bins
    int nBins = 5;
    double binEdges[nBins+1] = {0.0, 9.0, 20.0, 30.0, 40.0, 60.0};

    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 60, 0.0, 60.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 60, 0.0, 60.0);

    // HF Cross Sections
    TH1F *charmPt = new TH1F("charm_full","Heavy Flavour Events Cross-Section (Hardest Process);p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 30, 0.0, 60.0);
    TH1F *charmPtPart = new TH1F("charm_part","", 30, 0.0, 60.0);

    TH1F *bottomPt = new TH1F("bottom_full","", 30, 0.0, 60.0);
    TH1F *bottomPtPart = new TH1F("bottom_part","", 30, 0.0, 60.0);

    // HF Meson Cross Sections
    TH1F *DPt = new TH1F("D_full","D Meson Total Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 30, 0.0, 60.0);
    TH1F *DPtPart = new TH1F("D_part","", 30, 0.0, 60.0);

    TH1F *DHardPt = new TH1F("D_hard_full","", 30, 0.0, 60.0);
    TH1F *DHardPtPart = new TH1F("D_hard_part","", 30, 0.0, 60.0);

    TH1F *DShowerPt = new TH1F("D_shower_full","", 30, 0.0, 60.0);
    TH1F *DShowerPtPart = new TH1F("D_shower_part","", 30, 0.0, 60.0);
    
    TH1F *BPt = new TH1F("B_full","B Meson Total Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 30, 0.0, 60.0);
    TH1F *BPtPart = new TH1F("B_part","", 30, 0.0, 60.0);

    TH1F *BHardPt = new TH1F("B_hard_full","", 30, 0.0, 60.0);
    TH1F *BHardPtPart = new TH1F("B_hard_part","", 30, 0.0, 60.0);

    TH1F *BShowerPt = new TH1F("B_shower_full","", 30, 0.0, 60.0);
    TH1F *BShowerPtPart = new TH1F("B_shower_part","", 30, 0.0, 60.0);

    // Muon Decay Cross Sections
    TH1F *muPt = new TH1F("Decay Muons","Decay Muons Cross-Section;p_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 30, 0.0, 60.0);
    TH1F *muPtPart = new TH1F("mu_part","", 30, 0.0, 60.0);

    TH1F *muHFPt = new TH1F("mu_HF_hadron","", 30, 0.0, 60.0);
    TH1F *muHFPtPart = new TH1F("mu_HF_hadron_part","", 30, 0.0, 60.0);

    TH1F *muDPt = new TH1F("mu_D_hadron","", 30, 0.0, 60.0);
    TH1F *muDPtPart = new TH1F("mu_D_hadron_part","", 30, 0.0, 60.0);

    TH1F *muBPt = new TH1F("mu_B_hadron","", 30, 0.0, 60.0);
    TH1F *muBPtPart = new TH1F("mu_B_hadron_part","", 30, 0.0, 60.0);

    // TH1F *muOniumPt = new TH1F("mu_Onium_hadron","", 30, 0.0, 60.0);
    // TH1F *muOniumPtPart = new TH1F("mu_Onium_hadron_part","", 30, 0.0, 60.0);

    // Number of events to generate per bin.
    int N_events = 100000;

    for (int iBin = 0; iBin < nBins; ++iBin) {
        // set pythiaHard initialization variables
        if (iBin == 0) {
            pythiaHard.readString("HardQCD:all = off");
            pythiaHard.readString("SoftQCD:nonDiffractive = on");
        } else {
            pythiaHard.readString("HardQCD:all = on");
            // pythiaHard.readString("HardQCD:hardccbar = on");
            // pythiaHard.readString("HardQCD:hardbbbar = on");
            pythiaHard.readString("SoftQCD:nonDiffractive = off");
        }

        pythiaHard.readString("Beams:eCM = 5020.");
        pythiaHard.readString("Tune:pp = 14");
        pythiaHard.readString("411:onMode=off");
        pythiaHard.readString("411:onIfAny=13");
        pythiaHard.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythiaHard.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythiaHard.init();

        hardPtPart->Reset();
        charmPtPart->Reset();
        bottomPtPart->Reset();
        DPtPart->Reset();
        DHardPtPart->Reset();
        DShowerPtPart->Reset();
        BPtPart->Reset();
        BHardPtPart->Reset();
        BShowerPtPart->Reset();
        muPtPart->Reset();
        muHFPtPart->Reset();
        muDPtPart->Reset();
        muBPtPart->Reset();
        //muOniumPtPart->Reset();

        for (int iEvent = 0; iEvent < N_events; ++iEvent) {
            if (!pythiaHard.next()) continue;
            //cout << "====Start of new event====" << endl;

            double pTHat  = pythiaHard.info.pTHat();
            if (iBin == 0 && pythiaHard.info.isNonDiffractive()
                && pTHat > binEdges[1]) continue;

            hardPtPart->Fill(pTHat);

            bool charmAdded = false;
            bool bottomAdded = false;
            for (int i = 0; i < pythiaHard.event.size(); ++i) {
                int particleID = pythiaHard.event[i].id();
                int particleStatus = pythiaHard.event[i].status();

                // heavy flavours
                if ((abs(particleID) == 4) && (!charmAdded)) { // charm
                    if ((abs(particleStatus) == 23) || (abs(particleStatus) == 33)) { // outgoing from hardest process
                        charmPtPart->Fill(pythiaHard.event[i].pT());
                        charmAdded = true;
                        //decayStatus(i, pythiaHard.event);
                    }
                }

                if ((abs(particleID) == 5) && (!bottomAdded)) { // bottom
                    if ((abs(particleStatus) == 23) || (abs(particleStatus) == 33)) { // outgoing from hardest process
                        bottomPtPart->Fill(pythiaHard.event[i].pT());
                        bottomAdded = true;
                        //decayStatus(i, pythiaHard.event);
                    }
                }

                // hadrons
                if ((abs(particleStatus) >= 81) && (abs(particleStatus) <= 89)) { // primary hadrons produced by hadronization process
                    int charmStatus = abs(hadronHFStatus(i, pythiaHard.event));
                    bool containsCharm = true;
                    if ((charmStatus == 0) || (charmStatus >= 60)) containsCharm = false;

                    int bottomStatus = abs(hadronHFStatus(i, pythiaHard.event, 5));
                    bool containsBottom = true;
                    if ((bottomStatus == 0) || (bottomStatus >= 60)) containsBottom = false;

                    if (containsCharm || containsBottom) {
                        // look for muons decay
                        int muonIndex = muonDecay(i, pythiaHard.event);
                        if (muonIndex != -1) {
                            muHFPtPart->Fill(pythiaHard.event[muonIndex].pT());

                            if ((abs(particleID) >= 411) && (abs(particleID) <= 435)) { // D meson
                                muDPtPart->Fill(pythiaHard.event[muonIndex].pT());
                            }

                            else if ((abs(particleID) >= 511) && (abs(particleID) <= 545)) { // B meson
                                muBPtPart->Fill(pythiaHard.event[muonIndex].pT());
                            }

                            // else if (((abs(particleID) >= 441) && (abs(particleID) <= 445)) ||
                            //         ((abs(particleID) >= 551) && (abs(particleID) <= 557))) { // Quarkonia
                            //     muOniumPtPart->Fill(pythiaHard.event[muonIndex].pT());
                            // }
                        }

                        // look at D,B hadrons seperately
                        if ((abs(particleID) >= 411) && (abs(particleID) <= 435)) { // D meson
                            DPtPart->Fill(pythiaHard.event[i].pT());

                            if (charmStatus < 40) DHardPtPart->Fill(pythiaHard.event[i].pT());
                            else if (charmStatus < 60) DShowerPtPart->Fill(pythiaHard.event[i].pT());
                        }

                        else if ((abs(particleID) >= 511) && (abs(particleID) <= 545)) { // B meson
                            BPtPart->Fill(pythiaHard.event[i].pT());

                            if (bottomStatus < 40) BHardPtPart->Fill(pythiaHard.event[i].pT());
                            else if (bottomStatus < 60) BShowerPtPart->Fill(pythiaHard.event[i].pT());
                        }
                    }
                } 

                // final state muons
                if ((abs(particleID) == 13) && (pythiaHard.event[i].isFinal())) {
                    muPtPart->Fill(pythiaHard.event[i].pT());
                }

                // if (((abs(particleID) >= 441) && (abs(particleID) <= 445)) ||
                //         ((abs(particleID) >= 551) && (abs(particleID) <= 557))) { // Quarkonia
                //     // int muonIndex = muonDecay(i, pythiaHard.event);
                //     // if (muonIndex != -1) {
                //     //     cout << particleStatus << endl;   
                //     // }
                //     decayStatus(i, pythiaHard.event);
                //     cout << endl;
                // }
            }
        }

        // cross-section for the bin
        float luminocity_hard = N_events/(pythiaHard.info.sigmaGen()*pow(10,9));
        hardPtPart->Scale(1/luminocity_hard, "width");
        charmPtPart->Scale(1/luminocity_hard, "width");
        bottomPtPart->Scale(1/luminocity_hard, "width");
        DPtPart->Scale(1/luminocity_hard, "width");
        DHardPtPart->Scale(1/luminocity_hard, "width");
        DShowerPtPart->Scale(1/luminocity_hard, "width");
        BPtPart->Scale(1/luminocity_hard, "width");
        BHardPtPart->Scale(1/luminocity_hard, "width");
        BShowerPtPart->Scale(1/luminocity_hard, "width");
        muPtPart->Scale(1/luminocity_hard, "width");
        muHFPtPart->Scale(1/luminocity_hard, "width");
        muDPtPart->Scale(1/luminocity_hard, "width");
        muBPtPart->Scale(1/luminocity_hard, "width");
        //muOniumPtPart->Scale(1/luminocity_hard, "width");

        // add to final distribution
        hardPt->Add(hardPtPart);
        charmPt->Add(charmPtPart);
        bottomPt->Add(bottomPtPart);
        DPt->Add(DPtPart);
        DHardPt->Add(DHardPtPart);
        DShowerPt->Add(DShowerPtPart);
        BPt->Add(BPtPart);
        BHardPt->Add(BHardPtPart);
        BShowerPt->Add(BShowerPtPart);
        muPt->Add(muPtPart);
        muHFPt->Add(muHFPtPart);
        muDPt->Add(muDPtPart);
        muBPt->Add(muBPtPart);
        //muOniumPt->Add(muOniumPtPart);
    }

    //Plotting
    TCanvas *c1 = new TCanvas("total_sigma","total_sigma");

    hardPt->SetLineColor(1);
    hardPt->Draw();

    c1->Write();

    TCanvas *c2 = new TCanvas("mu_sigma","mu_sigma");

    muPt->SetLineColor(1);
    muPt->Draw();

    muHFPt->SetLineColor(2);
    muHFPt->Draw("SAME");

    muDPt->SetLineColor(3);
    muDPt->Draw("SAME");

    muBPt->SetLineColor(4);
    muBPt->Draw("SAME");

    // muOniumPt->SetLineColor(6);
    // muOniumPt->Draw("SAME");

    auto legend2 = new TLegend();
    legend2->AddEntry(muPt,"Final State #mu","l");
    legend2->AddEntry(muHFPt,"c,b->X->#mu","l");
    legend2->AddEntry(muDPt,"c->D->#mu","l");
    legend2->AddEntry(muBPt,"b->B->#mu","l");
    //legend2->AddEntry(muOniumPt,"Quarkonium->#mu","l");
    legend2->Draw("SAME");

    c2->Write();

    TCanvas *c3 = new TCanvas("HF_sigma","HF_sigma");

    charmPt->SetLineColor(1);
    charmPt->SetStats(0);
    charmPt->Draw();

    bottomPt->SetLineColor(2);
    bottomPt->SetStats(0);
    bottomPt->Draw("SAME");

    auto legend3 = new TLegend();
    legend3->AddEntry(charmPt,"Charm Total","l");
    legend3->AddEntry(bottomPt,"Bottom Total","l");
    legend3->Draw("SAME");

    c3->Write();

    TCanvas *c4 = new TCanvas("D_sigma","D_sigma");

    DPt->SetLineColor(1);
    DPt->SetStats(0);
    DPt->Draw();

    DHardPt->SetLineColor(2);
    DHardPt->SetStats(0);
    DHardPt->Draw("SAME");

    DShowerPt->SetLineColor(3);
    DShowerPt->SetStats(0);
    DShowerPt->Draw("SAME");

    auto legend4 = new TLegend();
    legend4->AddEntry(DPt,"D Total","l");
    legend4->AddEntry(DHardPt,"c(Hard+MPI)->D","l");
    legend4->AddEntry(DShowerPt,"c(Showers)->D","l");
    legend4->Draw("SAME");

    c4->Write();

    TCanvas *c5 = new TCanvas("B_sigma","B_sigma");

    BPt->SetLineColor(1);
    BPt->SetStats(0);
    BPt->Draw();

    BHardPt->SetLineColor(2);
    BHardPt->SetStats(0);
    BHardPt->Draw("SAME");

    BShowerPt->SetLineColor(3);
    BShowerPt->SetStats(0);
    BShowerPt->Draw("SAME");

    auto legend5 = new TLegend();
    legend5->AddEntry(BPt,"B Total","l");
    legend5->AddEntry(BHardPt,"b(Hard+MPI)->B","l");
    legend5->AddEntry(BShowerPt,"b(Showers)->B","l");
    legend5->Draw("SAME");

    c5->Write();

    delete outFile;

    return 0;
}