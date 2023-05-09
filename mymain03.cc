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
        cout << decayLevel << " (" << eventObj[particleIndex].name() << ")" << ":";  
        for (int momIndex: eventObj[particleIndex].motherList()) {
            cout << ' ' << eventObj[momIndex].name() << " (" << eventObj[momIndex].id() << ")";
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

// Checks if passed particle index has passed id as a mother
bool motherSearch(int particleIndex, int motherID, Pythia8::Event eventObj) {
    if ((particleIndex == 0) || (particleIndex == 1) || (particleIndex == 2)) return false;

    if (abs(eventObj[particleIndex].id()) == motherID) return true;

    bool hasMother = false;
    for (int momIndex: eventObj[particleIndex].motherList()) {
        if (motherSearch(momIndex, motherID, eventObj)) {
            hasMother = true;
        }
    }

    return hasMother;
}

int main() {
    Pythia pythia;

    // Collision settings/parameters
    pythia.readString("HardQCD:all = on");

    // Collision energy
    pythia.readString("Beams:eCM = 5020.");

    // Monash 13
    pythia.readString("Tune:pp = 14");

    // initialize
    pythia.init();

    // ROOT file for histograms
    TFile* outFile = new TFile("mymain03.root", "RECREATE");

    // TNtuple appraoch for cuts and histograms
    TNtuple* muInclusive = new TNtuple("mu_yield_inc", "mu", "pt:y:eta");
    TNtuple* muCharm = new TNtuple("mu_yield_c", "mu", "pt:y:eta");
    TNtuple* muCharmMeson = new TNtuple("mu_yield_d_meson", "mu", "pt:y:eta");
    TNtuple* muBottom = new TNtuple("mu_yield_b", "mu", "pt:y:eta");
    TNtuple* muOnium = new TNtuple("mu_yield_onium", "mu", "pt:y:eta");
    TNtuple* muLightMeson = new TNtuple("mu_yield_light_meson", "mu", "pt:y:eta");
    //TNtuple* muTau = new TNtuple("mu_yield_tau", "mu", "pt:y:eta");

    // Generate events
    int N_events = 300000;
    int muCountInc = 0;
    int muCountC = 0;
    int muCountB = 0;
    int muCountC_Meson = 0;
    int muCountB_Meson = 0;
    int muCountC_Baryon = 0;
    int muCountB_Baryon = 0;
    int muCountOnium = 0;
    int muCountLight = 0;
    int muCountS_Baryon = 0;
    int muCountTau = 0;
    int muBoson = 0;

    for (int iEvent = 0; iEvent < N_events; ++iEvent) {
        if (!pythia.next()) continue;

        int muCountEvent = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if ((pythia.event[i].id() == 13) || (pythia.event[i].id() == -13) && (pythia.event[i].isFinal())) { // is final state muon
                muCountEvent++;
                muCountInc++;
                muInclusive->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());

                // Determine if some c decay
                bool c_parent = motherSearch(i, 4, pythia.event);
                if (c_parent) {
                    muCountC++;
                    muCharm->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    // decayStatus(i, pythia.event, 0);
                    // cout << endl;
                }

                // Determine if some b decay
                bool b_parent = motherSearch(i, 5, pythia.event);
                if (b_parent) {
                    muCountB++;
                    muBottom->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    // decayStatus(i, pythia.event, 0);
                    // cout << endl;
                }

                // Determine decay channel
                int motherID = abs(pythia.event[pythia.event[i].mother1()].id());

                // Charmed Mesons (D mesons)
                if ((motherID >= 411) && (motherID <= 435)) {
                    if (c_parent) { // must also come from charm
                        muCountC_Meson++;
                        muCharmMeson->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    }
                    //decayStatus(i, pythia.event, 0);
                }

                // Charmed Baryons
                else if ((motherID >= 4112) && (motherID <= 4444)) {
                    muCountC_Baryon++;
                    //decayStatus(i, pythia.event, 0);
                }

                // C-Cbar quarkonium
                 else if ((motherID >= 441) && (motherID <= 445)) {
                    muCountOnium++;
                    muOnium->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    //decayStatus(i, pythia.event, 0);
                }

                // Bottom mesons
                else if ((motherID >= 511) && (motherID <= 545)) {
                    muCountB_Meson++;
                    //decayStatus(i, pythia.event, 0);
                }

                // Bottom Baryons
                else if ((motherID >= 5112) && (motherID <= 5554)) {
                    muCountB_Baryon++;
                    //decayStatus(i, pythia.event, 0);
                }

                // B-Bbar quarkonium
                else if ((motherID >= 551) && (motherID <= 557)) {
                    muCountOnium++;
                    muOnium->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    //decayStatus(i, pythia.event, 0);
                }

                // Light mesons: omega, eta, pi, K
                else if ((motherID >= 111) && (motherID <= 337)) {
                    muCountLight++;
                    muLightMeson->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    //decayStatus(i, pythia.event, 0);
                }

                // From tau (W?)
                else if (abs(motherID)==15) {
                    muCountTau++;
                    // muTau->Fill(pythia.event[i].pT(), pythia.event[i].y(), pythia.event[i].eta());
                    // decayStatus(i, pythia.event, 0);
                    // cout << endl;
                    // Try invariant mass for W
                    int index1 = -1;
                    int index2 = -1;
                    for (int daughterIndex: pythia.event[pythia.event[i].mother1()].daughterList()) {
                        if (abs(pythia.event[daughterIndex].id()) != 16) {
                            if (index1 == -1) {
                                index1 = daughterIndex;
                            } else {
                                index2 = daughterIndex;
                            }
                        }
                    }
                    double invariantMass = m(pythia.event[index1], pythia.event[index2]);
                    cout << "Invariant Mass: " << invariantMass << endl;
                }

                // Strange baryon
                else if ((motherID >= 3112) && (motherID <= 3334)) {
                    muCountS_Baryon++;
                    //decayStatus(i, pythia.event, 0);
                }

                // Other
                // else {
                // //    decayStatus(i, pythia.event, 0);
                // //    cout << endl;
                // //   pythia.event.list();
                // }
            }
        }

        // if (muCountEvent > 0) {
        //     cout << "Number of final state muons in event: " << muCountEvent << endl;
        //     pythia.event.list();
        //     cout << endl;
        // }
    }

    cout << "Total: " << muCountInc << endl;
    cout << "Charm: " << muCountC << endl;
    cout << "Bottom: " << muCountB << endl;
    cout << "Quarkonium: " << muCountOnium << endl;
    cout << "Tau: " << muCountTau << endl;
    cout << "Other: " << (muCountInc - muCountC - muCountB - muCountOnium - muCountTau) << endl;

    TCanvas *c1 = new TCanvas("c1","c1");

    const Int_t NBINS = 12;
    Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20};
    //const Int_t NBINS = 17;
    //Double_t edges[NBINS + 1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40};

    // Fill and draw histograms
    TH1F *muIncPt = new TH1F("mu_inc_pt","Inlusive muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muIncPt->SetLineColor(1);
    //muIncPt->SetStats(0);
    muInclusive->Draw("pt>>mu_inc_pt", "pt<20");

    TH1F *muCPt = new TH1F("mu_c_pt","Charm -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muCPt->SetLineColor(2);
    muCPt->SetStats(0);
    muCharm->Draw("pt>>mu_c_pt", "pt<20", "SAME");

    TH1F *muCMesonPt = new TH1F("mu_c_meson_pt","Charm -> D -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muCMesonPt->SetLineColor(46);
    muCMesonPt->SetStats(0);
    muCharmMeson->Draw("pt>>mu_c_meson_pt", "pt<20", "SAME");

    TH1F *muBPt = new TH1F("mu_b_pt","Bottom -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muBPt->SetLineColor(4);
    muBPt->SetStats(0);
    muBottom->Draw("pt>>mu_b_pt", "pt<20", "SAME");

    TH1F *muOniumPt = new TH1F("mu_onium_pt","Quarkonium -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muOniumPt->SetLineColor(6);
    muOniumPt->SetStats(0);
    muOnium->Draw("pt>>mu_onium_pt", "pt<20", "SAME");

    TH1F *muLightPt = new TH1F("mu_light_pt","Light Meson -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    muLightPt->SetLineColor(14);
    muLightPt->SetStats(0);
    muLightMeson->Draw("pt>>mu_light_pt", "pt<20", "SAME");

    // TH1F *muTauPt = new TH1F("mu_tau_pt","Tau -> muon decay pt;Momentum (GeV/c);dN/dpT", NBINS, edges);
    // muTauPt->SetLineColor(3);
    // muTauPt->SetStats(0);
    // muTau->Draw("pt>>mu_tau_pt", "pt<40", "SAME");

    auto legend = new TLegend();
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend->AddEntry(muIncPt,"Inclusive #mu yield","l");
    legend->AddEntry(muCPt,"#mu <- c","l");
    legend->AddEntry(muCMesonPt,"#mu <- D <- c","l");
    legend->AddEntry(muBPt,"#mu <- b","l");
    legend->AddEntry(muOniumPt,"#mu <- Quarkonium","l");
    legend->AddEntry(muLightPt,"#mu <- Light Meson","l");
    //legend->AddEntry(muTauPt,"#mu <- #tau","l");
    legend->Draw();

    c1->Write();

    delete outFile;

    return 0;
}