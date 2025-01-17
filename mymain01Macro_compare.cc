#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain01Macro_compare() {
    // Access Multiplicity distributions
    TFile* infile_central = TFile::Open("mymain01Macro_central.root", "READ");
    TFile* infile_central_CR = TFile::Open("mymain01Macro_central_CR_off.root", "READ");
    TFile* infile_central_MPI = TFile::Open("mymain01Macro_central_MPI_off.root", "READ");
    TFile* infile_forward = TFile::Open("mymain01Macro_forward.root", "READ");
    TFile* infile_V0C = TFile::Open("mymain01Macro_V0C.root", "READ");

    // raw charged particle multilicity 
    TH1D* mult_raw_central = (TH1D*) infile_central->Get("multiplicity_events_raw");
    TH1D* mult_raw_central_CR = (TH1D*) infile_central_CR->Get("multiplicity_events_raw");
    TH1D* mult_raw_central_MPI = (TH1D*) infile_central_MPI->Get("multiplicity_events_raw");
    TH1D* mult_raw_forward = (TH1D*) infile_forward->Get("multiplicity_events_raw");
    TH1D* mult_raw_V0C = (TH1D*) infile_V0C->Get("multiplicity_events_raw");

    // output the <N_ch> for the different systems
    std::cout << "Central Full Monash: <N_ch> = " << mult_raw_central->GetMean() << "; sigma = " << mult_raw_central->GetStdDev() << std::endl;
    std::cout << "Central CR off: <N_ch> = " << mult_raw_central_CR->GetMean() << "; sigma = " << mult_raw_central_CR->GetStdDev() << std::endl;
    std::cout << "Central MPI off: <N_ch> = " << mult_raw_central_MPI->GetMean() << "; sigma = " << mult_raw_central_MPI->GetStdDev() << std::endl;
    std::cout << "Forward: <N_ch> = " << mult_raw_forward->GetMean() << std::endl;
    std::cout << "V0C: <N_ch> = " << mult_raw_V0C->GetMean() << std::endl;

    // Output file
    TFile* outFile = new TFile("mymain01Macro_compare.root", "RECREATE");

    // Plot with varying model parameters (MPI+CR)
    TCanvas *canvasMultModels = new TCanvas("mult_models","mult_models");
    gPad->SetLogy();

    mult_raw_central_CR->SetMaximum(10);
    mult_raw_central_CR->SetStats(0);
    mult_raw_central_CR->SetLineWidth(3);
    mult_raw_central_CR->SetLineColor(2);
    mult_raw_central_CR->SetMarkerStyle(21);
    mult_raw_central_CR->SetMarkerColor(2);
    mult_raw_central_CR->SetMarkerSize(1.3);
    mult_raw_central_CR->Draw("SAME");

    mult_raw_central->SetStats(0);
    mult_raw_central->SetLineWidth(3);
    mult_raw_central->SetLineColor(3);
    mult_raw_central->SetMarkerStyle(22);
    mult_raw_central->SetMarkerColor(3);
    mult_raw_central->SetMarkerSize(1.3);
    mult_raw_central->Draw("SAME");

    mult_raw_central_MPI->SetStats(0);
    mult_raw_central_MPI->SetLineWidth(3);
    mult_raw_central_MPI->SetLineColor(4);
    mult_raw_central_MPI->SetMarkerStyle(20);
    mult_raw_central_MPI->SetMarkerColor(4);
    mult_raw_central_MPI->SetMarkerSize(1.3);
    mult_raw_central_MPI->Draw("SAME");

    auto legendModels = new TLegend();
    legendModels->AddEntry(mult_raw_central, "Full Monash","p");
    legendModels->AddEntry(mult_raw_central_CR, "CR off","p");
    legendModels->AddEntry(mult_raw_central_MPI, "MPI off","p");
    legendModels->Draw("SAME");

    auto labelModels = new TLatex();
    labelModels->DrawLatex(0.0, 0.0, "Pythia8 pp @ #sqrt{s} = 5.36 TeV, Minimum Bias (QCD)");
    labelModels->DrawLatex(0.0, 0.0, "N_{ch} > 0, |#eta_{ch}| < 1");
    labelModels->Draw("SAME");

    canvasMultModels->Write();

    // plot with different rapidity regions
    TCanvas *canvasMultRegions = new TCanvas("mult_regions","mult_regions");
    gPad->SetLogy();

    mult_raw_central->SetMinimum(0.00001);
    mult_raw_central->SetStats(0);
    mult_raw_central->SetLineWidth(3);
    mult_raw_central->SetLineColor(3);
    mult_raw_central->SetMarkerStyle(22);
    mult_raw_central->SetMarkerColor(3);
    mult_raw_central->SetMarkerSize(1.3);
    mult_raw_central->Draw("SAME");

    mult_raw_V0C->SetStats(0);
    mult_raw_V0C->SetLineWidth(3);
    mult_raw_V0C->SetLineColor(2);
    mult_raw_V0C->SetMarkerStyle(21);
    mult_raw_V0C->SetMarkerColor(2);
    mult_raw_V0C->SetMarkerSize(1.3);
    mult_raw_V0C->Draw("SAME");

    mult_raw_forward->SetStats(0);
    mult_raw_forward->SetLineWidth(3);
    mult_raw_forward->SetLineColor(4);
    mult_raw_forward->SetMarkerStyle(20);
    mult_raw_forward->SetMarkerColor(4);
    mult_raw_forward->SetMarkerSize(1.3);
    mult_raw_forward->Draw("SAME");

    auto legendRegions = new TLegend();
    legendRegions->AddEntry(mult_raw_central, "Central (|#eta_{ch}| < 1)","p");
    legendRegions->AddEntry(mult_raw_V0C, "V0C (1.7 < #eta_{ch} < 3.7)","p");
    legendRegions->AddEntry(mult_raw_forward, "Muon Arm (2.5 < #eta_{ch} < 4.0)","p");
    legendRegions->Draw("SAME");

    auto labelRegions = new TLatex();
    labelRegions->DrawLatex(0.0, 0.0, "Pythia8 pp @ #sqrt{s} = 5.36 TeV, Monash Tune");
    labelRegions->DrawLatex(0.0, 0.0, "Minimum Bias (QCD), N_{ch} > 0");
    labelRegions->Draw("SAME");

    canvasMultRegions->Write();

    delete outFile;
}