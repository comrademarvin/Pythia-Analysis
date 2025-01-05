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

    // raw charged particle multilicity 
    TH1D* mult_raw_central = (TH1D*) infile_central->Get("multiplicity_events_raw");
    TH1D* mult_raw_central_CR = (TH1D*) infile_central_CR->Get("multiplicity_events_raw");
    TH1D* mult_raw_central_MPI = (TH1D*) infile_central_MPI->Get("multiplicity_events_raw");

    // Output file
    TFile* outFile = new TFile("mymain01Macro_compare.root", "RECREATE");

    // Plot with varying model parameters (MPI+CR)
    TCanvas *canvasMultModels = new TCanvas("mult_models","mult_models");
    gPad->SetLogy();

    mult_raw_central_CR->SetStats(0);
    mult_raw_central_CR->SetLineWidth(3);
    mult_raw_central_CR->SetLineColor(2);
    mult_raw_central_CR->SetMarkerStyle(21);
    mult_raw_central_CR->SetMarkerColor(2);
    mult_raw_central_CR->Draw("SAME");

    mult_raw_central->SetStats(0);
    mult_raw_central->SetLineWidth(3);
    mult_raw_central->SetLineColor(3);
    mult_raw_central->SetMarkerStyle(22);
    mult_raw_central->SetMarkerColor(3);
    mult_raw_central->Draw("SAME");

    mult_raw_central_MPI->SetStats(0);
    mult_raw_central_MPI->SetLineWidth(3);
    mult_raw_central_MPI->SetLineColor(4);
    mult_raw_central_MPI->SetMarkerStyle(20);
    mult_raw_central_MPI->SetMarkerColor(4);
    mult_raw_central_MPI->Draw("SAME");

    auto legendModels = new TLegend();
    legendModels->AddEntry(mult_raw_central, "Full Monash","p");
    legendModels->AddEntry(mult_raw_central_CR, "CR off","p");
    legendModels->AddEntry(mult_raw_central_MPI, "MPI off","p");
    legendModels->Draw("SAME");

    canvasMultModels->Write();

    delete outFile;
}