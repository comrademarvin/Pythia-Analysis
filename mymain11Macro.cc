#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <TNtuple.h>

void mymain11Macro() {
    TFile *infile_W_plus = TFile::Open("mymain11_W+_100k.root", "READ");
    TFile *infile_W_minus = TFile::Open("mymain11_W-_100k.root", "READ");

    // Ratio plot of pT cross-section
    TH1F* W_plus_pt = (TH1F*)infile_W_plus->Get("W_pt_sigma");
    TH1F* W_minus_pt = (TH1F*)infile_W_minus->Get("W_pt_sigma");
}