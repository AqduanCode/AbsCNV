#include <iostream>
#include <fstream>
#include <vector>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TString.h>

void export_autosomal_read_depth(const char *root_path, const char *csv_path) {
    TFile *file = TFile::Open(root_path, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: cannot open ROOT file: " << root_path << std::endl;
        return;
    }

    TDirectory *dir = dynamic_cast<TDirectory*>(file->Get("bin_1000"));
    if (!dir) {
        std::cerr << "ERROR: directory bin_1000 not found in " << root_path << std::endl;
        file->Close();
        return;
    }

    std::ofstream out(csv_path);
    if (!out) {
        std::cerr << "ERROR: cannot create CSV: " << csv_path << std::endl;
        file->Close();
        return;
    }
    out << "BinLowEdge,BinUpEdge,Content,Error,Chromosome\n";

    std::vector<TString> chromosomes;
    for (int i = 1; i <= 22; ++i) chromosomes.push_back(Form("chr%d", i));

    int found = 0;
    for (const auto &chrom : chromosomes) {
        TString hist_name = "his_rd_p_" + chrom + "_1000";
        TH1D *hist = dynamic_cast<TH1D*>(dir->Get(hist_name));
        if (!hist) {
            std::cerr << "ERROR: histogram not found: " << hist_name << std::endl;
            continue;
        }
        ++found;
        const int n_bins = hist->GetNbinsX();
        for (int i = 1; i <= n_bins; ++i) {
            out << hist->GetXaxis()->GetBinLowEdge(i) << ','
                << hist->GetXaxis()->GetBinUpEdge(i) << ','
                << hist->GetBinContent(i) << ','
                << hist->GetBinError(i) << ','
                << chrom << '\n';
        }
    }

    out.close();
    file->Close();
    if (found != 22) {
        std::cerr << "ERROR: expected 22 chromosome histograms, found " << found << std::endl;
    } else {
        std::cout << "Exported all 22 chromosome histograms to " << csv_path << std::endl;
    }
}
