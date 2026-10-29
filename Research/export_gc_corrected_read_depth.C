// Parameterized I/O only; original his_rd_p_chr*_1000_GC selection is unchanged.
#include <iostream>
#include <vector>
#include <fstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TH1D.h>
#include <TKey.h>
#include <TString.h>

void export_gc_corrected_read_depth(const char *root_path, const char *csv_path) {
    // 打开ROOT文件
    TFile *file = TFile::Open(root_path);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // 获取包含直方图的目录
    TDirectory *dir = (TDirectory*)file->Get("bin_1000");
    if (!dir) {
        std::cerr << "Directory 'bin_1000' not found!" << std::endl;
        file->Close();
        return;
    }

    // 打开CSV文件准备写入数据
    std::ofstream outputFile(csv_path);
    outputFile << "BinLowEdge,BinUpEdge,Content,Error,Chromosome\n"; // 写入表头

    // 定义要匹配的染色体名称
    std::vector<TString> chromosomes = {
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19",
        "chr20", "chr21", "chr22", "chrX", "chrY"
    };

    // 遍历目录中的所有键
    TIter iter(dir->GetListOfKeys());
    TKey *key;
    
    while ((key = (TKey*)iter())) {
        TString name = key->GetName();

        // 检查直方图名称是否匹配所需的模式
        for (const auto& chr : chromosomes) {
            TString pattern = "his_rd_p_" + chr + "_1000_GC";
            if (name == pattern) {
                // 获取直方图
                TH1D *hist = (TH1D*)dir->Get(name);
                if (hist) {
                    // 遍历直方图的每个bin，将数据写入CSV文件
                    int nBins = hist->GetNbinsX();
                    for (int i = 1; i <= nBins; ++i) {
                        // 获取每个bin的下限和上限值
                        double binLowEdge = hist->GetXaxis()->GetBinLowEdge(i);
                        double binUpEdge = hist->GetXaxis()->GetBinUpEdge(i);
                        double readDepth = hist->GetBinContent(i);
                        double binError = hist->GetBinError(i);
                        outputFile << binLowEdge << "," << binUpEdge << "," << readDepth << "," << binError << "," << chr << "\n";
                    }
                }
                break;  // 如果匹配成功，跳出循环
            }
        }
    }

    // 关闭文件
    outputFile.close();
    file->Close();
}