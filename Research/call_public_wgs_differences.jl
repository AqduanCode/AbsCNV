# Paths adapted for the curated archive; analysis parameters and calls unchanged.
length(ARGS) == 2 || error("Usage: julia --project=REPO SCRIPT INPUT_DIR NEW_OUTPUT_DIR")
const REPO_DIR = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_DIR, "MyFunctions.jl"))
include(joinpath(REPO_DIR, "MyProcess.jl"))

using CSV
using DataFrames
using JLD2
using DelimitedFiles
using StatsBase

input_dir = abspath(ARGS[1])
output_dir = abspath(ARGS[2])
mkpath(output_dir)
file_names = ["FASR.csv","MCF7.csv","TAMR.csv"]

#input_dir = "D:/桌面/GMM方法论文/背景更新后的结果/真实下载的数据/WGBS结果/input/"
#output_dir = "D:/桌面/GMM方法论文/背景更新后的结果/真实下载的数据/WGBS结果/output/"
#file_names = ["FASR_WGBS.csv","MCF7_WGBS.csv","TAMR_WGBS.csv"]

# ============================================================
# 1. 读取并过滤、标准化数据
# ============================================================

combineData = process_genomic_data(
    input_dir;
    file_names = file_names,
    file_extension = ".csv",
)

metaDataDF = filtering_normalizing_seqDepth(
    combineData;
    initial_quantiles=(0.01, 0.99),
    normalized_limits=[10, 250],
)

CSV.write(
    joinpath(output_dir, "metaDataDF.csv"),
    metaDataDF,
)


# ============================================================
# 2. 整数CN约束GMM＋完整宽高斯背景
# ============================================================


cn_param = fitting_copy_number_gmm(
    metaDataDF;
    window_size = 50,
    gmm_window_size = 10,
    modal_peak_copy_number = 3,
    copy_number_range = 1:10,
    share_modal_peak_position = true,
    share_copy_number_spacing = true,
    shared_modal_peak_max_distance_fraction = 0.50,
    shared_spacing_tolerance_fraction = 0.20,
    include_background = true,
    spacing_tolerance_fraction=0.30,

    max_background_weight = 0.001,
    background_min_std_fraction = 2.0,
    background_max_std_fraction = 5.0,

    shared_std = true,
)

@save joinpath(output_dir, "cn_param.jld2") cn_param

metaData = Matrix(metaDataDF[:,4:end])

moving_average_histogram(metaData[:,3], cn_param[3].means, 50)
# ============================================================
# 3. 背景感知的基因组分段
# ============================================================

segment_result = genome_segment_background(
    metaDataDF,
    cn_param,
    0.001,
    10,
    10,
    min_background_probability = 0.90,
)

CNVscaffold_background = segment_result.CNVscaffold

writedlm(
    joinpath(output_dir, "CNVscaffold_background.csv"),
    CNVscaffold_background,
    ',',
)

@save joinpath(output_dir, "segment_result.jld2") segment_result


# ============================================================
# 4. 背景感知的片段CN调用
# ============================================================

copynumber_result =
    copynumber_determination_background(
        metaDataDF,
        CNVscaffold_background,
        cn_param,
        10;
        thresholdBFactor = 2.0,
        max_background_probability = 0.90,
    )


@save joinpath(
    output_dir,
    "copynumber_result.jld2"
) copynumber_result

CSV.write(
    joinpath(output_dir, "copynumber_background.csv"),
    copynumber_result.copy_number,
)

# ============================================================
# 5. 样本间CN比较
# ============================================================

comparison_result =
    copyNumber_comparison_background(
        metaDataDF,
        CNVscaffold_background,
        copynumber_result,
        cn_param,
        10;
        pvalue_theSame = 0.0005,
    )

@save joinpath(
    output_dir,
    "comparison_result.jld2"
) comparison_result

# ============================================================
# 6. 相邻CNV的RD差扩展
# ============================================================



filter_result =
    copyNumber_comparison_filter_fn_background(
        metaDataDF,
        CNVscaffold_background,
        copynumber_result,
        comparison_result;
        highRDthreshold = 20,
        iterative = true,
        max_iterations = 3,
    )

@save joinpath(
    output_dir,
    "filter_result.jld2"
) filter_result


# ============================================================
# 7. 绘图
# ============================================================
sample_pairs = [1 2; 1 3; 2 3;]
plot_result =
    plotCNVregionMulti_RD_background(
        metaDataDF,
        CNVscaffold_background,
        filter_result,
        sample_pairs;
        layout = (1, 3),
        plotsize = (1200, 500),
        positive_cnv_content = 1,
        copynumber_result = copynumber_result,
        chr = "chr7",
        merge_regions = true,
        min_gap = 1_000_000,
        max_gap = 3_000_000,
        diff_threshold = 15,
        gap_diff_threshold = 10,
        filtered_depth_diff_threshold = 20,
    )


CSV.write(
    joinpath(output_dir, "CNVresult_raw.csv"),
    plot_result.cnv_regions,
)

CSV.write(
    joinpath(output_dir, "CNVresult_filter.csv"),
    plot_result.filtered_cnv_regions,
)

CSV.write(
    joinpath(output_dir, "background_regions.csv"),
    plot_result.background_regions,
)

