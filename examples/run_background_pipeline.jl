# Run from the repository root:
# julia --project=. examples/run_background_pipeline.jl [input_dir] [output_dir] [file1.csv file2.csv ...]
# With no arguments, use the bundled synthetic two-sample data.
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")
const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPOSITORY_ROOT, "MyFunctions.jl"))
include(joinpath(REPOSITORY_ROOT, "MyProcess.jl"))
using JLD2

function main(args)
    input_dir = length(args) >= 1 ? abspath(args[1]) : joinpath(REPOSITORY_ROOT, "sample data")
    output_dir = length(args) >= 2 ? abspath(args[2]) : joinpath(REPOSITORY_ROOT, "output")
    file_names = length(args) >= 3 ? String.(args[3:end]) :
                 (isempty(args) ? ["Sample_A.csv", "Sample_B.csv"] : String[])
    mkpath(output_dir)

    # 1. Loading; the keyword API validates aligned coordinates and keeps full file stems.
    combineData = process_genomic_data(input_dir; file_names=file_names, file_extension=".csv")
    ncol(combineData) >= 5 || error("This pairwise example requires at least two samples.")
    # 2. Filtering and normalization. These are example settings, not universal cutoffs.
    metaDataDF = filtering_normalizing_seqDepth(combineData;
        initial_quantiles=(0.01, 0.99), normalized_limits=[10, 250])
    CSV.write(joinpath(output_dir, "metaDataDF.csv"), metaDataDF)
    sample_names = String.(names(metaDataDF)[4:end])
    CSV.write(joinpath(output_dir, "sample_index.csv"),
        DataFrame(sample_index=collect(eachindex(sample_names)), sample_name=sample_names))

    # 3. Integer-CN constrained GMM with a full-support broad Gaussian background.
    cn_param = fitting_copy_number_gmm(metaDataDF;
        window_size=50, gmm_window_size=10,
        min_peak_prominence=0.05,
        modal_peak_copy_number=2, copy_number_range=1:10,
        share_modal_peak_position=true, share_copy_number_spacing=true,
        shared_modal_peak_max_distance_fraction=0.50,
        shared_spacing_tolerance_fraction=0.20,
        include_background=true, max_background_weight=0.001,
        background_min_std_fraction=2.0, background_max_std_fraction=5.0,
        shared_std=true)
    JLD2.jldsave(joinpath(output_dir, "cn_param.jld2"); cn_param)

    # 4. Background-aware segmentation; extract the scaffold from the returned NamedTuple.
    segment_result = genome_segment_background(metaDataDF, cn_param, 0.001, 10, 10;
        min_background_probability=0.90)
    CNVscaffold_background = segment_result.CNVscaffold
    JLD2.jldsave(joinpath(output_dir, "segment_result.jld2"); segment_result)
    CSV.write(joinpath(output_dir, "CNVscaffold_background.csv"),
        DataFrame(start_row=CNVscaffold_background[:, 1],
                  end_row=CNVscaffold_background[:, 2],
                  row_span=CNVscaffold_background[:, 3]))
    isempty(CNVscaffold_background) && error(
        "No examination intervals were retained. Inspect data, detected peaks, and segmentation settings.")

    # 5. Copy-number calls; missing values are uncallable, not CN0.
    copynumber_result = copynumber_determination_background(
        metaDataDF, CNVscaffold_background, cn_param, 10;
        thresholdBFactor=2.0, max_background_probability=0.90)
    JLD2.jldsave(joinpath(output_dir, "copynumber_result.jld2"); copynumber_result)
    CSV.write(joinpath(output_dir, "copynumber_background.csv"), copynumber_result.copy_number)

    # 6. Pairwise comparison and optional adjacent-fragment RD extension.
    comparison_result = copyNumber_comparison_background(
        metaDataDF, CNVscaffold_background, copynumber_result, cn_param, 10;
        pvalue_theSame=0.005)
    JLD2.jldsave(joinpath(output_dir, "comparison_result.jld2"); comparison_result)
    filter_result = copyNumber_comparison_filter_fn_background(
        metaDataDF, CNVscaffold_background, copynumber_result, comparison_result;
        highRDthreshold=18, iterative=true, max_iterations=10)
    JLD2.jldsave(joinpath(output_dir, "filter_result.jld2"); filter_result)

    # 7. Plot all sample pairs; layout must have enough panels for the selected pairs.
    sample_count = length(sample_names)
    sample_pairs = reduce(vcat, [reshape([a, b], 1, 2)
        for a in 1:(sample_count - 1) for b in (a + 1):sample_count])
    plot_result = plotCNVregionMulti_RD_background(
        metaDataDF, CNVscaffold_background, filter_result, sample_pairs;
        copynumber_result=copynumber_result, layout=(size(sample_pairs, 1), 1),
        plotsize=(1200, 450 * size(sample_pairs, 1)),
        positive_cnv_content=0.80, chr="all", filtered_depth_diff_threshold=18)
    CSV.write(joinpath(output_dir, "CNVresult_raw.csv"), plot_result.cnv_regions)
    CSV.write(joinpath(output_dir, "CNVresult_filter.csv"), plot_result.filtered_cnv_regions)
    CSV.write(joinpath(output_dir, "background_regions.csv"), plot_result.background_regions)
    savefig(plot_result.plot, joinpath(output_dir, "CNVplot.png"))
    println("Finished. Output directory: ", output_dir)
    println("Retained bins: ", nrow(metaDataDF), "; intervals: ", size(CNVscaffold_background, 1))
    println("CNV regions: ", nrow(plot_result.cnv_regions),
            "; RD-filtered regions: ", nrow(plot_result.filtered_cnv_regions))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
