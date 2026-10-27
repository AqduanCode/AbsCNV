
using Pkg
using Random

required_packages = [
    "CSV", "DataFrames", "Plots",
    "StatsBase", "Distributions", "DelimitedFiles",
    "MultipleTesting"
    ]

for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))
        println("$pkg loaded")
    catch
        println("$pkg is not installed; installing...")
        Pkg.add(pkg)
        eval(Meta.parse("using $pkg"))
        println("$pkg installed and loaded")
    end
end

# Compatibility overloads retained for existing callers.
function process_genomic_data(temp_path::String, file_names::Vector{String})::DataFrame

    if isempty(temp_path)
        throw(ArgumentError("temp_path must be a non-empty string"))
    end

    if isempty(file_names)
        throw(ArgumentError("file_names must be a non-empty vector"))
    end

    dfs = Vector{DataFrame}(undef, length(file_names))
    for i in 1:length(file_names)
        full_path = joinpath(temp_path, file_names[i])
        dfs[i] = CSV.read(full_path, DataFrame)
    end

    chr_strings = [
        "chr1", "chr2", "chr3", "chr4", "chr5",
        "chr6", "chr7", "chr8", "chr9", "chr10",
        "chr11", "chr12", "chr13", "chr14", "chr15",
        "chr16", "chr17", "chr18", "chr19", "chr20",
        "chr21", "chr22"
    ]

    n_rows = size(dfs[1], 1)
    chr_flag = zeros(n_rows, 1)

    for i in 1:n_rows
        chr_flag[i] = (sum(dfs[1].Chromosome[i] .== chr_strings) > 0)
    end

    defined_indices = findall(chr_flag[:] .== 1)

    combineData = DataFrame(
        Chromosome = dfs[1][defined_indices, :Chromosome],
        BinLowEdge = dfs[1][defined_indices, :BinLowEdge],
        BinUpEdge = dfs[1][defined_indices, :BinUpEdge],
    )

    for i in 1:length(dfs)

        column_prefix = split(file_names[i], "_")[1]
        column_name = replace(column_prefix, "-" => "")
        column_symbol = Symbol(column_name)

        combineData[!, column_symbol] = dfs[i][defined_indices, :Content]
    end

    return combineData
end

function moving_average_histogram(depth::Vector, Xbar::Vector, windowSize::Int)

    if windowSize > length(depth)
        throw(ArgumentError("Window size cannot exceed vector length"))
    end

    if windowSize < 1
        throw(ArgumentError("Window size must be a positive integer"))
    end

    n = length(depth)
    moving_avg = Vector{Float64}(undef, n - windowSize + 1)

    for i in 1:(n - windowSize + 1)
        window = depth[i:(i + windowSize - 1)]
        moving_avg[i] = mean(window)
    end

    p = histogram(moving_avg,
                  alpha=0.7,
                  fillcolor=:lightblue)

    ylims_current = ylims()
    max_height = ylims_current[2]

    for x_val in Xbar

        plot!([x_val, x_val], [0, max_height],
              color=:red,
              linewidth=2,
              legend = false)
    end

    display(p)

end

# Merge and plot CNV regions.
function merge_cnv_region_groups(
    region_df::DataFrame,
    metaDataDF::DataFrame;
    merge_regions::Bool = true,
    min_gap::Real = 1_000_000,
    max_gap::Real = 3_000_000,
    diff_threshold::Real = 15,
    gap_diff_threshold::Real = 10,
)
    merge_regions || return copy(region_df)
    isempty(region_df) && return copy(region_df)
    min_gap < max_gap || throw(ArgumentError("min_gap must be less than max_gap"))

    required = [
        :pair_sampleA, :pair_sampleB, :pair_sampleA_name, :pair_sampleB_name,
        :chr, :start_pos, :end_pos, :length, :start_pos_metaData,
        :end_pos_metaData, :mean_depth_A, :mean_depth_B, :depth_diff,
    ]
    missing_columns = setdiff(required, Symbol.(names(region_df)))
    isempty(missing_columns) ||
        throw(ArgumentError("region_df is missing columns: $(join(string.(missing_columns), ", "))"))

    sample_count = ncol(metaDataDF) - 3
    chr_indices = Dict{String, Vector{Int}}()
    for (idx, chr_value) in enumerate(metaDataDF[!, 1])
        push!(get!(chr_indices, String(chr_value), Int[]), idx)
    end

    sorted_df = sort(region_df, [:pair_sampleA, :pair_sampleB, :chr, :start_pos])
    merged_df = similar(sorted_df, 0)

    for group in groupby(sorted_df, [:pair_sampleA, :pair_sampleB, :chr])
        sample_a = Int(group[1, :pair_sampleA])
        sample_b = Int(group[1, :pair_sampleB])
        (1 <= sample_a <= sample_count && 1 <= sample_b <= sample_count) ||
            throw(ArgumentError("Sample index is outside the metaDataDF range: ($sample_a, $sample_b)"))

        chr_name = String(group[1, :chr])
        meta_idx = get(chr_indices, chr_name, Int[])
        low_edges = Float64.(metaDataDF[meta_idx, 2])
        up_edges = Float64.(metaDataDF[meta_idx, 3])
        depth_a = Float64.(coalesce.(metaDataDF[meta_idx, 3 + sample_a], 0.0))
        depth_b = Float64.(coalesce.(metaDataDF[meta_idx, 3 + sample_b], 0.0))

        i = 1
        while i <= nrow(group)
            first_row = group[i, :]
            first_diff = Float64(coalesce(first_row.depth_diff, 0.0))

            if abs(first_diff) <= diff_threshold
                push!(merged_df, first_row)
                i += 1
                continue
            end

            cluster_end = Float64(first_row.end_pos)
            total_weight = Float64(first_row.length)
            sum_a = Float64(first_row.mean_depth_A) * total_weight
            sum_b = Float64(first_row.mean_depth_B) * total_weight
            last_i = i
            j = i + 1

            while j <= nrow(group)
                next_row = group[j, :]
                gap_start = cluster_end
                gap_end = Float64(next_row.start_pos)
                gap = gap_end - gap_start
                gap >= max_gap && break

                next_diff = Float64(coalesce(next_row.depth_diff, 0.0))

                same_direction =
                    (first_diff > 0 && next_diff > 0) ||
                    (first_diff < 0 && next_diff < 0)
                same_direction || break

                gap_supports_merge = true
                if gap > 0
                    overlap = min.(up_edges, gap_end) .- max.(low_edges, gap_start)
                    in_gap = (up_edges .> gap_start) .& (low_edges .< gap_end) .& (overlap .> 0)
                    if any(in_gap)
                        weights = overlap[in_gap]
                        mean_gap_a = sum(depth_a[in_gap] .* weights) / sum(weights)
                        mean_gap_b = sum(depth_b[in_gap] .* weights) / sum(weights)
                        gap_supports_merge = abs(mean_gap_a - mean_gap_b) > gap_diff_threshold
                    else
                        gap_supports_merge = false
                    end
                end

                if gap >= min_gap &&
                   (abs(next_diff) <= diff_threshold || !gap_supports_merge)
                    break
                end

                row_weight = Float64(next_row.length)
                total_weight += row_weight
                sum_a += Float64(next_row.mean_depth_A) * row_weight
                sum_b += Float64(next_row.mean_depth_B) * row_weight
                cluster_end = max(cluster_end, Float64(next_row.end_pos))
                last_i = j
                j += 1
            end

            if last_i == i
                push!(merged_df, first_row)
            else
                total_weight > 0 || throw(ArgumentError("The total length of regions to merge must be positive"))
                mean_a = sum_a / total_weight
                mean_b = sum_b / total_weight
                last_row = group[last_i, :]
                push!(merged_df, (
                    pair_sampleA = first_row.pair_sampleA,
                    pair_sampleB = first_row.pair_sampleB,
                    pair_sampleA_name = first_row.pair_sampleA_name,
                    pair_sampleB_name = first_row.pair_sampleB_name,
                    chr = first_row.chr,
                    start_pos = first_row.start_pos,
                    end_pos = last_row.end_pos,
                    length = last_row.end_pos - first_row.start_pos,
                    start_pos_metaData = first_row.start_pos_metaData,
                    end_pos_metaData = last_row.end_pos_metaData,
                    mean_depth_A = mean_a,
                    mean_depth_B = mean_b,
                    depth_diff = mean_b - mean_a,
                ))
            end
            i = last_i + 1
        end
    end

    return merged_df
end

function plotCNVregionMulti_RD(
    metaDataDF,
    CNVscaffold,
    pairedCNV,
    sample_pair;
    layout = (1, 1),
    xlims_pre = [0 0],
    ylims_pre = [0 0],
    plotsize = (1000, 500),
    positive_cnv_content = 0.8,
    chr = "all",
    merge_regions = false,
    merge_min_gap = 1_000_000,
    merge_max_gap = 3_000_000,
    merge_diff_threshold = 15,
    merge_gap_diff_threshold = 10,
    filtered_depth_diff_threshold = 20,
)
    try

        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be a DataFrame"))
        end

        if !isa(CNVscaffold, Matrix) || size(CNVscaffold, 2) != 3
            throw(ArgumentError("CNVscaffold must be an n-by-3 matrix"))
        end

        if !isa(pairedCNV, Array) || length(size(pairedCNV)) != 3
            throw(ArgumentError("pairedCNV must be a three-dimensional array"))
        end
        if !isa(sample_pair, Matrix) || size(sample_pair, 2) != 2
            throw(ArgumentError("sample_pair must be an n-by-2 matrix"))
        end

        num_pairs = size(sample_pair, 1)
        sample_names = String.(names(metaDataDF)[4:end])
        sample_count = length(sample_names)
        all((1 .<= sample_pair) .& (sample_pair .<= sample_count)) ||
            throw(ArgumentError("sample_pair contains an index outside 1:$sample_count "))

        plt = plot(layout = layout, size = plotsize)
        xlims_current = zeros(num_pairs, 2)
        ylims_current = zeros(num_pairs, 2)
        plotcnv_regionDF = DataFrame()
        for pair_i = 1:num_pairs

            cnv_regionsDF_samplepair = cnv_region(
                metaDataDF,
                CNVscaffold,
                pairedCNV[:, sample_pair[pair_i, 1], sample_pair[pair_i, 2]],
                positive_cnv_content
            )

            metaDataDF_samplepair = select(
                metaDataDF,
                1, 2, 3,
                3 + sample_pair[pair_i, 1],
                3 + sample_pair[pair_i, 2]
            )

            depths_metaData = Matrix(metaDataDF_samplepair[:, 4:end])

            num_fragments = size(CNVscaffold, 1)
            depths_mean = zeros(num_fragments, 2)
            fragments_chr = String[]
            fragments_chr_startposi = Int64[]
            fragments_chr_endposi = Int64[]

            for i = 1:num_fragments
                push!(fragments_chr, String(metaDataDF_samplepair[CNVscaffold[i, 1], 1]))
                push!(fragments_chr_startposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 1], 2]))
                push!(fragments_chr_endposi, Int64(metaDataDF_samplepair[CNVscaffold[i, 2], 3]))

                start_idx = max(1, CNVscaffold[i, 1])
                end_idx = min(size(depths_metaData, 1), CNVscaffold[i, 2])
                if start_idx <= end_idx
                    depths_mean[i, 1] = mean(depths_metaData[start_idx:end_idx, 1])
                    depths_mean[i, 2] = mean(depths_metaData[start_idx:end_idx, 2])
                else
                    depths_mean[i, 1] = 0
                    depths_mean[i, 2] = 0
                end
            end

            if !isempty(cnv_regionsDF_samplepair)
                mean_depth_A = Float64[]
                mean_depth_B = Float64[]
                depth_diff = Float64[]
                for j in 1:nrow(cnv_regionsDF_samplepair)
                    start_bin = cnv_regionsDF_samplepair[j, :start_pos_metaData]
                    end_bin = cnv_regionsDF_samplepair[j, :end_pos_metaData]
                    if start_bin <= end_bin && start_bin >= 1 && end_bin <= size(depths_metaData, 1)
                        mean_A = mean(depths_metaData[start_bin:end_bin, 1])
                        mean_B = mean(depths_metaData[start_bin:end_bin, 2])
                        push!(mean_depth_A, mean_A)
                        push!(mean_depth_B, mean_B)
                        push!(depth_diff, mean_B - mean_A)
                    else
                        push!(mean_depth_A, 0.0)
                        push!(mean_depth_B, 0.0)
                        push!(depth_diff, 0.0)
                    end
                end
                cnv_regionsDF_samplepair[!, :mean_depth_A] = mean_depth_A
                cnv_regionsDF_samplepair[!, :mean_depth_B] = mean_depth_B
                cnv_regionsDF_samplepair[!, :depth_diff] = depth_diff
            end

            sample1 = repeat([sample_pair[pair_i,1]], size(cnv_regionsDF_samplepair, 1))
            sample2 = repeat([sample_pair[pair_i,2]], size(cnv_regionsDF_samplepair, 1))
            sample1_name = repeat(
                [sample_names[sample_pair[pair_i, 1]]],
                size(cnv_regionsDF_samplepair, 1),
            )
            sample2_name = repeat(
                [sample_names[sample_pair[pair_i, 2]]],
                size(cnv_regionsDF_samplepair, 1),
            )
            temp_DF = DataFrame(
                pair_sampleA = sample1,
                pair_sampleB = sample2,
                pair_sampleA_name = sample1_name,
                pair_sampleB_name = sample2_name,
            )
            temp_DF = hcat(temp_DF, cnv_regionsDF_samplepair)
            plotcnv_regionDF = vcat(plotcnv_regionDF, temp_DF)

            if chr == "all"

                selected_idx = findall(fragments_chr_endposi .- fragments_chr_startposi .> 0)

                temp_x = 1:length(selected_idx)
                temp_y1 = depths_mean[selected_idx, 1]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y1,
                    linewidth = 0.5,
                    color = "red",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 1])"
                )

                temp_y2 = depths_mean[selected_idx, 2]
                plot!(
                    plt,
                    subplot = pair_i,
                    temp_x,
                    temp_y2,
                    linewidth = 0.5,
                    color = "blue",
                    legend = false,
                    label = "Sample $(sample_pair[pair_i, 2])"
                )

                if size(cnv_regionsDF_samplepair, 1) > 0

                end
            else

                selected_idx = findall((fragments_chr .== chr) .&& (fragments_chr_endposi .- fragments_chr_startposi .> 0))
                if !isempty(selected_idx)
                    temp_x = hcat(fragments_chr_startposi[selected_idx], fragments_chr_endposi[selected_idx])
                    temp_y1 = hcat(depths_mean[selected_idx, 1], depths_mean[selected_idx, 1])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y1',
                        linewidth = 1,
                        color = "red",
                        legend = false,
                        title = "Sample $(sample_pair[pair_i, :])"
                    )

                    temp_y2 = hcat(depths_mean[selected_idx, 2], depths_mean[selected_idx, 2])
                    plot!(
                        plt,
                        subplot = pair_i,
                        temp_x',
                        temp_y2',
                        linewidth = 1,
                        color = "blue",
                        legend = false,
                    )

                    if isa(cnv_regionsDF_samplepair, DataFrame) && size(cnv_regionsDF_samplepair, 1) > 0 && size(cnv_regionsDF_samplepair, 2) >= 3
                        chr_cnv_rows = cnv_regionsDF_samplepair[:, :chr] .== chr
                        if any(chr_cnv_rows)
                            selected_cnv_idx = findall(chr_cnv_rows)

                            valid_mask = abs.(cnv_regionsDF_samplepair[selected_cnv_idx, :depth_diff]) .>= 15
                            selected_cnv_idx = selected_cnv_idx[valid_mask]
                            if !isempty(selected_cnv_idx)
                                temp_xCNV = hcat(
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :start_pos],
                                    cnv_regionsDF_samplepair[selected_cnv_idx, :end_pos]
                                )
                                temp_yCNV = ones(length(selected_cnv_idx), 2) * 1
                                plot!(
                                    plt,
                                    subplot = pair_i,
                                    temp_xCNV',
                                    temp_yCNV',
                                    linewidth = 2,
                                    color = "black",
                                    legend = false
                                )
                            end
                        end
                    end
                else
                    @warn "No data found for chromosome $chr "
                end
            end
        end

        for pair_i = 1:num_pairs
            xlims_current[pair_i, :] = [xlims(plt[pair_i])[1] xlims(plt[pair_i])[2]]
            ylims_current[pair_i, :] = [ylims(plt[pair_i])[1] ylims(plt[pair_i])[2]]
        end
        min_xlims = minimum(xlims_current[:, 1])
        min_ylims = minimum(ylims_current[:, 1])
        max_ylims = maximum(ylims_current[:, 2])
        max_xlims = maximum(xlims_current[:, 2])

        if xlims_pre == [0 0]
            for pair_i = 1:num_pairs
                xlims!(plt[pair_i], min_xlims, max_xlims)
            end
        else
            for pair_i = 1:num_pairs
                xlims!(plt[pair_i], xlims_pre[1], xlims_pre[2])
            end
        end
        if ylims_pre == [0 0]
            for pair_i = 1:num_pairs
                ylims!(plt[pair_i], min_ylims, max_ylims)
            end
        else
            for pair_i = 1:num_pairs
                ylims!(plt[pair_i], ylims_pre[1], ylims_pre[2])
            end
        end

        display(plt)

       if hasproperty(plotcnv_regionDF, :depth_diff) && !isempty(plotcnv_regionDF)
          filter_source = merge_regions ? merge_cnv_region_groups(
            plotcnv_regionDF,
            metaDataDF;
            min_gap = merge_min_gap,
            max_gap = merge_max_gap,
            diff_threshold = merge_diff_threshold,
            gap_diff_threshold = merge_gap_diff_threshold,
          ) : plotcnv_regionDF
          valid_filter = abs.(coalesce.(filter_source.depth_diff, 0.0)) .>= filtered_depth_diff_threshold
          subplotcnv_regionDF = filter_source[valid_filter, :]
        else
            subplotcnv_regionDF = DataFrame()
        end

        return plotcnv_regionDF, subplotcnv_regionDF, plt

    catch e
        @error "An error occurred in plotCNVregionMulti: " * string(e)
        rethrow(e)
    end
end

# Background-aware plotting helpers.
function _cnv_plot_state_with_background_barriers(
    metaDataDF,
    CNVscaffold,
    state,
    positive_cnv_content,
)
    fragment_count, sample_count1, sample_count2 = size(state)
    fragment_count == size(CNVscaffold, 1) || throw(DimensionMismatch(
        "The first state-matrix dimension does not match the number of CNVscaffold fragments"))
    output = zeros(Int8, size(state))
    fragment_chr = String.(metaDataDF[CNVscaffold[:, 1], 1])
    fragment_lengths = CNVscaffold[:, 2] .- CNVscaffold[:, 1] .+ 1

    for sample1 in 1:sample_count1, sample2 in 1:sample_count2
        status = vec(state[:, sample1, sample2])
        block_start = 1
        while block_start <= fragment_count
            while block_start <= fragment_count && status[block_start] == -1
                block_start += 1
            end
            block_start > fragment_count && break
            block_end = block_start
            while block_end < fragment_count &&
                  status[block_end + 1] != -1 &&
                  fragment_chr[block_end + 1] == fragment_chr[block_start]
                block_end += 1
            end
            block_status = status[block_start:block_end]
            if any(block_status .== 1)
                block_indices = block_start:block_end
                fragment_coords = hcat(
                    fragment_chr[block_indices],
                    CNVscaffold[block_indices, 1],
                    CNVscaffold[block_indices, 2],
                    fragment_lengths[block_indices],
                )
                regions = synth_cnv(
                    fragment_coords,
                    block_status,
                    positive_cnv_content,
                )
                if !(regions == -1) && !isempty(regions)
                    for region in axes(regions, 1)
                        global_start = block_start + regions[region, 1] - 1
                        global_end = block_start + regions[region, 2] - 1
                        output[global_start:global_end, sample1, sample2] .= 1
                    end
                end
            end
            block_start = block_end + 1
        end
    end
    return output
end

function _cnv_background_regions(metaDataDF, CNVscaffold, state, sample_pair)
    output = DataFrame(
        pair_sampleA=Int[], pair_sampleB=Int[],
        pair_sampleA_name=String[], pair_sampleB_name=String[],
        chr=String[], start_pos=Int[], end_pos=Int[],
        start_fragment=Int[], end_fragment=Int[], fragment_count=Int[],
    )
    sample_names = String.(names(metaDataDF)[4:end])
    fragment_chr = String.(metaDataDF[CNVscaffold[:, 1], 1])
    for pair_index in axes(sample_pair, 1)
        sample1, sample2 = sample_pair[pair_index, :]
        status = vec(state[:, sample1, sample2])
        fragment = 1
        while fragment <= length(status)
            if status[fragment] != -1
                fragment += 1
                continue
            end
            region_start = fragment
            while fragment < length(status) && status[fragment + 1] == -1 &&
                  fragment_chr[fragment + 1] == fragment_chr[region_start]
                fragment += 1
            end
            region_end = fragment
            push!(output, (
                pair_sampleA=sample1,
                pair_sampleB=sample2,
                pair_sampleA_name=sample_names[sample1],
                pair_sampleB_name=sample_names[sample2],
                chr=fragment_chr[region_start],
                start_pos=Int(metaDataDF[CNVscaffold[region_start, 1], 2]),
                end_pos=Int(metaDataDF[CNVscaffold[region_end, 2], 3]),
                start_fragment=region_start,
                end_fragment=region_end,
                fragment_count=region_end - region_start + 1,
            ))
            fragment += 1
        end
    end
    return output
end

function _copy_number_call_matrix(copynumber_result)
    if hasproperty(copynumber_result, :copy_number_calls)
        calls = copynumber_result.copy_number_calls
        calls isa AbstractMatrix || throw(ArgumentError(
            "copynumber_result.copy_number_calls must be a matrix"))
        return calls
    end

    if hasproperty(copynumber_result, :copy_number) &&
       copynumber_result.copy_number isa AbstractMatrix
        return copynumber_result.copy_number
    end
    throw(ArgumentError(
        "copynumber_result does not contain the internal copy_number_calls matrix"))
end

function _cnv_region_weighted_different_copy_number_pair(
    copynumber_result,
    CNVscaffold,
    event_state,
    sample_A::Int,
    sample_B::Int,
    region_start::Int,
    region_end::Int,
)
    copy_number = _copy_number_call_matrix(copynumber_result)
    callable = copynumber_result.callable
    weights_by_pair = Dict{Tuple{Int, Int}, Int}()

    for fragment in axes(CNVscaffold, 1)
        fragment_start = CNVscaffold[fragment, 1]
        fragment_end = CNVscaffold[fragment, 2]
        fragment_end < region_start && continue
        fragment_start > region_end && continue
        event_state[fragment, sample_A, sample_B] == 1 || continue
        callable[sample_A, fragment] || continue
        callable[sample_B, fragment] || continue

        cn_A = copy_number[sample_A, fragment]
        cn_B = copy_number[sample_B, fragment]
        (ismissing(cn_A) || ismissing(cn_B)) && continue
        integer_cn_A = Int(cn_A)
        integer_cn_B = Int(cn_B)
        integer_cn_A != integer_cn_B || continue

        overlap_length = min(fragment_end, region_end) -
                         max(fragment_start, region_start) + 1
        overlap_length > 0 || continue
        pair = (integer_cn_A, integer_cn_B)
        weights_by_pair[pair] =
            get(weights_by_pair, pair, 0) + overlap_length
    end

    isempty(weights_by_pair) && return (missing, missing)

    return first(sort!(collect(keys(weights_by_pair));
                       by=pair -> (-weights_by_pair[pair], pair[1], pair[2])))
end

function _cnv_add_region_copy_numbers(
    regions::DataFrame,
    CNVscaffold,
    copynumber_result,
    event_state,
)
    output = copy(regions)
    copy_number_A = Union{Missing, Int}[missing for _ in 1:nrow(output)]
    copy_number_B = Union{Missing, Int}[missing for _ in 1:nrow(output)]

    if nrow(output) == 0
        output[!, :copy_number_A] = copy_number_A
        output[!, :copy_number_B] = copy_number_B
        return output
    end

    required_columns = (
        :pair_sampleA,
        :pair_sampleB,
        :start_pos_metaData,
        :end_pos_metaData,
    )
    all(column -> column in propertynames(output), required_columns) ||
        throw(ArgumentError(
            "CNV region results lack sample indices or metaData row coordinates required to attach copy numbers"))
    expected_state_size = (
        size(CNVscaffold, 1),
        size(copynumber_result.callable, 1),
        size(copynumber_result.callable, 1),
    )
    size(event_state) == expected_state_size || throw(DimensionMismatch(
        "event_state dimensions do not match CNVscaffold or the sample count"))

    for row in 1:nrow(output)
        sample_A = Int(output[row, :pair_sampleA])
        sample_B = Int(output[row, :pair_sampleB])
        region_start = Int(output[row, :start_pos_metaData])
        region_end = Int(output[row, :end_pos_metaData])
        event_pair = _cnv_region_weighted_different_copy_number_pair(
            copynumber_result,
            CNVscaffold,
            event_state,
            sample_A,
            sample_B,
            region_start,
            region_end,
        )
        copy_number_A[row], copy_number_B[row] = event_pair
    end

    output[!, :copy_number_A] = copy_number_A
    output[!, :copy_number_B] = copy_number_B
    has_different_copy_number = [
        !ismissing(copy_number_A[row]) &&
        !ismissing(copy_number_B[row]) &&
        copy_number_A[row] != copy_number_B[row]
        for row in 1:nrow(output)
    ]
    return output[has_different_copy_number, :]
end

function plotCNVregionMulti_RD_background(
    metaDataDF,
    CNVscaffold,
    comparison_or_filter_result,
    sample_pair;
    copynumber_result=nothing,
    merge_regions::Bool=false,
    min_gap::Real=1_000_000,
    max_gap::Real=3_000_000,
    diff_threshold::Real=15,
    gap_diff_threshold::Real=10,
    kwargs...,
)
    state = hasproperty(comparison_or_filter_result, :state) ?
            comparison_or_filter_result.state : comparison_or_filter_result
    state isa AbstractArray && ndims(state) == 3 || throw(ArgumentError(
        "The third argument must be a background-comparison result with a state field or a three-dimensional state array"))
    all(in((-1, 0, 1)), state) || throw(ArgumentError(
        "State values must be -1, 0, or 1"))
    sample_pair isa Matrix && size(sample_pair, 2) == 2 ||
        throw(ArgumentError("sample_pair must be an n-by-2 matrix"))
    isnothing(copynumber_result) && throw(ArgumentError(
        "Pass the result of copynumber_determination_background through the copynumber_result keyword"))
    hasproperty(copynumber_result, :copy_number) &&
        hasproperty(copynumber_result, :callable) || throw(ArgumentError(
        "copynumber_result must come from copynumber_determination_background"))
    sample_count = ncol(metaDataDF) - 3
    fragment_count = size(CNVscaffold, 1)
    copy_number_calls = _copy_number_call_matrix(copynumber_result)
    size(copy_number_calls) == (sample_count, fragment_count) ||
        throw(DimensionMismatch(
            "copynumber_result.copy_number_calls dimensions do not match the samples or CNVscaffold"))
    size(copynumber_result.callable) == (sample_count, fragment_count) ||
        throw(DimensionMismatch(
            "copynumber_result.callable dimensions do not match the samples or CNVscaffold"))
    event_state = hasproperty(comparison_or_filter_result, :original_state) ?
                  comparison_or_filter_result.original_state : state
    size(event_state) == size(state) || throw(DimensionMismatch(
        "original_state dimensions do not match state"))

    keyword_values = (; kwargs...)
    positive_cnv_content = hasproperty(keyword_values, :positive_cnv_content) ?
                           keyword_values.positive_cnv_content : 0.8
    0 <= positive_cnv_content <= 1 || throw(ArgumentError(
        "positive_cnv_content must be in [0, 1]"))
    min_gap < max_gap || throw(ArgumentError("min_gap must be less than max_gap"))

    plot_state = _cnv_plot_state_with_background_barriers(
        metaDataDF, CNVscaffold, state, positive_cnv_content)
    background_regions = _cnv_background_regions(
        metaDataDF, CNVscaffold, state, sample_pair)
    forwarded_keywords = merge(
        keyword_values,
        (
            positive_cnv_content=1.0,
            merge_regions=merge_regions,
            merge_min_gap=min_gap,
            merge_max_gap=max_gap,
            merge_diff_threshold=diff_threshold,
            merge_gap_diff_threshold=gap_diff_threshold,
        ),
    )
    cnv_regions, filtered_cnv_regions, plt = plotCNVregionMulti_RD(
        metaDataDF,
        CNVscaffold,
        plot_state,
        sample_pair;
        forwarded_keywords...,
    )
    cnv_regions = _cnv_add_region_copy_numbers(
        cnv_regions, CNVscaffold, copynumber_result, event_state)
    filtered_cnv_regions = _cnv_add_region_copy_numbers(
        filtered_cnv_regions, CNVscaffold, copynumber_result, event_state)

    selected_chr = hasproperty(keyword_values, :chr) ?
                   String(keyword_values.chr) : "all"
    for pair_index in axes(sample_pair, 1)
        sample1, sample2 = sample_pair[pair_index, :]
        selected_background = background_regions[
            (background_regions.pair_sampleA .== sample1) .&
            (background_regions.pair_sampleB .== sample2), :]
        for region in eachrow(selected_background)
            if selected_chr == "all"
                vspan!(plt, [region.start_fragment, region.end_fragment];
                       subplot=pair_index, color=:gray, alpha=0.15,
                       label=false)
            elseif region.chr == selected_chr
                vspan!(plt, [region.start_pos, region.end_pos];
                       subplot=pair_index, color=:gray, alpha=0.15,
                       label=false)
            end
        end
    end
    display(plt)

    return (
        cnv_regions=cnv_regions,
        filtered_cnv_regions=filtered_cnv_regions,
        background_regions=background_regions,
        plot=plt,
        plot_state=plot_state,
        original_state=Int8.(state),
        background_semantics=merge_regions ?
            :background_aware_state_then_gap_merge :
            :uncallable_merge_barrier,
    )
end

function cnv_region(metaDataDF, CNVscarffold, fragment_status, threshold_proportion)
    try

        if !isa(metaDataDF, DataFrame)
            throw(ArgumentError("metaDataDF must be a DataFrame"))
        end

        if !isa(CNVscarffold, Matrix) || size(CNVscarffold, 2) != 3
            throw(ArgumentError("CNVscarffold must be an n-by-3 matrix"))
        end

        if !isa(fragment_status, Vector)
            throw(ArgumentError("fragment_status must be a Vector"))
        end

        if !isa(threshold_proportion, Number) || threshold_proportion < 0 || threshold_proportion > 1
            throw(ArgumentError("threshold_proportion must be in [0, 1]"))
        end

        if size(metaDataDF, 2) < 3
            throw(DimensionMismatch("metaDataDF must contain at least three columns"))
        end

        if length(fragment_status) != size(CNVscarffold, 1)
            throw(DimensionMismatch("fragment_status length must match the number of CNVscarffold rows"))
        end

        fragment_coords = hcat(
            metaDataDF[CNVscarffold[:, 1], 1],
            CNVscarffold[:, 1],
            CNVscarffold[:, 2],
            CNVscarffold[:, 2] - CNVscarffold[:, 1] .+ 1
        )

        temp_regions = synth_cnv(fragment_coords, fragment_status, threshold_proportion)

        if temp_regions == -1
            @error "An error occurred while combining CNV fragments"
            return []
        end

        chr = String[]
        start_pos = Int64[]
        end_pos = Int64[]
        num_regions = size(temp_regions, 1)
        start_pos_metaData = Int64[]
        end_pos_metaData = Int64[]
        if !isempty(temp_regions)
            for i in 1:num_regions
                temp_chr = String(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 1])
                temp_start_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 1], 1], 2])
                temp_end_pos = Int64(metaDataDF[CNVscarffold[temp_regions[i, 2], 2], 3])
                if temp_end_pos >= temp_start_pos
                    push!(chr, temp_chr)
                    push!(start_pos, temp_start_pos)
                    push!(end_pos, temp_end_pos)
                    push!(start_pos_metaData, CNVscarffold[temp_regions[i, 1], 1])
                    push!(end_pos_metaData, CNVscarffold[temp_regions[i, 2], 2])
                end
            end
        end

        coordinate_cnv_regions = DataFrame()
        coordinate_cnv_regions[!,:chr] = chr
        coordinate_cnv_regions[!,:start_pos] = start_pos
        coordinate_cnv_regions[!,:end_pos] = end_pos
        coordinate_cnv_regions[!,:length] = end_pos .- start_pos
        coordinate_cnv_regions[!,:start_pos_metaData] = start_pos_metaData
        coordinate_cnv_regions[!,:end_pos_metaData] = end_pos_metaData

        return coordinate_cnv_regions

    catch e
        @error "An error occurred in cnv_region: " * string(e)
        rethrow(e)
    end
end

# Background-aware CNV calling and comparison pipeline.
function copyNumber_comparison_filter_fn_background(
    metaDataDF,
    CNVscaffold,
    copynumber_result,
    comparison_result;
    highRDthreshold::Real=100,
    iterative::Bool=true,
    max_iterations::Union{Nothing,Integer}=nothing,
)
    metaDataDF isa DataFrame || throw(ArgumentError(
        "metaDataDF must be a DataFrame"))
    CNVscaffold isa Matrix && size(CNVscaffold, 2) == 3 ||
        throw(ArgumentError("CNVscaffold must be an n-by-3 matrix"))
    highRDthreshold >= 0 || throw(ArgumentError(
        "highRDthreshold must be nonnegative"))
    if !isnothing(max_iterations)
        max_iterations isa Integer && max_iterations >= 1 ||
            throw(ArgumentError("max_iterations must be a positive integer or nothing"))
    end
    hasproperty(copynumber_result, :copy_number) &&
        hasproperty(copynumber_result, :callable) &&
        hasproperty(copynumber_result, :fragment_depth) ||
        throw(ArgumentError(
            "copynumber_result must come from copynumber_determination_background"))
    hasproperty(comparison_result, :state) &&
        hasproperty(comparison_result, :callable) ||
        throw(ArgumentError(
            "comparison_result must come from copyNumber_comparison_background"))

    fragment_count = size(CNVscaffold, 1)
    sample_count = ncol(metaDataDF) - 3
    sample_count >= 2 || throw(ArgumentError(
        "At least two samples are required for between-sample CN comparison filtering"))
    copy_number_calls = _copy_number_call_matrix(copynumber_result)
    size(copy_number_calls) ==
        (sample_count, fragment_count) || throw(DimensionMismatch(
            "copynumber_result.copy_number_calls dimensions do not match the sample or fragment count"))
    size(copynumber_result.callable) ==
        (sample_count, fragment_count) || throw(DimensionMismatch(
            "copynumber_result.callable dimensions do not match"))
    size(copynumber_result.fragment_depth) ==
        (fragment_count, sample_count) || throw(DimensionMismatch(
            "copynumber_result.fragment_depth dimensions do not match"))
    expected_state_size = (fragment_count, sample_count, sample_count)
    size(comparison_result.state) == expected_state_size ||
        throw(DimensionMismatch("comparison_result.state dimensions do not match"))
    size(comparison_result.callable) == expected_state_size ||
        throw(DimensionMismatch("comparison_result.callable dimensions do not match"))
    all(in((-1, 0, 1)), comparison_result.state) || throw(ArgumentError(
        "comparison_result.state values must be -1, 0, or 1"))

    original_state = Int8.(comparison_result.state)
    final_state = copy(original_state)
    promoted_by_rd = falses(expected_state_size)
    promotion_iteration = zeros(Int, expected_state_size)
    fragment_rd_difference = fill(NaN, expected_state_size)
    for sample1 in 1:sample_count, sample2 in 1:sample_count,
        fragment in 1:fragment_count
        depth1 = copynumber_result.fragment_depth[fragment, sample1]
        depth2 = copynumber_result.fragment_depth[fragment, sample2]
        if isfinite(depth1) && isfinite(depth2)
            fragment_rd_difference[fragment, sample1, sample2] =
                abs(depth1 - depth2)
        end
    end

    completed_iterations = 0
    maximum_rounds = isnothing(max_iterations) ?
                     max(fragment_count, 1) : Int(max_iterations)
    while completed_iterations < maximum_rounds
        completed_iterations += 1
        previous_state = copy(final_state)
        changed = false
        if fragment_count >= 3
            for sample1 in 1:(sample_count - 1),
                sample2 in (sample1 + 1):sample_count,
                fragment in 2:(fragment_count - 1)

                previous_state[fragment, sample1, sample2] == 0 || continue
                comparison_result.callable[fragment, sample1, sample2] || continue
                copynumber_result.callable[sample1, fragment] || continue
                copynumber_result.callable[sample2, fragment] || continue
                has_positive_neighbor =
                    previous_state[fragment - 1, sample1, sample2] == 1 ||
                    previous_state[fragment + 1, sample1, sample2] == 1
                has_positive_neighbor || continue
                rd_difference = fragment_rd_difference[
                    fragment, sample1, sample2]
                isfinite(rd_difference) &&
                    rd_difference > highRDthreshold || continue

                final_state[fragment, sample1, sample2] = 1
                final_state[fragment, sample2, sample1] = 1
                promoted_by_rd[fragment, sample1, sample2] = true
                promoted_by_rd[fragment, sample2, sample1] = true
                promotion_iteration[fragment, sample1, sample2] =
                    completed_iterations
                promotion_iteration[fragment, sample2, sample1] =
                    completed_iterations
                changed = true
            end
        end
        (!iterative || !changed) && break
    end

    return (
        state=final_state,
        original_state=original_state,
        callable=comparison_result.callable,
        promoted_by_rd=promoted_by_rd,
        promotion_iteration=promotion_iteration,
        fragment_rd_difference=fragment_rd_difference,
        iterations=completed_iterations,
        highRDthreshold=Float64(highRDthreshold),
        iterative=iterative,
        state_legend=Dict(Int8(-1) => :skipped,
                          Int8(0) => :no_difference,
                          Int8(1) => :difference),
    )
end

function copynumber_determination_background(
    metaDataDF,
    CNVscaffold,
    gmmParam,
    flanking_length;
    thresholdBFactor::Real=2.0,
    max_background_probability::Real=0.50,
)
    metaDataDF isa DataFrame || throw(ArgumentError(
        "metaDataDF must be a DataFrame"))
    CNVscaffold isa Matrix && size(CNVscaffold, 2) == 3 ||
        throw(ArgumentError("CNVscaffold must be an n-by-3 matrix"))
    gmmParam isa Vector || throw(ArgumentError("gmmParam must be a Vector"))
    flanking_length isa Integer && flanking_length > 0 ||
        throw(ArgumentError("flanking_length must be a positive integer"))
    thresholdBFactor > 0 || throw(ArgumentError(
        "thresholdBFactor must be positive"))
    0 <= max_background_probability <= 1 || throw(ArgumentError(
        "max_background_probability must be in [0, 1]"))
    ncol(metaDataDF) >= 4 || throw(DimensionMismatch(
        "metaDataDF must contain three coordinate columns and at least one RD column"))

    subData = Matrix{Float64}(metaDataDF[:, 4:end])
    num_rows, num_samples = size(subData)
    num_fragments = size(CNVscaffold, 1)
    length(gmmParam) == num_samples || throw(DimensionMismatch(
        "gmmParam length must match the number of samples"))
    flanking_length * 2 < num_rows || throw(ArgumentError(
        "flanking_length is too large for the data length"))
    if !isempty(CNVscaffold)
        minimum(CNVscaffold[:, 1:2]) >= 1 || throw(BoundsError(
            "CNVscaffold indices cannot be less than 1"))
        maximum(CNVscaffold[:, 1:2]) <= num_rows || throw(BoundsError(
            "CNVscaffold indices exceed the metaDataDF range"))
    end

    all_mean_sub = fill(NaN, num_rows, num_samples)
    for row in (flanking_length + 1):(num_rows - flanking_length)
        for sample in 1:num_samples
            all_mean_sub[row, sample] = mean(@view subData[
                (row - flanking_length):(row + flanking_length), sample])
        end
    end
    fragment_depth = fill(NaN, num_fragments, num_samples)
    for fragment in 1:num_fragments
        fragment_start = Int(CNVscaffold[fragment, 1])
        fragment_end = Int(CNVscaffold[fragment, 2])
        for sample in 1:num_samples
            values = @view all_mean_sub[fragment_start:fragment_end, sample]
            finite_values = values[isfinite.(values)]
            isempty(finite_values) ||
                (fragment_depth[fragment, sample] = mean(finite_values))
        end
    end

    max_components = maximum(length(parameter.means) for parameter in gmmParam)
    copy_number_calls = Matrix{Union{Missing,Int}}(
        missing, num_samples, num_fragments)
    callable = falses(num_samples, num_fragments)
    status = fill(:unclassified, num_samples, num_fragments)
    background_probability = fill(NaN, num_samples, num_fragments)
    maximum_cn_probability = fill(NaN, num_samples, num_fragments)
    maximum_conditional_cn_probability = fill(NaN, num_samples, num_fragments)
    nearest_cn_z = fill(NaN, num_samples, num_fragments)
    background_full_support = falses(num_samples, num_fragments)
    component_probabilities = zeros(num_samples, num_fragments, max_components)
    conditional_component_probabilities = zeros(
        num_samples, num_fragments, max_components)

    for sample in 1:num_samples
        parameter = gmmParam[sample]
        labels = hasproperty(parameter, :copy_labels) ?
                 Int.(parameter.copy_labels) :
                 collect(1:length(parameter.means))
        length(labels) == length(parameter.means) || throw(DimensionMismatch(
            "Sample $sample  copy_labels and means have different lengths"))
        for fragment in 1:num_fragments
            observation = fragment_depth[fragment, sample]
            if !isfinite(observation)
                status[sample, fragment] = :invalid_depth
                continue
            end
            posterior = _cn_state_posterior(observation, parameter)
            cn_probabilities = posterior.cn_probabilities
            conditional_cn_probabilities =
                posterior.conditional_cn_probabilities
            background_p = posterior.background_probability
            component_probabilities[sample, fragment,
                                    1:length(cn_probabilities)] .= cn_probabilities
            conditional_component_probabilities[
                sample, fragment, 1:length(conditional_cn_probabilities)] .=
                conditional_cn_probabilities
            background_probability[sample, fragment] = background_p
            maximum_p, component_index = findmax(cn_probabilities)
            maximum_conditional_p, conditional_component_index = findmax(
                conditional_cn_probabilities)
            maximum_cn_probability[sample, fragment] = maximum_p
            maximum_conditional_cn_probability[sample, fragment] =
                maximum_conditional_p
            nearest_cn_z[sample, fragment] = posterior.nearest_cn_z
            background_full_support[sample, fragment] =
                posterior.background_full_support

            bayes_factor = maximum_conditional_p /
                           max(1 - maximum_conditional_p, eps())

            if background_p >= max_background_probability
                status[sample, fragment] = :background
            elseif bayes_factor <= thresholdBFactor
                status[sample, fragment] = :uncertain_copy_number
            else
                copy_number_calls[sample, fragment] =
                    labels[conditional_component_index]
                callable[sample, fragment] = true
                status[sample, fragment] = :copy_number
            end
        end
    end

    sample_names = String.(names(metaDataDF)[4:end])
    length(unique(sample_names)) == length(sample_names) || throw(ArgumentError(
        "Input sample column names must be unique"))
    coordinate_column_names = Set([
        "chr", "start", "end", "scaffoldstart", "scaffoldend"])
    conflicting_names = intersect(Set(sample_names), coordinate_column_names)
    isempty(conflicting_names) || throw(ArgumentError(
        "Sample column names conflict with output coordinate columns: " *
        join(sort!(collect(conflicting_names)), ", ")))

    fragment_chr = Vector{String}(undef, num_fragments)
    genomic_start = Vector{Int}(undef, num_fragments)
    genomic_end = Vector{Int}(undef, num_fragments)
    scaffold_start = Int.(CNVscaffold[:, 1])
    scaffold_end = Int.(CNVscaffold[:, 2])
    for fragment in 1:num_fragments
        first_row = scaffold_start[fragment]
        last_row = scaffold_end[fragment]
        first_row <= last_row || throw(ArgumentError(
            "CNVscaffold fragment $fragment  has a start index greater than its end index"))
        first_chr = string(metaDataDF[first_row, 1])
        last_chr = string(metaDataDF[last_row, 1])
        first_chr == last_chr || throw(ArgumentError(
            "CNVscaffold fragment $fragment  crosses chromosome boundaries: " *
            "$first_chr -> $last_chr"))
        fragment_chr[fragment] = first_chr
        genomic_start[fragment] = Int(metaDataDF[first_row, 2])

        genomic_end[fragment] = Int(metaDataDF[last_row, 2])
    end

    copy_number_table = DataFrame(
        :chr => fragment_chr,
        :start => genomic_start,
        :end => genomic_end,
        :scaffoldstart => scaffold_start,
        :scaffoldend => scaffold_end,
    )
    for (sample, sample_name) in enumerate(sample_names)
        copy_number_table[!, sample_name] = vec(copy_number_calls[sample, :])
    end

    return (
        copy_number=copy_number_table,
        copy_number_calls=copy_number_calls,
        callable=callable,
        status=status,
        background_probability=background_probability,
        maximum_cn_probability=maximum_cn_probability,
        maximum_conditional_cn_probability=
            maximum_conditional_cn_probability,
        nearest_cn_z=nearest_cn_z,
        background_full_support=background_full_support,
        background_eligible=background_full_support,
        component_probabilities=component_probabilities,
        conditional_component_probabilities=
            conditional_component_probabilities,
        fragment_depth=fragment_depth,
        sample_names=sample_names,
    )
end

function _cn_component_index_for_label(parameter, copy_label::Int)
    labels = hasproperty(parameter, :copy_labels) ?
             Int.(parameter.copy_labels) :
             collect(1:length(parameter.means))
    return findfirst(==(copy_label), labels)
end

function _prob_same_absolute_cn_background(
    fragment_depths,
    parameter1,
    parameter2,
    window_size::Int,
)
    labels1 = hasproperty(parameter1, :copy_labels) ?
              Int.(parameter1.copy_labels) :
              collect(1:length(parameter1.means))
    labels2 = hasproperty(parameter2, :copy_labels) ?
              Int.(parameter2.copy_labels) :
              collect(1:length(parameter2.means))
    common_labels = sort(intersect(labels1, labels2))
    isempty(common_labels) && return (
        probability=NaN, common_labels=Int[], label_probabilities=Float64[])

    fragment_length = size(fragment_depths, 1)
    pseudo_length = sqrt(max(fragment_length / window_size, eps()))
    mean_within = vec(mean(fragment_depths, dims=1))
    variance_within = Float64[
        fragment_length > 1 ? var(@view fragment_depths[:, sample]) : 0.0
        for sample in 1:2
    ] ./ window_size

    label_probabilities = Float64[]
    prior_weights = Float64[]
    for label in common_labels
        index1 = _cn_component_index_for_label(parameter1, label)
        index2 = _cn_component_index_for_label(parameter2, label)
        expected_means = [parameter1.means[index1], parameter2.means[index2]]
        model_variances = [parameter1.stds[index1]^2,
                           parameter2.stds[index2]^2]
        between_variance = max.(1e-16, model_variances .- variance_within)
        statistic = sum((mean_within .- expected_means).^2 ./
                        (variance_within ./ pseudo_length .+
                         between_variance))
        push!(label_probabilities, ccdf(Chisq(2), statistic))
        push!(prior_weights,
              parameter1.weights[index1] * parameter2.weights[index2])
    end
    weight_sum = sum(prior_weights)
    probability = weight_sum > 0 && isfinite(weight_sum) ?
                  sum(label_probabilities .* prior_weights) / weight_sum :
                  mean(label_probabilities)
    return (
        probability=probability,
        common_labels=common_labels,
        label_probabilities=label_probabilities,
    )
end

function copyNumber_comparison_background(
    metaDataDF,
    CNVscaffold,
    copynumber_result,
    gmmParam,
    windowSize;
    pvalue_theSame::Real=0.05,
)
    metaDataDF isa DataFrame || throw(ArgumentError(
        "metaDataDF must be a DataFrame"))
    CNVscaffold isa Matrix && size(CNVscaffold, 2) == 3 ||
        throw(ArgumentError("CNVscaffold must be an n-by-3 matrix"))
    gmmParam isa Vector || throw(ArgumentError("gmmParam must be a Vector"))
    windowSize isa Integer && windowSize > 0 || throw(ArgumentError(
        "windowSize must be a positive integer"))
    0 < pvalue_theSame < 1 || throw(ArgumentError(
        "pvalue_theSame must be in (0, 1)"))
    hasproperty(copynumber_result, :copy_number) &&
        hasproperty(copynumber_result, :callable) || throw(ArgumentError(
            "copynumber_result must come from copynumber_determination_background"))

    subData = Matrix{Float64}(metaDataDF[:, 4:end])
    num_rows, num_samples = size(subData)
    num_fragments = size(CNVscaffold, 1)
    copy_number_calls = _copy_number_call_matrix(copynumber_result)
    size(copy_number_calls) == (num_samples, num_fragments) ||
        throw(DimensionMismatch("copynumber_result dimensions do not match"))
    length(gmmParam) == num_samples || throw(DimensionMismatch(
        "gmmParam length must match the number of samples"))

    state = fill(Int8(-1), num_fragments, num_samples, num_samples)
    probability = fill(NaN, num_fragments, num_samples, num_samples)
    callable_pair = falses(num_fragments, num_samples, num_samples)
    skip_reason = fill(:not_evaluated, num_fragments, num_samples, num_samples)
    common_copy_numbers = Array{Vector{Int}}(undef,
        num_fragments, num_samples, num_samples)

    for sample1 in 1:num_samples, sample2 in 1:num_samples,
        fragment in 1:num_fragments
        common_copy_numbers[fragment, sample1, sample2] = Int[]
        if sample1 == sample2
            if copynumber_result.callable[sample1, fragment]
                state[fragment, sample1, sample2] = 0
                probability[fragment, sample1, sample2] = 1.0
                callable_pair[fragment, sample1, sample2] = true
                skip_reason[fragment, sample1, sample2] = :self_comparison
            else
                skip_reason[fragment, sample1, sample2] = :uncallable_sample
            end
            continue
        end
        if !copynumber_result.callable[sample1, fragment] ||
           !copynumber_result.callable[sample2, fragment]
            skip_reason[fragment, sample1, sample2] = :uncallable_sample
            continue
        end

        cn1 = copy_number_calls[sample1, fragment]
        cn2 = copy_number_calls[sample2, fragment]
        if cn1 == cn2
            state[fragment, sample1, sample2] = 0
            probability[fragment, sample1, sample2] = 1.0
            callable_pair[fragment, sample1, sample2] = true
            skip_reason[fragment, sample1, sample2] = :same_called_copy_number
            continue
        end

        fragment_start = max(1, Int(CNVscaffold[fragment, 1]))
        fragment_end = min(num_rows, Int(CNVscaffold[fragment, 2]))
        test_result = _prob_same_absolute_cn_background(
            @view(subData[fragment_start:fragment_end, [sample1, sample2]]),
            gmmParam[sample1],
            gmmParam[sample2],
            Int(windowSize),
        )
        common_copy_numbers[fragment, sample1, sample2] =
            test_result.common_labels
        if isempty(test_result.common_labels) ||
           !isfinite(test_result.probability)
            skip_reason[fragment, sample1, sample2] =
                :no_shared_copy_number_model
            continue
        end
        probability[fragment, sample1, sample2] = test_result.probability
        callable_pair[fragment, sample1, sample2] = true
        state[fragment, sample1, sample2] =
            test_result.probability < pvalue_theSame ? Int8(1) : Int8(0)
        skip_reason[fragment, sample1, sample2] = :tested
    end

    return (
        state=state,
        probability=probability,
        callable=callable_pair,
        skip_reason=skip_reason,
        common_copy_numbers=common_copy_numbers,
        state_legend=Dict(Int8(-1) => :skipped,
                          Int8(0) => :no_difference,
                          Int8(1) => :difference),
    )
end

function genome_segment_background(
    inputDataDF,
    givenDistParam,
    breaks_threshold,
    flanking_length,
    minimum_length;
    min_background_probability::Real=0.50,
)
    inputDataDF isa DataFrame || throw(ArgumentError(
        "inputDataDF must be a DataFrame"))
    givenDistParam isa Vector || throw(ArgumentError(
        "givenDistParam must be a Vector"))
    0 <= breaks_threshold <= 1 || throw(ArgumentError(
        "breaks_threshold must be in [0, 1]"))
    0 <= min_background_probability <= 1 || throw(ArgumentError(
        "min_background_probability must be in [0, 1]"))
    flanking_length isa Integer && flanking_length > 0 ||
        throw(ArgumentError("flanking_length must be a positive integer"))
    minimum_length isa Integer && minimum_length > 0 ||
        throw(ArgumentError("minimum_length must be a positive integer"))
    ncol(inputDataDF) >= 4 || throw(DimensionMismatch(
        "inputDataDF must contain three coordinate columns and at least one RD column"))

    subData = Matrix{Float64}(inputDataDF[:, 4:end])
    row_count, sample_count = size(subData)
    flanking_length * 2 < row_count || throw(ArgumentError(
        "flanking_length is too large for the data length"))
    length(givenDistParam) == sample_count || throw(DimensionMismatch(
        "givenDistParam length must match the number of samples"))

    left_means = fill(NaN, row_count, sample_count)
    right_means = fill(NaN, row_count, sample_count)
    coincidence_probability = fill(NaN, row_count, sample_count)
    left_background_probability = fill(NaN, row_count, sample_count)
    right_background_probability = fill(NaN, row_count, sample_count)
    dominant_left_state = Matrix{Any}(undef, row_count, sample_count)
    dominant_right_state = Matrix{Any}(undef, row_count, sample_count)
    fill!(dominant_left_state, :not_evaluated)
    fill!(dominant_right_state, :not_evaluated)
    sample_breaks = falses(row_count, sample_count)
    boundary_type = fill(:none, row_count, sample_count)

    for row in (flanking_length + 1):(row_count - flanking_length)
        physical_gap = inputDataDF[row + flanking_length, 3] -
                       inputDataDF[row - flanking_length, 3] > 100 * 10^3
        for sample in 1:sample_count
            left_means[row, sample] = mean(@view subData[
                (row - flanking_length):(row - 1), sample])
            right_means[row, sample] = mean(@view subData[
                (row + 1):(row + flanking_length), sample])
            details = rate_coincident_background(
                [left_means[row, sample], right_means[row, sample]],
                givenDistParam[sample];
                min_background_probability=min_background_probability,
                return_details=true,
            )
            coincidence_probability[row, sample] = details.probability
            left_background_probability[row, sample] =
                details.background_probabilities[1]
            right_background_probability[row, sample] =
                details.background_probabilities[2]
            dominant_left_state[row, sample] = details.dominant_states[1]
            dominant_right_state[row, sample] = details.dominant_states[2]

            probability_break = details.probability < breaks_threshold
            sample_breaks[row, sample] = probability_break || physical_gap
            if physical_gap
                boundary_type[row, sample] = :physical_gap
            elseif probability_break
                left_is_background = details.dominant_states[1] === :background
                right_is_background = details.dominant_states[2] === :background
                if xor(left_is_background, right_is_background)
                    boundary_type[row, sample] = :background_transition
                elseif !left_is_background && !right_is_background &&
                       details.dominant_states[1] != details.dominant_states[2]
                    boundary_type[row, sample] = :integer_cn_transition
                else
                    boundary_type[row, sample] = :probabilistic_state_change
                end
            end
        end
    end

    break_support_count = vec(sum(sample_breaks, dims=2))

    break_runs = run_length_nonzero(break_support_count)
    scaffold = zeros(Int, 0, 3)
    if size(break_runs, 1) > 1
        run_regions = zeros(Int, size(break_runs, 1) - 1, 3)
        for run in 1:(size(break_runs, 1) - 1)
            region_start = Int(break_runs[run, 3])
            region_end = Int(break_runs[run + 1, 2])
            run_regions[run, :] .=
                (region_start, region_end, region_end - region_start)
        end
        retained = findall(run_regions[:, 3] .>= minimum_length)
        if !isempty(retained)
            candidate_scaffold = run_regions[retained, :]
            same_chromosome = inputDataDF[candidate_scaffold[:, 1], 1] .==
                              inputDataDF[candidate_scaffold[:, 2], 1]
            scaffold = candidate_scaffold[findall(same_chromosome), :]
        end
    end

    return (
        CNVscaffold=scaffold,
        sample_breaks=sample_breaks,
        break_support_count=break_support_count,
        boundary_type=boundary_type,
        coincidence_probability=coincidence_probability,
        left_background_probability=left_background_probability,
        right_background_probability=right_background_probability,
        min_background_probability=Float64(min_background_probability),
        dominant_left_state=dominant_left_state,
        dominant_right_state=dominant_right_state,
        left_means=left_means,
        right_means=right_means,
        state_semantics=(
            background_is_distinct_from_integer_cn=true,
            background_matches_background=true,
        ),
    )
end

# Legacy normalization overload retained for compatibility.
function filtering_normalizing_seqDepth(dataframe_raw_sequencingdepth, depth_filter1, depth_filter2)
    try

        if !isa(dataframe_raw_sequencingdepth, DataFrame)
            throw(ArgumentError("dataframe_raw_sequencingdepth must be a DataFrame"))
        end

        if (!isa(depth_filter1, Vector) || length(depth_filter1) != 2) && (isa(depth_filter2, Vector) || length(depth_filter2) != 2)
            throw(DimensionMismatch("depth_filter must be a two-element vector"))
        end

        if depth_filter1[1] >= depth_filter1[2]
            throw(ArgumentError("The minimum depth in depth_filter1 must be less than the maximum depth"))
        end

        if depth_filter2[1] >= depth_filter2[2]
            throw(ArgumentError("The minimum depth in depth_filter2 must be less than the maximum depth"))
        end

        if size(dataframe_raw_sequencingdepth, 2) < 4
            throw(DimensionMismatch("dataframe_raw_sequencingdepth must contain at least four columns"))
        end

        tempDepth = sum(Matrix(dataframe_raw_sequencingdepth[:, 4:end]), dims = 2)

        tempidx = findall((tempDepth[:] .> depth_filter1[1]) .& (tempDepth[:] .< depth_filter1[2]) .== 1)

        if isempty(tempidx)
            @warn "No fragments satisfy the depth filters"
            return DataFrame()
        end

        subCombineData = dataframe_raw_sequencingdepth[tempidx, :]

        tempMean = mean(Matrix(subCombineData[:, 4:end]), dims = 1)
        normalSubCombineData = hcat(subCombineData[:, 1:3], Float64.(subCombineData[:, 4:end]))

        for i in 1:size(normalSubCombineData, 1)

            if all(tempMean .> 0)
                normalSubCombineData[i, 4:end] = (Vector{Float64}(normalSubCombineData[i, 4:end])' ./ tempMean) .* 100
            else
                @warn "Mean depth contains zero; skipping normalization"
                break
            end
        end

        miniDepth = minimum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        maxiDepth = maximum(Matrix(normalSubCombineData[:,4:end]), dims = 2)
        tempidx = findall((miniDepth[:] .> depth_filter2[1]) .& (maxiDepth[:] .< depth_filter2[2]))

        return normalSubCombineData[tempidx,:]

    catch e
        @error "An error occurred in filtering_normalizing_seqDepth: " * string(e)
        rethrow(e)
    end
end
