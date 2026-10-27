
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions

# Data ingestion and normalization.
function process_genomic_data(
    input_directory::AbstractString;
    file_names::Vector{String}=String[],
    file_extension::AbstractString=".csv",
)::DataFrame
    isdir(input_directory) ||
        throw(ArgumentError("Input directory does not exist: $input_directory"))

    extension = startswith(file_extension, ".") ?
        lowercase(file_extension) : "." * lowercase(file_extension)
    selected_files = isempty(file_names) ?
        sort([name for name in readdir(input_directory)
              if isfile(joinpath(input_directory, name)) &&
                 lowercase(splitext(name)[2]) == extension]) : copy(file_names)
    isempty(selected_files) &&
        throw(ArgumentError("No data files were found in the input directory: $input_directory"))

    required_columns = [:Chromosome, :BinLowEdge, :BinUpEdge, :Content]
    tables = DataFrame[]
    sample_names = String[]
    for file_name in selected_files
        full_path = isabspath(file_name) ? file_name :
            joinpath(input_directory, file_name)
        isfile(full_path) ||
            throw(ArgumentError("Input file does not exist: $full_path"))
        table = CSV.read(full_path, DataFrame)
        missing_columns = setdiff(required_columns, Symbol.(names(table)))
        isempty(missing_columns) || throw(ArgumentError(
            "$file_name is missing required columns: " * join(string.(missing_columns), ", ")))
        push!(tables, table)
        push!(sample_names, splitext(basename(file_name))[1])
    end
    length(unique(sample_names)) == length(sample_names) ||
        throw(ArgumentError("Sample names generated from file names must be unique"))

    coordinates = select(tables[1], :Chromosome, :BinLowEdge, :BinUpEdge)
    for sample_idx in 2:length(tables)
        nrow(tables[sample_idx]) == nrow(coordinates) ||
            throw(DimensionMismatch("Input files have inconsistent bin counts"))
        same_coordinates =
            all(tables[sample_idx].Chromosome .== coordinates.Chromosome) &&
            all(tables[sample_idx].BinLowEdge .== coordinates.BinLowEdge) &&
            all(tables[sample_idx].BinUpEdge .== coordinates.BinUpEdge)
        same_coordinates ||
            throw(ArgumentError("Input files have inconsistent genomic coordinates"))
    end

    autosomes = Set("chr" .* string.(1:22))
    keep = in.(String.(coordinates.Chromosome), Ref(autosomes))
    combined_data = coordinates[keep, :]
    for (sample_name, table) in zip(sample_names, tables)
        combined_data[!, Symbol(sample_name)] = table[keep, :Content]
    end
    return combined_data
end

function filtering_normalizing_seqDepth(
    data::DataFrame;
    initial_quantiles::Tuple{<:Real,<:Real}=(0.005, 0.995),
    normalized_quantiles::Tuple{<:Real,<:Real}=(0.005, 0.995),
    normalized_limits=nothing,
)
    size(data, 2) >= 4 ||
        throw(DimensionMismatch("The input DataFrame must contain three coordinate columns and at least one RD column"))
    0 <= initial_quantiles[1] < initial_quantiles[2] <= 1 ||
        throw(ArgumentError("initial_quantiles must satisfy0 <= low < high <= 1"))
    0 <= normalized_quantiles[1] < normalized_quantiles[2] <= 1 ||
        throw(ArgumentError("normalized_quantiles must satisfy0 <= low < high <= 1"))
    if !isnothing(normalized_limits)
        normalized_limits isa Tuple || normalized_limits isa AbstractVector ||
            throw(ArgumentError(
                "normalized_limits must be a two-element Tuple or AbstractVector"))
        length(normalized_limits) == 2 || throw(DimensionMismatch(
            "normalized_limits must contain exactly two bounds"))
        all(value -> value isa Real && isfinite(value), normalized_limits) ||
            throw(ArgumentError("normalized_limits must contain two finite real values"))
        normalized_limits[1] < normalized_limits[2] || throw(ArgumentError(
            "normalized_limits must satisfylow < high"))
    end

    raw_matrix = Matrix{Float64}(data[:, 4:end])
    total_depth = vec(sum(raw_matrix, dims=2))
    finite_total = total_depth[isfinite.(total_depth)]
    isempty(finite_total) &&
        throw(ArgumentError("The input contains no finite total RD values"))
    raw_limits = Float64.(quantile(
        finite_total, collect(initial_quantiles)))
    raw_keep = isfinite.(total_depth) .&
               (total_depth .> raw_limits[1]) .&
               (total_depth .< raw_limits[2])
    any(raw_keep) ||
        throw(ArgumentError("No rows remain after initial_quantiles filtering"))

    filtered = copy(data[raw_keep, :])
    filtered_depth = Matrix{Float64}(filtered[:, 4:end])
    sample_means = vec(mean(filtered_depth, dims=1))
    all(isfinite.(sample_means) .& (sample_means .> 0)) ||
        throw(ArgumentError("Every sample must have a positive finite mean RD"))
    normalized_depth = filtered_depth ./ reshape(sample_means, 1, :) .* 100.0
    for (offset, column_name) in enumerate(names(filtered)[4:end])
        filtered[!, column_name] = normalized_depth[:, offset]
    end

    finite_rows = vec(all(isfinite.(normalized_depth), dims=2))
    any(finite_rows) ||
        throw(ArgumentError("No rows have finite RD values after normalization"))
    row_minima = vec(minimum(normalized_depth[finite_rows, :], dims=2))
    row_maxima = vec(maximum(normalized_depth[finite_rows, :], dims=2))
    resolved_normalized_limits = if isnothing(normalized_limits)
        Float64[
            quantile(row_minima, normalized_quantiles[1]),
            quantile(row_maxima, normalized_quantiles[2]),
        ]
    else
        Float64[normalized_limits[1], normalized_limits[2]]
    end
    normalized_row_minima = vec(minimum(normalized_depth, dims=2))
    normalized_row_maxima = vec(maximum(normalized_depth, dims=2))
    normalized_keep = finite_rows .&
                      (normalized_row_minima .> resolved_normalized_limits[1]) .&
                      (normalized_row_maxima .< resolved_normalized_limits[2])
    any(normalized_keep) ||
        throw(ArgumentError("No rows remain after normalized-range filtering"))
    return filtered[normalized_keep, :]
end

function filtering_normalizing_seqDepth(
    data::DataFrame,
    normalized_limits;
    initial_quantiles::Tuple{<:Real,<:Real}=(0.005, 0.995),
    normalized_quantiles::Tuple{<:Real,<:Real}=(0.005, 0.995),
)
    return filtering_normalizing_seqDepth(
        data;
        initial_quantiles=initial_quantiles,
        normalized_quantiles=normalized_quantiles,
        normalized_limits=normalized_limits,
    )
end

# Convert adjacent positive fragments into CNV regions.
function synth_cnv(fragment_coords, fragment_status, threshold_proportion)

    start_indices = Int64.(fragment_coords[:, 2])
    end_indices = Int64.(fragment_coords[:, 3])

    for i in eachindex(start_indices)
        if start_indices[i] > end_indices[i]
            println("Error: CNV fragment start index exceeds its end index; row: $i")
            return -1
        end

        if i > firstindex(start_indices) && start_indices[i] < start_indices[i - 1]
            println("Error: CNV fragments are not ordered by global metaDataDF row index; row: $i")
            return -1
        end
    end

    positive_indices = findall(fragment_status .== 1)

    if isempty(positive_indices)
        println("Error: no positive fragments")
        return -1
    end

    cnv_table = hcat(positive_indices, positive_indices)

    k = size(cnv_table, 1)
    if k <= 1

        return cnv_table
    end

    merged = true
    while merged
        merged = false

        for i in 2:size(cnv_table, 1)

            last_idx_fragment_coords = cnv_table[i-1, 1]
            current_idx_fragment_coords = cnv_table[i, 2]

            temp_lengths = fragment_coords[last_idx_fragment_coords:current_idx_fragment_coords, 4]
            possible_proportion = sum(temp_lengths .* (fragment_status[last_idx_fragment_coords:current_idx_fragment_coords] .== 1)) /
                                  sum(temp_lengths)

            if (possible_proportion >= threshold_proportion) &
                (fragment_coords[last_idx_fragment_coords, 1] == fragment_coords[current_idx_fragment_coords, 1])

                cnv_table[i, 1] = cnv_table[i - 1, 1]
                cnv_table[i - 1, 1] = -1
                cnv_table[i - 1, 2] = -1
                merged = true
            end
        end

        temp_cnvtable = hcat(cnv_table[findall(cnv_table[:, 1] .> 0), 1], cnv_table[findall(cnv_table[:, 2] .> 0), 2])
        cnv_table = temp_cnvtable
    end

    return cnv_table
end

# Integer copy-number constrained GMM helpers.
function _cn_resolve_modal_copy_number(specification, sample_index, sample_name)
    value = if specification isa Integer
        specification
    elseif specification isa AbstractVector
        1 <= sample_index <= length(specification) || throw(ArgumentError(
            "modal_peak_copy_number has no entry for sample $sample_index "))
        specification[sample_index]
    elseif specification isa AbstractDict
        if haskey(specification, sample_name)
            specification[sample_name]
        elseif haskey(specification, Symbol(sample_name))
            specification[Symbol(sample_name)]
        elseif haskey(specification, sample_index)
            specification[sample_index]
        else
            throw(ArgumentError(
                "modal_peak_copy_number dictionary has no entry for sample $sample_name "))
        end
    else
        throw(ArgumentError(
            "modal_peak_copy_number must be an integer, a vector, or a dictionary keyed by sample name"))
    end
    value isa Integer || throw(ArgumentError(
        "Sample  $sample_name  must have an integer copy number for its highest-density peak"))
    return Int(value)
end

function _cn_chromosome_ranges(chromosomes)
    isempty(chromosomes) && return UnitRange{Int}[]
    output = UnitRange{Int}[]
    range_start = 1
    for row_index in 2:length(chromosomes)
        if chromosomes[row_index] != chromosomes[row_index - 1]
            push!(output, range_start:(row_index - 1))
            range_start = row_index
        end
    end
    push!(output, range_start:length(chromosomes))
    return output
end

function _cn_window_means(values::AbstractVector{<:Real}, window_size::Int,
                          chromosome_ranges)
    output = Float64[]
    for chromosome_range in chromosome_ranges
        length(chromosome_range) >= window_size || continue
        chromosome_values = Float64.(@view values[chromosome_range])
        cumulative = vcat(0.0, cumsum(chromosome_values))
        append!(output,
                (cumulative[(window_size + 1):end] .-
                 cumulative[1:(end - window_size)]) ./ window_size)
    end
    return output
end

function _cn_smooth_histogram(counts::AbstractVector{<:Real}, window::Int)
    radius = fld(window, 2)
    return Float64[
        mean(@view counts[max(firstindex(counts), index - radius):min(lastindex(counts), index + radius)])
        for index in eachindex(counts)
    ]
end

function _cn_detect_histogram_peaks(
    counts::AbstractVector{<:Real};
    min_peak_prominence::Real,
    prominence_window_bins::Int,
    min_peak_distance_bins::Int,
    min_peak_mass::Real,
    max_components::Int,
)
    maximum_count = maximum(counts)
    maximum_count > 0 || return Int[]
    total_count = sum(counts)
    candidates = NamedTuple[]

    for index in eachindex(counts)

        counts[index] > 0 || continue
        left_neighbor = index == firstindex(counts) ? -Inf : counts[index - 1]
        right_neighbor = index == lastindex(counts) ? -Inf : counts[index + 1]
        counts[index] >= left_neighbor && counts[index] >= right_neighbor &&
            (counts[index] > left_neighbor || counts[index] > right_neighbor) ||
            continue

        left_range = max(firstindex(counts), index - prominence_window_bins):index
        right_range = index:min(lastindex(counts), index + prominence_window_bins)
        baseline = index == firstindex(counts) ? minimum(@view counts[right_range]) :
                   index == lastindex(counts) ? minimum(@view counts[left_range]) :
                   max(minimum(@view counts[left_range]),
                       minimum(@view counts[right_range]))
        prominence = (counts[index] - baseline) / maximum_count
        mass_range = max(firstindex(counts), index - min_peak_distance_bins):min(lastindex(counts), index + min_peak_distance_bins)
        local_mass = sum(@view counts[mass_range]) / total_count
        prominence >= min_peak_prominence || continue
        local_mass >= min_peak_mass || continue
        push!(candidates,
              (index=index, prominence=prominence, mass=local_mass))
    end

    ranked = sort(candidates; by=x -> (x.prominence, x.mass), rev=true)
    retained = Int[]
    for candidate in ranked
        all(abs(candidate.index - previous) >= min_peak_distance_bins
            for previous in retained) || continue
        push!(retained, candidate.index)
        length(retained) >= max_components && break
    end
    return sort(retained)
end

function _cn_assign_integer_lattice(
    peak_means::Vector{Float64},
    peak_heights::Vector{Float64},
    modal_copy_number::Int;
    copy_number_range,
    number_of_active_components,
    max_copy_jump::Int,
    spacing_tolerance_fraction::Real,
    missing_state_penalty::Real,
    modal_peak_index::Union{Nothing,Integer}=nothing,
    reference_spacing::Union{Nothing,Real}=nothing,
    reference_spacing_tolerance_fraction::Real=0.20,
)
    length(peak_means) >= 2 || throw(ArgumentError(
        "The integer copy-number grid requires at least two candidate peaks"))
    modal_peak_index = isnothing(modal_peak_index) ? argmax(peak_heights) :
                       Int(modal_peak_index)
    1 <= modal_peak_index <= length(peak_means) || throw(ArgumentError(
        "modal_peak_index is outside the candidate-peak range"))
    minimum_components = isnothing(number_of_active_components) ? 2 :
                         Int(number_of_active_components)
    maximum_components = isnothing(number_of_active_components) ?
                         length(peak_means) : Int(number_of_active_components)
    2 <= minimum_components <= length(peak_means) || throw(ArgumentError(
        "number_of_active_components must be between 2 and the number of candidate peaks"))

    best = nothing

    for subset_size in maximum_components:-1:minimum_components
        local_best = nothing
        candidate_count = length(peak_means)
        for mask in 0:((1 << candidate_count) - 1)
            count_ones(mask) == subset_size || continue
            (mask & (1 << (modal_peak_index - 1))) != 0 || continue
            selected_indices = [index for index in 1:candidate_count
                                if (mask & (1 << (index - 1))) != 0]
            selected_means = peak_means[selected_indices]
            selected_heights = peak_heights[selected_indices]
            selected_anchor = findfirst(==(modal_peak_index), selected_indices)
            anchor_mean = selected_means[selected_anchor]

            spacing_candidates = Float64[]
            for left in 1:(length(selected_means) - 1)
                for right in (left + 1):length(selected_means)
                    difference = selected_means[right] - selected_means[left]
                    for jump in 1:max_copy_jump
                        push!(spacing_candidates, difference / jump)
                    end
                end
            end

            for initial_spacing in spacing_candidates
                initial_spacing > 0 || continue
                if !isnothing(reference_spacing)
                    abs(initial_spacing - reference_spacing) /
                        reference_spacing <=
                        reference_spacing_tolerance_fraction || continue
                end
                relative_labels = round.(Int,
                    (selected_means .- anchor_mean) ./ initial_spacing)
                relative_labels[selected_anchor] == 0 || continue
                issorted(relative_labels) || continue
                length(unique(relative_labels)) == length(relative_labels) || continue
                all(1 .<= diff(relative_labels) .<= max_copy_jump) || continue
                copy_labels = modal_copy_number .+ relative_labels
                all(first(copy_number_range) .<= copy_labels .<=
                    last(copy_number_range)) || continue

                weights = selected_heights ./ max(sum(selected_heights), eps())
                x_bar = sum(weights .* relative_labels)
                y_bar = sum(weights .* selected_means)
                denominator = sum(weights .* (relative_labels .- x_bar).^2)
                denominator > eps() || continue
                fitted_spacing = sum(weights .* (relative_labels .- x_bar) .*
                                     (selected_means .- y_bar)) / denominator
                fitted_spacing > 0 || continue
                if !isnothing(reference_spacing)
                    abs(fitted_spacing - reference_spacing) /
                        reference_spacing <=
                        reference_spacing_tolerance_fraction || continue
                end
                fitted_anchor = y_bar - fitted_spacing * x_bar
                predicted = fitted_anchor .+ fitted_spacing .* relative_labels
                residuals = abs.(selected_means .- predicted) ./ fitted_spacing
                maximum(residuals) <= spacing_tolerance_fraction || continue
                missing_states = sum(max.(0, diff(copy_labels) .- 1))
                score = sum(weights .* residuals.^2) +
                        missing_state_penalty * missing_states
                proposal = (
                    score=score,
                    selected_indices=selected_indices,
                    observed_peak_means=selected_means,
                    observed_peak_heights=selected_heights,
                    modal_peak_observed_mean=anchor_mean,
                    modal_peak_copy_number=modal_copy_number,
                    copy_labels=Int.(copy_labels),
                    anchor_mean=fitted_anchor,
                    spacing=fitted_spacing,
                    means=fitted_anchor .+ fitted_spacing .* relative_labels,
                    maximum_spacing_residual=maximum(residuals),
                )
                if local_best === nothing || proposal.score < local_best.score
                    local_best = proposal
                end
            end
        end
        if local_best !== nothing
            best = local_best
            break
        end
    end

    best === nothing && throw(ArgumentError(
        "Candidate peaks cannot form an integer copy-number grid under the requested tolerance, copy-number range, and maximum jump"))
    return best
end

function _cn_state_posterior(observation::Real, parameter)
    isfinite(observation) || throw(ArgumentError("observation must be finite"))
    component_count = length(parameter.means)
    length(parameter.stds) == component_count &&
        length(parameter.weights) == component_count ||
        throw(DimensionMismatch("means, stds, and weights must have equal lengths"))
    cn_scores = parameter.weights .* pdf.(
        Normal.(parameter.means, max.(parameter.stds, eps())), observation)
    nearest_cn_z = minimum(
        abs(observation - parameter.means[k]) /
        max(parameter.stds[k], eps()) for k in 1:component_count)
    background_enabled = hasproperty(parameter, :background) &&
                         parameter.background.enabled
    background_score = background_enabled ?
        parameter.background.weight * pdf(
            Normal(parameter.background.mean,
                   max(parameter.background.std, eps())), observation) : 0.0
    total_score = sum(cn_scores) + background_score
    if !isfinite(total_score) || total_score <= 1e-300

        cn_probabilities = zeros(component_count)
        cn_probabilities[argmin(abs.(observation .- parameter.means))] = 1.0
        return (
            cn_probabilities=cn_probabilities,
            conditional_cn_probabilities=copy(cn_probabilities),
            background_probability=0.0,
            nearest_cn_z=nearest_cn_z,
            background_eligible=background_enabled,
            background_full_support=background_enabled,
        )
    end
    cn_probabilities = cn_scores ./ total_score
    background_probability = background_score / total_score
    cn_total = sum(cn_probabilities)
    conditional_cn_probabilities = cn_total > 0 ?
        cn_probabilities ./ cn_total : fill(1 / component_count, component_count)
    return (
        cn_probabilities=cn_probabilities,
        conditional_cn_probabilities=conditional_cn_probabilities,
        background_probability=background_probability,
        nearest_cn_z=nearest_cn_z,
        background_eligible=background_enabled,
        background_full_support=background_enabled,
    )
end

function _cn_constrained_em(
    data::Vector{Float64},
    lattice;
    shared_std::Bool,
    min_std_fraction::Real,
    max_std_fraction::Real,
    include_background::Bool,
    max_background_weight::Real,
    background_min_std_fraction::Real,
    background_max_std_fraction::Real,
    spacing_update_fraction::Real,
    max_iter::Int,
    tol::Real,
)
    labels = lattice.copy_labels
    relative_labels = labels .- lattice.modal_peak_copy_number
    component_count = length(labels)
    means = copy(lattice.means)
    initial_spacing = lattice.spacing
    spacing = initial_spacing
    anchor_mean = lattice.anchor_mean
    minimum_std = max(eps(), min_std_fraction * spacing)
    maximum_std = max(minimum_std, max_std_fraction * spacing)
    initial_std = clamp(std(data) / max(component_count, 1),
                        minimum_std, maximum_std)
    stds = fill(initial_std, component_count)
    background_weight = include_background ? min(0.02, max_background_weight) : 0.0
    cn_weight_total = 1 - background_weight
    peak_weights = lattice.observed_peak_heights
    weights = cn_weight_total .* peak_weights ./ sum(peak_weights)

    background_mean = median(data)
    background_std_lower = background_min_std_fraction * spacing
    background_std_upper = max(
        background_std_lower, background_max_std_fraction * spacing)
    background_std = clamp(std(data), background_std_lower,
                           background_std_upper)
    responsibilities = zeros(length(data), component_count)
    background_responsibility = include_background ? zeros(length(data)) : Float64[]
    converged = false
    iterations = 0

    for iteration in 1:max_iter
        iterations = iteration
        previous_means = copy(means)
        previous_stds = copy(stds)

        mixture_density = zeros(length(data))
        for component in 1:component_count
            responsibilities[:, component] .= weights[component] .* pdf.(
                Normal(means[component], max(stds[component], eps())), data)
            mixture_density .+= responsibilities[:, component]
        end
        if include_background
            background_responsibility .= background_weight .* pdf.(
                Normal(background_mean, max(background_std, eps())), data)
            mixture_density .+= background_responsibility
        end
        mixture_density .= max.(mixture_density, 1e-300)
        responsibilities ./= mixture_density
        include_background && (background_responsibility ./= mixture_density)

        effective_counts = vec(sum(responsibilities, dims=1))
        raw_centers = similar(means)
        for component in 1:component_count
            raw_centers[component] = effective_counts[component] > eps() ?
                sum(responsibilities[:, component] .* data) /
                effective_counts[component] : means[component]
        end

        regression_weights = effective_counts ./ max(sum(effective_counts), eps())
        x_bar = sum(regression_weights .* relative_labels)
        y_bar = sum(regression_weights .* raw_centers)
        denominator = sum(regression_weights .*
                          (relative_labels .- x_bar).^2)
        if denominator > eps()
            proposed_spacing = sum(regression_weights .*
                (relative_labels .- x_bar) .* (raw_centers .- y_bar)) /
                denominator
            spacing_lower = initial_spacing * (1 - spacing_update_fraction)
            spacing_upper = initial_spacing * (1 + spacing_update_fraction)
            spacing = clamp(proposed_spacing, spacing_lower, spacing_upper)
            anchor_mean = y_bar - spacing * x_bar
            means .= anchor_mean .+ spacing .* relative_labels
        end

        minimum_std = max(eps(), min_std_fraction * spacing)
        maximum_std = max(minimum_std, max_std_fraction * spacing)
        if shared_std
            variance_numerator = 0.0
            for component in 1:component_count
                variance_numerator += sum(responsibilities[:, component] .*
                    (data .- means[component]).^2)
            end
            shared_value = sqrt(variance_numerator /
                                max(sum(effective_counts), eps()))
            stds .= clamp(shared_value, minimum_std, maximum_std)
        else
            for component in 1:component_count
                component_std = sqrt(sum(responsibilities[:, component] .*
                    (data .- means[component]).^2) /
                    max(effective_counts[component], eps()))
                stds[component] = clamp(component_std,
                                        minimum_std, maximum_std)
            end
        end

        if include_background
            background_count = sum(background_responsibility)
            if background_count > eps()
                background_mean = sum(background_responsibility .* data) /
                                  background_count
                background_std = sqrt(sum(background_responsibility .*
                    (data .- background_mean).^2) / background_count)
                background_std_lower = background_min_std_fraction * spacing
                background_std_upper = max(
                    background_std_lower,
                    background_max_std_fraction * spacing)
                background_std = clamp(background_std, background_std_lower,
                                       background_std_upper)
            end
            minimum_background_weight = min(1e-6, max_background_weight)
            background_weight = clamp(background_count / length(data),
                                      minimum_background_weight,
                                      max_background_weight)
        end
        available_cn_weight = 1 - background_weight
        weights .= available_cn_weight .* effective_counts ./
                   max(sum(effective_counts), eps())

        change = max(maximum(abs.(means .- previous_means)),
                     maximum(abs.(stds .- previous_stds)))
        if change < tol
            converged = true
            break
        end
    end

    final_density = zeros(length(data))
    for component in 1:component_count
        final_density .+= weights[component] .* pdf.(
            Normal(means[component], max(stds[component], eps())), data)
    end
    if include_background
        final_density .+= background_weight .* pdf.(
            Normal(background_mean, max(background_std, eps())), data)
    end
    loglikelihood = sum(log.(max.(final_density, 1e-300)))

    all_copy_labels = collect(minimum(labels):maximum(labels))
    all_means = anchor_mean .+
        spacing .* (all_copy_labels .- lattice.modal_peak_copy_number)
    active = in.(all_copy_labels, Ref(labels))
    return (
        means=means,
        stds=stds,
        weights=weights,
        copy_labels=labels,
        all_copy_labels=all_copy_labels,
        all_means=all_means,
        active=active,
        modal_peak_copy_number=lattice.modal_peak_copy_number,
        modal_peak_mean=anchor_mean,
        spacing=spacing,
        background=(enabled=include_background,
                    mean=background_mean,
                    std=background_std,
                    weight=background_weight,
                    mode=:full_gaussian,
                    full_support=true,
                    min_std_fraction=Float64(background_min_std_fraction),
                    max_std_fraction=Float64(background_max_std_fraction)),
        loglikelihood=loglikelihood,
        converged=converged,
        iterations=iterations,
        maximum_initial_spacing_residual=
            lattice.maximum_spacing_residual,
        observed_peak_means=lattice.observed_peak_means,
        selected_peak_indices=lattice.selected_indices,
    )
end

# Fit the shared integer copy-number model.
function fitting_copy_number_gmm(
    metaDataDF::DataFrame;
    window_size::Integer=10,
    gmm_window_size::Union{Nothing,Integer}=nothing,
    modal_peak_copy_number=2,
    copy_number_range=0:10,
    number_of_active_components::Union{Nothing,Integer}=nothing,
    nbins::Integer=100,
    histogram_smoothing_bins::Integer=5,
    min_peak_prominence::Real=0.001,
    prominence_window_bins::Integer=20,
    min_peak_distance_bins::Integer=6,
    min_peak_mass::Real=0.001,
    max_components::Integer=8,
    max_copy_jump::Integer=3,
    spacing_tolerance_fraction::Real=0.20,
    missing_state_penalty::Real=0.02,
    share_modal_peak_position::Bool=false,
    share_copy_number_spacing::Bool=false,
    shared_modal_peak_max_distance_fraction::Real=0.50,
    shared_spacing_tolerance_fraction::Real=0.20,
    shared_std::Bool=true,
    min_std_fraction::Real=0.05,
    max_std_fraction::Real=0.45,
    include_background::Bool=true,
    max_background_weight::Real=0.10,
    background_min_std_fraction::Real=1.0,
    background_max_std_fraction::Real=5.0,
    spacing_update_fraction::Real=0.25,
    max_iter::Integer=400,
    tol::Real=1e-6,
)
    ncol(metaDataDF) >= 4 || throw(DimensionMismatch(
        "metaDataDF must contain three coordinate columns and at least one RD column"))
    window_size >= 1 || throw(ArgumentError("window_size must be a positive integer"))
    resolved_gmm_window_size = isnothing(gmm_window_size) ?
                               Int(window_size) : Int(gmm_window_size)
    resolved_gmm_window_size >= 1 || throw(ArgumentError(
        "gmm_window_size must be a positive integer or nothing"))
    nbins >= 20 || throw(ArgumentError("nbins must be at least 20"))
    histogram_smoothing_bins >= 1 || throw(ArgumentError(
        "histogram_smoothing_bins must be a positive integer"))
    0 <= min_peak_prominence <= 1 || throw(ArgumentError(
        "min_peak_prominence must be in [0, 1]"))
    0 <= min_peak_mass <= 1 || throw(ArgumentError(
        "min_peak_mass must be in [0, 1]"))
    0 <= spacing_tolerance_fraction < 1 || throw(ArgumentError(
        "spacing_tolerance_fraction must be in [0, 1)"))
    shared_modal_peak_max_distance_fraction > 0 || throw(ArgumentError(
        "shared_modal_peak_max_distance_fraction must be positive"))
    0 <= shared_spacing_tolerance_fraction < 1 || throw(ArgumentError(
        "shared_spacing_tolerance_fraction must be in [0, 1)"))
    0 < min_std_fraction <= max_std_fraction || throw(ArgumentError(
        "Standard-deviation fractions must satisfy 0 < min_std_fraction <= max_std_fraction"))
    0 <= max_background_weight < 1 || throw(ArgumentError(
        "max_background_weight must be in [0, 1)"))
    include_background && max_background_weight <= 0 && throw(ArgumentError(
        "max_background_weight must be positive when include_background=true"))
    0 < background_min_std_fraction <= background_max_std_fraction ||
        throw(ArgumentError(
            "Background standard-deviation fractions must satisfy 0 < background_min_std_fraction <= background_max_std_fraction"))
    0 <= spacing_update_fraction < 1 || throw(ArgumentError(
        "spacing_update_fraction must be in [0, 1)"))
    first(copy_number_range) <= last(copy_number_range) ||
        throw(ArgumentError("copy_number_range must be nonempty and increasing"))

    sample_names = String.(names(metaDataDF)[4:end])
    chromosomes = String.(metaDataDF[:, 1])
    chromosome_ranges = _cn_chromosome_ranges(chromosomes)
    depth_matrix = Matrix{Float64}(metaDataDF[:, 4:end])
    all(isfinite, depth_matrix) || throw(ArgumentError(
        "metaDataDF RD columns contain missing or non-finite values"))

    profiles = Any[]
    for sample_index in 1:length(sample_names)
        sample_name = sample_names[sample_index]
        peak_moving = _cn_window_means(depth_matrix[:, sample_index],
                                       Int(window_size), chromosome_ranges)
        isempty(peak_moving) && throw(ArgumentError(
            "Sample  $sample_name  has insufficient data for candidate-peak window means"))
        gmm_moving = resolved_gmm_window_size == Int(window_size) ?
                     peak_moving :
                     _cn_window_means(depth_matrix[:, sample_index],
                                      resolved_gmm_window_size,
                                      chromosome_ranges)
        isempty(gmm_moving) && throw(ArgumentError(
            "Sample  $sample_name  has insufficient data for GMM window means"))
        minimum_value, maximum_value = extrema(peak_moving)
        minimum_value < maximum_value || throw(ArgumentError(
            "Sample  $sample_name  has constant windowed RD values"))
        padding = 3 * (maximum_value - minimum_value) / nbins
        edges = collect(range(minimum_value - padding,
                              maximum_value + padding;
                              length=nbins + 1))
        centers = (edges[1:end-1] .+ edges[2:end]) ./ 2
        counts = Float64.(fit(Histogram, peak_moving, edges).weights)
        smoothed_counts = _cn_smooth_histogram(
            counts, Int(histogram_smoothing_bins))
        peak_indices = _cn_detect_histogram_peaks(
            smoothed_counts;
            min_peak_prominence=min_peak_prominence,
            prominence_window_bins=Int(prominence_window_bins),
            min_peak_distance_bins=Int(min_peak_distance_bins),
            min_peak_mass=min_peak_mass,
            max_components=Int(max_components),
        )
        length(peak_indices) >= 2 || throw(ArgumentError(
            "Sample  $sample_name  has fewer than two valid detected peaks"))
        peak_means = Float64.(centers[peak_indices])
        peak_heights = Float64.(smoothed_counts[peak_indices])
        modal_cn = _cn_resolve_modal_copy_number(
            modal_peak_copy_number, sample_index, sample_name)
        first(copy_number_range) <= modal_cn <= last(copy_number_range) ||
            throw(ArgumentError(
                "Sample  $sample_name  has highest peak CN=$modal_cn  outside copy_number_range"))
        push!(profiles, (
            sample_index=sample_index,
            sample_name=sample_name,
            peak_moving=peak_moving,
            gmm_moving=gmm_moving,
            peak_means=peak_means,
            peak_heights=peak_heights,
            own_modal_peak_index=argmax(peak_heights),
            modal_cn=modal_cn,
        ))
    end

    if share_modal_peak_position
        modal_labels = unique(profile.modal_cn for profile in profiles)
        length(modal_labels) == 1 || throw(ArgumentError(
            "When share_modal_peak_position=true, all samples must use the same modal_peak_copy_number"))
    end

    shared_modal_reference = share_modal_peak_position ? median(Float64[
        profile.peak_means[profile.own_modal_peak_index] for profile in profiles
    ]) : nothing

    anchor_indices = Int[]
    for profile in profiles
        if share_modal_peak_position
            push!(anchor_indices,
                  argmin(abs.(profile.peak_means .- shared_modal_reference)))
        else
            push!(anchor_indices, profile.own_modal_peak_index)
        end
    end

    preliminary_lattices = Any[]
    for (profile, anchor_index) in zip(profiles, anchor_indices)
        push!(preliminary_lattices, _cn_assign_integer_lattice(
            profile.peak_means,
            profile.peak_heights,
            profile.modal_cn;
            copy_number_range=copy_number_range,
            number_of_active_components=number_of_active_components,
            max_copy_jump=Int(max_copy_jump),
            spacing_tolerance_fraction=spacing_tolerance_fraction,
            missing_state_penalty=missing_state_penalty,
            modal_peak_index=anchor_index,
        ))
    end

    shared_spacing_reference = share_copy_number_spacing ?
        median(Float64[lattice.spacing for lattice in preliminary_lattices]) :
        nothing

    if share_modal_peak_position

        distance_scale = isnothing(shared_spacing_reference) ?
            median(Float64[lattice.spacing for lattice in preliminary_lattices]) :
            shared_spacing_reference
        for (profile, anchor_index) in zip(profiles, anchor_indices)
            distance = abs(profile.peak_means[anchor_index] - shared_modal_reference)
            distance <= shared_modal_peak_max_distance_fraction * distance_scale ||
                throw(ArgumentError(
                    "Sample  $(profile.sample_name)  has no candidate peak near the shared anchor $(round(shared_modal_reference; digits=3)) ; nearest-peak distance=$(round(distance; digits=3))"))
        end
    end

    results = Any[]
    for (profile, anchor_index, preliminary_lattice) in
        zip(profiles, anchor_indices, preliminary_lattices)
        lattice = share_copy_number_spacing ? _cn_assign_integer_lattice(
            profile.peak_means,
            profile.peak_heights,
            profile.modal_cn;
            copy_number_range=copy_number_range,
            number_of_active_components=number_of_active_components,
            max_copy_jump=Int(max_copy_jump),
            spacing_tolerance_fraction=spacing_tolerance_fraction,
            missing_state_penalty=missing_state_penalty,
            modal_peak_index=anchor_index,
            reference_spacing=shared_spacing_reference,
            reference_spacing_tolerance_fraction=
                shared_spacing_tolerance_fraction,
        ) : preliminary_lattice
        fitted = _cn_constrained_em(
            profile.gmm_moving,
            lattice;
            shared_std=shared_std,
            min_std_fraction=min_std_fraction,
            max_std_fraction=max_std_fraction,
            include_background=include_background,
            max_background_weight=max_background_weight,
            background_min_std_fraction=background_min_std_fraction,
            background_max_std_fraction=background_max_std_fraction,
            spacing_update_fraction=spacing_update_fraction,
            max_iter=Int(max_iter),
            tol=tol,
        )
        push!(results, merge(fitted, (
            sample_index=profile.sample_index,
            sample_name=profile.sample_name,
            window_size=Int(window_size),
            peak_detection_window_size=Int(window_size),
            gmm_window_size=resolved_gmm_window_size,
            candidate_peak_means=profile.peak_means,
            candidate_peak_heights=profile.peak_heights,
            own_highest_peak_mean=
                profile.peak_means[profile.own_modal_peak_index],
            selected_modal_anchor_peak_mean=profile.peak_means[anchor_index],
            share_modal_peak_position=share_modal_peak_position,
            shared_modal_peak_reference=shared_modal_reference,
            share_copy_number_spacing=share_copy_number_spacing,
            shared_spacing_reference=shared_spacing_reference,
            absolute_cn_anchor=share_modal_peak_position ?
                :shared_modal_peak_position : :sample_highest_peak,
        )))
    end
    return results
end

# Background-aware breakpoint scoring.
function rate_coincident_background(
    obs_data,
    givenParameters;
    min_background_probability::Real=0.50,
    return_details::Bool=false,
)
    0 <= min_background_probability <= 1 || throw(ArgumentError(
        "min_background_probability must be in [0, 1]"))
    observations = vec(Float64.(obs_data))
    length(observations) == 2 || throw(ArgumentError(
        "obs_data must contain exactly two flanking RD observations"))
    all(isfinite, observations) || throw(ArgumentError(
        "obs_data contains non-finite values"))
    component_count = length(givenParameters.means)
    length(givenParameters.stds) == component_count &&
        length(givenParameters.weights) == component_count ||
        throw(DimensionMismatch("means, stds, and weights must have equal lengths"))
    labels = hasproperty(givenParameters, :copy_labels) ?
             Int.(givenParameters.copy_labels) :
             collect(1:component_count)
    length(labels) == component_count || throw(DimensionMismatch(
        "copy_labels length does not match the number of integer-CN components"))

    background_enabled = hasproperty(givenParameters, :background) &&
                         givenParameters.background.enabled
    state_count = component_count + (background_enabled ? 1 : 0)
    posterior = zeros(state_count, 2)
    conditional_cn_posterior = zeros(component_count, 2)
    background_probabilities = zeros(2)
    nearest_cn_z = zeros(2)
    background_full_support = falses(2)
    dominant_states = Any[]
    for side in 1:2
        state_posterior = _cn_state_posterior(
            observations[side], givenParameters)
        posterior[1:component_count, side] .=
            state_posterior.cn_probabilities
        conditional_cn_posterior[:, side] .=
            state_posterior.conditional_cn_probabilities
        background_probabilities[side] =
            state_posterior.background_probability
        nearest_cn_z[side] = state_posterior.nearest_cn_z
        background_full_support[side] =
            state_posterior.background_full_support
        if background_enabled
            posterior[end, side] = state_posterior.background_probability
        end
        maximum_cn_probability, maximum_cn_index = findmax(
            state_posterior.cn_probabilities)
        background_is_dominant =
            background_enabled &&
            state_posterior.background_probability >= min_background_probability &&
            state_posterior.background_probability > maximum_cn_probability
        push!(dominant_states,
              background_is_dominant ? :background :
                  labels[maximum_cn_index])
    end

    left_background = dominant_states[1] === :background
    right_background = dominant_states[2] === :background
    coincidence_probability = if left_background && right_background
        1.0
    elseif xor(left_background, right_background)
        0.0
    else

        sum(conditional_cn_posterior[:, 1] .*
            conditional_cn_posterior[:, 2])
    end
    details = (
        probability=coincidence_probability,
        posterior=posterior,
        conditional_cn_posterior=conditional_cn_posterior,
        copy_labels=labels,
        background_enabled=background_enabled,
        background_probabilities=background_probabilities,
        min_background_probability=Float64(min_background_probability),
        dominant_states=dominant_states,
        nearest_cn_z=nearest_cn_z,
        background_full_support=background_full_support,
        background_eligible=background_full_support,
    )
    return return_details ? details : coincidence_probability
end

function run_length_nonzero(v)

    result = Matrix{Int}(undef, 0, 3)
    n = length(v)
    n == 0 && return result

    current_val = v[1]
    start_idx = 1

    for i in 2:n
        if v[i] != current_val

            if current_val != 0
                result = vcat(result, [current_val start_idx i-1])
            end
            current_val = v[i]
            start_idx = i
        end
    end

    if current_val != 0
        result = vcat(result, [current_val start_idx n])
    end

    return result
end
