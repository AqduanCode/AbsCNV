using CSV
using DataFrames
using DelimitedFiles
using Statistics
using Distributions

if length(ARGS) < 1
    error("Usage: julia generate_workflow_intermediates.jl <project_root>")
end

root = abspath(ARGS[1])
# The workflow illustration uses its original algorithm snapshot.
# Run in a separate Julia process from the current benchmark entries.
include(joinpath(@__DIR__, "workflow_functions_snapshot.jl"))
include(joinpath(@__DIR__, "workflow_process_snapshot.jl"))

sim_dir = joinpath(root, "00_simulation_data")
norm_dir = joinpath(root, "01_Fig1b_RD_normalization_smoothing")
gmm_dir = joinpath(root, "02_Fig1c_sample_specific_GMM")
bp_dir = joinpath(root, "03_Fig1d_breakpoints_shared_scaffold")
state_dir = joinpath(root, "04_Fig1e_posterior_state_calls")
pair_dir = joinpath(root, "05_Fig1f_pairwise_refinement")
summary_dir = joinpath(root, "06_summary")
for d in [norm_dir, gmm_dir, bp_dir, state_dir, pair_dir, summary_dir]
    mkpath(d)
end

samples = ["Sample_A", "Sample_B", "Sample_C"]
raw_df = CSV.read(joinpath(sim_dir, "simulated_raw_RD_20Mb_3samples.csv"), DataFrame)
truth_df = CSV.read(joinpath(sim_dir, "simulation_truth_bin_level.csv.gz"), DataFrame)

# Run the original AbsCNV filtering and per-sample normalization.
meta_df = filtering_normalizing_seqDepth(raw_df, [100, 2500], [15, 250])
CSV.write(joinpath(norm_dir, "normalized_RD_after_AbsCNV_filtering.csv.gz"), meta_df; compress=true)

meta = Matrix{Float64}(meta_df[:, 4:end])
n_bins_retained, nsample = size(meta)
window_size = 10
smooth10 = fill(NaN, n_bins_retained, nsample)
# This is exactly the forward 10-bin mean constructed inside fitting_gmm:
# window i contains rows i:(i + window_size - 1).
for i in 1:(n_bins_retained - window_size + 1), j in 1:nsample
    smooth10[i, j] = mean(meta[i:(i + window_size - 1), j])
end
smooth_df = select(meta_df, 1:3)
for (j, s) in enumerate(samples)
    smooth_df[!, Symbol(s * "_normalized_rd")] = meta[:, j]
    smooth_df[!, Symbol(s * "_smooth10_rd")] = smooth10[:, j]
end
CSV.write(joinpath(norm_dir, "normalized_and_smooth10_RD.csv.gz"), smooth_df; compress=true)

# Use the actual AbsCNV fitting_gmm function; four states are represented in this controlled simulation.
pre_mean = [50.0, 100.0, 150.0, 200.0]
esti_param = fitting_gmm(meta, pre_mean, window_size)
if length(esti_param) != nsample
    error("fitting_gmm returned $(length(esti_param)) parameter sets for $nsample samples")
end

gmm_records = DataFrame(sample=String[], state=Int[], mean=Float64[], std=Float64[], weight=Float64[])
for j in 1:nsample
    ord = sortperm(esti_param[j].means)
    means = esti_param[j].means[ord]
    stds = esti_param[j].stds[ord]
    weights = esti_param[j].weights[ord]
    esti_param[j] = (means=means, stds=stds, weights=weights)
    for k in eachindex(means)
        push!(gmm_records, (samples[j], k, means[k], stds[k], weights[k]))
    end
end
CSV.write(joinpath(gmm_dir, "gmm_parameters_from_fitting_gmm.csv"), gmm_records)

function posterior_matrix(values::Vector{Float64}, param)
    k = length(param.means)
    out = zeros(length(values), k)
    for i in eachindex(values)
        dens = [param.weights[z] * pdf(Normal(param.means[z], max(param.stds[z], 1e-6)), values[i]) for z in 1:k]
        den = sum(dens)
        out[i, :] .= den > 0 ? dens ./ den : fill(1 / k, k)
    end
    out
end

# Export left/right means and the same-state probability used by AbsCNV breakpoint detection.
flank = 10
bp_scores = select(meta_df, 1:3)
for j in 1:nsample
    left = fill(NaN, n_bins_retained)
    right = fill(NaN, n_bins_retained)
    coincident = fill(NaN, n_bins_retained)
    candidate = zeros(Int, n_bins_retained)
    for i in (flank + 1):(n_bins_retained - flank)
        left[i] = mean(meta[(i - flank):(i - 1), j])
        right[i] = mean(meta[(i + 1):(i + flank), j])
        coincident[i] = rate_coincident([left[i], right[i]], esti_param[j])
        candidate[i] = coincident[i] < 0.001 ? 1 : 0
    end
    bp_scores[!, Symbol(samples[j] * "_left_mean")] = left
    bp_scores[!, Symbol(samples[j] * "_right_mean")] = right
    bp_scores[!, Symbol(samples[j] * "_p_same")] = coincident
    bp_scores[!, Symbol(samples[j] * "_score")] = -log10.(clamp.(coincident, 1e-300, 1.0))
    bp_scores[!, Symbol(samples[j] * "_candidate")] = candidate
end
CSV.write(joinpath(bp_dir, "breakpoint_scores_from_AbsCNV.csv.gz"), bp_scores; compress=true)

# Reproduce and export the exact multi-sample breakpoint runs used by
# genome_segment. run_length_nonzero separates adjacent runs when the number
# of supporting samples changes, exactly as the original function does.
candidate_sum = zeros(Int, n_bins_retained)
for s in samples
    candidate_sum .+= bp_scores[!, Symbol(s * "_candidate")]
end
merged_runs = run_length_nonzero(candidate_sum)
merged_runs_df = DataFrame(
    run_id = 1:size(merged_runs, 1),
    supporting_sample_count = merged_runs[:, 1],
    start_row_1based = merged_runs[:, 2],
    end_row_1based = merged_runs[:, 3],
    start_bp = Int.(meta_df[merged_runs[:, 2], 2]),
    end_bp = Int.(meta_df[merged_runs[:, 3], 3])
)
CSV.write(joinpath(bp_dir, "merged_breakpoint_runs_from_genome_segment_logic.csv"), merged_runs_df)

# Run the actual AbsCNV shared-scaffold function.
scaffold = genome_segment(meta_df, esti_param, 0.001, 10, 10)
scaffold_df = DataFrame(
    segment_id = 1:size(scaffold, 1),
    start_row_1based = scaffold[:, 1],
    end_row_1based = scaffold[:, 2],
    length_bins = scaffold[:, 3],
    chr = string.(meta_df[scaffold[:, 1], 1]),
    start_bp = Int.(meta_df[scaffold[:, 1], 2]),
    end_bp = Int.(meta_df[scaffold[:, 2], 3])
)
CSV.write(joinpath(bp_dir, "shared_scaffold_from_genome_segment.csv"), scaffold_df)

# Audit the key invariant: no candidate breakpoint row may occur strictly
# inside a scaffold. Candidate rows at the two endpoints are boundaries and
# are therefore excluded from the internal count.
audit_df = DataFrame(
    segment_id = scaffold_df.segment_id,
    start_row_1based = scaffold_df.start_row_1based,
    end_row_1based = scaffold_df.end_row_1based,
    internal_candidate_rows = zeros(Int, nrow(scaffold_df))
)
for i in 1:nrow(scaffold_df)
    lo = scaffold_df.start_row_1based[i]
    hi = scaffold_df.end_row_1based[i]
    audit_df.internal_candidate_rows[i] =
        count(candidate_sum[(lo + 1):(hi - 1)] .> 0)
end
CSV.write(joinpath(bp_dir, "scaffold_internal_breakpoint_audit.csv"), audit_df)
if any(audit_df.internal_candidate_rows .> 0)
    error("A candidate breakpoint was found strictly inside an AbsCNV scaffold")
end

# Run actual AbsCNV interval-level posterior-odds state calling.
cn_info = copynumber_determination(meta_df, scaffold, esti_param, 10; thresholdBFactor=2.0)
call_df = copy(scaffold_df)
posterior_long = DataFrame(segment_id=Int[], sample=String[], state=Int[], posterior_probability=Float64[])
for j in 1:nsample
    call_df[!, Symbol(samples[j] * "_called_state")] = Int.(round.(cn_info[j, :, 1]))
    call_df[!, Symbol(samples[j] * "_max_posterior")] = vec(maximum(cn_info[j, :, 2:end], dims=2))
    for i in 1:size(scaffold, 1), k in 1:length(esti_param[j].means)
        push!(posterior_long, (i, samples[j], k, cn_info[j, i, k + 1]))
    end
end

# Add truth mode after coordinate matching.
for s in samples
    truth_col = Symbol(s * "_truth_effective_state")
    modes = Float64[]
    for i in 1:nrow(call_df)
        mask = (truth_df.chr .== call_df.chr[i]) .&
               (truth_df.binLowEdge .>= call_df.start_bp[i]) .&
               (truth_df.binUpEdge .<= call_df.end_bp[i])
        vals = truth_df[mask, truth_col]
        if isempty(vals)
            push!(modes, NaN)
        else
            counts = DataFrames.combine(groupby(DataFrame(v=vals), :v), DataFrames.nrow => :n)
            push!(modes, counts.v[argmax(counts.n)])
        end
    end
    call_df[!, Symbol(s * "_truth_mode")] = modes
end
CSV.write(joinpath(state_dir, "shared_segments_and_state_calls.csv"), call_df)
CSV.write(joinpath(state_dir, "state_posteriors_long.csv"), posterior_long)

# Run actual AbsCNV pairwise comparison and RD-neighbour refinement.
paired = copyNumber_comparison(meta_df, scaffold, cn_info, esti_param, 10;
    pvalue_theSame=0.05e-3, fdr_threshold=0.1e-3)
paired_refined = copyNumber_comparison_filter_fn(meta_df, scaffold, cn_info, paired;
    highRDthreshold=20)

pair_tables = DataFrame[]
for (label, a, b) in [("A_vs_B", 1, 2), ("A_vs_C", 1, 3), ("B_vs_C", 2, 3)]
    temp = select(scaffold_df, [:segment_id, :chr, :start_bp, :end_bp, :length_bins])
    temp.pair = fill(label, nrow(temp))
    temp.sample_1 = fill(samples[a], nrow(temp))
    temp.sample_2 = fill(samples[b], nrow(temp))
    temp.sample_1_state = Int.(round.(cn_info[a, :, 1]))
    temp.sample_2_state = Int.(round.(cn_info[b, :, 1]))
    temp.initial_positive = Int.(round.(paired[:, a, b]))
    temp.refined_positive = Int.(round.(paired_refined[:, a, b]))
    push!(pair_tables, temp)
end
pair_records = vcat(pair_tables...; cols=:union)
CSV.write(joinpath(pair_dir, "pairwise_segment_calls_from_AbsCNV.csv"), pair_records)

summary = DataFrame(
    metric = ["raw_bins", "retained_bins", "shared_segments", "no_calls", "pairwise_initial_positive", "pairwise_refined_positive"],
    value = [
        nrow(raw_df),
        nrow(meta_df),
        nrow(scaffold_df),
        sum([sum(call_df[!, Symbol(s * "_called_state")] .== -1) for s in samples]),
        sum(pair_records.initial_positive),
        sum(pair_records.refined_positive)
    ]
)
CSV.write(joinpath(summary_dir, "run_summary.csv"), summary)

println("AbsCNV analysis completed: $(nrow(meta_df)) retained bins, $(nrow(scaffold_df)) shared segments")
