# AbsCNV

<p align="center">
  <b>A Julia-based toolkit for copy number variation detection from sequencing depth data</b>
</p>

## Overview

**AbsCNV** estimates absolute integer copy numbers from aligned read-depth (RD) bins and detects copy-number differences between samples. The workflow fits an integer-CN constrained Gaussian mixture model (GMM) with an optional broad Gaussian background, constructs shared examination intervals, calls CN where supported, and performs pairwise comparisons.

Absolute CN labels depend on the supplied `modal_peak_copy_number` assumption. RD normalization alone does not establish ploidy or tumor purity. Check the modal-CN assignment against independent biological information when available.

### Key Features

- **Integer-CN constrained GMM**: automatic candidate-peak detection and EM fitting with component means on an integer-CN lattice.
- **Broad Gaussian background**: a full-support background component competes with integer-CN components over the entire RD domain.
- **Multi-sample alignment**: optional common reference for modal-peak selection and initial CN spacing, followed by sample-specific EM fits.
- **Background-aware inference**: explicit background/uncertain calls and skipped pairwise states.
- **Pairwise CN comparison**: compares shared absolute CN labels across samples.
- **Adjacent-region extension**: optionally promotes callable neighboring intervals using normalized RD differences.
- **Visualization and export**: CNV regions, representative CN pairs, skipped regions, and RD plots.

Pairwise comparison uses the reported probability threshold without a multiple-testing correction.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Step 1: Data Loading](#step-1-data-loading)
  - [Step 2: Filtering and Normalization](#step-2-filtering-and-normalization)
  - [Step 3: GMM Parameter Estimation](#step-3-gmm-parameter-estimation)
  - [Step 4: Genome Segmentation](#step-4-genome-segmentation)
  - [Step 5: Copy Number Determination](#step-5-copy-number-determination)
  - [Step 6: Pairwise Comparison](#step-6-pairwise-comparison)
  - [Step 7: Visualization](#step-7-visualization)
- [Input Format](#input-format)
- [Output Format](#output-format)
- [Parameters](#parameters)
- [Examples](#examples)
- [Citation](#citation)
- [Contact](#contact)

## Installation

### Requirements

- Julia 1.11 (tested with Julia 1.11.5).
- Packages: CSV, DataFrames, Plots, StatsBase, Distributions, DelimitedFiles, MultipleTesting, Statistics, Random, and JLD2 (for saved intermediate results).

### Setup

```bash
git clone https://github.com/AqduanCode/AbsCNV
cd AbsCNV
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

In Windows PowerShell, the same single-quoted `-e` argument works. Start Julia in the repository environment:

```bash
julia --project=.
```

```julia
include("MyFunctions.jl")
include("MyProcess.jl")
using JLD2
```

Load `MyFunctions.jl` first. `MyProcess.jl` tries to install its listed packages if loading fails; instantiate the project before including it.

### Project Structure

```text
AbsCNV/
├── README.md
├── MyFunctions.jl
├── MyProcess.jl
├── Project.toml
├── Manifest.toml
├── .gitignore
├── examples/
│   └── run_background_pipeline.jl
└── sample data/
    ├── Sample_A.csv
    ├── Sample_B.csv
    └── README.md
```

`Project.toml` declares the dependencies; `Manifest.toml` records the tested dependency resolution for Julia 1.11.5.

## Quick Start

Run the complete two-sample synthetic example from the repository root:

```bash
julia --project=. examples/run_background_pipeline.jl
```

The script creates `output/`, writes intermediate JLD2 files and CSV tables, and saves `CNVplot.png`. The data demonstrate execution and are not an accuracy benchmark.

To analyze your own aligned CSVs:

```bash
julia --project=. examples/run_background_pipeline.jl /path/to/data /path/to/results sample1.csv sample2.csv
```

Quote paths containing spaces. With an input directory and output directory but no filenames, the script loads all `.csv` files in the input directory, sorted by name. Select files explicitly if that directory contains other CSVs.

The core workflow is:

```julia
include("MyFunctions.jl")
include("MyProcess.jl")

combineData = process_genomic_data("sample data";
    file_names=["Sample_A.csv", "Sample_B.csv"])
metaDataDF = filtering_normalizing_seqDepth(combineData;
    initial_quantiles=(0.01, 0.99), normalized_limits=[10, 250])

cn_param = fitting_copy_number_gmm(metaDataDF;
    window_size=50, gmm_window_size=10,
    min_peak_prominence=0.05,
    modal_peak_copy_number=2, copy_number_range=1:10,
    share_modal_peak_position=true, share_copy_number_spacing=true,
    include_background=true, max_background_weight=0.001,
    background_min_std_fraction=2.0, background_max_std_fraction=5.0,
    shared_std=true)

segment_result = genome_segment_background(metaDataDF, cn_param, 0.001, 10, 10;
    min_background_probability=0.90)
CNVscaffold = segment_result.CNVscaffold
copynumber_result = copynumber_determination_background(
    metaDataDF, CNVscaffold, cn_param, 10;
    thresholdBFactor=2.0, max_background_probability=0.90)
comparison_result = copyNumber_comparison_background(
    metaDataDF, CNVscaffold, copynumber_result, cn_param, 10;
    pvalue_theSame=0.005)
filter_result = copyNumber_comparison_filter_fn_background(
    metaDataDF, CNVscaffold, copynumber_result, comparison_result;
    highRDthreshold=18, iterative=true, max_iterations=10)
plot_result = plotCNVregionMulti_RD_background(
    metaDataDF, CNVscaffold, filter_result, [1 2];
    copynumber_result=copynumber_result, layout=(1, 1),
    positive_cnv_content=0.80, chr="all", filtered_depth_diff_threshold=18)
CSV.write("CNVresult_raw.csv", plot_result.cnv_regions)
CSV.write("CNVresult_filter.csv", plot_result.filtered_cnv_regions)
CSV.write("background_regions.csv", plot_result.background_regions)
savefig(plot_result.plot, "CNVplot.png")
```

The example uses `min_peak_prominence=0.05` to suppress weak transition-related peaks in the synthetic data. Example settings differ from several function defaults documented below. Inspect peak diagnostics before applying them to biological data.

## Usage

### Step 1: Data Loading

```julia
input_dir = "/path/to/data"
file_names = ["sample1.csv", "sample2.csv", "sample3.csv"]
combineData = process_genomic_data(input_dir;
    file_names=file_names, file_extension=".csv")
```

| Parameter | Type / default | Description |
| --- | --- | --- |
| `input_directory` | `AbstractString` | Input directory |
| `file_names` | `Vector{String}`, empty | Selected files in sample order; empty means discover and sort matching files |
| `file_extension` | `AbstractString`, `".csv"` | File-discovery extension; contents are read as CSV |

Use the keyword API above. It checks required columns, equal bin counts, and identical coordinates in the same row order across samples. Full file stems become sample column names (e.g. `Sample_A.csv` becomes `Sample_A`), and only `chr1`–`chr22` are retained. The loader does not sort bins; provide chromosome-grouped, position-sorted inputs.

### Step 2: Filtering and Normalization

```julia
metaDataDF = filtering_normalizing_seqDepth(combineData;
    initial_quantiles=(0.01, 0.99), normalized_limits=[10, 250])
```

1. Sum raw RD across samples per bin and keep bins strictly inside the requested total-RD quantiles.
2. Divide each sample by its mean over retained bins and multiply by 100.
3. Keep bins for which every sample is strictly inside the normalized RD limits.

| Parameter | Function default | Description |
| --- | --- | --- |
| `initial_quantiles` | `(0.005, 0.995)` | Lower/upper quantiles of total raw RD |
| `normalized_quantiles` | `(0.005, 0.995)` | Derive lower limit from row minima and upper limit from row maxima |
| `normalized_limits` | `nothing` | Explicit two-value bounds override the quantile-derived limits |

Filtering removes rows jointly across samples. The mean is scaled before the final filtering pass, so the final retained sample mean need not equal exactly 100. All later scaffold indices refer to rows of this filtered `metaDataDF`.

### Step 3: GMM Parameter Estimation

```julia
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

for p in cn_param
    println((sample=p.sample_name, labels=p.copy_labels, means=p.means,
             spacing=p.spacing, background=p.background, converged=p.converged))
end
```

The function computes chromosome-wise sliding-window means, detects histogram peaks, assigns candidate peaks to an integer-CN lattice, then fits each sample. It requires at least two valid detected peaks per sample; a flat or effectively unimodal profile may fail.

Inspect `candidate_peak_means`, `candidate_peak_heights`, `copy_labels`, fitted means, and spacing together. Weak peaks arising from transition windows can affect lattice assignment. Tune peak-detection settings against the observed histogram and independent biological information; a fixed component count alone does not ensure correct CN assignments.

The fitted integer means satisfy `mean(CN) = modal_peak_mean + spacing * (CN - modal_peak_copy_number)`. Use `copy_labels` to map components to CN: the component index is not the absolute copy number, and labels can skip states. `copy_number_range` bounds candidate labels; it does not create an active Gaussian for every CN in that range.

`share_modal_peak_position=true` selects candidate anchors near the median of sample modal peaks and requires the same modal CN across samples. `share_copy_number_spacing=true` constrains initial lattice assignment around a median reference spacing. Each sample subsequently updates its own anchor and spacing in EM; final fitted means/spacings are not forced to be identical across samples.

When enabled, the background is one broad Gaussian with `mode=:full_gaussian`, not a separate CN state and not a tail-only component. Its weight is capped by `max_background_weight`; its standard deviation is bounded by the configured fractions of CN spacing. `sum(p.weights) + p.background.weight` is approximately 1.

```julia
# Optional diagnostic; the plot displays the chosen windowed RD histogram and fitted CN means.
metaData = Matrix(metaDataDF[:, 4:end])
moving_average_histogram(metaData[:, 1], cn_param[1].means, 50)
```

This diagnostic helper uses a continuous row-wise window; for multiple chromosomes, subset a chromosome before plotting to avoid windows across chromosome boundaries. The fitting function itself uses chromosome-wise windows.

### Step 4: Genome Segmentation

```julia
segment_result = genome_segment_background(metaDataDF, cn_param, 0.001, 10, 10;
    min_background_probability=0.90)
CNVscaffold = segment_result.CNVscaffold
```

| Parameter | Default / example | Description |
| --- | --- | --- |
| `breaks_threshold` | Required; example `0.001` | Break if left/right state coincidence score is below this value |
| `flanking_length` | Required; example `10` | Bins on each side, excluding the focal bin |
| `minimum_length` | Required; example `10` | Minimum retained `end_row - start_row` |
| `min_background_probability` | `0.50`; example `0.90` | Background posterior threshold for dominant-state assignment |

Background must also exceed every individual CN posterior to be the dominant state. Background–background has coincidence 1, background–CN has coincidence 0, and CN–CN uses the dot product of conditional CN posteriors. A physical coordinate span exceeding 100 kb across the flanking window also creates a break; this threshold is hard-coded in the implementation.

The returned NamedTuple includes `CNVscaffold`, `boundary_type`, `sample_breaks`, `break_support_count`, left/right background probabilities, and other diagnostics. The scaffold is a three-column integer matrix: start row, end row, and `end_row - start_row` (not base pairs, and not the inclusive bin count).

Intervals are constructed between successive runs of breakpoint support, with short intervals and cross-chromosome endpoints removed. They do not necessarily cover every retained bin or chromosome end. Check for an empty scaffold before proceeding.

### Step 5: Copy Number Determination

```julia
copynumber_result = copynumber_determination_background(
    metaDataDF, CNVscaffold, cn_param, 10;
    thresholdBFactor=2.0, max_background_probability=0.90)
CSV.write("copynumber_background.csv", copynumber_result.copy_number)
```

The function computes centered local RD means (using `2 * flanking_length + 1` bins), averages available values within each scaffold interval, and evaluates CN/background posteriors. A background posterior at or above the threshold yields an uncallable background status. Otherwise the largest conditional CN probability `p` must satisfy `p / (1 - p) > thresholdBFactor` for a CN call. Other intervals remain uncertain or invalid.

| Field | Shape / meaning |
| --- | --- |
| `copy_number` | DataFrame with coordinates, scaffold indices, and one CN column per sample |
| `copy_number_calls` | Samples × intervals; integer CN or `missing` |
| `callable`, `status` | Samples × intervals; valid-call mask and status symbols |
| `background_probability` | Samples × intervals; background posterior |
| `component_probabilities` | Samples × intervals × components; CN posterior including competition with background |
| `conditional_component_probabilities` | Same shape; probabilities conditioned on integer CN |
| `fragment_depth` | Intervals × samples; local-mean RD summary |

Missing calls are not CN0. The implementation excludes unavailable centered means at dataset edges; callers should inspect callable/status diagnostics.

### Step 6: Pairwise Comparison

```julia
comparison_result = copyNumber_comparison_background(
    metaDataDF, CNVscaffold, copynumber_result, cn_param, 10;
    pvalue_theSame=0.005)

filter_result = copyNumber_comparison_filter_fn_background(
    metaDataDF, CNVscaffold, copynumber_result, comparison_result;
    highRDthreshold=18, iterative=true, max_iterations=10)
```

The comparison returns `state`, `probability`, `callable`, `skip_reason`, and `common_copy_numbers`. Its arrays have dimensions **intervals × samples × samples**.

| State | Meaning |
| --- | --- |
| `-1` | Skipped / uncallable; not evidence of equal CN |
| `0` | No detected difference (including equal called CN) |
| `1` | Detected difference |

Equal callable CN labels return state 0 directly. For unequal calls, the function tests common absolute CN labels using weighted chi-square upper-tail probabilities and returns state 1 when the resulting score is below `pvalue_theSame`. If there are no shared modeled labels, the comparison is skipped. There is no FDR adjustment in this function.

The optional extension examines internal intervals with state 0. A callable interval is promoted to state 1 if a neighboring interval was positive in the preceding iteration and the absolute normalized RD difference exceeds `highRDthreshold`. Skipped intervals remain skipped. Endpoint intervals are not candidates for promotion. Inspect `promoted_by_rd`, `promotion_iteration`, and `original_state` to distinguish promoted from initial calls. Extension is a heuristic; reduced false negatives or improved accuracy is not guaranteed.

### Step 7: Visualization

```julia
sample_pairs = [1 2; 1 3; 2 3]  # Requires three loaded samples
plot_result = plotCNVregionMulti_RD_background(
    metaDataDF, CNVscaffold, filter_result, sample_pairs;
    copynumber_result=copynumber_result,
    layout=(3, 1), plotsize=(1200, 1350),
    positive_cnv_content=0.80, chr="all",
    filtered_depth_diff_threshold=18)

CSV.write("CNVresult_raw.csv", plot_result.cnv_regions)
CSV.write("CNVresult_filter.csv", plot_result.filtered_cnv_regions)
CSV.write("background_regions.csv", plot_result.background_regions)
savefig(plot_result.plot, "CNVplot.png")
```

For two samples use `[1 2]` and `layout=(1, 1)`. Sample indices follow the loaded file order. Pass the complete `copynumber_result` through the keyword argument; a CN matrix alone is insufficient.

| Parameter | Function default | Description |
| --- | --- | --- |
| `layout` | `(1, 1)` | Panel layout; must accommodate all selected pairs |
| `plotsize` | `(1000, 500)` | Figure size |
| `positive_cnv_content` | `0.8` | Minimum positive-fragment proportion during region construction |
| `chr` | `"all"` | Chromosome shown in the figure; `"all"` uses fragment indices on the x-axis |
| `xlims_pre`, `ylims_pre` | `[0 0]` | Plot limits; zero pair means automatic |
| `filtered_depth_diff_threshold` | `20` | Retain exported filtered regions with absolute RD difference at least this value |
| `merge_regions` | `false` | Enable additional distance/RD-based gap merging |
| `min_gap`, `max_gap` | `1_000_000`, `3_000_000` | Gap thresholds in base pairs for optional merging |
| `diff_threshold`, `gap_diff_threshold` | `15`, `10` | RD thresholds for optional merging |

With default `merge_regions=false`, skipped intervals act as region-construction barriers. Optional gap merging is a separate operation; do not assume it preserves those barriers. Gray spans mark skipped pairwise intervals, which can include background, uncertain CN, or no shared CN model.

The `chr` option controls visualization; exported region tables are computed for the supplied scaffold and selected pairs and are not restricted to that plotted chromosome. Subset the tables explicitly if required.

## Input Format

Each sample CSV must contain:

| Column | Name | Description |
| --- | --- | --- |
| 1 | `Chromosome` | Human autosome identifier `chr1`–`chr22` |
| 2 | `BinLowEdge` | Bin start position |
| 3 | `BinUpEdge` | Bin end position |
| 4 | `Content` | Finite numeric read depth |

```csv
Chromosome,BinLowEdge,BinUpEdge,Content
chr1,0,1000,100.0
chr1,1000,2000,102.5
chr1,2000,3000,98.7
```

Files must use identical bins in identical order, a consistent coordinate convention, and unique file stems. The loader does not convert coordinate conventions. The bundled synthetic inputs use zero-based starts and exclusive upper edges. Avoid sample names `chr`, `start`, `end`, `scaffoldstart`, and `scaffoldend`, which conflict with CN output columns.

## Output Format

### Complete Example Outputs

| File | Content |
| --- | --- |
| `metaDataDF.csv` | Filtered normalized RD matrix with three coordinate columns |
| `sample_index.csv` | Mapping of one-based sample indices to sample names |
| `cn_param.jld2` | Fitted CN/background models and fit diagnostics |
| `CNVscaffold_background.csv` | Headered scaffold table: `start_row`, `end_row`, `row_span` |
| `segment_result.jld2` | Scaffold and breakpoint diagnostics |
| `copynumber_background.csv` | Per-interval CN calls by sample |
| `copynumber_result.jld2` | CN calls, masks, statuses, posteriors, and RD summaries |
| `comparison_result.jld2` | Initial pairwise states and test diagnostics |
| `filter_result.jld2` | Extended states and promotion diagnostics |
| `CNVresult_raw.csv` | Constructed CNV regions before the final RD table filter |
| `CNVresult_filter.csv` | Regions after the absolute RD-difference table filter |
| `background_regions.csv` | Consecutive skipped pairwise intervals |
| `CNVplot.png` | Saved RD/CNV plot |

The example writes the scaffold using `CSV.write` with headers. Scaffold coordinates are row indices in filtered `metaDataDF`.

### Copy Number DataFrame

| Column | Meaning |
| --- | --- |
| `chr` | Chromosome |
| `start` | `BinLowEdge` of the interval's first retained row |
| `end` | `BinLowEdge` of the interval's last retained row, as implemented |
| `scaffoldstart`, `scaffoldend` | Inclusive row indices in filtered `metaDataDF` |
| Sample-name columns | Integer CN, or an empty CSV field for `missing` |

**Coordinate detail:** the CN table's `end` uses the last bin's lower edge, whereas CNV region tables use the last bin's `BinUpEdge`. Do not interpret these two end columns interchangeably. For a full last-bin endpoint in a derived table, use `metaDataDF.BinUpEdge[table.scaffoldend]` explicitly. Use the documented coordinate definitions when exporting to other genomic formats.

### CNV Regions DataFrame

| Column | Description |
| --- | --- |
| `pair_sampleA`, `pair_sampleB` | One-based sample indices |
| `pair_sampleA_name`, `pair_sampleB_name` | Sample names |
| `chr`, `start_pos`, `end_pos` | Chromosome, first-bin lower edge, last-bin upper edge |
| `length` | `end_pos - start_pos` |
| `start_pos_metaData`, `end_pos_metaData` | Row indices in filtered `metaDataDF` |
| `mean_depth_A`, `mean_depth_B` | Mean normalized RD within the region |
| `depth_diff` | `mean_depth_B - mean_depth_A` |
| `copy_number_A`, `copy_number_B` | Representative differing CN pair |

The representative CN pair is chosen by overlap bin count among callable intervals that had an initial positive comparison and different CN calls. When extension output is passed, this uses its `original_state`, not RD-promoted intervals. Region rows lacking such a differing CN pair are removed. These values do not imply every bin in a region has the same CN pair. Empty region tables may contain only the CN columns because of the underlying implementation; check `nrow` before assuming a full column schema.

Despite its filename, `background_regions.csv` contains all selected pairwise state `-1` runs, including non-background skip reasons. Its fields are sample indices/names, `chr`, `start_pos`, `end_pos`, `start_fragment`, `end_fragment`, and `fragment_count`. Inspect `comparison_result.skip_reason` and `copynumber_result.status` for causes.

```julia
using JLD2
cn_param = JLD2.load("output/cn_param.jld2", "cn_param")
comparison_result = JLD2.load("output/comparison_result.jld2", "comparison_result")
```

## Parameters

### Main Model Parameters

| Parameter | Function default | Complete example | Meaning |
| --- | --- | --- | --- |
| `window_size` | `10` | `50` | Sliding-window width for candidate-peak detection, in retained bins |
| `gmm_window_size` | `nothing` | `10` | EM fitting window; `nothing` reuses `window_size` |
| `min_peak_prominence` | `0.001` | `0.05` | Minimum local histogram prominence relative to the highest peak |
| `modal_peak_copy_number` | `2` | `2` | Assumed modal CN; scalar, per-sample vector, or dictionary keyed by sample name/symbol/index |
| `copy_number_range` | `0:10` | `1:10` | Allowed integer labels; example excludes CN0 |
| `number_of_active_components` | `nothing` | `nothing` | Optional required number of active CN components, excluding background |
| `share_modal_peak_position` | `false` | `true` | Use a common modal-peak reference in anchor selection |
| `share_copy_number_spacing` | `false` | `true` | Use a common initial-spacing reference |
| `shared_modal_peak_max_distance_fraction` | `0.50` | `0.50` | Maximum anchor distance in units of reference spacing |
| `shared_spacing_tolerance_fraction` | `0.20` | `0.20` | Initial lattice tolerance around shared spacing |
| `shared_std` | `true` | `true` | Equal CN-component SDs within a sample |
| `min_std_fraction`, `max_std_fraction` | `0.05`, `0.45` | Same | CN SD bounds relative to CN spacing |
| `include_background` | `true` | `true` | Add the broad Gaussian background |
| `max_background_weight` | `0.10` | `0.001` | Maximum background mixture weight |
| `background_min_std_fraction` | `1.0` | `2.0` | Minimum background SD / CN spacing |
| `background_max_std_fraction` | `5.0` | `5.0` | Maximum background SD / CN spacing |
| `spacing_update_fraction` | `0.25` | Same | Allowed EM spacing deviation from its initialization |
| `max_iter`, `tol` | `400`, `1e-6` | Same | EM iteration limit and mean/SD-change convergence tolerance |

Additional peak/lattice controls are `nbins=100`, `histogram_smoothing_bins=5`, `prominence_window_bins=20`, `min_peak_distance_bins=6`, `min_peak_mass=0.001`, `max_components=8`, `max_copy_jump=3`, `spacing_tolerance_fraction=0.20`, and `missing_state_penalty=0.02`.

### Downstream Thresholds

| Parameter | Function default | Complete example |
| --- | --- | --- |
| `breaks_threshold` | Required positional argument | `0.001` |
| `flanking_length`, `minimum_length` for segmentation | Required positional arguments | `10`, `10` |
| `min_background_probability` for segmentation | `0.50` | `0.90` |
| `thresholdBFactor` for CN calls | `2.0` | `2.0` |
| `max_background_probability` for CN calls | `0.50` | `0.90` |
| `pvalue_theSame` for pairwise comparison | `0.05` | `0.005` |
| `highRDthreshold` for extension | `100` | `18` |
| `iterative` | `true` | `true` |
| `max_iterations` | `nothing` (up to interval count) | `10` |
| `filtered_depth_diff_threshold` for exported regions | `20` | `18` |

Window lengths count retained rows, not a fixed genomic distance after filtering. Increase/decrease settings based on RD distributions, bin width, callable coverage, and independent validation. The background weight cap and posterior cutoffs control different quantities. A lower mixture-weight cap does not imply every observation has a low background posterior.

The supplied segmentation/CN-call routines compute flanking means over the combined row matrix; those means are not explicitly reset at chromosome or filtered-coordinate gaps, although cross-chromosome scaffold endpoints are removed. Review boundary-adjacent results for multi-chromosome data. The single-chromosome demo does not establish correctness of boundary behavior.

## Examples

### Example 1: Basic CNV Detection

```bash
julia --project=. examples/run_background_pipeline.jl
```

This runs two synthetic samples with CN-like depth blocks and all stages of the background-aware workflow. Read [the sample-data description](sample%20data/README.md) before interpreting results.

### Example 2: Specific Chromosome Analysis

After running Steps 1–6 interactively:

```julia
chr1_result = plotCNVregionMulti_RD_background(
    metaDataDF, CNVscaffold, filter_result, [1 2];
    copynumber_result=copynumber_result, layout=(1, 1),
    chr="chr1", xlims_pre=[1e7 1.3e7], positive_cnv_content=0.8)
chr1_table = filter(:chr => ==("chr1"), chr1_result.cnv_regions)
CSV.write("chr1_cnv_regions.csv", chr1_table)
```

For a single sample, Steps 1–5 provide CN calls when model fitting succeeds. Pairwise comparison, adjacent-pair extension, and pairwise plotting require multiple samples.

## Citation

For reproducible reporting, cite the [AbsCNV repository](https://github.com/AqduanCode/AbsCNV) and record the exact commit used for the analysis.

## Contact

For questions and bug reports, open a [GitHub issue](https://github.com/AqduanCode/AbsCNV/issues) or contact Helab at duananqi@fudan.edu.cn.
