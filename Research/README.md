# AbsCNV manuscript research code

Copy **only this `Research/` folder** into the existing AbsCNV repository root, beside `MyFunctions.jl`, `MyProcess.jl` and `Project.toml`. Existing software, examples and sample data are not duplicated in this add-on. `README_zh.md` is the Chinese version.

The companion dataset provides the inputs and saved results used by these analyses. Dataset paths below refer to its `datasets/` directory. Obtain the code version associated with that dataset and keep newly generated results in separate directories.

## Requirements and data preparation

Use Julia 1.11 (the supplied Manifest is resolved with 1.11.5), Python 3.9+ and R with `jsonlite`, `data.table`, `ggplot2`, `scales`, `systemfonts`, `ragg` and `svglite`. Install the Python requirements below. Arial is used for paper plots; install it or supply `--font-dir` where supported. Illustrator is needed only for the sampling-scheme JSX. External CNV tools, ART, alignment and ROOT tools are needed only to rerun raw-sequencing stages.

The archive stores identical RD inputs once. Single/paired A share one stored file, and the six length A inputs share another; the two groups differ and are never merged. **Run `prepare_working_data.py` before inference**: it creates a new working copy, decompresses text inputs and restores all eight logical input aliases by hard links or copies. The archive remains unchanged. Original compressed plotting snapshots and VCFs are retained. Use new output directories.

## Analyses in manuscript order

Paths in the table are relative to the prepared `datasets/` directory. The small `read_single_sample_predictions.py` and `paired_interval_metrics.py` modules are required dependencies of the evaluators.

| Paper component | Purpose-named entry points | Archived inputs | Expected outputs |
|---|---|---|---|
| Figure 1 | `plot_method_workflow.R` | `Figure1_AbsCNV workflow/` (seven saved illustrative tables) | Method panels b–f |
| Single benchmark | `call_single_sample_copy_numbers.jl`, `standardize_single_sample_predictions.py`, `evaluate_single_sample_predictions.py` | `single_sample/input`, `benchmark`, `native`, `output/copynumber_background.csv`, `othertools` | Single-sample CN calls; five unmerged tables; four-standard summary and optional match details |
| Paired benchmark | `call_paired_sample_differences.jl`, `standardize_paired_sample_predictions.py`, `evaluate_paired_sample_predictions.py` | `paired_samples/input`, `benchmark`, `raw_for_normalization`, `normalized_calls`, `output` | Paired CNV calls; native-CN-priority normalized predictions; four-standard summary |
| Figure 2 | `plot_benchmark_performance.py` | Single/paired `evaluation/*four_standard*summary.csv` | 175 × 113 mm PDF/PNG, 20 plotted metric rows |
| Six lengths / Figure S1 | `call_length_benchmark_differences.jl`, `evaluate_and_plot_length_performance.R` | `length_series/input`, `benchmark`, `output` | Six CNV tables; event/base metrics; F1-by-length figure |
| Public WGS / Figure 3 | `call_public_wgs_differences.jl`, `export_public_copy_number_peaks.jl`, `plot_public_copy_number_histograms.R`, `plot_public_pairwise_copy_number_profiles.R` | `public_wgs/input`, `output`, `histogram_current/cn_peaks.csv` | Three pairwise call sets; latest saved-GMM histograms; pairwise profiles |
| LUMPY / Table S5 | `evaluate_lumpy_support.py` | `public_wgs/LUMPY_SV`, `output/CNVresult_filter.csv` | Final directional support counts and proportions |
| U251 / Figure 4 | `call_clonal_lineage_differences.jl`, `plot_lineage_chromosome_events.py`, `evaluate_lineage_event_inheritance.py`, `summarize_lineage_specific_differences.py` | Ten `u251/input` CSVs; `u251/output` calls/CN background | Chromosome matrices and positional enlargements; inherited/same/opposite fractions; spatially distinct lengths |
| Figure S3 | `export_lineage_copy_number_peaks.jl`, `plot_lineage_copy_number_histograms.R` | U251 saved GMM and normalized RD; `peak_parameters` | Ten-sample histogram PDF |
| Figures S2 / S4 | `draw_clonal_sampling_scheme.jsx`, `assemble_supplementary_lineage_panels.py` | Fixed sampling scheme; generated U12/U13 PNGs | Editable sampling diagram; supplementary Word assembly |

## Recompute metrics and figures from the archived results

Replace `SOURCE_DATASETS`, `DATASETS` and `RESULTS` with actual paths. Commands are run from the AbsCNV repository root; create the RESULTS root before shell output redirection. The figure command uses the archived full-precision summaries. To plot newly recalculated summaries instead, provide `--single-metrics RESULTS/single/single_sample_four_standard_summary.csv --paired-metrics RESULTS/paired/four_standard_f1_summary.csv`.

```text
# Install from the AbsCNV repository root.
julia --project=. -e "using Pkg; Pkg.instantiate()"
python -m pip install -r Research/requirements.txt

# SOURCE_DATASETS = the extracted Zenodo/datasets directory.
python Research/prepare_working_data.py --source "SOURCE_DATASETS" --output "NEW_WORK_DATASETS"

# In the following commands DATASETS is NEW_WORK_DATASETS, RESULTS is a new output folder.
python Research/evaluate_single_sample_predictions.py --base "DATASETS/single_sample" --output-dir "RESULTS/single"
python Research/evaluate_paired_sample_predictions.py --dataset "DATASETS/paired_samples" --output-dir "RESULTS/paired"
python Research/plot_benchmark_performance.py --data-root "DATASETS" --output-dir "RESULTS/figure2"
Rscript Research/evaluate_and_plot_length_performance.R "DATASETS/length_series" "RESULTS/lengths"

Rscript Research/plot_method_workflow.R "DATASETS/Figure1_AbsCNV workflow" "RESULTS/method_panels"
julia --project=. Research/export_public_copy_number_peaks.jl "DATASETS/public_wgs/output/cn_param.jld2" "RESULTS/public_hist/cn_peaks.csv"
Rscript Research/plot_public_copy_number_histograms.R "DATASETS/public_wgs/output/metaDataDF.csv" "RESULTS/public_hist/cn_peaks.csv" "RESULTS/public_hist"
Rscript Research/plot_public_pairwise_copy_number_profiles.R "DATASETS/public_wgs/output" "RESULTS/public_pairs"
python Research/evaluate_lumpy_support.py "DATASETS/public_wgs" > "RESULTS/lumpy_support.json"

python Research/plot_lineage_chromosome_events.py --cnv "DATASETS/u251/output/CNVresult_filter.csv" --background "DATASETS/u251/output/copynumber_background.csv" --output "RESULTS/lineage_events"
python Research/evaluate_lineage_event_inheritance.py --input "DATASETS/u251/output" --output "RESULTS/inheritance"
python Research/summarize_lineage_specific_differences.py --input "DATASETS/u251/output" --output "RESULTS/lineage_specific"
julia --project=. Research/export_lineage_copy_number_peaks.jl "DATASETS/u251/output/cn_param.jld2" "RESULTS/u251_hist/cn_peak_parameters_background_updated.csv"
python Research/assemble_supplementary_lineage_panels.py --input-dir "RESULTS/lineage_events" --output "RESULTS/supplementary_lineages.docx"
```

For the ten-sample histogram set `CNV_PLOT_PROJECT_DIR=DATASETS/u251` and `CNV_PLOT_OUTPUT_DIR=RESULTS/u251_hist`, then run `Rscript Research/plot_lineage_copy_number_histograms.R`. The peak CSV must already exist in that output folder. In Bash use `export NAME=value`; in PowerShell use `$env:NAME="value"`. Run `draw_clonal_sampling_scheme.jsx` from Illustrator and choose its output folder. The supplementary assembly script produces DOCX; the archived PDF was exported from Word, not generated by that Python script.

## Rerun AbsCNV inference from RD

Inference parameters and computational calls are unchanged; paths and filenames are purpose-based. Run each branch into a new directory. Repeat the last command for 50kb, 100kb, 200kb, 500kb and 1000kb. The length output root gets an `output/LABEL/` tree.

```text
julia --project=. Research/call_single_sample_copy_numbers.jl "DATASETS/single_sample/input" "RESULTS/single_inference"
julia --project=. Research/call_paired_sample_differences.jl "DATASETS/paired_samples/input" "RESULTS/paired_inference"
julia --project=. Research/call_public_wgs_differences.jl "DATASETS/public_wgs/input" "RESULTS/public_inference"
julia --project=. Research/call_clonal_lineage_differences.jl "DATASETS/u251/input" "RESULTS/u251_inference"
julia --project=. Research/call_length_benchmark_differences.jl 20kb "DATASETS/length_series" "RESULTS/length_inference"
```

After new single-sample inference, point `standardize_single_sample_predictions.py --dataset` at a working single-sample tree containing the new `output/copynumber_background.csv` and native calls; write its results to that tree's `othertools/`. It preserves the original unmerged coordinate conventions. For paired normalization set `CNV_RAW_DIR` to `paired_samples/raw_for_normalization`, `CNV_EVAL_OUT` to a new directory and run `standardize_paired_sample_predictions.py`; the final `my_method_greedy` event rule is the default. The normalized_calls folder it generates can replace the working paired tree's normalized_calls. Do not overwrite archived inputs/results.

## Optional raw-sequencing reproduction

The following scripts preserve the tool settings and complete ART calls, but need the external references/software/resources specified at their top. Configure these paths for your machine. `ABSCNV_PROJECT` / `ABSCNV_BASE` override applicable historical project roots; shared length-A roots use `ABSCNV_SAMPLE_A_FASTQ` / `ABSCNV_SAMPLE_A_ALIGNMENT`, and length design sources use `ABSCNV_LENGTH_SOURCE`. Internal code calls resolve to this Research directory. These expensive stages were not rerun during packaging.

| Stage | Code | Required external inputs / products |
|---|---|---|
| Paired profiles / FASTQ | `build_paired_benchmark_truth.py`, `simulate_paired_benchmark_reads.py`, `generate_paired_benchmark_fastq.sh` | Archived gain/loss source + chromosome FASTAs → profiles, truth, ART FASTQ |
| Length profiles / FASTQ | `build_length_benchmark_profiles.py`, `simulate_length_benchmark_reads.py`, `generate_length_benchmark_fastq.sh`, `run_length_benchmark_experiments.sh`, `run_one_length_experiment.sh` | Archived master/profile sources + FASTAs → six independent experiments |
| Simulated alignment / RD | `align_simulated_reads_and_export_depth.sh`, `align_single_sample_benchmark.sh`, `align_length_benchmark_reads.sh`, `export_autosomal_read_depth.C` | FASTQ + reference and indexes → BAM/ROOT/1-kb RD |
| Single comparators | `run_cnvnator_single_sample.sh`, `run_freec_single_sample.sh`, `run_cony_single_sample.sh`, `build_cony_depth_windows.py`, `call_cony_single_sample.R`, `cony_compatible_runtime.R`, `prepare_rdxplorer_reference.py`, `run_rdxplorer_single_sample.sh`, `call_rdxplorer_chromosome.R` | BAM/reference + external tool installations → native predictions; CONY/RDXplorer baseline=2 |
| Paired comparators | `run_paired_sample_comparators.sh` | A/B BAM + external CNVkit/BIC-seq2/FREEC/ichorCNA resources → native calls |
| Public WGS | `download_public_wgs_reads.sh`, `convert_and_validate_public_fastq.sh`, `align_public_wgs_reads.sh`, `extract_public_wgs_read_depth.sh` | GSE118715 reads + reference → sorted BAM and non-GC RD |
| U251 WGBS | `process_wgbs_reads_and_export_depth.sh`, `export_gc_corrected_read_depth.C` | Raw WGBS + Bismark/split references → deduplicated BAM, ROOT and GC-corrected RD |

ART is **Q 2.5.8 (June 7, 2016)**; set `ABSCNV_ART` to the executable. Read length=100, fragment mean=200 and SD=10. Paired and length simulators retain their distinct seed rules. `build_paired_benchmark_truth.py` takes GAIN LOSS REF_DIR OUT_DIR; `build_length_benchmark_profiles.py` takes GAIN LOSS REF_DIR MASTER OUT_ROOT BASELINE. The required design sources are under paired_samples/source and length_series/simulation_sources. To reproduce archived RD-based analyses exactly, use the archived RD inputs; regenerating FASTQ is a separate expensive workflow.

CONY: `bash Research/run_cony_single_sample.sh BAM REF_FASTA NEW_OUTPUT`. RDXplorer: run `prepare_rdxplorer_reference.py --official-source RDXPLORER_3_2_SOURCE --reference REF_FASTA --output-dir NEW_RESOURCES`, then `bash Research/run_rdxplorer_single_sample.sh BAM REF_FASTA NEW_RESOURCES/patched NEW_OUTPUT`. Upstream RDXplorer is an external installation, not a missing internal Research module. CNVnator/FREEC shell entries retain the benchmark's project/reference/tool settings. U251 WGBS entry arguments are SAMPLE R1 R2 BISMARK_GENOME SPLIT_FASTA_DIR NEW_OUTPUT_DIR; `RUN_METHYLATION=1` additionally enables the method's methylation extraction, which is not required for CNV RD. The WGBS entry is reconstructed from the method plus recovered CNVnator commands, not claimed to be a recovered historical master script.

## Definitions, expected checks and provenance

- Single truth=372; unmerged predictions: AbsCNV 2385, CNVnator 321, FREEC 470, CONY 11954, RDXplorer 123180. AbsCNV base F1≈0.968965; many-to-many event F1≈0.822535; strict either-80% one-to-one F1≈0.169024.
- Paired truth=50; predictions: AbsCNV 51, BICseq2 59, CNVkit 34, FREEC 69, ichorCNA 52. AbsCNV base F1≈0.972737; many-to-many event F1≈0.970297. Native CN fields are prioritized for CNVkit/ichorCNA. Only the current AbsCNV branch is retained.
- Length-series base TP sums accepted one-to-one intersections; it is not the main benchmarks' global base-union metric. At 20kb, event F1≈0.857143 and base F1≈0.490566.
- Public calls: 85/102/44; final filtered directional LUMPY support: 59/85, 50/102, 27/44. At least 1 bp overlap is the support threshold, not an independent validation of exact breakpoints. Figure 3A uses the latest saved RD/GMM; no refit is performed by its exporter/plotter.
- U251 inheritance is an operational filtered-event overlap definition, not independent equality of every absolute CN state. Spatial non-overlap is not proof of a biologically de novo event; U1.2.2 same-direction overlap is 0.0574%.
- Public raw data: GSE118715; U251 accession reported in the manuscript: OMIX015037 (publication/access status must be stated at submission). Raw reads, BAMs, ROOTs and whole reference FASTAs are not duplicated here. No private reviewer URL is included.
- Inspected U251 BAM headers report Bismark 0.24.2 while the manuscript says 0.22.3; public WGS alignment records do not establish a separate deduplication step. These text/provenance discrepancies remain for author review. License decisions remain deferred as instructed.

Evaluation/normalization equivalence, alias restoration and syntax were checked after renaming. Large FASTQ/BAM and full inference jobs were not rerun for packaging. Final supplementary workbooks are author-prepared files; no nonexistent automatic workbook-builder is implied. Algorithmic helper functions retain their established function names; public script filenames describe their purpose.


## Generate workflow inputs, truth, and source designs

The workflow example uses an original algorithm snapshot, kept in `workflow_functions_snapshot.jl` and `workflow_process_snapshot.jl`. These two files are required by `generate_workflow_intermediates.jl`; run this entry in its own Julia process. They preserve the method illustration rather than replacing the current core software used for benchmarks.

```text
Rscript Research/simulate_workflow_read_depth.R "NEW_WORKFLOW" revised
julia --project=. Research/generate_workflow_intermediates.jl "NEW_WORKFLOW"
Rscript Research/plot_method_workflow.R "NEW_WORKFLOW" "NEW_RESULTS/workflow"
python Research/build_single_sample_truth.py --profile "DATASETS/single_sample/benchmark/Sample_A.absolute_CN_segments.tsv" --output "NEW_RESULTS/Sample_A.truth_events.tsv"
python Research/export_lumpy_support_summary.py --dataset "DATASETS/public_wgs" --output "NEW_RESULTS/filtered_LUMPY_summary.csv"
```

`build_length_benchmark_master.py GAIN LOSS REF_DIR NEW_MASTER_DIR` generates the original combined 120-event design, profiles and design summary. Its `Sample_B_vs_Sample_A.benchmark.tsv` supplies the MASTER argument of `build_length_benchmark_profiles.py`; the latter splits the design into six experiments. Both require the original chromosome FASTAs and gain/loss source files. The archived master table can be used directly when only rerunning the six experiments.

For the archived workflow snapshots, `plot_method_workflow.R` accepts `DATASETS/Figure1_AbsCNV workflow` directly. It also accepts the historical directories created by the workflow generator. The simulator defaults to the archived revised noise profile (seed 20260725). Figure 1 represents a controlled method illustration, separate from the genome-wide performance benchmarks.
