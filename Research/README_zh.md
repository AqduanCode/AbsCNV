# AbsCNV 论文研究代码

将本文件夹 **`Research/`单独复制到已有AbsCNV仓库根目录**，与`MyFunctions.jl`、`MyProcess.jl`和`Project.toml`并列。这里不重复放已有的软件主体、示例和测试数据。英文说明见`README.md`。

配套数据包含论文分析的输入和已保存结果。下文的 DATASETS 指准备后的 datasets 目录；使用与数据对应的代码版本，并将新结果写入独立目录。

## 环境与数据准备

使用Julia 1.11（所附Manifest基于1.11.5）、Python 3.9+，以及安装了`jsonlite`、`data.table`、`ggplot2`、`scales`、`systemfonts`、`ragg`、`svglite`的R。Python依赖按下方命令安装。图中文字使用Arial；必要时安装字体或使用支持的`--font-dir`参数。Illustrator仅用于采样示意图JSX。只有从原始测序重新运行时才需要外部CNV工具、ART、比对工具及ROOT环境。

归档中，相同RD只保存一份：单/配对A共用一个文件，六长度A共用另一个文件；两组数据并不相同，未互相合并。**推断前先运行`prepare_working_data.py`**，生成新的工作副本、解压文本并通过硬链接或复制恢复8个实验输入路径。归档原件不变；原压缩绘图快照与VCF保留。所有分析写入新输出目录。

## 按论文顺序的代码、输入与输出

表中数据路径相对于准备好的`datasets/`。`read_single_sample_predictions.py`与`paired_interval_metrics.py`是评价入口必需的辅助模块。

| 论文内容 | 按用途命名的入口 | 归档输入 | 主要输出 |
|---|---|---|---|
| Figure 1 | `plot_method_workflow.R` | `Figure1_AbsCNV workflow/`下7份已保存的示意表 | 方法示意b–f面板 |
| 单样本比较 | `call_single_sample_copy_numbers.jl`、`standardize_single_sample_predictions.py`、`evaluate_single_sample_predictions.py` | `single_sample/input`、`benchmark`、`native`、`output/copynumber_background.csv`、`othertools` | CN预测、五工具未合并结果、四口径评价及可再生成的匹配明细 |
| 配对比较 | `call_paired_sample_differences.jl`、`standardize_paired_sample_predictions.py`、`evaluate_paired_sample_predictions.py` | `paired_samples/input`、`benchmark`、`raw_for_normalization`、`normalized_calls`、`output` | 配对CNV、优先原生CN字段的标准化结果、四口径评价 |
| Figure 2 | `plot_benchmark_performance.py` | 单/配对`evaluation/*four_standard*summary.csv` | 175×113mm PDF/PNG、20行实际绘图指标 |
| 六组长度 / Figure S1 | `call_length_benchmark_differences.jl`、`evaluate_and_plot_length_performance.R` | `length_series/input`、`benchmark`、`output` | 六组最终CNV、事件/碱基指标、长度性能图 |
| 公开WGS / Figure 3 | `call_public_wgs_differences.jl`、`export_public_copy_number_peaks.jl`、`plot_public_copy_number_histograms.R`、`plot_public_pairwise_copy_number_profiles.R` | `public_wgs/input`、`output`、`histogram_current/cn_peaks.csv` | 三组两两结果、最新模型峰直方图、两两CNV轨迹 |
| LUMPY / Table S5 | `evaluate_lumpy_support.py` | `public_wgs/LUMPY_SV`、`output/CNVresult_filter.csv` | 最终方向一致的支持数量与比例 |
| U251 / Figure 4 | `call_clonal_lineage_differences.jl`、`plot_lineage_chromosome_events.py`、`evaluate_lineage_event_inheritance.py`、`summarize_lineage_specific_differences.py` | `u251/input`十份CSV及`output`中的CNV/CN背景 | 染色体矩阵与位置放大图、继承/同向/反向比例、位置特异长度 |
| Figure S3 | `export_lineage_copy_number_peaks.jl`、`plot_lineage_copy_number_histograms.R` | U251模型、归一化RD、`peak_parameters` | 十样本直方图PDF |
| Figure S2 / S4 | `draw_clonal_sampling_scheme.jsx`、`assemble_supplementary_lineage_panels.py` | 固定采样设计、生成的U12/U13 PNG | 可编辑谱系示意、补充图Word拼版 |

## 从归档结果重算评价与绘图

将以下`SOURCE_DATASETS`、`DATASETS`、`RESULTS`替换为实际路径。全部命令从AbsCNV仓库根目录执行；使用输出重定向前先创建RESULTS根目录。图2默认读归档的全精度评价表；如果要绘制本次重算的指标，使用`--single-metrics RESULTS/single/single_sample_four_standard_summary.csv --paired-metrics RESULTS/paired/four_standard_f1_summary.csv`。

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

十样本直方图：设置`CNV_PLOT_PROJECT_DIR=DATASETS/u251`及`CNV_PLOT_OUTPUT_DIR=RESULTS/u251_hist`，再执行`Rscript Research/plot_lineage_copy_number_histograms.R`。峰参数CSV须先生成在该输出目录。Bash用`export NAME=value`设置环境变量，PowerShell用`$env:NAME="value"`。在Illustrator运行`draw_clonal_sampling_scheme.jsx`后选择保存目录。补充图拼版脚本生成DOCX；归档PDF来自Word排版导出，不是该Python脚本自动生成。

## 从RD重新运行AbsCNV

只调整了入口路径及代码文件名，推断参数和计算调用不变。各实验使用新输出目录。最后一条命令依次换为50kb、100kb、200kb、500kb、1000kb；长度实验输出根目录下会生成`output/LABEL/`。

```text
julia --project=. Research/call_single_sample_copy_numbers.jl "DATASETS/single_sample/input" "RESULTS/single_inference"
julia --project=. Research/call_paired_sample_differences.jl "DATASETS/paired_samples/input" "RESULTS/paired_inference"
julia --project=. Research/call_public_wgs_differences.jl "DATASETS/public_wgs/input" "RESULTS/public_inference"
julia --project=. Research/call_clonal_lineage_differences.jl "DATASETS/u251/input" "RESULTS/u251_inference"
julia --project=. Research/call_length_benchmark_differences.jl 20kb "DATASETS/length_series" "RESULTS/length_inference"
```

新单样本推断完成后，将新的`output/copynumber_background.csv`和原生工具结果组织到工作single_sample目录，再运行`standardize_single_sample_predictions.py --dataset WORK_SINGLE_SAMPLE --output-dir WORK_SINGLE_SAMPLE/othertools`；坐标转换和未合并事件口径与原分析一致。配对标准化时，设置`CNV_RAW_DIR`为`paired_samples/raw_for_normalization`，`CNV_EVAL_OUT`为新输出目录，再运行`standardize_paired_sample_predictions.py`。默认事件口径是最终采用的`my_method_greedy`，生成的normalized_calls可用于新的工作配对数据目录。不要覆盖归档原件。

## 可选：从原始测序重新处理

以下代码保留工具设置和完整ART调用，需要配置代码顶部列出的外部软件、参考和资源路径。部分入口可用`ABSCNV_PROJECT`或`ABSCNV_BASE`覆盖项目根目录；共享长度A的来源由`ABSCNV_SAMPLE_A_FASTQ`、`ABSCNV_SAMPLE_A_ALIGNMENT`控制，长度设计源由`ABSCNV_LENGTH_SOURCE`指定。内部脚本调用均指向本Research目录。整理时没有重跑这些昂贵流程。

| 阶段 | 代码 | 外部输入与产物 |
|---|---|---|
| 配对模拟 | `build_paired_benchmark_truth.py`、`simulate_paired_benchmark_reads.py`、`generate_paired_benchmark_fastq.sh` | gain/loss源、染色体FASTA→profile、真值、ART FASTQ |
| 六长度模拟 | `build_length_benchmark_profiles.py`、`simulate_length_benchmark_reads.py`、`generate_length_benchmark_fastq.sh`、`run_length_benchmark_experiments.sh`、`run_one_length_experiment.sh` | 主真值和profile源、FASTA→六组独立实验 |
| 模拟比对与RD | `align_simulated_reads_and_export_depth.sh`、`align_single_sample_benchmark.sh`、`align_length_benchmark_reads.sh`、`export_autosomal_read_depth.C` | FASTQ、参考与索引→BAM、ROOT、1-kb RD |
| 单样本比较工具 | `run_cnvnator_single_sample.sh`、`run_freec_single_sample.sh`、`run_cony_single_sample.sh`、`build_cony_depth_windows.py`、`call_cony_single_sample.R`、`cony_compatible_runtime.R`、`prepare_rdxplorer_reference.py`、`run_rdxplorer_single_sample.sh`、`call_rdxplorer_chromosome.R` | BAM、参考、外部工具→原生预测；CONY/RDXplorer基线为2 |
| 配对比较工具 | `run_paired_sample_comparators.sh` | A/B BAM及四工具资源→原生预测 |
| 公开WGS | `download_public_wgs_reads.sh`、`convert_and_validate_public_fastq.sh`、`align_public_wgs_reads.sh`、`extract_public_wgs_read_depth.sh` | GSE118715读段和参考→排序BAM、非GC的RD |
| U251 WGBS | `process_wgbs_reads_and_export_depth.sh`、`export_gc_corrected_read_depth.C` | 原始WGBS及Bismark/拆分参考→去重BAM、ROOT、GC校正RD |

ART版本为**Q 2.5.8（2016-06-07）**；通过`ABSCNV_ART`指定程序。读长100、插入片段均值200/标准差10，配对版和长度版各自的随机种子规则均保留。`build_paired_benchmark_truth.py`参数为GAIN LOSS REF_DIR OUT_DIR；`build_length_benchmark_profiles.py`参数为GAIN LOSS REF_DIR MASTER OUT_ROOT BASELINE。设计材料见paired_samples/source及length_series/simulation_sources。要精确复算已归档的RD分析，请使用已归档RD；重新生成FASTQ属于另一条高成本流程。

CONY入口：`bash Research/run_cony_single_sample.sh BAM REF_FASTA NEW_OUTPUT`。RDXplorer先执行`prepare_rdxplorer_reference.py --official-source RDXPLORER_3_2_SOURCE --reference REF_FASTA --output-dir NEW_RESOURCES`，再执行`bash Research/run_rdxplorer_single_sample.sh BAM REF_FASTA NEW_RESOURCES/patched NEW_OUTPUT`；RDXplorer原始工具需要外部安装，不是漏收的Research内部代码。CNVnator/FREEC的shell入口保留基准所用项目、参考和工具设置。U251 WGBS入口参数为SAMPLE R1 R2 BISMARK_GENOME SPLIT_FASTA_DIR NEW_OUTPUT_DIR；`RUN_METHYLATION=1`才另行执行甲基化提取，该步骤不参与CNV RD计算。这个WGBS入口按方法描述和找回的CNVnator流程重新整理，不冒充已找回完整历史master脚本。

## 口径、预期结果与来源

- 单样本真值372个；AbsCNV、CNVnator、FREEC、CONY、RDXplorer分别为2385、321、470、11954、123180个未合并事件。AbsCNV碱基F1约0.968965，多对多事件F1约0.822535，either-80%严格一对一F1约0.169024。
- 配对真值50个；AbsCNV、BICseq2、CNVkit、FREEC、ichorCNA分别为51、59、34、69、52个事件。AbsCNV碱基F1约0.972737，多对多事件F1约0.970297。CNVkit/ichorCNA优先使用原生CN字段；仅保留当前AbsCNV分支。
- 长度实验碱基TP只累加已接受的一对一匹配交集，与主比较的全局碱基并集口径不同。20kb事件F1约0.857143，碱基F1约0.490566。
- 公开数据CNV为85/102/44个；最终LUMPY方向一致支持为59/85、50/102、27/44。支持阈值为至少重叠1bp，不是精确断点独立验证。Figure3A使用最新保存的RD/GMM，参数导出和绘图不重新拟合。
- U251继承按筛选后事件的区间重叠作操作性划分，不是对每个碱基绝对CN相等的独立检验；位置不重叠也不证明生物学新发。U1.2.2同向变化比例为0.0574%。
- 公开原始数据为GSE118715；稿件报告的U251登录号为OMIX015037，投稿时需说明当时的公开/访问状态。这里不重复归档巨大原始读段、BAM、ROOT或参考FASTA；不包含审稿人私有链接。
- 已查U251 BAM头部为Bismark0.24.2，而论文为0.22.3；公开WGS的现有比对记录未证明单独去重步骤。相关正文/来源差异仍由作者核对。许可证事项按要求暂缓。

重命名后检查了评价/标准化等价性、输入别名还原和语法。没有为整理而重跑大型FASTQ/BAM或完整推断任务。Zenodo中的`analysis_metadata.json`保存版本和来源事实，校验清单固定数据与源码快照。附表是作者已整理的最终成品，没有虚构不存在的自动附表生成程序。辅助函数保留既有函数名，公开的脚本文件名按代码用途命名。


## 流程图输入、单样本真值与长度主设计

- `simulate_workflow_read_depth.R NEW_WORKFLOW revised`：按原始随机种子和 revised 噪声参数生成三样本20 Mb RD与真值。
- `generate_workflow_intermediates.jl NEW_WORKFLOW`：调用同目录的 `workflow_functions_snapshot.jl` 和 `workflow_process_snapshot.jl`，生成 Figure 1 的七份绘图输入。两个函数快照保留原方法示例的算法版本；请单独启动 Julia 进程运行，与当前基准分析所用核心代码区分。
- `build_single_sample_truth.py --profile PROFILE --output NEW_TRUTH`：使用原始真值读取及方向合并规则，从绝对CN profile恢复单样本事件真值。
- `build_length_benchmark_master.py GAIN LOSS REF_DIR NEW_MASTER_DIR`：生成120事件的原始合并设计，输出主真值、A/B profile及设计摘要；再由 `build_length_benchmark_profiles.py` 拆分六组。
- `export_lumpy_support_summary.py --dataset PUBLIC_WGS_DATASET --output NEW_CSV`：调用现有LUMPY评价程序，按原公式导出正文支持率汇总CSV。

R模拟脚本需要 data.table 和 jsonlite。运行命令在仓库根目录前缀加 `Rscript Research/`、`julia --project=. Research/` 或 `python Research/`。`plot_method_workflow.R` 已兼容归档的 Figure1_AbsCNV workflow 目录和历史带面板字母的目录名，无需手动复制目录。
