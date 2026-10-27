# Synthetic example data

`Sample_A.csv` and `Sample_B.csv` are deterministic demonstration inputs, not biological observations.

Each file contains 2,400 consecutive 1 kb bins on `chr1`, using zero-based starts and exclusive upper edges. There are eight 300-bin depth blocks:

| Sample | Block means |
| --- | --- |
| Sample_A | 50, 100, 50, 100, 150, 100, 150, 100 |
| Sample_B | 100, 50, 50, 100, 100, 150, 150, 100 |

For zero-based position `i` within each block, RD is `block_mean + 4 * sin(2*pi*i/37)`, rounded to six decimal places. The main depth peak is 100, assigned CN2 by the demonstration's assumption; 50 and 150 are CN1-like and CN3-like depths.

The two profiles have equal numbers of each depth block and different block ordering. They exercise aligned data loading, peak fitting, segmentation, pairwise calling, and export. They do not establish sensitivity, specificity, statistical calibration, background-detection accuracy, or behavior at chromosome boundaries. Background fitting is enabled but this dataset has no dedicated broad-background interval.

Run from the repository root:

```bash
julia --project=. examples/run_background_pipeline.jl
```
