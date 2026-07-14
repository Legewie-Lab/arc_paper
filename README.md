# RON kinetic-model fitting and knockdown analysis code

This repository contains MATLAB code and input data used for RON kinetic-model fitting, empirical error-model estimation, downstream RON model analysis, and RBP knockdown prediction analysis.

The code was organized into separate analysis folders:

- `RON/`: RON one-step, two-step, and control model fitting scripts for HPC, plus downstream fitting analysis.
- `KD/`: knockdown prediction analysis scripts.
- `error model/`: empirical error-model estimation scripts.
- `Fitting.m`: a local example of the older RON fitting workflow. The main RON fitting reported here was performed on HPC with the scripts in `RON/`.

Each analysis folder contains its own `Data/` directory with the input files required by that part of the workflow.

## Software

The code was developed and run with:

- MATLAB R2023a.
- Data2Dynamics (D2D) / `arFramework3`.
- Java enabled in MATLAB, required by D2D.
- MATLAB Optimization Toolbox, used for nonlinear least-squares fitting such as `lsqnonlin`.
- MATLAB Parallel Computing Toolbox, used/recommended for `parfor` sections and D2D parallel compilation on HPC.

The scripts also use standard MATLAB table, plotting, and file I/O functions such as `readtable`, `writetable`, `savefig`, and `print`.

## Repository Layout

```text
.
├── Fitting.m
├── README.md
├── RON/
│   ├── Data/
│   ├── Setup_RON_onestep_compile.m
│   ├── Setup_RON_twostep_compile.m
│   ├── Setup_RON_control_compile.m
│   ├── Run_RON_onestep_LHS.m
│   ├── Run_RON_twostep_LHS.m
│   ├── Run_RON_control_LHS.m
│   └── analyse_fitting.m
├── KD/
│   ├── Data/
│   ├── kd_fitting_prediction_chi.m
│   └── KD_fitting_prediction_RMSE.m
└── error model/
    ├── Data/
    └── error_model.m
```

## Folder Overview

`RON/` contains the D2D compile and LHS fitting scripts for the RON models:

- `Setup_RON_onestep_compile.m`: compiles the one-step RON model on the HPC/Linux system.
- `Setup_RON_twostep_compile.m`: compiles the two-step RON model on the HPC/Linux system.
- `Setup_RON_control_compile.m`: compiles the two-step control model on the HPC/Linux system.
- `Run_RON_onestep_LHS.m`: runs one-step model fitting as a SLURM array job.
- `Run_RON_twostep_LHS.m`: runs two-step model fitting as a SLURM array job.
- `Run_RON_control_LHS.m`: runs control model fitting as a SLURM array job.
- `analyse_fitting.m`: post-processes fitted one-step and two-step model objects and generates comparison plots and fitted perturbation factors.
- `Data/`: RON model definition/data files and related inputs.

`KD/` contains the RBP knockdown prediction analyses:

- `kd_fitting_prediction_chi.m`: knockdown prediction analysis using the chi-square style objective.
- `KD_fitting_prediction_RMSE.m`: knockdown prediction analysis using the RMSE objective.
- `Data/`: input data and intermediate results used by the KD analysis.

`error model/` contains the empirical error-model workflow:

- `error_model.m`: estimates empirical isoform standard-deviation/error models from replicate data, generates per-isoform and global `sigma(mu)` functions, and writes D2D input/error tables.
- `Data/`: replicate input data and generated/required error-model inputs.

`Fitting.m` is retained as a local example script for RON fitting. It is not the main HPC workflow, but it documents the older local D2D setup and can be adapted for local tests.

## Path Configuration

The HPC scripts in `RON/` assume the following project layout on the cluster:

```text
$HOME/HPC/d2d-master/arFramework3
$HOME/HPC/D2D_1
```

The run scripts also require `RESULT_DIR` to be set by the SLURM job or shell environment. `SLURM_ARRAY_TASK_ID` is used when available to seed and label array jobs.

`Fitting.m` avoids hard-coded personal paths and can be configured with environment variables:

- `D2D_FRAMEWORK_DIR`: path to D2D `arFramework3`.
- `RON_D2D1_DIR`: path to the main D2D project folder.
- `RON_D2D2_DIR`: path to the second D2D project folder used in refinement sections.
- `RON_REFERENCE_DIR`: path to reference fitting results used for bounds or starting values.

`error model/error_model.m` writes the main error-model output to the current directory by default. Set `ERROR_MODEL_OUTPUT_DIR` to write it elsewhere.

## Typical Workflow

1. Estimate or update the empirical error model:
   ```matlab
   cd('error model')
   addpath('Data')
   error_model
   ```

2. Compile RON models on the HPC system:
   ```matlab
   cd('RON')
   Setup_RON_onestep_compile
   Setup_RON_twostep_compile
   Setup_RON_control_compile
   ```

3. Submit SLURM array jobs for LHS fitting. Each job should set `RESULT_DIR`, then run one of:
   ```matlab
   Run_RON_onestep_LHS
   Run_RON_twostep_LHS
   Run_RON_control_LHS
   ```

4. Analyze fitted RON models:
   ```matlab
   cd('RON')
   addpath('Data')
   analyse_fitting
   ```

5. Run knockdown prediction analyses:
   ```matlab
   cd('KD')
   addpath('Data')
   kd_fitting_prediction_chi
   KD_fitting_prediction_RMSE
   ```

## Figures Generated by This Code

The following table maps figures in the published article [*Coordinated alternative splicing decisions via stepwise exon definition*](https://doi.org/10.1093/nar/gkag464) (Nucleic Acids Research, 2026) to the corresponding MATLAB scripts and plotting sections in this repository.

| Published figure | MATLAB source | Relevant section or generated output |
|---|---|---|
| Fig. 2B and Fig. 2D | `RON/analyse_fitting.m` | Section `Plot 1.1: Model vs Data (exon 2 only & exon 1 3)`. The one-step and two-step exon-2 fits are the two panels generated as `Plot1.1_Model_vs_Data_Exon2.svg` and `.fig`. |
| Fig. 2C | `RON/analyse_fitting.m` | Section `Plot 2.1: PSI vs Efficiency comparison (Data + One-step + Two-step)`, generating `Plot2_PSI_vs_Efficiency_exon2.svg` and `.fig`. |
| Fig. 5A | `KD/kd_fitting_prediction_chi.m` | In the `data & Plot` section, the control-versus-KD isoform scatter plots are generated by the loop using `data_co`, `contro_matx`, `ihek`, and `ikd`. The PUF60 and SMU1 panels were selected for the main figure. |
| Fig. 5B | `KD/kd_fitting_prediction_chi.m` | Section explicitly labelled `Figure 5B`, which interleaves fitted and predicted normalized chi-square values and plots the model-selection heatmap. The underlying values are saved as `Figure5B_data.mat`. |
| Fig. 5C | `KD/kd_fitting_prediction_chi.m` | Section labelled `Figure 5C`, using the selected PUF60 and SMU1 models (`j=[2,4]`) to plot predicted versus experimental isoform frequencies. |
| Fig. 5D | `KD/kd_fitting_prediction_chi.m` | Section labelled `Figure 5D: PSI vs. Efficiency`, again using PUF60 and SMU1 (`j=[2,4]`). |
| Supplementary Fig. S6A | `error model/error_model.m` | Section `error model -- alphas` plots the effect of different alpha values on the error-model shape. The empirical symmetric-error visualization is also supported by the deviation-distribution sections later in the same script. |
| Supplementary Fig. S6B | `error model/error_model.m` | Sections beginning with `Load table` generate replicate-minus-mean deviations for individual isoforms and the combined `All Isoforms` plot. |
| Supplementary Fig. S6C-G | `error model/error_model.m` | The error-model calibration and plotting block with `Figure_ErrorModel_sigma_vs_mu` generates the five isoform-specific mean-versus-standard-deviation fits. Outputs include `Figure_ErrorModel_sigma_vs_mu.fig`, `Figure_ErrorModel_sigma_vs_mu.svg`, and `ErrorModel_sigma_vs_mu.svg`. |
| Supplementary Fig. S8A and S8C | `RON/analyse_fitting.m` | The `Model_vs_Data_Exon1_3` block plots one-step and two-step model fits for mutations affecting the outer exons, generating `Plot1.2_Model_vs_Data_Exon1_3.svg` and `.fig`. |
| Supplementary Fig. S9, top row | `KD/kd_fitting_prediction_chi.m` | In the `data & Plot` section, the loop plotting `contro_matx` against each RBP-KD dataset generates the measured control-versus-KD isoform shifts. |
| Supplementary Fig. S9, middle row | `KD/kd_fitting_prediction_chi.m` | Section labelled `Figure S5 (middel): scatter plot of Model predictions vs. experimental values`. Despite the older section label, this code corresponds to the final Supplementary Fig. S9 middle row. |
| Supplementary Fig. S9, bottom row | `KD/kd_fitting_prediction_chi.m` | Section labelled `Figure S5: PSI vs Efficiency (Top + btm)`, particularly the loop using `best_sort_pi`, generates the measured and predicted PSI-efficiency plots. Despite the older section label, these panels correspond to the final Supplementary Fig. S9 bottom row. |

The manuscript figures were assembled from these MATLAB outputs together with panel labels, legends, model schematics, and other layout elements added during final figure preparation. Therefore, some final multi-panel figures do not correspond to a single exported MATLAB file.

## Notes

- Some scripts expect intermediate `.mat` files produced by earlier workflow steps or by HPC fitting runs.
- Error parameters are not fitted in the D2D LHS scripts (`ar.config.fiterrors = 0`); measurement uncertainty is supplied by the precomputed empirical error model.
- Random seeds in the SLURM fitting scripts are offset by `SLURM_ARRAY_TASK_ID`, so each array task starts from a different seed.
- The folder name `error model/` contains a space. In MATLAB commands, wrap it in quotes as shown above.
