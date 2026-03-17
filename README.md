# RPOE_MRI

**R**esearch **P**rogram **O**f **E**scellence (**RPOE**) — 7 Tesla MRI Processing and Analysis

This repository contains the full pipeline for processing, analyzing, and visualizing 7T MRI data from the RPOE study, including structural (T1w/T2w) and functional (resting-state and task-based) MRI modalities. Analyses span anatomical volumetrics, functional signal characterization, task-fMRI design and activation modeling, functional connectivity, network correspondence, and mitochondrial mapping — all linked to behavioral and cognitive phenotypes (language, IQ).

---

## Repository Structure

```
RPOE_MRI/
├── src/
│   ├── processing-pipelines/       # Shell scripts for raw data processing (DICOM → analysis-ready)
│   └── *.R / *.py / *.m            # Analysis and visualization scripts
└── figs/                           # Output figures
```

---

## Processing Pipelines

Located in [`src/processing-pipelines/`](src/processing-pipelines/), these Bash scripts handle all upstream data processing steps and are designed to run on an HPC cluster (SLURM-based job submission via `99_pipelines-submission.sh`).

| Script | Description |
|--------|-------------|
| `01_01_DICOM-download.sh` | Download raw DICOM files from the scanner/archive |
| `01_02_DICOM-conversion.sh` | Convert DICOMs to NIfTI format |
| `02_anat-T1_V2.sh` | T1-weighted anatomical preprocessing (bias field correction, skull stripping, segmentation) |
| `03_anat-T2_V2.sh` | T2-weighted anatomical preprocessing |
| `04_anat-MNI-registration.sh` | Nonlinear registration of anatomical images to MNI standard space |
| `05_anat-T2-T1-ratio.sh` | Compute voxelwise T2/T1 ratio maps (myelin proxy) |
| `06_anat-3Dprint.sh` | Generate 3D-printable surface meshes from anatomical data |
| `07_func-preprocessing.sh` | Functional MRI preprocessing (slice timing, motion correction, smoothing, normalization) |
| `08_01_func-psvc-AFNI-subject-level-analysis.sh` | Subject-level task-fMRI GLM analysis using AFNI (PSVC task) |
| `08_02_func-psvc-AFNI-group-level-analysis.sh` | Group-level task-fMRI analysis using AFNI |
| `08_03_func-psvc-AFNI-WAT-resid-to-lang-analysis.sh` | Residualized WAT-task fMRI signal mapped to language phenotypes |
| `08_04_func-rest-AFNI-to-lang-iq-PCs.sh` | Resting-state fMRI AFNI analysis regressed on language/IQ principal components |
| `99_pipelines-submission.sh` | HPC job submission wrapper for all pipeline scripts |

---

## Analysis Scripts

Located in [`src/`](src/), these R, Python, and MATLAB scripts handle downstream statistical analysis and visualization.

### Participant Metadata & QC

| Script | Description |
|--------|-------------|
| `00_participants-MRI-metadata-and-cleaning.R` | Load and clean participant metadata; apply MRI QC filters and exclusion criteria |

### Task-fMRI (PSVC) Behavioral Log Processing

| Script | Description |
|--------|-------------|
| `01_processing-PSVC-tfMRI-log-files.R` | Parse and clean PSVC task behavioral log files from scanner runs |
| `02_PSVC-tfMRI-log-corr-w-other-tests.R` | Correlate PSVC task performance with external cognitive/behavioral assessments |

### Structural MRI

| Script | Description |
|--------|-------------|
| `03_anat-volumetrics-cleaning.R` | Clean and QC anatomical volumetric outputs (cortical thickness, subcortical volumes) |

### Functional MRI Analyses

| Script | Description |
|--------|-------------|
| `04_00_fMRI-signal.R` | Compute and summarize raw fMRI signal metrics across ROIs |
| `04_01_PSVC-fMRI-task-design.R` | Build PSVC task design matrices for GLM modeling |
| `04_02_PSVC-AFNI-analysis.R` | Post-process AFNI GLM outputs; extract beta maps and contrasts |
| `04_04_NCT-run.R` | Run network correspondence toolbox (NCT) analysis on functional connectivity matrices |
| `04_04_NCT.py` | Python helper for NCT computations |
| `04_05_mito-maps.R` | Generate and analyze mitochondrial function proxy maps |
| `04_06_task-falff-euc-or-dotprod-to-others-and-themselves.R` | Compute Euclidean distance and dot product similarity between fALFF maps |
| `04_08_functional-connectivity.R` | Whole-brain and ROI-based functional connectivity analysis |
| `04_08_lang-iq-PCs-to-MRI-AFNI-prep.R` | Prepare language and IQ PCs as AFNI regressors for fMRI analyses |
| `04_08_lang-iq-PCs-to-MRI-SPM-prep.R` | Prepare language and IQ PCs for SPM second-level analyses |
| `04_09_spm_falff_to_lang_iq_pcs.m` | SPM MATLAB script: regress fALFF maps on language/IQ PCs |
| `04_09_spm_reho_to_lang_iq_pcs.m` | SPM MATLAB script: regress ReHo maps on language/IQ PCs |
| `04_09_spm_task_resids_to_lang_iq_pcs.m` | SPM MATLAB script: regress task-fMRI residuals on language/IQ PCs |

### Visualization

| Script | Description |
|--------|-------------|
| `05_viz-participant-proportion-difference-than-bg.py` | Visualize participant-level fMRI signal proportion difference relative to background |
| `99_aes-viz.R` | Shared aesthetic/theme functions for R-based figures (ggplot2) |
| `99_viz.py` | Shared visualization utilities for Python-based figures (matplotlib/nilearn) |

---

## Data

This repository contains **code only**. Raw MRI data, processed NIfTI files, and participant-level outputs are not included due to data use agreements and IRB restrictions. Access to the RPOE dataset may be requested through the appropriate data governance process.

---

## Citation

If you use this code, please cite the associated publication (forthcoming) or contact the author.

---

## Author

**Muhammad Elsadany**  
PhD Candidate, Interdisciplinary Graduate Program in Genetics (Computational Genetics Track)  
Department of Psychiatry, University of Iowa  
[GitHub](https://github.com/melsadany)

---

## License

This repository is for academic and research use. Please contact the author before reuse or redistribution.
