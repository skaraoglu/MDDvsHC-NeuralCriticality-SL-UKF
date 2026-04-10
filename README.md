<div align="center">

# Cross-Condition Stuart-Landau Estimation in Major Depressive Disorder

### Multi-Method Bifurcation Parameter Comparison Across Healthy Controls & MDD

[![R](https://img.shields.io/badge/R-≥4.2-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-≥3.9-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-Academic_Use-lightgrey)]()
[![Status](https://img.shields.io/badge/Status-Analysis_Complete-brightgreen)]()
[![Atlas](https://img.shields.io/badge/Parcellation-216_ROI-blue)]()
[![Subjects](https://img.shields.io/badge/N-49_(30_HC_+_19_MDD)-orange)]()

---

*If the bifurcation parameter depends on how you estimate it,  
what exactly is each method measuring — and does depression change it?*

</div>

---

## Overview

This repository contains the analysis pipeline for a **cross-condition comparison** of whole-brain dynamics between healthy controls (HC) and unmedicated Major Depressive Disorder (MDD), using three independent estimation approaches applied to the **Stuart-Landau oscillator** — the normal form of a supercritical Hopf bifurcation.

Chapter 4 established that MDD resting-state dynamics are deeply subcritical ($a \approx -0.29$) and that neurofeedback can perturb the bifurcation parameter within-subject. This chapter asks two follow-up questions: (1) do HC dynamics occupy a different regime, and (2) do different estimation methods agree on where that regime lies?

The central finding is a **clean dissociation**: all three methods converge on a null HC–MDD difference at the whole-brain level, but they place the dynamical operating point at fundamentally different locations ($a \approx -0.26$ to $+0.01$). Time-domain methods (UKF, spectral) detect neurofeedback treatment effects; the spatial correlation method (SC-FC matching) does not.

---

## Study Design

<table>
<tr>
<td width="50%">

**Healthy Controls — HNU1 Dataset**
- 30 subjects, 10 sessions each
- Consortium for Reliability and Reproducibility (CoRR)
- Hangzhou Normal University test-retest
- Multi-session averaging for trait-level estimates

**MDD Cohort — Neurofeedback Dataset**
- 19 paired subjects (from 23 enrolled)
- Unmedicated MDD (DSM-IV-TR)
- Double-blind, sham-controlled rtfMRI-NF
- Active: left amygdala upregulation
- Sham: left intraparietal sulcus

</td>
<td width="50%">

**Acquisition**
- MDD: Siemens 3T, TR = 2.0 s, 260 volumes
- HC: 3T, TR = 2.0 s, 300 volumes
- Resting-state fMRI pre- and post-NF (MDD)
- 10 resting-state sessions (HC)

**Parcellation**
- Primary: Schaefer-200 + Melbourne-16 subcortical (216 ROIs)
- Validation: Harvard-Oxford 110-ROI (adaptive sphere radius)

**Preprocessing**
- MDD: AFNI pipeline (motion, nuisance, bandpass 0.01–0.10 Hz)
- HC: LIBR-matched confound regression with per-subject calibration

</td>
</tr>
</table>

---

## Three Estimation Approaches

The core methodological contribution is a systematic comparison of three operationalizations of the Hopf control parameter on identical data:

| Approach | Domain | What It Estimates | Typical $a$ |
|----------|--------|-------------------|-------------|
| **UKF State-Space** | Time | Per-TR decay rate via Unscented Kalman Filter | $\approx -0.26$ |
| **Spectral Lorentzian** | Frequency | BOLD power spectrum peak width via NLS fitting | $\approx -0.16$ |
| **SC-FC Matching** | Spatial correlation | Global coupling optimized to match empirical FC | $\approx +0.01$ |

The Stuart-Landau equation in complex form:

$$\dot{z} = (a + i\omega)\,z \ - \ |z|^2\,z \ + \ \eta(t)$$

| Parameter | Meaning | Estimated by |
|-----------|---------|-------------|
| $a$ | Distance from critical point | UKF / Spectral / SC-FC |
| $\omega$ | Natural oscillation frequency | Hilbert instantaneous phase / Spectral peak |
| $K$ | Inter-regional coupling | Tested → non-identifiable at TR = 2 s |

---

## Pipeline Architecture

```
┌──────────────────────────────────────────────────────────────────────────┐
│                     parcellate_hc_hnu1_v3.ipynb                          │
│  HC NIfTI  ──▸  Confound regression  ──▸  Atlas parcellation  ──▸  CSVs  │
└──────────────────────────────────┬───────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                    MDD-HC_analysis_v4.ipynb  (R kernel)                  │
│                                                                          │
│  ┌──────────────────┐  ┌──────────────────┐  ┌───────────────────────┐   │
│  │  Path 1: UKF     │  │  Path 2: Spectral│  │  Path 3: SC-FC        │   │
│  │  Multi-session   │  │  Lorentzian NLS  │  │  FC-matching global   │   │
│  │  averaging       │  │  peak-width fit  │  │  coupling optim.      │   │
│  │  (71,928 fits)   │  │                  │  │  (216 × 216 FC)       │   │
│  └────────┬─────────┘  └────────┬─────────┘  └──────────┬────────────┘   │
│           │                     │                        │               │
│           ▼                     ▼                        ▼               │
│  ┌──────────────────────────────────────────────────────────────────┐    │
│  │  Cross-condition: HC vs MDD  ·  Cross-method: UKF vs Spec vs FC  │    │
│  │  Circuit-specific: DMN/limbic ROI analysis                       │    │
│  │  Treatment replication: NF effects with corrected group labels   │    │
│  └──────────────────────────────────────────────────────────────────┘    │
└──────────────────────────────────────────────────────────────────────────┘
                                   │
                                   ▼
┌──────────────────────────────────────────────────────────────────────────┐
│                    ch5_supplement_corrected.ipynb                        │
│  NF treatment effects  ──▸  Corrected group labels (participants.tsv)    │
│  SC-FC circuit analysis ──▸  Permutation cluster test                    │
└──────────────────────────────────────────────────────────────────────────┘
```

---

## Key Results

<table>
<tr>
<td width="50%">

**Cross-Condition (HC vs MDD)**
- Whole-brain null: UKF $d = 0.13$, $p = 0.68$
- Spectral null: $d = -0.14$, $p > 0.6$
- All network-level tests: Bonferroni $p > 1.0$
- TOST equivalence confirmed

**Method Dissociation**
- UKF: $a \approx -0.26$ (deeply subcritical)
- Spectral: $a \approx -0.16$ (moderately subcritical)
- SC-FC: $a \approx +0.01$ (near-critical)
- Cross-method $r$ = 0.37 (HC), 0.53 (MDD)

</td>
<td width="50%">

**Treatment Effects (Corrected)**
- UKF depression-circuit: $p \approx 0.027$ (large effect)
- SC-FC treatment effect: collapsed to $\approx 0$ post-correction
- Time-domain detects NF effects; SC-FC does not

**Circuit-Level Analysis**
- Permutation cluster test: HC–MDD $p \approx 0.012$
- DMN/limbic regions drive the difference
- 5 suggestive ROIs in depression circuitry

**Reliability (HC, 10 sessions)**
- ICC improves from 0.06 (single) to 0.25 (averaged)
- $2.8\times$ variance reduction from multi-session averaging
- Recommendation: $\geq 5$ sessions for stable estimates

</td>
</tr>
</table>

---

## Repository Structure

```
├── R/
│   ├── sl_models.R                       # Stuart-Landau ODE definitions
│   ├── ukf_engine.R                      # Unscented Kalman Filter core
│   ├── optim.R                           # Iterative & L-BFGS-B optimization
│   ├── preprocessing.R                   # Smoothing & signal conditioning
│   └── constants.R                       # UKF tuning constants
├── data/
│   ├── source/
│   │   ├── processed rest scans/         # MDD AFNI BRIK/HEAD (rest1)
│   │   ├── processed rest2 scans/        # MDD AFNI BRIK/HEAD (rest2)
│   │   └── participants.tsv              # Authoritative group assignments
│   ├── HNU1/                             # HC raw NIfTI (30 subjects × 10 sessions)
│   │   ├── 0025427/session_1..10/rest_1/
│   │   └── ...0025456/
│   └── parcellated/
│       ├── 219roi/rest1/                 # MDD parcellated time series
│       ├── 219roi/rest2/
│       └── hc_hnu1/                      # HC parcellated time series
├── atlases/
│   ├── Schaefer2018_200Parcels_*.nii.gz
│   └── Tian_Subcortex_S1_3T_2009cAsym.nii.gz
```

---

## Critical Methodological Notes

**Group Label Integrity** — A corrected notebook (`ch5_supplement_corrected.ipynb`) was built after discovering that supplement analyses were assigning active/sham neurofeedback labels via a hardcoded fallback rather than the authoritative `participants.tsv` file. The corrected pipeline matches subjects by exam ID (E####) using the exact procedure from the Chapter 4 pipeline. This correction collapsed the SC-FC treatment effect to near zero while preserving the UKF circuit-level effect.

**Per-Subject R_SCALE Calibration** — The UKF observation noise scale is calibrated per subject to match each individual's BOLD signal variance, eliminating a confound that could produce spurious group differences from acquisition or preprocessing differences.

**Atlas Affine Robustness** — NIfTI affine extraction uses a three-level fallback chain (sform → qform → manual pixdim/qoffset construction with MNI bounding-box validation) to handle inconsistent atlas headers across datasets.

**Multi-Session Averaging** — HC estimates average across 10 sessions ($\sqrt{10} \approx 3.2\times$ noise reduction); MDD estimates average rest1 + rest2 ($\sqrt{2}$ improvement). This produces trait-level estimates suitable for between-group comparison.

---

## Requirements

<table>
<tr>
<td>

**R packages**
```
pracma, MASS, Matrix, dplyr,
tidyr, ggplot2, scales, glmnet,
igraph, parallel, zoo, psych,
equivalence
```

</td>
<td>

**Python packages**
```
nibabel, nilearn, numpy,
pandas, scipy, tqdm,
matplotlib
```

</td>
</tr>
</table>

**System:** R ≥ 4.2 · Python ≥ 3.9 · AFNI (MDD preprocessing only)

---

## Quick Start

```bash
# 1. Clone and set up
git clone ...
cd ...

# 2. Place source data
#    data/source/processed rest scans/      (MDD rest1 BRIK/HEAD)
#    data/source/processed rest2 scans/     (MDD rest2 BRIK/HEAD)
#    data/source/participants.tsv            (group assignments — CRITICAL)
#    data/HNU1/0025427..0025456/             (HC NIfTI, 10 sessions each)
#    atlases/Tian_Subcortex_S1_3T_2009cAsym.nii.gz

# 3. Parcellate HC data (~5 hours)
jupyter execute parcellate_hc_hnu1_v3.ipynb

# 4. Run cross-condition analysis (~18 hours, 71,928 UKF fits)
jupyter execute MDD-HC_analysis_v4.ipynb

# Results → results/
```

---

## Citation

If you use this pipeline or build on this work, please cite:

---

<div align="center">

*Built with the [Stuart-Landau](https://en.wikipedia.org/wiki/Stuart%E2%80%93Landau_equation) normal form · [Unscented Kalman Filter](https://github.com/insilico/UKF) · [MDD-fMRI_with_SL-UKF](https://github.com/skaraoglu/MDD-fMRI_with_SL-UKF) · [Schaefer 2018](https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal) + [Melbourne Subcortex](https://github.com/yetianmed/subcortex) · [HNU1/CoRR](https://fcon_1000.projects.nitrc.org/indi/CoRR/html/hnu_1.html)*

</div>
