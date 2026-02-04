# SC CTSI – Diabetes Prevention Program Management and Training
## Linear Mixed-Effects Model (LMM) Analysis for Knowledge Scores

This project contains the analytic workflow used to evaluate changes in diabetes knowledge following a structured training program, with a focus on whether learning and retention differed across sequential implementation phases. The primary analysis uses a **linear mixed-effects model (LMM)** to account for repeated measurements within participants and unbalanced follow-up.

This folder is intended to support **manuscript development** and **reproducible analysis**, with manuscript-ready figures, tables, and diagnostics generated directly from the R script.

---

## Repository Navigation

- **R analysis script:** `Projects/DPPM/LMM.R`  
  https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/LMM.R

- **De-identified long-format dataset:** `Projects/DPPM/Data/LMM_long.csv`  
  https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/Data/LMM_long.csv

- **Results (figures/tables):** `Projects/DPPM/Results/`  
  https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/Results  
  Current output:
  - `LMM_Results.png`  
    https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/Results/LMM_Results.png

---

## Project Structure (Quick View)

- `Projects/DPPM/`
  - `LMM.R` — end-to-end analysis script (data prep → modeling → figures/tables)
  - `Data/`
    - `LMM_long.csv` — de-identified long-format participant dataset
  - `Results/`
    - `LMM_Results.png` — exported summary figure (more outputs forthcoming)
  - `README.md` — project documentation

---

## Data Format

The dataset is stored in **long format**, where **each row represents one participant at one timepoint**.

**Expected columns:**
- `ID` — participant identifier (de-identified)
- `Phase` — implementation phase (`Phase 1`, `Phase 2`, `Phase 3`)
- `Timepoint` — assessment timing (`pretraining`, `posttraining`, `followup`)
- `Score` — knowledge score as **proportion correct** (0–1)

---

## Analytic Overview

### Study Design
- Longitudinal repeated-measures design
- Knowledge assessed at:
  - **Pre-training**
  - **Post-training**
  - **1-Month Follow-up**
- Participants enrolled across three sequential **implementation phases**

### Why a Linear Mixed-Effects Model (LMM)?
The LMM was selected because it:
- Accounts for **correlated repeated measures** within participants
- Accommodates **unbalanced data** and **missing follow-up** without dropping participants
- Supports phase-specific comparisons via a **Timepoint × Phase interaction**
- Produces interpretable, manuscript-ready estimates (e.g., adjusted means and contrasts)

---

## Model Specification

The primary analysis uses a linear mixed-effects model with fixed effects for **timepoint**, **implementation phase**, and their interaction, and a **participant-level random intercept** to account for repeated measurements.

**Model form:**

Score_ij = β0 + β1(Post-training)ij + β2(1-Month Follow-up)ij + β3(Phase 2)ij + β4(Phase 3)ij + β5(Post-training × Phase 2)ij + β6(Follow-up × Phase 2)ij + β7(Post-training × Phase 3)ij + β8(Follow-up × Phase 3)ij + u_i + ε_ij  

where:

- **i** indexes participants and **j** indexes timepoints  
- **β0** is the mean baseline knowledge score for Phase 1 at pre-training  
- **β1–β2** represent within-phase changes over time (Phase 1 reference)  
- **β3–β4** represent baseline differences by phase  
- **β5–β8** represent phase-specific differences in knowledge change over time  
- **u_i** is a participant-level random intercept (assumed ~ N(0, σ²_u))  
- **ε_ij** is the residual error term (assumed ~ N(0, σ²))

**Reference categories (for interpretation):**
- Reference phase: **Phase 1**
- Reference timepoint: **Pre-training**

---

## What `LMM.R` Produces

### 1) Descriptive Visualizations
- Frequency distributions of knowledge scores:
  - **Pooled by timepoint** (across phases)
  - **Stratified by Phase × Timepoint** (counts with free y-scale for readability)
- Optional density overlays for exploratory distribution checks

### 2) Model Diagnostics (Assumption Checks)
- Residuals vs. fitted values
- Normal Q–Q plot of residuals

### 3) Model-Adjusted Estimates
- Estimated marginal means (LS-means) for timepoint **within each phase**
- 95% confidence intervals
- Manuscript-ready profile plot of adjusted means over time

### 4) SAS-Style Fixed Effects Table
A “Solution for Fixed Effects”-style output table including:
- β coefficients
- Standard errors
- Degrees of freedom (Satterthwaite)
- *p*-values
- 95% confidence intervals

---

## Reproducibility Notes

- Factor reference levels are explicitly set to support consistent interpretation:
  - `Timepoint`: `pretraining` is the reference
  - `Phase`: `Phase 1` is the reference
- Manuscript labels are defined once and reused consistently for plots
- Missing observations are handled using all available data (consistent with LMM assumptions)

---

## Intended Use

This project is designed to:
- Support peer-reviewed manuscript preparation
- Provide transparent documentation of analytic decisions
- Serve as a reusable template for similar longitudinal training evaluations
- Facilitate collaboration across evaluation, research, and implementation teams

---

## Contributions / Updates

Additional exported figures and tables will be added as manuscript development progresses. Suggestions and extensions are welcome.