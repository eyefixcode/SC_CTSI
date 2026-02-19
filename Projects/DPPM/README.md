# SC CTSI – Diabetes Prevention Program Management and Training
## Linear Mixed-Effects Model (LMM) Analysis for Knowledge and Confidence Outcomes

This project contains the analytic workflow used to evaluate changes in diabetes knowledge and self-efficacy (confidence) following a structured training program, with a focus on whether learning and retention differed across sequential implementation phases. The primary analyses use a **linear mixed-effects model (LMM)** to account for repeated measurements within participants and unbalanced follow-up.

This folder is intended to support **manuscript development** and **reproducible analysis**, with manuscript-ready figures, tables, and diagnostics generated directly from the R scripts.

---

## Repository Navigation

### R Analysis Scripts

- **Knowledge analysis script:** `Projects/DPPM/R_scripts/Knowledge.R`  
  https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/R_scripts/Knowledge.R

- **Confidence analysis script:** `Projects/DPPM/R_scripts/Confidence.R`  
  https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/R_scripts/Confidence.R

---

### Data

- **Knowledge dataset (long format):**  
  `Projects/DPPM/Data/LMM_long.csv`

- **Confidence dataset (long format):**  
  `Projects/DPPM/Data/Confidence_long.csv`

---

### Results (Figures)

📂 Results directory:  
https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/Results

#### Knowledge Outputs
- `knowledge_box_by_phase_time.png`
- `knowledge_emm_by_phase.png`
- `knowledge_forest_fixed_effects.png`
- `knowledge_residuals_vs_fitted.png`

#### Confidence Outputs
- `confidence_box_by_phase_time.png`
- `confidence_emm_by_phase.png`
- `confidence_forest_fixed_effects.png`
- `confidence_residuals_vs_fitted.png`

All figures are generated directly from their respective R scripts to ensure reproducibility.

---

## Project Structure (Quick View)

- `Projects/DPPM/`
  - `R_scripts/`
    - `Knowledge.R` — end-to-end knowledge analysis (data prep → modeling → figures/tables)
    - `Confidence.R` — parallel analysis pipeline for confidence outcomes
  - `Data/`
    - `LMM_long.csv`
    - `Confidence_long.csv`
  - `Results/`
    - Exported PNG figures (knowledge and confidence)
  - `README.md` — project documentation

---

## Data Format

The datasets are stored in **long format**, where **each row represents one participant at one timepoint**.

### Knowledge Dataset Columns
- `ID` — participant identifier (de-identified)
- `Phase` — implementation phase (`Phase 1`, `Phase 2`, `Phase 3`)
- `Timepoint` — assessment timing (`pretraining`, `posttraining`, `followup`)
- `Score` — knowledge score as **proportion correct** (0–1)

### Confidence Dataset Columns
- `ID` — participant identifier
- `Phase` — implementation phase
- `Timepoint` — assessment timing (`pretraining`, `posttraining`)
- `Confidence` — confidence score (0–10 scale)

---

## Analytic Overview

### Study Design
- Longitudinal repeated-measures design
- Participants enrolled across three sequential **implementation phases**

**Knowledge assessed at:**
- Pre-training
- Post-training
- 1-Month Follow-up

**Confidence assessed at:**
- Pre-training
- Post-training

---

## Why a Linear Mixed-Effects Model (LMM)?

The LMM was selected because it:

- Accounts for **correlated repeated measures** within participants
- Accommodates **unbalanced data** and **missing follow-up** without dropping participants
- Supports phase-specific comparisons via a **Timepoint × Phase interaction**
- Produces interpretable, manuscript-ready estimates (e.g., adjusted means and contrasts)

Models were fit using the `lme4` package in R, with Satterthwaite-approximated degrees of freedom and p-values provided by `lmerTest`.

---

## Model Specification – Knowledge

The primary knowledge analysis uses a linear mixed-effects model with fixed effects for **timepoint**, **implementation phase**, and their interaction, and a **participant-level random intercept** to account for repeated measurements.

**Model form:**

Score ij = β0 + β1(Post-training)ij + β2(1-Month Follow-up)ij + β3(Phase 2)ij + β4(Phase 3)ij + β5(Post-training × Phase 2)ij + β6(Follow-up × Phase 2)ij + β7(Post-training × Phase 3)ij + β8(Follow-up × Phase 3)ij + ui + εij  

where:

- **i** indexes participants and **j** indexes timepoints  
- **β0** is the mean baseline knowledge score for Phase 1 at pre-training  
- **β1–β2** represent within-phase changes over time (Phase 1 reference)  
- **β3–β4** represent baseline differences by phase  
- **β5–β8** represent phase-specific differences in knowledge change over time  
- **ui** is a participant-level random intercept (assumed ~ N(0, σ²_u))  
- **εij** is the residual error term (assumed ~ N(0, σ²))

**Reference categories (for interpretation):**
- Reference phase: **Phase 1**
- Reference timepoint: **Pre-training**

---

## Model Specification – Confidence

The confidence model follows the same framework but includes only two timepoints (Pre and Post).

**Model form:**

Confidence ij = β0 + β1(Post)ij + β2(Phase 2)ij + β3(Phase 3)ij + β4(Post × Phase 2)ij + β5(Post × Phase 3)ij + ui + εij  

**Reference categories (for interpretation):**
- Reference phase: **Phase 1**
- Reference timepoint: **Pre-training**

---

## What the Scripts Produce

### 1) Descriptive Visualizations
- Boxplot + jittered individual values by Phase × Timepoint

### 2) Model Diagnostics (Assumption Checks)
- Residuals vs. fitted values
- Normal Q–Q plot of residuals

### 3) Model-Adjusted Estimates
- Estimated marginal means (LS-means) for timepoint within each phase
- 95% confidence intervals
- Manuscript-ready profile plot of adjusted means over time

### 4) Fixed Effects Tables
- β coefficients
- Standard errors
- Degrees of freedom (Satterthwaite)
- p-values
- 95% confidence intervals
- Forest plot of fixed-effect estimates

---

## Intended Use

This project is designed to:
- Support peer-reviewed manuscript preparation
- Provide transparent documentation of analytic decisions
- Serve as a reusable template for similar longitudinal training evaluations
- Facilitate collaboration across evaluation, research, and implementation teams

---

## Contributions / Updates

All figures are generated directly from the analysis scripts (`Knowledge.R` and `Confidence.R`) to ensure reproducibility and analytic transparency.

Additional exported figures and refinements may be added as manuscript development progresses.