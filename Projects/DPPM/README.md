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

Score ij = β0 + β1(Post-training)ij + β2(1-Month Follow-up)ij + β3(Phase 2)ij + β4(Phase 3)ij + β5(Post-training × Phase 2)ij + β6(Follow-up × Phase 2)ij + β7(Post-training × Phase 3)ij + β8(Follow-up × Phase 3)ij + ui + εij  

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

## Results Outputs

📂 Results directory:  
https://github.com/eyefixcode/SC_CTSI/tree/main/Projects/DPPM/Results

### Descriptive Distributions

- **`Histogram_Freq_By_Phase_Timepoint.png`**  
  Frequency histograms of knowledge scores stratified by **Phase × Timepoint**.  
  Panels use free y-axis scaling to accommodate unequal sample sizes and highlight distributional shape, ceiling effects, and phase-specific patterns.

- **`Histogram_Pooled_By_Phase.png`**  
  Knowledge score distributions pooled across timepoints and stratified by phase.  
  Intended as a high-level descriptive summary of overall score distributions by implementation phase.

- **`Density_Distribution_PhaseTimepoint.png`**  
  Density-based representations of score distributions by phase and timepoint.  
  These plots emphasize distributional shape rather than counts and are primarily intended for exploratory or supplemental use.

---

### Model-Adjusted Results

- **`KnowledgeScoreOverTimebyPhase.png`**  
  Estimated marginal means (LS-means) of knowledge scores over time within each phase, derived from the linear mixed-effects model.  
  Points represent adjusted means and error bars indicate 95% confidence intervals.  
  This figure is the primary visualization supporting phase-specific learning and retention effects.

- **`LMM_Results.png`**  
  Summary visualization of fixed-effect estimates from the linear mixed-effects model.  
  Designed to align directly with manuscript tables and reported coefficients.

---

### Model Diagnostics

- **`ResidualPlot.png`**  
  Residuals versus fitted values plot used to assess linearity and homoscedasticity assumptions of the linear mixed-effects model.

- **`QQ_plot.png`**  
  Normal Q–Q plot of model residuals used to assess approximate normality.

---

### Notes

- All figures are generated directly from the analysis script (`LMM.R`) to ensure reproducibility.

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