###############################################
# Linear Mixed-Effects Model (LMM) Analysis
# Diabetes Management Training – Knowledge Scores
#
# Purpose:
#   1) Describe and visualize knowledge score distributions (pre, post, follow-up)
#   2) Fit an LMM to test whether change over time differs by implementation Phase
#   3) Produce manuscript-ready outputs (diagnostics, EMM plot, fixed-effects table)
#
# Model (conceptual):
#   Score_ij = β0 + β1(Post) + β2(Follow) + β3(Phase2) + β4(Phase3)
#            + β5(Post×Phase2) + β6(Follow×Phase2)
#            + β7(Post×Phase3) + β8(Follow×Phase3)
#            + u_i + ε_ij
#   where u_i is a participant-specific random intercept.
###############################################

# -----------------------------
# 0) Setup
# -----------------------------
# Working directory: set to your local folder containing LMM_long.csv
# Tip: if sharing scripts, consider using here::here() to avoid hardcoding.
setwd("/Users/bryce/Desktop/SC_CTSI/Diabetes Prevention Management")

# -----------------------------
# 1) Libraries
# -----------------------------
# lme4       : fits linear mixed-effects models
# lmerTest   : adds Satterthwaite df + p-values for fixed effects (SAS-like)
# dplyr      : data manipulation
# readr      : CSV import
# emmeans    : estimated marginal means (LS-means) + contrasts
# ggplot2    : plotting
# scales     : percent formatting
# broom.mixed: tidy model output to tables
library(lme4)
library(lmerTest)
library(dplyr)
library(readr)
library(emmeans)
library(ggplot2)
library(scales)
library(broom.mixed)

# -----------------------------
# 2) Read data (long format)
# -----------------------------
# Expected columns:
#   ID        : participant identifier
#   Phase     : Phase 1 / Phase 2 / Phase 3 (between-subject factor)
#   Timepoint : pretraining / posttraining / followup (within-subject factor)
#   Score     : knowledge score as proportion correct (0–1)
lmm_data <- read_csv("LMM_long.csv")

# Quick integrity checks
head(lmm_data)
str(lmm_data)

# -----------------------------
# 3) Standardize factor coding + set manuscript labels ONCE
# -----------------------------
# IMPORTANT: Factor levels below define reference categories for interpretation:
#   - Timepoint reference = pretraining
#   - Phase reference     = Phase 1
#
# Also create one set of manuscript-ready labels to use consistently across ALL plots.
time_levels  <- c("pretraining", "posttraining", "followup")
time_labels  <- c("Pre-training", "Post-training", "1-Month Follow-up")

phase_levels <- c("Phase 1", "Phase 2", "Phase 3")
phase_labels <- c("Phase 1", "Phase 2", "Phase 3")

lmm_data <- lmm_data %>%
  mutate(
    ID = factor(ID),
    Timepoint = factor(Timepoint, levels = time_levels),
    Phase     = factor(Phase, levels = phase_levels),
    Score     = as.numeric(Score),
    
    # Manuscript-friendly display labels (use these for plotting only)
    Timepoint_label = factor(Timepoint, levels = time_levels, labels = time_labels),
    Phase_label     = factor(Phase,     levels = phase_levels, labels = phase_labels)
  )

# Convenience object for plotting (drop missing scores so ggplot doesn't warn)
plot_data <- lmm_data %>% filter(!is.na(Score))

# -----------------------------
# 4) Missingness check (descriptive)
# -----------------------------
# Shows attrition / missing data by timepoint.
table(lmm_data$Timepoint, is.na(lmm_data$Score))

# Optional: number of observations per participant (how many timepoints each contributes)
# lmm_data %>% count(ID) %>% count(n)

# -----------------------------
# 5) Fit the Linear Mixed-Effects Model
# -----------------------------
# Fixed effects:
#   Timepoint + Phase + Timepoint×Phase
# Random effects:
#   (1 | ID) participant-specific intercept (accounts for repeated measures)
#
# Missing data:
#   na.omit removes rows with missing Score; participants can still contribute partial data.
lmm <- lmer(
  Score ~ Timepoint * Phase + (1 | ID),
  data = lmm_data,
  na.action = na.omit
)

# Model summary with df + p-values (via lmerTest)
summary(lmm)

# -----------------------------
# 6) Descriptive distribution plots (manuscript-oriented)
# -----------------------------

# 6a) Pooled across phases: distribution by timepoint (3 panels; counts)
# Purpose: simple “orientation” figure for readers (like Dupuis et al.).
ggplot(plot_data, aes(x = Score)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "grey85",
    color = "black"
  ) +
  facet_wrap(~ Timepoint_label, nrow = 1) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Distribution of Knowledge Scores by Timepoint (Pooled Across Phases)",
    x = "Knowledge Score (% Correct)",
    y = "Frequency"
  ) +
  theme_classic(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))

# 6b) Stratified by phase and timepoint: counts with free y-scale
# Purpose: show distribution shape by phase/timepoint without Phase 2 dominating visually.
# NOTE: y-axis varies by panel to accommodate different sample sizes.
ggplot(plot_data, aes(x = Score)) +
  geom_histogram(
    binwidth = 0.05,
    fill = "grey85",
    color = "black"
  ) +
  facet_grid(Phase_label ~ Timepoint_label, scales = "free_y") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "Frequency Distributions of Knowledge Scores by Phase and Timepoint",
    x = "Knowledge Score (% Correct)",
    y = "Number of Participants"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

# 6c) Optional: density overlays (diagnostic / exploratory; usually Supplement-only)
# Purpose: compare distribution shapes (not counts) within each timepoint.
ggplot(plot_data, aes(x = Score, fill = Phase_label)) +
  geom_histogram(aes(y = after_stat(density)),
                 binwidth = 0.05, alpha = 0.5, position = "identity") +
  geom_density(color = "black", linewidth = 1) +
  facet_wrap(~ Timepoint_label) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Density Distributions of Knowledge Scores by Phase and Timepoint",
    x = "Knowledge Score (% Correct)",
    y = "Density",
    fill = "Phase"
  ) +
  theme_minimal(base_size = 13)

# -----------------------------
# 7) Model diagnostics (assumptions checks)
# -----------------------------
# 7a) Residuals vs fitted: check for obvious patterns / heteroscedasticity
plot(lmm, which = 1)

# 7b) Normal Q–Q plot of residuals: check approximate normality
qqnorm(resid(lmm))
qqline(resid(lmm))

# -----------------------------
# 8) Estimated Marginal Means (LS-means) + plot (model-adjusted)
# -----------------------------
# EMMs are a clean way to visualize and report adjusted means per timepoint within each phase.
emm <- emmeans(lmm, ~ Timepoint | Phase)

emm_df <- as.data.frame(emm) %>%
  mutate(
    # Apply manuscript labels in the output dataframe so plots read cleanly
    Timepoint_label = factor(Timepoint, levels = time_levels, labels = time_labels),
    Phase_label     = factor(Phase,     levels = phase_levels, labels = phase_labels)
  )

head(emm_df)

ggplot(emm_df, aes(x = Timepoint_label, y = emmean, group = Phase_label, color = Phase_label)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Estimated Knowledge Scores Over Time by Phase (Model-Adjusted)",
    x = "Timepoint",
    y = "Estimated Mean Knowledge Score (% Correct)",
    color = "Phase"
  ) +
  theme_minimal(base_size = 14)

# Useful for exporting or inspecting the EMM table
emm_df

# -----------------------------
# 9) SAS-style fixed effects table (Solution for Fixed Effects)
# -----------------------------
# This mirrors PROC MIXED-style reporting: Estimate, SE, df, t, p, and 95% CI.
fixef_tbl <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Effect   = term,
    Estimate = estimate,
    SE       = std.error,
    DF       = df,
    tValue   = statistic,
    pValue   = p.value,
    CI_Lower = conf.low,
    CI_Upper = conf.high
  )

# Rounded version for manuscript tables
fixef_tbl_fmt <- fixef_tbl %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

fixef_tbl_fmt

# -----------------------------
# Optional next steps (if needed)
# -----------------------------
# 1) Type III tests (SAS-like):
# anova(lmm, type = 3)
#
# 2) Planned contrasts within phase (Pre→Post, Pre→Follow, Post→Follow):
# emm_within <- emmeans(lmm, ~ Timepoint | Phase)
# pairs(emm_within, adjust = "tukey")
#
# 3) Export fixed effects table:
# write.csv(fixef_tbl_fmt, "FixedEffectsTable.csv", row.names = FALSE)
###############################################