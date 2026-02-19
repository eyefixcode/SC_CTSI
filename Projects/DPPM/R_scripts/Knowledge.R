###############################################
# Linear Mixed-Effects Model (LMM) Analysis
# Diabetes Management Training – Knowledge Scores
#
# Goal (manuscript-ready workflow; consistent with Confidence analysis):
#   1) Describe knowledge score distributions over time (Pre, Post, Follow-up), overall and by Phase
#   2) Fit an LMM to test whether change over time differs by implementation Phase
#   3) Produce publishable outputs:
#        - Missingness summary
#        - Descriptive plots (box + jitter)
#        - Model diagnostics (residuals vs fitted; Q–Q)
#        - EMM (model-adjusted means) plot + EMM table
#        - Fixed-effects tables (full, reduced, selected)
#        - Forest plot of fixed effects (APA/clinical style)
#
# Key difference vs Confidence:
#   - Knowledge has three timepoints (Pre, Post, Follow-up), so the model includes follow-up terms.
#
# Reference categories (interpretation):
#   - Timepoint reference = Pre-training
#   - Phase reference     = Phase 1
###############################################

# -----------------------------
# 0) Setup
# -----------------------------
setwd("/Users/bryce/Desktop/SC_CTSI/Diabetes Prevention Management")

# -----------------------------
# 1) Libraries
# -----------------------------
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

# Quick integrity checks (QA)
head(lmm_data)
str(lmm_data)

# -----------------------------
# 3) Standardize factor coding + manuscript labels ONCE
# -----------------------------
time_levels_know <- c("pretraining", "posttraining", "followup")
time_labels_know <- c("Pre-training", "Post-training", "1-Month Follow-up")

phase_levels <- c("Phase 1", "Phase 2", "Phase 3")
phase_labels <- c("Phase 1", "Phase 2", "Phase 3")

lmm_data <- lmm_data %>%
  mutate(
    ID        = factor(ID),
    Timepoint = factor(Timepoint, levels = time_levels_know),  # Pre-training reference
    Phase     = factor(Phase, levels = phase_levels),          # Phase 1 reference
    Score     = as.numeric(Score),
    
    # Plot-only labels
    Timepoint_label = factor(Timepoint, levels = time_levels_know, labels = time_labels_know),
    Phase_label     = factor(Phase,     levels = phase_levels,    labels = phase_labels)
  )

plot_data_know <- lmm_data %>% filter(!is.na(Score))

# -----------------------------
# 4) Missingness check (descriptive)
# -----------------------------
table(lmm_data$Timepoint, is.na(lmm_data$Score))

# -----------------------------
# 5) Fit candidate LMMs and test interaction (ML; not REML)
# -----------------------------
# Use ML (REML = FALSE) to compare models that differ in fixed effects.
lmm_full <- lmer(
  Score ~ Timepoint * Phase + (1 | ID),
  data = lmm_data,
  na.action = na.omit,
  REML = FALSE
)

lmm_reduced <- lmer(
  Score ~ Timepoint + Phase + (1 | ID),
  data = lmm_data,
  na.action = na.omit,
  REML = FALSE
)

# Likelihood ratio test (LRT) for interaction contribution
interaction_test_knowledge <- anova(lmm_reduced, lmm_full)
interaction_test_knowledge

# Select final model based on interaction LRT p-value
p_int_know <- interaction_test_knowledge$`Pr(>Chisq)`[2]
lmm <- if (!is.na(p_int_know) && p_int_know < 0.05) lmm_full else lmm_reduced

cat(
  "\nFinal model selected:",
  if (identical(lmm, lmm_full)) "FULL (with interaction)" else "REDUCED (no interaction)",
  "\nInteraction LRT p-value:", p_int_know, "\n\n"
)

summary(lmm)

# -----------------------------
# 5b) Main-effects-only model (no interaction)
# -----------------------------
# This model provides estimates assuming no Timepoint × Phase interaction.
# Useful for reporting “unadjusted for interaction” effects.

summary(lmm_reduced)

fixef_tbl_main_effects <- broom.mixed::tidy(
  lmm_reduced,
  effects = "fixed",
  conf.int = TRUE
) %>%
  transmute(
    Effect   = term,
    Estimate = round(estimate, 3),
    SE       = round(std.error, 3),
    DF       = round(df, 3),
    tValue   = round(statistic, 3),
    pValue   = round(p.value, 3),
    CI_Lower = round(conf.low, 3),
    CI_Upper = round(conf.high, 3)
  )

fixef_tbl_main_effects

# -----------------------------
# 6) Descriptive distribution plot (box + jitter)
# -----------------------------
# Unadjusted visualization: shows medians/IQR and individual values by phase and timepoint.
# Tweaks for manuscript style:
#   - Reduced point size and alpha to limit overplotting (especially Phase 2).
#   - Keep boxplot prominent.
p_know_box <- ggplot(plot_data_know, aes(x = Timepoint_label, y = Score)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 1) +
  geom_jitter(width = 0.10, alpha = 0.35, size = 1.4) +
  facet_wrap(~ Phase_label, nrow = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Knowledge Score Distributions by Timepoint and Phase",
    x = "",
    y = "Knowledge Score (% Correct)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_know_box

# -----------------------------
# 7) Model diagnostics (assumptions checks)
# -----------------------------
# Replace base diagnostic scatter with ggplot version for consistent export and readability.
diag_df_know <- data.frame(
  fitted = fitted(lmm),
  resid  = resid(lmm, type = "pearson")
)

p_know_resid <- ggplot(diag_df_know, aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_point(size = 1.6, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  labs(
    title = "Residuals vs Fitted (Knowledge Model)",
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_know_resid

# Q–Q plot (base R acceptable; title clarified for manuscript consistency)
qqnorm(resid(lmm), main = "Normal Q–Q Plot of Residuals (Knowledge Model)")
qqline(resid(lmm))

# -----------------------------
# 8) Estimated Marginal Means (EMMs) + plot (model-adjusted)
# -----------------------------
emm_know <- emmeans(lmm, ~ Timepoint | Phase)

emm_know_df <- as.data.frame(emm_know) %>%
  mutate(
    Timepoint_label = factor(Timepoint, levels = time_levels_know, labels = time_labels_know),
    Phase_label     = factor(Phase,     levels = phase_levels,     labels = phase_labels)
  )

emm_know_df

# Manuscript-style EMM plot:
#   - Avoid reliance on color; use linetype + shape with black ink.
#   - Classic theme; legend moved to bottom.
p_know_emm <- ggplot(
  emm_know_df,
  aes(
    x = Timepoint_label, y = emmean,
    group = Phase_label,
    linetype = Phase_label,
    shape = Phase_label
  )
) +
  geom_line(linewidth = 0.9, color = "black") +
  geom_point(size = 2.6, color = "black") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.7, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title = "Estimated Knowledge Scores Over Time by Phase (Model-Adjusted)",
    x = "Timepoint",
    y = "Estimated Mean Knowledge Score (% Correct)",
    linetype = "Phase",
    shape = "Phase"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p_know_emm

# -----------------------------
# 9) Fixed effects tables (FULL, REDUCED, SELECTED)
# -----------------------------
know_fixef_tbl_full <- broom.mixed::tidy(lmm_full, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Effect   = term,
    Estimate = estimate,
    SE       = std.error,
    DF       = df,
    tValue   = statistic,
    pValue   = p.value,
    CI_Lower = conf.low,
    CI_Upper = conf.high
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

know_fixef_tbl_reduced <- broom.mixed::tidy(lmm_reduced, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Effect   = term,
    Estimate = estimate,
    SE       = std.error,
    DF       = df,
    tValue   = statistic,
    pValue   = p.value,
    CI_Lower = conf.low,
    CI_Upper = conf.high
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

know_fixef_tbl_selected <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Effect   = term,
    Estimate = estimate,
    SE       = std.error,
    DF       = df,
    tValue   = statistic,
    pValue   = p.value,
    CI_Lower = conf.low,
    CI_Upper = conf.high
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

know_fixef_tbl_full
know_fixef_tbl_reduced
know_fixef_tbl_selected

# -----------------------------
# 10) Forest plot of fixed effects (APA/clinical style)
# -----------------------------
# Key change for consistency:
#   - Build the forest plot from the SELECTED model (lmm), not always the full model.
#   - Term mapping adapts automatically depending on whether interaction terms exist.
pp <- 100  # express effects as percentage points

forest_data_know <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    estimate  = estimate  * pp,
    conf.low  = conf.low  * pp,
    conf.high = conf.high * pp,
    
    term_clean = case_when(
      term == "PhasePhase 2" ~ "Phase 2 vs Phase 1 (reference)",
      term == "PhasePhase 3" ~ "Phase 3 vs Phase 1 (reference)",
      term == "Timepointposttraining" ~ "Post-training vs Pre-training",
      term == "Timepointfollowup"     ~ "Follow-up vs Pre-training",
      term == "Timepointposttraining:PhasePhase 2" ~ "Post-training × Phase 2",
      term == "Timepointposttraining:PhasePhase 3" ~ "Post-training × Phase 3",
      term == "Timepointfollowup:PhasePhase 2"     ~ "Follow-up × Phase 2",
      term == "Timepointfollowup:PhasePhase 3"     ~ "Follow-up × Phase 3",
      TRUE ~ term
    ),
    
    block = ifelse(grepl("×", term_clean), "Interaction terms", "Main effects"),
    block = factor(block, levels = c("Main effects", "Interaction terms"))
  )

# Order only the terms that are present (prevents errors if reduced model was selected)
preferred_order_know <- c(
  "Phase 2 vs Phase 1 (reference)",
  "Phase 3 vs Phase 1 (reference)",
  "Post-training vs Pre-training",
  "Follow-up vs Pre-training",
  "Post-training × Phase 2",
  "Post-training × Phase 3",
  "Follow-up × Phase 2",
  "Follow-up × Phase 3"
)

present_terms <- preferred_order_know[preferred_order_know %in% forest_data_know$term_clean]
forest_data_know <- forest_data_know %>%
  mutate(term_clean = factor(term_clean, levels = rev(present_terms)))

# Symmetric x-limits around 0; clean tick spacing
lim_know <- max(abs(c(forest_data_know$conf.low, forest_data_know$conf.high)), na.rm = TRUE)
lim_know <- ceiling(lim_know / 10) * 10

p_know_forest <- ggplot(forest_data_know, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(
    limits = c(-lim_know, lim_know),
    breaks = seq(-lim_know, lim_know, by = 10)
  ) +
  labs(
    title = "Fixed-Effects Estimates for Knowledge Scores",
    x = "Estimate (percentage-point change in knowledge score)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.text.y = element_text(size = 11, lineheight = 1.1),
    axis.text.x = element_text(size = 11),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(face = "plain", size = 9),
    plot.margin = margin(10, 30, 10, 10)
  )

p_know_forest

# -----------------------------
# Optional exports (tables)
# -----------------------------
# write.csv(know_fixef_tbl_full,     "Knowledge_FixedEffectsTable_FULL.csv",     row.names = FALSE)
# write.csv(know_fixef_tbl_reduced,  "Knowledge_FixedEffectsTable_REDUCED.csv",  row.names = FALSE)
# write.csv(know_fixef_tbl_selected, "Knowledge_FixedEffectsTable_SELECTED.csv", row.names = FALSE)

# -----------------------------
# 11) Export key plots to PNG (standardized)
# -----------------------------
fig_dir <- "figures_png_knowledge"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

plots_to_save <- list(
  "knowledge_box_by_phase_time"    = p_know_box,
  "knowledge_residuals_vs_fitted"  = p_know_resid,
  "knowledge_emm_by_phase"         = p_know_emm,
  "knowledge_forest_fixed_effects" = p_know_forest
)

for (nm in names(plots_to_save)) {
  ggsave(
    filename = file.path(fig_dir, paste0(nm, ".png")),
    plot     = plots_to_save[[nm]],
    device   = "png",
    width    = 7,
    height   = 4.5,
    units    = "in",
    dpi      = 300,
    bg       = "white"
  )
}

message("Saved ", length(plots_to_save), " plots to: ", normalizePath(fig_dir))
###############################################
