###############################################
# Linear Mixed-Effects Model (LMM) Analysis
# Diabetes Management Training – Confidence Scores
#
# Goal (manuscript-ready workflow; consistent with Knowledge analysis):
#   1) Describe confidence score distributions over time (Pre, Post), overall and by Phase
#   2) Fit an LMM to test whether change over time differs by implementation Phase
#   3) Produce publishable outputs:
#        - Missingness summary
#        - Descriptive plots (box + jitter)
#        - Model diagnostics (residuals vs fitted; Q–Q)
#        - EMM (model-adjusted means) plot + EMM table
#        - Fixed-effects tables (full, reduced, selected)
#        - Forest plot of fixed effects (APA/clinical style)
#
# Key difference vs Knowledge:
#   - Confidence has only two timepoints (Pre, Post), so the model has fewer time terms.
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
library(broom.mixed)

# -----------------------------
# 2) Read data (long format)
# -----------------------------
# Expected columns:
#   ID        : participant identifier
#   Phase     : Phase 1 / Phase 2 / Phase 3 (between-subject factor)
#   Timepoint : pretraining / posttraining (within-subject factor)
#   Confidence: confidence score (e.g., 0–10 scale)
conf_data <- read_csv("Confidence_long.csv")

# Quick integrity checks (QA)
head(conf_data)
str(conf_data)

# -----------------------------
# 3) Standardize factor coding + manuscript labels ONCE
# -----------------------------
time_levels_conf <- c("pretraining", "posttraining")
time_labels_conf <- c("Pre-training", "Post-training")

phase_levels <- c("Phase 1", "Phase 2", "Phase 3")
phase_labels <- c("Phase 1", "Phase 2", "Phase 3")

conf_data <- conf_data %>%
  mutate(
    ID         = factor(ID),
    Timepoint  = factor(Timepoint, levels = time_levels_conf), # Pre-training reference
    Phase      = factor(Phase, levels = phase_levels),         # Phase 1 reference
    Confidence = as.numeric(Confidence),
    
    # Plot-only labels
    Timepoint_label = factor(Timepoint, levels = time_levels_conf, labels = time_labels_conf),
    Phase_label     = factor(Phase,     levels = phase_levels,    labels = phase_labels)
  )

plot_data_conf <- conf_data %>% filter(!is.na(Confidence))

# -----------------------------
# 4) Missingness check (descriptive)
# -----------------------------
table(conf_data$Timepoint, is.na(conf_data$Confidence))

# -----------------------------
# 5) Fit candidate LMMs and test interaction (ML; not REML)
# -----------------------------
# Use ML (REML = FALSE) to compare models that differ in fixed effects.
conf_lmm_full <- lmer(
  Confidence ~ Timepoint * Phase + (1 | ID),
  data = conf_data,
  na.action = na.omit,
  REML = FALSE
)

conf_lmm_reduced <- lmer(
  Confidence ~ Timepoint + Phase + (1 | ID),
  data = conf_data,
  na.action = na.omit,
  REML = FALSE
)

# Likelihood ratio test (LRT) for interaction contribution
interaction_test_confidence <- anova(conf_lmm_reduced, conf_lmm_full)
interaction_test_confidence

# Select final model based on interaction LRT p-value
p_int_conf <- interaction_test_confidence$`Pr(>Chisq)`[2]
conf_lmm <- if (!is.na(p_int_conf) && p_int_conf < 0.05) conf_lmm_full else conf_lmm_reduced

cat(
  "\nFinal model selected:",
  if (identical(conf_lmm, conf_lmm_full)) "FULL (with interaction)" else "REDUCED (no interaction)",
  "\nInteraction LRT p-value:", p_int_conf, "\n\n"
)

summary(conf_lmm)

# -----------------------------
# 6) Descriptive distribution plot (box + jitter)
# -----------------------------
# Unadjusted visualization: shows medians/IQR and individual values by phase and timepoint.
# Tweaks for manuscript style:
#   - Reduced point size and alpha to limit overplotting (especially Phase 2).
#   - Keep boxplot prominent.
p_conf_box <- ggplot(plot_data_conf, aes(x = Timepoint_label, y = Confidence)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 1) +
  geom_jitter(width = 0.10, alpha = 0.35, size = 1.4) +
  facet_wrap(~ Phase_label, nrow = 1) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  labs(
    title = "Confidence Score Distributions by Timepoint and Phase",
    x = "",
    y = "Confidence Score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

p_conf_box

# -----------------------------
# 7) Model diagnostics (assumptions checks)
# -----------------------------
# Replace base diagnostic scatter with ggplot version for consistent export and readability.
diag_df_conf <- data.frame(
  fitted = fitted(conf_lmm),
  resid  = resid(conf_lmm, type = "pearson")
)

p_conf_resid <- ggplot(diag_df_conf, aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_point(size = 1.6, alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.8) +
  labs(
    title = "Residuals vs Fitted (Confidence Model)",
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_conf_resid

# Q–Q plot (base R is acceptable; title clarified for manuscript consistency)
qqnorm(resid(conf_lmm), main = "Normal Q–Q Plot of Residuals (Confidence Model)")
qqline(resid(conf_lmm))

# -----------------------------
# 8) Estimated Marginal Means (EMMs) + plot (model-adjusted)
# -----------------------------
emm_conf <- emmeans(conf_lmm, ~ Timepoint | Phase)

emm_conf_df <- as.data.frame(emm_conf) %>%
  mutate(
    Timepoint_label = factor(Timepoint, levels = time_levels_conf, labels = time_labels_conf),
    Phase_label     = factor(Phase,     levels = phase_levels,     labels = phase_labels)
  )

emm_conf_df

# Manuscript-style EMM plot:
#   - Avoid reliance on color; use linetype + shape with black ink.
#   - Classic theme; legend moved to bottom.
p_conf_emm <- ggplot(
  emm_conf_df,
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
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  labs(
    title = "Estimated Confidence Scores Over Time by Phase (Model-Adjusted)",
    x = "Timepoint",
    y = "Estimated Mean Confidence Score",
    linetype = "Phase",
    shape = "Phase"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p_conf_emm

# -----------------------------
# 9) Fixed effects tables (FULL, REDUCED, SELECTED)
# -----------------------------
conf_fixef_tbl_full <- broom.mixed::tidy(conf_lmm_full, effects = "fixed", conf.int = TRUE) %>%
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

conf_fixef_tbl_reduced <- broom.mixed::tidy(conf_lmm_reduced, effects = "fixed", conf.int = TRUE) %>%
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

conf_fixef_tbl_selected <- broom.mixed::tidy(conf_lmm, effects = "fixed", conf.int = TRUE) %>%
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

conf_fixef_tbl_full
conf_fixef_tbl_reduced
conf_fixef_tbl_selected

# -----------------------------
# 10) Forest plot of fixed effects (APA/clinical style)
# -----------------------------
# Key change for consistency:
#   - Build the forest plot from the SELECTED model (conf_lmm), not always the full model.
#   - Term mapping adapts automatically depending on whether interaction terms exist.

forest_data_conf <- broom.mixed::tidy(conf_lmm, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = case_when(
      term == "PhasePhase 2" ~ "Phase 2 vs Phase 1 (baseline)",
      term == "PhasePhase 3" ~ "Phase 3 vs Phase 1 (baseline)",
      term == "Timepointposttraining" ~ "Post vs Pre",
      term == "Timepointposttraining:PhasePhase 2" ~ "Post × Phase 2",
      term == "Timepointposttraining:PhasePhase 3" ~ "Post × Phase 3",
      TRUE ~ term
    ),
    block = ifelse(grepl("×", term_clean), "Interaction terms", "Main effects"),
    block = factor(block, levels = c("Main effects", "Interaction terms"))
  )

# Order only the terms that are present (prevents errors if reduced model was selected)
preferred_order_conf <- c(
  "Phase 2 vs Phase 1 (baseline)",
  "Phase 3 vs Phase 1 (baseline)",
  "Post vs Pre",
  "Post × Phase 2",
  "Post × Phase 3"
)

present_terms <- preferred_order_conf[preferred_order_conf %in% forest_data_conf$term_clean]
forest_data_conf <- forest_data_conf %>%
  mutate(term_clean = factor(term_clean, levels = rev(present_terms)))

# Symmetric x-limits around 0; clean tick spacing
lim_conf <- max(abs(c(forest_data_conf$conf.low, forest_data_conf$conf.high)), na.rm = TRUE)
lim_conf <- ceiling(lim_conf)

p_conf_forest <- ggplot(forest_data_conf, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(
    limits = c(-lim_conf, lim_conf),
    breaks = seq(-lim_conf, lim_conf, by = 2)
  ) +
  labs(
    title = "Fixed-Effects Estimates for Confidence Scores",
    x = "Estimate (change in confidence score)",
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

p_conf_forest

# -----------------------------
# Optional exports (tables)
# -----------------------------
# write.csv(conf_fixef_tbl_full,     "Confidence_FixedEffectsTable_FULL.csv",     row.names = FALSE)
# write.csv(conf_fixef_tbl_reduced,  "Confidence_FixedEffectsTable_REDUCED.csv",  row.names = FALSE)
# write.csv(conf_fixef_tbl_selected, "Confidence_FixedEffectsTable_SELECTED.csv", row.names = FALSE)

# -----------------------------
# 11) Export key plots to PNG (standardized)
# -----------------------------
fig_dir <- "figures_png_confidence"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

plots_to_save <- list(
  "confidence_box_by_phase_time"    = p_conf_box,
  "confidence_residuals_vs_fitted"  = p_conf_resid,
  "confidence_emm_by_phase"         = p_conf_emm,
  "confidence_forest_fixed_effects" = p_conf_forest
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