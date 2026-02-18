###############################################
# Linear Mixed-Effects Model (LMM) Analysis
# Diabetes Management Training – Knowledge Scores
#
# Goal (manuscript-ready workflow):
#   1) Describe score distributions over time (Pre, Post, Follow-up), overall and by Phase
#   2) Fit an LMM to test whether change over time differs by implementation Phase
#   3) Produce publishable visuals + tables:
#        - Descriptive plots
#        - Model diagnostics
#        - EMM (model-adjusted means) plot + EMM table
#        - Fixed-effects tables
#        - Forest plot of fixed effects (APA/clinical style)
#
# Conceptual model:
#   Score_ij = β0 + β(Timepoint) + β(Phase) + β(Timepoint×Phase) + u_i + ε_ij
#   where u_i is a participant-level random intercept (accounts for repeated measures).
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
# 2) Read data
# -----------------------------
lmm_data <- read_csv("LMM_long.csv")

# -----------------------------
# 3) Factor coding
# -----------------------------
time_levels  <- c("pretraining", "posttraining", "followup")
time_labels  <- c("Pre-training", "Post-training", "1-Month Follow-up")

phase_levels <- c("Phase 1", "Phase 2", "Phase 3")
phase_labels <- c("Phase 1", "Phase 2", "Phase 3")

lmm_data <- lmm_data %>%
  mutate(
    ID        = factor(ID),
    Timepoint = factor(Timepoint, levels = time_levels),
    Phase     = factor(Phase, levels = phase_levels),
    Score     = as.numeric(Score),
    Timepoint_label = factor(Timepoint, levels = time_levels, labels = time_labels),
    Phase_label     = factor(Phase, levels = phase_levels, labels = phase_labels)
  )

plot_data <- lmm_data %>% filter(!is.na(Score))

# -----------------------------
# 4) Missingness
# -----------------------------
table(lmm_data$Timepoint, is.na(lmm_data$Score))

# -----------------------------
# 5) Model fitting (ML for LRT)
# -----------------------------
lmm_full <- lmer(
  Score ~ Timepoint * Phase + (1 | ID),
  data = lmm_data,
  REML = FALSE,
  na.action = na.omit
)

lmm_reduced <- lmer(
  Score ~ Timepoint + Phase + (1 | ID),
  data = lmm_data,
  REML = FALSE,
  na.action = na.omit
)

interaction_test_knowledge <- anova(lmm_reduced, lmm_full)
p_int <- interaction_test_knowledge$`Pr(>Chisq)`[2]

lmm <- if (!is.na(p_int) && p_int < 0.05) lmm_full else lmm_reduced

summary(lmm)

# -----------------------------
# 6) Descriptive Box + Jitter (publication style)
# -----------------------------
p_know_box <- ggplot(plot_data, aes(x = Timepoint_label, y = Score)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 1) +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.4) +
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
# 7) Diagnostics (ggplot-based, consistent with Confidence)
# -----------------------------
diag_df_know <- data.frame(
  fitted = fitted(lmm),
  resid  = resid(lmm, type = "pearson")
)

p_know_resid <- ggplot(diag_df_know, aes(x = fitted, y = resid)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_point(size = 1.6, alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.8, se = FALSE, linewidth = 0.8) +
  labs(
    title = "Residuals vs Fitted (Knowledge Model)",
    x = "Fitted values",
    y = "Pearson residuals"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_know_resid

qqnorm(resid(lmm), main = "Normal Q–Q Plot of Residuals (Knowledge Model)")
qqline(resid(lmm))

# -----------------------------
# 8) EMM Plot (black-only, manuscript style)
# -----------------------------
emm <- emmeans(lmm, ~ Timepoint | Phase)

emm_df <- as.data.frame(emm) %>%
  mutate(
    Timepoint_label = factor(Timepoint, levels = time_levels, labels = time_labels),
    Phase_label     = factor(Phase, levels = phase_levels, labels = phase_labels)
  )

p_know_emm <- ggplot(
  emm_df,
  aes(
    x = Timepoint_label,
    y = emmean,
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
    legend.position = "bottom"
  )

p_know_emm

# -----------------------------
# 9) Fixed Effects Tables
# -----------------------------
fixef_tbl_selected <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
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

fixef_tbl_selected

# -----------------------------
# 10) Forest Plot (built from SELECTED model)
# -----------------------------
pp <- 100

forest_data_know <- broom.mixed::tidy(lmm, effects = "fixed", conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    estimate  = estimate  * pp,
    conf.low  = conf.low  * pp,
    conf.high = conf.high * pp,
    term_clean = case_when(
      term == "PhasePhase 2" ~ "Phase 2 vs Phase 1 (baseline)",
      term == "PhasePhase 3" ~ "Phase 3 vs Phase 1 (baseline)",
      term == "Timepointposttraining" ~ "Post vs Pre",
      term == "Timepointfollowup" ~ "Follow-up vs Pre",
      term == "Timepointposttraining:PhasePhase 2" ~ "Post × Phase 2",
      term == "Timepointposttraining:PhasePhase 3" ~ "Post × Phase 3",
      term == "Timepointfollowup:PhasePhase 2" ~ "Follow-up × Phase 2",
      term == "Timepointfollowup:PhasePhase 3" ~ "Follow-up × Phase 3",
      TRUE ~ term
    ),
    block = ifelse(grepl("×", term_clean), "Interaction terms", "Main effects"),
    block = factor(block, levels = c("Main effects", "Interaction terms"))
  )

preferred_order_know <- c(
  "Phase 2 vs Phase 1 (baseline)",
  "Phase 3 vs Phase 1 (baseline)",
  "Post vs Pre",
  "Follow-up vs Pre",
  "Post × Phase 2",
  "Post × Phase 3",
  "Follow-up × Phase 2",
  "Follow-up × Phase 3"
)

present_terms <- preferred_order_know[preferred_order_know %in% forest_data_know$term_clean]

forest_data_know <- forest_data_know %>%
  mutate(term_clean = factor(term_clean, levels = rev(present_terms)))

lim <- max(abs(c(forest_data_know$conf.low, forest_data_know$conf.high)), na.rm = TRUE)
lim <- ceiling(lim / 10) * 10

p_know_forest <- ggplot(forest_data_know, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18, linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(-lim, lim), breaks = seq(-lim, lim, by = 10)) +
  labs(
    title = "Fixed-Effects Estimates for Knowledge Scores",
    x = "Estimate (percentage-point change in knowledge score)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 9)
  )

p_know_forest


# -----------------------------
# 11) Export plots to PNG
# -----------------------------
fig_dir <- "figures_png_knowledge"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

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

message("Saved ", length(plots_to_save), " knowledge plots to: ", normalizePath(fig_dir))