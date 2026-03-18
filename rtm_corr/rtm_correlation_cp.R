# =========================================================
# Correlation analysis from raw CP data
# - Start from trial-level data
# - Clean practice trials by excluding Range == 0
# - Aggregate ACC by Subject x Group x ExperimentName x Regression
# - Pivot wide
# - Compute Magnitude = biasp - biasm
# - Run Pearson correlations
# - Plot CP and TD side-by-side with race-specific annotations
# =========================================================

library(tidyverse)

# ---------------------------------------------------------
# Helpers
# Assumes these already exist in your project
# ---------------------------------------------------------
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

base_family <- getOption("project.base_family", "serif")


# ---------------------------------------------------------
# Load raw data
# ---------------------------------------------------------
df <- read.csv(data_path("full_data_cp.csv"))

# ---------------------------------------------------------
# 1) Clean data
# Practice / null trials are Range == 0, so exclude them
# ---------------------------------------------------------
df_clean <- df %>%
  filter(Range != 0)

# ---------------------------------------------------------
# 2) Aggregate to subject-level condition means
# This gives 4 rows per subject:
#   2 ExperimentName levels x 2 Regression levels
# ---------------------------------------------------------
df_agg <- df_clean %>%
  group_by(Subject, Group, ExperimentName, Regression) %>%
  summarise(
    mean_acc = mean(ACC, na.rm = TRUE),
    n_trials = n(),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 3) Pivot wide
# Final correlation table: 2 rows per subject
#   one for Asian, one for Caucasian
# ---------------------------------------------------------
df_corr <- df_agg %>%
  select(Subject, Group, ExperimentName, Regression, mean_acc, n_trials) %>%
  pivot_wider(
    names_from  = Regression,
    values_from = c(mean_acc, n_trials),
    names_sep   = "_"
  ) %>%
  mutate(
    biasM     = mean_acc_biasm,
    biasP     = mean_acc_biasp,
    Magnitude = biasP - biasM
  ) %>%
  select(
    Subject, Group, ExperimentName,
    biasM, biasP, Magnitude,
    n_trials_biasm, n_trials_biasp
  ) %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  )

# ---------------------------------------------------------
# 4) Correlation statistics
# Pearson correlation within each Group x ExperimentName
# ---------------------------------------------------------
corr_stats <- df_corr %>%
  group_by(Group, ExperimentName) %>%
  summarise(
    n = sum(complete.cases(biasM, Magnitude)),
    r = cor(biasM, Magnitude, method = "pearson", use = "complete.obs"),
    p = cor.test(biasM, Magnitude, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    race_label = recode(
      ExperimentName,
      "Caucasian" = "Own-Race",
      "Asian" = "Other-Race"
    ),
    sig_stars = case_when(
      p < .001 ~ "***",
      p < .01  ~ "**",
      p < .05  ~ "*",
      TRUE     ~ ""
    ),
    p_fmt = if_else(p < .001, "< .001", paste0("= ", sprintf("%.3f", p))),
    label = paste0(
      race_label,
      ": r = ", sprintf("%.3f", r),
      ", p ", p_fmt, sig_stars
    ),
    annotation_rank = match(as.character(ExperimentName), c("Caucasian", "Asian"))
  )

print(corr_stats)

# ---------------------------------------------------------
# 5) Annotation positions — fixed top-left per panel
# ---------------------------------------------------------
corr_stats_plot <- corr_stats %>%
  mutate(
    x_pos = 0.2,
    y_pos = 0.4 - (annotation_rank - 1) * 0.05
  )

# ---------------------------------------------------------
# 6) Plot
# ---------------------------------------------------------
p <- ggplot(
  df_corr,
  aes(x = biasM, y = Magnitude)
) +
  facet_grid(cols = vars(Group)) +
  geom_point(
    aes(fill = ExperimentName),
    shape = 21,
    color = "white",
    stroke = 0.25,
    size = 2.2,
    alpha = 0.62
  ) +
  geom_smooth(
    aes(color = ExperimentName, fill = ExperimentName),
    method = "lm",
    se = TRUE,
    alpha = 0.20,
    linewidth = 1.6
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5,
    alpha = 0.5
  ) +
  geom_text(
    data = corr_stats_plot,
    aes(x = x_pos, y = y_pos, label = label),
    color = "#4F4F4F",
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 4.2,
    family = base_family,
    fontface = "bold"
  ) +
  labs(
    x = "Bias- accuracy",
    y = "RTM magnitude (Bias+ - Bias-)"
  ) +
  theme_pub(base_family = base_family) +
  scale_colour_manual(
    values = c(
      "Caucasian" = "#A6D8C3",
      "Asian"     = "#F6B48F"
    ),
    labels = c(
      "Caucasian" = "Own-Race",
      "Asian"     = "Other-Race"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Caucasian" = "#A6D8C3",
      "Asian"     = "#F6B48F"
    ),
    labels = c(
      "Caucasian" = "Own-Race",
      "Asian"     = "Other-Race"
    ),
    guide = "none"
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

if (interactive()) {
  print(p)
}

# ---------------------------------------------------------
# 7) Save plot
# ---------------------------------------------------------
ggsave(
  filename = plot_path("Correlation_biasM_vs_magnitude_CP.png"),
  plot     = p,
  width    = 10,
  height   = 6,
  dpi      = 300,
  device   = ragg::agg_png
)

