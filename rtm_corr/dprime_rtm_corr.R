# =========================================================
# Correlation analysis from raw CP data using d-prime
# - Start from trial-level data
# - Merge per-condition d-prime values
# - Clean practice trials by excluding Range == 0
# - Aggregate ACC by Subject x Group x ExperimentName x Regression
# - Compute Magnitude = biasp - biasm per subject x race
# - Use biasM d-prime only per subject x race
# - Run Pearson correlations between biasM d-prime and RTM magnitude
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
# Load data
# ---------------------------------------------------------
df <- read.csv(data_path("full_data_cp.csv"))
dprime_df <- read.csv(data_path("dprime_results.csv"))

# ---------------------------------------------------------
# 1) Merge per-condition d-prime onto the trial-level data
# Each subject has a different d-prime in each race x regression condition.
# ---------------------------------------------------------
merge_keys <- c("Subject", "Regression", "ExperimentName")

df_with_dprime <- df %>%
  left_join(
    dprime_df %>%
      select(all_of(merge_keys), Dprime),
    by = merge_keys
  )

# ---------------------------------------------------------
# 2) Clean data
# Practice / null trials are Range == 0, so exclude them
# ---------------------------------------------------------
df_clean <- df_with_dprime %>%
  filter(Range != 0)

# ---------------------------------------------------------
# 3) Aggregate to subject-level condition means
# This gives 4 rows per subject:
#   2 ExperimentName levels x 2 Regression levels
# ---------------------------------------------------------
df_agg <- df_clean %>%
  group_by(Subject, Group, ExperimentName, Regression) %>%
  summarise(
    mean_acc = mean(ACC, na.rm = TRUE),
    Dprime = mean(Dprime, na.rm = TRUE),
    n_trials = n(),
    .groups = "drop"
  )

# ---------------------------------------------------------
# 4) Pivot wide to one row per subject x race
# RTM magnitude is computed per subject x race.
# d-prime is taken from biasM only.
# ---------------------------------------------------------
df_corr <- df_agg %>%
  select(Subject, Group, ExperimentName, Regression, mean_acc, Dprime, n_trials) %>%
  pivot_wider(
    names_from  = Regression,
    values_from = c(mean_acc, Dprime, n_trials),
    names_sep   = "_"
  ) %>%
  mutate(
    biasM     = mean_acc_biasm,
    biasP     = mean_acc_biasp,
    Magnitude = biasP - biasM,
    Dprime    = Dprime_biasm
  ) %>%
  select(
    Subject, Group, ExperimentName,
    Dprime, Magnitude,
    biasM, biasP,
    Dprime_biasm, Dprime_biasp,
    n_trials_biasm, n_trials_biasp
  ) %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  )

# ---------------------------------------------------------
# 5) Correlation statistics
# Pearson correlation within each Group x ExperimentName
# using biasM d-prime only
# ---------------------------------------------------------
corr_stats <- df_corr %>%
  group_by(Group, ExperimentName) %>%
  summarise(
    n = sum(complete.cases(Dprime, Magnitude)),
    r = cor(Dprime, Magnitude, method = "pearson", use = "complete.obs"),
    p = cor.test(Dprime, Magnitude, method = "pearson")$p.value,
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
# 6) Annotation positions — scaled to d-prime range
# ---------------------------------------------------------
x_range <- range(df_corr$Dprime, na.rm = TRUE)
y_range <- range(df_corr$Magnitude, na.rm = TRUE)
x_span <- diff(x_range)
y_span <- diff(y_range)

if (!is.finite(x_span) || x_span == 0) {
  x_span <- 1
}

if (!is.finite(y_span) || y_span == 0) {
  y_span <- 1
}

corr_stats_plot <- corr_stats %>%
  mutate(
    x_pos = x_range[1] + 0.04 * x_span,
    y_pos = y_range[2] - (0.06 + (annotation_rank - 1) * 0.12) * y_span
  )

# ---------------------------------------------------------
# 7) Plot
# ---------------------------------------------------------
p <- ggplot(
  df_corr,
  aes(x = Dprime, y = Magnitude)
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
    x = "d-prime (Bias- trials)",
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
  scale_x_continuous(
    expand = expansion(mult = c(0.06, 0.14))
  ) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

if (interactive()) {
  print(p)
}

# ---------------------------------------------------------
# 8) Save plot
# ---------------------------------------------------------
ggsave(
  filename = plot_path("Correlation_dprime_biasM_vs_magnitude_CP.png"),
  plot     = p,
  width    = 10,
  height   = 6,
  dpi      = 300,
  device   = ragg::agg_png
)