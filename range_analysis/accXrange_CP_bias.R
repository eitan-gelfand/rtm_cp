# Accuracy-by-range plot with Race as the plotted contrast.
# This figure shows how discrimination accuracy changes with morph distance,
# separately by Group and Bias condition, after averaging within subject.
#
# ──────────────────────────────────────────────────────────────
# Libraries
# ──────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(scales)

# ── Shared helpers ─────────────────────────────────────────────
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

# ──────────────────────────────────────────────────────────────
# 1) Load trial-level data and average within subject
# Range == 0 trials are excluded because they are not part of the main
# discrimination analysis.
# ──────────────────────────────────────────────────────────────
df_summary <- read.csv(data_path("full_data_cp.csv")) %>%
  filter(Range != 0) %>%
  group_by(Subject, ExperimentName, Range, Regression, Group) %>%
  summarize(mean_ACC = mean(ACC, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  ) %>%
  mutate(
    Regression = recode(Regression, "biasP" = "biasp", "biasM" = "biasm"),
    Regression = factor(Regression, levels = c("biasm", "biasp"))
  )

# ──────────────────────────────────────────────────────────────
# 2) Compute group-level means and standard errors
# ──────────────────────────────────────────────────────────────
stats_df_all <- df_summary %>%
  group_by(Group, ExperimentName, Range, Regression) %>%
  summarize(mean_val = mean(mean_ACC), se = sd(mean_ACC)/sqrt(n()), .groups = "drop")

# ──────────────────────────────────────────────────────────────
# 3) Plot Race differences across morph range
# ──────────────────────────────────────────────────────────────
p_all <- ggplot() +
  geom_line(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, group = ExperimentName, color = ExperimentName),
    linewidth = 1.2
  ) +
  geom_point(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, color = ExperimentName),
    size = 3
  ) +
  geom_errorbar(
    data = stats_df_all,
    aes(x = factor(Range), ymin = mean_val - se, ymax = mean_val + se, color = ExperimentName),
    width = 0.1, linewidth = 0.8
  ) +
  facet_grid(Regression ~ Group ,
             labeller = labeller(Regression = c("biasp" = "Bias+", "biasm" = "Bias−"))) +
  scale_x_discrete(
    breaks = as.character(unique(df_summary$Range)),
    labels = as.character(unique(df_summary$Range))
  ) +
  scale_colour_manual(
    values = c("Caucasian" = "#A6D8C3", "Asian" = "#F6B48F"),
    labels = c("Caucasian" = "Own-Race", "Asian" = "Other-Race"),
    name   = "Race"
  ) +
  labs(x = "Difference in morph levels", y = "Accuracy") +
  theme_pub()

# ──────────────────────────────────────────────────────────────
# 4) Save the figure
# ──────────────────────────────────────────────────────────────
ggsave(
  filename = plot_path("acc_range_CP_bias.png"),
  plot = p_all, width = 8, height = 6, dpi = 300,
  device = ragg::agg_png
)
