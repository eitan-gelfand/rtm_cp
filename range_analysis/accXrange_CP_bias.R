# ──────────────────────────────────────────────────────────────
# Libraries
# ──────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(showtext)
library(scales)

font_add("Times New Roman", regular = "times.ttf")
showtext_opts(dpi = 300)
showtext_auto()

# ──────────────────────────────────────────────────────────────
# 1) Load & summarize
# ──────────────────────────────────────────────────────────────
df_summary <- read.csv("C:/Users/eitas/OneDrive/University/Semester D/Galia's Lab/R/data/full_data_cp_experiment.csv") %>%
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
# 2) Means & SEs
# ──────────────────────────────────────────────────────────────
stats_df_all <- df_summary %>%
  group_by(Group, ExperimentName, Range, Regression) %>%
  summarize(mean_val = mean(mean_ACC), se = sd(mean_ACC)/sqrt(n()), .groups = "drop")

# ──────────────────────────────────────────────────────────────
# 3) Plot
# ──────────────────────────────────────────────────────────────
p_all <- ggplot() +
  geom_line(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, group = ExperimentName, color = ExperimentName),
    size = 1.2
  ) +
  geom_point(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, color = ExperimentName),
    size = 3
  ) +
  geom_errorbar(
    data = stats_df_all,
    aes(x = factor(Range), ymin = mean_val - se, ymax = mean_val + se, color = ExperimentName),
    width = 0.1, size = 0.8
  ) +
  facet_grid(Group ~ Regression,
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
# 4) Save
# ──────────────────────────────────────────────────────────────
ggsave(
  filename = "C:/Users/eitas/OneDrive/University/Semester D/Galia's Lab/R/CP_TD_experiment/CP_experiment_plots/acc_range_CP_bias.png",
  plot = p_all, width = 8, height = 6, dpi = 300
)
