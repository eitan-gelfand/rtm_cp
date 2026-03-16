# ──────────────────────────────────────────────────────────────
# Libraries
# ──────────────────────────────────────────────────────────────
library(ggplot2)
library(dplyr)
library(showtext)
library(scales)

# ── Shared helpers ─────────────────────────────────────────────
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

# ──────────────────────────────────────────────────────────────
# 1) Load & summarize
# ──────────────────────────────────────────────────────────────
df_summary <- read.csv(data_path("full_data_cp.csv")) %>%
  filter(Range != 0) %>%
  group_by(Subject, ExperimentName, Range, Regression, Group) %>%
  summarize(mean_ACC = mean(ACC, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),  # adapt to your CP data
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  ) %>%
  mutate(
    Regression = recode(Regression, "biasP" = "biasp", "biasM" = "biasm"),
    Regression = factor(Regression, levels = c("biasp", "biasm"))
  )

# ──────────────────────────────────────────────────────────────
# 2) Means & SEs
# ──────────────────────────────────────────────────────────────
stats_df_all <- df_summary %>%
  group_by(Group, ExperimentName, Range, Regression) %>%
  summarize(
    mean_val = mean(mean_ACC),
    se       = sd(mean_ACC) / sqrt(n()),
    .groups  = "drop"
  )

# ──────────────────────────────────────────────────────────────
# 3) Plot
# ──────────────────────────────────────────────────────────────
p_all <- ggplot() +
  # geom_jitter(
  #   data = df_summary,
  #   aes(x = factor(Range), y = mean_ACC, color = Regression),
  #   size = 1.5, alpha = 0.5,
  #   position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0)
  # ) +
  geom_line(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, group = Regression, color = Regression),
    size = 1.2
  ) +
  geom_point(
    data = stats_df_all,
    aes(x = factor(Range), y = mean_val, color = Regression),
    size = 3
  ) +
  geom_errorbar(
    data = stats_df_all,
    aes(x = factor(Range), ymin = mean_val - se, ymax = mean_val + se, color = Regression),
    width = 0.1, size = 0.8
  ) +
  facet_grid( ExperimentName ~ Group,
             labeller = labeller(ExperimentName = c("Asian" = "Other-Race", "Caucasian" = "Own-Race"))) +
  scale_x_discrete(
    breaks = as.character(unique(df_summary$Range)),
    labels = as.character(unique(df_summary$Range))
  ) +
  scale_colour_manual(
    values = c(biasp = "#7FB3D5", biasm = "#F08080"),
    labels = c(biasp = "Bias+", biasm = "Bias−")
  ) +
  labs(x = "Difference in morph levels", y = "Accuracy") +
  theme_pub()

# ──────────────────────────────────────────────────────────────
# 4) Save
# ──────────────────────────────────────────────────────────────
ggsave(
  filename = plot_path("acc_range_CP_race.png"),
  plot = p_all, width = 8, height = 6, dpi = 300
)
