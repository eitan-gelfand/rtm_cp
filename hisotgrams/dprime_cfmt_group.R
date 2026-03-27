# ──────────────────────────────────────────────────────────────
# Libraries
# ──────────────────────────────────────────────────────────────
library(dplyr)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(tibble)
library(scales)

# ── Shared helpers ─────────────────────────────────────────────
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

# ──────────────────────────────────────────────────────────────
# 1) Load CP Dprime dataset with CFMT grouping
# ──────────────────────────────────────────────────────────────
df <- read.csv(data_path("dprime_results_with_range_and_group.csv"))

# ──────────────────────────────────────────────────────────────
# 2) Aggregate over Range (subject-level histogram)
# ──────────────────────────────────────────────────────────────
cfmt_levels <- sort(unique(stats::na.omit(as.character(df$CFMTgroup))))
preferred_cfmt_levels <- c("L", "H")
cfmt_levels <- c(preferred_cfmt_levels[preferred_cfmt_levels %in% cfmt_levels],
                 setdiff(cfmt_levels, preferred_cfmt_levels))

df_grouped <- df %>%
  filter(!is.na(CFMTgroup), CFMTgroup != "") %>%
  group_by(Subject, ExperimentName, Regression) %>%
  summarize(
    Dprime    = mean(Dprime, na.rm = TRUE),
    CR        = mean(CR, na.rm = TRUE),
    CFMTgroup = first(CFMTgroup),
    .groups   = "drop"
  ) %>%
  mutate(
    CFMTgroup = factor(CFMTgroup, levels = cfmt_levels),
    ExperimentName = factor(
      ExperimentName,
      levels = c("Caucasian", "Asian")
    ),
    Regression = recode(Regression, "biasP" = "biasp", "biasM" = "biasm"),
    Regression = factor(Regression, levels = c("biasm", "biasp"))
  )

# ──────────────────────────────────────────────────────────────
# 3) Summary for plotting
# ──────────────────────────────────────────────────────────────
# number of within-subject cells per participant:
# ExperimentName (2) × Regression (2) = 4
k_within <- nlevels(df_grouped$ExperimentName) * nlevels(df_grouped$Regression)
morey_correction <- sqrt(k_within / (k_within - 1))

# normalize within each subject, separately inside each between-subject CFMT group
df_grouped_ws <- df_grouped %>%
  group_by(CFMTgroup, Subject) %>%
  mutate(subject_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(CFMTgroup) %>%
  mutate(group_grand_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Dprime_norm = Dprime - subject_mean + group_grand_mean
  )

df_summary <- df_grouped_ws %>%
  group_by(ExperimentName, Regression, CFMTgroup) %>%
  summarize(
    meanDprime = mean(Dprime, na.rm = TRUE),
    sd_norm    = sd(Dprime_norm, na.rm = TRUE),
    n          = n_distinct(Subject),
    seDprime   = (sd_norm / sqrt(n)) * morey_correction,
    .groups    = "drop"
  )

# ──────────────────────────────────────────────────────────────
# 4) Base histogram
# ──────────────────────────────────────────────────────────────
p_cp_cfmt <- ggplot(df_summary, aes(x = CFMTgroup, y = meanDprime)) +
  geom_col(
    aes(fill = Regression),
    position = position_dodge(width = 0.8),
    width = 0.8,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = meanDprime - seDprime,
      ymax = meanDprime + seDprime,
      group = Regression
    ),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(
    ~ ExperimentName,
    labeller = labeller(
      ExperimentName = c(
        "Asian" = "Other-Race",
        "Caucasian" = "Own-Race"
      )
    )
  ) +
  labs(x = "CFMT Group", y = "d-prime") +
  scale_fill_manual(
    name = "Regression",
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias-"),
    drop = FALSE
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  theme_pub()

# ──────────────────────────────────────────────────────────────
# 5) CFMT-group comparison within each Race panel, separately by Regression
# ──────────────────────────────────────────────────────────────
stat_test_cfmt <- df_grouped %>%
  group_by(ExperimentName, Regression) %>%
  t_test(Dprime ~ CFMTgroup, paired = FALSE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif)) %>%
  filter(p.signif != "ns") %>%
  add_xy_position(x = "CFMTgroup", dodge = 0.8)

y_range <- diff(range(df_summary$meanDprime, na.rm = TRUE))
if (!is.finite(y_range) || y_range <= 0) y_range <- 1

cfmt_base_y <- df_summary %>%
  group_by(ExperimentName, Regression) %>%
  summarize(
    bracket_base = max(meanDprime + seDprime, na.rm = TRUE),
    .groups = "drop"
  )

stat_test_cfmt <- stat_test_cfmt %>%
  left_join(cfmt_base_y, by = c("ExperimentName", "Regression")) %>%
  mutate(
    y.position = bracket_base + 0.08 * y_range +
      (as.numeric(Regression) - 1) * 0.03 * y_range
  )

# ──────────────────────────────────────────────────────────────
# 6) Add significance brackets
# ──────────────────────────────────────────────────────────────
final_plot_cp_cfmt <- p_cp_cfmt +
  stat_pvalue_manual(
    stat_test_cfmt,
    label = "p.signif",
    bracket.size = 0.8,
    tip.length = 0.01,
    size = 6
  )

if (interactive()) {
  print(final_plot_cp_cfmt)
}

# ──────────────────────────────────────────────────────────────
# 7) Save
# ──────────────────────────────────────────────────────────────
ggsave(
  plot_path("dprime_cfmt_group.png"),
  final_plot_cp_cfmt,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)
