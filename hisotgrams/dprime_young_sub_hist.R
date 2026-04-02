# Subject-level d-prime summary plot for the younger subgroup only.
# This figure mirrors the main d-prime histogram but restricts the sample to
# the younger-age subset used elsewhere in the project (Age <= 42).
#
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
# 1) Load d-prime results with Range-level values
# ──────────────────────────────────────────────────────────────
df <- read.csv(data_path("dprime_results_with_range.csv"))
# Restrict to the younger subgroup used in the manuscript's age-based follow-up.
df <- df %>% filter(Age <= 42)
# ──────────────────────────────────────────────────────────────
# 2) Aggregate over Range to obtain one value per subject x condition
# ──────────────────────────────────────────────────────────────
df_grouped <- df %>%
  group_by(Subject, ExperimentName, Regression) %>%
  summarize(
    Dprime = mean(Dprime, na.rm = TRUE),
    CR     = mean(CR, na.rm = TRUE),
    Group  = first(Group),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, 
                            levels = c("Caucasian", "Asian")),
    Regression = recode(Regression, "biasP" = "biasp", "biasM" = "biasm"),
    Regression = factor(Regression, levels = c("biasm", "biasp"))
  )

# ──────────────────────────────────────────────────────────────
# 3) Summary for plotting
# Error bars use a within-subject Morey correction because Race and Regression
# vary within subject whereas Group varies between subjects.
# ──────────────────────────────────────────────────────────────

# number of within-subject cells per participant:
# ExperimentName (2) × Regression (2) = 4
k_within <- nlevels(df_grouped$ExperimentName) * nlevels(df_grouped$Regression)
morey_correction <- sqrt(k_within / (k_within - 1))

# normalize within each subject, separately inside each between-subject Group
df_grouped_ws <- df_grouped %>%
  group_by(Group, Subject) %>%
  mutate(subject_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Group) %>%
  mutate(group_grand_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Dprime_norm = Dprime - subject_mean + group_grand_mean
  )

df_summary <- df_grouped_ws %>%
  group_by(ExperimentName, Regression, Group) %>%
  summarize(
    meanDprime = mean(Dprime, na.rm = TRUE),                 # raw mean for bar height
    sd_norm    = sd(Dprime_norm, na.rm = TRUE),             # SD after within-subject normalization
    n          = n_distinct(Subject),
    seDprime   = (sd_norm / sqrt(n)) * morey_correction,    # Morey-corrected within-subject SE
    .groups    = "drop"
  )

y_range <- diff(range(df_summary$meanDprime, na.rm = TRUE))
if (!is.finite(y_range) || y_range <= 0) y_range <- 1

within_offset_step <- 0.03 * y_range
between_offset_base <- 0.10 * y_range
between_offset_step <- 0.03 * y_range
dodge_width <- 0.8
n_reg <- nlevels(df_grouped$Regression)
group_levels <- levels(df_grouped$Group)
x_td <- which(group_levels == "TD")
x_cp <- which(group_levels == "CP")


# ──────────────────────────────────────────────────────────────
# 4) Base histogram
# ──────────────────────────────────────────────────────────────
p_cp <- ggplot(df_summary, aes(x = Group, y = meanDprime)) +
  geom_col(
    aes(fill = Regression),
    position = position_dodge(width = 0.8),
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin  = meanDprime - seDprime,
      ymax  = meanDprime + seDprime,
      group = Regression
    ),
    width    = 0.2,
    position = position_dodge(0.8)
  ) +
  facet_wrap(
    ~ ExperimentName,
    labeller = labeller(
      ExperimentName = c("Asian" = "Other-Race",
                         "Caucasian" = "Own-Race")
    )
  ) +
  labs(x = "Group", y = "d-prime") +
  scale_fill_manual(
    name   = "Regression",
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias-")
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  theme_pub()

# ──────────────────────────────────────────────────────────────
# 5) Bias+ vs Bias- comparison within each Group x Race panel
# ──────────────────────────────────────────────────────────────
stat_test_reg <- df_grouped %>%
  group_by(ExperimentName, Group) %>%
  t_test(Dprime ~ Regression, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif)) %>%
  filter(p.signif != "ns") %>%
  add_xy_position(x = "Group", dodge = 0.8)

reg_base_y <- df_summary %>%
  group_by(ExperimentName, Group) %>%
  summarize(
    reg_base = max(meanDprime + seDprime, na.rm = TRUE),
    reg_se_top = max(seDprime, na.rm = TRUE),
    .groups = "drop"
  )

stat_test_reg <- stat_test_reg %>%
  left_join(reg_base_y, by = c("ExperimentName", "Group")) %>%
  mutate(
    within_gap = pmax(0.02 * y_range, 0.35 * reg_se_top),
    y.position = reg_base + within_gap +
      (as.numeric(Group) - 1) * within_offset_step
  )

# ──────────────────────────────────────────────────────────────
# 6) TD vs CP comparison within each Race x Bias condition
# ──────────────────────────────────────────────────────────────
stat_test_group <- df_grouped %>%
  group_by(ExperimentName, Regression) %>%
  t_test(Dprime ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif)) %>%
  filter(p.signif != "ns") %>%
  add_xy_position(x = "Group", dodge = 0.8)

group_base_y <- df_summary %>%
  group_by(ExperimentName, Regression) %>%
  summarize(group_base = max(meanDprime + seDprime, na.rm = TRUE), .groups = "drop")

stat_test_group <- stat_test_group %>%
  left_join(group_base_y, by = c("ExperimentName", "Regression")) %>%
  mutate(y.position = group_base + between_offset_base +
           (as.numeric(Regression) - 1) * between_offset_step) %>%
  mutate(
    reg_index = as.numeric(Regression),
    x_shift = (reg_index - (n_reg + 1) / 2) * (dodge_width / n_reg),
    xmin = x_td + x_shift,
    xmax = x_cp + x_shift
  ) %>%
  select(-reg_index, -x_shift)

# ──────────────────────────────────────────────────────────────
# 7) Add significance brackets to the final plot
# ──────────────────────────────────────────────────────────────
final_plot_cp <- p_cp +
  stat_pvalue_manual(
    stat_test_reg,
    label        = "p.signif",
    bracket.size = 0.8,
    tip.length   = 0.01,
    size         = 6
  ) +
  stat_pvalue_manual(
    stat_test_group,
    label        = "p.signif",
    bracket.size = 0.8,
    tip.length   = 0.01,
    size         = 6
  )

if (interactive()) {
  print(final_plot_cp)
}

# ──────────────────────────────────────────────────────────────
# 8) Save
# ──────────────────────────────────────────────────────────────
ggsave(
  plot_path("dprime_histogram_young.png"),
  final_plot_cp,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)
