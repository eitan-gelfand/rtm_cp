# Lambda histogram for Bias- and Bias+ trials in young and full samples.
# Input metric: lambda from weibull_metric_results.csv.
# Output plot: output/plots/bias_split_lambda_young_and_all.png

library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(scales)

.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

df <- read.csv(data_path("weibull_metric_results.csv"))

age_lookup <- read.csv(data_path("dprime_results_with_range.csv")) %>%
  group_by(Subject) %>%
  summarize(Age = first(Age), .groups = "drop")

df_bias <- df %>%
  left_join(age_lookup, by = "Subject") %>%
  mutate(
    Regression = recode(Regression, "biasM" = "biasm", "biasP" = "biasp"),
    Regression = factor(Regression, levels = c("biasm", "biasp")),
    BiasLabel = recode(
      as.character(Regression),
      "biasm" = "Bias- trials",
      "biasp" = "Bias+ trials"
    ),
    BiasLabel = factor(BiasLabel, levels = c("Bias- trials", "Bias+ trials"))
  )

df_plot_source <- bind_rows(
  df_bias %>%
    filter(Age <= 42) %>%
    mutate(Sample = "Young (Age <= 42)"),
  df_bias %>%
    mutate(Sample = "All ages")
) %>%
  mutate(
    Sample = factor(Sample, levels = c("Young (Age <= 42)", "All ages"))
  )

# Average to obtain one value per subject x bias x race x sample.
df_grouped <- df_plot_source %>%
  group_by(Sample, BiasLabel, Regression, Subject, ExperimentName) %>%
  summarize(
    lambda = mean(lambda, na.rm = TRUE),
    Age = first(Age),
    Group = first(Group),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  )

# Race is a within-subject factor and Group is between-subject, so error bars
# use within-subject normalization separately inside each Sample x Bias x Group cell.
k_within <- nlevels(df_grouped$ExperimentName)
morey_correction <- sqrt(k_within / (k_within - 1))

df_grouped_ws <- df_grouped %>%
  group_by(Sample, BiasLabel, Group, Subject) %>%
  mutate(subject_mean = mean(lambda, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample, BiasLabel, Group) %>%
  mutate(group_grand_mean = mean(lambda, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(lambda_norm = lambda - subject_mean + group_grand_mean)

df_summary <- df_grouped_ws %>%
  group_by(Sample, BiasLabel, Group, ExperimentName) %>%
  summarize(
    meanlambda = mean(lambda, na.rm = TRUE),
    sd_norm = sd(lambda_norm, na.rm = TRUE),
    n = n_distinct(Subject),
    selambda = (sd_norm / sqrt(n)) * morey_correction,
    .groups = "drop"
  )

dodge_width <- 0.8
bar_width <- 0.65

p <- ggplot(df_summary, aes(x = Group, y = meanlambda)) +
  geom_col(
    aes(fill = ExperimentName),
    position = position_dodge(width = dodge_width),
    width = bar_width,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = meanlambda - selambda,
      ymax = meanlambda + selambda,
      group = ExperimentName
    ),
    width = 0.2,
    position = position_dodge(dodge_width)
  ) +
  facet_grid(BiasLabel ~ Sample) +
  labs(x = "Group", y = "lambda") +
  scale_fill_manual(
    name = "Race",
    values = c("Caucasian" = "#A6D8C3", "Asian" = "#F6B48F"),
    labels = c("Caucasian" = "Own-Race", "Asian" = "Other-Race")
  ) +
  scale_y_continuous(
    labels = number_format(accuracy = 0.01),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08))
  ) +
  theme_pub()

# Paired race comparisons within each sample x bias x group.
stat_test_race_all <- df_grouped %>%
  group_by(Sample, BiasLabel, Group) %>%
  t_test(lambda ~ ExperimentName, paired = TRUE) %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif)) %>%
  add_xy_position(x = "Group", dodge = dodge_width)

race_base_y <- df_summary %>%
  group_by(Sample, BiasLabel, Group) %>%
  summarize(y.position = max(meanlambda + selambda, na.rm = TRUE) + 0.08,
            .groups = "drop")

stat_test_race_all <- stat_test_race_all %>%
  select(-y.position) %>%
  left_join(race_base_y, by = c("Sample", "BiasLabel", "Group"))

cat("\nPaired race comparisons within each sample x bias x group:\n")
print(
  stat_test_race_all %>%
    select(Sample, BiasLabel, Group, group1, group2, n1, n2, statistic, df, p, p.signif),
  n = Inf,
  width = Inf
)

stat_test_race <- stat_test_race_all %>%
  filter(p.signif != "ns")

# Independent TD vs CP comparisons within each sample x bias x race.
stat_test_group <- df_grouped %>%
  group_by(Sample, BiasLabel, ExperimentName) %>%
  t_test(lambda ~ Group, paired = FALSE) %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif))

facet_y <- df_summary %>%
  group_by(Sample, BiasLabel) %>%
  summarize(
    facet_top = max(meanlambda + selambda, na.rm = TRUE),
    facet_range = diff(range(meanlambda + selambda, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(facet_range = if_else(is.finite(facet_range) & facet_range > 0, facet_range, 1))

group_base_y <- df_summary %>%
  group_by(Sample, BiasLabel, ExperimentName) %>%
  summarize(group_base = max(meanlambda + selambda, na.rm = TRUE), .groups = "drop")

stat_test_group <- stat_test_group %>%
  left_join(facet_y, by = c("Sample", "BiasLabel")) %>%
  left_join(group_base_y, by = c("Sample", "BiasLabel", "ExperimentName")) %>%
  mutate(
    race_index = as.numeric(ExperimentName),
    x_shift = (race_index - (nlevels(ExperimentName) + 1) / 2) *
      (dodge_width / nlevels(ExperimentName)),
    xmin = as.numeric(factor(group1, levels = levels(df_grouped$Group))) + x_shift,
    xmax = as.numeric(factor(group2, levels = levels(df_grouped$Group))) + x_shift,
    y.position = pmax(group_base, facet_top) + 0.16 * facet_range +
      (race_index - 1) * 0.08 * facet_range
  ) %>%
  select(-facet_top, -facet_range, -group_base, -race_index, -x_shift)

cat("\nIndependent TD vs CP comparisons within each sample x bias x race:\n")
print(
  stat_test_group %>%
    select(Sample, BiasLabel, ExperimentName, group1, group2, n1, n2, statistic, df, p, p.signif),
  n = Inf,
  width = Inf
)

stat_test_group <- stat_test_group %>%
  filter(p.signif != "ns")

final_plot <- p +
  stat_pvalue_manual(
    stat_test_race,
    label = "p.signif",
    bracket.size = 0.6,
    tip.length = 0.01,
    size = 5
  ) +
  stat_pvalue_manual(
    stat_test_group,
    label = "p.signif",
    bracket.size = 0.6,
    tip.length = 0.01,
    size = 5
  )

if (interactive()) {
  print(final_plot)
}

ggsave(
  plot_path("bias_split_lambda_young_and_all.png"),
  final_plot,
  width = 9,
  height = 8.5,
  dpi = 300,
  device = ragg::agg_png
)
