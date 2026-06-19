# D-prime histogram for Bias- trials in young and full samples.
# Input metric: Dprime from dprime_results_with_range.csv.
# Output plot: output/plots/biasm_dprime_young_and_all.png

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

df <- read.csv(data_path("dprime_results_with_range.csv"))

df_biasm <- df %>%
  mutate(
    Regression = recode(Regression, "biasM" = "biasm", "biasP" = "biasp")
  ) %>%
  filter(Regression == "biasm")

df_plot_source <- bind_rows(
  df_biasm %>%
    filter(Age <= 42) %>%
    mutate(Sample = "Young (Age <= 42)"),
  df_biasm %>%
    mutate(Sample = "All ages")
) %>%
  mutate(
    Sample = factor(Sample, levels = c("Young (Age <= 42)", "All ages"))
  )

# Average over Range to obtain one value per subject x race x sample.
df_grouped <- df_plot_source %>%
  group_by(Sample, Subject, ExperimentName) %>%
  summarize(
    Dprime = mean(Dprime, na.rm = TRUE),
    Age = first(Age),
    Group = first(Group),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian", "Asian"))
  )

# Race is a within-subject factor and Group is between-subject, so error bars
# use within-subject normalization separately inside each Sample x Group cell.
k_within <- nlevels(df_grouped$ExperimentName)
morey_correction <- sqrt(k_within / (k_within - 1))

df_grouped_ws <- df_grouped %>%
  group_by(Sample, Group, Subject) %>%
  mutate(subject_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Sample, Group) %>%
  mutate(group_grand_mean = mean(Dprime, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Dprime_norm = Dprime - subject_mean + group_grand_mean)

df_summary <- df_grouped_ws %>%
  group_by(Sample, Group, ExperimentName) %>%
  summarize(
    meanDprime = mean(Dprime, na.rm = TRUE),
    sd_norm = sd(Dprime_norm, na.rm = TRUE),
    n = n_distinct(Subject),
    seDprime = (sd_norm / sqrt(n)) * morey_correction,
    .groups = "drop"
  )

dodge_width <- 0.8
bar_width <- 0.65

p <- ggplot(df_summary, aes(x = Group, y = meanDprime)) +
  geom_col(
    aes(fill = ExperimentName),
    position = position_dodge(width = dodge_width),
    width = bar_width,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = meanDprime - seDprime,
      ymax = meanDprime + seDprime,
      group = ExperimentName
    ),
    width = 0.2,
    position = position_dodge(dodge_width)
  ) +
  facet_wrap(~ Sample, nrow = 1) +
  labs(x = "Group", y = "D-prime (Bias- trials)") +
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

# Paired race comparisons within each sample x group.
stat_test_race_all <- df_grouped %>%
  group_by(Sample, Group) %>%
  t_test(Dprime ~ ExperimentName, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  mutate(p.signif = ifelse(nchar(p.adj.signif) > 3, "***", p.adj.signif)) %>%
  add_xy_position(x = "Group", dodge = dodge_width)

race_base_y <- df_summary %>%
  group_by(Sample, Group) %>%
  summarize(y.position = max(meanDprime + seDprime, na.rm = TRUE) + 0.08,
            .groups = "drop")

stat_test_race_all <- stat_test_race_all %>%
  select(-y.position) %>%
  left_join(race_base_y, by = c("Sample", "Group"))

cat("\nPaired race comparisons within each sample x group:\n")
print(
  stat_test_race_all %>%
    select(Sample, Group, group1, group2, n1, n2, statistic, df, p, p.adj, p.signif),
  n = Inf,
  width = Inf
)

stat_test_race <- stat_test_race_all %>%
  filter(p.signif != "ns")

# Independent TD vs CP comparisons within each sample x race.
stat_test_group <- df_grouped %>%
  group_by(Sample, ExperimentName) %>%
  t_test(Dprime ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  mutate(p.signif = ifelse(nchar(p.adj.signif) > 3, "***", p.adj.signif))

cat("\nIndependent TD vs CP comparisons within each sample x race:\n")
print(
  stat_test_group %>%
    select(Sample, ExperimentName, group1, group2, n1, n2, statistic, df, p, p.adj, p.signif),
  n = Inf,
  width = Inf
)

final_plot <- p +
  stat_pvalue_manual(
    stat_test_race,
    label = "p.signif",
    bracket.size = 0.6,
    tip.length = 0.01,
    size = 5
  )

if (interactive()) {
  print(final_plot)
}

ggsave(
  plot_path("biasm_dprime_young_and_all.png"),
  final_plot,
  width = 9,
  height = 5.5,
  dpi = 300,
  device = ragg::agg_png
)
