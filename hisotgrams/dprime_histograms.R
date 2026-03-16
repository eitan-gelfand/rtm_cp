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
# 1) Load CP Dprime dataset
# ──────────────────────────────────────────────────────────────
df <- read.csv(data_path("dprime_results_with_range.csv"))
  
# ──────────────────────────────────────────────────────────────
# 2) Aggregate over Range (since histogram uses subject-level)
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
# ──────────────────────────────────────────────────────────────
df_summary <- df_grouped %>%
  group_by(ExperimentName, Regression, Group) %>%
  summarize(
    meanDprime = mean(Dprime, na.rm = TRUE),
    sdDprime   = sd(Dprime,  na.rm = TRUE),
    n          = n(),
    seDprime   = sdDprime / sqrt(n),
    .groups    = "drop"
  )


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
  labs(x = "Group", y = "D-prime") +
  scale_fill_manual(
    name   = "Regression",
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias-")
  ) +
  scale_y_continuous(labels = number_format(accuracy = 0.01)) +
  theme_pub()

# ──────────────────────────────────────────────────────────────
# 5) Bias+ vs Bias– (paired test)
# ──────────────────────────────────────────────────────────────
stat_test_reg <- df_grouped %>%
  group_by(ExperimentName) %>%
  t_test(Dprime ~ Regression, paired = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif)) %>%
  add_xy_position(x = "Group", dodge = 0.5)

# Place regression brackets above bars
df_max_reg <- df_summary %>%
  group_by(ExperimentName) %>%
  summarize(max_bar = max(meanDprime), .groups = "drop")

stat_test_reg <- stat_test_reg %>%
  left_join(df_max_reg, by = "ExperimentName") %>%
  mutate(y.position = max_bar + 0.15)

# ──────────────────────────────────────────────────────────────
# 6) TD vs CP comparison
# ──────────────────────────────────────────────────────────────
stat_test_group <- df_grouped %>%
  group_by(ExperimentName, Regression) %>%
  t_test(Dprime ~ Group, paired = FALSE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p") %>%
  mutate(p.signif = ifelse(nchar(p.signif) > 3, "***", p.signif))

# get y positions
df_max_group <- df_summary %>%
  group_by(ExperimentName, Regression) %>%
  summarize(max_bar = max(meanDprime), .groups = "drop")

stat_test_group <- stat_test_group %>%
  left_join(df_max_group,
            by = c("ExperimentName", "Regression")) %>%
  mutate(y.position = max_bar + 0.35)

# ──────────────────────────────────────────────────────────────
# 7) Add significance brackets to final CP plot
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
  plot_path("dprime_histogram.png"),
  final_plot_cp,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)
