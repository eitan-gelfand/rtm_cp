# Analysis: CP/TD updating-model pipeline (last trail)
# Input:  data/updating_data_cp.csv
# Output: output/updating_beta_results.csv
#         output/last_trail_updating_analysis.png

library(lme4)
library(ggplot2)
library(dplyr)

.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

center_scale <- function(x) as.numeric(scale(x, scale = FALSE))
z_scale <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mu) / (2 * s))
}

race_label_map <- c("Caucasian" = "Own-Race", "Asian" = "Other-Race")
panel_label_map <- c(
  "TD_L" = "TD with low CFMT score",
  "TD_H" = "TD with high CFMT score",
  "CP" = "CP"
)

extract_term <- function(model_std, term_regex) {
  coef_tbl <- coef(summary(model_std))
  rn <- rownames(coef_tbl)
  idx <- grep(term_regex, rn)

  if (length(idx) == 0) {
    stop(paste0("Could not find coefficient matching regex: ", term_regex))
  }

  row <- coef_tbl[idx[1], , drop = FALSE]
  list(beta = as.numeric(row[1, "Estimate"]), se = as.numeric(row[1, "Std. Error"]))
}

fit_condition <- function(df_all, group_name, race_name, cfmt_level = NA_character_) {
  if (is.na(cfmt_level)) {
    df <- df_all %>% filter(.data$Group == group_name, .data$ExperimentName == race_name)
    panel_key <- "CP"
  } else {
    df <- df_all %>%
      filter(
        .data$Group == group_name,
        .data$ExperimentName == race_name,
        .data$CFMTgroup == cfmt_level
      )
    panel_key <- paste0("TD_", cfmt_level)
  }

  if (nrow(df) == 0) {
    stop(paste("No rows for condition:", group_name, race_name, cfmt_level))
  }

  df$c_New_tinf <- center_scale(df$New_tinf)
  df$c_New_t1 <- center_scale(df$New_t1)
  df$c_Range <- center_scale(df$Range)

  df$z_New_tinf <- z_scale(df$c_New_tinf)
  df$z_New_t1 <- z_scale(df$c_New_t1)
  df$z_Range <- z_scale(df$c_Range)

  model0 <- glmer(
    ACC ~ c_New_tinf + c_Range + (1 | Subject),
    data = df,
    family = binomial(link = "probit"),
    nAGQ = 0
  )

  model1 <- glmer(
    ACC ~ c_New_tinf + c_New_t1 + c_Range + (1 | Subject),
    data = df,
    family = binomial(link = "probit"),
    nAGQ = 0
  )

  model0_std <- glmer(
    ACC ~ z_New_tinf + z_Range + (1 | Subject),
    data = df,
    family = binomial(link = "probit"),
    nAGQ = 0
  )

  model1_std <- glmer(
    ACC ~ z_New_tinf + z_New_t1 + z_Range + (1 | Subject),
    data = df,
    family = binomial(link = "probit"),
    nAGQ = 0
  )

  cat("\n==============================\n")
  cat("Condition:", group_name, "|", race_name, "| CFMT:", ifelse(is.na(cfmt_level), "NA", cfmt_level), "\n")
  print(anova(model0, model1, test = "Chisq"))

  b0_tinf <- extract_term(model0_std, "z_New_tinf")
  b1_tinf <- extract_term(model1_std, "z_New_tinf")
  b1_t1 <- extract_term(model1_std, "z_New_t1")

  race_out <- unname(race_label_map[[race_name]])
  panel_out <- unname(panel_label_map[[panel_key]])

  rbind(
    data.frame(
      Group = group_name,
      CFMTgroup = ifelse(is.na(cfmt_level), "", cfmt_level),
      ExperimentName = race_name,
      Panel = panel_out,
      Race = race_out,
      Point = "D0_tinf",
      Label = "t-inf model",
      Beta = b0_tinf$beta,
      SE = b0_tinf$se,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Group = group_name,
      CFMTgroup = ifelse(is.na(cfmt_level), "", cfmt_level),
      ExperimentName = race_name,
      Panel = panel_out,
      Race = race_out,
      Point = "D1_t1",
      Label = "t-1",
      Beta = b1_t1$beta,
      SE = b1_t1$se,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Group = group_name,
      CFMTgroup = ifelse(is.na(cfmt_level), "", cfmt_level),
      ExperimentName = race_name,
      Panel = panel_out,
      Race = race_out,
      Point = "D1_tinf",
      Label = "t-inf",
      Beta = b1_tinf$beta,
      SE = b1_tinf$se,
      stringsAsFactors = FALSE
    )
  )
}

data <- read.csv(data_path("updating_data_cp.csv"))

data$Range <- as.numeric(data$Range)
data$ACC <- as.numeric(data$ACC)
data$New_t1 <- as.numeric(data$New_t1)
data$New_tinf <- as.numeric(data$New_tinf)

data <- subset(
  data,
  Range > 0 & !is.na(ACC) & !is.na(New_t1) & !is.na(New_tinf) & !is.na(Subject)
)

conditions <- data.frame(
  Group = c("TD", "TD", "TD", "TD", "CP", "CP"),
  ExperimentName = c("Caucasian", "Asian", "Caucasian", "Asian", "Caucasian", "Asian"),
  CFMTgroup = c("L", "L", "H", "H", NA, NA),
  stringsAsFactors = FALSE
)

results_list <- lapply(seq_len(nrow(conditions)), function(i) {
  fit_condition(
    df_all = data,
    group_name = conditions$Group[i],
    race_name = conditions$ExperimentName[i],
    cfmt_level = conditions$CFMTgroup[i]
  )
})

results <- bind_rows(results_list) %>%
  mutate(
    Panel = factor(Panel, levels = c("TD with low CFMT score", "TD with high CFMT score", "CP")),
    Race = factor(Race, levels = c("Own-Race", "Other-Race")),
    Point = factor(Point, levels = c("D0_tinf", "D1_t1", "D1_tinf")),
    Label = factor(Label, levels = c("t-inf model", "t-1", "t-inf"))
  ) %>%
  arrange(Panel, Race, Point)

write.csv(results, output_path("updating_beta_results.csv"), row.names = FALSE)
cat("Saved results table to:", output_path("updating_beta_results.csv"), "\n")

pal_race <- c("Own-Race" = "#A6D8C3", "Other-Race" = "#F6B48F")

p <- ggplot(results, aes(x = Point, y = Beta, color = Race)) +
  geom_errorbar(aes(ymin = Beta - SE, ymax = Beta + SE), width = 0.1, linewidth = 0.8) +
  geom_line(
    data = results %>% filter(.data$Point != "D0_tinf"),
    aes(group = Race),
    linewidth = 0.9
  ) +
  geom_point(size = 2.8) +
  facet_wrap(~ Panel, nrow = 1) +
  scale_x_discrete(labels = c("D0_tinf" = "t-inf model", "D1_t1" = "t-1", "D1_tinf" = "t-inf")) +
  scale_color_manual(values = pal_race, name = "Race") +
  labs(x = NULL, y = "Beta") +
  coord_cartesian(ylim = c(-0.1, 0.55)) +
  theme_pub() +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 12)
  )

ggsave(
  filename = output_path("last_trail_updating_analysis.png"),
  plot = p,
  width = 10,
  height = 4,
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("Saved plot to:", output_path("last_trail_updating_analysis.png"), "\n")