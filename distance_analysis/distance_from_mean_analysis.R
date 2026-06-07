# Distance-from-mean control analyses for the CP RTM revision.
#
# Goal:
#   Address the concern that the Bias+/Bias- contrast may be driven by the
#   second face's distance from the morph midpoint rather than by RTM.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(emmeans)
})

.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))

output_file <- file.path(.PROJECT_ROOT, "distance_analysis", "distance_from_mean_analysis_output.txt")
if (file.exists(output_file)) file.remove(output_file)

write_report <- function(...) {
  text <- paste(..., sep = "")
  cat(text, "\n", sep = "")
  cat(text, "\n", sep = "", file = output_file, append = TRUE)
}

write_block <- function(x = "") {
  lines <- capture.output(print(x))
  cat(lines, sep = "\n")
  cat("\n")
  cat(lines, sep = "\n", file = output_file, append = TRUE)
  cat("\n", file = output_file, append = TRUE)
}

write_section <- function(title) {
  rule <- paste(rep("=", nchar(title)), collapse = "")
  write_report("\n", title, "\n", rule)
}

print_model_summary <- function(model, title) {
  write_section(title)
  write_report("\nModel summary")
  cat(capture.output(summary(model)), sep = "\n")
  cat("\n")
  cat(capture.output(summary(model)), sep = "\n", file = output_file, append = TRUE)
  cat("\n", file = output_file, append = TRUE)

  write_report("\nType-III fixed-effect tests")
  if (requireNamespace("car", quietly = TRUE)) {
    write_block(car::Anova(model, type = 3))
  } else {
    write_report("Package 'car' is not installed; skipping Type-III tests.")
  }
}

options(contrasts = c("contr.sum", "contr.poly"))

write_section("Distance-from-mean control analyses")
write_report("Input file: ", data_path("updating_data_cp.csv"))
write_report("Output file: ", output_file)

raw_df <- read.csv(data_path("updating_data_cp.csv"))

write_section("Trial counts before cleaning")
write_report("Total rows: ", nrow(raw_df))
write_report("Rows with RT <= 250: ", sum(raw_df$RT <= 250, na.rm = TRUE))
write_report("Rows with missing RT: ", sum(is.na(raw_df$RT)))

df <- raw_df %>%
  filter(!is.na(RT), RT > 250) %>%
  mutate(
    Subject = factor(Subject),
    Group = factor(Group, levels = c("TD", "CP")),
    ExperimentName = factor(
      ExperimentName,
      levels = c("Caucasian", "Asian"),
      labels = c("Own-Race", "Other-Race")
    ),
    Condition = factor(Condition, levels = c("Same", "Different")),
    Regression = recode(Regression, "biasP" = "biasp", "biasM" = "biasm"),
    Regression = factor(Regression, levels = c("biasm", "biasp")),
    F1_distance_50 = abs(Value1 - 50),
    F2_distance_50 = abs(Value2 - 50),
    F1_distance_50_c = as.numeric(scale(F1_distance_50, center = TRUE, scale = TRUE)),
    F2_distance_50_c = as.numeric(scale(F2_distance_50, center = TRUE, scale = TRUE)),
    Range_c = as.numeric(scale(Range, center = TRUE, scale = TRUE)),
    Age_c = as.numeric(scale(Age, center = TRUE, scale = FALSE))
  )

same_df <- df %>% filter(Condition == "Same")
different_df <- df %>% filter(Condition == "Different", !is.na(Regression))

write_section("Trial counts after RT cleaning")
write_report("Total rows retained: ", nrow(df))
write_report("Same rows retained: ", nrow(same_df))
write_report("Different rows retained: ", nrow(different_df))

write_report("\nRows by Condition x Group x Race")
write_block(
  df %>%
    group_by(Condition, Group, ExperimentName) %>%
    summarise(n = n(), .groups = "drop")
)

write_report("\nDifferent-trial rows by Bias x Group x Race")
write_block(
  different_df %>%
    group_by(Regression, Group, ExperimentName) %>%
    summarise(n = n(), .groups = "drop")
)

write_section("Distance summaries")
write_report("Same trials: F1 distance from midpoint by Group x Race")
write_block(
  same_df %>%
    group_by(Group, ExperimentName) %>%
    summarise(
      n = n(),
      mean_F1_distance_50 = mean(F1_distance_50, na.rm = TRUE),
      sd_F1_distance_50 = sd(F1_distance_50, na.rm = TRUE),
      .groups = "drop"
    )
)

write_report("\nDifferent trials: F2 distance from midpoint by Bias x Group x Race")
write_block(
  different_df %>%
    group_by(Regression, Group, ExperimentName) %>%
    summarise(
      n = n(),
      mean_F2_distance_50 = mean(F2_distance_50, na.rm = TRUE),
      sd_F2_distance_50 = sd(F2_distance_50, na.rm = TRUE),
      .groups = "drop"
    )
)

write_report("\nPoint-biserial correlation: F2 distance x Bias+ indicator")
bias_indicator <- ifelse(different_df$Regression == "biasp", 1, 0)
write_report(
  "cor(F2_distance_50, Bias+ indicator) = ",
  round(cor(different_df$F2_distance_50, bias_indicator, use = "complete.obs"), 4)
)

same_model <- glmmTMB(
  ACC ~ F1_distance_50_c * Group * ExperimentName + Age_c + (1 | Subject),
  data = same_df,
  family = binomial(link = "probit")
)

different_model <- glmmTMB(
  ACC ~ Range_c + Regression * Group * ExperimentName + F2_distance_50_c + Age_c + (1 | Subject),
  data = different_df,
  family = binomial(link = "probit")
)

print_model_summary(
  same_model,
  "Same-trial model: accuracy as a function of first-face distance from midpoint"
)

print_model_summary(
  different_model,
  "Different-trial model: Bias effect while controlling second-face distance from midpoint"
)

write_section("Estimated Bias contrasts from different-trial model")
write_report("Bias+ vs Bias- within Group x Race, controlling Range, F2 distance, and Age.")
write_block(
  pairs(emmeans(different_model, ~ Regression | Group * ExperimentName))
)

write_section("Done")
write_report("Analysis complete.")
