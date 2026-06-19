# =========================================================
# POC mixed linear models for d-prime and RTM magnitude
# - Uses bias- d-prime only, matching dprime_rtm_corr.R
# - Computes Magnitude = bias+ accuracy - bias- accuracy
# - Tests the full Dprime x Group x Race interaction
# - Prints subject-level and exploratory trial-weighted summaries
# =========================================================

required_packages <- c("dplyr", "tidyr", "glmmTMB", "emmeans")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(tidyr)
library(glmmTMB)
library(emmeans)

.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))

report_path <- output_path("rtm_dprime_lmm_poc_summary.txt")

emit <- function(...) {
  text <- paste(..., sep = "")
  cat(text, "\n", sep = "")
  cat(text, "\n", sep = "", file = report_path, append = TRUE)
}

emit_object <- function(object) {
  lines <- capture.output(print(object))
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  cat(paste(lines, collapse = "\n"), "\n", sep = "", file = report_path, append = TRUE)
}

emit_section <- function(title) {
  emit("")
  emit("=========================================================")
  emit(title)
  emit("=========================================================")
}

safe_p <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) {
    return(NA_real_)
  }
  cor.test(x[ok], y[ok], method = "pearson")$p.value
}

fit_and_report <- function(data, label) {
  emit_section(label)
  emit("Rows: ", nrow(data))
  emit("Subjects: ", dplyr::n_distinct(data$Subject))

  full_model <- glmmTMB(
    Magnitude ~ Dprime_z * Group * Race + (1 | Subject),
    data = data,
    family = gaussian(),
    REML = FALSE
  )

  no_three_way_model <- glmmTMB(
    Magnitude ~ (Dprime_z + Group + Race)^2 + (1 | Subject),
    data = data,
    family = gaussian(),
    REML = FALSE
  )

  emit_section(paste(label, "- full model summary"))
  emit_object(summary(full_model))

  emit_section(paste(label, "- likelihood-ratio test for 3-way interaction"))
  emit("Full:      Magnitude ~ Dprime_z * Group * Race + (1 | Subject)")
  emit("Reduced:   Magnitude ~ (Dprime_z + Group + Race)^2 + (1 | Subject)")
  emit_object(anova(no_three_way_model, full_model))

  emit_section(paste(label, "- estimated Dprime slopes by Group x Race"))
  slopes <- emtrends(full_model, specs = ~ Group * Race, var = "Dprime_z")
  emit_object(summary(slopes))

  emit_section(paste(label, "- pairwise slope comparisons"))
  emit_object(pairs(slopes))

  invisible(full_model)
}

if (file.exists(report_path)) {
  invisible(file.remove(report_path))
}

df <- read.csv(data_path("full_data_cp.csv"))
dprime_df <- read.csv(data_path("dprime_results.csv"))

df_clean <- df %>%
  filter(Range != 0)

dprime_biasm <- dprime_df %>%
  filter(Regression == "biasm") %>%
  group_by(Subject, ExperimentName) %>%
  summarise(
    Dprime = mean(Dprime, na.rm = TRUE),
    .groups = "drop"
  )

subject_condition_acc <- df_clean %>%
  group_by(Subject, Group, ExperimentName, Regression) %>%
  summarise(
    mean_acc = mean(ACC, na.rm = TRUE),
    n_trials = n(),
    .groups = "drop"
  )

subject_model_data <- subject_condition_acc %>%
  pivot_wider(
    names_from = Regression,
    values_from = c(mean_acc, n_trials),
    names_sep = "_"
  ) %>%
  mutate(
    biasM = mean_acc_biasm,
    biasP = mean_acc_biasp,
    Magnitude = biasP - biasM
  ) %>%
  inner_join(dprime_biasm, by = c("Subject", "ExperimentName")) %>%
  mutate(
    Subject = factor(Subject),
    Group = factor(Group, levels = c("TD", "CP")),
    Race = factor(
      ExperimentName,
      levels = c("Caucasian", "Asian"),
      labels = c("Own-Race", "Other-Race")
    )
  ) %>%
  filter(complete.cases(Magnitude, Dprime, Group, Race)) %>%
  mutate(
    Dprime_z = as.numeric(scale(Dprime))
  )

trial_weighted_data <- df_clean %>%
  mutate(
    Subject = factor(Subject),
    Group = factor(Group, levels = c("TD", "CP")),
    Race = factor(
      ExperimentName,
      levels = c("Caucasian", "Asian"),
      labels = c("Own-Race", "Other-Race")
    )
  ) %>%
  inner_join(
    subject_model_data %>%
      select(Subject, ExperimentName, Dprime, Dprime_z, Magnitude),
    by = c("Subject", "ExperimentName")
  ) %>%
  filter(complete.cases(Magnitude, Dprime_z, Group, Race))

cell_summary <- subject_model_data %>%
  group_by(Group, Race) %>%
  summarise(
    n_subjects = n_distinct(Subject),
    mean_dprime = mean(Dprime, na.rm = TRUE),
    sd_dprime = sd(Dprime, na.rm = TRUE),
    mean_magnitude = mean(Magnitude, na.rm = TRUE),
    sd_magnitude = sd(Magnitude, na.rm = TRUE),
    r = cor(Dprime, Magnitude, use = "complete.obs"),
    p = safe_p(Dprime, Magnitude),
    .groups = "drop"
  )

emit_section("Cell summary from subject-level data")
emit("Race labels: Caucasian = Own-Race, Asian = Other-Race")
emit("Dprime predictor: bias- trials only")
emit("Magnitude outcome: mean ACC bias+ - mean ACC bias-")
emit_object(as.data.frame(cell_summary))

subject_model <- fit_and_report(
  subject_model_data,
  "SUBJECT-LEVEL LMM"
)

emit_section("Exploratory trial-weighted note")
emit(
  "The trial-weighted model repeats each subject x race Magnitude value ",
  "across all non-practice trials. Treat it as a POC weighting check, ",
  "not the clean inferential model."
)

trial_weighted_model <- fit_and_report(
  trial_weighted_data,
  "EXPLORATORY TRIAL-WEIGHTED LMM"
)

emit_section("Output")
emit("Wrote summary to: ", report_path)
