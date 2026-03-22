# Load required packages
library(glmmTMB)
library(dplyr)

# ── Shared helpers ─────────────────────────────────────────────
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))

# ──────────────────────────────────────────────
# Bootstrap controls
# ──────────────────────────────────────────────

# Keep the main controls together so scaling from a quick debug run
# to a larger bootstrap only requires changing a few values here.
bootstrap_seed       <- 20260321L
bootstrap_iterations <- 1000L
bootstrap_age_step   <- 1
bootstrap_save_every <- 4L
bootstrap_age_grid   <- seq(19, 78, by = bootstrap_age_step)
bootstrap_file_stem  <- "glmm_cp_peak_bootstrap_debug"

# ──────────────────────────────────────────────
# Read and preprocess data
# ──────────────────────────────────────────────

df <- read.csv(data_path("full_data_cp.csv")) %>%
  filter(Range != 0)

df <- df %>%
  mutate(
    Regression = factor(Regression, levels = c("biasm", "biasp")),
    ExperimentName = factor(
      ExperimentName,
      levels = c("Caucasian", "Asian"),
      labels = c("Own-Race", "Other-Race")
    ),
    Group = factor(Group, levels = c("TD", "CP")),
    Subject = factor(Subject)
  )

df$Range_c <- scale(df$Range, center = TRUE, scale = TRUE)

# ──────────────────────────────────────────────
# Fit the same official GLMM
# ──────────────────────────────────────────────

m_simple <- glmmTMB(
  ACC ~
    (poly(Age, 2) + Regression + ExperimentName + Group)^3 +
    Range_c +
    (1 | Subject),
  data   = df,
  family = binomial(link = "probit")
)

summary(m_simple)

# ──────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────

# Subject-level cluster bootstrap within Group. This is the preferred
# approach here because the random intercept is defined at the Subject level.
# When a subject is sampled multiple times, each copy gets a unique Subject ID
# so the refit treats duplicated clusters as separate bootstrap units.
subject_lookup <- df %>%
  distinct(Subject, Group) %>%
  arrange(Group, Subject)

subjects_by_group <- split(
  as.character(subject_lookup$Subject),
  as.character(subject_lookup$Group)
)

resample_subject_clusters <- function(data, subject_ids_by_group) {
  sampled_subjects <- lapply(
    subject_ids_by_group,
    function(ids) sample(ids, size = length(ids), replace = TRUE)
  )

  boot_chunks <- list()
  chunk_index <- 1L

  for (group_name in names(sampled_subjects)) {
    group_draws <- sampled_subjects[[group_name]]

    for (draw_index in seq_along(group_draws)) {
      subject_id <- group_draws[[draw_index]]

      subject_rows <- data %>%
        filter(as.character(Subject) == subject_id)

      subject_rows$Subject <- paste0(group_name, "_boot_", draw_index, "_", subject_id)
      boot_chunks[[chunk_index]] <- subject_rows
      chunk_index <- chunk_index + 1L
    }
  }

  boot_df <- bind_rows(boot_chunks)
  boot_df$Subject <- factor(boot_df$Subject)
  boot_df
}

make_bootstrap_prediction_grid <- function(data, age_grid) {
  pred_grid <- expand.grid(
    Age = age_grid,
    Regression = levels(data$Regression),
    ExperimentName = levels(data$ExperimentName),
    Group = levels(data$Group),
    Range_c = 0,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  pred_grid$Regression <- factor(pred_grid$Regression, levels = levels(data$Regression))
  pred_grid$ExperimentName <- factor(
    pred_grid$ExperimentName,
    levels = levels(data$ExperimentName)
  )
  pred_grid$Group <- factor(pred_grid$Group, levels = levels(data$Group))

  pred_grid
}

extract_group_peaks <- function(pred_grid) {
  pred_grid %>%
    group_by(Group, Age) %>%
    summarise(
      peak_pred = mean(pred),
      .groups = "drop"
    ) %>%
    arrange(Group, Age) %>%
    group_by(Group) %>%
    slice_max(order_by = peak_pred, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(peak_age = Age)
}

extract_group_race_peaks <- function(pred_grid) {
  pred_grid %>%
    group_by(Group, ExperimentName, Age) %>%
    summarise(
      peak_pred = mean(pred),
      .groups = "drop"
    ) %>%
    arrange(Group, ExperimentName, Age) %>%
    group_by(Group, ExperimentName) %>%
    slice_max(order_by = peak_pred, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(peak_age = Age)
}

extract_peak_differences <- function(peaks_df) {
  overall_td <- peaks_df %>%
    filter(peak_scope == "overall_group", Group == "TD") %>%
    transmute(td_peak_age = peak_age, td_peak_pred = peak_pred)

  overall_cp <- peaks_df %>%
    filter(peak_scope == "overall_group", Group == "CP") %>%
    transmute(cp_peak_age = peak_age, cp_peak_pred = peak_pred)

  if (nrow(overall_td) == 1 && nrow(overall_cp) == 1) {
    overall_diff <- bind_cols(overall_td, overall_cp) %>%
      mutate(
        comparison_scope = "overall_group",
        ExperimentName = NA_character_,
        cp_minus_td_peak_age = cp_peak_age - td_peak_age,
        cp_gt_td = cp_peak_age > td_peak_age,
        .before = 1
      )
  } else {
    overall_diff <- data.frame(
      comparison_scope = character(),
      ExperimentName = character(),
      td_peak_age = numeric(),
      td_peak_pred = numeric(),
      cp_peak_age = numeric(),
      cp_peak_pred = numeric(),
      cp_minus_td_peak_age = numeric(),
      cp_gt_td = logical(),
      stringsAsFactors = FALSE
    )
  }

  race_td <- peaks_df %>%
    filter(peak_scope == "group_race", Group == "TD") %>%
    transmute(
      ExperimentName = as.character(ExperimentName),
      td_peak_age = peak_age,
      td_peak_pred = peak_pred
    )

  race_cp <- peaks_df %>%
    filter(peak_scope == "group_race", Group == "CP") %>%
    transmute(
      ExperimentName = as.character(ExperimentName),
      cp_peak_age = peak_age,
      cp_peak_pred = peak_pred
    )

  race_diff <- merge(race_td, race_cp, by = "ExperimentName") %>%
    mutate(
      comparison_scope = "group_race",
      cp_minus_td_peak_age = cp_peak_age - td_peak_age,
      cp_gt_td = cp_peak_age > td_peak_age,
      .before = 1
    )

  bind_rows(overall_diff, race_diff) %>%
    select(
      comparison_scope,
      ExperimentName,
      td_peak_age,
      td_peak_pred,
      cp_peak_age,
      cp_peak_pred,
      cp_minus_td_peak_age,
      cp_gt_td
    )
}

percentile_stat <- function(x, prob) {
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }

  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

summarise_peak_results <- function(results_df, grouping_vars) {
  results_df %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      successful_iterations = n(),
      mean_peak_age = mean(peak_age, na.rm = TRUE),
      sd_peak_age = stats::sd(peak_age, na.rm = TRUE),
      median_peak_age = stats::median(peak_age, na.rm = TRUE),
      peak_age_p2_5 = percentile_stat(peak_age, 0.025),
      peak_age_p97_5 = percentile_stat(peak_age, 0.975),
      mean_peak_pred = mean(peak_pred, na.rm = TRUE),
      sd_peak_pred = stats::sd(peak_pred, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_difference_results <- function(differences_df) {
  differences_df %>%
    group_by(comparison_scope, ExperimentName) %>%
    summarise(
      successful_iterations = n(),
      mean_cp_minus_td_peak_age = mean(cp_minus_td_peak_age, na.rm = TRUE),
      sd_cp_minus_td_peak_age = stats::sd(cp_minus_td_peak_age, na.rm = TRUE),
      median_cp_minus_td_peak_age = stats::median(cp_minus_td_peak_age, na.rm = TRUE),
      p2_5_cp_minus_td_peak_age = percentile_stat(cp_minus_td_peak_age, 0.025),
      p97_5_cp_minus_td_peak_age = percentile_stat(cp_minus_td_peak_age, 0.975),
      prop_cp_peak_gt_td_peak = mean(cp_gt_td, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      comparison_label = ifelse(
        comparison_scope == "overall_group",
        "Overall",
        as.character(ExperimentName)
      ),
      .after = ExperimentName
    )
}

save_bootstrap_checkpoint <- function(results_df, differences_df, failures_df, file_stem) {
  checkpoint_rds <- output_path(paste0(file_stem, "_checkpoint.rds"))
  checkpoint_csv <- output_path(paste0(file_stem, "_results_checkpoint.csv"))
  diff_csv       <- output_path(paste0(file_stem, "_peak_differences_checkpoint.csv"))
  failures_csv   <- output_path(paste0(file_stem, "_failures_checkpoint.csv"))

  saveRDS(
    list(
      results = results_df,
      peak_differences = differences_df,
      failures = failures_df
    ),
    checkpoint_rds
  )
  write.csv(results_df, checkpoint_csv, row.names = FALSE)
  write.csv(differences_df, diff_csv, row.names = FALSE)
  write.csv(failures_df, failures_csv, row.names = FALSE)
}

empty_results_df <- function() {
  data.frame(
    iteration = integer(),
    peak_scope = character(),
    Group = character(),
    ExperimentName = character(),
    peak_age = numeric(),
    peak_pred = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_differences_df <- function() {
  data.frame(
    iteration = integer(),
    comparison_scope = character(),
    ExperimentName = character(),
    td_peak_age = numeric(),
    td_peak_pred = numeric(),
    cp_peak_age = numeric(),
    cp_peak_pred = numeric(),
    cp_minus_td_peak_age = numeric(),
    cp_gt_td = logical(),
    stringsAsFactors = FALSE
  )
}

empty_failures_df <- function() {
  data.frame(
    iteration = integer(),
    error_message = character(),
    stringsAsFactors = FALSE
  )
}

empty_group_summary_df <- function() {
  data.frame(
    Group = character(),
    successful_iterations = integer(),
    mean_peak_age = numeric(),
    sd_peak_age = numeric(),
    median_peak_age = numeric(),
    peak_age_p2_5 = numeric(),
    peak_age_p97_5 = numeric(),
    mean_peak_pred = numeric(),
    sd_peak_pred = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_group_race_summary_df <- function() {
  data.frame(
    Group = character(),
    ExperimentName = character(),
    successful_iterations = integer(),
    mean_peak_age = numeric(),
    sd_peak_age = numeric(),
    median_peak_age = numeric(),
    peak_age_p2_5 = numeric(),
    peak_age_p97_5 = numeric(),
    mean_peak_pred = numeric(),
    sd_peak_pred = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_difference_summary_df <- function() {
  data.frame(
    comparison_scope = character(),
    ExperimentName = character(),
    comparison_label = character(),
    successful_iterations = integer(),
    mean_cp_minus_td_peak_age = numeric(),
    sd_cp_minus_td_peak_age = numeric(),
    median_cp_minus_td_peak_age = numeric(),
    p2_5_cp_minus_td_peak_age = numeric(),
    p97_5_cp_minus_td_peak_age = numeric(),
    prop_cp_peak_gt_td_peak = numeric(),
    stringsAsFactors = FALSE
  )
}

format_mean_sd <- function(mean_value, sd_value) {
  if (length(mean_value) == 0 || length(sd_value) == 0) {
    return("NA")
  }

  if (any(is.na(c(mean_value, sd_value)))) {
    return("NA")
  }

  sprintf("%.2f +/- %.2f", mean_value, sd_value)
}

format_diff_ci <- function(mean_value, lower_value, upper_value) {
  if (length(mean_value) == 0 || length(lower_value) == 0 || length(upper_value) == 0) {
    return("NA")
  }

  if (any(is.na(c(mean_value, lower_value, upper_value)))) {
    return("NA")
  }

  sprintf("%.2f (95%% percentile CI %.2f to %.2f)", mean_value, lower_value, upper_value)
}

format_proportion <- function(value) {
  if (length(value) == 0) {
    return("NA")
  }

  if (is.na(value)) {
    return("NA")
  }

  sprintf("%.3f", value)
}

# ──────────────────────────────────────────────
# Bootstrap loop
# ──────────────────────────────────────────────

set.seed(bootstrap_seed)

bootstrap_results_list <- vector("list", bootstrap_iterations)
bootstrap_differences_list <- vector("list", bootstrap_iterations)
bootstrap_failures_list <- vector("list", bootstrap_iterations)
successful_bootstrap_fits <- 0L
bootstrap_start_time <- Sys.time()

cat(
  sprintf(
    "\nBOOTSTRAP PEAK-AGE RUN\nSeed: %d | Iterations: %d | Age grid: %s to %s by %s\n",
    bootstrap_seed,
    bootstrap_iterations,
    min(bootstrap_age_grid),
    max(bootstrap_age_grid),
    bootstrap_age_step
  )
)

for (iter in seq_len(bootstrap_iterations)) {
  iter_start_time <- Sys.time()
  cat(sprintf("\n[Bootstrap] Iteration %d/%d started\n", iter, bootstrap_iterations))

  iter_result <- tryCatch(
    {
      boot_df <- resample_subject_clusters(df, subjects_by_group)
      boot_fit <- update(m_simple, data = boot_df)

      boot_pred_grid <- make_bootstrap_prediction_grid(boot_df, bootstrap_age_grid)
      boot_pred_grid$pred <- predict(
        boot_fit,
        newdata = boot_pred_grid,
        type = "response",
        re.form = NA
      )

      boot_group_peaks <- extract_group_peaks(boot_pred_grid)
      boot_group_race_peaks <- extract_group_race_peaks(boot_pred_grid)

      if (!all(levels(df$Group) %in% as.character(boot_group_peaks$Group))) {
        stop("Peak extraction did not return both groups.")
      }

      expected_group_race <- expand.grid(
        Group = levels(df$Group),
        ExperimentName = levels(df$ExperimentName),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )

      observed_group_race <- boot_group_race_peaks %>%
        transmute(
          Group = as.character(Group),
          ExperimentName = as.character(ExperimentName)
        )

      missing_group_race <- expected_group_race %>%
        anti_join(observed_group_race, by = c("Group", "ExperimentName"))

      if (nrow(missing_group_race) > 0) {
        stop("Peak extraction did not return all Group x Race combinations.")
      }

      boot_peaks <- bind_rows(
        boot_group_peaks %>%
          mutate(
            peak_scope = "overall_group",
            ExperimentName = NA_character_,
            .before = 1
          ),
        boot_group_race_peaks %>%
          mutate(
            peak_scope = "group_race",
            .before = 1
          )
      ) %>%
        mutate(iteration = iter, .before = 1)

      boot_differences <- extract_peak_differences(boot_peaks) %>%
        mutate(iteration = iter, .before = 1)

      race_difference_levels <- boot_differences %>%
        filter(comparison_scope == "group_race") %>%
        pull(ExperimentName) %>%
        as.character()

      if (
        sum(boot_differences$comparison_scope == "overall_group") != 1L ||
        !all(levels(df$ExperimentName) %in% race_difference_levels)
      ) {
        stop("Peak-difference extraction did not return the expected comparisons.")
      }

      list(
        status = "success",
        peaks = boot_peaks,
        peak_differences = boot_differences,
        error_message = NA_character_
      )
    },
    error = function(e) {
      list(
        status = "failure",
        peaks = NULL,
        error_message = conditionMessage(e)
      )
    }
  )

  iter_elapsed_seconds <- as.numeric(
    difftime(Sys.time(), iter_start_time, units = "secs")
  )

  if (identical(iter_result$status, "success")) {
    successful_bootstrap_fits <- successful_bootstrap_fits + 1L
    bootstrap_results_list[[iter]] <- iter_result$peaks
    bootstrap_differences_list[[iter]] <- iter_result$peak_differences

    cat(
      sprintf(
        "[Bootstrap] Iteration %d/%d success | elapsed %.1f s | successful fits %d\n",
        iter,
        bootstrap_iterations,
        iter_elapsed_seconds,
        successful_bootstrap_fits
      )
    )
  } else {
    bootstrap_failures_list[[iter]] <- data.frame(
      iteration = iter,
      error_message = iter_result$error_message,
      stringsAsFactors = FALSE
    )

    cat(
      sprintf(
        "[Bootstrap] Iteration %d/%d failure | elapsed %.1f s | successful fits %d\n",
        iter,
        bootstrap_iterations,
        iter_elapsed_seconds,
        successful_bootstrap_fits
      )
    )
    cat(sprintf("[Bootstrap] Error: %s\n", iter_result$error_message))
  }

  if (iter %% bootstrap_save_every == 0L || iter == bootstrap_iterations) {
    checkpoint_results <- bind_rows(bootstrap_results_list)
    checkpoint_differences <- bind_rows(bootstrap_differences_list)
    checkpoint_failures <- bind_rows(bootstrap_failures_list)

    if (nrow(checkpoint_results) == 0) {
      checkpoint_results <- empty_results_df()
    }

    if (nrow(checkpoint_differences) == 0) {
      checkpoint_differences <- empty_differences_df()
    }

    if (nrow(checkpoint_failures) == 0) {
      checkpoint_failures <- empty_failures_df()
    }

    save_bootstrap_checkpoint(
      results_df = checkpoint_results,
      differences_df = checkpoint_differences,
      failures_df = checkpoint_failures,
      file_stem = bootstrap_file_stem
    )

    cat(sprintf("[Bootstrap] Checkpoint saved after iteration %d\n", iter))
  }
}

# ──────────────────────────────────────────────
# Collect and summarise results
# ──────────────────────────────────────────────

bootstrap_results <- bind_rows(bootstrap_results_list)
bootstrap_peak_differences <- bind_rows(bootstrap_differences_list)
bootstrap_failures <- bind_rows(bootstrap_failures_list)

if (nrow(bootstrap_results) == 0) {
  bootstrap_results <- empty_results_df()
}

if (nrow(bootstrap_peak_differences) == 0) {
  bootstrap_peak_differences <- empty_differences_df()
}

if (nrow(bootstrap_failures) == 0) {
  bootstrap_failures <- empty_failures_df()
}

bootstrap_results <- bootstrap_results %>%
  arrange(iteration, peak_scope, ExperimentName, Group)

bootstrap_peak_differences <- bootstrap_peak_differences %>%
  arrange(iteration, comparison_scope, ExperimentName)

bootstrap_group_results <- bootstrap_results %>%
  filter(peak_scope == "overall_group")

bootstrap_group_race_results <- bootstrap_results %>%
  filter(peak_scope == "group_race")

if (nrow(bootstrap_group_results) > 0) {
  bootstrap_summary_by_group <- summarise_peak_results(
    results_df = bootstrap_group_results,
    grouping_vars = "Group"
  ) %>%
    arrange(match(Group, levels(df$Group)))
} else {
  bootstrap_summary_by_group <- empty_group_summary_df()
}

if (nrow(bootstrap_group_race_results) > 0) {
  bootstrap_summary_by_group_race <- summarise_peak_results(
    results_df = bootstrap_group_race_results,
    grouping_vars = c("Group", "ExperimentName")
  ) %>%
    arrange(
      match(Group, levels(df$Group)),
      match(ExperimentName, levels(df$ExperimentName))
    )
} else {
  bootstrap_summary_by_group_race <- empty_group_race_summary_df()
}

if (nrow(bootstrap_peak_differences) > 0) {
  bootstrap_difference_summary <- summarise_difference_results(
    bootstrap_peak_differences
  ) %>%
    arrange(
      match(comparison_scope, c("overall_group", "group_race")),
      match(ExperimentName, levels(df$ExperimentName))
    )
} else {
  bootstrap_difference_summary <- empty_difference_summary_df()
}

bootstrap_total_elapsed <- as.numeric(
  difftime(Sys.time(), bootstrap_start_time, units = "secs")
)

bootstrap_summary <- list(
  seed = bootstrap_seed,
  iterations_requested = bootstrap_iterations,
  age_grid = bootstrap_age_grid,
  successful_iterations = successful_bootstrap_fits,
  failed_iterations = nrow(bootstrap_failures),
  by_group = bootstrap_summary_by_group,
  by_group_race = bootstrap_summary_by_group_race,
  peak_age_difference = bootstrap_difference_summary,
  total_elapsed_seconds = bootstrap_total_elapsed
)

# ──────────────────────────────────────────────
# Save outputs
# ──────────────────────────────────────────────

bootstrap_results_rds <- output_path(paste0(bootstrap_file_stem, "_results.rds"))
bootstrap_results_csv <- output_path(paste0(bootstrap_file_stem, "_results.csv"))
bootstrap_differences_csv <- output_path(paste0(bootstrap_file_stem, "_peak_differences.csv"))
bootstrap_failures_csv <- output_path(paste0(bootstrap_file_stem, "_failures.csv"))
bootstrap_summary_rds <- output_path(paste0(bootstrap_file_stem, "_summary.rds"))
bootstrap_summary_csv <- output_path(paste0(bootstrap_file_stem, "_summary_by_group.csv"))
bootstrap_group_race_summary_csv <- output_path(paste0(bootstrap_file_stem, "_summary_by_group_race.csv"))
bootstrap_difference_csv <- output_path(paste0(bootstrap_file_stem, "_peak_age_difference_summary.csv"))

saveRDS(
  list(
    results = bootstrap_results,
    peak_differences = bootstrap_peak_differences,
    failures = bootstrap_failures,
    summary = bootstrap_summary,
    by_group = bootstrap_summary_by_group,
    by_group_race = bootstrap_summary_by_group_race,
    peak_age_difference = bootstrap_difference_summary
  ),
  bootstrap_results_rds
)

write.csv(bootstrap_results, bootstrap_results_csv, row.names = FALSE)
write.csv(bootstrap_peak_differences, bootstrap_differences_csv, row.names = FALSE)
write.csv(bootstrap_failures, bootstrap_failures_csv, row.names = FALSE)
saveRDS(bootstrap_summary, bootstrap_summary_rds)
write.csv(bootstrap_summary_by_group, bootstrap_summary_csv, row.names = FALSE)
write.csv(bootstrap_summary_by_group_race, bootstrap_group_race_summary_csv, row.names = FALSE)
write.csv(bootstrap_difference_summary, bootstrap_difference_csv, row.names = FALSE)

# ──────────────────────────────────────────────
# Print final summary
# ──────────────────────────────────────────────

cat("\nBOOTSTRAP PEAK-AGE SUMMARY\n")
cat(sprintf("Successful iterations: %d\n", successful_bootstrap_fits))
cat(sprintf("Failed iterations: %d\n", nrow(bootstrap_failures)))

td_summary <- bootstrap_summary_by_group %>% filter(Group == "TD")
cp_summary <- bootstrap_summary_by_group %>% filter(Group == "CP")
overall_difference_summary <- bootstrap_difference_summary %>%
  filter(comparison_scope == "overall_group")
own_race_difference_summary <- bootstrap_difference_summary %>%
  filter(comparison_scope == "group_race", ExperimentName == "Own-Race")
other_race_difference_summary <- bootstrap_difference_summary %>%
  filter(comparison_scope == "group_race", ExperimentName == "Other-Race")

if (nrow(td_summary) == 1) {
  cat(
    sprintf(
      "TD peak age: %s\n",
      format_mean_sd(td_summary$mean_peak_age, td_summary$sd_peak_age)
    )
  )
}

if (nrow(cp_summary) == 1) {
  cat(
    sprintf(
      "CP peak age: %s\n",
      format_mean_sd(cp_summary$mean_peak_age, cp_summary$sd_peak_age)
    )
  )
}

if (nrow(overall_difference_summary) == 1) {
  cat(
    sprintf(
      "Overall CP - TD peak age difference: %s\n",
      format_diff_ci(
        overall_difference_summary$mean_cp_minus_td_peak_age,
        overall_difference_summary$p2_5_cp_minus_td_peak_age,
        overall_difference_summary$p97_5_cp_minus_td_peak_age
      )
    )
  )
}

if (nrow(own_race_difference_summary) == 1) {
  cat(
    sprintf(
      "Own-Race CP - TD peak age difference: %s\n",
      format_diff_ci(
        own_race_difference_summary$mean_cp_minus_td_peak_age,
        own_race_difference_summary$p2_5_cp_minus_td_peak_age,
        own_race_difference_summary$p97_5_cp_minus_td_peak_age
      )
    )
  )
}

if (nrow(other_race_difference_summary) == 1) {
  cat(
    sprintf(
      "Other-Race CP - TD peak age difference: %s\n",
      format_diff_ci(
        other_race_difference_summary$mean_cp_minus_td_peak_age,
        other_race_difference_summary$p2_5_cp_minus_td_peak_age,
        other_race_difference_summary$p97_5_cp_minus_td_peak_age
      )
    )
  )
}

cat(
  sprintf(
    "Proportion with CP > TD peak age | Overall: %s | Own-Race: %s | Other-Race: %s\n",
    format_proportion(overall_difference_summary$prop_cp_peak_gt_td_peak),
    format_proportion(own_race_difference_summary$prop_cp_peak_gt_td_peak),
    format_proportion(other_race_difference_summary$prop_cp_peak_gt_td_peak)
  )
)
cat(sprintf("Bootstrap total elapsed time: %.1f s\n", bootstrap_total_elapsed))
cat("Saved outputs:\n")
cat(sprintf("  Raw results RDS: %s\n", bootstrap_results_rds))
cat(sprintf("  Raw peaks CSV: %s\n", bootstrap_results_csv))
cat(sprintf("  Raw differences CSV: %s\n", bootstrap_differences_csv))
cat(sprintf("  Failures CSV: %s\n", bootstrap_failures_csv))
cat(sprintf("  Group summary CSV: %s\n", bootstrap_summary_csv))
cat(sprintf("  Group x Race summary CSV: %s\n", bootstrap_group_race_summary_csv))
cat(sprintf("  Difference summary CSV: %s\n", bootstrap_difference_csv))
