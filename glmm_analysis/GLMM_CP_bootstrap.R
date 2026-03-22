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
bootstrap_iterations <- 10L
bootstrap_age_step   <- 1
bootstrap_save_every <- 2L
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

save_bootstrap_checkpoint <- function(results_df, failures_df, file_stem) {
  checkpoint_rds <- output_path(paste0(file_stem, "_checkpoint.rds"))
  checkpoint_csv <- output_path(paste0(file_stem, "_results_checkpoint.csv"))
  failures_csv   <- output_path(paste0(file_stem, "_failures_checkpoint.csv"))

  saveRDS(
    list(results = results_df, failures = failures_df),
    checkpoint_rds
  )
  write.csv(results_df, checkpoint_csv, row.names = FALSE)
  write.csv(failures_df, failures_csv, row.names = FALSE)
}

empty_results_df <- function() {
  data.frame(
    iteration = integer(),
    Group = character(),
    peak_age = numeric(),
    peak_pred = numeric(),
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

# ──────────────────────────────────────────────
# Bootstrap loop
# ──────────────────────────────────────────────

set.seed(bootstrap_seed)

bootstrap_results_list <- vector("list", bootstrap_iterations)
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

      boot_peaks <- extract_group_peaks(boot_pred_grid) %>%
        mutate(iteration = iter, .before = 1)

      if (!all(levels(df$Group) %in% as.character(boot_peaks$Group))) {
        stop("Peak extraction did not return both groups.")
      }

      list(
        status = "success",
        peaks = boot_peaks,
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
    checkpoint_failures <- bind_rows(bootstrap_failures_list)

    if (nrow(checkpoint_results) == 0) {
      checkpoint_results <- empty_results_df()
    }

    if (nrow(checkpoint_failures) == 0) {
      checkpoint_failures <- empty_failures_df()
    }

    save_bootstrap_checkpoint(
      results_df = checkpoint_results,
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
bootstrap_failures <- bind_rows(bootstrap_failures_list)

if (nrow(bootstrap_results) == 0) {
  bootstrap_results <- empty_results_df()
}

if (nrow(bootstrap_failures) == 0) {
  bootstrap_failures <- empty_failures_df()
}

bootstrap_results <- bootstrap_results %>%
  arrange(iteration, Group)

bootstrap_summary_by_group <- bootstrap_results %>%
  group_by(Group) %>%
  summarise(
    successful_iterations = n(),
    mean_peak_age = mean(peak_age),
    sd_peak_age = sd(peak_age),
    mean_peak_pred = mean(peak_pred),
    sd_peak_pred = sd(peak_pred),
    .groups = "drop"
  )

td_bootstrap <- bootstrap_results %>%
  filter(Group == "TD") %>%
  select(iteration, td_peak_age = peak_age, td_peak_pred = peak_pred)

cp_bootstrap <- bootstrap_results %>%
  filter(Group == "CP") %>%
  select(iteration, cp_peak_age = peak_age, cp_peak_pred = peak_pred)

bootstrap_peak_differences <- merge(td_bootstrap, cp_bootstrap, by = "iteration")

if (nrow(bootstrap_peak_differences) > 0) {
  bootstrap_peak_differences$cp_minus_td_peak_age <-
    bootstrap_peak_differences$cp_peak_age - bootstrap_peak_differences$td_peak_age

  bootstrap_difference_summary <- data.frame(
    successful_iterations = nrow(bootstrap_peak_differences),
    mean_cp_minus_td_peak_age = mean(bootstrap_peak_differences$cp_minus_td_peak_age),
    sd_cp_minus_td_peak_age = sd(bootstrap_peak_differences$cp_minus_td_peak_age)
  )
} else {
  bootstrap_difference_summary <- data.frame(
    successful_iterations = 0L,
    mean_cp_minus_td_peak_age = NA_real_,
    sd_cp_minus_td_peak_age = NA_real_
  )
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
  peak_age_difference = bootstrap_difference_summary,
  total_elapsed_seconds = bootstrap_total_elapsed
)

# ──────────────────────────────────────────────
# Save outputs
# ──────────────────────────────────────────────

bootstrap_results_rds <- output_path(paste0(bootstrap_file_stem, "_results.rds"))
bootstrap_results_csv <- output_path(paste0(bootstrap_file_stem, "_results.csv"))
bootstrap_failures_csv <- output_path(paste0(bootstrap_file_stem, "_failures.csv"))
bootstrap_summary_rds <- output_path(paste0(bootstrap_file_stem, "_summary.rds"))
bootstrap_summary_csv <- output_path(paste0(bootstrap_file_stem, "_summary_by_group.csv"))
bootstrap_difference_csv <- output_path(paste0(bootstrap_file_stem, "_peak_age_difference.csv"))

saveRDS(
  list(
    results = bootstrap_results,
    failures = bootstrap_failures,
    summary = bootstrap_summary,
    peak_age_difference = bootstrap_peak_differences
  ),
  bootstrap_results_rds
)

write.csv(bootstrap_results, bootstrap_results_csv, row.names = FALSE)
write.csv(bootstrap_failures, bootstrap_failures_csv, row.names = FALSE)
saveRDS(bootstrap_summary, bootstrap_summary_rds)
write.csv(bootstrap_summary_by_group, bootstrap_summary_csv, row.names = FALSE)
write.csv(bootstrap_difference_summary, bootstrap_difference_csv, row.names = FALSE)

# ──────────────────────────────────────────────
# Print final summary
# ──────────────────────────────────────────────

cat("\nBOOTSTRAP PEAK-AGE SUMMARY\n")
cat(sprintf("Successful iterations: %d\n", successful_bootstrap_fits))
cat(sprintf("Failed iterations: %d\n", nrow(bootstrap_failures)))

td_summary <- bootstrap_summary_by_group %>% filter(Group == "TD")
cp_summary <- bootstrap_summary_by_group %>% filter(Group == "CP")

if (nrow(td_summary) == 1) {
  cat(
    sprintf(
      "TD peak age: mean = %.2f, SD = %.2f\n",
      td_summary$mean_peak_age,
      td_summary$sd_peak_age
    )
  )
}

if (nrow(cp_summary) == 1) {
  cat(
    sprintf(
      "CP peak age: mean = %.2f, SD = %.2f\n",
      cp_summary$mean_peak_age,
      cp_summary$sd_peak_age
    )
  )
}

cat(
  sprintf(
    "Peak age difference (CP - TD): mean = %.2f\n",
    bootstrap_difference_summary$mean_cp_minus_td_peak_age
  )
)
cat(sprintf("Bootstrap total elapsed time: %.1f s\n", bootstrap_total_elapsed))
cat(sprintf("Saved bootstrap results to:\n%s\n", bootstrap_results_rds))
