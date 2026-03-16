# Load required packages
library(glmmTMB)
library(dplyr)
library(emmeans)
library(DHARMa)
library(marginaleffects)
library(interactions)
library(ggeffects)
library(ggplot2)
library(broom.mixed)
library(sjPlot)

# ── Shared helpers ─────────────────────────────────────────────
.root <- if (file.exists("R/paths.R")) "." else ".."
source(file.path(.root, "R/paths.R"))
source(file.path(.root, "R/theme_pub.R"))
source(file.path(.root, "R/setup_fonts.R"))
setup_fonts()

# 1. Read & preprocess
df <- read.csv(data_path("full_data_cp.csv")) %>%
  filter(Range != 0)


df <- df %>%
  mutate(
    Regression     = factor(Regression,     levels = c("biasm","biasp")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian","Asian"), labels = c("Own-Race", "Other-Race")),
    Group     = factor(Group,     levels = c("TD","CP")),
    Subject        = factor(Subject),
  )

# center & scale Range in one step
df$Range_c <- scale(df$Range, center = TRUE, scale = TRUE)


m_simple <- glmmTMB(
  ACC ~ 
    poly(Age, 2) * Regression  * ExperimentName * Group 
  + Range_c
  + (1 | Subject),                 
  data   = df,
  family = binomial(link = "probit")
)


# ──────────────────────────────────────────────
# Predictions from GLMM model
# ──────────────────────────────────────────────
p <- plot_model(
  m_simple,
  type = "pred",
  terms = c("Age [19:79]", "Regression", "ExperimentName", "Group"),
  transform = "response",
  grid.type = "facet"
)

# ──────────────────────────────────────────────
# Adjust colors and labels
# ──────────────────────────────────────────────
final_plot <- p +
  labs(
    title = NULL,
    y = "Accuracy",
    x = "Age"
  ) +
  scale_color_manual(
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias−")
  ) +
  scale_fill_manual(
    values = c(biasp = "#7FB3D5", biasm = "#F08080"),
    labels = c(biasp = "Bias+", biasm = "Bias−")
  ) +
  scale_y_continuous(
    limits = c(0.15, 1),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  guides(fill = "none") +
  scale_x_continuous(
    limits = c(19, 79),
    breaks = seq(20, 80, by = 20)
  ) +
  theme_pub() 

# ──────────────────────────────────────────────
# Save plot
# ──────────────────────────────────────────────
ggsave(
  filename = plot_path("glmm_cp_.png"),
  plot = final_plot,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)




#-------LOSO test

# ---- metrics ----
brier <- function(y, p) mean((y - p)^2)
logloss <- function(y, p, eps = 1e-15) {
  p <- pmin(pmax(p, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

# ---- LOSO helper that also returns Group ----
loso_once <- function(test_id) {
  train_ids <- setdiff(unique(df$Subject), test_id)
  train_df  <- df %>% dplyr::filter(Subject %in% train_ids)
  test_df   <- df %>% dplyr::filter(Subject == test_id)
  
  # Get the held-out subject's Group value
  gvals <- unique(test_df$Group)
  
  fit <- update(m_simple, data = train_df)
  
  # population-level predictions (no RE)
  p <- predict(fit, newdata = test_df, type = "response", re.form = NA)
  
  data.frame(
    Subject = test_id,
    Group   = gvals,
    Brier   = brier(test_df$ACC, p),
    LogLoss = logloss(test_df$ACC, p),
    stringsAsFactors = FALSE
  )
}

# # ---- run LOSO across all subjects ----
# set.seed(1)
# subjects <- unique(df$Subject)
# loso_results <- do.call(rbind, lapply(subjects, loso_once))

# # nice ordering for readability
# # loso_results <- loso_results %>% dplyr::arrange(Group, Subject)

# # ---- print (full and preview) ----
# print(head(loso_results, 10))   # quick peek
# print(loso_results)             # full per-subject table with Group

# # ---- save to CSV ----
# out_path <- output_path("GLMM_LOSO_results.csv")
# write.csv(loso_results, out_path, row.names = FALSE)
# cat("\nSaved per-subject metrics with Group to:\n", out_path, "\n")

# # ---- (optional) group-wise summary, if you want it too ----
# # dplyr::summarise will give you per-Group means & SEs
# loso_by_group <- loso_results %>%
#   dplyr::group_by(Group) %>%
#   dplyr::summarise(
#     n = dplyr::n(),
#     Brier_mean   = mean(Brier),
#     Brier_se     = sd(Brier)/sqrt(n),
#     LogLoss_mean = mean(LogLoss),
#     LogLoss_se   = sd(LogLoss)/sqrt(n),
#     .groups = "drop"
#   )
# print(loso_by_group)


# loso_overall <- loso_results %>%
#   dplyr::summarise(
#     n = dplyr::n(),
#     Brier_mean = mean(Brier),  Brier_se = sd(Brier)/sqrt(n),
#     LogLoss_mean = mean(LogLoss), LogLoss_se = sd(LogLoss)/sqrt(n)
#   )

# cat(sprintf("LOSO (n=%d): Brier = %.3f (SE = %.3f), LogLoss = %.3f (SE = %.3f)\n",
#             loso_overall$n, loso_overall$Brier_mean, loso_overall$Brier_se,
#             loso_overall$LogLoss_mean, loso_overall$LogLoss_se))
