# Load required packages
library(glmmTMB)
library(dplyr)
library(emmeans)
library(ggplot2)
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
df$Age_c <- scale(df$Age, center = TRUE, scale = FALSE)


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
# Predictions from GLMM model
# ──────────────────────────────────────────────
p <- plot_model(
  m_simple,
  type = "pred",
  terms = c("Age [19:79]", "Regression", "Group", "ExperimentName"),
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


# ──────────────────────────────────────────────
# Simple effects analysis
# ──────────────────────────────────────────────

# 2. Follow-up model
# better interpret the Age predictor by including seperatly linear and quadratic terms 

m_followup <- glmmTMB(
  ACC ~ (Age_c + I(Age_c^2)) * Regression * ExperimentName * Group +
    Range_c + (1 | Subject),
  data   = df,
  family = binomial(link = "probit")
)

# 3. Simple categorical effects within each group
race_by_group <- emmeans(m_followup, pairwise ~ ExperimentName | Group)
reg_by_group  <- emmeans(m_followup, pairwise ~ Regression | Group)

# 4. Simple linear age effects within each group
age_race_by_group <- emtrends(m_followup, pairwise ~ ExperimentName | Group, var = "Age_c")
age_reg_by_group  <- emtrends(m_followup, pairwise ~ Regression | Group, var = "Age_c")



# 6. Print outputs
cat("\nRACE EFFECT WITHIN EACH GROUP\n")
print(race_by_group$contrasts)

cat("\nREGRESSION EFFECT WITHIN EACH GROUP\n")
print(reg_by_group$contrasts)

cat("\nLINEAR AGE × RACE WITHIN EACH GROUP\n")
print(age_race_by_group$contrasts)

cat("\nLINEAR AGE × REGRESSION WITHIN EACH GROUP\n")
print(age_reg_by_group$contrasts)



# ──────────────────────────────────────────────
# Peak age from full GLMM predictions
# Overall by Group + by Group×Race×Regression
# ──────────────────────────────────────────────


# 1. Dense prediction grid
pred_grid <- expand.grid(
  Age = seq(19, 78, by = 0.1),
  Regression = levels(df$Regression),
  ExperimentName = levels(df$ExperimentName),
  Group = levels(df$Group),
  Range_c = 0
)

# 2. Predict from full model
pred_grid$pred <- predict(
  m_simple,
  newdata = pred_grid,
  type = "response",
  re.form = NA
)

# ── Output 1: one peak per whole group ───────────────────────
peak_input_group <- pred_grid %>%
  group_by(Group, Age) %>%
  summarise(
    pred_mean = mean(pred),
    .groups = "drop"
  )

peak_by_group <- peak_input_group %>%
  group_by(Group) %>%
  slice_max(order_by = pred_mean, n = 1, with_ties = FALSE) %>%
  ungroup()

cat("\nPEAK BY WHOLE GROUP\n")
print(peak_by_group)

# ── Output 2: one peak per Group × Race × Regression ─────────
peak_by_condition <- pred_grid %>%
  group_by(Group, ExperimentName, Regression) %>%
  slice_max(order_by = pred, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Group, ExperimentName, Regression)

cat("\nPEAK BY GROUP × RACE × REGRESSION\n")
print(peak_by_condition)


# ── Output 3: one peak per Group × Race ─────────
peak_by_group_race <- pred_grid %>%
  group_by(Group, ExperimentName, Age) %>%
  summarise(
    pred_mean = mean(pred),
    .groups = "drop"
  ) %>%
  group_by(Group, ExperimentName) %>%
  slice_max(order_by = pred_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Group, ExperimentName)

cat("\nPEAK BY GROUP × RACE\n")
print(peak_by_group_race)




# =========================================================
# Mean ages at subject level
# =========================================================

subject_age <- df %>%
  group_by(Subject, Group) %>%
  summarise(age = mean(Age, na.rm = TRUE), .groups = "drop")

young_mean_age <- subject_age %>%
  filter(age <= 42) %>%
  summarise(mean_age = mean(age, na.rm = TRUE)) %>%
  pull(mean_age)

old_mean_age <- subject_age %>%
  filter(age > 42) %>%
  summarise(mean_age = mean(age, na.rm = TRUE)) %>%
  pull(mean_age)

print(young_mean_age)
print(old_mean_age)

range_lock <- 0


# =========================================================
# YOUNG-LOCKED BETWEEN-GROUP CONTRASTS
# =========================================================

# 1) Group contrast within Race
emm_group_by_race <- emmeans(
  m_simple,
  ~ Group | ExperimentName,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_group_by_race <- pairs(emm_group_by_race, adjust = "holm")
contrast_group_by_race


# 2) Group contrast within Regression
emm_group_by_reg <- emmeans(
  m_simple,
  ~ Group | Regression,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_group_by_reg <- pairs(emm_group_by_reg, adjust = "holm")
contrast_group_by_reg


# 3) Group contrast within Race × Regression
emm_group_by_race_reg <- emmeans(
  m_simple,
  ~ Group | ExperimentName * Regression,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_group_by_race_reg <- pairs(emm_group_by_race_reg, adjust = "holm")
contrast_group_by_race_reg


# =========================================================
# YOUNG-LOCKED WITHIN-GROUP CONTRASTS
# =========================================================

# 1) Race contrast within Group
emm_race_by_group <- emmeans(
  m_simple,
  ~ ExperimentName | Group,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_race_by_group <- pairs(emm_race_by_group, adjust = "holm")
contrast_race_by_group


# 2) Regression contrast within Group
emm_reg_by_group <- emmeans(
  m_simple,
  ~ Regression | Group,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_reg_by_group <- pairs(emm_reg_by_group, adjust = "holm")
contrast_reg_by_group


# 3) Race contrast within Group × Regression
emm_race_by_group_reg <- emmeans(
  m_simple,
  ~ ExperimentName | Group * Regression,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_race_by_group_reg <- pairs(emm_race_by_group_reg, adjust = "holm")
contrast_race_by_group_reg


# 4) Regression contrast within Group × Race
emm_reg_by_group_race <- emmeans(
  m_simple,
  ~ Regression | Group * ExperimentName,
  at = list(Age = young_mean_age, Range_c = range_lock),
  type = "response"
)

contrast_reg_by_group_race <- pairs(emm_reg_by_group_race, adjust = "holm")
contrast_reg_by_group_race


# =========================================================
# AGE CONTRASTS: YOUNG MEAN vs OLD MEAN
# =========================================================
# This compares the two locked ages directly,
# within Group / Race / Regression conditions

# 1) Age contrast within Group
emm_age_by_group <- emmeans(
  m_simple,
  ~ Age | Group,
  at = list(
    Age = c(young_mean_age, old_mean_age),
    Range_c = range_lock
  ),
  cov.reduce = FALSE,
  type = "response"
)

contrast_age_by_group <- pairs(emm_age_by_group, adjust = "holm")
contrast_age_by_group


# 2) Age contrast within Group × Race
emm_age_by_group_race <- emmeans(
  m_simple,
  ~ Age | Group * ExperimentName,
  at = list(
    Age = c(young_mean_age, old_mean_age),
    Range_c = range_lock
  ),
  cov.reduce = FALSE,
  type = "response"
)

contrast_age_by_group_race <- pairs(emm_age_by_group_race, adjust = "holm")
contrast_age_by_group_race


# 3) Age contrast within Group × Regression
emm_age_by_group_reg <- emmeans(
  m_simple,
  ~ Age | Group * Regression,
  at = list(
    Age = c(young_mean_age, old_mean_age),
    Range_c = range_lock
  ),
  cov.reduce = FALSE,
  type = "response"
)

contrast_age_by_group_reg <- pairs(emm_age_by_group_reg, adjust = "holm")
contrast_age_by_group_reg


# 4) Age contrast within Group × Race × Regression
emm_age_by_group_race_reg <- emmeans(
  m_simple,
  ~ Age | Group * ExperimentName * Regression,
  at = list(
    Age = c(young_mean_age, old_mean_age),
    Range_c = range_lock
  ),
  cov.reduce = FALSE,
  type = "response"
)

contrast_age_by_group_race_reg <- pairs(emm_age_by_group_race_reg, adjust = "holm")
contrast_age_by_group_race_reg




# ──────────────────────────────────────────────
# Nested Model Comparison
# ──────────────────────────────────────────────

m_4_way <- glmmTMB(
  ACC ~ 
    poly(Age, 2) * Regression  * ExperimentName * Group 
  + Range_c
  + (1 | Subject),                 
  data   = df,
  family = binomial(link = "probit")
)
summary(m_4_way)
anova(m_4_way, m_simple)


m_2_way <- glmmTMB(
  ACC ~
    (poly(Age, 2) + Regression + ExperimentName + Group)^2 +
    Range_c +
    (1 | Subject),
  data   = df,
  family = binomial(link = "probit")
)
anova(m_2_way, m_simple)


# ──────────────────────────────────────────────
# LOSO analysis (Leave-One-Subject-Out)
# ──────────────────────────────────────────────
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

# ---- run LOSO across all subjects ----
set.seed(1)
subjects <- unique(df$Subject)
loso_results <- do.call(rbind, lapply(subjects, loso_once))

# nice ordering for readability
# loso_results <- loso_results %>% dplyr::arrange(Group, Subject)

# ---- print (full and preview) ----
print(head(loso_results, 10))   # quick peek
print(loso_results)             # full per-subject table with Group

# ---- save to CSV ----
out_path <- output_path("GLMM_LOSO_results.csv")
write.csv(loso_results, out_path, row.names = FALSE)
cat("\nSaved per-subject metrics with Group to:\n", out_path, "\n")

# ---- (optional) group-wise summary, if you want it too ----
# dplyr::summarise will give you per-Group means & SEs
loso_by_group <- loso_results %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(
    n = dplyr::n(),
    Brier_mean   = mean(Brier),
    Brier_se     = sd(Brier)/sqrt(n),
    LogLoss_mean = mean(LogLoss),
    LogLoss_se   = sd(LogLoss)/sqrt(n),
    .groups = "drop"
  )
print(loso_by_group)


loso_overall <- loso_results %>%
  dplyr::summarise(
    n = dplyr::n(),
    Brier_mean = mean(Brier),  Brier_se = sd(Brier)/sqrt(n),
    LogLoss_mean = mean(LogLoss), LogLoss_se = sd(LogLoss)/sqrt(n)
  )

cat(sprintf("LOSO (n=%d): Brier = %.3f (SE = %.3f), LogLoss = %.3f (SE = %.3f)\n",
            loso_overall$n, loso_overall$Brier_mean, loso_overall$Brier_se,
            loso_overall$LogLoss_mean, loso_overall$LogLoss_se))
