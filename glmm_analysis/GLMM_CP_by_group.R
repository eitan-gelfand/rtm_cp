# ──────────────────────────────────────────────
# Load required packages
# ──────────────────────────────────────────────
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

# ──────────────────────────────────────────────
# 1. Read & preprocess
# ──────────────────────────────────────────────
df <- read.csv(data_path("full_data_cp.csv")) %>%
  filter(Range != 0) %>%
  mutate(
    Regression     = factor(Regression,     levels = c("biasm","biasp")),
    ExperimentName = factor(ExperimentName, levels = c("Caucasian","Asian"), 
                            labels = c("Own-Race", "Other-Race")),
    Group          = factor(Group, levels = c("TD","CP")),
    Subject        = factor(Subject)
  )

# Center & scale Range
df$Range_c <- scale(df$Range, center = TRUE, scale = TRUE)

# Create categorical variable of Age
df$AgeGroup <- ifelse(df$Age < 40, "Younger", "Older")

# ──────────────────────────────────────────────
# 2. Split datasets
# ──────────────────────────────────────────────
df_TD <- df %>% filter(Group == "TD")
df_CP <- df %>% filter(Group == "CP")

# ──────────────────────────────────────────────
# 3. TD model
# ──────────────────────────────────────────────
m_TD <- glmmTMB(
  ACC ~ poly(Age, 2) * Regression * ExperimentName + Range_c + (1 | Subject),
  data = df_TD,
  family = binomial(link = "probit")
)

summary(m_TD)

# Plot TD model predictions
p_TD <- plot_model(
  m_TD,
  type = "pred",
  terms = c("Age [19:79]", "Regression", "ExperimentName"),
  transform = "response",
  grid.type = "facet"
) +
  labs(title = "Typical Development (TD)",
       y = "Accuracy",
       x = "Age") +
  scale_color_manual(
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias−")
  ) +
  theme_pub()

ggsave(
  filename = plot_path("glmm_TD_by_group.png"),
  plot = p_TD,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)

# ──────────────────────────────────────────────
# 4. CP model
# ──────────────────────────────────────────────
m_CP <- glmmTMB(
  ACC ~ 
  poly(Age, 2) * Regression 
  + poly(Age, 2) * ExperimentName
  + ExperimentName : Regression
  + Range_c + (1 | Subject),
  data = df_CP,
  family = binomial(link = "probit")
)

summary(m_CP)

# Plot CP model predictions
p_CP <- plot_model(
  m_CP,
  type = "pred",
  terms = c("Age [19:79]", "Regression", "ExperimentName"),
  transform = "response",
  grid.type = "facet"
) +
  labs(title = "Congenital Prosopagnosia (CP)",
       y = "Accuracy",
       x = "Age") +
  scale_color_manual(
    values = c("biasp" = "#7FB3D5", "biasm" = "#F08080"),
    labels = c("biasp" = "Bias+", "biasm" = "Bias−")
  ) +
  theme_pub()

ggsave(
  filename = plot_path("glmm_CP_by_group.png"),
  plot = p_CP,
  width = 8,
  height = 6,
  dpi = 300,
  device = ragg::agg_png
)


