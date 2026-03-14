# R/theme_pub.R
# Shared publication-style ggplot theme used across analysis scripts.

# Returns a consistent black-and-white base theme with project defaults.
theme_pub <- function(base_size = 14, base_family = "Times New Roman") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      axis.text        = element_text(size = base_size - 2),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title       = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      strip.text       = element_text(size = base_size, face = "plain"),
      strip.background = element_rect(fill = "gray90", linewidth = 0.1),
      legend.position  = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.4)
    )
}