#' Configure a safe serif font for plot rendering.
#'
#' Chooses the first available Times-like system serif so plots render
#' consistently across operating systems without relying on showtext.
#'
#' @param dpi Unused. Kept for compatibility with existing scripts.
#' @return NULL (side effects only)
setup_fonts <- function(dpi = 300) {
  options(project.base_family = "serif")

  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    return(invisible(NULL))
  }

  serif_candidates <- c(
    "Times New Roman",
    "Times",
    "Baskerville",
    "Georgia",
    "Liberation Serif",
    "Nimbus Roman",
    "DejaVu Serif"
  )

  for (family in serif_candidates) {
    font_match <- tryCatch(systemfonts::match_fonts(family), error = function(...) NULL)
    if (!is.null(font_match) && nrow(font_match) > 0 && nzchar(font_match$path[[1]])) {
      options(project.base_family = family)
      break
    }
  }

  invisible(NULL)
}