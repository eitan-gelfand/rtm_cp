#' Configure Font for Plot Rendering
#'
#' Sets up Times New Roman font for high-quality plot exports.
#' Automatically falls back to serif if Times New Roman is not available.
#' Compatible across different systems (macOS, Windows, Linux).
#'
#' @param dpi Dots per inch for plot rendering (default: 300)
#' @return NULL (called for side effects)
#' @examples
#' setup_fonts()
#' setup_fonts(dpi = 600)
setup_fonts <- function(dpi = 300) {
  if (!requireNamespace("showtext", quietly = TRUE) ||
      !requireNamespace("sysfonts", quietly = TRUE)) {
    return(invisible(NULL))
  }
  
  # Define font path (macOS system location)
  font_path <- "/System/Library/Fonts/Supplemental/Times New Roman.ttf"
  
  # Add font with fallback to serif if not found
  font_target <- if (file.exists(font_path)) font_path else "serif"
  try(
    sysfonts::font_add("Times New Roman", regular = font_target),
    silent = TRUE
  )
  
  # Configure showtext settings
  showtext::showtext_opts(dpi = dpi)
  showtext::showtext_auto()
  
  invisible(NULL)
}
