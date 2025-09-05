#' Plot Case vs Control Methylation Fractions with Density-Weighted Points
#'
#' This function generates a scatter plot of methylation fractions in case vs control samples
#' from a CH3 SQLite database table, with optional density shading and threshold lines.
#'
#' @param ch3_db A path to the CH3 SQLite database, or an open database connection object.
#' @param table A string specifying the table name within the database, containing `mh_frac_case` and `mh_frac_control` columns.
#' @param thresh A numeric vector of length 1 or 2 giving the positive and negative threshold values to draw as dashed red lines above and below the identity line (default is `c(0.2, 0.3)`). If `NULL`, no threshold lines are drawn.
#' @param palette A string specifying the viridis color palette for density shading (default is `"turbo"`). Other options include `"viridis"`, `"plasma"`, `"cividis"`, etc.
#' @param size A numeric value controlling point size in the plot (default is `0.8`).
#'
#' @return Invisibly returns the closed database connection object. The plot is printed to the current graphics device.
#'
#' @details
#' This function visualizes the relationship between `mh_frac_case` and `mh_frac_control` methylation values
#' by plotting them on a scatter plot. It optionally enhances visualization with local point density using
#' the `ggpointdensity` package if available.
#'
#' Features:
#' \itemize{
#'   \item If \strong{ggpointdensity} is installed, point density is used to color points.
#'   \item If not installed, a fallback to plain transparent points is used.
#'   \item A solid black identity line (`y = x`) is added for reference.
#'   \item Optional threshold lines at ±`thresh` distance from the identity line are shown as red dashed lines.
#' }
#'
#' The plot is limited to the [0, 1] range for both axes and uses a clean classic theme.
#'
#' @examples
#' \dontrun{
#' plot_ch3_diff_density("my_data/ch3.sqlite", table = "chr5", thresh = c(0.25, 0.35))
#' }
#'
#' @importFrom DBI dbExistsTable
#' @importFrom dplyr tbl select filter collect
#' @importFrom ggplot2 ggplot aes geom_abline geom_point coord_equal labs theme_classic theme_minimal
#' @importFrom ggplot2 scale_colour_viridis_c guide_colorbar
#' @export

plot_ch3_diff_density <- function(ch3_db,
                                  table,
                                  thresh   = c(0.2, 0.3),
                                  palette  = "turbo",
                                  n_kde    = 250)
{
  
  ch3_db <- .ch3helper_connectDB(ch3_db)
  on.exit(ch3_db <- .ch3helper_closeDB(ch3_db), add = TRUE)
  
  if (!DBI::dbExistsTable(ch3_db$con, table))
    stop(table, " table does not exist.")
  
  df <- dplyr::tbl(ch3_db$con, table) |>
    dplyr::select(mh_frac_case, mh_frac_control) |>
    dplyr::filter(!is.na(mh_frac_case), !is.na(mh_frac_control)) |>
    dplyr::collect()
  
  if (nrow(df) < 20) stop("Too few points to draw a density surface.")
  
  # ── KDE surface  ---------------------------------------------------------
  kd <- MASS::kde2d(df$mh_frac_case, df$mh_frac_control,
                    n    = n_kde,
                    lims = c(0, 1, 0, 1))
  
  surf <- expand.grid(x = kd$x, y = kd$y)
  surf$z <- as.vector(kd$z)
  
  # ── base plot  -----------------------------------------------------------
  p <- ggplot2::ggplot(surf, ggplot2::aes(x, y, fill = z)) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(
      option = palette, name = "Density", guide = "colourbar"
    ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "black")
  
  p <- p +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      title = "Case vs Control Mod Fractions",
      x = "mh_frac_case",
      y = "mh_frac_control"
    ) +
    ggplot2::theme_minimal()
  
  print(p)
  invisible(ch3_db)
}
