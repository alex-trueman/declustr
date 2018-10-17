#' Calculate decluster weights for 2D samples using tessellations.
#'
#' @param x Numeric vector. x-axis coordinates with length equal to `y`.
#' @param y Numeric vector. x-axis coordinates with length equal to `x`.
#' @param mask Numeric dataframe. x and y coordinates of grid mask, This
#'     limits the area of tessellation and hence the magnitude of edge weights.
#' @param expand_mask Scalar number (default `0`) of units to exapand mask.
#'     Sometimes when samples are on the mask boundary it is necessary to
#'     slightly expand the mask so that edge samples are included in the
#'     weighting process. Usually only a small value (e.g., 1 unit) is needed.
#' @param normalize Boolean (default `TRUE`). Normalize the weights so that they
#'     sum to 1. Otherwise weights sum to the area of the mask.
#' @param estdomain Boolean (default `TRUE`). Estimate the spatial domain using
#'     a Ripley-Rasson estimate (`spatstat::ripras`) and intersect with the input
#'     mask to create a window for tesselation. This approach may reduce edge
#'     weights where the input mask is extrapolated well beyond the limits of
#'     the data.
#'
#' @return A labelled list with a dataframe of positive weights optionally
#'    normalised such that the sum of weights is 1. The weights are initially
#'    calculated as the ratio of the tesellated area of influence of the sample
#'    divided by the total mask area. Duplicate samples are assigned a zero
#'    weight. The list also contains the point pattern object.
#' @export
#'
#' @importFrom spatstat area dirichletAreas dirichletWeights expand.owin
#'     intersect.owin owin ppp ripras
#'
#' @examples
polydeclust2d <- function(x, y, mask, expand_mask = 0, normalize = TRUE,
                          estdomain = TRUE) {

  # Create mask-only ppp for distance functions.
  points_mask_only <- ppp(x, y, window = owin(mask = mask))

  # No point declustering if n = 1. Don't need to expand mask either.
  if (length(x) == 1) {
    w <- owin(mask = mask)
    points <- ppp(x, y, window = w)
    if (normalize) {
      return(list(weights = data.frame(weight = 1), ppp = points,
             ppp_mask = points_mask_only))
    } else {
      return(list(weights = data.frame(weight = area(points)), ppp = points,
             ppp_mask = points_mask_only))
    }
  }

  # Make mask to control tessellation.
  if (estdomain) {
    w0 <- intersect.owin(owin(mask = mask), ripras(data.frame(x = x, y = y)))
  } else {
    w0 <- owin(mask = mask)
  }
  w <- expand.owin(w0, distance = expand_mask)

  # Calculate weights.
  points <- ppp(x, y, window = w)
  weights <- dirichletWeights(points)
  if (normalize) weights <- weights/sum(dirichletAreas(points))

  return(list(weights = data.frame(weight = weights), ppp = points,
              ppp_mask = points_mask_only))
}

#' Create a plotable tesselation dataframe from a point pattern.
#'
#' @param ppp 2D point pattern object created by package `spatstat`.
#' @param raster_out Boolean (default `FALSE`). Output the tesselation as a point
#'     grid rather than as a path-type dataframe.
#'
#' @return A dataframe with x and y coordinates and group variable for each
#'     tessellation path. May be plotted with `ggplot2` using `geom_path` (or
#'     `geom_raster` in the case of `raster_out = TRUE`.
#' @export
#'
#' @importFrom dplyr arrange mutate select
#' @importFrom ggplot2 fortify
#' @importFrom magrittr %>%
#' @importFrom raster rasterFromXYZ rasterToPolygons
#' @importFrom rlang .data
#' @importFrom spatstat as.data.frame.owin as.data.frame.tess dirichlet npoints
#'     Window
#'
#' @examples
plotable_tess <- function(ppp, raster_out = FALSE) {

  # If only one point there will be an error with `dirichlet`, so just rasterize
  # the window of the ppp.
  if (npoints(ppp) == 1) {
    df_raster <- as.data.frame.owin(Window(ppp)) %>%
      mutate(group = 1) %>%
      select(.data$x, .data$y, .data$group)
  } else {
    # Convert ppp tessellation to ratser dataframe.
    df_raster <- as.data.frame.tess(dirichlet(ppp)) %>%
      arrange(.data$Tile, .data$x, .data$y) %>%
      mutate(group = .data$Tile) %>%
      select(.data$x, .data$y, .data$group)
  }

  if (raster_out) return(df_raster)

  rast <- rasterFromXYZ(df_raster)
  polyg <- rasterToPolygons(rast, dissolve = TRUE)
  df_path <- fortify(polyg, group) %>%
    mutate(x = .data$long, y = .data$lat) %>%
    select(.data$group, .data$x, .data$y)

  return(df_path)
}

#' Calculate 2D point sample spacing.
#'
#' @param points A point pattern object from package `spatstat`.
#' @param nk Scalar number of nearest points to use.
#'
#' @return A dataframe with columns `x`, `y`, and `distance` defining a grid of
#'     points (as defined by the `points` window). This can be used in `ggplot`
#'     with `geom_raster` to vizualize point spacing. Distance is the mean
#'     distance to the `kn` nearest samples converted to a square grid sapcing:
#'         `distance = 2 * sqrt(mean_distance^2 / 2)``
#'     If there are less than `kn` points in the point pattern the distance
#'     is calculated on the available samples.
#' @export
#'
#' @importFrom dplyr left_join mutate rename select
#' @importFrom magrittr %>% %<>%
#' @importFrom rlang := !! .data
#' @importFrom spatstat as.data.frame.im as.im distfun
#'
#' @examples
point_spacing_2d <- function(points, nk = 4) {
  value_cols <- c(3:6)
  for (i in 1:nk) {
    dfun <- as.data.frame.im(as.im(distfun(points, k = i)))
    if (i == 1) {
      dmap <- dfun %>%
        rename(k1 = .data$value)
    } else {
      # If there are less than kn samples in the point pattern break.
      if (is.infinite(dfun$value[1])) {
        value_cols <- c(3:(2 + (i - 1)))
        break
      }
      colname <- paste0("k", i)
      dmap %<>%
        left_join(dfun, by = c("x", "y")) %>%
        rename(`:=`(!!colname, .data$value))
    }
  }
    dmap %<>%
      mutate(distance = 2 * sqrt((rowMeans(select(.data, c(value_cols)),
                                           na.rm = TRUE))^2/2)) %>%
      select(.data$x, .data$y, .data$distance)

    return(dmap)
}
