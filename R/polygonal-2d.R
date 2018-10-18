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
#'    weight. The list also contains the point pattern object used for the
#'    tessellation and a point pattern object with only the mask (no ripras).
#' @export
#'
#' @importFrom spatstat area dirichletAreas dirichletWeights expand.owin
#'     intersect.owin owin ppp ripras
#'
#' @examples
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask)
#' samples_dec <- cbind(samples, dec$weights)
polydeclust2d <- function(x, y, mask, expand_mask = 0, normalize = TRUE,
                          estdomain = TRUE) {

  # Create a mask for the whole domain as defined by the `mask` argument.
  domain_mask <- owin(mask = mask)

  # Create mask-only ppp for distance function, etc.
  points_mask_only <- ppp(x, y, window = domain_mask)

  # No point declustering if n = 1. Don't need to expand mask either.
  # Just create  the required objects and set the weight.
  if (length(x) == 1) {
    points <- ppp(x, y, window = domain_mask)
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
    w0 <- intersect.owin(domain_mask, ripras(data.frame(x = x, y = y)))
  } else {
    w0 <- domain_mask
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
#' @param points 2D point pattern object created by package `spatstat`.
#'
#' @return A simple features (sf) polygons dataframe.  May be plotted with
#'     `ggplot2` using `geom_sf`.
#' @export
#'
#' @importFrom dplyr arrange bind_rows mutate select
#' @importFrom magrittr %>% %<>%
#' @importFrom raster rasterFromXYZ rasterToPolygons
#' @importFrom rlang .data
#' @importFrom sf st_as_sf
#' @importFrom spatstat as.data.frame.owin as.data.frame.tess dirichlet npoints
#'     Window
#'
#' @examples
#' library(ggplot2)
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask)
#' samples_dec <- cbind(samples, dec$weights)
#' ptess <- plotable_tess(dec$ppp)
#' ggplot() +
#'   geom_raster(data = mask, aes(x, y), fill = "lightblue") +
#'   geom_sf(data = ptess, fill = NA) +
#'   geom_point(data = samples_dec, aes(x, y, size = weight))
plotable_tess <- function(points) {

  # If only one point there will be an error with `dirichlet`, so just
  # return the polygon of the whole domain.
  if (npoints(points) == 1) {
    # Make a dataframe of the window of the ppp.
    df_tess <- as.data.frame.owin(Window(points))
    # If there are multiple regions there will be an id column.
    if("id" %in% colnames(df_tess)) {
      df_tess %<>%
        mutate(group = as.numeric(.data$id)) %>%
        select(.data$x, .data$y, .data$group)
    } else {
      df_tess %<>%
        mutate(group = 1) %>%
        select(.data$x, .data$y, .data$group)
    }
  } else {
    df_tess <- as.data.frame.tess(dirichlet(points)) %>%
      mutate(group = as.numeric(.data$Tile)) %>%
      select(.data$x, .data$y, .data$group)
  }

  # Convert ppp tessellation to dataframe.
  sf_tess <- df_tess %>%
    rasterFromXYZ() %>%
    rasterToPolygons(dissolve = TRUE) %>%
    st_as_sf()

  return(sf_tess)
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
#' library(ggplot2)
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask)
#' spacing <- point_spacing_2d(dec$ppp_mask)
#' ggplot() +
#'   geom_raster(data = spacing, aes(x, y, fill = distance)) +
#'   geom_point(data = samples, aes(x, y)) +
#'   scale_fill_viridis_c()
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
        rename(!!colname := .data$value)
    }
  }
    dmap %<>%
      mutate(distance = 2 * sqrt((rowMeans(select(., c(value_cols)),
                                           na.rm = TRUE))^2/2)) %>%
      select(.data$x, .data$y, .data$distance)

    return(dmap)
}
