#' Calculate Decluster Weights for 2D Samples Using Tessellations.
#'
#' \code{polydeclust2d} calcuates decluster weights for 2D point sample data
#' using Voronoi tesselations of the point data. Weights are calculated as the
#' area of influence of each sample divided by the total area of the domain and
#' then normalized so that the sum of all weights is 1. This normalization is
#' the default but can be turned off with the \code{normalize} argument.
#'
#' @param x,y Numeric vectors. x and y axis coordinates with equal length.
#'   Duplicate samples are assigned a zero weight.
#' @param mask Numeric dataframe. x and y coordinates of grid mask. This limits
#'   the area of tessellation and hence the magnitude of edge weights. The
#'   dataframe must have only two columns with x and y coordinates in that order
#'   (the name of the columns is not important).
#' @param expand_mask Scalar number (default 0). Units to exapand the mask.
#'   Sometimes when samples are on the mask boundary it is necessary to slightly
#'   expand the mask so that edge samples are included in the weighting process.
#'   Usually only a small value (e.g., 1 unit) is needed.
#' @param normalize Boolean (default \code{TRUE}). Normalize the weights so that
#'   they sum to 1. Otherwise weights sum to the area of the mask.
#' @param estdomain Boolean (default \code{TRUE}). Estimate the spatial domain
#'   using a Ripley-Rasson estimate (\code{\link{ripras}}) and
#'   intersect with the input mask to create a window for tesselation. This
#'   approach may reduce edge weights where the input mask is extrapolated well
#'   beyond the limits of the data.
#' @return A labelled list with a named vector of positive weights optionally
#'   normalized such that the sum of weights is 1. The list also contains the
#'   point pattern object used for the tessellation and a point pattern object
#'   with only the mask (no ripras).
#' @export
#' @importFrom assertthat are_equal assert_that is.scalar not_empty noNA
#' @importFrom spatstat area dirichletAreas dirichletWeights expand.owin
#'   inside.owin intersect.owin npoints owin ppp ripras
#' @examples
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask)
#' samples_dec <- cbind(samples, dec$weights)
polydeclust2d <- function(x, y, mask, expand_mask = 0, normalize = TRUE,
                          estdomain = TRUE) {

  # Assertion of function arguments.
  # `x` and `y` are non-empty, numeric vectors of equal length with no NA values.
  assert_that(is.numeric(x))
  assert_that(is.vector(x))
  assert_that(not_empty(x))
  assert_that(noNA(x))
  assert_that(is.numeric(y))
  assert_that(is.vector(y))
  assert_that(not_empty(y))
  assert_that(noNA(y))
  assert_that(are_equal(length(x), length(y)),
    msg = "x and y vectors must have equal length")
  # `mask` is a dataframe with two numeric columns and no NA values.
  assert_that(is.data.frame(mask))
  assert_that(are_equal(ncol(mask), 2), msg = "`mask` must have two columns with
      x and y coordinates, respectively.")
  assert_that(is.numeric(mask[[1]]))
  assert_that(is.numeric(mask[[2]]))
  # `expand_mask` must be a scalar number >= 0.
  assert_that(is.scalar(expand_mask))
  assert_that(is.numeric(expand_mask))
  assert_that(expand_mask >= 0, msg = "`expand_mask` must be >= 0")
  # `normalize` and `estdomain` are booleans.
  assert_that(is.logical(normalize))
  assert_that(is.logical(estdomain))

  # Create a mask for the whole domain as defined by the `mask` argument.
  domain_mask <- owin(mask = mask)

  # Create mask-only ppp for distance function, etc.
  points_mask_only <- ppp(x, y, window = domain_mask)

  # No point declustering if n = 1. Don't need to expand mask either.
  # Just create  the required objects and set the weight.
  if (length(x) == 1) {
    points <- ppp(x, y, window = domain_mask)
    poutside <- c()
    if (normalize) {
      return(list(weights = data.frame(weight = 1), ppp = points,
             ppp_mask = points_mask_only))
    } else {
      return(list(weights = data.frame(weight = area(points)), ppp = points,
             ppp_mask = points_mask_only))
    }
  }

  # Make window to control tessellation.
  if (estdomain) {
    w0 <- intersect.owin(domain_mask, ripras(data.frame(x = x, y = y)))
  } else {
    w0 <- domain_mask
  }
  w <- expand.owin(w0, distance = expand_mask)

  # Identify points outside of defined window.
  poutside <- inside.owin(x, y, w = w)

  # Calculate weights.
  points <- ppp(x, y, window = w)
  weights <- c(dirichletWeights(points))
  if (normalize) weights <- weights/sum(dirichletAreas(points))

  # Assertions on return values.
  # `weights` is a numeric vector with no NA values having equal length to the
  # input x and y vectors.
  assert_that(is.numeric(weights))
  assert_that(is.vector(weights))
  assert_that(noNA(weights))
  assert_that(are_equal(length(weights), length(x)))
  # Objects `points` and `points_mask_only` should be of class `ppp` and should
  # have the same number of points as the input x and y vectors.
  assert_that(class(points) == "ppp", msg = "retuned point pattern is not of
    class 'ppp'")
  assert_that(class(points_mask_only) == "ppp", msg = "retuned point pattern is
    not of class 'ppp'")
  assert_that(are_equal(length(x), npoints(points)))
  assert_that(are_equal(length(x), npoints(points_mask_only)))

  return(list(weights = data.frame(weight = weights), ppp = points,
              ppp_mask = points_mask_only, outside_point = poutside))
}

#' Create Plotable Tesselation from a Point Pattern.
#'
#' \code{plotable_tess} returns a dataframe of points defining the edges of the
#' tessellation constructed by \code{\link{polydeclust2d}}. These points can be
#' used to visualize the tessellations.
#'
#' @param points 2D point pattern object created by package
#'   \code{\link{spatstat}}. It must have a window defined using a
#'   binary mask.
#' @return A simple features (\code{sf}) polygons dataframe.  May be plotted
#'   with \code{\link{geom_sf}}.
#' @export
#' @importFrom assertthat assert_that
#' @importFrom dplyr arrange bind_rows mutate select
#' @importFrom magrittr %>% %<>%
#' @importFrom raster rasterFromXYZ rasterToPolygons
#' @importFrom rlang .data
#' @importFrom sf st_as_sf
#' @importFrom spatstat as.data.frame.owin as.data.frame.tess dirichlet is.mask
#'   npoints Window
#' @examples
#' library(ggplot2)
#'
#' # Tessellation with a mask only may produce excessive edge weights.
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask, estdomain = FALSE)
#' samples_dec <- cbind(samples, dec$weights)
#' ptess <- plotable_tess(dec$ppp)
#' ggplot() +
#'   geom_raster(data = mask, aes(x, y), fill = "lightblue") +
#'   geom_sf(data = ptess, fill = NA) +
#'   geom_point(data = samples_dec, aes(x, y, size = weight))
#'
#' # Using the `ripras` option reduces excessive edge weights.
#' dec <- polydeclust2d(samples$x, samples$y, mask = mask)
#' samples_dec <- cbind(samples, dec$weights)
#' ptess <- plotable_tess(dec$ppp)
#' ggplot() +
#'   geom_raster(data = mask, aes(x, y), fill = "lightblue") +
#'   geom_sf(data = ptess, fill = NA) +
#'   geom_point(data = samples_dec, aes(x, y, size = weight))
plotable_tess <- function(points) {

  # Assertions on input arguments.
  # `point` is a point pattern object with >= 1 point using a binary mask.
  assert_that(class(points) == "ppp",
    msg = "retuned point pattern is not of class 'ppp'")
  assert_that(is.mask(Window(points)),
    msg = "point pattern must be defined with binary mask")
  assert_that(npoints(points) >= 1,
    msg = "point pattern must have >= 1 points")

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

  # Assertions on return value.
  # `sf_tess` should be of class 'sf'
  assert_that("sf" %in% class(sf_tess),
    msg = "return value is not of class 'sf'")

  return(sf_tess)
}

#' Calculate 2D Point Sample Spacing.
#'
#' \code{point_spacing_2d} computes the avererage distance to the four nearest
#' samples in the mask domain. The number of nearest samples can be adjusted.
#' The returned object is a dataframe of regular grid points (as defined by the
#' input point pattern object) with column of distances. Distance is the mean
#' distance to the \code{kn} nearest samples converted to a square grid sapcing:
#' \code{distance = 2 * sqrt(mean_distance^2 / 2)}. If there are less than
#' \code{kn} points in the point pattern, the distance is calculated on the
#' available samples.
#'
#' @param points A point pattern object from package
#'   \code{\link{spatstat}}.
#' @param nk Scalar number (default 4) of nearest points to use (>= 1).
#'
#' @return A dataframe with columns 'x', 'y', and 'distance' defining a grid of
#'   points (as defined by the \code{points} window). This can be used with
#'   \code{\link{geom_raster}} to vizualize point spacing.
#' @export
#'
#' @importFrom assertthat assert_that is.scalar
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
point_spacing_2d <- function(points, nk = 4L) {

  # Assertions on input arguments.
  # `point` is a point pattern object with >= 1 point using a binary mask.
  assert_that(class(points) == "ppp",
  msg = "retuned point pattern is not of class 'ppp'")
  assert_that(is.mask(Window(points)),
  msg = "point pattern must be defined with binary mask")
  assert_that(npoints(points) >= 1,
    msg = "point pattern must have >= 1 points")
  # `nk` is a scalar whole number >= 1.
  assert_that(is.scalar(nk))
  assert_that(
    (function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol)(nk),
    msg = "`nk` must be a whole number")
  assert_that(nk >= 1, msg = "`nk` must be >= 1")

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
