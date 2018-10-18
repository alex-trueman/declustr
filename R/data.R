#' 2D Point Sample Data.
#'
#' A dataset containing sample point data for a 2D (xy) domain. The samples
#'     record a concentration of metal the domain thickness and the accumulation
#'     (metal concentration × thickness).
#'
#' @docType data
#'
#' @usage data(samples)
#'
#' @format A data frame with 44 rows and 6 variables:
#' \describe{
#'   \item{id}{unique point ID}
#'   \item{x}{x coordinate in metres}
#'   \item{y}{y coordinate in metres}
#'   \item{value}{value metal concentration in ppm}
#'   \item{thickness}{domain thickness in metres}
#'   \item{accumulation}{metal accumulation in ppm·m}
#' }
"samples"

#' 2D Sample Domain Mask.
#'
#' A dataset containing 2D (x and y) coordinates of a 1 × 1 metre regular grid
#'     that defines the sampled domain associated with the `samples` dataset.
#'
#' @docType data
#'
#' @usage data(mask)
#'
#' @format A data frame with 21,601 rows and 2 variables:
#' \describe{
#'   \item{x}{x coordinate in metres}
#'   \item{y}{y coordinate in metres}
#' }
"mask"
