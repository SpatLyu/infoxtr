.surd_ts = \(data, target, agent, lag = 1, bin = 5, method = "equal",
             max.combs = 3, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method))
}

.surd_lattice = \(data, target, agent, lag = 1, bin = 5, method = "equal", 
                  max.combs = 3, threads = 1, base = 2.0, normalize = TRUE, nb = NULL) {
  if (is.null(nb)) nb = sdsfun::spdep_nb(data)
  mat = .convert2mat(data)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method, nb))
}

.surd_grid = \(data, target, agent, lag = 1, bin = 5, method = "equal",
               max.combs = 3, threads = 1, base = 2.0, normalize = TRUE) {
  mat = .convert2mat(data)
  return(RcppSURD(mat, target, agent, lag, bin, max.combs, 
                  threads, base, normalize, method, NULL, terra::nrow(data[[1]])))
}

#' SURD
#' 
#' Synergistic-Unique-Redundant Decomposition of causality
#'
#' @inheritParams te
#' @param lag (optional) Lag of the agent variables.
#' @param bin (optional) Number of discretization bins.
#' @param method (optional) Discretization method. One of
#'   `"sd"`, `"equal"`, `"geometric"`, `"quantile"`,
#'   `"natural("jenks")"`, or `"headtail"("headtails")`.
#' @param max.combs (optional) Maximum combination order.
#' @param threads (optional) Number of threads used.
#' @param nb (optional) Neighbours list.
#'
#' @return A list.
#' \describe{
#'   \item{unique}{Unique information contributions per variable.}
#'   \item{synergistic}{Synergistic information components by agent combinations.}
#'   \item{redundant}{Redundant information shared by agent subsets.}
#'   \item{mutual_info}{Mutual information measures for each combination.}
#'   \item{info_leak}{Information leak ratio.}
#' }
#'
#' @export
#' @name surd
#' @aliases surd,data.frame-method
#' @references
#' Martinez-Sanchez, A., Arranz, G. & Lozano-Duran, A. Decomposing causality into its synergistic, unique, and redundant components. Nat Commun 15, 9296 (2024).
#'
#' @examples
#' columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
#' \donttest{
#' tryCatch(
#'   infoxtr::surd(columbus, "hoval", c("inc", "crime")),
#'   error = \(e) message("Skipping Python-dependent example: ", e$message)
#' )
#' }
methods::setMethod("surd", "data.frame", .surd_ts)

#' @rdname surd
methods::setMethod("surd", "sf", .surd_lattice)

#' @rdname surd
methods::setMethod("surd", "SpatRaster", .surd_grid)
