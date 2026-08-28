#' Information Imbalance
#'
#' @param mx Numeric matrix of hypothesised driving variable measurements.
#' @param my Numeric matrix of hypothesised response variable measurements.
#' @param alpha (optional) Scaling parameter weighting the putative driver measurements.
#' @param lib (optional) Library indices.
#' @param pred (optional) Prediction indices.
#' @param h (optional) Prediction horizon.
#' @param k (optional) Number of nearest neighbors when estimating ranks.
#' @param threads (optional) Number of parallel threads.
#' @param method (optional) Distance measure to be used: `"euclidean"`, `"manhattan"`, or `maximum"`.
#'
#' @returns A numeric vector.
#' @export
#'
#' @examples
#' set.seed(42)
#' mx = embed(rnorm(100), 3)
#' my = embed(rnorm(100), 3)
#' infoxtr::infoimbalance(mx, my)
#' 
infoimbalance = \(mx, my, alpha = seq(0,1,0.1), lib = NULL, pred = NULL,
                  h = 1, k = 3, threads = 1, method = "euclidean") {
  if (is.null(lib)) lib = seq_len(nrow(mx))
  if (is.null(pred)) pred = lib
  return(RcppInfoImbalance(as.matrix(mx), as.matrix(my), abs(alpha), abs(lib), abs(pred),
                           abs(h), abs(k), abs(threads), method))
}
